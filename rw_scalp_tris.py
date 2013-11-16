import ipdb

def get_scalp_tris(in_file, original_T1, final_volume, mesh_file, mesh_id, electrode_name_file, basename="scalp_triangles_", session_number=5):
    import scipy.io as sio
    import os.path as op
    import numpy as np
    import nibabel as nb
    struct = sio.loadmat(in_file, squeeze_me=True)
    electrode_names = np.loadtxt(electrode_name_file, dtype=str)
    electrodes = struct['ELECTRODESF']
    position = electrodes[session_number - 1]

    in_t1 = nb.load(original_T1)
    t1_data = in_t1.get_data()
    t1_header = in_t1.get_data()
    out_files = []
    for pt_id, data in enumerate(position):
        x = -float(data[0]) + 112
        y = float(data[2]) - 127.5
        z = float(data[1]) - 127.5
        point = [x, y, z]
        out_basename = op.abspath(basename + str(session_number) + "_" + electrode_names[pt_id])
        out_file = get_elements_near_point(
            point, mesh_file, mesh_id, pt_id+6000, pt_id+5000, out_basename)
        print('Writing electrodes to {ec} to {f}'.format(ec=electrode_names[pt_id], f=out_file))
        out_files.append(out_file)
    return out_files


def get_elements_near_point(point, mesh_filename, mesh_id, new_physical_surface_id, 
        new_elementary_entity_id, out_basename):
    import numpy as np
    from scipy.spatial.distance import cdist
    import os.path as op
    import time
    start_time = time.time()
    mesh_file = open(mesh_filename, 'r')

    while True:
        line = mesh_file.readline()
        if '$Nodes' in line:
            line = mesh_file.readline()
            number_of_nodes = int(line)
            vertices = []

            for i in xrange(0, number_of_nodes):
                line = mesh_file.readline()
                node_data = line.split()
                vert = node_data
                vertices.append(vert)

            vertices = np.array((vertices))
            dist = cdist([point], vertices[:, 1:], 'euclidean')
            min_dist = np.min(dist)
            closest = np.where(dist == min_dist)[1]
            closest_idx = vertices[closest][0][0]
            #print(vertices[closest])
            #print(min_dist)

        elif '$Elements' in line:
            line = mesh_file.readline()
            number_of_elements = int(line)
            polygons = []
            for i in xrange(0, number_of_elements):
                # -- If all elements were quads, each line has 10 numbers.
                #    The first one is a tag, and the last 4 are node numbers.
                line = mesh_file.readline()
                elem_data = line.split()
                if int(elem_data[3]) == mesh_id:
                    polygons.append(elem_data)

            polygons = np.array((polygons))
            poly_indices = polygons[:, -3:]
            rows = np.where(poly_indices == closest_idx)

            # This is the ID of all the polygons connected to the closest
            # vertex
            connected_polygons = polygons[rows[0]][:, 0]
            connected_element_data = polygons[rows[0]]
            break

    mesh_file.close()
    elapsed_time = time.time() - start_time
    verts_to_keep = np.unique(polygons[rows[0]][:, -3:])
    new_idx = 1
    remapping = {}
    old_nodes = []
    new_nodes = []
    new_elements = []
    for vert_id in verts_to_keep:
        idx = np.where(vertices[:, 0] == vert_id)
        new_node = vertices[idx]
        old_node = new_node.copy()
        new_node[0, 0] = new_idx
        remapping[vert_id] = new_idx
        new_idx += 1
        new_nodes.append(new_node)
        old_nodes.append(old_node)

    new_idx = 1
    for elem_id in connected_polygons:
        idx = np.where(polygons[:, 0] == elem_id)
        new_elem = polygons[idx]
        new_elem[0, 0] = new_idx
        new_elem[0, -3] = remapping[new_elem[0, -3]]
        new_elem[0, -2] = remapping[new_elem[0, -2]]
        new_elem[0, -1] = remapping[new_elem[0, -1]]
        new_idx += 1
        new_elements.append(new_elem)

    #print("Old Nodes")
    #print(old_nodes)
    #print("New Nodes")
    #print(new_nodes)

    #print("Old Elements")
    #print(connected_element_data)
    #print("New Elements")
    #print(new_elements)
    out_file = op.abspath(out_basename + ".msh")

    f = open(out_file, 'w')
    f.write("$MeshFormat\n")
    f.write("2.2 0 8\n")
    f.write("$EndMeshFormat\n")
    f.write("$Nodes\n")
    f.write(str(len(new_nodes)) + "\n")

    for idx in xrange(0, len(new_nodes)):
        line = new_nodes[idx][0]
        line = line.tolist()
        f.write(" ".join(line) + "\n")

    f.write("$EndNodes\n")
    f.write("$Elements\n")
    f.write(str(len(new_elements)) + "\n")
    print(new_physical_surface_id)
    print(new_elementary_entity_id)
    for idx in xrange(0, len(new_elements)):
        line = np.empty((8), dtype=int)
        line[0] = idx+1 #element number
        line[1] = 2 # element type (triangle)
        line[2] = 2 # number of tags
        line[3] = new_physical_surface_id # tag
        line[4] = new_elementary_entity_id # tag
        line[5:] = new_elements[idx][0][-3:]
        line = line.astype(str, copy=False)
        strline = line.tolist()
        import ipdb
        #ipdb.set_trace()
        f.write(" ".join(strline) + "\n")
    f.write("$EndElements\n")
    f.close()
    #print(elapsed_time)
    return out_file

in_file = "/Users/erik/Dropbox/Analysis/TMSEEG/TMS007/TMS/SESS_CHAR_TMS007.mat"
electrode_name_file = "/Users/erik/Dropbox/Analysis/TMSEEG/ElectrodeNames.txt"

original_T1 = "/Users/erik/Dropbox/Analysis/TMSEEG/TMS007/20121204_134045t1mpragetms20121113s002a001.nii"
final_volume = "/Users/erik/Dropbox/Analysis/TMSEEG/TMS007/orig_out_flirt.nii"
mesh_file = "/Users/erik/Dropbox/Analysis/TMSEEG/ForwardProblem/TMS007_gmsh.msh"
mesh_id = 2005  # Scalp

get_scalp_tris(in_file, original_T1, final_volume, mesh_file, mesh_id, electrode_name_file)
#get_elements_near_point(
#    [58.00, -51.40, 82.70], "/Users/erik/Dropbox/Analysis/TMSEEG/ForwardProblem/skin.msh", 0, "scalp_tri")

#get_elements_near_point(
#    [3.90, -116.60,  5.40], "/Users/erik/Dropbox/Analysis/TMSEEG/ForwardProblem/skin.msh", 0, "scalp_tri")