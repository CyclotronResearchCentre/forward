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
    points = []
    new_phys_ids = []
    new_elem_ids = []
    for pt_id, data in enumerate(position):
        x = -float(data[0]) + 112
        y = float(data[2]) - 127.5
        z = float(data[1]) - 127.5
        point = [x, y, z]
        points.append(point)
        new_phys_ids.append(pt_id+5000)
        new_elem_ids.append(pt_id+6000)

    points = np.array(points)
    out_basename = op.abspath(basename + str(session_number) + "_" + "elec_tris")
    out_file = get_elements_near_points(
        points, mesh_file, mesh_id, new_phys_ids, new_elem_ids, out_basename)
    print('Writing electrode to {f}'.format(ec=electrode_names[pt_id], f=out_file))
    return out_file


def get_elements_near_points(points, mesh_filename, mesh_id, new_phys_ids, 
        new_elem_ids, out_basename):
    import numpy as np
    from scipy.spatial.distance import cdist
    import os.path as op
    import time
    start_time = time.time()
    mesh_file = open(mesh_filename, 'r')
    out_file = op.abspath(out_basename + ".msh")

    f = open(out_file, 'w')
    while True:
        line = mesh_file.readline()
        f.write(line)

        if line == '$Nodes\n':
            ipdb.set_trace()
            line = mesh_file.readline()
            number_of_nodes = int(line)
            f.write(line)
            vertices = []

            for i in xrange(0, number_of_nodes):
                line = mesh_file.readline()
                f.write(line)

                node_data = line.split()
                vert = node_data
                vertices.append(vert)

            vertices = np.array((vertices))

            dist = cdist(points, vertices[:, 1:], 'euclidean')

            min_dist_by_elec = dist.min(axis=1)

            closest_nodes = []
            for idx, min_dist in enumerate(min_dist_by_elec):
                closest = np.where(dist[idx] == min_dist)[0][0]
                closest_idx = vertices[closest][0]
                closest_nodes.append(closest_idx)

            closest_nodes = np.array(closest_nodes)
            print(closest_nodes)

        elif line == '$Elements\n':
            line = mesh_file.readline()
            f.write(line)

            number_of_elements = int(line)
            for i in xrange(0, number_of_elements):
                # -- If all elements were quads, each line has 10 numbers.
                #    The first one is a tag, and the last 4 are node numbers.
                line = mesh_file.readline()
                elem_data = line.split()
                if int(elem_data[3]) == mesh_id:
                    nodes = np.array(elem_data[5:])
                    mask = np.in1d(nodes,closest_nodes)
                    if np.sum(mask) >= 1:
                        which_elec = np.where(closest_nodes == nodes[mask][0])[0]
                        elem_data[3] = str(new_phys_ids[which_elec])
                        #elem_data[4] = str(new_elem_ids[which_elec])

                    rw_line = " ".join(elem_data) + "\n"
                    f.write(rw_line)
                else:
                    f.write(line)
            break

    mesh_file.close()
    f.write("$EndElements\n")
    f.close()

    elapsed_time = time.time() - start_time
    print(elapsed_time)
    return out_file

in_file = "/Users/erik/Dropbox/Analysis/TMSEEG/TMS007/TMS/SESS_CHAR_TMS007.mat"
electrode_name_file = "/Users/erik/Dropbox/Analysis/TMSEEG/ElectrodeNames.txt"

original_T1 = "/Users/erik/Dropbox/Analysis/TMSEEG/TMS007/20121204_134045t1mpragetms20121113s002a001.nii"
final_volume = "/Users/erik/Dropbox/Analysis/TMSEEG/TMS007/orig_out_flirt.nii"
#mesh_file = "/Users/erik/Dropbox/Analysis/TMSEEG/ForwardProblem/skinmsh.msh"
mesh_file = "/Users/erik/Dropbox/Analysis/TMSEEG/ForwardProblem/TMS007_gmsh.msh"
mesh_id = 1005  # Scalp

get_scalp_tris(in_file, original_T1, final_volume, mesh_file, mesh_id, electrode_name_file)