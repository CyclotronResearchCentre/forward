import numpy as np
import nibabel as nb
import ipdb
import shutil

def get_num_nodes_from_elm_type(elm_type):
    mapping = [[1, 2], #2-node line.
               [2, 3], #3-node triangle.
               [3, 4], #4-node quadrangle.
               [4, 4], #4-node tetrahedron.
               [5, 8], #8-node hexahedron.
               [6, 6], #6-node prism.
               [7, 5], #5-node pyramid.
               [8, 3], #3-node second order line (2 nodes associated with the vertices and 1 with the edge).
               [9, 6], #6-node second order triangle (3 nodes associated with the vertices and 3 with the edges).
               [10, 9], #9-node second order quadrangle (4 nodes associated with the vertices, 4 with the edges and 1 with the face).
               [11, 10], #10-node second order tetrahedron (4 nodes associated with the vertices and 6 with the edges).
               [12, 27], #27-node second order hexahedron (8 nodes associated with the vertices, 12 with the edges, 6 with the faces and 1 with the volume).
               [13, 18], #18-node second order prism (6 nodes associated with the vertices, 9 with the edges and 3 with the quadrangular faces).
               [14, 14], #14-node second order pyramid (5 nodes associated with the vertices, 8 with the edges and 1 with the quadrangular face).
               [15, 1], #1-node point.
               [16, 8], #8-node second order quadrangle (4 nodes associated with the vertices and 4 with the edges).
               [17, 20], #20-node second order hexahedron (8 nodes associated with the ver- tices and 12 with the edges).
               [18, 15], #15-node second order prism (6 nodes associated with the vertices and 9 with the edges).
               [19, 13], #13-node second order pyramid (5 nodes associated with the vertices and 8 with the edges).
               [20, 9], #9-node third order incomplete triangle (3 nodes associated with the vertices, 6 with the edges)
               [21, 10], #10-node third order triangle (3 nodes associated with the vertices, 6 with the edges, 1 with the face)
               [22, 12], #12-node fourth order incomplete triangle (3 nodes associated with the vertices, 9 with the edges)
               [23, 15], #15-node fourth order triangle (3 nodes associated with the vertices, 9 with the edges, 3 with the face)
               [24, 15], #15-node fifth order incomplete triangle (3 nodes associated with the vertices, 12 with the edges)
               [25, 21], #21-node fifth order complete triangle (3 nodes associated with the vertices, 12 with the edges, 6 with the face)
               [26, 4], #4-node third order edge (2 nodes associated with the vertices, 2 internal to the edge)
               [27, 5], #5-node fourth order edge (2 nodes associated with the vertices, 3 internal to the edge)
               [28, 6], #6-node fifth order edge (2 nodes associated with the vertices, 4 internal to the edge)
               [29, 20], #20-node third order tetrahedron (4 nodes associated with the vertices, 12 with the edges, 4 with the faces)
               [30, 35], #35-node fourth order tetrahedron (4 nodes associated with the vertices, 18 with the edges, 12 with the faces, 1 in the volume)
               [31, 56] #56-node fifth order tetrahedron (4 nodes associated with the vertices, 24 with the edges, 24 with the faces, 4 in the volume)
               ]

    number_of_nodes = mapping[int(elm_type)-1][1]
    return number_of_nodes


def read_mesh(mesh_filename, elements_to_consider):
    print("Reading mesh file: %s" % mesh_filename)
    print("Mesh ids to consider: %s" % elements_to_consider)

    mesh_file = open(mesh_filename, 'r')
    while True:
        line = mesh_file.readline()
        if '$Nodes' in line:
            line = mesh_file.readline()
            number_of_nodes = int(line)
            print("%d nodes in mesh" % number_of_nodes)
            vertex_list = []

            for i in xrange(0, number_of_nodes):
                line = mesh_file.readline()
                node_data = line.split()
                vertex_list.append(node_data)

            vertex_list = np.array((vertex_list)).astype(float)
            print("Done reading nodes")

        elif '$Elements' in line:
            line = mesh_file.readline()
            number_of_elements = int(line)
            print("%d elements in mesh" % number_of_elements)
            polygons = []
            for i in xrange(0, number_of_elements):
                # -- If all elements were quads, each line has 10 numbers.
                #    The first one is a tag, and the last 4 are node numbers.
                line = mesh_file.readline()
                elem_data = line.split()
                if int(elem_data[3]) in elements_to_consider:
                    polygons.append(elem_data)

            polygons = np.array((polygons))
            poly_indices = polygons[:, -3:]
            print("Done reading elements")
            break
    
    # Loop through and assign points to each polygon, save as a dictionary
    mesh_data = []
    num_polygons = len(polygons)
    print("%d polygons found with mesh IDs: %s" % (num_polygons, elements_to_consider))
    for idx, polygon in enumerate(polygons):
        poly_data = {}
        poly_data["element_id"] = int(polygon[0])
        poly_data["number_of_points"] = get_num_nodes_from_elm_type(polygon[1])
        poly_data["node_ids"] = polygon[-poly_data["number_of_points"]:].astype(int)
        poly_data["node_locations"] = vertex_list[poly_data["node_ids"]-1][:,1:]
        poly_data["centroid"] = np.mean(poly_data["node_locations"],0)

        mesh_data.append(poly_data)
        print("%3.3f%%" % (float(idx)/num_polygons*100.0))

    return mesh_data

def include_gmsh_tensor_elements(mesh_file, tensor_file, mask_file, mask_threshold=0.3):
    import numpy as np
    import nibabel as nb
    from nipype.utils.filemanip import split_filename
    import os.path as op

    # Load 4D (6 volume: xx, xy, xz, yy, yz, zz) conductivity tensor image
    tensor_image = nb.load(tensor_file)
    tensor_data = np.flipud(tensor_image.get_data())
    
    # Load mask (usually fractional anisotropy) image
    mask_image = nb.load(mask_file)
    mask_data = np.flipud(mask_image.get_data())
    header = tensor_image.get_header()

    # Make sure the files share the same (3D) dimensions before continuing
    assert(np.shape(tensor_data)[0:3] == np.shape(mask_data)[0:3])

    # Define various constants
    node_id = 0
    element_id = 10000000
    elements_to_consider = [1001] #Use only white matter
    vx, vy, vz = header.get_zooms()[0:3]
    max_x, max_y, max_z = np.shape(tensor_data)[0:3]
    halfx, halfy, halfz = np.array((vx*max_x, vy*max_y, vz*max_z))/2.0
    elements = []
    tensors = {}
    node_info = []
    node_positions = {}

    mesh_data = read_mesh(mesh_file, elements_to_consider)

    # Create the output mesh file
    path, name, ext = split_filename(mesh_file)
    out_file = op.abspath(name + "_cond.msh")
    print('Copying current mesh file to %s' % out_file)
    shutil.copyfile(mesh_file, out_file)

    f = open(out_file,'a') #Append to the end of the file
    print('Appending Conductivity tensors to %s' % out_file)


    # Write the tag information to the file:
    num_polygons = len(mesh_data)
    f.write('$ElementData\n')
    str_tag = '"Conductivity"'
    n_real_tags = str(1)
    timestep = 0.0001

    f.write('1\n') #Num String tags
    f.write(str_tag + '\n')
    f.write('1\n') #Num Real tags
    f.write('%f\n' % timestep)

    #Three integer tags: timestep, num field components, num elements
    f.write('3\n') #Three int tags
    f.write('0\n') #Time step index
    f.write('9\n') #Num field components

    # Get the centroid of all white matter elements
    # Find out which voxel they lie inside
    print("Getting tensor for each element")
    #ipdb.set_trace()
    nonzero = 0
    elem_list = []
    for idx, poly in enumerate(mesh_data):
        i = np.round((poly['centroid'][0]+halfx)/vx).astype(int)
        j = np.round((poly['centroid'][1]+halfy)/vy).astype(int)
        k = np.round((poly['centroid'][2]+halfz)/vz).astype(int)
        #print("%d %d %d" % (i,j,k))
        poly["tensor_triform"] = tensor_data[i,j,k]
        if not all(poly["tensor_triform"] == 0):
            elementdata_str = ('%d %e %e %e %e %e %e %e %e %e\n' % (poly["element_id"], 
                poly["tensor_triform"][0], poly["tensor_triform"][1], poly["tensor_triform"][2],
                poly["tensor_triform"][1], poly["tensor_triform"][3], poly["tensor_triform"][4],
                poly["tensor_triform"][3], poly["tensor_triform"][4], poly["tensor_triform"][5]))
            elem_list.append(elementdata_str)
            nonzero += 1
        print("%3.3f%%" % (float(idx)/num_polygons*100.0))
        #print(nonzero)


    f.write('%d\n' % nonzero) #Num nonzero field components
    for elementdata_str in elem_list:
        f.write(elementdata_str)

    f.write('$EndElementData\n')
    f.close()

    print("Finished writing to %s" % out_file)



    ## Nifti from gmsh
    # Create empty nifti volume (e.g. 256 conformed, 1mm isotropic)
    # Find centroid of all elements
    # Find voxels they lie in, assign e_field value appropriately
    # (many volumes, almost 180 x1 y1 z1 x2 y2 z2)
    return out_file

import os.path as op
mesh_file = op.abspath("ForwardProblem/TMS007_running.msh")
tensor_file = op.abspath("cond_tensor.nii")
mask_file = op.abspath("fa_instructspace.nii")
out_mesh_file = include_gmsh_tensor_elements(mesh_file, tensor_file, mask_file)