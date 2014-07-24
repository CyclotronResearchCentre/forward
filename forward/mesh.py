import numpy as np
import logging
logger = logging.getLogger('electrode')

def read_mesh(mesh_filename, elements_to_consider=None):
    from nipype import logging
    iflogger = logging.getLogger('interface')

    iflogger.info("Reading mesh file: %s" % mesh_filename)
    iflogger.info("Mesh ids to consider: %s" % elements_to_consider)

    elem_ids_in_file = []

    mesh_file = open(mesh_filename, 'r')
    while True:
        line = mesh_file.readline()
        if '$Nodes' in line:
            line = mesh_file.readline()
            number_of_nodes = int(line)
            iflogger.info("%d nodes in mesh" % number_of_nodes)
            vertex_list = []

            for i in xrange(0, number_of_nodes):
                line = mesh_file.readline()
                node_data = line.split()
                vertex_list.append(node_data[1:])

            vertex_list = np.array((vertex_list)).astype(float)
            iflogger.info("Done reading nodes")

        elif '$Elements' in line:
            line = mesh_file.readline()
            number_of_elements = int(line)
            iflogger.info("%d elements in mesh" % number_of_elements)
            polygons = []
            for i in xrange(0, number_of_elements):
                # -- If all elements were quads, each line has 10 numbers.
                #    The first one is a tag, and the last 4 are node numbers.
                line = mesh_file.readline()
                elem_data = line.split()
                if isinstance(elements_to_consider, int):
                    elements_to_consider = [elements_to_consider]

                if elements_to_consider is None:
                    polygons.append(np.array(elem_data))
                    elem_ids_in_file.append(int(elem_data[3]))

                elif int(elem_data[3]) in elements_to_consider:
                    polygons.append(np.array(elem_data))

            polygons = np.array((polygons))
            iflogger.info("Done reading elements")
            break
    
    # Loop through and assign points to each polygon, save as a dictionary
    mesh_data = []
    num_polygons = len(polygons)
    if elements_to_consider is None:
        iflogger.info("%d polygons found with mesh IDs: %s" % (num_polygons, elem_ids_in_file))
    else:
        iflogger.info("%d polygons found with mesh IDs: %s" % (num_polygons, elements_to_consider))

    elem = []
    labels = []
    for idx, polygon in enumerate(polygons):
        poly_data = {}
        poly_data["element_id"] = int(polygon[0])
        poly_data["phys_id"] = int(polygon[3])
        poly_data["number_of_points"] = get_num_nodes_from_elm_type(polygon[1])
        poly_data["node_ids"] = polygon[-poly_data["number_of_points"]:].astype(int)
        poly_data["node_locations"] = vertex_list[poly_data["node_ids"]-1][:,:]
        poly_data["centroid"] = np.mean(poly_data["node_locations"],0)

        mesh_data.append(poly_data)
        elem.append(poly_data["node_ids"])
        labels.append(poly_data["phys_id"])
        #nodes             n. of nodes*3 field with the position of the nodes
        #elements          n. of elements*8 field with the elements
        #labels            n. of elements vector with the elements labels
        #iflogger.info("%3.3f%%" % (float(idx)/num_polygons*100.0))

    return mesh_data, vertex_list, elem, labels


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


def get_num_nodes_elements(mesh_filename):
    mesh_file = open(mesh_filename, 'r')
    while True:
        line = mesh_file.readline()

        if line == '$Nodes\n':
            line = mesh_file.readline()
            number_of_nodes = int(line)

        elif line == '$Elements\n':
            line = mesh_file.readline()
            number_of_elements = int(line)
            break
    return number_of_nodes, number_of_elements


def mask_from_labels_fn(in_file, label_values):
    import os.path as op
    import numpy as np
    import nibabel as nb
    from nipype.utils.filemanip import split_filename

    _, name, ext = split_filename(in_file)
    in_file = nb.load(in_file)
    in_data = in_file.get_data()
    new_data = np.zeros(in_data.shape)
    for label in label_values:
        new_data[in_data == label] = 1

    new_image = nb.Nifti1Image(data=new_data,
                               affine=in_file.get_affine(), header=in_file.get_header())
    out_file = op.abspath(name + "_mask" + ext)
    nb.save(new_image, out_file)
    return out_file


def check_intersecting_fn(mesh1, mesh2):
    import subprocess
    intersecting = False
    args = ['meshfix']
    args.append(mesh1)
    args.append(mesh2)
    args.extend(["--shells", "2", "--no-clean", "--intersect"])

    try:
        output = subprocess.check_output(args)
        intersecting = True
    except subprocess.CalledProcessError:
        # No intersections
        intersecting = False
    return intersecting


def cut_inner_fn(outer_mesh, inner_mesh):
    import subprocess
    from nipype.utils.filemanip import split_filename
    import os.path as op
    path, name, ext = split_filename(outer_mesh)
    cut_mesh = op.join(path, name + "_CI.off")
    args = ['meshfix']
    args.extend([outer_mesh, inner_mesh])
    args.extend(["-a", "2.0", "--shells", "2", "--cut-inner", "0"])
    args.extend(["-o", cut_mesh])
    try:
        subprocess.call(args)
    except:
        print("Something went wrong")
    return cut_mesh


def cut_outer_fn(inner_mesh, outer_mesh):
    import subprocess
    from nipype.utils.filemanip import split_filename
    import os.path as op
    path, name, ext = split_filename(inner_mesh)
    cut_mesh = op.join(path, name + "_CO.off")
    args = ['meshfix']
    args.extend([inner_mesh, outer_mesh])
    args.extend(["-a", "2.0", "--shells", "2", "--cut-outer", "0"])
    args.extend(["-o", cut_mesh])
    try:
        subprocess.call(args)
    except:
        print("Something went wrong")
    return cut_mesh


def decouple_outout_fn(outer_mesh, inner_mesh):
    import subprocess
    from nipype.utils.filemanip import split_filename
    import os.path as op
    path, name, ext = split_filename(outer_mesh)
    cut_mesh = op.join(path, name + "_DOO.off")
    args = ['meshfix']
    args.extend([outer_mesh, inner_mesh])
    args.extend(["-a", "2.0", "--shells", "2", "--decouple-outout", "0"])
    args.extend(["-o", cut_mesh])
    try:
        subprocess.call(args)
    except:
        print("Something went wrong")
    return cut_mesh


def decouple_inin_fn(inner_mesh, outer_mesh):
    import subprocess
    from nipype.utils.filemanip import split_filename
    import os.path as op
    path, name, ext = split_filename(inner_mesh)
    cut_mesh = op.join(path, name + "_DII.off")
    args = ['meshfix']
    args.extend([inner_mesh, outer_mesh])
    args.extend(["-a", "2.0", "--shells", "2", "--decouple-inin", "0"])
    args.extend(["-o", cut_mesh])
    try:
        subprocess.call(args)
    except:
        print("Something went wrong")
    return cut_mesh

def decouple_outin_fn(inner_mesh, outer_mesh):
    import subprocess
    from subprocess import CalledProcessError
    import os.path as op
    path, name, ext = split_filename(inner_mesh)
    cut_mesh = op.join(path, name + "_DOI.off")
    args = ['meshfix']
    args.extend([inner_mesh, outer_mesh])
    args.extend(["-a", "2.0", "--shells", "2", "--decouple-outin", "0"])
    args.extend(["-o", cut_mesh])
    try:
        subprocess.call(args)
    except:
        print("Something went wrong")
    return cut_mesh


def clean_mesh_fn(mesh_file):
    import subprocess
    from nipype.utils.filemanip import split_filename
    import os.path as op
    path, name, ext = split_filename(mesh_file)
    cleaned1 = op.join(path, name + "_c1.off")
    cleaned2 = op.join(path, name + "_c2.off")
    args1 = ['meshfix']
    args1.append(mesh_file)
    args1.extend(["-a", "2.0", "-u", "1", "-q"])
    args1.extend(["-o", cleaned1])

    args2 = ['meshfix']
    args2.append(cleaned1)
    args2.extend(["-a", "2.0", "-q"])
    args2.extend(["-o", cleaned2])
    try:
        subprocess.call(args1)
        subprocess.call(args2)
    except:
        print("Something went wrong")
    return cleaned2


def decouple_and_cut_inner_fn(outer_mesh, inner_mesh):
    outer_mesh = decouple_outout_fn(outer_mesh, inner_mesh)
    outer_mesh = cut_inner_fn(outer_mesh, inner_mesh)
    return outer_mesh


def decouple_surfaces_fn(outer_mesh, inner_mesh):
    # This function loops until outer_mesh and inner_mesh
    # have no intersections.
    # At each iteration it pushes outer_mesh out, cuts parts of
    # the inner_mesh away, and cleans the mesh.
    # The pushing out is also iterative.
    #
    from forward.mesh import (
        iter_push_out_fn, check_intersecting_fn,
        cut_inner_fn, clean_mesh_fn)

    for i in range(0, 2):
        intersections = check_intersecting_fn(outer_mesh, inner_mesh)
        if intersections:
            outer_mesh = iter_push_out_fn(outer_mesh, inner_mesh, iterations=3)
            outer_mesh = cut_inner_fn(outer_mesh, inner_mesh)
            outer_mesh = clean_mesh_fn(outer_mesh)
        else:
            break
    return outer_mesh


def iter_remove_throats_fn(outer_mesh, inner_mesh, iterations=5):
    # This function loops until there are no
    # intersections between two surfaces.
    # It cuts away the inner surface and
    # uses decouple_outout:
    #
    # "Treat 1st file as outer, 2nd file as inner component.
    ## "Resolve overlaps by moving outers triangles outwards."
    ## "Constrain the min distance between the components > d."
    #
    # and cut-inner:
    #
    ## "Remove triangles of 1st that are inside  of the 2nd shell."
    ## "Dilate 2nd by d; Fill holes and keep only 1st afterwards."
    #
    # at each iteration.
    from forward.mesh import (check_intersecting_fn,
                                           decouple_and_cut_inner_fn)
    for i in range(0, iterations):
        intersections = check_intersecting_fn(outer_mesh, inner_mesh)
        if intersections:
            outer_mesh = decouple_and_cut_inner_fn(outer_mesh, inner_mesh)
        else:
            break
    return outer_mesh


def iter_push_out_fn(outer_mesh, inner_mesh, iterations=5):
    # This function loops until there are no
    # intersections between two surfaces.
    # It cuts away the inner surface and
    # uses decouple_outout:
    #
    # "Treat 1st file as outer, 2nd file as inner component.
    ## "Resolve overlaps by moving outers triangles outwards."
    ## "Constrain the min distance between the components > d."
    #
    # at each iteration.
    from forward.mesh import (check_intersecting_fn,
                                           decouple_outout_fn)
    for i in range(0, iterations):
        intersecting = check_intersecting_fn(outer_mesh, inner_mesh)
        if intersecting:
            outer_mesh = decouple_outout_fn(outer_mesh, inner_mesh)
        else:
            break
    return outer_mesh


def iter_decoupling_fn(outer_mesh, inner_mesh):
    from forward.mesh import (check_intersecting_fn,
                                           cut_inner_fn, decouple_outout_fn, clean_mesh_fn)

    intersections = check_intersecting_fn(outer_mesh, inner_mesh)
    while(intersections):
        for i in range(0, 3):
            intersections = check_intersecting_fn(outer_mesh, inner_mesh)
            if intersections:
                outer_mesh = decouple_outout_fn(outer_mesh, inner_mesh)
            else:
                break
            outer_mesh = cut_inner_fn(outer_mesh, inner_mesh)
            outer_mesh = clean_mesh_fn(outer_mesh)
    return outer_mesh


def remove_spikes_fn(outer_mesh, inner_mesh):
    from forward.mesh import (check_intersecting_fn,
                                           cut_outer_fn, decouple_inin_fn, clean_mesh_fn)

    intersections = check_intersecting_fn(outer_mesh, inner_mesh)
    if intersections:
        inner_mesh = cut_outer_fn(inner_mesh, outer_mesh)
        inner_mesh = clean_mesh_fn(inner_mesh)

        for i in xrange(0, 3):
            intersections = check_intersecting_fn(outer_mesh, inner_mesh)
            if intersections:
                inner_mesh = decouple_inin_fn(inner_mesh, outer_mesh)
                inner_mesh = clean_mesh_fn(inner_mesh)
            else:
                break
    return inner_mesh


def decouple_ventricles_fn(ventricles, white_matter):
    from forward.mesh import (check_intersecting_fn,
                                           cut_outer_fn, decouple_inin_fn, clean_mesh_fn)

    intersections = check_intersecting_fn(white_matter, ventricles)
    while(intersections):
        ventricles = decouple_inin_fn(ventricles, white_matter)
        ventricles = cut_outer_fn(ventricles, white_matter)
        ventricles = clean_mesh_fn(ventricles)
        intersections = check_intersecting_fn(white_matter, ventricles)
    return ventricles

def decouple_input_from_GM_fn(mesh_file, gray_matter):
    from forward.mesh import (check_intersecting_fn,
                                           cut_outer_fn, decouple_outin_fn, clean_mesh_fn)

    intersections = check_intersecting_fn(gray_matter, mesh_file)
    while(intersections):
        mesh_file = decouple_outin_fn(mesh_file, gray_matter)
        mesh_file = cut_outer_fn(mesh_file, gray_matter)
        mesh_file = clean_mesh_fn(mesh_file)
        intersections = check_intersecting_fn(gray_matter, mesh_file)
    return mesh_file


def decouple_outout_cutin_fn(outer_mesh, inner_mesh):
    from forward.mesh import (check_intersecting_fn,
                                           cut_inner_fn, decouple_outout_fn, clean_mesh_fn)

    intersections = check_intersecting_fn(outer_mesh, inner_mesh)
    while(intersections):
        outer_mesh = decouple_outout_fn(outer_mesh, inner_mesh)
        outer_mesh = cut_inner_fn(outer_mesh, inner_mesh)
        outer_mesh = clean_mesh_fn(outer_mesh)
        intersections = check_intersecting_fn(outer_mesh, inner_mesh)
    return outer_mesh


def get_closest_element_index(point, mesh_filename, mesh_id):
    import numpy as np
    from scipy.spatial.distance import cdist
    import os.path as op
    import time
    start_time = time.time()
    mesh_file = open(mesh_filename, 'r')
    closest_element_idx = 0
    n_brain_elem = 0
    nodes_of_potential_elements = []
    potential_elem_ids = []
    while True:
        line = mesh_file.readline()

        if line == '$Nodes\n':
            line = mesh_file.readline()
            number_of_nodes = int(line)
            vertices = []

            for i in xrange(0, number_of_nodes):
                line = mesh_file.readline()
                node_data = line.split()
                vert = node_data
                vertices.append(vert)

            vertices = np.array((vertices))

            dist = cdist(point, vertices[:, 1:], 'euclidean')

            min_dist_by_elec = dist.min(axis=1)

            closest_nodes = []
            for idx, min_dist in enumerate(min_dist_by_elec):
                closest = np.where(dist[idx] == min_dist)[0][0]
                closest_idx = vertices[closest][0]
                closest_nodes.append(closest_idx)

            closest_nodes = np.array(closest_nodes)
            logger.info("Closest Node IDs")
            logger.info(closest_nodes)

        elif line == '$Elements\n':
            line = mesh_file.readline()

            number_of_elements = int(line)
            for i in xrange(0, number_of_elements):
                line = mesh_file.readline()
                elem_data = line.split()
                if int(elem_data[3]) == mesh_id:
                    n_brain_elem = n_brain_elem + 1
                    nodes = np.array(elem_data[5:])
                    mask = np.in1d(nodes, closest_nodes)
                    if np.sum(mask) >= 1:
                        # Needs to be reworked to give the single closest element!
                        nodes_of_potential_elements.append(nodes)
                        potential_elem_ids.append(n_brain_elem)
                        closest_element_idx = n_brain_elem

        elif line == '$EndElementData\n' or len(line) == 0:
            break

    #get_closest_element(vertices, nodes_of_potential_elements)
    mesh_file.close()
    elapsed_time = time.time() - start_time
    print(elapsed_time)
    return closest_element_idx, number_of_elements, n_brain_elem

def get_closest_element_to_point(mesh_filename, elements_to_consider, point):
    from nipype import logging
    from scipy.spatial.distance import cdist
    iflogger = logging.getLogger('interface')

    iflogger.info("Reading mesh file: %s" % mesh_filename)
    iflogger.info("Mesh ids to consider: %s" % elements_to_consider)

    mesh_file = open(mesh_filename, 'r')
    while True:
        line = mesh_file.readline()
        if '$Nodes' in line:
            line = mesh_file.readline()
            number_of_nodes = int(line)
            iflogger.info("%d nodes in mesh" % number_of_nodes)
            vertex_list = []

            for i in xrange(0, number_of_nodes):
                line = mesh_file.readline()
                node_data = line.split()
                vertex_list.append(node_data)

            vertex_list = np.array((vertex_list)).astype(float)
            iflogger.info("Done reading nodes")

        elif '$Elements' in line:
            line = mesh_file.readline()
            number_of_elements = int(line)
            iflogger.info("%d elements in mesh" % number_of_elements)
            polygons = []
            for i in xrange(0, number_of_elements):
                # -- If all elements were quads, each line has 10 numbers.
                #    The first one is a tag, and the last 4 are node numbers.
                line = mesh_file.readline()
                elem_data = line.split()
                if isinstance(elements_to_consider, int):
                    elements_to_consider = [elements_to_consider]
                if int(elem_data[3]) in elements_to_consider:
                    polygons.append(np.array(elem_data))

            polygons = np.array((polygons))
            iflogger.info("Done reading elements")
            break
    
    # Loop through and assign points to each polygon, save as a dictionary
    mesh_data = []
    num_polygons = len(polygons)
    iflogger.info("%d polygons found with mesh IDs: %s" % (num_polygons, elements_to_consider))
    min_dist = 200
    for idx, polygon in enumerate(polygons):
        poly_data = {}
        poly_data["element_id"] = int(polygon[0])
        poly_data["number_of_points"] = get_num_nodes_from_elm_type(polygon[1])
        poly_data["node_ids"] = polygon[-poly_data["number_of_points"]:].astype(int)
        poly_data["node_locations"] = vertex_list[poly_data["node_ids"]-1][:,1:]
        poly_data["centroid"] = np.mean(poly_data["node_locations"],0)

        poly_data["distance_to_pt"] = dist = cdist(point, np.array([poly_data["centroid"]]), 'euclidean')
        if dist <= min_dist:
            min_dist = dist
            closest_element_idx = poly_data["element_id"]
            lf_idx = idx
            closest_element_data = poly_data
            centroid = poly_data["centroid"]


        mesh_data.append(poly_data)
        #iflogger.info("%3.3f%%" % (float(idx)/num_polygons*100.0))

    return mesh_data, closest_element_idx, centroid, closest_element_data, lf_idx