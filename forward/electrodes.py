import logging
logger = logging.getLogger('electrode')

def rewrite_mesh_with_electrodes(electrode_position_file, electrode_name_file, mesh_file, mesh_id):
    import os.path as op
    import numpy as np
    from nipype.utils.filemanip import split_filename
    _, name, _ = split_filename(mesh_file)
    out_basename = op.abspath(name + "_elec")

    try:
        positions = np.loadtxt(electrode_position_file, dtype=float)
    except ValueError:
        positions = np.loadtxt(electrode_position_file, dtype=float, delimiter=",")

    electrode_names = np.loadtxt(electrode_name_file, dtype=str)
    num_electrodes = len(positions)
    new_phys_ids = range(5000, 5000 + num_electrodes)
    new_elem_ids = range(6000, 6000 + num_electrodes)

    out_file = get_closest_element_to_points(
        positions, mesh_file, mesh_id, new_phys_ids, new_elem_ids, out_basename)
    print('Writing electrodes to %s' % out_file)
    write_electrode_labels(positions, electrode_names)
    return out_file


def CRC_get_electrode_location(in_file, electrode_name_file, session_number=1):
    '''
    This function is only useful if you're using the TMS/EEG system at the CRC
    and are dealing with the associated MATLAB files
    '''
    import scipy.io as sio
    import os.path as op
    import numpy as np
    struct = sio.loadmat(in_file, squeeze_me=True)
    electrodes = struct['ELECTRODESF']
    position = electrodes[session_number - 1]
    points = []
    for pt_id, data in enumerate(position):
        x = -float(data[0]) + 112
        y = float(data[2]) - 127.5
        z = float(data[1]) - 127.5
        point = [x, y, z]
        points.append(point)

    points = np.array(points)
    out_file = op.abspath("ElectrodePositions_sess" + str(session_number) + ".txt")
    np.savetxt(out_file, points)
    return out_file

def write_electrode_labels(points, electrode_names, out_file="ElectrodeLabels.geo"):
    import os.path as op
    f = open(op.abspath(out_file), 'w')
    f.write('View "Electrode Names" {\n')
    for pt_id, pt_data in enumerate(points):
        pt_str = ('T3(%5.2f, %5.2f, %5.2f,0){"%s"};' % (
            pt_data[0], pt_data[1], pt_data[2], electrode_names[pt_id]))
        f.write(pt_str + '\n')
    f.write('};\n')
    f.close()
    print('Writing electrode labels to %s' % out_file)
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
            logger.info("Closest Node IDs")
            logger.info(closest_nodes)

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
                    mask = np.in1d(nodes, closest_nodes)
                    if np.sum(mask) >= 1:
                        which_elec = np.where(
                            closest_nodes == nodes[mask][0])[0][0]
                        elem_data[3] = str(new_phys_ids[which_elec])
                        elem_data[4] = str(new_elem_ids[which_elec])

                    rw_line = " ".join(elem_data) + "\n"
                    f.write(rw_line)
                else:
                    f.write(line)
        elif line == '$EndElementData\n' or len(line) == 0:
            break

    mesh_file.close()
    f.close()
    elapsed_time = time.time() - start_time
    print(elapsed_time)
    return out_file

def get_closest_element_to_points(points, mesh_filename, elements_to_consider, new_phys_ids,
                             new_elem_ids, out_basename):
    import numpy as np
    from scipy.spatial.distance import cdist
    import os.path as op
    import time
    from nipype import logging
    from forward.mesh import get_num_nodes_from_elm_type

    iflogger = logging.getLogger('interface')

    iflogger.info("Reading mesh file: %s" % mesh_filename)
    iflogger.info("Mesh ids to consider: %s" % elements_to_consider)

    start_time = time.time()
    mesh_file = open(mesh_filename, 'r')
    out_file = op.abspath(out_basename + ".msh")

    f = open(out_file, 'w')
    while True:
        line = mesh_file.readline()
        f.write(line)

        if line == '$Nodes\n':
            line = mesh_file.readline()
            f.write(line)
            number_of_nodes = int(line)
            iflogger.info("%d nodes in mesh" % number_of_nodes)
            vertex_list = []

            for i in xrange(0, number_of_nodes):
                line = mesh_file.readline()
                f.write(line)
                node_data = line.split()
                vertex_list.append(node_data)

            vertex_list = np.array((vertex_list)).astype(float)
            iflogger.info("Done reading nodes")

        elif line == '$Elements\n':
            line = mesh_file.readline()
            f.write(line)
            number_of_elements = int(line)
            iflogger.info("%d elements in mesh" % number_of_elements)
            elem_lines = []
            centroids = []
            ids = []
            for i in xrange(0, number_of_elements):
                # -- If all elements were quads, each line has 10 numbers.
                #    The first one is a tag, and the last 4 are node numbers.
                line = mesh_file.readline()
                elem_lines.append(line)
                elem_data = line.split()
                if isinstance(elements_to_consider, int):
                    elements_to_consider = [elements_to_consider]
                if int(elem_data[3]) in elements_to_consider:
                    #polygons.append(np.array(elem_data))
                    polygon = np.array(elem_data)
                    poly_data = {}

                    poly_data["element_id"] = int(polygon[0])
                    poly_data["number_of_points"] = get_num_nodes_from_elm_type(polygon[1])
                    poly_data["node_ids"] = polygon[-poly_data["number_of_points"]:].astype(int)
                    poly_data["node_locations"] = vertex_list[poly_data["node_ids"]-1][:,1:]
                    poly_data["centroid"] = np.mean(poly_data["node_locations"],0)
                    centroids.append(poly_data["centroid"])
                    ids.append(poly_data["element_id"])

            centroids = np.array(centroids)
            dist = cdist(points, centroids, 'euclidean')
            min_by_elec = np.min(dist,1)
            closest_elem_ids = []
            for idx in xrange(0,dist.shape[0]):
                min_idx = np.where(dist[idx]==min_by_elec[idx])
                closest_elem_ids.append(ids[min_idx[0][0]])

            # Rewrite the file with the changed element
            for line in elem_lines:
                elem_data = line.split()
                if int(elem_data[0]) in closest_elem_ids:
                    which_elec = closest_elem_ids.index(int(elem_data[0]))
                    elem_data[3] = str(new_phys_ids[which_elec])
                    elem_data[4] = str(new_elem_ids[which_elec])
                    line = " ".join(elem_data) + "\n"

                f.write(line)

        elif line == '$EndElementData\n' or len(line) == 0:
            break

    mesh_file.close()
    f.close()
    elapsed_time = time.time() - start_time
    print(elapsed_time)
    return out_file
