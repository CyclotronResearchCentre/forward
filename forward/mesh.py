import numpy as np

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