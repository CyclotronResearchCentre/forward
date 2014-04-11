def read_mesh_elem_data(mesh_filename, view_name="Cost function"):
    import numpy as np

    from nipype import logging
    iflogger = logging.getLogger('interface')

    iflogger.info("Reading mesh file: %s" % mesh_filename)

    mesh_file = open(mesh_filename, 'r')
    while True:
        line = mesh_file.readline()
        if '$ElementData' in line:
            line = mesh_file.readline()
            name = mesh_file.readline()
            if name == '"' + view_name + '"\n':
                line = mesh_file.readline()
                line = mesh_file.readline()
                line = mesh_file.readline()
                line = mesh_file.readline()
                line = mesh_file.readline()
                number_of_elements = int(mesh_file.readline().replace("\n",""))

                elem_data = []
                iflogger.info("%d elements in element data" % number_of_elements)
                for i in xrange(0, number_of_elements):
                    line = mesh_file.readline()
                    line_data = line.split()
                    polygon = {}
                    polygon["element_id"] = int(line_data[0])
                    polygon["data"] = float(line_data[1])
                    elem_data.append(polygon)

                iflogger.info("Done reading element data")
                break
    
    # Loop through and assign points to each polygon, save as a dictionary
    num_polygons = len(elem_data)
    iflogger.info("%d polygons found" % num_polygons)

    return elem_data

def subtract_meshes(mesh1, mesh2, view_name1, view_name2, out_file="Mesh1-Mesh2.msh"):
    '''
    Compares two mesh files by their element data

    e.g. cost functions calculated from leadfields with isotropic
    white matter conductivity and anisotropic conductivity
    '''
    import os.path as op

    from nipype import logging    
    iflogger = logging.getLogger('interface')

    print("Reading mesh 1: %s" % mesh1)
    mesh_data1 = read_mesh_elem_data(mesh1, view_name1)

    print("Reading mesh 2: %s" % mesh2)
    mesh_data2 = read_mesh_elem_data(mesh2, view_name2)

    # Create the output mesh file
    out_file = op.abspath(out_file)

    f = open(out_file,'a') #Append to the end of the file
    iflogger.info('Writing difference file to %s' % out_file)
    import ipdb
    ipdb.set_trace()
    elem_data_diff_list = []
    for idx, val in enumerate(mesh_data1):
        d1 = mesh_data1[idx]["data"]
        d2 = mesh_data2[idx]["data"]
        diff = d1 - d2

        diff_str = ('%d %e \n' % (mesh_data1[idx]["element_id"], diff))
        elem_data_diff_list.append(diff_str)


    # Write the first part of the mesh file
    mesh_file = open(mesh1, 'r')
    while True:
        line = mesh_file.readline()
        f.write(line)
        if '$ElementData' in line:
            break

    str_tag = '"Difference %s - %s"' % (mesh1,mesh2)
    timestep = 0.0001
    f.write('1\n') #Num String tags
    f.write(str_tag + '\n')
    f.write('1\n') #Num Real tags
    f.write('%f\n' % timestep)
    #Three integer tags: timestep, num field components, num elements
    f.write('3\n') #Three int tags
    f.write('0\n') #Time step index
    f.write('1\n') #Num field components
    f.write('%d\n' % len(mesh_data1)) #Num nonzero field components
    for elementdata_str in elem_data_diff_list:
        f.write(elementdata_str)
    f.write('$EndElementData\n')
    f.close()

    iflogger.info("Finished writing to %s" % out_file)
    return out_file