import nipype.interfaces.utility as util     # utility
import nipype.pipeline.engine as pe          # pypeline engine

def dipoles_to_geo(in_file, row=0, dipole_name="Dipole1"):
    import os.path as op
    import numpy as np
    from nipype.utils.filemanip import split_filename
    dipole_data = np.loadtxt(in_file)
    _, name, _ = split_filename(in_file)
    out_file = op.abspath(name + "_" + dipole_name + ".geo")
    f = open(out_file, 'w')
    print('Writing dipole .geo file to {f}'.format(f=out_file))
    f.write('View "Dipole locations" {\n')

    if len(dipole_data.shape) == 1:
        dipole_data = [dipole_data]

    dipole_data = [dipole_data[int(row)]]

    for data in dipole_data:
        x = float(data[2])
        y = float(data[3])
        z = float(data[4])
        q_x = float(data[6])
        q_y = float(data[7])
        q_z = float(data[8])
        pt_str = (
            'VP(%5.2f, %5.2f, %5.2f){%5.2f, %5.2f, %5.2f};' % (x, y, z, q_x, q_y, q_z))
        f.write(pt_str + '\n')

    f.write('};\n')
    f.close()
    return out_file

def cost_function_mapping(potential, leadfield, mesh_file, mesh_id):
    import os.path as op
    import numpy as np
    import h5py
    import shutil
    from forward.mesh import read_mesh
    from nipype.utils.filemanip import split_filename
    from nipype import logging
    iflogger = logging.getLogger('interface')

    data_name = "leadfield"
    lf_data_file = h5py.File(leadfield, "r")
    lf_data = lf_data_file.get(data_name)
    leadfield_matrix = lf_data.value
    _, lf_name, _ = split_filename(leadfield)
    V = np.load(potential)


    # ---- Resave the mesh with the cost function as element data ---- #
    mesh_data, _, _, _ = read_mesh(mesh_file, [mesh_id])

    # Create the output mesh file
    path, name, ext = split_filename(mesh_file)
    cost_mesh_file = op.abspath(name + "_cost.msh")

    iflogger.info('Copying current mesh file to %s' % cost_mesh_file)
    shutil.copyfile(mesh_file, cost_mesh_file)

    f = open(cost_mesh_file,'a') #Append to the end of the file
    iflogger.info('Appending cost function scalars to %s' % cost_mesh_file)

    # Write the tag information to the file:
    f.write('$ElementData\n')
    str_tag = '"Cost function %s"' % (leadfield)
    timestep = 0.0001

    f.write('1\n') #Num String tags
    f.write(str_tag + '\n')
    f.write('1\n') #Num Real tags
    f.write('%f\n' % timestep)

    #Three integer tags: timestep, num field components, num elements
    f.write('3\n') #Three int tags
    f.write('0\n') #Time step index
    f.write('1\n') #Num field components
    elem_list = []
    elem_idx_list = []
    cost_list = []

    assert(len(mesh_data)*3 == leadfield_matrix.shape[0])

    for lf_idx, poly in enumerate(mesh_data):
        L = np.transpose(leadfield_matrix[lf_idx * 3:lf_idx * 3 + 3])        
        I = np.eye(len(V))
        Linv = np.linalg.pinv(L)
        val = np.dot(V.T, I - np.dot(L,Linv))
        R = np.dot(val,V)
        cost = np.linalg.norm(R)
        #cost = np.sqrt(R / np.linalg.norm(V))
        cost_list.append(cost)

        elem_idx_list.append(poly["element_id"])

        elementdata_str = ('%d %e\n' % (poly["element_id"], cost))
        elem_list.append(elementdata_str)
        #iflogger.info("%3.3f%%" % (float(lf_idx)/num_polygons*100.0))

    import ipdb
    ipdb.set_trace()
    minidx = cost_list.index(np.min(cost_list))
    maxidx = cost_list.index(np.max(cost_list))
    print(elem_idx_list[minidx], elem_idx_list[maxidx])

    f.write('%d\n' % len(mesh_data)) #Num nonzero field components
    for elementdata_str in elem_list:
        f.write(elementdata_str)

    f.write('$EndElementData\n')
    f.close()

    iflogger.info("Finished writing to %s" % cost_mesh_file)
    cost_geo_file = ""
    return cost_mesh_file, cost_geo_file


def calc_surface_potential_fn(dipole_file, dipole_row, leadfield, mesh_file, mesh_id):
    import os.path as op
    import numpy as np
    import h5py
    from forward.mesh import get_closest_element_to_point

    from nipype import logging
    iflogger = logging.getLogger('interface')

    dipole_data = np.loadtxt(dipole_file)
    if len(dipole_data.shape) == 1:
        dipole_data = [dipole_data]
    dip = dipole_data[int(dipole_row)]
    
    x, y, z = [float(i) for i in dip[2:5]]
    q_x, q_y, q_z = [float(j) for j in dip[6:9]]
    _, element_idx, centroid, element_data, lf_idx = get_closest_element_to_point(mesh_file, mesh_id, [[x, y, z]])

    lf_data_file = h5py.File(leadfield, "r")
    lf_data = lf_data_file.get("leadfield")
    leadfield_matrix = lf_data.value

    L = np.transpose(leadfield_matrix[lf_idx * 3:lf_idx * 3 + 3])
    J = np.array([q_x,q_y,q_z])
    potential = np.dot(L,J)
    out_potential = op.abspath("potential.npy")
    np.save(out_potential,potential)
    return out_potential


def create_cost_function_workflow(name):

    #Inputs - dipole (x,y,z)
    # - leadfield
    # - mesh file
    # - conductivity tensor included?

    workflow = pe.Workflow(name=name)
    inputnode = pe.Node(
        interface=util.IdentityInterface(fields=["dipole_row", "dipole_file", "mesh_file", "mesh_id", "leadfield"]), name="inputnode")

    outputnode = pe.Node(
        interface=util.IdentityInterface(fields=["mesh_file", "dipole_geo_file", "potential_from_dipole"]), name="outputnode")

    dipoles_to_geo_interface = util.Function(
        input_names=["in_file", "row", "dipole_name"],
        output_names=["out_file"], function=dipoles_to_geo)

    write_dipole_geo_file = pe.Node(interface=dipoles_to_geo_interface,
        name="write_dipole_geo_file")

    #Get potential on sensors from "measured" data (from dipole forward simulation)
    calc_surface_potential_interface = util.Function(
        input_names=["dipole_file", "dipole_row", "leadfield", "mesh_file", "mesh_id"],
        output_names=["out_potential"], function=calc_surface_potential_fn)

    calc_surface_potential = pe.Node(interface=calc_surface_potential_interface,
        name="calc_surface_potential")


    #For every element in gray matter, compute two-norm misfit between potential created by
    #element and the "measured potential", save as elementdata in mesh (or scalar point in .geo)

    cost_function_interface = util.Function(
        input_names=["potential", "leadfield", "mesh_file", "mesh_id"],
        output_names=["mesh_file", "geo_file"], function=cost_function_mapping)

    map_cost_function = pe.Node(interface=cost_function_interface,
        name="map_cost_function")

    workflow.connect(
        [(inputnode, calc_surface_potential, [("dipole_file", "dipole_file")])])
    workflow.connect(
        [(inputnode, calc_surface_potential, [("dipole_row", "dipole_row")])])
    workflow.connect(
        [(inputnode, calc_surface_potential, [("leadfield", "leadfield")])])
    workflow.connect(
        [(inputnode, calc_surface_potential, [("mesh_file", "mesh_file")])])
    workflow.connect(
        [(inputnode, calc_surface_potential, [("mesh_id", "mesh_id")])])
    workflow.connect(
        [(calc_surface_potential, map_cost_function, [("out_potential", "potential")])])


    workflow.connect(
        [(inputnode, map_cost_function, [("leadfield", "leadfield")])])
    workflow.connect(
        [(inputnode, map_cost_function, [("mesh_id", "mesh_id")])])
    workflow.connect(
        [(inputnode, map_cost_function, [("mesh_file", "mesh_file")])])

    workflow.connect(
        [(inputnode, write_dipole_geo_file, [("dipole_file", "in_file")])])
    workflow.connect(
        [(inputnode, write_dipole_geo_file, [("dipole_row", "row")])])

    workflow.connect(
        [(write_dipole_geo_file, outputnode, [("out_file", "dipole_geo_file")])])
    workflow.connect(
        [(calc_surface_potential, outputnode, [("out_potential", "potential_from_dipole")])])
    workflow.connect(
        [(map_cost_function, outputnode, [("mesh_file", "mesh_file")])])
    return workflow