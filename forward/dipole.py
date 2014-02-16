import nipype.interfaces.utility as util     # utility
import nipype.pipeline.engine as pe          # pypeline engine
from .leadfield import get_electrode_indices
from nipype.interfaces.gmsh import GetDP


def write_dipole_pro_file(mesh_file, conductivity_tensor_included, ground_index, electrode_indices, electrode_name_file, orig_pro_file=None, dipole_name="Dipole1"):
    import os
    import os.path as op
    import numpy as np

    if orig_pro_file is None:
        orig_pro_file = op.join(os.environ["FWD_DIR"], "etc/eeg_dipole.pro")

    problem_file = op.abspath("eeg_forward_dipole_" + str(dipole_name) + ".pro")

    elec_names = np.loadtxt(electrode_name_file, dtype=str)
    electrode_namelist = elec_names.tolist()
    sensor_idlist = [x+5000 for x in electrode_indices]
    
    sensor_regions = {}
    for val in electrode_indices:
        sensor_regions[electrode_namelist[val]] = sensor_idlist[val]

    ground_region = ground_index + 5000

    original = open(orig_pro_file, 'r')
    modified = open(problem_file, 'w')
    while True:
        line = original.readline()
        if len(line) == 0:
            break # EOF
        elif line.rfind('// **This is replaced by dipole.py**') >= 0:
            modified.write("  /* ----- The lines below were replaced by dipole.py ----- */ \n")
            for elec_id in sensor_regions.keys():
                modified.write("  %s = Region[%i];\n" % (elec_id, sensor_regions[elec_id]))
            line = "  /* ----- The lines above were replaced by dipole.py ----- */\n"
        elif line.rfind('Sink = Region[') >= 0:
            line = "  Sink = Region[%d];\n" % ground_region
        elif line.rfind('Electrodes = Region[{Sink, Source}];') >= 0:
            line = "  Electrodes = Region[{%s}];\n" % ",".join(sensor_regions.keys())
        elif line.rfind('sigma[WhiteMatter_Cerebellum]') >= 0 and conductivity_tensor_included:
            line = "    sigma[WhiteMatter_Cerebellum] = TensorField[XYZ[]] #1 ? #1 : 0.33 ;//: 0.33*1e-3 ;\n"
        elif line.rfind('sigma[WhiteMatter_Cerebellum]') >= 0 and not conductivity_tensor_included:
            line = "    sigma[WhiteMatter_Cerebellum] = 0.33;\n"
        modified.write(line)

    modified.close()
    original.close()
    return problem_file


def rewrite_mesh_with_dipole(dipole_file, mesh_file, mesh_id, row=0):
    import os.path as op
    import numpy as np
    from nipype.utils.filemanip import split_filename
    from forward.electrodes import get_closest_element_to_points
    _, name, _ = split_filename(mesh_file)
    out_basename = op.abspath(name + "_dipole")

    dipole_data = np.loadtxt(dipole_file)
    dip = dipole_data[int(row)]
    
    x = float(dip[2])
    y = float(dip[3])
    z = float(dip[4])
    #q_x = float(dip[6])
    #q_y = float(dip[7])
    #q_z = float(dip[8])

    positions = [[x, y, z]]

    new_phys_ids = [8000]
    new_elem_ids = [9000]

    out_file = get_closest_element_to_points(
        positions, mesh_file, mesh_id, new_phys_ids, new_elem_ids, out_basename)
    return out_file


def dipoles_to_geo(in_file, row=None, dipole_name="Dipole1"):
    import os.path as op
    import numpy as np
    from nipype.utils.filemanip import split_filename
    dipole_data = np.loadtxt(in_file)
    _, name, _ = split_filename(in_file)
    out_file = op.abspath(name + "_" + dipole_name + ".geo")
    f = open(out_file, 'w')
    print('Writing dipole .geo file to {f}'.format(f=out_file))
    f.write('View "Dipole locations" {\n')

    if row is not None:
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


def cost_function_mapping(table_files, leadfield, mesh_file, ground_index):
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

    voltage_filename = "v_elec"

    if isinstance(table_files, str):
        table_files = [table_files]

    for in_file in table_files:
        _, name, _ = split_filename(in_file)
        if name.rfind(voltage_filename) >= 0:
            potential = np.loadtxt(in_file)

    if not potential.any():
        raise Exception

    import ipdb
    ipdb.set_trace()

    potential_from_dipole = potential[:,8]

    dipole_data = read_mesh(mesh_file, 8000)

    # Two-norm misfit
    #np.linalg.norm(potential_from_dipole - dipole_data[0]['element_id'])

    # ---- Resave the mesh with the cost function as element data ---- #
    elements_to_consider = [1002]
    mesh_data = read_mesh(mesh_file, elements_to_consider)

    # Create the output mesh file
    path, name, ext = split_filename(mesh_file)
    cost_mesh_file = op.abspath(name + "_cost.msh")

    iflogger.info('Copying current mesh file to %s' % cost_mesh_file)
    shutil.copyfile(mesh_file, cost_mesh_file)

    f = open(cost_mesh_file,'a') #Append to the end of the file
    iflogger.info('Appending root mean squared error scalars to %s' % cost_mesh_file)

    # Write the tag information to the file:
    num_polygons = len(mesh_data)
    f.write('$ElementData\n')
    str_tag = '"Cost function %s %s"' % (leadfield)
    timestep = 0.0001

    f.write('1\n') #Num String tags
    f.write(str_tag + '\n')
    f.write('1\n') #Num Real tags
    f.write('%f\n' % timestep)

    #Three integer tags: timestep, num field components, num elements
    f.write('3\n') #Three int tags
    f.write('0\n') #Time step index
    f.write('1\n') #Num field components

    nonzero = 0
    cost = []
    elem_list = []
    for idx, poly in enumerate(mesh_data):
        elementdata_str = ('%d %e\n' % (poly["element_id"], cost[idx]))
        elem_list.append(elementdata_str)
        nonzero += 1
        iflogger.info("%3.3f%%" % (float(idx)/num_polygons*100.0))

    f.write('%d\n' % nonzero) #Num nonzero field components
    for elementdata_str in elem_list:
        f.write(elementdata_str)

    f.write('$EndElementData\n')
    f.close()

    iflogger.info("Finished writing to %s" % cost_mesh_file)


    cost_geo_file = ""
    return cost_mesh_file, cost_geo_file



def create_cost_function_workflow(name, conductivity_tensor_included=False):

    #Inputs - dipole (x,y,z)
    # - leadfield
    # - mesh file
    # - conductivity tensor included?

    #Write dipole as geo file for display
    #Write dipole into mesh as ID 8000
    #Run dipole model

    '''
    Define the problem file
    '''

    dipole_row = 0

    workflow = pe.Workflow(name=name)
    inputnode = pe.Node(
        interface=util.IdentityInterface(fields=["dipole_file", "mesh_file", "leadfield", "electrode_name_file",
            "marker_list", "ground_electrode", "conductivity_tensor_included"]), name="inputnode")
    outputnode = pe.Node(
        interface=util.IdentityInterface(fields=["mesh_file", "cost_function_geo"]), name="outputnode")

    get_elec_indices_interface = util.Function(
        input_names=["electrode_name_file", "marker_list", "ground_electrode"],
        output_names=["electrode_indices", "ground_index"], function=get_electrode_indices)

    get_elec_indices = pe.Node(
        interface=get_elec_indices_interface, name="get_elec_indices")


    # Adding dipole to the mesh file
    add_dipole_to_mesh_interface = util.Function(
        input_names=["dipole_file", "mesh_file", "mesh_id", "row"],
        output_names=["out_file"], function=rewrite_mesh_with_dipole)

    add_dipole_to_mesh = pe.Node(interface=add_dipole_to_mesh_interface,
        name="add_dipole_to_mesh")
    add_dipole_to_mesh.inputs.mesh_id = 1002
    add_dipole_to_mesh.inputs.row = dipole_row

    dipoles_to_geo_interface = util.Function(
        input_names=["in_file", "row", "dipole_name"],
        output_names=["out_file"], function=dipoles_to_geo)

    write_dipole_geo_file = pe.Node(interface=dipoles_to_geo_interface,
        name="write_dipole_geo_file")

    # Define interface to rewrite the Problem for each electrode
    write_problem_file_interface = util.Function(
        input_names=["mesh_file", "conductivity_tensor_included", "dipole_name", "ground_index", "electrode_indices", "electrode_name_file", "orig_pro_file"],
        output_names=["problem_file"], function=write_dipole_pro_file)

    write_problem_file = pe.Node(interface=write_problem_file_interface,
                                    name="write_problem_file")

    write_problem_file.inputs.dipole_name = "Dipole1_RightHemi" # For the dipole source

    # Define MapNode to run the forward model M times
    run_forward_model = pe.Node(interface=GetDP(), name="run_forward_model")
    run_forward_model.inputs.solve = "Electrostatics"
    run_forward_model.inputs.binary_output_files = True #for debugging
    
    run_forward_model.inputs.out_table_filenames = ["v_elec", "e_brain"]
    run_forward_model.inputs.out_pos_filenames = ["v", "V_sink", "I_sink", "R_sink"]


    #Get potential on sensors from "measured" data (from dipole forward simulation)    
    #For every element in gray matter, compute two-norm misfit between potential created by
    #element and the "measured potential", save as elementdata in mesh (or scalar point in .geo)

    cost_function_interface = util.Function(
        input_names=["table_files", "leadfield", "mesh_file", "ground_index"],
        output_names=["mesh_file", "geo_file"], function=cost_function_mapping)

    map_cost_function = pe.Node(interface=cost_function_interface,
        name="map_cost_function")

    map_cost_function.inputs.ground_index = 5008


    # Start connecting the inputnode and the electrode identification function
    workflow.connect(
        [(inputnode, get_elec_indices, [("electrode_name_file", "electrode_name_file")])])
    workflow.connect(
        [(inputnode, get_elec_indices, [("marker_list", "marker_list")])])
    workflow.connect(
        [(inputnode, get_elec_indices, [("ground_electrode", "ground_electrode")])])

    workflow.connect(
        [(inputnode, add_dipole_to_mesh, [("mesh_file", "mesh_file")])])
    workflow.connect(
        [(inputnode, add_dipole_to_mesh, [("dipole_file", "dipole_file")])])

    workflow.connect(
        [(inputnode, write_dipole_geo_file, [("dipole_file", "in_file")])])

    workflow.connect(
        [(add_dipole_to_mesh, write_problem_file, [("out_file", "mesh_file")])])
    workflow.connect(
        [(get_elec_indices, write_problem_file, [("ground_index", "ground_index")])])
    workflow.connect(
        [(get_elec_indices, write_problem_file, [("electrode_indices", "electrode_indices")])])
    workflow.connect(
        [(inputnode, write_problem_file, [("conductivity_tensor_included", "conductivity_tensor_included")])])
    workflow.connect(
        [(inputnode, write_problem_file, [("electrode_name_file", "electrode_name_file")])])

    workflow.connect(
        [(write_problem_file, run_forward_model, [("problem_file", "problem_file")])])
    workflow.connect(
        [(add_dipole_to_mesh, run_forward_model, [("out_file", "mesh_file")])])
    if conductivity_tensor_included:
        workflow.connect(
            [(add_dipole_to_mesh, run_forward_model, [("out_file", "gmsh_read_file")])])

    workflow.connect(
        [(inputnode, map_cost_function, [("leadfield", "leadfield")])])
    workflow.connect(
        [(get_elec_indices, map_cost_function, [("ground_index", "ground_index")])])
    workflow.connect(
        [(add_dipole_to_mesh, map_cost_function, [("out_file", "mesh_file")])])
    workflow.connect(
        [(run_forward_model, map_cost_function, [("table_files", "table_files")])])

    workflow.connect(
        [(map_cost_function, outputnode, [("mesh_file", "mesh_file")])])
    workflow.connect(
        [(map_cost_function, outputnode, [("geo_file", "cost_function_geo")])])
    return workflow