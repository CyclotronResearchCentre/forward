import numpy as np
import os.path as op
import subprocess
from forward.mesh import get_num_nodes_elements
import nipype.interfaces.io as nio           # Data i/o
import nipype.interfaces.utility as util     # utility
import nipype.pipeline.engine as pe          # pypeline engine
from nipype.interfaces.gmsh import GetDP


def check_potential(table_files):
    import numpy as np
    from nipype.utils.filemanip import split_filename
    import h5py
    import os.path as op

    '''
    Assess sanity check data to make sure potential at source / sink electrodes
    are 1 and -1 Volts, respectively.
    '''

    field_filename = "e_brain"
    voltage_filename = "v_elec"

    if isinstance(table_files, str):
        table_files = [table_files]

    field_file = ''
    voltage_table_file = ''

    for in_file in table_files:
        _, name, _ = split_filename(in_file)
        if name.rfind(field_filename) >= 0:
            field_file = in_file
            electric_field = np.loadtxt(field_file)
            out_hdf5_file = op.abspath(name + ".hdf5")

        elif name.rfind(voltage_filename) >= 0:
            voltage_table_file = in_file
            potential = np.loadtxt(voltage_table_file)

    if not voltage_table_file or not field_file:
        raise Exception

    ###############
    # Placeholder for now
    sane = True
    if not sane:
        raise Exception("Sanity check on potential failed")
    ###############

    f = h5py.File(out_hdf5_file, "w")
    dset = f.create_dataset("e_field", electric_field.shape, dtype="f", compression="lzf")
    dset[...] = electric_field
    f.close()
    return out_hdf5_file


def write_pro_file(mesh_file, conductivity_tensor_included, source_index, ground_index, electrode_name_file, orig_pro_file=None):
    import numpy as np
    import os
    import os.path as op

    if orig_pro_file is None:
        orig_pro_file = op.join(os.environ["FWD_DIR"], "etc/eeg_forward.pro")

    elec_names = np.loadtxt(electrode_name_file, dtype=str)
    electrode_namelist = elec_names.tolist()
    problem_file = op.abspath("eeg_forward_src_" + electrode_namelist[source_index] + ".pro")

    original = open(orig_pro_file, 'r')
    modified = open(problem_file, 'w')

    while True:
        line = original.readline()
        if len(line) == 0:
            break # EOF
        
        if line.rfind('Anode = Region[') >= 0:
            line = '  Anode = Region[' + str(ground_index+5000) + '];\n'
        elif line.rfind('Cathode = Region[') >= 0:
            line = '  Cathode = Region[' + str(source_index+5000) + '];\n'
        elif line.rfind('sigma[WhiteMatter_Cerebellum]') >= 0 and conductivity_tensor_included:
            line = "    sigma[WhiteMatter_Cerebellum] = TensorField[XYZ[]] #1 ? #1 : 0.33 ;//: 0.33*1e-3 ;\n"
        elif line.rfind('sigma[WhiteMatter_Cerebellum]') >= 0 and not conductivity_tensor_included:
            line = "    sigma[WhiteMatter_Cerebellum] = 0.33;\n"
        modified.write(line)

    modified.close()
    original.close()
    return problem_file


def combine_leadfield_rows(row_data_files, source_indices):
    import h5py
    from nipype.utils.filemanip import split_filename
    import os.path as op
    import numpy as np

    data_name = "e_field"
    '''
    Read electric field results files and append into matrix
    '''
    #field_file = 'e_eeg_forward.txt'
    #electric_field = np.loadtxt(field_file)
    '''
    The number M is the total non-ground recording sites
    '''
    M = len(source_indices)

    '''
    The number N is the total elements in the mesh file
    '''
    field_data_file = h5py.File(row_data_files[0], "r")
    data = field_data_file.get(data_name)
    electric_field = data.value
    print(np.shape(electric_field[0]))
    N = np.shape(electric_field[0])

    '''
    Pre-allocate the Lead Field matrix, which is N rows x M columns
    '''
    N = np.shape(electric_field)[0] * 3

    '''
    Write a single row in the lead field matrix (L_e)
    X,Y,Z elements of the field are flattened, so each row is
    e.g.:
            [E1x E1y E1z E2x E2y E2z]
    '''
    name ="leadfield"
    out_filename = op.abspath(name + ".hdf5")
    out_leadfield_file = h5py.File(out_filename, "w")
    dset = out_leadfield_file.create_dataset("leadfield", (N, M), dtype="f", compression="lzf")
    for index, electric_field_file in enumerate(row_data_files):
        print("Reading field file: %s" % electric_field_file)
        field_data_file = h5py.File(electric_field_file, "r")
        data = field_data_file.get(data_name)
        electric_field = data.value
        src_id = source_indices[index]
        dset[:,src_id-1] = np.ravel(electric_field[:, -3:])

    out_leadfield_file.close()
    print("Saved leadfield matrix as %s" % out_filename)
    return out_filename


def get_electrode_indices(electrode_name_file, marker_list, ground_electrode="IZ"):
    import numpy as np
    elec_names = np.loadtxt(electrode_name_file, dtype=str)
    electrode_namelist = elec_names.tolist()
    elec_indices = range(0, len(electrode_namelist))
    for marker in marker_list:
        try:
            mark_ind = electrode_namelist.index(marker)
        except ValueError:
            raise Exception("Marker electrode not found")
        elec_indices.remove(mark_ind)

    ground_index = electrode_namelist.index(ground_electrode)
    elec_indices.remove(ground_index)
    return elec_indices, ground_index


def create_forward_model_workflow(name, conductivity_tensor_included=False):
    '''
    Define the problem file
    '''

    workflow = pe.Workflow(name=name)
    inputnode = pe.Node(
        interface=util.IdentityInterface(fields=["mesh_file", "electrode_name_file",
            "marker_list", "ground_electrode", "conductivity_tensor_included"]), name="inputnode")
    outputnode = pe.Node(
        interface=util.IdentityInterface(fields=["leadfield"]), name="outputnode")

    get_elec_indices_interface = util.Function(
        input_names=["electrode_name_file", "marker_list", "ground_electrode"],
        output_names=["electrode_indices", "ground_index"], function=get_electrode_indices)

    get_elec_indices = pe.Node(
        interface=get_elec_indices_interface, name="get_elec_indices")

    # Define interface to rewrite the Problem for each electrode
    write_problem_file_interface = util.Function(
        input_names=["mesh_file", "conductivity_tensor_included", "source_index", "ground_index", "electrode_name_file", "orig_pro_file"],
        output_names=["problem_file"], function=write_pro_file)

    write_problem_file = pe.MapNode(interface=write_problem_file_interface,
                                    name="write_problem_file", iterfield=["source_index"])

    # Define MapNode to run the forward model M times
    run_forward_model = pe.MapNode(interface=GetDP(), name="run_forward_model",
                                   iterfield=["problem_file"])
    run_forward_model.inputs.solve = "Electrostatics"
    #run_forward_model.inputs.binary_output_files = True #for debugging
    run_forward_model.inputs.binary_output_files = False #for debugging
    run_forward_model.inputs.out_table_filenames = ["v_elec", "e_brain"]
    run_forward_model.inputs.out_pos_filenames = ["v", "j", "e"]

    # Sanity check for results of forward model for each electrode
    sanity_check_voltage_interface = util.Function(
        input_names=["table_files"],
        output_names=["out_hdf5_file"], function=check_potential)

    run_sanity_check = pe.MapNode(interface=sanity_check_voltage_interface,
                                  iterfield=["table_files"], name="run_sanity_check")

    # Define interface to combine all resulting rows of the forward model
    create_leadfield_interface = util.Function(
        input_names=["row_data_files", "source_indices"],
        output_names=["leadfield_matrix_file"], function=combine_leadfield_rows)

    create_leadfield = pe.Node(
        interface=create_leadfield_interface, name="create_leadfield")

    # Start connecting the inputnode and the electrode identification function
    workflow.connect(
        [(inputnode, get_elec_indices, [("electrode_name_file", "electrode_name_file")])])
    workflow.connect(
        [(inputnode, get_elec_indices, [("marker_list", "marker_list")])])
    workflow.connect(
        [(inputnode, get_elec_indices, [("ground_electrode", "ground_electrode")])])

    # This connection spawns multiple instances because write_problem_file is
    # a MapNode with iterfield source_index
    workflow.connect(
        [(get_elec_indices, write_problem_file, [("electrode_indices", "source_index")])])
    workflow.connect(
        [(get_elec_indices, write_problem_file, [("ground_index", "ground_index")])])
    workflow.connect(
        [(inputnode, write_problem_file, [("electrode_name_file", "electrode_name_file")])])

    workflow.connect(
        [(inputnode, write_problem_file, [("mesh_file", "mesh_file")])])
    workflow.connect(
        [(inputnode, write_problem_file, [("conductivity_tensor_included", "conductivity_tensor_included")])])
    workflow.connect(
        [(write_problem_file, run_forward_model, [("problem_file", "problem_file")])])
    workflow.connect(
        [(inputnode, run_forward_model, [("mesh_file", "mesh_file")])])
    if conductivity_tensor_included:
        workflow.connect(
            [(inputnode, run_forward_model, [("mesh_file", "gmsh_read_file")])])

    workflow.connect(
        [(run_forward_model, run_sanity_check, [("table_files", "table_files")])])

    # Now all the MapNodes are combined again to provide a list of
    # row_data_files to the create_leadfield interface
    workflow.connect(
        [(run_sanity_check, create_leadfield, [("out_hdf5_file", "row_data_files")])])
    workflow.connect(
        [(get_elec_indices, create_leadfield, [("electrode_indices", "source_indices")])])
    workflow.connect(
        [(create_leadfield, outputnode, [("leadfield_matrix_file", "leadfield")])])
    return workflow
