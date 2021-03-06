import nipype.interfaces.utility as util     # utility
import nipype.pipeline.engine as pe          # pypeline engine
from nipype.interfaces.gmsh import GetDP


def check_potential(table_files):
    import numpy as np
    from nipype.utils.filemanip import split_filename
    import h5py
    import os.path as op

    '''
    Sanity check data to make sure potential at source / sink electrodes
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
        
        if line.rfind('Sink = Region[') >= 0:
            line = '  Sink = Region[' + str(ground_index+5000) + '];\n'
        elif line.rfind('Source = Region[') >= 0:
            line = '  Source = Region[' + str(source_index+5000) + '];\n'
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
    run_forward_model.inputs.binary_output_files = True #for debugging
    
    run_forward_model.inputs.out_table_filenames = ["v_elec", "e_brain"]
    run_forward_model.inputs.out_pos_filenames = ["v_elec", "e_brain"]

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

def compare_leadfields(leadfield1, leadfield2, mesh_file, write_mesh=True):
    '''
    Compares two leadfield matrices and outputs the difference as
    scalars attached to a mesh file. The mesh file must have the 
    same number of elements as the leadfield's second dimension

    FEM reciprocity leadfield is M x N where M is the number of sensors
    and N is the number of elements. Mesh file must have N elements.

    e.g. isotropic white matter conductivity
    vs. anisotropic conductivity from estimated diffusion tensors
    '''
    import os.path as op
    import h5py
    import numpy as np
    import shutil
    
    from forward.mesh import read_mesh

    from nipype.utils.filemanip import split_filename

    from nipype import logging    
    iflogger = logging.getLogger('interface')

    data_name = "leadfield"

    print("Reading leadfield 1: %s" % leadfield1)
    lf1_data_file = h5py.File(leadfield1, "r")
    lf1_data = lf1_data_file.get(data_name)
    leadfield_matrix1 = lf1_data.value

    print("Reading leadfield 2: %s" % leadfield2)
    lf2_data_file = h5py.File(leadfield2, "r")
    lf2_data = lf2_data_file.get(data_name)
    leadfield_matrix2 = lf2_data.value


    path, name, ext = split_filename(mesh_file)

    if write_mesh:
        # Electric field elements are only saved in the gray matter
        #elements_to_consider = [1001] #For Sphere
        elements_to_consider = [1002] #For head models
        mesh_data, _, _, _ = read_mesh(mesh_file, elements_to_consider)
        # Create the output mesh file
        rms_mesh_file = op.abspath(name + "_rmse.msh")
        iflogger.info('Copying current mesh file to %s' % rms_mesh_file)
        shutil.copyfile(mesh_file, rms_mesh_file)
        f = open(rms_mesh_file,'a') #Append to the end of the file
        iflogger.info('Appending root mean squared error scalars to %s' % rms_mesh_file)
    else:
        rms_mesh_file = None

    out_rms_hdf5_file_x = op.abspath(name + "_rmse_x.hdf5")
    out_rms_hdf5_file_y = op.abspath(name + "_rmse_y.hdf5")
    out_rms_hdf5_file_z = op.abspath(name + "_rmse_z.hdf5")
    out_rms_hdf5_file_avg = op.abspath(name + "_rmse_avg.hdf5")

    rms_hdf5_file_x = h5py.File(out_rms_hdf5_file_x, "w")
    rms_hdf5_file_y = h5py.File(out_rms_hdf5_file_y, "w")
    rms_hdf5_file_z = h5py.File(out_rms_hdf5_file_z, "w")
    rms_hdf5_file_avg = h5py.File(out_rms_hdf5_file_avg, "w")


    # Write the tag information to the file:
    if write_mesh:
        num_polygons = len(mesh_data)
    _, name1, _ = split_filename(leadfield1)
    _, name2, _ = split_filename(leadfield2)

    #For testing
    #noise = np.random.normal(0,1,leadfield_matrix2.shape)
    #leadfield_matrix2 = leadfield_matrix2 + noise

    # Reshape the leadfields so they are M x 3 x N vectors (x, y, z direction of electric field)
    leadfield_mesh1 = np.reshape(leadfield_matrix1, (leadfield_matrix1.shape[0]/3,3,leadfield_matrix1.shape[1]))
    leadfield_mesh2 = np.reshape(leadfield_matrix2, (leadfield_matrix2.shape[0]/3,3,leadfield_matrix2.shape[1]))

    #Check that the dimensions are appropriate
    try:
        if write_mesh:
            assert(len(mesh_data) == leadfield_mesh1.shape[0] == leadfield_mesh2.shape[0])
        else:
            assert(leadfield_mesh1.shape[0] == leadfield_mesh2.shape[0])
    except AssertionError:
        iflogger.error("Lead fields could not be compared because the " \
            "number of elements in the files do not match")
        if write_mesh:
            iflogger.error("Elements in %s: %d" % (mesh_file, len(mesh_data)))
        iflogger.error("Elements in leadfield 1 %s: %d" % (leadfield1, leadfield_mesh1.shape[0]))
        iflogger.error("Elements in leadfield 2 %s: %d" % (leadfield2, leadfield_mesh2.shape[0]))
        raise Exception

    # Get the mean squared difference between electric field vectors by sensor and element
    diff = leadfield_mesh2 - leadfield_mesh1

    # Gets the norm of the difference for X, Y, and Z directions
    norm_x_diff = np.linalg.norm(diff[:,0,:],axis=1)
    norm_y_diff = np.linalg.norm(diff[:,1,:],axis=1)
    norm_z_diff = np.linalg.norm(diff[:,2,:],axis=1)

    norm_lf1_x = np.linalg.norm(leadfield_mesh1[:,0,:], axis=1)
    norm_lf1_y = np.linalg.norm(leadfield_mesh1[:,1,:], axis=1)
    norm_lf1_z = np.linalg.norm(leadfield_mesh1[:,2,:], axis=1)

    rmse_x = norm_x_diff / norm_lf1_x
    rmse_y = norm_y_diff / norm_lf1_y
    rmse_z = norm_z_diff / norm_lf1_z

    rmse = (rmse_x + rmse_y + rmse_z) / 3.
    assert(leadfield_mesh2.shape[0] == rmse.shape[0])

    if write_mesh:
        rmse_avg_list = []
        rmse_x_list = []
        rmse_y_list = []
        rmse_z_list = []
        for idx, poly in enumerate(mesh_data):
            rmse_avg_str = ('%d %e \n' % (poly["element_id"], rmse[idx]))
            rmse_avg_list.append(rmse_avg_str)
            rmse_x_str = ('%d %e \n' % (poly["element_id"], rmse_x[idx]))
            rmse_x_list.append(rmse_x_str)
            rmse_y_str = ('%d %e \n' % (poly["element_id"], rmse_y[idx]))
            rmse_y_list.append(rmse_y_str)
            rmse_z_str  = ('%d %e \n' % (poly["element_id"], rmse_z[idx]))        
            rmse_z_list.append(rmse_z_str)
            iflogger.info("%3.3f%%" % (float(idx)/num_polygons*100.0))


        # ----- Average RMSE ----- #
        f.write('$ElementData\n')
        str_tag = '"Average Root Mean Squared Error"'
        timestep = 0.0001
        f.write('1\n') #Num String tags
        f.write(str_tag + '\n')
        f.write('1\n') #Num Real tags
        f.write('%f\n' % timestep)
        #Three integer tags: timestep, num field components, num elements
        f.write('3\n') #Three int tags
        f.write('0\n') #Time step index
        f.write('1\n') #Num field components
        f.write('%d\n' % len(mesh_data)) #Num nonzero field components
        for elementdata_str in rmse_avg_list:
            f.write(elementdata_str)
        f.write('$EndElementData\n')


        # ----- RMSE X ----- #
        f.write('$ElementData\n')
        str_tag = '"Root Mean Squared Error X-direction"'
        timestep = 0.0002
        f.write('1\n') #Num String tags
        f.write(str_tag + '\n')
        f.write('1\n') #Num Real tags
        f.write('%f\n' % timestep)
        #Three integer tags: timestep, num field components, num elements
        f.write('3\n') #Three int tags
        f.write('0\n') #Time step index
        f.write('1\n') #Num field components
        f.write('%d\n' % len(mesh_data)) #Num nonzero field components
        for elementdata_str in rmse_x_list:
            f.write(elementdata_str)
        f.write('$EndElementData\n')

        # ----- RMSE X ----- #
        f.write('$ElementData\n')
        str_tag = '"Root Mean Squared Error Y-direction"'
        timestep = 0.0003
        f.write('1\n') #Num String tags
        f.write(str_tag + '\n')
        f.write('1\n') #Num Real tags
        f.write('%f\n' % timestep)
        #Three integer tags: timestep, num field components, num elements
        f.write('3\n') #Three int tags
        f.write('0\n') #Time step index
        f.write('1\n') #Num field components
        f.write('%d\n' % len(mesh_data)) #Num nonzero field components
        for elementdata_str in rmse_y_list:
            f.write(elementdata_str)
        f.write('$EndElementData\n')

        # ----- RMSE X ----- #
        f.write('$ElementData\n')
        str_tag = '"Root Mean Squared Error Z-direction"'
        timestep = 0.0004
        f.write('1\n') #Num String tags
        f.write(str_tag + '\n')
        f.write('1\n') #Num Real tags
        f.write('%f\n' % timestep)
        #Three integer tags: timestep, num field components, num elements
        f.write('3\n') #Three int tags
        f.write('0\n') #Time step index
        f.write('1\n') #Num field components
        f.write('%d\n' % len(mesh_data)) #Num nonzero field components
        for elementdata_str in rmse_z_list:
            f.write(elementdata_str)
        f.write('$EndElementData\n')

        f.close()

        iflogger.info("Finished writing to %s" % rms_mesh_file)

    ## Save RMSE data to an HDF5 file
    dset_x = rms_hdf5_file_x.create_dataset("rmse_x", data=rmse_x)
    dset_x[...] = rmse_x
    dset_y = rms_hdf5_file_y.create_dataset("rmse_y", data=rmse_y)
    dset_y[...] = rmse_y
    dset_z = rms_hdf5_file_z.create_dataset("rmse_z", data=rmse_z)
    dset_z[...] = rmse_z
    dset = rms_hdf5_file_avg.create_dataset("rmse_avg", data=rmse)
    dset[...] = rmse

    rms_hdf5_file_x.close()
    rms_hdf5_file_y.close()
    rms_hdf5_file_z.close()
    rms_hdf5_file_avg.close()
    print("Saved RMSE-X matrix as %s" % out_rms_hdf5_file_x)
    print("Saved RMSE-Y matrix as %s" % out_rms_hdf5_file_y)
    print("Saved RMSE-Z matrix as %s" % out_rms_hdf5_file_z)
    print("Saved RMSE-Avg matrix as %s" % out_rms_hdf5_file_avg)
    ###
    return rms_mesh_file, out_rms_hdf5_file_avg, out_rms_hdf5_file_x, out_rms_hdf5_file_y, out_rms_hdf5_file_z