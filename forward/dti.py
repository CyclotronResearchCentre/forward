import nipype.interfaces.utility as util     # utility
import nipype.pipeline.engine as pe          # pypeline engine
import nipype.interfaces.fsl as fsl
import nipype.interfaces.dipy as dipy
import shutil

from nipype import logging
iflogger = logging.getLogger('interface')

fsl.FSLCommand.set_default_output_type('NIFTI')


def split_tissue_classes_fn(in_files):
    for in_file in in_files:
        if 'seg_0' in in_file:
            gray_matter = in_file
        elif 'seg_1' in in_file:
            white_matter = in_file
    return gray_matter, white_matter


def calc_scale_factor_fn(white_matter, gray_matter):
    import numpy as np
    conductivity_WM = 0.126  # [S/m]
    conductivity_GM = 0.276  # [S/m]

    scale_factor = ((white_matter * conductivity_WM + gray_matter * conductivity_GM) /
                   ((np.power(white_matter, 2) + np.power(gray_matter, 2)) * 1000))
    return scale_factor


def create_get_scale_factor_workflow(name):
    inputnode = pe.Node(
        util.IdentityInterface(
            fields=["tissue_class_files", "L1", "L2", "L3"]),
        name="inputnode")
    split_tissue_classes_interface = util.Function(input_names=["in_files"],
                                                   output_names=[
                                                       "gray_matter", "white_matter"],
                                                   function=split_tissue_classes_fn)
    split_tissue_classes = pe.Node(
        interface=split_tissue_classes_interface, name='split_tissue_classes')

    calc_scale_factor_interface = util.Function(
        input_names=["gray_matter", "white_matter"],
        output_names=["scale_factor"],
        function=calc_scale_factor_fn)
    calc_scale_factor = pe.Node(
        interface=calc_scale_factor_interface, name='calc_scale_factor')

    mult_EV1_GM = pe.Node(
        interface=fsl.MultiImageMaths(), name='mult_EV1_GM')
    mult_EV1_GM.inputs.op_string = "-mul %s -mul 1000"
    mult_EV2_GM = mult_EV1_GM.clone("mult_EV2_GM")
    mult_EV3_GM = mult_EV1_GM.clone("mult_EV3_GM")
    GM_EV_merge = pe.Node(interface=util.Merge(2), name='GM_EV_merge')

    mult_EV1_WM = mult_EV1_GM.clone("mult_EV1_WM")
    mult_EV2_WM = mult_EV1_GM.clone("mult_EV2_WM")
    mult_EV3_WM = mult_EV1_GM.clone("mult_EV3_WM")
    WM_EV_merge = pe.Node(interface=util.Merge(2), name='WM_EV_merge')

    multiply_GM = pe.Node(
        interface=fsl.MultiImageMaths(), name='multiply_GM')
    multiply_GM.inputs.op_string = "-mul %s -mul %s"
    multiply_WM = multiply_GM.clone("multiply_WM")

    GM_stats = pe.Node(interface=fsl.ImageStats(), name='GM_stats')
    GM_stats.inputs.op_string = '-M'
    WM_stats = GM_stats.clone('WM_stats')

    workflow = pe.Workflow(name=name)

    workflow.connect(
        [(inputnode, split_tissue_classes, [("tissue_class_files", "in_files")])])

    workflow.connect([(inputnode, mult_EV1_GM, [("L1", "in_file")])])
    workflow.connect(
        [(split_tissue_classes, mult_EV1_GM, [("gray_matter", "operand_files")])])

    workflow.connect([(inputnode, mult_EV2_GM, [("L2", "in_file")])])
    workflow.connect(
        [(split_tissue_classes, mult_EV2_GM, [("gray_matter", "operand_files")])])

    workflow.connect([(inputnode, mult_EV3_GM, [("L3", "in_file")])])
    workflow.connect(
        [(split_tissue_classes, mult_EV3_GM, [("gray_matter", "operand_files")])])

    workflow.connect([(inputnode, mult_EV1_WM, [("L1", "in_file")])])
    workflow.connect(
        [(split_tissue_classes, mult_EV1_WM, [("white_matter", "operand_files")])])

    workflow.connect([(inputnode, mult_EV2_WM, [("L2", "in_file")])])
    workflow.connect(
        [(split_tissue_classes, mult_EV2_WM, [("white_matter", "operand_files")])])

    workflow.connect([(inputnode, mult_EV3_WM, [("L3", "in_file")])])
    workflow.connect(
        [(split_tissue_classes, mult_EV3_WM, [("white_matter", "operand_files")])])

    workflow.connect([(mult_EV1_GM, multiply_GM, [("out_file", "in_file")])])
    workflow.connect([(mult_EV2_GM, GM_EV_merge, [("out_file", "in1")])])
    workflow.connect([(mult_EV3_GM, GM_EV_merge, [("out_file", "in2")])])
    workflow.connect([(GM_EV_merge, multiply_GM, [("out", "operand_files")])])

    workflow.connect([(mult_EV1_WM, multiply_WM, [("out_file", "in_file")])])
    workflow.connect([(mult_EV2_WM, WM_EV_merge, [("out_file", "in1")])])
    workflow.connect([(mult_EV3_WM, WM_EV_merge, [("out_file", "in2")])])
    workflow.connect([(WM_EV_merge, multiply_WM, [("out", "operand_files")])])

    workflow.connect([(multiply_GM, GM_stats, [("out_file", "in_file")])])
    workflow.connect([(multiply_WM, WM_stats, [("out_file", "in_file")])])

    workflow.connect(
        [(WM_stats, calc_scale_factor, [("out_stat", "white_matter")])])
    workflow.connect(
        [(GM_stats, calc_scale_factor, [("out_stat", "gray_matter")])])

    return workflow


def include_gmsh_tensor_elements(mesh_file, tensor_file, mask_file, mask_threshold=0.3):
    import numpy as np
    import nibabel as nb
    from nipype.utils.filemanip import split_filename
    import os.path as op
    from forward.mesh import read_mesh
    from nipype import logging
    import shutil
    iflogger = logging.getLogger('interface')

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
    elements_to_consider = [1001] #Use only white matter
    vx, vy, vz = header.get_zooms()[0:3]
    max_x, max_y, max_z = np.shape(tensor_data)[0:3]
    halfx, halfy, halfz = np.array((vx*max_x, vy*max_y, vz*max_z))/2.0

    mesh_data = read_mesh(mesh_file, elements_to_consider)

    # Create the output mesh file
    path, name, ext = split_filename(mesh_file)
    out_file = op.abspath(name + "_cond.msh")
    iflogger.info('Copying current mesh file to %s' % out_file)
    shutil.copyfile(mesh_file, out_file)

    f = open(out_file,'a') #Append to the end of the file
    iflogger.info('Appending Conductivity tensors to %s' % out_file)

    # Write the tag information to the file:
    num_polygons = len(mesh_data)
    f.write('$ElementData\n')
    str_tag = '"Conductivity"'
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
    iflogger.info("Getting tensor for each element")
    #ipdb.set_trace()
    nonzero = 0
    elem_list = []
    for idx, poly in enumerate(mesh_data):
        i = np.round((poly['centroid'][0]+halfx)/vx).astype(int)
        j = np.round((poly['centroid'][1]+halfy)/vy).astype(int)
        k = np.round((poly['centroid'][2]+halfz)/vz).astype(int)
        poly["tensor_triform"] = tensor_data[i,j,k]
        if not (all(poly["tensor_triform"] == 0) and mask_data[i,j,k] >= mask_threshold):
            elementdata_str = ('%d %e %e %e %e %e %e %e %e %e\n' % (poly["element_id"], 
                poly["tensor_triform"][0], poly["tensor_triform"][1], poly["tensor_triform"][2],
                poly["tensor_triform"][1], poly["tensor_triform"][3], poly["tensor_triform"][4],
                poly["tensor_triform"][3], poly["tensor_triform"][4], poly["tensor_triform"][5]))
            elem_list.append(elementdata_str)
            nonzero += 1
        iflogger.info("%3.3f%%" % (float(idx)/num_polygons*100.0))

    f.write('%d\n' % nonzero) #Num nonzero field components
    for elementdata_str in elem_list:
        f.write(elementdata_str)

    f.write('$EndElementData\n')
    f.close()

    iflogger.info("Finished writing to %s" % out_file)
    return out_file


def create_conductivity_tensor_mesh_workflow(name="add_conductivity"):
    inputnode = pe.Node(interface=util.IdentityInterface(
        fields=["dwi", "bvecs", "bvals", "struct", "mesh_file"]), name="inputnode")
    outputnode = pe.Node(
        interface=util.IdentityInterface(fields=["mesh_file"]), name="outputnode")

    bet_T1 = pe.Node(interface=fsl.BET(), name='bet_T1')

    bet_b0 = pe.Node(interface=fsl.BET(), name='bet_b0')
    bet_b0.inputs.frac = 0.2
    bet_b0.inputs.mask = True

    dtifit = pe.Node(interface=fsl.DTIFit(), name='dtifit')
    dtifit.inputs.save_tensor = True

    affine_register = pe.Node(interface=fsl.FLIRT(), name='affine_register')
    affine_register.inputs.cost = ('normmi')
    affine_register.inputs.dof = 6

    nonlinear_warp = pe.Node(interface=fsl.FNIRT(), name='nonlinear_warp')
    #nonlinear_warp.inputs.field_file = True
    nonlinear_warp.inputs.fieldcoeff_file = "FA_T1_warp.nii"

    register_tensor = pe.Node(interface=fsl.VecReg(), name='register_tensor')

    decompose_tensor = pe.Node(
        interface=fsl.DecomposeTensor(), name='decompose_tensor')

    fast_segment = pe.Node(interface=fsl.FAST(), name='fast_segment')
    fast_segment.inputs.segments = True

    conductivity_mapping = pe.Node(
        interface=dipy.EstimateConductivity(), name='conductivity_mapping')

    add_tensors_to_mesh_interface = util.Function(
        input_names=["mesh_file", "tensor_file", "mask_file", "mask_threshold"],
        output_names=["out_file"],
        function=include_gmsh_tensor_elements)
    add_conductivity_tensor_to_mesh = pe.Node(
        interface=add_tensors_to_mesh_interface, name='add_conductivity_tensor_to_mesh')

    erode_FA_map = pe.Node(interface=fsl.BET(), name="erode_FA_map")
    erode_FA_map.inputs.mask = True

    """
    Define the Workflow
    -------------------
    """

    workflow = pe.Workflow(name='workflow')
    scale_factor_workflow = create_get_scale_factor_workflow(
        "scale_factor_workflow")

    workflow.connect([(inputnode, dtifit, [('dwi', 'dwi'),
                    ('bvals', 'bvals'),
                    ('bvecs', 'bvecs'),
    ])])

    workflow.connect([(inputnode, bet_b0, [('dwi', 'in_file')])])
    workflow.connect([(bet_b0, dtifit, [('mask_file', 'mask')])])

    workflow.connect([(inputnode, bet_T1, [('struct', 'in_file')])])
    workflow.connect([(bet_T1, affine_register, [('out_file', 'reference')])])

    workflow.connect([(dtifit, affine_register, [('FA', 'in_file')])])
    workflow.connect(
        [(affine_register, nonlinear_warp, [('out_file', 'in_file')])])
    workflow.connect([(bet_T1, nonlinear_warp, [('out_file', 'ref_file')])])

    workflow.connect([(dtifit, register_tensor, [('tensor', 'in_file')])])
    workflow.connect([(bet_T1, register_tensor, [('out_file', 'ref_vol')])])
    workflow.connect(
        [(affine_register, nonlinear_warp, [('out_matrix_file', 'affine_file')])])
    workflow.connect(
        [(nonlinear_warp, register_tensor, [('fieldcoeff_file', 'warp_field')])])
    workflow.connect(
        [(register_tensor, decompose_tensor, [('out_file', 'in_file')])])

    workflow.connect([(bet_T1, fast_segment, [("out_file", "in_files")])])

    workflow.connect(
        [(fast_segment, scale_factor_workflow, [("tissue_class_files", "inputnode.tissue_class_files")])])
    workflow.connect(
        [(decompose_tensor, scale_factor_workflow, [('L1', 'inputnode.L1'),
                    ('L2', 'inputnode.L2'),
                    ('L3', 'inputnode.L3'),
          ])])

    workflow.connect([(scale_factor_workflow, conductivity_mapping,
                       [("calc_scale_factor.scale_factor", "eigenvalue_scaling_factor")])])

    workflow.connect([(register_tensor, conductivity_mapping,
                       [("out_file", "in_file")])])

    workflow.connect(
        [(conductivity_mapping, add_conductivity_tensor_to_mesh, [("out_file", "tensor_file")])])
    workflow.connect(
        [(nonlinear_warp, add_conductivity_tensor_to_mesh, [("warped_file", "mask_file")])])
    workflow.connect(
        [(inputnode, add_conductivity_tensor_to_mesh, [("mesh_file", "mesh_file")])])


    workflow.connect(
        [(add_conductivity_tensor_to_mesh, outputnode, [("out_file", "mesh_file")])])
    return workflow
