import nipype.interfaces.io as nio           # Data i/o
import nipype.interfaces.utility as util     # utility
import nipype.pipeline.engine as pe          # pypeline engine
import nipype.interfaces.fsl as fsl
import nipype.interfaces.dipy as dipy
import os.path as op                     # system functions
import ipdb

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
    conductivity_WM = 0.126 # [S/m]
    conductivity_GM = 0.276 # [S/m]

    scale_factor = ((white_matter * conductivity_WM + gray_matter * conductivity_GM) / 
        ((np.power(white_matter,2) + np.power(gray_matter,2))*1000))
    return scale_factor

def create_get_scale_factor_workflow(name):
    inputnode = pe.Node(util.IdentityInterface(fields=["tissue_class_files", "L1", "L2", "L3"]),
                               name="inputnode")
    split_tissue_classes_interface = util.Function(input_names=["in_files"],
                                output_names=["gray_matter", "white_matter"],
                                function=split_tissue_classes_fn)
    split_tissue_classes = pe.Node(interface=split_tissue_classes_interface, name='split_tissue_classes')

    calc_scale_factor_interface = util.Function(input_names=["gray_matter", "white_matter"],
                                output_names=["scale_factor"],
                                function=calc_scale_factor_fn)
    calc_scale_factor = pe.Node(interface=calc_scale_factor_interface, name='calc_scale_factor')

    mult_EV1_GM = pe.Node(interface=fsl.MultiImageMaths(), name = 'mult_EV1_GM')
    mult_EV1_GM.inputs.op_string = "-mul %s -mul 1000"
    mult_EV2_GM = mult_EV1_GM.clone("mult_EV2_GM")
    mult_EV3_GM = mult_EV1_GM.clone("mult_EV3_GM")
    GM_EV_merge = pe.Node(interface=util.Merge(2), name='GM_EV_merge')

    mult_EV1_WM = mult_EV1_GM.clone("mult_EV1_WM")
    mult_EV2_WM = mult_EV1_GM.clone("mult_EV2_WM")
    mult_EV3_WM = mult_EV1_GM.clone("mult_EV3_WM")
    WM_EV_merge = pe.Node(interface=util.Merge(2), name='WM_EV_merge')

    multiply_GM = pe.Node(interface=fsl.MultiImageMaths(), name = 'multiply_GM')
    multiply_GM.inputs.op_string = "-mul %s -mul %s"
    multiply_WM = multiply_GM.clone("multiply_WM")

    GM_stats = pe.Node(interface=fsl.ImageStats(), name = 'GM_stats')
    GM_stats.inputs.op_string = '-M'
    WM_stats= GM_stats.clone('WM_stats')

    workflow = pe.Workflow(name=name)

    workflow.connect([(inputnode, split_tissue_classes,[("tissue_class_files","in_files")])])

    workflow.connect([(inputnode, mult_EV1_GM,[("L1","in_file")])])
    workflow.connect([(split_tissue_classes, mult_EV1_GM,[("gray_matter","operand_files")])])

    workflow.connect([(inputnode, mult_EV2_GM,[("L2","in_file")])])
    workflow.connect([(split_tissue_classes, mult_EV2_GM,[("gray_matter","operand_files")])])

    workflow.connect([(inputnode, mult_EV3_GM,[("L3","in_file")])])
    workflow.connect([(split_tissue_classes, mult_EV3_GM,[("gray_matter","operand_files")])])

    workflow.connect([(inputnode, mult_EV1_WM,[("L1","in_file")])])
    workflow.connect([(split_tissue_classes, mult_EV1_WM,[("white_matter","operand_files")])])

    workflow.connect([(inputnode, mult_EV2_WM,[("L2","in_file")])])
    workflow.connect([(split_tissue_classes, mult_EV2_WM,[("white_matter","operand_files")])])

    workflow.connect([(inputnode, mult_EV3_WM,[("L3","in_file")])])
    workflow.connect([(split_tissue_classes, mult_EV3_WM,[("white_matter","operand_files")])])

    workflow.connect([(mult_EV1_GM, multiply_GM,[("out_file","in_file")])])
    workflow.connect([(mult_EV2_GM, GM_EV_merge,[("out_file","in1")])])
    workflow.connect([(mult_EV3_GM, GM_EV_merge,[("out_file","in2")])])
    workflow.connect([(GM_EV_merge, multiply_GM,[("out","operand_files")])])

    workflow.connect([(mult_EV1_WM, multiply_WM,[("out_file","in_file")])])
    workflow.connect([(mult_EV2_WM, WM_EV_merge,[("out_file","in1")])])
    workflow.connect([(mult_EV3_WM, WM_EV_merge,[("out_file","in2")])])
    workflow.connect([(WM_EV_merge, multiply_WM,[("out","operand_files")])])

    workflow.connect([(multiply_GM, GM_stats,[("out_file","in_file")])])
    workflow.connect([(multiply_WM, WM_stats,[("out_file","in_file")])])

    workflow.connect([(WM_stats, calc_scale_factor,[("out_stat","white_matter")])])
    workflow.connect([(GM_stats, calc_scale_factor,[("out_stat","gray_matter")])])

    return workflow


def nifti_tensors_to_gmsh(in_file, mask_file, fa_file, threshold=0.5):
    import numpy as np
    import nibabel as nb
    from nipype.utils.filemanip import split_filename
    import os.path as op
    tensor_image = nb.load(in_file)
    mask_image = nb.load(mask_file)
    fa_image = nb.load(fa_file)
    path, name, ext = split_filename(in_file)
    out_file = op.abspath(name + '.geo')
    f = open(out_file,'w')
    print 'Writing tensors to {f}'.format(f=out_file)
    data = tensor_image.get_data()
    mask_data = mask_image.get_data()
    fa_data = fa_image.get_data()
    header = tensor_image.get_header()
    zooms = header.get_zooms()
    vx = zooms[0]
    vy = zooms[1]
    vz = zooms[2]
    halfx = vx*np.shape(data)[0]/2.0
    halfy = vy*np.shape(data)[1]/2.0
    halfz = vz*np.shape(data)[2]/2.0
    max_x = np.shape(data)[0]
    max_y = np.shape(data)[1]
    max_z = np.shape(data)[2]
    view_str = 'View "Diffusion tensor field Fractional Anisotropy >= %f"{' % (threshold)
    f.write(view_str + '\n')
    for x in range(0,max_x):
        for y in range(0,max_y):
            for z in range(0,max_z):
                if mask_data[x,y,z] > 0 and fa_data[x,y,z] >= threshold:
                    tensor = np.zeros((3,3))
                    tensor[0,0] = data[x,y,z,0]
                    tensor[1,0] = data[x,y,z,1]
                    tensor[1,1] = data[x,y,z,2]
                    tensor[2,0] = data[x,y,z,3]
                    tensor[2,1] = data[x,y,z,4]
                    tensor[2,2] = data[x,y,z,5]
                    tensor = tensor + tensor.T
                    t_list = list(tensor.flatten())
                    pt_str = ('TP(%3.2f,%3.2f,%3.2f){%e,%e,%e,%e,%e,%e,%e,%e,%e};' % (x*vx-halfx, y*vy-halfy, z*vz-halfz, t_list[0], 
                    t_list[1], t_list[2], t_list[3], t_list[4], t_list[5], t_list[6], t_list[7], t_list[8]))
                    f.write(pt_str + '\n')
    closing_str = '};'
    f.write(closing_str + '\n')
    f.close()
    return out_file

def nifti_to_gmsh_tensor_elements(in_file, mask_file, fa_file, mask_threshold=0.1, fa_threshold=0.3):
    import numpy as np
    import nibabel as nb
    from nipype.utils.filemanip import split_filename
    import os.path as op
    tensor_image = nb.load(in_file)
    mask_image = nb.load(mask_file)
    fa_image = nb.load(fa_file)
    path, name, ext = split_filename(in_file)
    out_file = op.abspath(name + "_elem.msh")

    f = open(out_file,'w')
    print 'Writing tensors to {f}'.format(f=out_file)
    data = tensor_image.get_data()
    fa_data = fa_image.get_data()
    mask_data = mask_image.get_data()
    data = np.flipud(data)
    mask_data = np.flipud(mask_data)
    fa_data = np.flipud(fa_data)
    header = tensor_image.get_header()
    zooms = header.get_zooms()
    elements = []
    tensors = {}
    node_positions = {}

    node_id = 0
    element_id = 10000000
    vx = zooms[0]
    vy = zooms[1]
    vz = zooms[2]
    halfx = vx*np.shape(data)[0]/2.0
    halfy = vy*np.shape(data)[1]/2.0
    halfz = vz*np.shape(data)[2]/2.0
    max_x = np.shape(data)[0]
    max_y = np.shape(data)[1]
    max_z = np.shape(data)[2]

    print 'Writing the header for the file...'
    intro_list = ['$MeshFormat']
    intro_list.append('2.0 0 8')
    intro_list.append('$EndMeshFormat')
    intro_list.append('$Nodes')
    for intro_str in intro_list:
        f.write(intro_str + '\n')

    node_info = []
    nodes = np.zeros((np.shape(data)[0:3]))
    print 'Writing the position of each node...'
    print np.shape(data)
    for x in range(0,max_x):
        for y in range(0,max_y):
            for z in range(0,max_z):
                node_id += 1
                nodes[x,y,z] = node_id
                node_str = ('%d %f %f %f' % (node_id, x*vx-halfx, y*vy-halfy, z*vz-halfz))
                node_info.append(node_str)

    n_node_str = ('%d' % (len(node_info)))
    f.write(n_node_str + '\n')
    for node_str in node_info:
        f.write(node_str + '\n')

    print 'Calculating the tensor for each element...'
    for x in xrange(0,np.shape(data)[0]-1):
        for y in xrange(0,np.shape(data)[1]-1):
            for z in xrange(0,np.shape(data)[2]-1):
                if mask_data[x,y,z] >= mask_threshold and fa_data[x,y,z] >= fa_threshold:
                    tensor = np.zeros((3,3))
                    tensor[0,0] = data[x,y,z,0]
                    tensor[0,1] = data[x,y,z,1]
                    tensor[1,0] = data[x,y,z,1]
                    tensor[0,2] = data[x,y,z,2]
                    tensor[2,0] = data[x,y,z,2]
                    tensor[1,1] = data[x,y,z,3]
                    tensor[1,2] = data[x,y,z,4]
                    tensor[2,1] = data[x,y,z,4]
                    tensor[2,2] = data[x,y,z,5]
                    element_id += 1
                    elements.append(element_id)
                    try:
                        node_positions[element_id] = [nodes[x,y,z], nodes[x+1,y,z], nodes[x+1,y+1,z], nodes[x,y+1,z], nodes[x,y,z+1], nodes[x+1,y,z+1], nodes[x+1,y+1,z+1], nodes[x,y+1,z+1]]
                        tensors[element_id] = tensor
                    except IndexError:
                        continue
                
    f.write('$EndNodes\n')

    print 'Write the Elements block'
    f.write('$Elements\n')
    n_elements = ('%f' % (len(node_positions.keys())))
    f.write(n_elements + '\n')
    block_type = 5 # Hexahedron aka Cube
    number_of_tags = 2
    tag_1 = 6
    material_tag = 100
    for element_id in elements:
        try:
            node = node_positions[element_id]
            element_str = ('%d %d %d %d %d %d %d %d %d %d %d %d %d' % (element_id, block_type, number_of_tags, tag_1, material_tag,
            node[0], node[1], node[2], node[3], node[4], node[5], node[6], node[7]))
            f.write(element_str + '\n')
        except KeyError:
            continue
    f.write('$EndElements\n')
    print 'End of the Elements Block'

    f.write('$ElementData\n')
    n_str_tags = str(1)
    str_tag = '"Conductivity"'

    f.write(n_str_tags + '\n')
    f.write(str_tag + '\n')

    n_real_tags = str(1)
    real_tag = str(10)

    f.write(n_real_tags + '\n')
    f.write(real_tag + '\n')

    n_int_tags = str(3)
    int_tags = []
    int_tags.append(3)
    int_tags.append(0)
    int_tags.append(9) # 9 for tensor
    for tag in int_tags:
        f.write(str(tag) + '\n')
        
    elementdata_list = []
    print 'Writing the element data block...'
    for element_id in elements:
        tensor = tensors[element_id]
        t_list = list(tensor.flatten())
        elementdata_str = ('%d %e %e %e %e %e %e %e %e %e' % (element_id, t_list[0], 
        t_list[1], t_list[2], t_list[3], t_list[4], t_list[5], t_list[6], t_list[7], t_list[8]))
        elementdata_list.append(elementdata_str)
        
    n_elementdata = str(len(elementdata_list))
    f.write('%d\n' % n_elementdata)   
    for elementdata_str in elementdata_list:
        f.write(elementdata_str + '\n')
    f.write('$EndElementData\n')
    f.close()

    return out_file


data_dir = op.abspath(op.curdir)
subject_list = ['TMS007']

infosource = pe.Node(interface=util.IdentityInterface(
    fields=['subject_id']), name="infosource")
infosource.iterables = ('subject_id', subject_list)

info = dict(dwi=[['subject_id', 'DWI_1000']],
            bvecs=[['subject_id', 'grad_1000']],
            bvals=[['subject_id', 'bval_1000']])

datasource = pe.Node(interface=nio.DataGrabber(infields=['subject_id'],
                                               outfields=info.keys()),
                     name='datasource')

datasource.inputs.sort_filelist = True
datasource.inputs.template = "%s/%s"
datasource.inputs.base_directory = data_dir
datasource.inputs.field_template = dict(
    dwi='%s/*%s.nii', bvecs='%s/*%s', bvals='%s/*%s')
datasource.inputs.template_args = info

inputnode = pe.Node(interface=util.IdentityInterface(
    fields=["dwi", "bvecs", "bvals", "struct"]), name="inputnode")

inputnode.inputs.struct = "/Users/erik/Dropbox/Analysis/TMSEEG/TMS007/orig_out_flirt.nii"

bet_T1 = pe.Node(interface=fsl.BET(), name = 'bet_T1')

bet_b0 = pe.Node(interface=fsl.BET(), name = 'bet_b0')
bet_b0.inputs.frac = 0.2
bet_b0.inputs.mask = True

dtifit = pe.Node(interface=fsl.DTIFit(), name = 'dtifit')
dtifit.inputs.save_tensor = True

affine_register = pe.Node(interface=fsl.FLIRT(), name = 'affine_register')
affine_register.inputs.cost = ('normmi')
affine_register.inputs.dof = 6

nonlinear_warp = pe.Node(interface=fsl.FNIRT(), name = 'nonlinear_warp')

register_tensor = pe.Node(interface=fsl.VecReg(), name = 'register_tensor')

decompose_tensor = pe.Node(interface=fsl.DecomposeTensor(), name = 'decompose_tensor')

fast_segment = pe.Node(interface=fsl.FAST(), name = 'fast_segment')
fast_segment.inputs.segments = True

conductivity_mapping = pe.Node(interface=dipy.EstimateConductivity(), name = 'conductivity_mapping')

# FractionalAnisotropyToGmsh_interface = util.Function(input_names=["in_file"],
#                             output_names=["out_file"],
#                             function=FAtoGmsh)
# fa_to_gmsh = pe.Node(interface=FractionalAnisotropyToGmsh_interface, name='fa_to_gmsh')

niftiTensorsToGmsh_interface = util.Function(input_names=["in_file", "mask_file", "fa_file"],
                             output_names=["out_file"],
                             function=nifti_tensors_to_gmsh)
tensors_to_gmsh_for_viz = pe.Node(interface=niftiTensorsToGmsh_interface, name='tensors_to_gmsh_for_viz')

niftiTensorsToGmshElementData_interface = util.Function(input_names=["in_file", "mask_file", "fa_file"],
                             output_names=["out_file"],
                             function=nifti_to_gmsh_tensor_elements)
tensors_to_gmsh_for_GetDP = pe.Node(interface=niftiTensorsToGmshElementData_interface, name='tensors_to_gmsh_for_GetDP')

erode_FA_map = pe.Node(interface=fsl.BET(), name="erode_FA_map")
erode_FA_map.inputs.mask = True

"""
Define the Workflow
-------------------
"""


preproc = pe.Workflow(name='preproc')
preproc.base_dir = op.abspath('preproc')

preproc.connect([(infosource, datasource, [("subject_id", "subject_id")])])

preproc.connect([(datasource, inputnode, [('dwi', 'dwi'),
               ('bvals', 'bvals'),
    ('bvecs', 'bvecs')
])])

scale_factor_workflow = create_get_scale_factor_workflow("scale_factor_workflow")

preproc.connect([(inputnode, dtifit, [('dwi', 'dwi'),
               ('bvals', 'bvals'),
               ('bvecs', 'bvecs'),
                ])])


preproc.connect([(inputnode, bet_b0,[('dwi','in_file')])])
preproc.connect([(bet_b0, dtifit,[('mask_file','mask')])])

preproc.connect([(inputnode, bet_T1,[('struct','in_file')])])
preproc.connect([(bet_T1, affine_register,[('out_file','reference')])])

preproc.connect([(dtifit, affine_register,[('FA','in_file')])])
preproc.connect([(affine_register, nonlinear_warp,[('out_file','in_file')])])
preproc.connect([(bet_T1, nonlinear_warp,[('out_file','ref_file')])])

preproc.connect([(dtifit, register_tensor,[('tensor','in_file')])])
preproc.connect([(bet_T1, register_tensor,[('out_file','ref_vol')])])
preproc.connect([(affine_register, register_tensor,[('out_matrix_file','affine_mat')])])
preproc.connect([(nonlinear_warp, register_tensor,[('field_file','warp_field')])])
preproc.connect([(register_tensor, decompose_tensor,[('out_file','in_file')])])

preproc.connect([(bet_T1, fast_segment,[("out_file","in_files")])])

preproc.connect([(fast_segment, scale_factor_workflow,[("tissue_class_files","inputnode.tissue_class_files")])])
preproc.connect([(decompose_tensor, scale_factor_workflow, [('L1', 'inputnode.L1'),
                                                            ('L2', 'inputnode.L2'),
                                                            ('L3', 'inputnode.L3'),
                ])])

preproc.connect([(scale_factor_workflow, conductivity_mapping,
    [("calc_scale_factor.scale_factor","eigenvalue_scaling_factor")])])

preproc.connect([(register_tensor, conductivity_mapping,
    [("out_file","in_file")])])

preproc.connect([(conductivity_mapping, tensors_to_gmsh_for_viz,[("out_file","in_file")])])
preproc.connect([(nonlinear_warp, tensors_to_gmsh_for_viz,[("warped_file","fa_file")])])
preproc.connect([(nonlinear_warp, erode_FA_map,[("warped_file","in_file")])])
preproc.connect([(erode_FA_map, tensors_to_gmsh_for_viz,[("mask_file","mask_file")])])

preproc.connect([(conductivity_mapping, tensors_to_gmsh_for_GetDP,[("out_file","in_file")])])
preproc.connect([(nonlinear_warp, tensors_to_gmsh_for_GetDP,[("warped_file","fa_file")])])
preproc.connect([(nonlinear_warp, tensors_to_gmsh_for_GetDP,[("warped_file","mask_file")])])


preproc.config['execution'] = {'remove_unnecessary_outputs': 'False',
                                   'hash_method': 'timestamp'}

if __name__ == '__main__':
    preproc.write_graph()
    preproc.run()

