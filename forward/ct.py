import nipype.interfaces.io as nio           # Data i/o
import nipype.interfaces.utility as util     # utility
import nipype.pipeline.engine as pe          # pypeline engine
import nipype.interfaces.freesurfer as fs
import nipype.interfaces.fsl as fsl
import nipype.interfaces.ants as ants
import os
import os.path as op
from nipype import logging
iflogger = logging.getLogger('interface')

def create_CT_MRI_coregistration_wf(name="coreg_wf"):
    workflow = pe.Workflow(name=name)
    inputnode = pe.Node(
        interface=util.IdentityInterface(fields=["T1", "CT"]), name="inputnode")
    outputnode = pe.Node(
        interface=util.IdentityInterface(fields=["conformed_T1", "conformed_CT"]), name="outputnode")


    conform_T1 = pe.Node(interface=fs.MRIConvert(), name="conform_T1")
    conform_T1.inputs.conform = True
    
    fix_CT_centroid = pe.Node(interface=fs.MRIConvert(), name="fix_CT_centroid")
    fix_CT_centroid.inputs.no_change = True
    fix_CT_centroid.inputs.conform = True

    workflow.connect(
        [(inputnode, conform_T1, [("T1", "in_file")])])
    workflow.connect(
        [(conform_T1, outputnode, [("out_file", "conformed_T1")])])

    workflow.connect(
        [(inputnode, fix_CT_centroid, [("CT", "in_file")])])
    workflow.connect(
        [(conform_T1, fix_CT_centroid, [("out_file", "reslice_like")])])
    workflow.connect(
        [(fix_CT_centroid, outputnode, [("out_file", "conformed_CT")])])

def create_CT_skin_seg_wf(name="skin_wf"):
    workflow = pe.Workflow(name=name)
    inputnode = pe.Node(
        interface=util.IdentityInterface(fields=["CT"]), name="inputnode")
    outputnode = pe.Node(
        interface=util.IdentityInterface(fields=["skin_mask"]), name="outputnode")

    init_high_thresh = 150
    skin_low = -705
    skin_high = -8
    
    step1_remove_high_val = pe.Node(interface=fsl.ImageMaths(), name="step1_remove_high_val")
    step1_remove_high_val.inputs.op_string = "-uthr %f -fmedian" % init_high_thresh

    step2_seg_skin = pe.Node(interface=fsl.ImageMaths(), name="step2_seg_skin")
    step2_seg_skin.inputs.op_string = "-thr %f -uthr %f" % (skin_low, skin_high)

    workflow.connect(
        [(inputnode, step1_remove_high_val, [("CT", "in_file")])])
    workflow.connect(
        [(step1_remove_high_val, step2_seg_skin, [("out_file", "in_file")])])
    workflow.connect(
        [(step2_seg_skin, outputnode, [("out_file", "skin_mask")])])
    return workflow



def preprocess_CT_wf(name="preproc_CT"):
    workflow = pe.Workflow(name=name)
    inputnode = pe.Node(
        interface=util.IdentityInterface(fields=["CT", "T1"]), name="inputnode")
    outputnode = pe.Node(
        interface=util.IdentityInterface(fields=["conformed_T1", "conformed_CT"]), name="outputnode")

    fixed_image = op.abspath("conform_MNI_MoutBija_20090622.nii")
    moving_image = op.abspath("conf_OrigCT.nii")

    conform_T1 = pe.Node(interface=fs.MRIConvert(), name="conform_T1")
    conform_T1.inputs.conform = True

    conform_CT = pe.Node(interface=fs.MRIConvert(), name="conform_CT")
    conform_CT.inputs.no_change = True
    conform_CT.inputs.conform = True

    register_CT = pe.Node(interface=ants.Registration(), name="register_CT")
    register_CT.inputs.fixed_image = [fixed_image]
    register_CT.inputs.moving_image = [moving_image]
    register_CT.inputs.output_transform_prefix = "Test2_output"
    register_CT.inputs.transforms = ['Translation', 'Rigid']
    register_CT.inputs.transform_parameters = [(0.1,), (0.1,)]
    register_CT.inputs.number_of_iterations = ([[10000, 111110, 11110]]*2)
    register_CT.inputs.dimension = 3
    register_CT.inputs.write_composite_transform = True
    register_CT.inputs.collapse_output_transforms = False
    register_CT.inputs.metric = ['Mattes'] * 2
    register_CT.inputs.metric_weight = [1] * 2
    register_CT.inputs.radius_or_number_of_bins = [32] * 2
    register_CT.inputs.sampling_strategy = ['Regular'] * 2
    register_CT.inputs.sampling_percentage = [0.3] * 2
    register_CT.inputs.convergence_threshold = [1.e-8] * 2
    register_CT.inputs.convergence_window_size = [20] * 2
    register_CT.inputs.smoothing_sigmas = [[4, 2, 1]] * 2
    register_CT.inputs.sigma_units = ['vox'] * 2
    register_CT.inputs.shrink_factors = [[6, 4, 2]] + [[3, 2, 1]]
    register_CT.inputs.use_estimate_learning_rate_once = [True] * 2
    register_CT.inputs.use_histogram_matching = [False] * 2
    register_CT.inputs.initial_moving_transform_com = True
    register_CT.inputs.output_warped_image = 'conf_CT_fixed.nii'

    workflow.connect(
        [(inputnode, conform_CT, [("CT", "in_file")])])
    workflow.connect(
        [(conform_CT, register_CT, [("out_file", "moving_image")])])

    workflow.connect(
        [(inputnode, conform_T1, [("T1", "in_file")])])
    workflow.connect(
        [(conform_T1, outputnode, [("out_file", "conformed_T1")])])
    workflow.connect(
        [(inputnode, register_CT, [("T1", "fixed_image")])])
    workflow.connect(
        [(register_CT, outputnode, [("out_file", "skin_mask")])])
    return workflow


def create_CT_seg_wf(name="CT_seg_wf"):
    workflow = pe.Workflow(name=name)
    inputnode = pe.Node(
        interface=util.IdentityInterface(fields=["CT"]), name="inputnode")
    outputnode = pe.Node(
        interface=util.IdentityInterface(fields=["skin_mask", "skull_mask",
            "spongiform_mask", "instrument_mask", "implant_mask", "volumes"]), name="outputnode")

    skull_low = 200
    skull_high = 2000

    instrument_low = 1700

    implant_low = 155
    implant_high = 195

    sponge_low = 500
    sponge_high = 700
    
    skull_seg = pe.Node(interface=fsl.ImageMaths(), name="skull_seg")
    skull_seg.inputs.op_string = "-thr %f -uthr %f" % (skull_low, skull_high)
    skull_seg.inputs.out_file = "skull.nii.gz"

    spongiform_seg = pe.Node(interface=fsl.ImageMaths(), name="spongiform_seg")
    spongiform_seg.inputs.op_string = "-thr %f -uthr %f -fmedian" % (sponge_low, sponge_high)
    spongiform_seg.inputs.out_file = "spongiform_bone.nii.gz"

    instrument_seg = pe.Node(interface=fsl.ImageMaths(), name="instrument_seg")
    instrument_seg.inputs.op_string = "-thr %f" % instrument_low
    instrument_seg.inputs.out_file = "instruments.nii.gz"

    implant_seg = pe.Node(interface=fsl.ImageMaths(), name="implant_seg")
    implant_seg.inputs.op_string = "-thr %f -uthr %f -fmedian -fmedian" % (implant_low, implant_high)
    implant_seg.inputs.out_file = "implant.nii.gz"

    init_high_thresh = 150
    skin_low = -705
    skin_high = -8
    
    tissue_low = -1000

    create_head_mask = pe.Node(interface=fsl.ImageMaths(), name="create_head_mask")
    create_head_mask.inputs.op_string = "-thr %f -bin -fillh26 -dilM -dilM -dilM -dilM" % tissue_low

    step1_remove_high_val = pe.Node(interface=fsl.ImageMaths(), name="step1_remove_high_val")
    step1_remove_high_val.inputs.op_string = "-uthr %f -fmedian" % init_high_thresh

    step2_seg_skin = pe.Node(interface=fsl.ImageMaths(), name="step2_seg_skin")
    step2_seg_skin.inputs.op_string = "-thr %f -uthr %f" % (skin_low, skin_high)

    create_skin_mask = pe.Node(interface=fsl.MultiImageMaths(), name="create_skin_mask")
    create_skin_mask.inputs.op_string = "-mul %s"
    create_skin_mask.inputs.out_file = "skin.nii.gz"

    
    merge_outputs = pe.Node(
        interface=util.Merge(5), name='merge_for_dilate_subcortical_mask')

    workflow.connect(
        [(inputnode, skull_seg, [("CT", "in_file")])])
    workflow.connect(
        [(inputnode, spongiform_seg, [("CT", "in_file")])])
    workflow.connect(
        [(inputnode, implant_seg, [("CT", "in_file")])])
    workflow.connect(
        [(inputnode, instrument_seg, [("CT", "in_file")])])
    workflow.connect(
        [(inputnode, step1_remove_high_val, [("CT", "in_file")])])
    workflow.connect(
        [(inputnode, create_head_mask, [("CT", "in_file")])])
    workflow.connect(
        [(step1_remove_high_val, step2_seg_skin, [("out_file", "in_file")])])
    workflow.connect(
        [(step2_seg_skin, create_skin_mask, [("out_file", "in_file")])])
    workflow.connect(
        [(create_head_mask, create_skin_mask, [("out_file", "operand_files")])])
    workflow.connect(
        [(create_skin_mask, outputnode, [("out_file", "skin_mask")])])
    workflow.connect(
        [(skull_seg, outputnode, [("out_file", "skull_mask")])])
    workflow.connect(
        [(spongiform_seg, outputnode, [("out_file", "spongiform_mask")])])
    workflow.connect(
        [(instrument_seg, outputnode, [("out_file", "instrument_mask")])])
    workflow.connect(
        [(implant_seg, outputnode, [("out_file", "implant_mask")])])
    workflow.connect(
        [(create_skin_mask, merge_outputs, [("out_file", "in1")])])
    workflow.connect(
        [(skull_seg, merge_outputs, [("out_file", "in2")])])
    workflow.connect(
        [(instrument_seg, merge_outputs, [("out_file", "in3")])])
    workflow.connect(
        [(implant_seg, merge_outputs, [("out_file", "in4")])])
    workflow.connect(
        [(create_head_mask, merge_outputs, [("out_file", "in5")])])
    workflow.connect(
        [(merge_outputs, outputnode, [("out", "volumes")])])

    return workflow
