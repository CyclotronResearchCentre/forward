import nipype.interfaces.io as nio           # Data i/o
import nipype.interfaces.utility as util     # utility
import nipype.pipeline.engine as pe          # pypeline engine
import nipype.interfaces.fsl as fsl
import nipype.interfaces.freesurfer as fs
import nipype.interfaces.meshfix as mf
from nipype.interfaces.gmsh import Gmsh
import os
import os.path as op                     # system functions
import ipdb
from nipype.interfaces.meshfix import (iter_push_out_fn,
                                       iter_decoupling_fn, iter_remove_throats_fn, remove_spikes_fn,
                                       decouple_surfaces_fn, mask_from_labels_fn,
                                       decouple_ventricles_fn, decouple_input_from_GM_fn,
                                       decouple_outout_cutin_fn)

fsl.FSLCommand.set_default_output_type('NIFTI')


def volume_mesh_fn(subject_id, gm, wm, csf, skull, skin, cerebellum, ventricles):
    import os.path as op
    out_file = op.abspath(subject_id + ".geo")
    f = open(out_file, 'w')
    f.write('Mesh.Algorithm3D=4; //1=delaunay (tetgen) and 4=frontal (netgen)\n')
    f.write('Mesh.Optimize=1;\n')
    f.write('Mesh.OptimizeNetgen=1;\n')
    f.write('Merge "%s"; // first surface\n' % wm)
    f.write('Merge "%s";\n' % gm)
    f.write('Merge "%s";\n' % csf)
    f.write('Merge "%s";\n' % skull)
    f.write('Merge "%s";\n' % skin)
    f.write('Merge "%s";\n' % cerebellum)
    f.write('Merge "%s";\n' % ventricles)
    f.write('Surface Loop(1) = {1}; // surface number on rhs; 1: wm.stl\n')
    f.write('Surface Loop(2) = {2}; // gm.stl\n')
    f.write('Surface Loop(3) = {3}; // csf.stl\n')
    f.write('Surface Loop(4) = {4}; // skull.stl\n')
    f.write('Surface Loop(5) = {5}; // skin.stl\n')
    f.write('Surface Loop(6) = {6}; // cerebellum.stl\n')
    f.write('Surface Loop(7) = {7}; // ventricles.stl\n')
    f.write('Volume(1) = {1, 7};    // surface loop numbers on rhs; e.g., volume between surface loops 1 (wm) and 7 (ventricles) is sm \n')
    f.write('Volume(2) = {1, 2};    // gm volume\n')
    f.write('Volume(3) = {2, 6, 3}; // csf volume: outside gm, cerebellum; inside csf\n')
    f.write('Volume(4) = {3, 4};    // skull volume\n')
    f.write('Volume(5) = {4, 5};    // skin volume\n')
    f.write('Volume(6) = {6};       // cerebellum volume\n')
    f.write('Volume(7) = {7};       // ventricle volume\n')
    f.write('Physical Surface(1) = {1,6}; // merge wm and cerebellum; LHS: target surface region number, RHS: surface number (i.e. from merge ...)\n')
    f.write('Physical Surface(2) = {2}; \n')
    f.write('Physical Surface(3) = {3,7}; // merge csf and ventricles \n')
    f.write('Physical Surface(4) = {4};\n')
    f.write('Physical Surface(5) = {5};\n')
    f.write('Physical Volume(1) = {1,6}; // merge wm and cerebellum; LHS: target volume region number, RHS: volume number\n')
    f.write('Physical Volume(2) = {2}; \n')
    f.write('Physical Volume(3) = {3,7}; // merge csf and ventricles\n')
    f.write('Physical Volume(4) = {4};\n')
    f.write('Physical Volume(5) = {5};\n')
    f.close()
    return out_file


def create_prepare_surface_wf(name, number_of_vertices):
    inputnode = pe.Node(interface=util.IdentityInterface(
        fields=["in_file", "optic_rad"]), name="inputnode")
    outputnode = pe.Node(
        interface=util.IdentityInterface(fields=["out_file"]), name="outputnode")

    initial_clean = pe.Node(interface=mf.MeshFix(), name='initial_clean')
    initial_clean.inputs.quiet_mode = True
    initial_clean.inputs.epsilon_angle = 2
    initial_clean.inputs.remove_handles = True

    restrict_vertices = initial_clean.clone("restrict_vertices")
    restrict_vertices.inputs.uniform_remeshing_steps = 5
    restrict_vertices.inputs.uniform_remeshing_vertices = number_of_vertices

    second_smooth = initial_clean.clone("second_smooth")
    second_smooth.inputs.uniform_remeshing_steps = 1

    finalize_surf = initial_clean.clone("finalize_surf")

    cut_away_optic_rad = initial_clean.clone("cut_away_optic_rad")
    cut_away_optic_rad.inputs.number_of_biggest_shells = 2
    cut_away_optic_rad.inputs.cut_inner = 0

    hemi_surf_cleanup = pe.Workflow(name=name)
    hemi_surf_cleanup.connect(
        [(inputnode, initial_clean, [('in_file', 'in_file1')])])
    hemi_surf_cleanup.connect(
        [(initial_clean, restrict_vertices, [('mesh_file', 'in_file1')])])
    hemi_surf_cleanup.connect(
        [(restrict_vertices, second_smooth, [('mesh_file', 'in_file1')])])
    hemi_surf_cleanup.connect(
        [(second_smooth, finalize_surf, [('mesh_file', 'in_file1')])])
    hemi_surf_cleanup.connect(
        [(finalize_surf, cut_away_optic_rad, [('mesh_file', 'in_file1')])])
    hemi_surf_cleanup.connect(
        [(inputnode, cut_away_optic_rad, [('optic_rad', 'in_file2')])])
    hemi_surf_cleanup.connect(
        [(cut_away_optic_rad, outputnode, [('mesh_file', 'out_file')])])
    return hemi_surf_cleanup


def create_iterative_cleanup_wf(name):
    inputnode = pe.Node(
        interface=util.IdentityInterface(fields=["in_file"]), name="inputnode")
    outputnode = pe.Node(
        interface=util.IdentityInterface(fields=["out_file"]), name="outputnode")

    clean1 = pe.Node(interface=mf.MeshFix(), name='clean1')
    clean1.inputs.quiet_mode = True
    clean1.inputs.epsilon_angle = 2
    clean1.inputs.uniform_remeshing_steps = 1

    clean2 = clean1.clone("clean2")
    clean3 = clean1.clone("clean3")
    clean4 = clean2.clone("clean4")

    clean_wf = pe.Workflow(name=name + '_clean_wf')
    clean_wf.connect([(inputnode, clean1, [('in_file', 'in_file1')])])
    clean_wf.connect([(clean1, clean2, [('mesh_file', 'in_file1')])])
    clean_wf.connect([(clean2, clean3, [('mesh_file', 'in_file1')])])
    clean_wf.connect([(clean3, clean4, [('mesh_file', 'in_file1')])])
    clean_wf.connect([(clean4, outputnode, [('mesh_file', 'out_file')])])
    return clean_wf

    

def create_tess_constrain_smooth_wf(name):
    inputnode = pe.Node(
        interface=util.IdentityInterface(fields=["in_file"]), name="inputnode")
    outputnode = pe.Node(
        interface=util.IdentityInterface(fields=["out_file"]), name="outputnode")

    workflow = pe.Workflow(name=name + '_workflow')

    mask_to_fsspace = pe.Node(
            interface=fsl.ApplyXfm(), name='mask_to_fsspace')
    mask_to_fsspace.inputs.in_matrix_file = "/Users/erik/Dropbox/Analysis/TMSEEG/fs2fsl_conform_inverse.mat"

    mask_to_mgz = pe.Node(interface=fs.MRIConvert(), name='mask_to_mgz')
    mask_to_mgz.inputs.out_type = "mgz"
    mask_to_mgz.inputs.conform = True

    tessellate = pe.Node(
        interface=fs.MRIMarchingCubes(), name='tessellate')
    tessellate.inputs.label_value = 255
    tessellate.inputs.out_file = "csf_FS.stl"

    smooth1 = pe.Node(interface=mf.MeshFix(), name='smooth1')
    smooth1.inputs.laplacian_smoothing_steps = 5
    smooth1.inputs.save_as_stl = True

    n_vertices = 60000
    constrain_vertices = pe.Node(
        interface=mf.MeshFix(), name='constrain_vertices')
    constrain_vertices.inputs.epsilon_angle = 2.0
    constrain_vertices.inputs.uniform_remeshing_steps = 5
    constrain_vertices.inputs.uniform_remeshing_vertices = n_vertices

    smooth2 = pe.Node(interface=mf.MeshFix(), name='smooth2')
    smooth2.inputs.epsilon_angle = 2.0
    smooth2.inputs.laplacian_smoothing_steps = 1

    workflow.connect([(inputnode, mask_to_fsspace, [("in_file", "in_file")])])
    workflow.connect([(inputnode, mask_to_fsspace, [("in_file", "reference")])])
    workflow.connect([(mask_to_fsspace, mask_to_mgz, [("out_file", "in_file")])])

    workflow.connect([(mask_to_mgz, tessellate, [("out_file", "in_file")])])
    workflow.connect([(tessellate, smooth1, [("surface", "in_file1")])])

    workflow.connect([(smooth1, constrain_vertices, [("mesh_file", "in_file1")])])
    workflow.connect([(constrain_vertices, smooth2, [("mesh_file", "in_file1")])])
    workflow.connect([(smooth2, outputnode, [("mesh_file", "out_file")])])
    return workflow

def create_stl_fsmesh_floodfilled_wf(name):
    inputnode = pe.Node(
        interface=util.IdentityInterface(fields=["in_file", "t1_fsl_space"]), name="inputnode")
    outputnode = pe.Node(
        interface=util.IdentityInterface(fields=["stl_mesh", "filled_fsl_space"]), name="outputnode")

    workflow = pe.Workflow(name=name + '_workflow')
    input_to_STL = pe.Node(interface=mf.MeshFix(), name='input_to_STL')
    input_to_STL.inputs.save_as_stl = True
    input_to_STL.inputs.out_filename = "input.stl"

    input_to_fsmesh = pe.Node(interface=mf.MeshFix(), name='input_to_fsmesh')
    input_to_fsmesh.inputs.dont_clean = True
    input_to_fsmesh.inputs.save_as_freesurfer_mesh = True

    unityxform = pe.Node(
        interface=fs.TransformSurface2Talairach(), name='unityxform')
    unityxform.inputs.transform_file = "/Users/erik/Dropbox/Analysis/TMSEEG/unity.xfm"

    floodfill = pe.Node(interface=fs.FloodfillSurface(), name='floodfill')
    floodfill.inputs.conform_before_writing = True

    filled_to_fsl_space = pe.Node(interface=fsl.ApplyXfm(), name='filled_to_fsl_space')
    filled_to_fsl_space.inputs.in_matrix_file = "/Users/erik/Dropbox/Analysis/TMSEEG/fs2fsl_conform.mat"

    workflow.connect([(inputnode, input_to_STL, [("in_file", "in_file1")])])
    workflow.connect([(input_to_STL, outputnode, [("mesh_file", "stl_mesh")])])
    workflow.connect([(input_to_STL, input_to_fsmesh, [("mesh_file", "in_file1")])])
    workflow.connect([(input_to_fsmesh, unityxform, [("mesh_file", "in_file")])])
    workflow.connect([(inputnode, unityxform, [("t1_fsl_space", "destination_volume")])])
    workflow.connect([(inputnode, unityxform, [("t1_fsl_space", "source_volume")])])
    workflow.connect([(unityxform, floodfill, [("surface", "in_file")])])
    workflow.connect(
        [(floodfill, filled_to_fsl_space, [("out_file", "in_file")])])
    workflow.connect([(inputnode, filled_to_fsl_space, [("t1_fsl_space", "reference")])])
    workflow.connect([(filled_to_fsl_space, outputnode, [("out_file", "filled_fsl_space")])])
    return workflow

def brain_workflow(name):
    inputnode = pe.Node(
        interface=util.IdentityInterface(fields=["subject_id", "subjects_dir"]), name="inputnode")
    outputnode = pe.Node(
        interface=util.IdentityInterface(fields=["gray_matter_volume", "white_matter_volume",
                                                 "white_matter_surface", "gray_matter_surface",
                                                 "t1_fsl_space", "nu", "labels_fsl_space",
                                                 "Conform2MNI", "MNI2Conform"]), name="outputnode")

    workflow = pe.Workflow(name=name + "_workflow")

    FreeSurferSource = pe.Node(
        interface=nio.FreeSurferSource(), name='fssource')

    FreeSurferSourceLH = pe.Node(
        interface=nio.FreeSurferSource(), name='fssourceLH')
    FreeSurferSourceLH.inputs.hemi = 'lh'
    FreeSurferSourceRH = pe.Node(
        interface=nio.FreeSurferSource(), name='fssourceRH')
    FreeSurferSourceRH.inputs.hemi = 'rh'

    workflow.connect(
        [(inputnode, FreeSurferSource, [("subject_id", "subject_id")])])
    workflow.connect(
        [(inputnode, FreeSurferSourceLH, [("subject_id", "subject_id")])])
    workflow.connect(
        [(inputnode, FreeSurferSourceRH, [("subject_id", "subject_id")])])

    workflow.connect(
        [(inputnode, FreeSurferSource, [("subjects_dir", "subjects_dir")])])
    workflow.connect(
        [(inputnode, FreeSurferSourceLH, [("subjects_dir", "subjects_dir")])])
    workflow.connect(
        [(inputnode, FreeSurferSourceRH, [("subjects_dir", "subjects_dir")])])

    # Convert surface meshes to STL
    lh_wm_to_STL = pe.Node(interface=fs.MRIsConvert(), name='lh_wm_to_STL')
    lh_wm_to_STL.inputs.out_datatype = "stl"
    rh_wm_to_STL = lh_wm_to_STL.clone('rh_wm_to_STL')
    lh_gm_to_STL = lh_wm_to_STL.clone('lh_gm_to_STL')
    rh_gm_to_STL = lh_wm_to_STL.clone('rh_gm_to_STL')

    workflow.connect(
        [(FreeSurferSourceLH, lh_wm_to_STL, [("white", "in_file")])])
    workflow.connect(
        [(FreeSurferSourceLH, lh_gm_to_STL, [("pial", "in_file")])])

    workflow.connect(
        [(FreeSurferSourceRH, rh_wm_to_STL, [("white", "in_file")])])
    workflow.connect(
        [(FreeSurferSourceRH, rh_gm_to_STL, [("pial", "in_file")])])

    # Transform T1 from conformed Freesurfer space to FSL
    orig_mgz_to_nii = pe.Node(
        interface=fs.MRIConvert(), name='orig_mgz_to_nii')
    nu_mgz_to_nii = orig_mgz_to_nii.clone("nu_mgz_to_nii")
    t1_to_fsl_space = pe.Node(interface=fsl.ApplyXfm(), name='t1_to_fsl_space')
    t1_to_fsl_space.inputs.in_matrix_file = "/Users/erik/Dropbox/Analysis/TMSEEG/fs2fsl_conform.mat"

    workflow.connect(
        [(FreeSurferSource, orig_mgz_to_nii, [("orig", "in_file")])])
    workflow.connect(
        [(FreeSurferSource, nu_mgz_to_nii, [("nu", "in_file")])])
    workflow.connect(
        [(nu_mgz_to_nii, outputnode, [("out_file", "nu")])])    
    workflow.connect(
        [(orig_mgz_to_nii, t1_to_fsl_space, [("out_file", "in_file")])])
    workflow.connect(
        [(orig_mgz_to_nii, t1_to_fsl_space, [("out_file", "reference")])])

    # Mask T1 and get transformations conformed space <-> MNI space
    mask_T1 = pe.Node(interface=fsl.ImageMaths(), name="mask_T1")
    mask_T1.inputs.op_string = "-bin"

    get_conform_to_MNI = pe.Node(interface=fsl.FLIRT(), name='get_conform_to_MNI')
    get_conform_to_MNI.inputs.cost = "mutualinfo"
    get_conform_to_MNI.inputs.reference = op.join(os.environ['FSL_DIR'],"data","standard", "MNI152_T1_1mm_brain.nii.gz")

    get_MNI_to_conform = pe.Node(
        interface=fsl.ConvertXFM(), name='get_MNI_to_conform')
    get_MNI_to_conform.inputs.invert_xfm = True

    workflow.connect([(t1_to_fsl_space, mask_T1, [("out_file", "in_file")])])
    workflow.connect([(t1_to_fsl_space, get_conform_to_MNI, [("out_file", "in_file")])])
    workflow.connect([(get_conform_to_MNI, get_MNI_to_conform, [("out_matrix_file", "in_file")])])

    # Label mask from conformed Freesurfer space to FSL
    labels_mgz_to_nii = orig_mgz_to_nii.clone("labels_mgz_to_nii")
    labels_to_fsl_space = t1_to_fsl_space.clone('labels_to_fsl_space')

    workflow.connect(
        [(FreeSurferSource, labels_mgz_to_nii, [("aseg", "in_file")])])
    workflow.connect(
        [(labels_mgz_to_nii, labels_to_fsl_space, [("out_file", "in_file")])])
    workflow.connect(
        [(t1_to_fsl_space, labels_to_fsl_space, [("out_file", "reference")])])

    # Create corpus callosum mask
    isolate_CC_labels = pe.Node(
        interface=fsl.ImageMaths(), name="isolate_CC_labels")
    isolate_CC_labels.inputs.op_string = "-thr 251 -uthr 255 -bin"
    apply_kernel_to_CC = pe.Node(
        interface=fsl.MultiImageMaths(), name="apply_kernel_to_CC")
    apply_kernel_to_CC.inputs.op_string = "-kernel file %s -fmean"
    apply_kernel_to_CC.output_datatype = "float"
    apply_kernel_to_CC.inputs.operand_files = "/Users/erik/Dropbox/Analysis/TMSEEG/CC_filterkernel.nii.gz"

    finalize_callosal_mask = pe.Node(
        interface=fsl.ImageMaths(), name="finalize_callosal_mask")
    finalize_callosal_mask.inputs.op_string = "-thr 0.1 -bin"

    workflow.connect(
        [(labels_to_fsl_space, isolate_CC_labels, [("out_file", "in_file")])])
    workflow.connect(
        [(isolate_CC_labels, apply_kernel_to_CC, [("out_file", "in_file")])])
    workflow.connect(
        [(apply_kernel_to_CC, finalize_callosal_mask, [("out_file", "in_file")])])

    # Create mask with shifted central CSF slab (to mask the fornix)
    isolate_CSF_slab = pe.Node(
        interface=fsl.ImageMaths(), name="isolate_CSF_slab")
    isolate_CSF_slab.inputs.op_string = "-thr 24 -uthr 24 -dilM -bin"
    CSF_yshift = pe.Node(interface=fsl.ApplyXfm(), name='CSF_yshift')
    CSF_yshift.inputs.in_matrix_file = "/Users/erik/Dropbox/Analysis/TMSEEG/yshift.mat"

    workflow.connect(
        [(labels_to_fsl_space, isolate_CSF_slab, [("out_file", "in_file")])])
    workflow.connect(
        [(isolate_CSF_slab, CSF_yshift, [("out_file", "in_file")])])
    workflow.connect(
        [(isolate_CSF_slab, CSF_yshift, [("out_file", "reference")])])

    # Create large subcortical mask (for joining the hemispheres)
    create_mask_from_labels_interface = util.Function(
        input_names=["in_file", "label_values"],
        output_names=["out_file"],
        function=mask_from_labels_fn)
    create_subcortical_mask = pe.Node(
        interface=create_mask_from_labels_interface, name="create_subcortical_mask")
    create_subcortical_mask.inputs.label_values = [
        4, 10, 14, 24, 28, 31, 43, 49, 60, 63]
    dilate_subcortical_mask = pe.Node(
        interface=fsl.MultiImageMaths(), name="dilate_subcortical_mask")
    dilate_subcortical_mask.inputs.op_string = "-add %s -add %s -dilM -ero -bin"

    subcort_mask_to_fsspace = pe.Node(
        interface=fsl.ApplyXfm(), name='subcort_mask_to_fsspace')
    subcort_mask_to_fsspace.inputs.in_matrix_file = "/Users/erik/Dropbox/Analysis/TMSEEG/fs2fsl_conform_inverse.mat"
    subcort_to_mgz = pe.Node(interface=fs.MRIConvert(), name='subcort_to_mgz')
    subcort_to_mgz.inputs.out_type = "mgz"
    subcort_to_mgz.inputs.conform = True
    merge_for_dilate_subcortical_mask = pe.Node(
        interface=util.Merge(2), name='merge_for_dilate_subcortical_mask')

    workflow.connect(
        [(labels_to_fsl_space, create_subcortical_mask, [("out_file", "in_file")])])
    workflow.connect(
        [(create_subcortical_mask, dilate_subcortical_mask, [("out_file", "in_file")])])

    workflow.connect(
        [(finalize_callosal_mask, merge_for_dilate_subcortical_mask, [("out_file", "in1")])])
    workflow.connect(
        [(CSF_yshift, merge_for_dilate_subcortical_mask, [("out_file", "in2")])])
    workflow.connect(
        [(merge_for_dilate_subcortical_mask, dilate_subcortical_mask, [("out", "operand_files")])])

    workflow.connect(
        [(dilate_subcortical_mask, subcort_mask_to_fsspace, [("out_file", "in_file")])])
    workflow.connect(
        [(dilate_subcortical_mask, subcort_mask_to_fsspace, [("out_file", "reference")])])
    workflow.connect(
        [(subcort_mask_to_fsspace, subcort_to_mgz, [("out_file", "in_file")])])

    # Create subcortical surface, smooth and convert to STL
    tessellate_subcort = pe.Node(
        interface=fs.MRIMarchingCubes(), name='tessellate_subcort')
    tessellate_subcort.inputs.label_value = 255
    tessellate_subcort.inputs.out_file = "subcortical_FS.stl"
    smooth_subcort = pe.Node(interface=mf.MeshFix(), name='smooth_subcort')
    smooth_subcort.inputs.laplacian_smoothing_steps = 5
    smooth_subcort.inputs.save_as_stl = True

    workflow.connect(
        [(subcort_to_mgz, tessellate_subcort, [("out_file", "in_file")])])
    workflow.connect(
        [(tessellate_subcort, smooth_subcort, [("surface", "in_file1")])])

    # Create optic radiation mask and smoothed surface
    mask_optic_rad = pe.Node(interface=fsl.ImageMaths(), name="mask_optic_rad")
    mask_optic_rad.inputs.op_string = "-thr 85 -uthr 85 -bin"
    optic_rad_mask_to_fsspace = subcort_mask_to_fsspace.clone(
        "optic_rad_mask_to_fsspace")
    optic_rad_to_mgz = subcort_to_mgz.clone("optic_rad_to_mgz")
    tesselate_optic_rad = tessellate_subcort.clone("tesselate_optic_rad")
    smooth_optic_rad = smooth_subcort.clone("smooth_optic_rad")

    workflow.connect(
        [(labels_to_fsl_space, mask_optic_rad, [("out_file", "in_file")])])
    workflow.connect(
        [(mask_optic_rad, optic_rad_mask_to_fsspace, [("out_file", "in_file")])])
    workflow.connect(
        [(mask_optic_rad, optic_rad_mask_to_fsspace, [("out_file", "reference")])])
    workflow.connect(
        [(optic_rad_mask_to_fsspace, optic_rad_to_mgz, [("out_file", "in_file")])])
    workflow.connect(
        [(optic_rad_to_mgz, tesselate_optic_rad, [("out_file", "in_file")])])
    workflow.connect(
        [(tesselate_optic_rad, smooth_optic_rad, [("surface", "in_file1")])])

    # Mesh fixing tasks for subcortical masks
    clean_subcort = pe.Node(interface=mf.MeshFix(), name='clean_subcort')
    clean_subcort.inputs.epsilon_angle = 2.0
    clean_subcort.inputs.remove_handles = True
    clean_subcort.inputs.quiet_mode = True

    workflow.connect(
        [(smooth_subcort, clean_subcort, [("mesh_file", "in_file1")])])

    n_vertices = 60000
    constrain_subcort_vertices = pe.Node(
        interface=mf.MeshFix(), name='constrain_subcort_vertices')
    constrain_subcort_vertices.inputs.epsilon_angle = 2.0
    constrain_subcort_vertices.inputs.uniform_remeshing_steps = 5
    constrain_subcort_vertices.inputs.uniform_remeshing_vertices = n_vertices / 8
    constrain_subcort_vertices.inputs.quiet_mode = True

    workflow.connect(
        [(clean_subcort, constrain_subcort_vertices, [("mesh_file", "in_file1")])])

    # Create large, small, and very small subcortical surfaces
    subcort_small = pe.Node(interface=mf.MeshFix(), name='subcort_small')
    subcort_small.inputs.dilation = -1
    subcort_small.inputs.epsilon_angle = 2.0
    subcort_small.inputs.quiet_mode = True

    subcort_verysmall = subcort_small.clone("subcort_verysmall")
    subcort_verysmall.inputs.dilation = -2

    subcort_large = subcort_small.clone("subcort_large")
    subcort_large.inputs.dilation = 1

    workflow.connect(
        [(constrain_subcort_vertices, subcort_small, [("mesh_file", "in_file1")])])
    workflow.connect(
        [(constrain_subcort_vertices, subcort_verysmall, [("mesh_file", "in_file1")])])
    workflow.connect(
        [(constrain_subcort_vertices, subcort_large, [("mesh_file", "in_file1")])])

    smooth_subcort_small = pe.Node(
        interface=mf.MeshFix(), name='smooth_subcort_small')
    smooth_subcort_small.inputs.epsilon_angle = 2.0
    smooth_subcort_small.inputs.uniform_remeshing_steps = 5
    smooth_subcort_small.inputs.quiet_mode = True

    workflow.connect(
        [(subcort_small, smooth_subcort_small, [("mesh_file", "in_file1")])])

    # Clean optic radiation surface
    clean_optic_rad = clean_subcort.clone("clean_optic_rad")
    smooth_optic_rad2 = smooth_subcort_small.clone("smooth_optic_rad2")
    optic_rad_verysmall = subcort_small.clone("optic_rad_verysmall")
    optic_rad_verysmall.inputs.dilation = -2

    workflow.connect(
        [(smooth_optic_rad, clean_optic_rad, [("mesh_file", "in_file1")])])
    workflow.connect(
        [(clean_optic_rad, smooth_optic_rad2, [("mesh_file", "in_file1")])])
    workflow.connect(
        [(smooth_optic_rad2, optic_rad_verysmall, [("mesh_file", "in_file1")])])

    # Prepare white matter surfaces
    prepare_wm_lh_surf = create_prepare_surface_wf(
        "prepare_wm_lh_surf", n_vertices)
    prepare_wm_rh_surf = create_prepare_surface_wf(
        "prepare_wm_rh_surf", n_vertices)

    workflow.connect(
        [(lh_wm_to_STL, prepare_wm_lh_surf, [("converted", "inputnode.in_file")])])
    workflow.connect(
        [(rh_wm_to_STL, prepare_wm_rh_surf, [("converted", "inputnode.in_file")])])
    workflow.connect(
        [(optic_rad_verysmall, prepare_wm_lh_surf, [("mesh_file", "inputnode.optic_rad")])])
    workflow.connect(
        [(optic_rad_verysmall, prepare_wm_rh_surf, [("mesh_file", "inputnode.optic_rad")])])

    # Prepare gray matter surfaces
    prepare_gm_lh_surf = create_prepare_surface_wf(
        "prepare_gm_lh_surf", 115 * n_vertices / 100)
    prepare_gm_rh_surf = create_prepare_surface_wf(
        "prepare_gm_rh_surf", 115 * n_vertices / 100)

    workflow.connect(
        [(lh_gm_to_STL, prepare_gm_lh_surf, [("converted", "inputnode.in_file")])])
    workflow.connect(
        [(rh_gm_to_STL, prepare_gm_rh_surf, [("converted", "inputnode.in_file")])])
    workflow.connect(
        [(optic_rad_verysmall, prepare_gm_lh_surf, [("mesh_file", "inputnode.optic_rad")])])
    workflow.connect(
        [(optic_rad_verysmall, prepare_gm_rh_surf, [("mesh_file", "inputnode.optic_rad")])])

    # Join the white matter hemispheres
    join_wm_with_subcort_small = pe.Node(
        interface=mf.MeshFix(), name='join_wm_with_subcort_small')
    join_wm_with_subcort_small.inputs.quiet_mode = True
    join_wm_with_subcort_small.inputs.epsilon_angle = 2
    join_wm_with_subcort_small.inputs.remove_handles = True
    join_wm_with_subcort_small.inputs.join_overlapping_largest_components = True
    join_wm_with_subcort_small.inputs.number_of_biggest_shells = 2

    join_wm_hemis = join_wm_with_subcort_small.clone("join_wm_hemis")

    workflow.connect([(prepare_wm_lh_surf, join_wm_with_subcort_small,
                       [("outputnode.out_file", "in_file1")])])
    workflow.connect(
        [(subcort_small, join_wm_with_subcort_small, [("mesh_file", "in_file2")])])

    workflow.connect(
        [(join_wm_with_subcort_small, join_wm_hemis, [("mesh_file", "in_file1")])])
    workflow.connect([(prepare_wm_rh_surf, join_wm_hemis,
                       [("outputnode.out_file", "in_file2")])])

    remove_throats_interface = util.Function(
        input_names=["outer_mesh", "inner_mesh", "iterations"],
        output_names=["outer_mesh"],
        function=iter_remove_throats_fn)

    iter_push_out_interface = util.Function(
        input_names=["outer_mesh", "inner_mesh", "iterations"],
        output_names=["outer_mesh"],
        function=iter_push_out_fn)

    iter_remove_throats_wm = pe.Node(
        interface=remove_throats_interface, name="iter_remove_throats_wm")
    iter_push_out_wm = pe.Node(
        interface=iter_push_out_interface, name="iter_push_out_wm")

    workflow.connect(
        [(join_wm_hemis, iter_remove_throats_wm, [("mesh_file", "outer_mesh")])])
    workflow.connect(
        [(subcort_small, iter_remove_throats_wm, [("mesh_file", "inner_mesh")])])
    workflow.connect(
        [(iter_remove_throats_wm, iter_push_out_wm, [("outer_mesh", "outer_mesh")])])
    workflow.connect(
        [(constrain_subcort_vertices, iter_push_out_wm, [("mesh_file", "inner_mesh")])])

    # Join the gray matter hemispheres

    join_gm_with_subcort_small = join_wm_hemis.clone(
        'join_gm_with_subcort_small')
    join_gm_hemis = join_wm_hemis.clone('join_gm_hemis')

    workflow.connect([(prepare_gm_lh_surf, join_gm_with_subcort_small,
                       [("outputnode.out_file", "in_file1")])])
    workflow.connect(
        [(subcort_small, join_gm_with_subcort_small, [("mesh_file", "in_file2")])])

    workflow.connect(
        [(join_gm_with_subcort_small, join_gm_hemis, [("mesh_file", "in_file1")])])
    workflow.connect([(prepare_gm_rh_surf, join_gm_hemis,
                       [("outputnode.out_file", "in_file2")])])

    iter_remove_throats_gm = pe.Node(
        interface=remove_throats_interface, name="iter_remove_throats_gm")
    iter_push_out_gm = pe.Node(
        interface=iter_push_out_interface, name="iter_push_out_gm")

    workflow.connect(
        [(join_gm_hemis, iter_remove_throats_gm, [("mesh_file", "outer_mesh")])])
    workflow.connect(
        [(subcort_small, iter_remove_throats_gm, [("mesh_file", "inner_mesh")])])
    workflow.connect(
        [(iter_remove_throats_gm, iter_push_out_gm, [("outer_mesh", "outer_mesh")])])
    workflow.connect(
        [(subcort_large, iter_push_out_gm, [("mesh_file", "inner_mesh")])])

    # Decoupling final GM from final WM

    clean_gm_before_decouple = create_iterative_cleanup_wf("gm_cleanup_iter")
    clean_wm_before_decouple = create_iterative_cleanup_wf("wm_cleanup_iter")

    workflow.connect(
        [(iter_push_out_wm, clean_wm_before_decouple, [("outer_mesh", "inputnode.in_file")])])
    workflow.connect(
        [(iter_push_out_gm, clean_gm_before_decouple, [("outer_mesh", "inputnode.in_file")])])

    decouple_surfaces_interface = util.Function(
        input_names=["outer_mesh", "inner_mesh"],
        output_names=["outer_mesh"],
        function=decouple_surfaces_fn)
    decouple_wm_gm = pe.Node(
        interface=decouple_surfaces_interface, name="decouple_wm_gm")

    workflow.connect(
        [(clean_wm_before_decouple, decouple_wm_gm, [("outputnode.out_file", "inner_mesh")])])
    workflow.connect(
        [(clean_gm_before_decouple, decouple_wm_gm, [("outputnode.out_file", "outer_mesh")])])

    # Fine tune for volume meshing
    # Gray matter
    finetune_gm = pe.Node(interface=mf.MeshFix(), name='finetune_gm')
    finetune_gm.inputs.quiet_mode = True
    finetune_gm.inputs.epsilon_angle = 2
    finetune_gm.inputs.number_of_biggest_shells = 2
    finetune_gm.inputs.finetuning_distance = 0.2
    finetune_gm.inputs.finetuning_substeps = 4
    finetune_gm.inputs.finetuning_outwards = True

    postfinetune_smooth_gm = pe.Node(
        interface=mf.MeshFix(), name='postfinetune_smooth_gm')
    postfinetune_smooth_gm.inputs.quiet_mode = True
    postfinetune_smooth_gm.inputs.epsilon_angle = 2
    postfinetune_smooth_gm.inputs.uniform_remeshing_steps = 1

    postfinetune_cleanup_gm = postfinetune_smooth_gm.clone(
        "postfinetune_cleanup_gm")
    postfinetune_cleanup_gm.inputs.uniform_remeshing_steps = 1

    workflow.connect(
        [(decouple_wm_gm, finetune_gm, [("outer_mesh", "in_file1")])])
    workflow.connect(
        [(finetune_gm, postfinetune_smooth_gm, [("mesh_file", "in_file1")])])
    workflow.connect(
        [(postfinetune_smooth_gm, postfinetune_cleanup_gm, [("mesh_file", "in_file1")])])

    # White matter
    finetune_wm = pe.Node(interface=mf.MeshFix(), name='finetune_wm')
    finetune_wm.inputs.quiet_mode = True
    finetune_wm.inputs.epsilon_angle = 2
    finetune_wm.inputs.number_of_biggest_shells = 2
    finetune_wm.inputs.finetuning_distance = 0.2
    finetune_wm.inputs.finetuning_substeps = 4
    finetune_wm.inputs.finetuning_inwards = True

    postfinetune_smooth_wm = postfinetune_smooth_gm.clone(
        "postfinetune_smooth_wm")
    postfinetune_cleanup_wm = postfinetune_cleanup_gm.clone(
        "postfinetune_cleanup_wm")

    workflow.connect(
        [(clean_wm_before_decouple, finetune_wm, [("outputnode.out_file", "in_file1")])])
    workflow.connect(
        [(finetune_wm, postfinetune_smooth_wm, [("mesh_file", "in_file1")])])
    workflow.connect(
        [(postfinetune_smooth_wm, postfinetune_cleanup_wm, [("mesh_file", "in_file1")])])

    remove_spikes_interface = util.Function(
        input_names=["outer_mesh", "inner_mesh"],
        output_names=["inner_mesh"],
        function=remove_spikes_fn)
    remove_spikes = pe.Node(
        interface=remove_spikes_interface, name="remove_spikes")

    workflow.connect(
        [(postfinetune_cleanup_wm, remove_spikes, [("mesh_file", "inner_mesh")])])
    workflow.connect(
        [(postfinetune_cleanup_gm, remove_spikes, [("mesh_file", "outer_mesh")])])

    # Final decoupling

    iter_decoupling_interface = util.Function(
        input_names=["outer_mesh", "inner_mesh"],
        output_names=["outer_mesh"],
        function=iter_decoupling_fn)
    final_GM_push_out = pe.Node(
        interface=iter_decoupling_interface, name="final_GM_push_out")

    workflow.connect(
        [(remove_spikes, final_GM_push_out, [("inner_mesh", "inner_mesh")])])
    workflow.connect(
        [(postfinetune_cleanup_gm, final_GM_push_out, [("mesh_file", "outer_mesh")])])

    # Save surfaces as stl and fsmesh
    wm_stl_floodfill_workflow = create_stl_fsmesh_floodfilled_wf(
        "wm_stl_floodfill")
    wm_stl_floodfill_workflow.inputs.input_to_STL.out_filename = "wm.stl"
    wm_stl_floodfill_workflow.inputs.filled_to_fsl_space.out_file = "wm.nii"

    gm_stl_floodfill_workflow = create_stl_fsmesh_floodfilled_wf(
        "gm_stl_floodfill")
    gm_stl_floodfill_workflow.inputs.input_to_STL.out_filename = "gm.stl"
    gm_stl_floodfill_workflow.inputs.filled_to_fsl_space.out_file = "gm.nii"

    workflow.connect(
        [(remove_spikes, wm_stl_floodfill_workflow, [("inner_mesh", "inputnode.in_file")])])
    workflow.connect(
        [(t1_to_fsl_space, wm_stl_floodfill_workflow, [("out_file", "inputnode.t1_fsl_space")])])

    workflow.connect(
        [(final_GM_push_out, gm_stl_floodfill_workflow, [("outer_mesh", "inputnode.in_file")])])
    workflow.connect(
        [(t1_to_fsl_space, gm_stl_floodfill_workflow, [("out_file", "inputnode.t1_fsl_space")])])


    workflow.connect(
        [(gm_stl_floodfill_workflow, outputnode, [("outputnode.filled_fsl_space", "gray_matter_volume")])])
    workflow.connect(
        [(gm_stl_floodfill_workflow, outputnode, [("outputnode.stl_mesh", "gray_matter_surface")])])

    workflow.connect(
        [(wm_stl_floodfill_workflow, outputnode, [("outputnode.filled_fsl_space", "white_matter_volume")])])
    workflow.connect(
        [(wm_stl_floodfill_workflow, outputnode, [("outputnode.stl_mesh", "white_matter_surface")])])


    workflow.connect([(get_conform_to_MNI, outputnode, [("out_file", "Conform2MNI")])])
    workflow.connect([(get_MNI_to_conform, outputnode, [("out_file", "MNI2Conform")])])

    workflow.connect(
        [(t1_to_fsl_space, outputnode, [("out_file", "t1_fsl_space")])])
    workflow.connect(
        [(labels_to_fsl_space, outputnode, [("out_file", "labels_fsl_space")])])
    return workflow


'''
Ventricle workflow
'''


def ventricle_workflow(name):
    inputnode = pe.Node(
        interface=util.IdentityInterface(
            fields=["t1_fsl_space", "segmentation", "white_matter_surface",
                    "white_matter_volume"]), name="inputnode")
    outputnode = pe.Node(
        interface=util.IdentityInterface(fields=["ventricle_surface", "ventricle_volume"]), name="outputnode")

    create_mask_from_ventricle_labels_interface = util.Function(
        input_names=["in_file", "label_values"],
        output_names=["out_file"],
        function=mask_from_labels_fn)
    create_ventricle_mask = pe.Node(
        interface=create_mask_from_ventricle_labels_interface, name="create_ventricle_mask")
    create_ventricle_mask.inputs.label_values = [4, 5, 14, 24, 31, 43, 44, 63]

    shrink_wm_volume = pe.Node(
        interface=fsl.ImageMaths(), name="shrink_wm_volume")
    shrink_wm_volume.inputs.op_string = "-ero -bin"

    # clean mask and ensure that ventricles are fully inside shrunken WM
    # (volume decoupling)
    shrink_ventricle_mask1 = pe.Node(
        interface=fsl.ImageMaths(), name="shrink_ventricle_mask1")
    shrink_ventricle_mask1.inputs.op_string = "-dilM -ero -bin"

    shrink_ventricle_mask2 = pe.Node(
        interface=fsl.MultiImageMaths(), name="shrink_ventricle_mask2")
    shrink_ventricle_mask2.inputs.op_string = "-mul %s -bin"

    ventricle_tess = create_tess_constrain_smooth_wf("ventricle_tess")

    remesh_ventricles = pe.Node(
        interface=mf.MeshFix(), name='remesh_ventricles')
    remesh_ventricles.inputs.epsilon_angle = 2.0
    remesh_ventricles.inputs.uniform_remeshing_steps = 5

    fix_ventricles1 = pe.Node(interface=mf.MeshFix(), name='fix_ventricles1')
    fix_ventricles1.inputs.epsilon_angle = 2.0

    fix_ventricles2 = pe.Node(interface=mf.MeshFix(), name='fix_ventricles2')
    fix_ventricles2.inputs.epsilon_angle = 2.0
    fix_ventricles2.inputs.uniform_remeshing_steps = 1

    fix_ventricles3 = fix_ventricles1.clone("fix_ventricles3")

    decouple_ventricles_interface = util.Function(
        input_names=["ventricles", "white_matter"],
        output_names=["ventricles"], function=decouple_ventricles_fn)

    decouple_ventricles = pe.Node(
        interface=decouple_ventricles_interface, name="decouple_ventricles")

    ventricle_stl_floodfill_workflow = create_stl_fsmesh_floodfilled_wf(
    "ventricle_stl_floodfill")
    ventricle_stl_floodfill_workflow.inputs.input_to_STL.out_filename = "ventricles.stl"
    ventricle_stl_floodfill_workflow.inputs.filled_to_fsl_space.out_file = "ventricles.nii"

    workflow = pe.Workflow(name=name + '_workflow')
    workflow.connect(
        [(inputnode, create_ventricle_mask, [("segmentation", "in_file")])])
    workflow.connect(
        [(inputnode, shrink_wm_volume, [("white_matter_volume", "in_file")])])
    workflow.connect(
        [(create_ventricle_mask, shrink_ventricle_mask1, [("out_file", "in_file")])])
    workflow.connect(
        [(shrink_ventricle_mask1, shrink_ventricle_mask2, [("out_file", "in_file")])])
    workflow.connect(
        [(shrink_wm_volume, shrink_ventricle_mask2, [("out_file", "operand_files")])])
    workflow.connect(
        [(shrink_ventricle_mask2, ventricle_tess, [("out_file", "inputnode.in_file")])])
    workflow.connect(
        [(ventricle_tess, remesh_ventricles, [("outputnode.out_file", "in_file1")])])
    workflow.connect(
        [(remesh_ventricles, fix_ventricles1, [("mesh_file", "in_file1")])])
    workflow.connect(
        [(fix_ventricles1, fix_ventricles2, [("mesh_file", "in_file1")])])
    workflow.connect(
        [(fix_ventricles2, fix_ventricles3, [("mesh_file", "in_file1")])])
    workflow.connect(
        [(fix_ventricles3, decouple_ventricles, [("mesh_file", "ventricles")])])
    workflow.connect(
        [(inputnode, decouple_ventricles, [("white_matter_surface", "white_matter")])])
    workflow.connect(
        [(decouple_ventricles, ventricle_stl_floodfill_workflow, [("ventricles", "inputnode.in_file")])])
    workflow.connect(
        [(inputnode, ventricle_stl_floodfill_workflow, [("t1_fsl_space", "inputnode.t1_fsl_space")])])
    workflow.connect(
        [(ventricle_stl_floodfill_workflow, outputnode, [("outputnode.stl_mesh", "ventricle_surface")])])
    workflow.connect(
        [(ventricle_stl_floodfill_workflow, outputnode, [("outputnode.filled_fsl_space", "ventricle_volume")])])
    return workflow


'''
End of ventricle workflow
'''

'''
Cerebellum workflow
'''


def cerebellum_workflow(name):
    inputnode = pe.Node(
        interface=util.IdentityInterface(
            fields=["t1_fsl_space", "segmentation", "gray_matter_surface",
                    "gray_matter_volume"]), name="inputnode")
    outputnode = pe.Node(
        interface=util.IdentityInterface(fields=["cerebellum_surface", "cerebellum_volume", "gray_matter_large"]), name="outputnode")

    create_mask_from_cerebellum_labels_interface = util.Function(
        input_names=["in_file", "label_values"],
        output_names=["out_file"],
        function=mask_from_labels_fn)
    create_cerebellum_mask = pe.Node(
        interface=create_mask_from_cerebellum_labels_interface, name="create_cerebellum_mask")
    create_cerebellum_mask.inputs.label_values = [7, 8, 16, 46, 47]

    enlarge_gm_volume = pe.Node(interface=fsl.ImageMaths(), name="enlarge_gm")
    enlarge_gm_volume.inputs.op_string = "-dilM -bin"

    # ensure that mask is fully outside enlarged GM (volume decoupling)
    shrink_cerebellum_mask1 = pe.Node(
        interface=fsl.MultiImageMaths(), name="shrink_cerebellum_mask1")
    shrink_cerebellum_mask1.inputs.op_string = "-sub %s -bin"

    shrink_cerebellum_mask2 = pe.Node(
        interface=fsl.ImageMaths(), name="shrink_cerebellum_mask2")
    shrink_cerebellum_mask2.inputs.op_string = "-ero -dilM -bin"

    cerebellum_tess = create_tess_constrain_smooth_wf("cerebellum_tess")

    remesh_cerebellum = pe.Node(
        interface=mf.MeshFix(), name='remesh_cerebellum')
    remesh_cerebellum.inputs.epsilon_angle = 2.0
    remesh_cerebellum.inputs.uniform_remeshing_steps = 5

    fix_cerebellum1 = pe.Node(interface=mf.MeshFix(), name='fix_cerebellum1')
    fix_cerebellum1.inputs.epsilon_angle = 2.0

    fix_cerebellum2 = pe.Node(interface=mf.MeshFix(), name='fix_cerebellum2')
    fix_cerebellum2.inputs.epsilon_angle = 2.0
    fix_cerebellum2.inputs.uniform_remeshing_steps = 1

    fix_cerebellum3 = fix_cerebellum1.clone("fix_cerebellum3")

    decouple_cerebellum_interface = util.Function(
        input_names=["mesh_file", "gray_matter"],
        output_names=["mesh_file"], function=decouple_input_from_GM_fn)

    decouple_cerebellum = pe.Node(
        interface=decouple_cerebellum_interface, name="decouple_cerebellum")

    cb_stl_floodfill_workflow = create_stl_fsmesh_floodfilled_wf(
        "cb_stl_floodfill")
    cb_stl_floodfill_workflow.inputs.input_to_STL.out_filename = "cerebellum.stl"
    cb_stl_floodfill_workflow.inputs.filled_to_fsl_space.out_file = "cerebellum.nii"

    workflow = pe.Workflow(name=name + '_workflow')
    workflow.connect(
        [(inputnode, create_cerebellum_mask, [("segmentation", "in_file")])])
    workflow.connect(
        [(inputnode, enlarge_gm_volume, [("gray_matter_volume", "in_file")])])
    workflow.connect(
        [(create_cerebellum_mask, shrink_cerebellum_mask1, [("out_file", "in_file")])])
    workflow.connect(
        [(enlarge_gm_volume, shrink_cerebellum_mask1, [("out_file", "operand_files")])])
    workflow.connect(
        [(shrink_cerebellum_mask1, shrink_cerebellum_mask2, [("out_file", "in_file")])])
    workflow.connect(
        [(shrink_cerebellum_mask2, cerebellum_tess, [("out_file", "inputnode.in_file")])])
    workflow.connect(
        [(cerebellum_tess, remesh_cerebellum, [("outputnode.out_file", "in_file1")])])
    workflow.connect(
        [(remesh_cerebellum, fix_cerebellum1, [("mesh_file", "in_file1")])])
    workflow.connect(
        [(fix_cerebellum1, fix_cerebellum2, [("mesh_file", "in_file1")])])
    workflow.connect(
        [(fix_cerebellum2, fix_cerebellum3, [("mesh_file", "in_file1")])])
    workflow.connect(
        [(fix_cerebellum3, decouple_cerebellum, [("mesh_file", "mesh_file")])])
    workflow.connect(
        [(inputnode, decouple_cerebellum, [("gray_matter_surface", "gray_matter")])])
    workflow.connect(
        [(decouple_cerebellum, cb_stl_floodfill_workflow, [("mesh_file", "inputnode.in_file")])])
    workflow.connect(
        [(inputnode, cb_stl_floodfill_workflow, [("t1_fsl_space", "inputnode.t1_fsl_space")])])
    workflow.connect(
        [(cb_stl_floodfill_workflow, outputnode, [("outputnode.stl_mesh", "cerebellum_surface")])])
    workflow.connect(
        [(cb_stl_floodfill_workflow, outputnode, [("outputnode.filled_fsl_space", "cerebellum_volume")])])
    workflow.connect(
        [(enlarge_gm_volume, outputnode, [("out_file", "gray_matter_large")])])
    return workflow


'''
End of cerebellum workflow
'''

'''
Cerebrospinal Fluid workflow
'''


def cerebrospinal_fluid_workflow(name):
    inputnode = pe.Node(
        interface=util.IdentityInterface(fields=["nu", "cerebellum_volume", "gray_matter_large", "t1_fsl_space", "Conform2MNI",
                                                 "gray_matter_surface", "cerebellum_surface"]), name="inputnode")
    outputnode = pe.Node(
        interface=util.IdentityInterface(fields=["csf_surface", "csf_volume", "NU_bet_meshfile"]), name="outputnode")

    nu_to_MNI = pe.Node(interface=fsl.ApplyXfm(), name='nu_to_MNI')
    nu_to_MNI.inputs.in_matrix_file = "/Users/erik/Dropbox/Analysis/TMSEEG/fs2fsl_conform.mat"

    get_initial_csf_surface = pe.Node(
        interface=fsl.BET(), name='get_initial_csf_surface')
    get_initial_csf_surface.inputs.mesh = True
    get_initial_csf_surface.inputs.mask = True
    get_initial_csf_surface.inputs.vertical_gradient = -0.2

    create_CSF_surface = pe.Node(
        interface=fsl.BETSurface(), name='create_CSF_surface')
    create_CSF_surface.inputs.mask = True
    create_CSF_surface.inputs.t1_only = True
    create_CSF_surface.inputs.outline = True

    # create enlarged csf volume and add it to enlarged gm volume
    create_enlarged_brainmask_with_CB = pe.Node(
        interface=fsl.MultiImageMaths(), name="create_enlarged_brainmask_with_CB")
    create_enlarged_brainmask_with_CB.inputs.op_string = "-dilM -add %s -bin"

    # ensure that CSF mask is fully outside brain_large (volume decoupling)
    csf_mask_outside_brain = pe.Node(
        interface=fsl.MultiImageMaths(), name="csf_mask_outside_brain")
    csf_mask_outside_brain.inputs.op_string = "-add %s -bin"

    csf_tess = create_tess_constrain_smooth_wf("csf_tess")

    decouple_outout_cutin_interface = util.Function(input_names=["outer_mesh", "inner_mesh"],
                                           output_names=["outer_mesh"], function=decouple_outout_cutin_fn)

    decouple_csf_from_gm = pe.Node(
        interface=decouple_outout_cutin_interface, name="decouple_csf_from_gm")
    decouple_csf_from_cerebellum = pe.Node(
        interface=decouple_outout_cutin_interface, name="decouple_csf_from_cerebellum")

    csf_stl_floodfill_workflow = create_stl_fsmesh_floodfilled_wf(
        "csf_stl_floodfill")
    csf_stl_floodfill_workflow.inputs.input_to_STL.out_filename = "csf.stl"
    csf_stl_floodfill_workflow.inputs.filled_to_fsl_space.out_file = "csf.nii"

    workflow = pe.Workflow(name=name + '_workflow')

    workflow.connect([(inputnode, nu_to_MNI, [("nu", "in_file")])])
    workflow.connect([(inputnode, nu_to_MNI, [("t1_fsl_space", "reference")])])
    workflow.connect(
            [(nu_to_MNI, get_initial_csf_surface, [("out_file", "in_file")])])
    workflow.connect(
        [(get_initial_csf_surface, create_CSF_surface, [("out_file", "t1_file")])])
    workflow.connect(
        [(get_initial_csf_surface, create_CSF_surface, [("meshfile", "meshfile")])])
    workflow.connect(
        [(get_initial_csf_surface, outputnode, [("meshfile", "NU_bet_meshfile")])])
    workflow.connect(
        [(inputnode, create_CSF_surface, [("Conform2MNI", "t1_to_standard_matrix")])])

    workflow.connect(
        [(inputnode, create_enlarged_brainmask_with_CB, [("cerebellum_volume", "in_file")])])
    workflow.connect(
        [(inputnode, create_enlarged_brainmask_with_CB, [("gray_matter_large", "operand_files")])])

    workflow.connect(
        [(create_enlarged_brainmask_with_CB, csf_mask_outside_brain, [("out_file", "operand_files")])])
    workflow.connect(
        [(create_CSF_surface, csf_mask_outside_brain, [("inskull_mask_file", "in_file")])])
    workflow.connect(
        [(csf_mask_outside_brain, csf_tess, [("out_file", "inputnode.in_file")])])
    workflow.connect(
        [(csf_tess, decouple_csf_from_gm, [("outputnode.out_file", "outer_mesh")])])
    workflow.connect(
        [(inputnode, decouple_csf_from_gm, [("gray_matter_surface", "inner_mesh")])])
    workflow.connect(
        [(inputnode, decouple_csf_from_cerebellum, [("cerebellum_surface", "inner_mesh")])])
    workflow.connect(
        [(decouple_csf_from_gm, decouple_csf_from_cerebellum, [("outer_mesh", "outer_mesh")])])
    workflow.connect(
        [(decouple_csf_from_cerebellum, csf_stl_floodfill_workflow, [("outer_mesh", "inputnode.in_file")])])
    workflow.connect(
        [(inputnode, csf_stl_floodfill_workflow, [("t1_fsl_space", "inputnode.t1_fsl_space")])])
    workflow.connect(
        [(csf_stl_floodfill_workflow, outputnode, [("outputnode.stl_mesh", "csf_surface")])])
    workflow.connect(
        [(csf_stl_floodfill_workflow, outputnode, [("outputnode.filled_fsl_space", "csf_volume")])])
    return workflow

'''
End of Cerebrospinal Fluid workflow
'''

'''
Skull workflow
'''


def skull_workflow(name):
    inputnode = pe.Node(
        interface=util.IdentityInterface(fields=["csf_volume", "csf_surface",
                                                 "t1_fsl_space","Conform2MNI",
                                                 "NU_bet_meshfile"]), name="inputnode")
    outputnode = pe.Node(
        interface=util.IdentityInterface(fields=["skull_volume", "skull_surface"]), name="outputnode")

    create_skullmask = pe.Node(
        interface=fsl.BETSurface(), name='create_skullmask')
    create_skullmask.inputs.mask = True
    create_skullmask.inputs.outline = True
    create_skullmask.inputs.t1_only = True

    # create enlarged csf volume and add it to enlarged gm volume
    create_enlarged_CSF = pe.Node(
        interface=fsl.ImageMaths(), name="create_enlarged_CSF")
    create_enlarged_CSF.inputs.op_string = "-dilM -bin"

    # ensure that skull mask is fully outside csf_large (volume decoupling)
    skull_outside_csf = pe.Node(
        interface=fsl.MultiImageMaths(), name="skull_outside_csf")
    skull_outside_csf.inputs.op_string = "-add %s -bin"

    n_vertices = 60000
    skull_tess = create_tess_constrain_smooth_wf("skull_tess")
    skull_tess.inputs.constrain_vertices.uniform_remeshing_vertices = n_vertices / 3

    decouple_outout_cutin_interface = util.Function(input_names=["outer_mesh", "inner_mesh"],
                                           output_names=["outer_mesh"], function=decouple_outout_cutin_fn)

    # Uses decouple-outout and cut-inner, same as the CSF function,
    # so the same interface is used.
    decouple_skull_from_csf = pe.Node(
        interface=decouple_outout_cutin_interface, name="decouple_skull_from_csf")

    skull_stl_floodfill_workflow = create_stl_fsmesh_floodfilled_wf(
        "skull_stl_floodfill")
    skull_stl_floodfill_workflow.inputs.input_to_STL.out_filename = "skull.stl"
    skull_stl_floodfill_workflow.inputs.filled_to_fsl_space.out_file = "skull.nii"
    skull_stl_floodfill_workflow.inputs.floodfill.resolution = 0.5

    workflow = pe.Workflow(name=name + '_workflow')
    workflow.connect(
        [(inputnode, create_enlarged_CSF, [("csf_volume", "in_file")])])

    workflow.connect(
        [(inputnode, create_skullmask, [("t1_fsl_space", "t1_file")])])
    workflow.connect(
        [(inputnode, create_skullmask, [("NU_bet_meshfile", "meshfile")])])

    workflow.connect(
        [(inputnode, create_skullmask, [("Conform2MNI", "t1_to_standard_matrix")])])

    workflow.connect(
        [(create_skullmask, skull_outside_csf, [("outskull_mask_file", "in_file")])])
    workflow.connect(
        [(create_enlarged_CSF, skull_outside_csf, [("out_file", "operand_files")])])
    workflow.connect(
        [(skull_outside_csf, skull_tess, [("out_file", "inputnode.in_file")])])

    ## These are reversed on purpose!
    workflow.connect(
        [(skull_tess, decouple_skull_from_csf, [("outputnode.out_file", "outer_mesh")])])
    workflow.connect(
        [(inputnode, decouple_skull_from_csf, [("csf_surface", "inner_mesh")])])
    workflow.connect(
        [(decouple_skull_from_csf, skull_stl_floodfill_workflow, [("outer_mesh", "inputnode.in_file")])])
    ## We are decoupling the skull from the csf!

    workflow.connect(
        [(inputnode, skull_stl_floodfill_workflow, [("t1_fsl_space", "inputnode.t1_fsl_space")])])
    workflow.connect(
        [(skull_stl_floodfill_workflow, outputnode, [("outputnode.stl_mesh", "skull_surface")])])
    workflow.connect(
        [(skull_stl_floodfill_workflow, outputnode, [("outputnode.filled_fsl_space", "skull_volume")])])
    return workflow
'''
End of Skull workflow
'''


'''
Skin workflow
'''


def skin_workflow(name):
    inputnode = pe.Node(
        interface=util.IdentityInterface(fields=["t1_fsl_space", "Conform2MNI",
                "skull_volume", "skull_surface","NU_bet_meshfile"]), name="inputnode")
    outputnode = pe.Node(
        interface=util.IdentityInterface(fields=["skin_surface", "skin_volume"]), name="outputnode")

    # get top slice of skin mask
    # w=$(( `fslstats $TMP_DIR/skin_raw.nii.gz -w | awk '{print $5}'` + `fslstats $TMP_DIR/skin_raw.nii.gz -w | awk '{print $6}'` - 1 ))
    # fslmaths $TMP_DIR/T1_conform -mul $TMP_DIR/skin_raw.nii.gz -roi 0 256 0 256 $w 1 0 1 $TMP_DIR/skin_top_slice.nii.gz
    # get mean intensity in top slice of T1_conform
    # skin_threshold=`fslstats $TMP_DIR/skin_top_slice.nii.gz -M`

    create_skin_mask = pe.Node(
        interface=fsl.BETSurface(), name='create_skin_mask')
    create_skin_mask.inputs.mask = True
    create_skin_mask.inputs.outline = True
    create_skin_mask.inputs.t1_only = True

    # create enlarged csf volume and add it to enlarged gm volume
    create_enlarged_skull = pe.Node(
        interface=fsl.ImageMaths(), name="create_enlarged_skull")
    create_enlarged_skull.inputs.op_string = "-dilM -bin"

    # ensure that skin mask is fully outside csf_large (volume decoupling)
    skin_outside_skull = pe.Node(
        interface=fsl.MultiImageMaths(), name="skin_outside_skull")
    skin_outside_skull.inputs.op_string = "-add %s -bin"

    n_vertices = 60000
    skin_tess = create_tess_constrain_smooth_wf("skin_tess")
    skin_tess.inputs.constrain_vertices.uniform_remeshing_vertices = n_vertices / 4


    decouple_outout_cutin_interface = util.Function(input_names=["outer_mesh", "inner_mesh"],
                                           output_names=["outer_mesh"], function=decouple_outout_cutin_fn)
    # Uses decouple-outout and cut-inner, same as the CSF function,
    # so the same interface is used.
    decouple_skin_from_skull = pe.Node(
        interface=decouple_outout_cutin_interface, name="decouple_skin_from_skull")

    skin_stl_floodfill_workflow = create_stl_fsmesh_floodfilled_wf(
        "skin_stl_floodfill")
    skin_stl_floodfill_workflow.inputs.input_to_STL.out_filename = "skin.stl"
    skin_stl_floodfill_workflow.inputs.filled_to_fsl_space.out_file = "skin.nii"
    skin_stl_floodfill_workflow.inputs.floodfill.resolution = 0.5

    workflow = pe.Workflow(name=name + '_workflow')
    workflow.connect(
        [(inputnode, create_skin_mask, [("Conform2MNI", "t1_to_standard_matrix")])])
    workflow.connect(
        [(inputnode, create_skin_mask, [("t1_fsl_space", "t1_file")])])
    workflow.connect(
        [(inputnode, create_skin_mask, [("NU_bet_meshfile", "meshfile")])])
    workflow.connect(
        [(inputnode, create_enlarged_skull, [("skull_volume", "in_file")])])
    workflow.connect(
        [(create_skin_mask, skin_outside_skull, [("outskin_mask_file", "in_file")])])
    workflow.connect(
        [(create_enlarged_skull, skin_outside_skull, [("out_file", "operand_files")])])
    workflow.connect(
        [(skin_outside_skull, skin_tess, [("out_file", "inputnode.in_file")])])
    workflow.connect(
        [(skin_tess, decouple_skin_from_skull, [("outputnode.out_file", "outer_mesh")])])
    workflow.connect(
        [(inputnode, decouple_skin_from_skull, [("skull_surface", "inner_mesh")])])
    workflow.connect(
        [(decouple_skin_from_skull, skin_stl_floodfill_workflow, [("outer_mesh", "inputnode.in_file")])])
    workflow.connect(
        [(inputnode, skin_stl_floodfill_workflow, [("t1_fsl_space", "inputnode.t1_fsl_space")])])
    workflow.connect(
        [(skin_stl_floodfill_workflow, outputnode, [("outputnode.stl_mesh", "skin_surface")])])
    workflow.connect(
        [(skin_stl_floodfill_workflow, outputnode, [("outputnode.filled_fsl_space", "skin_volume")])])
    return workflow
'''
End of Skin workflow
'''


'''
Main structural preprocessing workflow
'''

def create_structural_mesh_workflow(name="structural_mesh"):
    workflow = pe.Workflow(name=name)
    inputnode = pe.Node(
        interface=util.IdentityInterface(fields=["subject_id", "subjects_dir"]), name="inputnode")
    outputnode = pe.Node(
        interface=util.IdentityInterface(fields=["volume_mesh", "surfaces", "volumes"]), name="outputnode")

    brain_wf = brain_workflow("brain")
    ventricle_wf = ventricle_workflow("ventricle")
    cerebellum_wf = cerebellum_workflow("cerebellum")
    cerebrospinal_fluid_wf = cerebrospinal_fluid_workflow("cerebrospinal_fluid")
    skull_wf = skull_workflow("skull")
    skin_wf = skin_workflow("skin")

    create_volume_mesh_script_interface = util.Function(input_names=["subject_id", "gm", "wm", "csf","skull", "skin", "cerebellum", "ventricles"],
                                           output_names=["out_file"], function=volume_mesh_fn)

    create_volume_mesh_script = pe.Node(interface=create_volume_mesh_script_interface, name="create_volume_mesh_script")
    run_volume_meshing = pe.Node(interface=Gmsh(), name="run_volume_meshing")
    run_volume_meshing.inputs.mesh_generation_dimension = '3'

    mergeSurfaces = pe.Node(
        interface=util.Merge(7), name='mergeSurfaces')
    mergeVolumes = pe.Node(
        interface=util.Merge(7), name='mergeVolumes')

    workflow.connect(
        [(inputnode, brain_wf, [("subject_id", "inputnode.subject_id")])])
    workflow.connect(
        [(inputnode, brain_wf, [("subjects_dir", "inputnode.subjects_dir")])])

    workflow.connect(
        [(brain_wf, ventricle_wf, [("outputnode.labels_fsl_space", "inputnode.segmentation")])])
    workflow.connect(
        [(brain_wf, ventricle_wf, [("outputnode.white_matter_volume", "inputnode.white_matter_volume")])])
    workflow.connect(
        [(brain_wf, ventricle_wf, [("outputnode.white_matter_surface", "inputnode.white_matter_surface")])])
    workflow.connect(
        [(brain_wf, ventricle_wf, [("outputnode.t1_fsl_space", "inputnode.t1_fsl_space")])])

    workflow.connect(
        [(brain_wf, cerebellum_wf, [("outputnode.labels_fsl_space", "inputnode.segmentation")])])
    workflow.connect(
        [(brain_wf, cerebellum_wf, [("outputnode.gray_matter_volume", "inputnode.gray_matter_volume")])])
    workflow.connect(
        [(brain_wf, cerebellum_wf, [("outputnode.gray_matter_surface", "inputnode.gray_matter_surface")])])
    workflow.connect(
        [(brain_wf, cerebellum_wf, [("outputnode.t1_fsl_space", "inputnode.t1_fsl_space")])])


    workflow.connect(
        [(brain_wf, cerebrospinal_fluid_wf, [("outputnode.t1_fsl_space", "inputnode.t1_fsl_space")])])
    workflow.connect(
        [(brain_wf, cerebrospinal_fluid_wf, [("outputnode.nu", "inputnode.nu")])])
    workflow.connect(
        [(brain_wf, cerebrospinal_fluid_wf, [("outputnode.Conform2MNI", "inputnode.Conform2MNI")])])
    workflow.connect(
        [(brain_wf, cerebrospinal_fluid_wf, [("outputnode.gray_matter_surface", "inputnode.gray_matter_surface")])])
    workflow.connect(
        [(cerebellum_wf, cerebrospinal_fluid_wf, [("outputnode.cerebellum_volume", "inputnode.cerebellum_volume")])])
    workflow.connect(
        [(cerebellum_wf, cerebrospinal_fluid_wf, [("outputnode.cerebellum_surface", "inputnode.cerebellum_surface")])])
    workflow.connect(
        [(cerebellum_wf, cerebrospinal_fluid_wf, [("outputnode.gray_matter_large", "inputnode.gray_matter_large")])])

    workflow.connect(
        [(brain_wf, skull_wf, [("outputnode.t1_fsl_space", "inputnode.t1_fsl_space")])])
    workflow.connect(
        [(brain_wf, skull_wf, [("outputnode.Conform2MNI", "inputnode.Conform2MNI")])])
    workflow.connect(
        [(cerebrospinal_fluid_wf, skull_wf, [("outputnode.csf_volume", "inputnode.csf_volume")])])
    workflow.connect(
        [(cerebrospinal_fluid_wf, skull_wf, [("outputnode.csf_surface", "inputnode.csf_surface")])])
    workflow.connect(
        [(cerebrospinal_fluid_wf, skull_wf, [("outputnode.NU_bet_meshfile", "inputnode.NU_bet_meshfile")])])


    workflow.connect(
        [(brain_wf, skin_wf, [("outputnode.t1_fsl_space", "inputnode.t1_fsl_space")])])
    workflow.connect(
        [(brain_wf, skin_wf, [("outputnode.Conform2MNI", "inputnode.Conform2MNI")])])
    workflow.connect(
        [(skull_wf, skin_wf, [("outputnode.skull_volume", "inputnode.skull_volume")])])
    workflow.connect(
        [(skull_wf, skin_wf, [("outputnode.skull_surface", "inputnode.skull_surface")])])
    workflow.connect(
        [(cerebrospinal_fluid_wf, skin_wf, [("outputnode.NU_bet_meshfile", "inputnode.NU_bet_meshfile")])])

    workflow.connect(
        [(brain_wf, mergeSurfaces, [("outputnode.gray_matter_surface", "in1")])])
    workflow.connect(
        [(brain_wf, mergeSurfaces, [("outputnode.white_matter_surface", "in2")])])
    workflow.connect(
        [(ventricle_wf, mergeSurfaces, [("outputnode.ventricle_surface", "in3")])])
    workflow.connect(
        [(cerebellum_wf, mergeSurfaces, [("outputnode.cerebellum_surface", "in4")])])
    workflow.connect(
        [(cerebrospinal_fluid_wf, mergeSurfaces, [("outputnode.csf_surface", "in5")])])
    workflow.connect(
        [(skull_wf, mergeSurfaces, [("outputnode.skull_surface", "in6")])])
    workflow.connect(
        [(skin_wf, mergeSurfaces, [("outputnode.skin_surface", "in7")])])

    workflow.connect(
        [(inputnode, create_volume_mesh_script, [("subject_id", "subject_id")])])
    workflow.connect(
        [(brain_wf, create_volume_mesh_script, [("outputnode.gray_matter_surface", "gm")])])
    workflow.connect(
        [(brain_wf, create_volume_mesh_script, [("outputnode.white_matter_surface", "wm")])])
    workflow.connect(
        [(ventricle_wf, create_volume_mesh_script, [("outputnode.ventricle_surface", "ventricles")])])
    workflow.connect(
        [(cerebellum_wf, create_volume_mesh_script, [("outputnode.cerebellum_surface", "cerebellum")])])
    workflow.connect(
        [(cerebrospinal_fluid_wf, create_volume_mesh_script, [("outputnode.csf_surface", "csf")])])
    workflow.connect(
        [(skull_wf, create_volume_mesh_script, [("outputnode.skull_surface", "skull")])])
    workflow.connect(
        [(skin_wf, create_volume_mesh_script, [("outputnode.skin_surface", "skin")])])

    workflow.connect(
        [(create_volume_mesh_script, run_volume_meshing, [("out_file", "in_files")])])

    workflow.connect(
        [(brain_wf, mergeVolumes, [("outputnode.gray_matter_volume", "in1")])])
    workflow.connect(
        [(brain_wf, mergeVolumes, [("outputnode.white_matter_volume", "in2")])])
    workflow.connect(
        [(ventricle_wf, mergeVolumes, [("outputnode.ventricle_volume", "in3")])])
    workflow.connect(
        [(cerebellum_wf, mergeVolumes, [("outputnode.cerebellum_volume", "in4")])])
    workflow.connect(
        [(cerebrospinal_fluid_wf, mergeVolumes, [("outputnode.csf_volume", "in5")])])
    workflow.connect(
        [(skull_wf, mergeVolumes, [("outputnode.skull_volume", "in6")])])
    workflow.connect(
        [(skin_wf, mergeVolumes, [("outputnode.skin_volume", "in7")])])

    workflow.connect([(mergeSurfaces, outputnode, [("out", "surfaces")])])
    workflow.connect([(mergeVolumes, outputnode, [("out", "volumes")])])
    workflow.connect([(run_volume_meshing, outputnode, [("mesh_file", "volume_mesh")])])
    return workflow

'''
End of main structural preprocessing workflow
'''
