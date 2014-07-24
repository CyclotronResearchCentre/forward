import os.path as op
import nipype.interfaces.io as nio           # Data i/o
import nipype.interfaces.utility as util     # utility
import nipype.interfaces.fsl as fsl
import nipype.pipeline.engine as pe          # pypeline engine
from forward.struct import create_structural_mesh_workflow

from forward.datasets import sample
data_path = sample.data_path(name="simnibs")
almi5_path = sample.data_path(name="almi5")

subjects_dir = op.abspath(op.curdir)
# Normally something like:
# subjects_dir = os.environ["SUBJECTS_DIR"]

fsl.FSLCommand.set_default_output_type('NIFTI_GZ')

# The fat-suppressed T1 should be run through freesurfer!
# See http://simnibs.de/documentation/mri_sequences

subject_list = ['almi5']
name = 'FourImage_'
infosource = pe.Node(interface=util.IdentityInterface(
    fields=['subject_id', 'subjects_dir']), name="infosource")
infosource.iterables = ('subject_id', subject_list)
infosource.inputs.subjects_dir = subjects_dir

info = dict(t2=[['subject_id', '_T2']],
            t2fs=[['subject_id', '_T2fs']],
            t1=[['subject_id', '_T1']])

datasource = pe.Node(interface=nio.DataGrabber(infields=['subject_id'],
                                               outfields=info.keys()),
                     name='datasource')

datasource.inputs.sort_filelist = True
datasource.inputs.template = "%s/%s"
datasource.inputs.base_directory = data_path
datasource.inputs.field_template = dict(
    t2=op.join(data_path,'mesh/%s%s.nii.gz'),
    t2fs=op.join(data_path,'mesh/%s%s.nii.gz'),
    t1=op.join(data_path,'mesh/%s%s.nii.gz'))
datasource.inputs.template_args = info


datasink = pe.Node(interface=nio.DataSink(),
                   name="datasink")
datasink.inputs.base_directory = op.abspath(name + 'structural_datasink')
datasink.inputs.container = 'subject'

preproc = create_structural_mesh_workflow(include_t1=True, include_t2=True, include_t2fs=True)

struct_proc = pe.Workflow(name="struct_proc")
struct_proc.base_dir = op.abspath(name + 'struct_proc')
struct_proc.connect([
                    (infosource,preproc,[('subject_id', 'inputnode.subject_id')]),
                    (infosource,preproc,[('subjects_dir', 'inputnode.subjects_dir')]),
                    (infosource,datasource,[('subject_id', 'subject_id')]),
                ])

struct_proc.connect([
                    (datasource,preproc,[('t2', 'inputnode.t2_file')]),
                    (datasource,preproc,[('t2fs', 'inputnode.t2fs_file')]),
                    (datasource,preproc,[('t1', 'inputnode.t1_file')]),
                ])

struct_proc.connect([
                    (preproc,datasink,[('outputnode.volumes', 'volumes')]),
                    (preproc,datasink,[('outputnode.surfaces', 'surfaces')]),
                    (preproc,datasink,[('outputnode.volume_mesh', 'volume_mesh')]),
                    (preproc,datasink,[('outputnode.t1_fsl_space', 't1_fsl_space')]),
                ])

struct_proc.connect([(infosource, datasink, [("subject_id", "subject_id")])])


if __name__ == '__main__':
    struct_proc.write_graph()
    struct_proc.run()
    #struct_proc.run(plugin='MultiProc', plugin_args={'n_procs' : 4}, updatehash=False)

