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
name = 'SingleImage_'
infosource = pe.Node(interface=util.IdentityInterface(
    fields=['subject_id', 'subjects_dir']), name="infosource")
infosource.iterables = ('subject_id', subject_list)
infosource.inputs.subjects_dir = subjects_dir

datasink = pe.Node(interface=nio.DataSink(),
                   name="datasink")
datasink.inputs.base_directory = op.abspath(name + 'structural_datasink')
datasink.inputs.container = 'subject'

preproc = create_structural_mesh_workflow()

struct_proc = pe.Workflow(name="struct_proc")
struct_proc.base_dir = op.abspath(name + 'struct_proc')
struct_proc.connect([
                    (infosource,preproc,[('subject_id', 'inputnode.subject_id')]),
                    (infosource,preproc,[('subjects_dir', 'inputnode.subjects_dir')]),
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

