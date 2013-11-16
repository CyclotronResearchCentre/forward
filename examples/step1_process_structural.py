import os
import os.path as op
import nipype.interfaces.io as nio           # Data i/o
import nipype.interfaces.utility as util     # utility
import nipype.pipeline.engine as pe          # pypeline engine
from forward.struct import create_structural_mesh_workflow

data_dir = op.abspath(op.curdir)
subject_list = ['fsaverage']

infosource = pe.Node(interface=util.IdentityInterface(
    fields=['subject_id', 'subjects_dir']), name="infosource")
infosource.iterables = ('subject_id', subject_list)
infosource.inputs.subjects_dir = os.environ["SUBJECTS_DIR"]

datasink = pe.Node(interface=nio.DataSink(),
                   name="datasink")
datasink.inputs.base_directory = op.abspath('structural_datasink')
datasink.inputs.container = 'subject'

preproc = create_structural_mesh_workflow()

struct_proc = pe.Workflow(name="struct_proc")
struct_proc.base_dir = op.abspath('struct_proc')
struct_proc.connect([
                    (infosource,preproc,[('subject_id', 'inputnode.subject_id')]),
                    (infosource,preproc,[('subjects_dir', 'inputnode.subjects_dir')]),
                ])

struct_proc.connect([
                    (preproc,datasink,[('outputnode.volumes', 'volumes')]),
                    (preproc,datasink,[('outputnode.surfaces', 'surfaces')]),
                    (preproc,datasink,[('outputnode.surfaces', 'volume_mesh')]),
                ])

struct_proc.connect([(infosource, datasink, [("subject_id", "subject_id")])])


if __name__ == '__main__':
    struct_proc.write_graph()
    #struct_proc.run()
