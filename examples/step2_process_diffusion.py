import os
import os.path as op
import nipype.interfaces.io as nio           # Data i/o
import nipype.interfaces.utility as util     # utility
import nipype.pipeline.engine as pe          # pypeline engine
from forward.dti import create_conductivity_tensor_mesh_workflow

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

preproc = create_conductivity_tensor_mesh_workflow()

preproc.inputs.inputnode.mesh_file = "TMS007_running.msh"
preproc.inputs.inputnode.struct = "/Users/erik/Dropbox/Analysis/TMSEEG/TMS007/orig_out_flirt.nii"

datasink = pe.Node(interface=nio.DataSink(),
                   name="datasink")
datasink.inputs.base_directory = op.abspath('structural_datasink')
datasink.inputs.container = 'subject'

dti_proc = pe.Workflow(name="dti_proc")
dti_proc.base_dir = os.path.abspath('dti_proc')
dti_proc.connect([
                    (infosource,datasource,[('subject_id', 'subject_id')]),
                    (datasource,preproc,[('dwi','inputnode.dwi'),
                                               ('bvals','inputnode.bvals'),
                                               ('bvecs','inputnode.bvecs'),
                                               ])
                ])
dti_proc.connect([(preproc, datasink, [("outputnode.mesh_file", "mesh_file")])])
dti_proc.connect([(infosource, datasink, [("subject_id", "subject_id")])])

if __name__ == '__main__':
    dti_proc.write_graph()
    #dti_proc.run()
