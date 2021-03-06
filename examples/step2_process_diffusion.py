import os
import os.path as op
import nipype.interfaces.io as nio           # Data i/o
import nipype.interfaces.utility as util     # utility
import nipype.pipeline.engine as pe          # pypeline engine
from forward.dti import create_conductivity_tensor_mesh_workflow

from forward.datasets import sample
data_path = sample.data_path()

subject_list = ['TMS007']

infosource = pe.Node(interface=util.IdentityInterface(
    fields=['subject_id']), name="infosource")
infosource.iterables = ('subject_id', subject_list)

info = dict(dwi=[['subject_id', 'DWI_1000']],
            bvecs=[['subject_id', 'grad_1000_invZ']],
            bvals=[['subject_id', 'bval_1000']],
            mesh_file=[['subject_id', 'subject_id']],
            struct=[['subject_id', 'T1mprage']],
            t1_fsl_space=[['subject_id', 'orig']])

datasource = pe.Node(interface=nio.DataGrabber(infields=['subject_id'],
                                               outfields=info.keys()),
                     name='datasource')

datasource.inputs.sort_filelist = True
datasource.inputs.template = "%s/%s"
datasource.inputs.base_directory = data_path
datasource.inputs.field_template = dict(
    dwi='%s/*%s.nii.gz', bvecs='%s/*%s', bvals='%s/*%s',
    mesh_file='../structural_datasink/subject/volume_mesh/*%s/%s*.msh',
    t1_fsl_space='../structural_datasink/subject/t1_fsl_space/*%s/%s*.nii.gz',
    struct=op.join(data_path,'%s/%s*.nii'))
datasource.inputs.template_args = info

preproc = create_conductivity_tensor_mesh_workflow()

datasink = pe.Node(interface=nio.DataSink(),
                   name="datasink")
datasink.inputs.base_directory = op.abspath('diffusion_datasink')
datasink.inputs.container = 'subject'

dti_proc = pe.Workflow(name="dti_proc")
dti_proc.base_dir = os.path.abspath('dti_proc')
dti_proc.connect([
                    (infosource,datasource,[('subject_id', 'subject_id')]),
                    (datasource,preproc,[('dwi','inputnode.dwi'),
                                               ('bvals','inputnode.bvals'),
                                               ('bvecs','inputnode.bvecs'),
                                               ('mesh_file','inputnode.mesh_file'),
                                               ('struct','inputnode.struct'),
                                               ('t1_fsl_space','inputnode.t1_fsl_space'),
                                               ])
                ])
dti_proc.connect([(preproc, datasink, [("outputnode.mesh_file", "mesh_file")])])
dti_proc.connect([(preproc, datasink, [("outputnode.diff_V1", "diff_tensor")])])
dti_proc.connect([(preproc, datasink, [("outputnode.cond_V1", "cond_tensor")])])
dti_proc.connect([(preproc, datasink, [("outputnode.mean_conductivity", "mean_conductivity")])])
dti_proc.connect([(preproc, datasink, [("outputnode.fa_t1_space", "fa_t1_space")])])
dti_proc.connect([(infosource, datasink, [("subject_id", "subject_id")])])

if __name__ == '__main__':
    dti_proc.write_graph()

    import time
    start = time.time()
    dti_proc.run()
    end = time.time()
    print(start)
    print(end)
    print(end-start)
    #dti_proc.run(plugin='MultiProc', plugin_args={'n_procs' : 4})
