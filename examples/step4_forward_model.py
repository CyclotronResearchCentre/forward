'''
The fourth step is to run the actual forward modelling and generate
a lead field matrix. This is done with the concept of reciprocity.
'''


import os
import os.path as op
import nipype.interfaces.io as nio           # Data i/o
import nipype.interfaces.utility as util     # utility
import nipype.pipeline.engine as pe          # pypeline engine
from forward.leadfield import create_forward_model_workflow

from forward.datasets import sample
data_path = sample.data_path()

'''
Define the markers, ground electrode, and any bad electrodes to ignore
'''

markers = ['LeftEar', 'RightEar', 'NZ']
ground_electrode = 'IZ'

'''
Define the mesh file to be used. This should be the volume mesh 
output of the structural preprocessing workflow.
'''
conductivity_tensor_included = True

'''
Define the subjects
'''

subject_list = ['TMS007']


'''
Create the nodes and forward modelling workflow
'''

infosource = pe.Node(interface=util.IdentityInterface(
    fields=['subject_id']), name="infosource")
infosource.iterables = ('subject_id', subject_list)

info = dict(mesh_file=[['subject_id']],
            electrode_name_file=[['subject_id']])

datasource = pe.Node(interface=nio.DataGrabber(infields=['subject_id'],
                                               outfields=info.keys()),
                     name='datasource')

datasource.inputs.sort_filelist = True
datasource.inputs.template = "%s"
datasource.inputs.base_directory = data_path
datasource.inputs.field_template = dict(
    mesh_file='../%s*.msh', electrode_name_file='%s/ElectrodeNames.txt')
datasource.inputs.template_args = info

fwd = create_forward_model_workflow("forward", conductivity_tensor_included)
fwd.inputs.inputnode.marker_list = markers
fwd.inputs.inputnode.ground_electrode = ground_electrode
fwd.inputs.inputnode.conductivity_tensor_included = conductivity_tensor_included

datasink = pe.Node(interface=nio.DataSink(),
                   name="datasink")
datasink.inputs.base_directory = op.abspath('forward_datasink')
datasink.inputs.container = 'subject'

fwd_proc = pe.Workflow(name="fwd_proc")
fwd_proc.base_dir = os.path.abspath('fwd_proc')
fwd_proc.connect([
                    (infosource,datasource,[('subject_id', 'subject_id')]),
                    (datasource,fwd,[('mesh_file','inputnode.mesh_file'),
                                               ('electrode_name_file','inputnode.electrode_name_file'),
                                               ])
                ])
fwd_proc.connect([(fwd, datasink, [("outputnode.leadfield", "leadfield")])])
fwd_proc.connect([(infosource, datasink, [("subject_id", "subject_id")])])

if __name__ == '__main__':
    fwd_proc.write_graph(graph2use="exec")
    fwd_proc.run()
    #fwd_proc.run(plugin='MultiProc', plugin_args={'n_procs' : 4})
