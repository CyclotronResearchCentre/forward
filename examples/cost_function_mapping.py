'''
Here we use the previously derived leadfield matrix and an example dipole
to map the cost function in the gray matter elements. 
'''


import os
import os.path as op
import nipype.interfaces.io as nio           # Data i/o
import nipype.interfaces.utility as util     # utility
import nipype.pipeline.engine as pe          # pypeline engine
from forward.dipole import create_cost_function_workflow

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
            electrode_name_file=[['subject_id']],
            dipole_file=[['subject_id']],
            leadfield=[['subject_id']])

datasource = pe.Node(interface=nio.DataGrabber(infields=['subject_id'],
                                               outfields=info.keys()),
                     name='datasource')

datasource.inputs.sort_filelist = True
datasource.inputs.template = "%s"
datasource.inputs.base_directory = data_path
datasource.inputs.field_template = dict(
    mesh_file='../%s_gmsh_cond_elec.msh', electrode_name_file='%s/ElectrodeNames.txt',
    dipole_file="%s/*-lh.dip", leadfield="../forward_datasink/leadfield/%s_leadfield.hdf5")
datasource.inputs.template_args = info

cost = create_cost_function_workflow("dipole_cost", conductivity_tensor_included)
cost.inputs.inputnode.marker_list = markers
cost.inputs.inputnode.ground_electrode = ground_electrode
cost.inputs.inputnode.conductivity_tensor_included = conductivity_tensor_included

datasink = pe.Node(interface=nio.DataSink(),
                   name="datasink")
datasink.inputs.base_directory = op.abspath('forward_datasink')
datasink.inputs.container = 'subject'

cost_proc = pe.Workflow(name="cost_proc")
cost_proc.base_dir = os.path.abspath('cost_proc')
cost_proc.connect([
                    (infosource, datasource,[('subject_id', 'subject_id')]),
                    (datasource, cost,[('mesh_file','inputnode.mesh_file'),
                                               ('electrode_name_file','inputnode.electrode_name_file'),
                                               ('dipole_file','inputnode.dipole_file'),
                                               ('leadfield','inputnode.leadfield'),
                                               ])
                ])
cost_proc.connect([(cost, datasink, [("outputnode.mesh_file", "mesh_file")])])
cost_proc.connect([(cost, datasink, [("outputnode.cost_function_geo", "cost_function_geo")])])
cost_proc.connect([(infosource, datasink, [("subject_id", "subject_id")])])

if __name__ == '__main__':
    cost_proc.write_graph(graph2use="exec")
    cost_proc.run()
    #cost_proc.run(plugin='MultiProc', plugin_args={'n_procs' : 4})
