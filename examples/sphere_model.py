import os
import os.path as op
import nipype.interfaces.io as nio           # Data i/o
import nipype.interfaces.utility as util     # utility
import nipype.pipeline.engine as pe          # pypeline engine
from forward.leadfield import create_forward_model_workflow
from forward.sphere import create_4_shell_model


'''
For the sphere we have no markers
'''

markers = []
ground_electrode = 'vertex001'
radii=[85, 88, 92, 100]

electrode_location_file = op.join(os.environ["FWD_DIR"], "etc", "icosahedron42.txt")

'''
Create the nodes and forward modelling workflow
'''
create_sphere_model_interface = util.Function(input_names=["electrode_location_file", "radii", "out_file"],
                                       output_names=["mesh_file", "electrode_name_file", "electrode_location_file"], function=create_4_shell_model)

create_sphere_model = pe.Node(interface=create_sphere_model_interface, name="create_sphere_model")
create_sphere_model.inputs.radii = radii
create_sphere_model.inputs.electrode_location_file = electrode_location_file
create_sphere_model.inputs.out_file = "4shell.msh"

fwd = create_forward_model_workflow("sphere", conductivity_tensor_included=False)
fwd.inputs.inputnode.marker_list = markers
fwd.inputs.inputnode.ground_electrode = ground_electrode
fwd.inputs.inputnode.conductivity_tensor_included = False
fwd.inputs.write_problem_file.orig_pro_file = op.join(os.environ["FWD_DIR"], "etc/eeg_forward_sphere.pro")
#fwd.inputs.run_forward_model.binary_output_files = False #for debugging
fwd.inputs.run_forward_model.binary_output_files = True

datasink = pe.Node(interface=nio.DataSink(),
                   name="datasink")
datasink.inputs.base_directory = op.abspath('sphere_datasink')

sphere_proc = pe.Workflow(name="sphere_proc")
sphere_proc.base_dir = os.path.abspath('sphere_proc')
sphere_proc.connect([
                    (create_sphere_model,fwd,[('electrode_name_file', 'inputnode.electrode_name_file')]),
                    (create_sphere_model,fwd,[('mesh_file', 'inputnode.mesh_file')]),
                ])
sphere_proc.connect([(create_sphere_model, datasink, [("mesh_file", "mesh_file")])])
sphere_proc.connect([(create_sphere_model, datasink, [("electrode_name_file", "electrode_name_file")])])
sphere_proc.connect([(create_sphere_model, datasink, [("electrode_location_file", "electrode_location_file")])])
sphere_proc.connect([(fwd, datasink, [("outputnode.leadfield", "leadfield")])])

if __name__ == '__main__':
    sphere_proc.write_graph(graph2use="exec")
    #sphere_proc.run()
    sphere_proc.run(plugin='MultiProc', plugin_args={'n_procs' : 4})
