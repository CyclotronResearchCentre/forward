import os
import os.path as op

from forward.sphere import create_4_shell_model
from forward.leadfield import create_forward_model_workflow

mesh_file, electrode_name_file = create_4_shell_model(radii=[85, 88, 92, 100], out_file='4shell.msh')

fwd = create_forward_model_workflow("sphere", conductivity_tensor_included=False)
fwd.base_dir = "sphere_proc"
fwd.inputs.inputnode.marker_list = []
fwd.inputs.inputnode.ground_electrode = "vertex001"
fwd.inputs.inputnode.conductivity_tensor_included = False
fwd.inputs.inputnode.electrode_name_file = electrode_name_file 
fwd.inputs.inputnode.mesh_file = mesh_file
fwd.inputs.write_problem_file.orig_pro_file = op.join(os.environ["FWD_DIR"], "etc/eeg_forward_sphere.pro")
fwd.write_graph(graph2use="exec")
fwd.run()