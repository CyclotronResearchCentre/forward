import os
import os.path as op
import nipype.interfaces.io as nio           # Data i/o
import nipype.pipeline.engine as pe          # pypeline engine
from forward.dipole import create_cost_function_workflow


'''
For the sphere we have no markers
'''

markers = []
ground_electrode = 'vertex001'
radii=[85, 88, 92, 100]

electrode_location_file = op.join(os.environ["FWD_DIR"], "etc", "icosahedron42.txt")
dipole_location_file = op.join(os.environ["FWD_DIR"], "etc", "sphere_dipole.dip")
leadfield_file = op.join(os.environ["FWD_DIR"], "examples", "sphere_datasink", "leadfield", "leadfield.hdf5")
mesh_file = op.join(os.environ["FWD_DIR"], "examples", "sphere_datasink", "mesh_file", "4shell_elec.msh")
electrode_name_file = op.join(os.environ["FWD_DIR"], "examples", "sphere_datasink", "electrode_name_file", "icosahedron42.txt_names.txt")

dipole_row = 0

'''
Create the nodes and forward modelling workflow
'''

cost = create_cost_function_workflow("sphere")
cost.inputs.inputnode.dipole_file = dipole_location_file
cost.inputs.inputnode.leadfield = leadfield_file
cost.inputs.inputnode.dipole_row = dipole_row
cost.inputs.inputnode.mesh_id = 1001
cost.inputs.inputnode.mesh_file = mesh_file

datasink = pe.Node(interface=nio.DataSink(),
                   name="datasink")
datasink.inputs.base_directory = op.abspath('sphere_cost_datasink')

sphere_proc = pe.Workflow(name="sphere_cost_proc")
sphere_proc.base_dir = os.path.abspath('sphere_cost_proc')
sphere_proc.connect([(cost, datasink, [("outputnode.mesh_file", "cost_mesh_file")])])
sphere_proc.connect([(cost, datasink, [("outputnode.potential_from_dipole", "potential_from_dipole")])])

if __name__ == '__main__':
    sphere_proc.write_graph(graph2use="exec")
    sphere_proc.run()
    #sphere_proc.run(plugin='MultiProc', plugin_args={'n_procs' : 4})
