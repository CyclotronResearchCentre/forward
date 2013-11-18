'''
The third step is to rewrite elements of the volume mesh with new element IDs.
These element IDs will be used to create the source and sink for the forward
modelling step. Elements are identified by their proximity the location of the
EEG electrodes provided.
'''

from forward.electrodes import rewrite_mesh_with_electrodes
import os.path as op
import glob

from forward.datasets import sample
data_path = sample.data_path()

subject_id = 'TMS007'

'''
The files required can be found in the ForwardSample data package,
and as the outputs of the diffusion pipeline (or structural)
'''

electrode_position_file = op.join(data_path, subject_id, "ElectrodePositions.txt")
electrode_name_file = op.join(data_path, subject_id, "ElectrodeNames.txt")

# If you want to use the structural mesh, rather than the mesh including the conductivity
# tensors, you should uncomment the following line, and comment the one below it.
#mesh_path = op.abspath("structural_datasink/subject/volume_mesh/")
mesh_path = op.abspath("diffusion_datasink/subject/mesh_file/")

assert(op.exists(mesh_path))
search_string = "/*%s/%s*.msh" % (subject_id, subject_id)
try:
    mesh_file = glob.glob(mesh_path + search_string)[0]
except IndexError:
    raise IOError("Mesh file not found. You must run step1 and/or step2 \
        before trying to add electrodes to the mesh.")

'''
We use the element ID for the scalp volume in the Gmsh mesh
'''
mesh_id = 5

out_file = rewrite_mesh_with_electrodes(electrode_position_file, electrode_name_file,
    mesh_file, mesh_id)