import os.path as op
from forward.meshmath import subtract_meshes

mesh1 = op.abspath("cost_proc/TMS007_gmsh_cond_elec_cost_ANISO.msh")
mesh2 = op.abspath("cost_proc/TMS007_gmsh_cond_elec_cost_ISO.msh")
view_name1 = "Cost function /Users/erik/Dropbox/Analysis/TMSEEG/forward/examples/LeadfieldSample/TMS007_aniso_leadfield.hdf5"
view_name2 = "Cost function /Users/erik/Dropbox/Analysis/TMSEEG/forward/examples/LeadfieldSample/TMS007_iso_leadfield.hdf5"

subtract_meshes(mesh1, mesh2, view_name1, view_name2, out_file="MeshMath.msh")