import os.path as op
from forward.datasets import sample
data_path = sample.data_path(name="leadfield")

from forward.leadfield import compare_leadfields

lf1 = op.join(data_path, "TMS007_iso_leadfield.hdf5")
lf2 = op.join(data_path, "TMS007_aniso_leadfield.hdf5")
mesh_file = op.join(data_path, "TMS007_gmsh_cond_elec.msh")
compare_leadfields(lf2, lf1, mesh_file)