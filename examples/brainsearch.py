import os
import os.path as op
from forward.search import single_dipole_search
import numpy as np

from forward.datasets import sample
data_path = sample.data_path(name="leadfield")


conductivity_tensor_included = False
potential_file = op.join(os.environ["FWD_DIR"], "examples", 'forward_datasink','subject','potential_from_dipole','_subject_id_TMS007','potential.npy')

if conductivity_tensor_included:
    leadfield = op.join(data_path, "TMS007_aniso_leadfield.hdf5")
else:
    leadfield = op.join(data_path, "TMS007_iso_leadfield.hdf5")

mesh_file = op.join(os.environ["FWD_DIR"], "examples", "sphere_cost_datasink", "cost_mesh_file", "4shell_elec_dipole_cost.msh")

elements_to_consider = [1002]
electrode_potential = np.load(potential_file)
element_idx, geo_file = single_dipole_search(electrode_potential, leadfield, mesh_file, elements_to_consider)