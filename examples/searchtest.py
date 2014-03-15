import os
import os.path as op
from forward.search import single_dipole_search
import numpy as np
import time

potential_file = op.join(os.environ["FWD_DIR"], "examples", "sphere_cost_datasink", "potential_from_dipole", "potential.npy")
leadfield_file = op.join(os.environ["FWD_DIR"], "examples", "sphere_datasink", "leadfield", "leadfield.hdf5")
mesh_file = op.join(os.environ["FWD_DIR"], "examples", "sphere_datasink", "mesh_file", "4shell_elec.msh")

elements_to_consider = [1001]
electrode_potential = np.load(potential_file)
start_time = time.time()
xopt, qopt, element_data, geo_file = single_dipole_search(electrode_potential, leadfield_file, mesh_file, elements_to_consider)
end_time = time.time()

print("Finished in %f seconds" % (end_time - start_time))
print(xopt)
print(qopt)
print(geo_file)