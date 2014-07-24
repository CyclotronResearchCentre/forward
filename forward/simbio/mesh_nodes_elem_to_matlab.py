import numpy as np
from forward.mesh import read_mesh
n_electrodes = 42
elem_to_consider = [1000, 1002, 1003, 1004]
elem_to_consider.extend(range(5000, 5000+n_electrodes))
mesh_data, nodes, elem, labels = read_mesh('4shell_elec.msh', elem_to_consider)
np.savetxt('nodes.txt', nodes)
np.savetxt('elem.txt', elem)
np.savetxt('labels.txt', labels)