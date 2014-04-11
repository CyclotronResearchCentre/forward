import os
import os.path as op
import numpy as np
from forward.electrodes import write_electrode_labels

fwd_dir = os.environ['FWD_DIR']

electrode_location_file = op.join(fwd_dir, 'etc', 'ElectrodeLocations_BESA256.txt')
electrode_name_file = op.join(fwd_dir, 'etc', 'ElectrodeNames_BESA256.txt')

points = np.loadtxt(electrode_location_file, delimiter=',')
electrode_names = np.loadtxt(electrode_name_file, dtype='str')

# This is a good time to move the electrodes around
# and get them to match your skull. e.g.:
# Scale by 10x
#points = points * 10
# Move values 5 mm to the right, 25 mm backward, 
# and 20 mm upward
#points[:,0] = points[:,0] + 5
#points[:,1] = points[:,1] - 25
#points[:,2] = points[:,2] + 20

write_electrode_labels(points, electrode_names, out_file="BESA256.geo")