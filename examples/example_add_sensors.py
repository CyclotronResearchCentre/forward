import os
import os.path as op
import glob
from forward.tdcs import add_electrodes

subject_id = 'TMS007'
fwd_dir = os.environ['FWD_DIR']
path = op.join(fwd_dir,'examples','%s*_mask.msh' % subject_id)
mesh_filename = glob.glob(path)[0]

add_electrodes(mesh_filename)