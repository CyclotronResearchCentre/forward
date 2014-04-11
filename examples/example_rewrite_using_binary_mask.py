import os
import os.path as op
import glob
from forward.tdcs import rewrite_mesh_from_binary_mask

from forward.datasets import sample
data_path = sample.data_path()

subject_id = "TMS007"
mesh_file = glob.glob('structural_datasink/subject/volume_mesh/*/*.msh')[0]
mask_file = op.join(data_path, subject_id, "DLPFC_mask.nii.gz")
mesh_id = 1002

# Binary mask should be drawn on the conformed T1 from FreeSurfer

print mask_file
new_mesh = rewrite_mesh_from_binary_mask(mask_file, mesh_file, mesh_id)


