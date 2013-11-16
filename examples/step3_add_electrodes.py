'''
The third step is to rewrite elements of the volume mesh with new element IDs.
These element IDs will be used to create the source and sink for the forward
modelling step. Elements are identified by their proximity the location of the
EEG electrodes provided.
'''

from forward.electrodes import get_scalp_tris

in_file = "/Users/erik/Dropbox/Analysis/TMSEEG/TMS007/TMS/SESS_CHAR_TMS007.mat"
electrode_name_file = "/Users/erik/Dropbox/Analysis/TMSEEG/ElectrodeNames.txt"

original_T1 = "/Users/erik/Dropbox/Analysis/TMSEEG/TMS007/20121204_134045t1mpragetms20121113s002a001.nii"
final_volume = "/Users/erik/Dropbox/Analysis/TMSEEG/TMS007/orig_out_flirt.nii"
#mesh_file = "/Users/erik/Dropbox/Analysis/TMSEEG/ForwardProblem/skinmsh.msh"
mesh_file = "/Users/erik/Dropbox/Analysis/TMSEEG/ForwardProblem/TMS007_gmsh.msh"
mesh_id = 1005  # Scalp

get_scalp_tris(in_file, original_T1, final_volume, mesh_file, mesh_id, electrode_name_file)