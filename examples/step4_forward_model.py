'''
The fourth step is to run the actual forward modelling and generate
a lead field matrix. This is done with the concept of reciprocity.
'''

from forward.leadfield import create_forward_model_workflow


'''
Define the markers, ground electrode, and any bad electrodes to ignore
'''

markers = ['LeftEar', 'RightEar', 'NZ']
ground_electrode = 'IZ'
bad_electrodes = []

'''
Define the mesh file to be used. This should be the volume mesh 
output of the structural preprocessing workflow.
'''
mesh_file = "../etc/TMS007_running.msh"
conductivity_tensor_included = True

'''
Define the problem file
'''
problem_file = "../etc/eeg_forward.pro"


'''
Get the electrode names, remove markers and bad channels
'''
elec_names = np.loadtxt("ElectrodeNames.txt", dtype=str)
electrode_namelist = elec_names.tolist()