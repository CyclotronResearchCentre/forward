import numpy as np
import os.path as op
import subprocess
from forward.mesh import get_num_nodes_elements

def check_potential(potential):
    sane = False
    sane = True
    return sane

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
mesh_file = "/Users/erik/Dropbox/Analysis/TMSEEG/ForwardProblem/TMS007_running.msh"


'''
Define the conductivity tensor file to be incorportated. This should 
be the output of the diffusion preprocessing workflow.
'''
conductivity_tensor_file = "/Users/erik/Dropbox/Analysis/TMSEEG/ForwardProblem/dtifit__tensor_vreg_conductivity_elem.pos"

'''
Define the problem file
'''
problem_file = "../etc/eeg_forward.pro"


'''
Get the electrode names, remove markers and bad channels
'''
elec_names = np.loadtxt("ElectrodeNames.txt", dtype=str)
electrode_namelist = elec_names.tolist()

for marker in markers:
    try:
        electrode_namelist.remove(marker)
    except ValueError:
        continue

for bad_electrodes in bad_electrodes:
    try:
        electrode_namelist.remove(marker)
    except ValueError:
        continue
n_electrodes = len(electrode_namelist)
ground_idx = electrode_namelist.index(ground_electrode)
electrode_namelist.remove(ground_electrode)
src_electrodes = electrode_namelist

'''
The number M is the total non-ground recording sites
'''
M = len(src_electrodes)

'''
The number N is the total elements in the mesh file
'''
num_nodes, num_elements = get_num_nodes_elements(mesh_file)
N = num_elements

'''
Loop through the electrodes and solve the forward problem each time
'''

for src_idx, source_electrode in enumerate(src_electrodes):

    print("Index %d of %d " % (src_idx, len(src_electrodes)))
    print("Electrode %s" % source_electrode)

    '''
    Write the eeg_forward.dat file with the new source ID
    '''
    data_file = open('eeg_forward.dat', 'w')
    data_file.write("Source = %d; // %s \n" %
                    (5000 + src_idx, source_electrode))
    data_file.write("Sink = %d; // %s \n" %
                    (5000 + ground_idx, ground_electrode))

    '''
	Write the eeg_forward.dat file with the new source ID
	'''
    data_file.write('tensor_file = "%s";\n' % conductivity_tensor_file)
    data_file.close()

    #
    # Run the forward model for the selected source electrode:
    #
    subprocess.call(["getdp", problem_file, "-msh", mesh_file, "-bin",
                     "-solve", "Electrostatics"])

    '''
	Assess sanity check data to make sure potential at source / sink electrodes
	are 1 and -1 Volts, respectively.
	'''

    potential_file = 'v_electrodes.txt'
    potential = np.loadtxt(potential_file)
    if not check_potential(potential):
        raise Exception("Sanity check on potential failed")
        break

    '''
	Read electric field results file
	'''
    field_file = 'e_eeg_forward.txt'
    electric_field = np.loadtxt(field_file)

    '''
	Pre-allocated the Lead Field matrix, which is N rows x M columns
	'''
    if src_idx == 0:
        N = np.shape(electric_field)[0]*3
        L_e = np.empty((N, M))

    '''
    Write a single row in the lead field matrix (L_e)
    X,Y,Z elements of the field are flattened, so each row is
    e.g.:
    		[E1x E1y E1z E2x E2y E2z]

    Iteratively save over the file in MATLAB (.mat) and Numpy (.npy)
    file format.
    '''
    L_e[:, src_idx] = np.ravel(electric_field[:,-3:])
    #mdict = {}
    #mdict["L_e"] = L_e

    #ipdb.set_trace()
    #sio.savemat("leadfield.mat", mdict)
    #np.save("leadfield.npy", L_e)
    print("Saved over leadfield matrix")

ipdb.set_trace()

    #break
