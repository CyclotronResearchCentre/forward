import os
import os.path as op
import numpy as np
from numpy.linalg import norm
import ipdb
import nipype.pipeline.engine as pe          # pypeline engine
from forward.analytical import FourShellAnalyticalModel
from matplotlib.pyplot import plot, show, legend
import matplotlib.pyplot as plt

def rereference(array, ground_idx):
    array = array - array[ground_idx]
    array = np.delete(array, (ground_idx), axis=0)
    return array

def calculate_rdm_and_mags(lf_test, radii, cond, probe_dipole=[0,0,70], n=42, ground_idx=0):
    name = "".join([str(x) for x in probe_dipole])
    fwd_dir = os.environ["FWD_DIR"]
    out_directory = op.join(fwd_dir, "examples", "dipole" + name)
    if not op.exists(out_directory):
        os.makedirs(out_directory)

    four = pe.Node(interface=FourShellAnalyticalModel(), name="four")
    four.base_dir = out_directory
    four.inputs.probe_dipole = probe_dipole
    four.inputs.script = "pyscript.m"
    four.inputs.sphere_radii = radii
    four.inputs.shell_conductivity = cond
    four.inputs.icosahedron_sides = n
    four.inputs.fieldtrip_path = "/Developer/fieldtrip"
    analytical_solution = op.join(out_directory, "four", "analytical.txt")
    if not op.exists(analytical_solution):
        four.run()
    
    # Re-reference to ground electrode
    lf_sphere = np.loadtxt(analytical_solution, delimiter=",")
    lf_sphere = rereference(lf_sphere, ground_idx)
    lf_sphere = lf_sphere[:,2]

    openmeeg_solution = op.join(out_directory, "four", "openmeeg.txt")
    if op.exists(openmeeg_solution):
        lf_openmeeg = np.loadtxt(openmeeg_solution, delimiter=",")
        lf_openmeeg = rereference(lf_openmeeg, ground_idx)
        lf_openmeeg = lf_openmeeg[:,2]

        rdm_OMEEG = norm( lf_openmeeg / norm(lf_openmeeg) - lf_sphere / norm(lf_sphere) )
        mag_OMEEG = np.divide(norm(lf_openmeeg), norm(lf_sphere))

    lf_test = rereference(lf_test, ground_idx)
    #ipdb.set_trace()
    #lf_sphere = lf_sphere[:,2]
    #assert(lf_test.shape == lf_sphere.shape)
    rdm_test = norm( lf_test / norm(lf_test) - lf_sphere / norm(lf_sphere) )
    mag_test = np.divide(norm(lf_test),norm(lf_sphere))

    return rdm_test, mag_test, rdm_OMEEG, mag_OMEEG


import scipy.io as sio
data = sio.loadmat("SimbioResults.mat")
data = data['d']

# To fix units
data = data / 1000 

radii=[85,88,92,100]

# Realistic conductivity values in S/m
cond=[0.33, 1.79, 0.0042, 0.33]

n=42
#ipdb.set_trace()
results = []
for idx, lf_test in enumerate(data):
    rdm_test, mag_test, rdm_OMEEG, mag_OMEEG = calculate_rdm_and_mags(
        lf_test, radii, cond, probe_dipole=[0,0,idx], n=42, ground_idx=0)

    print("%3.2f %3.2f %3.2f %3.2f\n" % (rdm_test, mag_test, rdm_OMEEG, mag_OMEEG))
    results.append([rdm_test, mag_test, rdm_OMEEG, mag_OMEEG])

np.save("SimBioRDM_MAG.npy", results)