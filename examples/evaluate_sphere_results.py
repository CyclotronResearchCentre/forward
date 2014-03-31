import os
import os.path as op
import h5py
import numpy as np
from numpy.linalg import norm
import ipdb
import nipype.pipeline.engine as pe          # pypeline engine
from forward.mesh import get_closest_element_to_point
from forward.analytical import FourShellAnalyticalModel
from matplotlib.pyplot import plot, show, legend
import matplotlib.pyplot as plt

lf_filename = "leadfield.hdf5"
leadfield_file = op.join(
    op.curdir, "sphere_datasink", "leadfield", lf_filename)

data_name = "leadfield"

leadfield_hdf5 = h5py.File(leadfield_file, "r")
leadfield_data = leadfield_hdf5.get(data_name)
leadfield_matrix = leadfield_data.value

filename = "lf_4shell_analytical_42.txt"
analytical_solution = op.join(os.environ["FWD_DIR"], "etc", filename)
lf_sphere = np.loadtxt(analytical_solution, delimiter=",")

# Set the position of the probe dipole
mesh_filename = "4shell_elec.msh"
elec_filename = "icosahedron42.txt_names.txt"
elec_loc_fname = "electrode_locations.txt"
ground_electrode = 'vertex001'
elec_names = np.loadtxt(elec_filename, dtype=str)
elec_names = elec_names.tolist()
ground_idx = elec_names.index(ground_electrode)

mesh_file = op.join(
    op.curdir, "sphere_datasink", "mesh_file", mesh_filename)
electrode_name_file = op.join(
    op.curdir, "sphere_datasink", "electrode_name_file", elec_filename)
electrode_location_file = op.join(
    op.curdir, "sphere_datasink", "electrode_location_file", elec_loc_fname)
mesh_id = 1001


def get_rdm_and_mags(mesh_file, mesh_id, probe_dipole, radii, cond, n=42):
    four = pe.Node(interface=FourShellAnalyticalModel(), name="four")
    name = str(probe_dipole[0]).replace(" ","")
    name = name.replace('[','')
    name = name.replace(']','')
    out_directory = op.abspath("dipole" + name)
    if not op.exists(out_directory):
        os.makedirs(out_directory)
    four.base_dir = out_directory
    four.inputs.probe_dipole = probe_dipole[0].tolist()
    four.inputs.script = "pyscript.m"
    four.inputs.sphere_radii = radii
    four.inputs.shell_conductivity = cond
    four.inputs.icosahedron_sides = n
    four.inputs.fieldtrip_path = "/Developer/fieldtrip"
    four.run()
    analytical_solution = op.join(out_directory, "four", "analytical.txt")
    # Re-reference to ground electrode
    lf_sphere = np.loadtxt(analytical_solution, delimiter=",")
    lf_sphere = lf_sphere - lf_sphere[ground_idx]
    lf_sphere = np.delete(lf_sphere, (ground_idx), axis=0)


    mesh_data, closest_element_idx, centroid, closest_element_data, lf_idx = get_closest_element_to_point(mesh_file, mesh_id, probe_dipole)
    lf_getdp = np.transpose(leadfield_matrix[lf_idx * 3:lf_idx * 3 + 3])

    openmeeg_solution = op.join(out_directory, "four", "openmeeg.txt")
    lf_openmeeg = np.loadtxt(openmeeg_solution, delimiter=",")
    
    try:
        lf_openmeeg = lf_openmeeg - lf_openmeeg[ground_idx]
        lf_openmeeg = np.delete(lf_openmeeg, (ground_idx), axis=0)

        rdm_openmeeg = norm( lf_openmeeg / norm(lf_openmeeg) - lf_sphere / norm(lf_sphere) )
        mag_openmeeg = np.divide(norm(lf_openmeeg),norm(lf_sphere))
    except IndexError:
        rdm_openmeeg = None
        mag_openmeeg = None

    rdms = np.empty((np.shape(lf_getdp)[1], 1))
    rows = np.shape(lf_getdp)[1]
    for idx in xrange(0, rows):
        rdms[idx, :] = norm(lf_getdp[:, idx] / norm(lf_getdp[:, idx]) -
                            lf_sphere[:, idx] / norm(lf_sphere[:, idx]))

    rdm_singlevalue = norm( lf_getdp / norm(lf_getdp) - lf_sphere / norm(lf_sphere) )
    

    mags = np.divide(np.sqrt(np.sum(lf_getdp ** 2, 0)),
                     np.sqrt(np.sum(lf_sphere ** 2, 0)) )
    mag_singlevalue = np.divide(norm(lf_getdp),norm(lf_sphere))

    print("RDMs %s" % str(rdms))
    print rdm_singlevalue
    print("MAGs %s" % str(mags))
    print mag_singlevalue
    return (rdms, mags, rdm_singlevalue, mag_singlevalue, rdm_openmeeg,
        mag_openmeeg, lf_getdp, lf_sphere)

results = []
radii=[85,88,92,100]
#radii = (np.array(radii)/2).tolist()

# Realistic conductivity values in S/m
cond=[0.33, 1.79, 0.0042, 0.33]
#cond=[1, 1/20, 1/80, 1]

n=42
#loop = False
loop = True

probe_dipole = np.array([[0, 0, 70]])
(rdms, mags, rdm_sv, mag_sv,
    rdm_OMEEG, mag_OMEEG, lf_getdp, lf_sphere) = get_rdm_and_mags(mesh_file, mesh_id, probe_dipole, radii, cond, n)


if loop:
    print("Looping through dipole positions:")
    distances = [10, 20, 40, 60, 70, 75, 80, 82, 84, 85]
    #distances = (np.array(distances)/2).tolist()

    distances.reverse()
    norms = []
    print(distances)
    for idx in distances:
        print(idx)
        probe_dipole = np.array([[0, 0, idx]])
        rdms, mags, rdm_sv, mag_sv, rdm_OMEEG, mag_OMEEG, lf_getdp, lf_sphere = get_rdm_and_mags(mesh_file, mesh_id, probe_dipole, radii, cond, n)
        results.append((idx, rdm_sv, mag_sv, rdm_OMEEG, mag_OMEEG))

    res = np.array(results)
    max_r = np.max(radii)

    np.save("SphereResults.npy", res)
    # Plot RDMs
    plot(res[:,0], res[:,1], '.', label="GetDP")
    plot(res[:,0], res[:,3], '+', label="OpenMEEG")
    plt.xlabel('Dipole position (from 0 in Z dir, max radius = %3.0f)' % max_r, fontsize=14)
    plt.ylabel('Relative difference measure (RDM)', fontsize=14)
    for r in radii:
        plt.axvline(x=r,ymin=0,ymax=1,color="r")
    legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
            ncol=2, mode="expand", borderaxespad=0.)
    show()
    plt.axis([0.0, 0.6, 0, 100])
    ax = plt.gca()
    ax.set_autoscale_on(False)


    # Plot MAGs
    plot(res[:,0], res[:,2], '.', label="GetDP")
    plot(res[:,0], res[:,4], '+', label="OpenMEEG")
    plt.xlabel('Dipole position (from 0 in Z, max radius = %3.0f)' % max_r, fontsize=14)
    plt.ylabel('Magnification errors (MAG)', fontsize=14)
    for r in radii:
        plt.axvline(x=r,ymin=0,ymax=1,color="r")
    legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
            ncol=2, mode="expand", borderaxespad=0.)
    show()
    plt.axis([0.0, 6, 0, 1])
    ax = plt.gca()
    ax.set_autoscale_on(False)

def write_gmsh_pos_file(scalars, electrode_location_file, ground_idx, scalar_name="RDM"):
    out_file = scalar_name + ".pos"
    points = np.loadtxt(electrode_location_file, delimiter=",")
    points = np.delete(points, (ground_idx), axis=0)
    f = open(op.abspath(out_file), "w")
    f.write('View "%s" {\n' % scalar_name)
    for idx, scalar in enumerate(scalars):
        towrite = 'SP(%3.3f, %3.3f, %3.3f){%3.3f};' % (
            points[idx][0], points[idx][1], points[idx][2], scalar)
        f.write(towrite + "\n")
    f.write("};\n")
    f.close()
    return out_file

#write_gmsh_pos_file(rdms[:, 0].tolist(), electrode_location_file, ground_idx)


## For testing purposes
# x_sph = plot(lf_sphere[:,0], label="X analytical")
# x_getdp = plot(lf_getdp[:,0], label="X GetDP")
# legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
#        ncol=2, mode="expand", borderaxespad=0.)
# show()

# y_sph = plot(lf_sphere[:,1], label="Y analytical")
# y_getdp = plot(lf_getdp[:,1], label="Y GetDP")
# legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
#        ncol=2, mode="expand", borderaxespad=0.)
# show()

# z_sph = plot(lf_sphere[:,2], label="Z analytical")
# z_getdp = plot(lf_getdp[:,2], label="Z GetDP")
# legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
#        ncol=2, mode="expand", borderaxespad=0.)
# show()