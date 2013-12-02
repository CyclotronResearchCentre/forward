import os
import os.path as op
import h5py
import numpy as np
from numpy.linalg import norm
import ipdb
from forward.mesh import get_closest_element_to_point


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
probe_dipole = np.array([[0, 0, 70]])
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
mesh_data, closest_element_idx, centroid, closest_element_data, lf_idx = get_closest_element_to_point(mesh_file, mesh_id, probe_dipole)
lf_getdp = np.transpose(leadfield_matrix[lf_idx * 3:lf_idx * 3 + 3])

# Re-reference to ground electrode
#lf_sphere = lf_sphere - lf_sphere[ground_idx]
print(lf_sphere[ground_idx])
lf_sphere = lf_sphere - lf_sphere[ground_idx]
lf_sphere = np.delete(lf_sphere, (ground_idx), axis=0)

lf_getdp = lf_getdp * np.max(lf_sphere,0)/np.max(lf_getdp,0)

rdms = np.empty((np.shape(lf_getdp)[1], 1))
rows = np.shape(lf_getdp)[1]
for idx in xrange(0, rows):
    rdms[idx, :] = norm(lf_getdp[:, idx] / norm(lf_getdp[:, idx]) -
                        lf_sphere[:, idx] / norm(lf_sphere[:, idx]))

mags = np.divide(np.sqrt(np.sum(lf_getdp ** 2, 0)),
                 np.sqrt(np.sum(lf_sphere ** 2, 0)) )

print("RDMs %s" % str(rdms))
print("Mean RDMs %s" % str(np.mean(rdms)))
print("MAGs %s" % str(mags))
ipdb.set_trace()

from matplotlib.pyplot import plot, show, legend
x_sph = plot(lf_sphere[:,0], label="X analytical")
x_getdp = plot(lf_getdp[:,0], label="X GetDP")
legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
       ncol=2, mode="expand", borderaxespad=0.)
show()

y_sph = plot(lf_sphere[:,1], label="Y analytical")
y_getdp = plot(lf_getdp[:,1], label="Y GetDP")
legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
       ncol=2, mode="expand", borderaxespad=0.)
show()

z_sph = plot(lf_sphere[:,2], label="Z analytical")
z_getdp = plot(lf_getdp[:,2], label="Z GetDP")
legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
       ncol=2, mode="expand", borderaxespad=0.)
show()

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
