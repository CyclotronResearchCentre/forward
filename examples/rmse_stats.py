import h5py
import numpy as np

# Get all the values
x = h5py.File("TMS007_gmsh_cond_elec_rmse_x.hdf5","r")
h5_x = x.get("rmse_x")
rmse_x = h5_x.value

y = h5py.File("TMS007_gmsh_cond_elec_rmse_y.hdf5","r")
h5_y = y.get("rmse_y")
rmse_y = h5_y.value

z = h5py.File("TMS007_gmsh_cond_elec_rmse_z.hdf5","r")
h5_z = z.get("rmse_z")
rmse_z = h5_z.value

avg = h5py.File("TMS007_gmsh_cond_elec_rmse_avg.hdf5","r")
h5_avg = avg.get("rmse_avg")
rmse_avg = h5_avg.value

# Write some basic values:
print("RMSE-Avg: %3.3f (%3.3f)" % (np.mean(rmse_avg), np.std(rmse_avg)))
print("RMSE-X: %3.3f (%3.3f)" % (np.mean(rmse_x), np.std(rmse_x)))
print("RMSE-Y: %3.3f (%3.3f)" % (np.mean(rmse_y), np.std(rmse_y)))
print("RMSE-Z: %3.3f (%3.3f)\n" % (np.mean(rmse_z), np.std(rmse_z)))
