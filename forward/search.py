import os.path as op
import scipy.optimize as so
import numpy as np
import h5py
from nipype.utils.filemanip import split_filename
from scipy.spatial.distance import cdist

def write_search_points(locations, geo_file):
    print('Writing search point .geo file to {f}'.format(f=geo_file))
    f = open(geo_file, "w")
    f.write('View "Simplex Search Points" {\n')
    for idx, location in enumerate(locations):
        x,y,z = location
        pt_str = (
            'SP(%5.2f, %5.2f, %5.2f){%5i};' % (x, y, z, idx))
        f.write(pt_str + '\n')
    f.write('};\n')
    f.write('View "Simplex Search Point Names" {\n')
    for idx, location in enumerate(locations):
        x,y,z = location
        pt_str = (
            'T3(%5.2f, %5.2f, %5.2f,0){"%s"};' % (x, y, z, "Point" + str(idx)))
        f.write(pt_str + '\n')
    f.write('};\n')
    f.close()
    return geo_file

def random_initial_guess(centroids):
    C = np.array(centroids)
    mx = np.max(C,0)
    mn = np.min(C,0)
    p_range = np.abs(mx) + np.abs(mn)
    p0 = np.random.rand(3)*p_range - np.abs(mx)
    bounds = [(mn[0],mx[0]), (mn[1],mx[1]), (mn[2],mx[2])]
    return p0, bounds

def calculate_element_cost(p, leadfield_matrix, V, centroids, element_ids):
    p = np.array([p])
    distances = cdist(p, centroids)
    mindist = np.min(distances)
    lf_idx = np.where(distances == mindist)[1]
    L = np.transpose(leadfield_matrix[lf_idx * 3:lf_idx * 3 + 3])        
    I = np.eye(len(V))
    Linv = np.linalg.pinv(L)
    val = np.dot(V.T, I - np.dot(L,Linv))
    R = np.dot(val,V)
    #cost = np.linalg.norm(R)
    cost = np.sqrt(R / np.linalg.norm(V))
    if np.isnan(cost):
        cost = 0
    print(p,cost)
    return cost

def single_dipole_search(electrode_potential, leadfield, mesh_file, elements_to_consider=[1002]):
    from forward.mesh import read_mesh, get_closest_element_to_point
    from nipype import logging
    iflogger = logging.getLogger('interface')

    geo_file = op.abspath("search_points.geo")
    data_name = "leadfield"
    lf_data_file = h5py.File(leadfield, "r")
    lf_data = lf_data_file.get(data_name)
    leadfield_matrix = lf_data.value
    _, lf_name, _ = split_filename(leadfield)

    bounds = (0,leadfield_matrix.shape[0]/3)
    mesh_data = read_mesh(mesh_file, elements_to_consider)

    centroids = []
    element_ids = []
    for poly in mesh_data:
        element_ids.append(poly["element_id"])
        centroids.append(poly["centroid"])
    

    p0, bounds = random_initial_guess(centroids)
    #p0 = np.array([0,0,60])

    mw = {}
    args = (leadfield_matrix, electrode_potential, centroids, element_ids)
    mw['args'] = args
    #mw['bounds'] = bounds

    #res = so.basinhopping(calculate_element_cost, p0, T=0.01, stepsize = 3, minimizer_kwargs=mw, disp=True, niter_success=20)

    # Perform downhill search, minimizing cost function
    xopt, fopt, n_iter, funcalls, _, allvecs = so.fmin(calculate_element_cost, p0, ftol=0.00000001, args=(leadfield_matrix, electrode_potential, centroids, element_ids), xtol = 1, disp=True, full_output=True, retall=True)

    write_search_points(allvecs, geo_file) 

    # Need to minimize cost by changing values of orientation and magnitude:
    _, element_idx, centroid, element_data, lf_idx = get_closest_element_to_point(mesh_file, elements_to_consider, np.array([xopt]))

    L = np.transpose(leadfield_matrix[lf_idx * 3:lf_idx * 3 + 3])
    V = electrode_potential
    qopt, _, _, _ = np.linalg.lstsq(L, V)
    return xopt, qopt, element_data, geo_file