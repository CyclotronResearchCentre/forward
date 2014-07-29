import os
import os.path as op
import numpy as np
from forward.electrodes import write_electrode_labels, rewrite_mesh_with_electrodes
from nipype import logging
iflogger = logging.getLogger('interface')


def add_electrodes(mesh_filename):
    fwd_dir = os.environ['FWD_DIR']

    electrode_location_file = op.join(fwd_dir, 'etc', 'ElectrodeLocations_TMS007.txt')
    electrode_name_file = op.join(fwd_dir, 'etc', 'ElectrodeNames_BESA256.txt')

    points = np.loadtxt(electrode_location_file, delimiter=',')
    electrode_names = np.loadtxt(electrode_name_file, dtype='str')

    write_electrode_labels(points, electrode_names, out_file="BESA256.geo")
    mesh_id = 1005

    out_file = rewrite_mesh_with_electrodes(electrode_location_file, electrode_name_file, mesh_filename, mesh_id)
    return out_file

def define_tDCS_optimization(mesh_filename, region_ID = 1015, maximize = True, move_anode=True, move_cathode=False):
    # Optimization function
    # Specify location limits for anode and cathode?
    # - e.g. a set of points on the scalp?
    # - Select anode and cathode
    # - Run simulation (brute force?)
    # - Extract current density in specified region using ID

    # Return best placement for anode and cathode
    return

# Specify region to optimize current for
# - e.g. binary mask of the region
# - rewrite mesh with new region ID
# - specify that new region ID has same conductivity as old region ID

def rewrite_mesh_from_binary_mask(mask_file, mesh_file, mesh_id, new_phys_id=1015):
    # Takes about 8 minutes
    import numpy as np
    import nibabel as nb
    import os.path as op
    from forward.mesh import read_mesh
    from nipype import logging

    iflogger = logging.getLogger('interface')

    mask_img = nb.load(mask_file)
    mask_data = mask_img.get_data()
    mask_affine = mask_img.get_affine()
    mask_header = mask_img.get_header()
    #tensor_data = np.flipud(tensor_image.get_data())

    # Define various constants
    elements_to_consider = [1002] #Use only white matter
    vx, vy, vz = mask_header.get_zooms()[0:3]
    max_x, max_y, max_z = np.shape(mask_data)[0:3]
    halfx, halfy, halfz = np.array((vx*max_x, vy*max_y, vz*max_z))/2.0

    mesh_data, _, _, _ = read_mesh(mesh_file, elements_to_consider)

    elem_list = []
    for idx, poly in enumerate(mesh_data):
        i = np.round((poly['centroid'][0]+halfx)/vx).astype(int)
        j = np.round((poly['centroid'][1]+halfy)/vy).astype(int)
        k = np.round((poly['centroid'][2]+halfz)/vz).astype(int)
        if mask_data[i,j,k]:
            elem_list.append(poly['element_id'])
        #iflogger.info("%3.3f%%" % (float(idx)/num_polygons*100.0))
    new_elem_id = new_phys_id+2000
    out_file = swap_element_ids(mesh_file, elem_list, new_elem_id, new_phys_id)
    print('Writing binary masked mesh file to %s' % out_file)
    return out_file

def swap_element_ids(mesh_filename, elem_list, new_elem_id, new_phys_id):
    import time
    import os.path as op
    from nipype.utils.filemanip import split_filename
    from nipype import logging
    iflogger = logging.getLogger('interface')

    iflogger.info("Reading mesh file: %s" % mesh_filename)

    start_time = time.time()
    mesh_file = open(mesh_filename, 'r')
    _, name, _ = split_filename(mesh_filename)
    out_file = op.abspath(name + "_mask.msh")

    f = open(out_file, 'w')
    while True:
        line = mesh_file.readline()
        f.write(line)

        if line == '$Nodes\n':
            line = mesh_file.readline()
            f.write(line)
            number_of_nodes = int(line)
            iflogger.info("%d nodes in mesh" % number_of_nodes)

            for i in xrange(0, number_of_nodes):
                line = mesh_file.readline()
                f.write(line)

        elif line == '$Elements\n':
            line = mesh_file.readline()
            f.write(line)
            number_of_elements = int(line)
            iflogger.info("%d elements in mesh" % number_of_elements)
            elem_lines = []
            for i in xrange(0, number_of_elements):
                # -- If all elements were quads, each line has 10 numbers.
                #    The first one is a tag, and the last 4 are node numbers.
                line = mesh_file.readline()
                elem_lines.append(line)
                elem_data = line.split()
                if int(elem_data[0]) in elem_list:
                    elem_data[3] = str(new_phys_id)
                    elem_data[4] = str(new_elem_id)
                line = " ".join(elem_data) + "\n"
                f.write(line)

        elif line == '$EndElementData\n' or len(line) == 0:
            break

    mesh_file.close()
    f.close()
    elapsed_time = time.time() - start_time
    print(elapsed_time)
    return out_file
