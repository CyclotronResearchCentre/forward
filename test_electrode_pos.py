def get_electrode_positions(in_file, original_T1, final_volume, electrode_name_file, out_basename="electrodes_", session_number=5):
    import scipy.io as sio
    import os
    import os.path as op
    import numpy as np
    import nibabel as nb
    struct = sio.loadmat(in_file, squeeze_me=True)
    electrode_names = np.loadtxt(electrode_name_file, dtype=str)
    electrodes = struct['ELECTRODESF']
    position = electrodes[session_number - 1]
    out_file = op.abspath(out_basename + str(session_number) + ".geo")

    in_t1 = nb.load(original_T1)
    t1_data = in_t1.get_data()
    t1_header = in_t1.get_data()
    f = open(out_file, 'w')
    print 'Writing points to {f}'.format(f=out_file)
    for pt_id, data in enumerate(position):
        x = -float(data[0])+112
        y = float(data[2])-127.5
        z = float(data[1])-127.5
        pt_str = ('Point(%d) = {%5.2f,%5.2f,%5.2f};' % (pt_id, x, y, z))
        f.write(pt_str + '\n')

    f.write('View "Electrode Names" {\n')
    for pt_id, data in enumerate(position):
        x = -float(data[0])+112
        y = float(data[2])-127.5
        z = float(data[1])-127.5
        pt_str = ('T3(%5.2f, %5.2f, %5.2f,0){"%s"};' % (
            x, y, z, electrode_names[pt_id]))
        f.write(pt_str + '\n')
    f.write('};\n')

    f.close()
    return out_file


def get_tms_coil_position(in_file, original_T1, out_basename="tmscoil_", session_number=1):
    import scipy.io as sio
    import os
    import os.path as op
    import numpy as np
    struct = sio.loadmat(in_file, squeeze_me=True)
    data = struct['SESS_MEAN']
    coildata_all = data[:, 2:11]
    # (Coil LOC X,Y,Z), (Coil Normal X,Y,Z), (COIL Dir X,Y,Z)
    data = coildata_all[session_number - 1]
    out_file = op.abspath(out_basename + str(session_number) + ".geo")

    # Y and Z are reversed in this function because we are using
    # a NexStim TMS EEG system, and they use a really stupid coordinate
    # system. Same reason why we need to involve the MRI dimensions
    f = open(out_file, 'w')

    print('Writing coil location, normal, and direction to {f}'.format(f=out_file))
    x = -float(data[0])+112
    y = float(data[2])-128
    z = float(data[1])-128
    pt_str = ('Point(1000) = {%5.2f,%5.2f,%5.2f};' % (x, y, z))
    f.write(pt_str + '\n')
    f.write('View "TMS Coil Location" {\n')
    pt_str = ('T3(%5.2f, %5.2f, %5.2f,0){"Coil Location"};' % (x, y, z))
    f.write(pt_str + '\n')
    f.write('};\n')

    n_x = float(data[3])
    n_y = float(data[4])
    n_z = float(data[5])

    f.write('View "TMS Coil Normal" {\n')
    pt_str = (
        'VP(%5.2f, %5.2f, %5.2f){%5.2f, %5.2f, %5.2f};' % (x, y, z, n_x, n_y, n_z))
    f.write(pt_str + '\n')
    f.write('};\n')

    d_x = float(data[6])
    d_y = float(data[7])
    d_z = float(data[8])
    f.write('View "TMS Coil Direction" {\n')
    pt_str = (
        'VP(%5.2f, %5.2f, %5.2f){%5.2f, %5.2f, %5.2f};' % (x, y, z, d_x, d_y, d_z))
    f.write(pt_str + '\n')
    f.write('};\n')

    f.close()
    return out_file

import scipy.io as sio
in_file = "/Users/erik/Dropbox/Analysis/TMSEEG/TMS007/TMS/SESS_CHAR_TMS007.mat"
electrode_name_file = "/Users/erik/Dropbox/Analysis/TMSEEG/ElectrodeNames.txt"

original_T1 = "/Users/erik/Dropbox/Analysis/TMSEEG/TMS007/20121204_134045t1mpragetms20121113s002a001.nii"
final_volume = "/Users/erik/Dropbox/Analysis/TMSEEG/TMS007/orig_out_flirt.nii"
get_electrode_positions(in_file, original_T1, final_volume, electrode_name_file)
get_tms_coil_position(in_file, original_T1)
