def dipoles_to_geo(in_file):
    import os.path as op
    import numpy as np
    from nipype.utils.filemanip import split_filename
    dipole_data = np.loadtxt(in_file)
    _, name, _ = split_filename(in_file)
    out_file = op.abspath(name + ".geo")
    f = open(out_file, 'w')
    print('Writing dipole sources to {f}'.format(f=out_file))
    f.write('View "Dipole locations" {\n')
    for data in dipole_data:
        x = float(data[2])
        y = float(data[3])
        z = float(data[4])
        q_x = float(data[6])
        q_y = float(data[7])
        q_z = float(data[8])
        pt_str = (
            'VP(%5.2f, %5.2f, %5.2f){%5.2f, %5.2f, %5.2f};' % (x, y, z, q_x, q_y, q_z))
        f.write(pt_str + '\n')
    f.write('};\n')
    return out_file

def dipoles_to_dat(in_file):
    import os.path as op
    import numpy as np
    from nipype.utils.filemanip import split_filename
    dipole_data = np.loadtxt(in_file)
    _, name, _ = split_filename(in_file)
    out_file = op.abspath(name + ".dat")
    f = open(out_file, 'w')
    print('Writing dipole sources to {f}'.format(f=out_file))
    f.write('View "Dipole locations" {\n')
    for data in dipole_data:
        x = float(data[2])
        y = float(data[3])
        z = float(data[4])
        q_x = float(data[6])
        q_y = float(data[7])
        q_z = float(data[8])
        pt_str = (
            'VP(%5.2f, %5.2f, %5.2f){%5.2f, %5.2f, %5.2f};' % (x, y, z, q_x, q_y, q_z))
        f.write(pt_str + '\n')
    f.write('};\n')
    return out_file

dipoles_to_geo("/Applications/freesurfer/subjects/TMS007/bem/TMS007-7-rh.dip")
dipoles_to_geo("/Applications/freesurfer/subjects/TMS007/bem/TMS007-7-lh.dip")
