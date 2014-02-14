## Current source from dipole
## Voltage sink at sink electrode on scalp



def write_dipole_pro_file(mesh_file, conductivity_tensor_included, in_file, row=0, orig_pro_file=None, dipole_name="Dipole1"):
    import numpy as np
    import os
    import os.path as op

    if orig_pro_file is None:
        orig_pro_file = op.join(os.environ["FWD_DIR"], "etc/eeg_dipole.pro")

    problem_file = op.abspath("eeg_forward_dipole_" + str(dipole_name) + ".pro")

    original = open(orig_pro_file, 'r')
    modified = open(problem_file, 'w')

    while True:
        line = original.readline()
        if len(line) == 0:
            break # EOF
        
        if line.rfind('Sink = Region[') >= 0:
            line = '  Sink = Region[' + str(ground_index+5000) + '];\n'


        elif line.rfind('sigma[WhiteMatter_Cerebellum]') >= 0 and conductivity_tensor_included:
            line = "    sigma[WhiteMatter_Cerebellum] = TensorField[XYZ[]] #1 ? #1 : 0.33 ;//: 0.33*1e-3 ;\n"
        elif line.rfind('sigma[WhiteMatter_Cerebellum]') >= 0 and not conductivity_tensor_included:
            line = "    sigma[WhiteMatter_Cerebellum] = 0.33;\n"
        modified.write(line)

    modified.close()
    original.close()
    return problem_file


def rewrite_mesh_with_dipole(dipole_file, mesh_file, mesh_id, row=0):
    import os.path as op
    import numpy as np
    from nipype.utils.filemanip import split_filename
    from forward.electrodes import get_closest_element_to_points
    _, name, _ = split_filename(mesh_file)
    out_basename = op.abspath(name + "_dipole")

    dipole_data = np.loadtxt(dipole_file)
    dip = dipole_data[int(row)]
    
    import ipdb
    ipdb.set_trace()

    x = float(dip[2])
    y = float(dip[3])
    z = float(dip[4])
    q_x = float(dip[6])
    q_y = float(dip[7])
    q_z = float(dip[8])

    positions = [[x, y, z]]

    new_phys_ids = [8000]
    new_elem_ids = [9000]

    out_file = get_closest_element_to_points(
        positions, mesh_file, mesh_id, new_phys_ids, new_elem_ids, out_basename)
    return out_file


def dipoles_to_geo(in_file, row=None, dipole_name="Dipole1"):
    import os.path as op
    import numpy as np
    from nipype.utils.filemanip import split_filename
    dipole_data = np.loadtxt(in_file)
    _, name, _ = split_filename(in_file)
    out_file = op.abspath(name + "_" + dipole_name + ".geo")
    f = open(out_file, 'w')
    print('Writing dipole .geo file to {f}'.format(f=out_file))
    f.write('View "Dipole locations" {\n')

    if row is not None:
        dipole_data = [dipole_data[int(row)]]

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
    f.close()
    return out_file

