import numpy as np
import os.path as op
import subprocess
from forward.electrodes import rewrite_mesh_with_electrodes


def create_sphere(radius, vol_id, lc=10, out_file="newsphere.msh"):
    f = open(op.abspath(out_file), "w")
    f.write("lc = %f;\n" % lc)
    f.write("R = %f;\n" % radius)
    f.write("x_center = 0;\n")
    f.write("y_center = 0;\n")
    f.write("z_center = 0;\n")
    f.write("Point(1) = {x_center, y_center, z_center, lc};\n")
    f.write("Point(2) = {x_center - R, y_center, z_center, lc};\n")
    f.write("Point(4) = {x_center, y_center - R, z_center, lc};\n")
    f.write("Point(5) = {x_center + R, y_center, z_center, lc};\n")
    f.write("Point(8) = {x_center, y_center, z_center - R, lc};\n")
    f.write("Point(11) = {x_center, y_center + R, z_center, lc};\n")
    f.write("Point(14) = {x_center, y_center, z_center + R, lc};\n")
    f.write("Circle (1) = {2, 1, 4} Plane{0, 0, 1};\n")
    f.write("Circle (2) = {4, 1, 5} Plane{0, 0, 1};\n")
    f.write("Circle (3) = {2, 1, 8} Plane{0, 0, 1};\n")
    f.write("Circle (4) = {4, 1, 8} Plane{0, 0, 1};\n")
    f.write("Circle (6) = {2, 1, 11} Plane{0, 0, 1};\n")
    f.write("Circle (7) = {8, 1, 11} Plane{0, 0, 1};\n")
    f.write("Circle (9) = {2, 1, 14} Plane{0, 0, 1};\n")
    f.write("Circle (10) = {11, 1, 14} Plane{0, 0, 1};\n")
    f.write("Circle (13) = {14, 1, 4} Plane{0, 0, 1};\n")
    f.write("Circle (15) = {8, 1, 5} Plane{0, 0, 1};\n")
    f.write("Circle (18) = {11, 1, 5} Plane{0, 0, 1};\n")
    f.write("Circle (21) = {14, 1, 5} Plane{0, 0, 1};\n")
    f.write("Line Loop (1000005) = {1, 4, -3};\n")
    f.write("Ruled Surface (5) = {1000005};\n")
    f.write("Line Loop (1000008) = {3, 7, -6};\n")
    f.write("Ruled Surface (8) = {1000008};\n")
    f.write("Line Loop (1000011) = {6, 10, -9};\n")
    f.write("Ruled Surface (11) = {1000011};\n")
    f.write("Line Loop (1000014) = {9, 13, -1};\n")
    f.write("Ruled Surface (14) = {1000014};\n")
    f.write("Line Loop (1000017) = {-15, -4, 2};\n")
    f.write("Ruled Surface (17) = {1000017};\n")
    f.write("Line Loop (1000020) = {-18, -7, 15};\n")
    f.write("Ruled Surface (20) = {1000020};\n")
    f.write("Line Loop (1000023) = {-21, -10, 18};\n")
    f.write("Ruled Surface (23) = {1000023};\n")
    f.write("Line Loop (1000026) = {-2, -13, 21};\n")
    f.write("Ruled Surface (26) = {1000026};\n")
    f.write("Surface Loop(%d) = { 5, 8, 11, 14, 17, 20, 23, 26 };\n" % vol_id)
    f.close()
    return out_file


def mesh_2D(in_file):
    out_file = op.basename(in_file).replace('.geo', '.stl')
    call_list = ["gmsh", in_file, "-v", "0", "-2", "-format", "stl",
                 "-o", out_file]
    print(" ".join(call_list))
    subprocess.call(call_list)
    return out_file


def merge_and_diff(in_files, ids_outside_inward, out_file):
    merge_script = op.abspath("merge.geo")
    f = open(merge_script, "w")
    f.write("Mesh.Algorithm3D = 4;\n")
    f.write("Mesh.Optimize = 1;\n")
    f.write("Mesh.OptimizeNetgen = 1;\n")
    for in_file in in_files:
        f.write('Merge "%s";\n' % in_file)

    f.write("Surface Loop (1) = {1}; //brain\n")
    f.write("Surface Loop (2) = {2}; //csf\n")
    f.write("Surface Loop (3) = {3}; //skull\n")
    f.write("Surface Loop (4) = {4}; //skin\n")

    f.write("Volume(4) = { 3,4 }; //skin volume\n")
    f.write("Volume(3) = { 2,3 }; //skull volume\n")
    f.write("Volume(2) = { 1,2 }; //csf volume\n")
    f.write("Volume(1) = { 1 }; //brain volume\n")

    f.write("Physical Surface(2001) = { 1 };\n")
    f.write("Physical Surface(2002) = { 2 };\n")
    f.write("Physical Surface(2003) = { 3 };\n")
    f.write("Physical Surface(2004) = { 4 };\n")

    f.write("Physical Volume(1004) = { 4 }; //skin volume\n")
    f.write("Physical Volume(1003) = { 3 }; //skull volume\n")
    f.write("Physical Volume(1002) = { 2 }; //csf volume\n")
    f.write("Physical Volume(1001) = { 1 }; //brain volume\n")

    f.close()
    out_file = op.abspath(out_file)
    call_list = ["gmsh", merge_script, "-v", "0", "-3",
                 "-format", "msh", "-o", out_file]
    print(" ".join(call_list))
    subprocess.call(call_list)

    return out_file


def add_sensors(electrode_location_file, mesh_file, radii):
    points = np.loadtxt(electrode_location_file, delimiter=",")
    points = np.max(radii) * points
    new_electrode_location_file = op.abspath("electrode_locations.txt")
    np.savetxt(new_electrode_location_file, points, delimiter=",")
    elec_name_file = op.abspath(
        op.basename(electrode_location_file) + "_names.txt")
    f = open(elec_name_file, "w")
    for idx in range(0, len(points)):
        name = "vertex%03d" % (idx + 1)
        f.write("%s\n" % name)
    f.close()

    print('Rewriting mesh with sensors')
    print(new_electrode_location_file, elec_name_file, mesh_file)
    out_file = rewrite_mesh_with_electrodes(
        new_electrode_location_file,
        electrode_name_file=elec_name_file,
        mesh_file=mesh_file,
        mesh_id=1004)
    return op.abspath(out_file), elec_name_file, new_electrode_location_file


def create_4_shell_model(electrode_location_file=None, radii=[85, 88, 92, 100], out_file='4shell.msh'):
    from forward.sphere import (create_sphere, mesh_2D,
                                merge_and_diff, add_sensors)
    assert(len(radii) == 4)

    characteristic_length_brain = 3
    characteristic_length_csf = 7
    characteristic_length_skull = 7
    characteristic_length_skin = 7

    brain = create_sphere(
        radii[0], 1, characteristic_length_brain, out_file="brain_sphere.geo")
    brain_msh = mesh_2D(brain)

    csf = create_sphere(radii[1], 2, characteristic_length_csf,
                        out_file="csf_sphere.geo")
    csf_msh = mesh_2D(csf)

    skull = create_sphere(
        radii[2], 3, characteristic_length_skull, out_file="skull_sphere.geo")
    skull_msh = mesh_2D(skull)

    skin = create_sphere(
        radii[3], 4, characteristic_length_skin, out_file="skin_sphere.geo")
    skin_msh = mesh_2D(skin)

    mesh_without_electrodes = merge_and_diff([brain_msh, csf_msh, skull_msh, skin_msh],
                               ids_outside_inward=[4, 3, 2, 1], out_file=out_file)
    if electrode_location_file is not None:
        mesh_file, electrode_name_file, electrode_location_file = add_sensors(
            electrode_location_file, mesh_without_electrodes, radii)
        return mesh_file, electrode_name_file, electrode_location_file, mesh_without_electrodes
    else:
        return mesh_without_electrodes, None, None, None

