from forward.analytical import FourShellAnalyticalModel

a = FourShellAnalyticalModel()
a.inputs.probe_dipole = [0,0,70]
a.inputs.sphere_radii = [85,88,92,100]
a.inputs.shell_conductivity = [1, 1/20.0, 1/80.0, 1]
a.inputs.icosahedron_sides = 42
a.inputs.fieldtrip_path = "/Developer/fieldtrip"
a.run()