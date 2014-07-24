% First run mesh_nodes_elem_to_matlab.py
load('nodes.txt')
load('elem.txt')
load('labels.txt')
write_vista_mesh('4shell.v', nodes, elem, labels)

% Write example dipole positioned at 70 mm
% and oriented outwards along Z direction
current = [0, 0 , 70, 0, 0, 1];
write_dip_sourcesimulation(current, 'current.txt');

% Use the same electrodes as OpenM/EEG and
% Fieldtrip analytical solution
r = [82,85,9,100];
[pnt, tri] = icosahedron42;
sens.elecpos = max(r) * pnt;
write_elc(sens.elecpos, 'elec.txt')

% ./ipm_linux_opt -i sourcesimulation -p blurred_vista_4layer.par -h yourheadmesh.v -sens EEG -s yourelectrodes.elc -dip yourdipoles.dip -fwd FEM -o results.msr
