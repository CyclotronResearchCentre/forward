import os
import os.path as op
from nipype.interfaces.base import File, Directory, traits
from nipype.interfaces.matlab import MatlabCommand, MatlabInputSpec


class FourShellAnalyticalModelInputSpec(MatlabInputSpec):
    probe_dipole_location = traits.List(
        traits.Float, exists=True, minlen=3, maxlen=3, mandatory=True)
    probe_dipole_orientation = traits.List([0,0,1],
        traits.Float, exists=True, minlen=3, maxlen=3, usedefault=True, mandatory=True)
    sphere_radii = traits.List(
        traits.Float, exists=True, minlen=4, maxlen=4, mandatory=True)
    shell_conductivity = traits.List(
        traits.Float, exists=True, minlen=4, maxlen=4, mandatory=True)
    icosahedron_sides = traits.Enum([42, 162, 642], mandatory=True)
    number_of_analytical_terms = traits.Int(60, usedefault=True, mandatory=True)
    fieldtrip_path = Directory(exists=True, desc='Fieldtrip directory')
    out_analytical_file = File('analytical.txt', usedefault=True)
    out_openmeeg_file = File('openmeeg.txt', usedefault=True)

class FourShellAnalyticalModelOutputSpec(MatlabInputSpec):
    out_analytical_file = File(exists=True)
    out_openmeeg_file = File(exists=True)
    matlab_output = traits.Str()


class FourShellAnalyticalModel(MatlabCommand):

    """ Runs a four-shell spherical model using OpenMEEG and
    Fieldtrip's analytical solution

    Examples
    --------

    >>> four = FourShellAnalyticalModel()
    """
    input_spec = FourShellAnalyticalModelInputSpec
    output_spec = FourShellAnalyticalModelOutputSpec

    def _my_script(self):

        probe_location_list = [str(i) for i in self.inputs.probe_dipole_location]
        probe_orientation_list = [str(i) for i in self.inputs.probe_dipole_orientation]
        radii_list = [str(i) for i in self.inputs.sphere_radii]
        cond_list = [str(i) for i in self.inputs.shell_conductivity]

        fieldtrip_path = op.abspath(self.inputs.fieldtrip_path)
        private_fwd_path = op.join(os.environ["FWD_DIR"],"etc")
        probe_dipole_location = "[" + " ".join(probe_location_list) + ']'
        probe_dipole_orientation = "[" + " ".join(probe_orientation_list) + ']'
        sphere_radii = "[" + " ".join(radii_list) + ']'
        shell_conductivity = "[" + " ".join(cond_list) + ']'
        icosahedron_sides = self.inputs.icosahedron_sides
        num_terms = int(self.inputs.number_of_analytical_terms)
        out_analytical_file = op.abspath(self.inputs.out_analytical_file)
        out_openmeeg_file = op.abspath(self.inputs.out_openmeeg_file)

        script = """

    addpath(genpath('%s'));
    addpath(genpath('%s'));

    pos = %s;
    moment = %s;
    num_terms = %d;

    r = %s;
    c = %s;
    [pnt, tri] = icosahedron%d;
    
    sens.elecpos = max(r) * pnt;
    sens.label = {};
    nsens = size(sens.elecpos,1);
    for ii=1:nsens
        sens.label{ii} = sprintf('vertex%%03d', ii);
    end

    vol_sphere.r = r;
    vol_sphere.cond = c;
    vol_sphere.t = eeg_leadfield4_prepare(vol_sphere, num_terms);
    lf_sphere = ft_compute_leadfield(pos, sens, vol_sphere);
    dlmwrite('%s', lf_sphere)

    %% Create OpenMEEG solution
    %% Create a triangulated mesh, the first boundary is inside

    vol = [];
    for ii=1:length(r)
        vol.bnd(ii).pnt = pnt * r(ii);
        vol.bnd(ii).tri = tri;
    end

    cfg.method = 'openmeeg';
    cfg.conductivity = c;
    vol1 = ft_prepare_bemmodel(cfg, vol);
    vol_bem = vol1;
    cfg.vol = vol_bem;
    cfg.grid.pos = pos;
    cfg.grid.mom = moment;
    cfg.elec = sens;

    grid = ft_prepare_leadfield(cfg);
    lf_openmeeg = grid.leadfield{1};
    dlmwrite('%s', lf_openmeeg)

        """ % (fieldtrip_path, private_fwd_path, probe_dipole_location, probe_dipole_orientation, num_terms, 
            sphere_radii, shell_conductivity, icosahedron_sides, out_analytical_file, out_openmeeg_file)
        return script

    def run(self, **inputs):
        # inject your script
        self.inputs.script = self._my_script()
        results = super(MatlabCommand, self).run(**inputs)
        stdout = results.runtime.stdout
        # attach stdout to outputs to access matlab results
        results.outputs.matlab_output = stdout
        return results

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs['out_analytical_file'] = op.abspath(self.inputs.out_analytical_file)
        outputs['out_openmeeg_file'] = op.abspath(self.inputs.out_openmeeg_file)
        return outputs
