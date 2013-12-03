import os
import os.path as op
from nipype.interfaces.base import File, Directory, traits
from nipype.interfaces.matlab import MatlabCommand, MatlabInputSpec


class FourShellAnalyticalModelInputSpec(MatlabInputSpec):
    probe_dipole = traits.List(
        traits.Float, exists=True, minlen=3, maxlen=3, mandatory=True)
    sphere_radii = traits.List(
        traits.Float, exists=True, minlen=4, maxlen=4, mandatory=True)
    shell_conductivity = traits.List(
        traits.Float, exists=True, minlen=4, maxlen=4, mandatory=True)
    icosahedron_sides = traits.Enum([42, 162, 642], mandatory=True)
    fieldtrip_path = Directory(exists=True, desc='Fieldtrip directory')
    out_file = File('analytical.txt', usedefault=True)


class FourShellAnalyticalModelOutputSpec(MatlabInputSpec):
    out_file = File(exists=True)
    matlab_output = traits.Str()


class FourShellAnalyticalModel(MatlabCommand):

    """ Basic Hello World that displays Hello <name> in MATLAB

    Returns
    -------

    matlab_output : capture of matlab output which may be
                    parsed by user to get computation results

    Examples
    --------

    >>> four = FourShellAnalyticalModel()
    """
    input_spec = FourShellAnalyticalModelInputSpec
    output_spec = FourShellAnalyticalModelOutputSpec

    def _my_script(self):

        probe_list = [str(i) for i in self.inputs.probe_dipole]
        radii_list = [str(i) for i in self.inputs.sphere_radii]
        cond_list = [str(i) for i in self.inputs.shell_conductivity]

        fieldtrip_path = op.abspath(self.inputs.fieldtrip_path)
        probe_dipole = "[" + " ".join(probe_list) + ']'
        sphere_radii = "[" + " ".join(radii_list) + ']'
        shell_conductivity = "[" + " ".join(cond_list) + ']'
        icosahedron_sides = self.inputs.icosahedron_sides
        out_file = op.abspath(self.inputs.out_file)

        script = """

    addpath(genpath('%s'));

    pos = %s;
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
    vol_sphere.c = c;
    lf_sphere = ft_compute_leadfield(pos, sens, vol_sphere);
    dlmwrite('%s', lf_sphere)
    
        """ % (fieldtrip_path, probe_dipole, sphere_radii, shell_conductivity, icosahedron_sides, out_file)
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
        outputs['out_file'] = os.path.abspath(self.inputs.out_file)
        return outputs
