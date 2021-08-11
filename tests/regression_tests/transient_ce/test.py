from difflib import unified_diff
import filecmp
import glob
import h5py
import numpy as np
import os

import openmc
import openmc.kinetics as kinetics
import openmc.mgxs

from tests.testing_harness import TestHarness, colorize


class TransientTestHarness(TestHarness):
    """Specialized TestHarness for running OpenMC transient tests."""

    def __init__(self, solver, log_name, model=None, inputs_true=None):
        self.solver = solver
        self._log_name = log_name

        if model is None:
            self._model = pwr_core()
        else:
            self._model = model
        self._model.plots = []

        self.inputs_true = "inputs_true.dat" if not inputs_true else inputs_true

    def execute_test(self):
        """Don't call _run_openmc as OpenMC will be called through solver.py for
        transient tests.

        """
        try:
            self._build_inputs()
            inputs = self._get_inputs()
            self._write_inputs(inputs)
            self._compare_inputs()
            self._run_transient()
            self._test_output_created()
            results = self._get_results()
            self._write_results(results)
            self._compare_results()
        finally:
            self._cleanup()

    def _build_inputs(self):
        """Write input XML files."""
        self._model.export_to_xml()

    def _get_inputs(self):
        """Return a hash digest of the input XML files."""
        xmls = ['geometry.xml', 'materials.xml', 'settings.xml',
                'tallies.xml', 'plots.xml']
        return ''.join([open(fname).read() for fname in xmls
                        if os.path.exists(fname)])

    def _write_inputs(self, input_digest):
        """Write the digest of the input XMLs to an ASCII file."""
        with open('inputs_test.dat', 'w') as fh:
            fh.write(input_digest)

    def _compare_inputs(self):
        """Make sure the current inputs agree with the _true standard."""
        compare = filecmp.cmp('inputs_test.dat', self.inputs_true)
        if not compare:
            expected = open(self.inputs_true, 'r').readlines()
            actual = open('inputs_test.dat', 'r').readlines()
            diff = unified_diff(expected, actual, self.inputs_true,
                                'inputs_test.dat')
            print('Input differences:')
            print(''.join(colorize(diff)))
            os.rename('inputs_test.dat', 'inputs_error.dat')
        assert compare, 'Input files are broken.'

    def _cleanup(self):
        """Delete XMLs, statepoints, tally, and test files."""
        super()._cleanup()
        output = ['materials.xml', 'geometry.xml', 'settings.xml',
                  'tallies.xml', 'plots.xml', 'inputs_test.dat']
        output += glob.glob('*.h5')
        for f in output:
            if os.path.exists(f):
                os.remove(f)

    def _run_transient(self):
        self.solver.solve()

    def _test_output_created(self):
        """Make sure statepoint.* and tallies.out have been created."""
        logfile = glob.glob(self._log_name)
        assert len(logfile) == 1, 'Either multiple or no log files' \
            ' exist.'
        assert logfile[0].endswith('h5'), \
            'Log file is not a HDF5 file.'

    def _get_results(self, hash_output=False):
        """Digest info in the log file and return as a string."""
        # Read the log file.
        logfile = glob.glob(self._log_name)[0]
        with h5py.File(logfile, 'r') as lg:
            outstr = ''
            outstr += 'k_crit:\n'
            form = '{0:12.6E}\n'
            outstr += form.format(lg.attrs['k_crit'])

            # Write out power data
            results = list(lg['power'])
            results = ['{0:12.6E}'.format(x) for x in results]
            outstr += 'power:\n'
            outstr += '\n'.join(results) + '\n'

        return outstr

def test_transient():
    """Test runs a continuous energy adibatic transient calculation
    using OpenMC kinetics.

    The transient consists of a prescribed water density drop.

    """
    model = openmc.model.Model()

    uo2 = openmc.Material(name='UO2')
    uo2.set_density('g/cm3', 10.29769)
    uo2.add_element('U', 1., enrichment=2.4)
    uo2.add_element('O', 2.)

    zircaloy = openmc.Material(name='Zircaloy')
    zircaloy.set_density('g/cm3', 6.55)
    zircaloy.add_nuclide('Zr90', 1.0)

    borated_water = openmc.Material(name='Moderator')
    borated_water.set_density('g/cm3', 0.740582) 
    borated_water.add_element('B', 2.7800E-5, 'ao')
    borated_water.add_element('H', 2*3.3500E-2, 'ao')
    borated_water.add_element('O', 3.3500E-2, 'ao')
    borated_water.add_s_alpha_beta('c_H_in_H2O')

    model.materials.extend([uo2, zircaloy, borated_water])

    fuel_or = openmc.ZCylinder(r=0.39218)
    clad_or = openmc.ZCylinder(r=0.45720)

    pitch = 1.26
    box = openmc.rectangular_prism(pitch, pitch, boundary_type='reflective')
    bottom = openmc.ZPlane(z0=-10, boundary_type='vacuum')
    top = openmc.ZPlane(z0= 10, boundary_type='vacuum')

    fuel = openmc.Cell(fill=uo2, region=-fuel_or & +bottom & -top)
    clad = openmc.Cell(fill=zircaloy, region=+fuel_or & -clad_or & +bottom & -top)
    water = openmc.Cell(fill=borated_water, region=+clad_or & box & +bottom & -top)

    model.geometry = openmc.Geometry([fuel, clad, water])
    
    model.settings.batches = 100
    model.settings.inactive = 10
    model.settings.particles = 100
    model.settings.output = {'tallies': False}

    lower_left = [-0.63, -0.63, -10]
    upper_right = [0.63, 0.63, 10]
    uniform_dist = openmc.stats.Box(lower_left, upper_right, only_fissionable=True)
    model.settings.source = openmc.source.Source(space=uniform_dist)

    sourcepoint = dict()
    sourcepoint['batches'] = [model.settings.batches]
    sourcepoint['separate'] = True
    sourcepoint['write'] = True
    model.settings.sourcepoint = sourcepoint

    full_pin_cell_mesh = openmc.RegularMesh()
    full_pin_cell_mesh.dimension = [1,1,20]
    full_pin_cell_mesh.lower_left = lower_left
    full_pin_cell_mesh.upper_right = upper_right

    energy_groups = openmc.mgxs.EnergyGroups()
    energy_groups.group_edges = [0, 0.625, 20.e6]

    one_group = openmc.mgxs.EnergyGroups()
    one_group.group_edges = [0, 20.e6]

    t_outer = np.arange(0., 1.0, 5.e-1)
    clock = openmc.kinetics.Clock(dt_inner=1.e-2, t_outer=t_outer)

    transient = {}
    for material in model.materials:
        MatChange = {material.name: {},}
        transient.update(MatChange)
        for t in t_outer:
            time_dict = {
            t: {
                'density' : {},
                'temperature' : {},
                }
            }
            transient[material.name].update(time_dict)

    for material in model.materials:
        if material.name == 'Moderator':
            transient[material.name][0]['density'] = material.density
            transient[material.name][0.5]['density'] = material.density*0.9
            for t in t_outer:
                transient[material.name][t]['temperature'] = material.temperature
        else: 
            for t in t_outer:
                transient[material.name][t]['density'] = material.density
                transient[material.name][t]['temperature'] = material.temperature

    # Initialize solver object 
    solver = openmc.kinetics.Solver()
    solver.amplitude_mesh  = full_pin_cell_mesh
    solver.shape_mesh      = full_pin_cell_mesh
    solver.tally_mesh      = full_pin_cell_mesh
    solver.one_group       = one_group
    solver.energy_groups   = energy_groups
    solver.fine_groups     = energy_groups
    solver.tally_groups    = energy_groups
    solver.geometry        = model.geometry
    solver.settings        = model.settings
    solver.materials       = model.materials
    solver.transient       = transient
    solver.outer_tolerance = np.inf
    solver.initial_power   = 1.0
    solver.core_volume     = 1.0
    solver.method          = 'ADIABATIC'
    solver.multi_group     = False
    solver.clock           = clock
    solver.run_kwargs      = {'threads': None, 'mpi_args': None}
    solver.min_outer_iters = 1
    use_pcmfd              = True
    use_agd                = True
    use_pregenerated_sps   = False
    solver.chi_delayed_by_delayed_group = True
    solver.chi_delayed_by_mesh          = True
    solver.log_file_name   = 'log_file.h5'

    # Initialize and run solver object
    harness = TransientTestHarness(solver, 'log_file.h5', model=model)
    harness.main()
