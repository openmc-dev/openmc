from difflib import unified_diff
import filecmp
import glob
import h5py
import numpy as np
import os

import openmc
import openmc.kinetics as kinetics
import openmc.mgxs

from tests.regression_tests import config
from tests.testing_harness import PyAPITestHarness, colorize


class TransientTestHarness(PyAPITestHarness):
    """Specialized TestHarness for running OpenMC transient tests."""

    def __init__(self, solver, log_name, model=None, inputs_true=None, statepoint_name=None):
        super().__init__(statepoint_name, model, inputs_true)
        self.solver = solver
        self._log_name = log_name

    def _cleanup(self):
        """Delete XMLs, statepoints, tally, and test files."""
        super()._cleanup()
        output = ['materials.xml', 'geometry.xml', 'settings.xml',
                  'tallies.xml', 'plots.xml', 'inputs_test.dat']
        output += glob.glob('*.h5')
        for f in output:
            if os.path.exists(f):
                os.remove(f)

    def _run_openmc(self):
        self.solver.solve()

    def _test_output_created(self):
        """Make sure statepoint.* and tallies.out have been created."""
        logfile = glob.glob(self._log_name)
        assert len(logfile) == 1, 'Either multiple or no log files' \
            ' exist.'
        assert logfile[0].endswith('h5'), 'Log file is not a HDF5 file.'

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
    top = openmc.ZPlane(z0=10, boundary_type='vacuum')

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

    model.settings.sourcepoint = {
        'batches': [model.settings.batches],
        'separate': True,
        'write': True
    }

    full_pin_cell_mesh = openmc.RegularMesh()
    full_pin_cell_mesh.dimension = [1, 1, 20]
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
        transient[material.name] = {}
        for t in t_outer:
            transient[material.name][t] = {
                'density': material.density,
                'temperature': material.temperature
            }
        if material.name == 'Water':
            transient[material.name][0.5]['density'] = material.density*0.9

    for cell in model.geometry.get_all_cells().values():
        transient[cell.name] = {}
        for t in t_outer:
            transient[cell.name][t] = {'translation': None}

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
    solver.method          = 'adiabatic'
    solver.multi_group     = False
    solver.clock           = clock
    solver.min_outer_iters = 1
    use_pcmfd              = True
    use_agd                = True
    use_pregenerated_sps   = False
    solver.chi_delayed_by_delayed_group = True
    solver.chi_delayed_by_mesh          = True
    solver.log_file_name   = 'log_file.h5'

    if config['mpi']:
        mpi_args = [config['mpiexec'], '-n', config['mpi_np']]
    else:
        mpi_args = None
    solver.run_kwargs = {'threads': None, 'mpi_args': mpi_args}
    
    # Initialize and run solver object
    harness = TransientTestHarness(solver, 'log_file.h5', model=model)
    harness.main()
