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

os.environ['OPENMC_MG_CROSS_SECTIONS'] = 'mgxs.h5'


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

    groups = openmc.mgxs.EnergyGroups(group_edges=[
        0, 0.0635, 10.0, 1.0e2, 1.0e3, 0.5e6, 1.0e6, 20.0e6])

    uo2_xsdata = openmc.XSdata('UO2', groups)
    uo2_xsdata.order = 0
    uo2_xsdata.num_delayed_groups = 8
    uo2_xsdata.set_total(
        [0.1779492, 0.3298048, 0.4803882, 0.5543674, 0.3118013, 0.3951678,
        0.5644058])
    uo2_xsdata.set_absorption([8.0248E-03, 3.7174E-03, 2.6769E-02, 9.6236E-02,
        3.0020E-02, 1.1126E-01, 2.8278E-01])
    scatter_matrix = np.array(
        [[[0.1275370, 0.0423780, 0.0000094, 0.0000000, 0.0000000, 0.0000000, 0.0000000],
        [0.0000000, 0.3244560, 0.0016314, 0.0000000, 0.0000000, 0.0000000, 0.0000000],
        [0.0000000, 0.0000000, 0.4509400, 0.0026792, 0.0000000, 0.0000000, 0.0000000],
        [0.0000000, 0.0000000, 0.0000000, 0.4525650, 0.0055664, 0.0000000, 0.0000000],
        [0.0000000, 0.0000000, 0.0000000, 0.0001253, 0.2714010, 0.0102550, 0.0000000],
        [0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0012968, 0.2658020, 0.0168090],
        [0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0085458, 0.2730800]]])
    scatter_matrix = np.rollaxis(scatter_matrix, 0, 3)
    uo2_xsdata.set_scatter_matrix(scatter_matrix)
    uo2_xsdata.set_inverse_velocity([1./2.23466E+09, 1./5.07347E+08, 1./3.86595E+07,
        1./5.13931E+06, 1./1.67734E+06, 1./7.28603E+05, 1./2.92902E+05])
    uo2_xsdata.set_fission([7.21206E-03, 8.19301E-04, 6.45320E-03,
        1.85648E-02, 1.78084E-02, 8.30348E-02,
        2.16004E-01])
    uo2_xsdata.set_nu_fission([2.005998E-02, 2.027303E-03, 1.570599E-02,
        4.518301E-02, 4.334208E-02, 2.020901E-01,
        5.257105E-01])
    uo2_xsdata.set_kappa_fission([7.21206E-03, 8.19301E-04, 6.45320E-03,
        1.85648E-02, 1.78084E-02, 8.30348E-02,
        2.16004E-01])
    uo2_xsdata.set_chi([5.8791E-01, 4.1176E-01, 3.3906E-04, 1.1761E-07, 0.0000E+00,
        0.0000E+00, 0.0000E+00])
    uo2_xsdata.set_chi_delayed([[0.00075, 0.98512, 0.01413, 0.00000, 0.00000, 0.00000, 0.00000],
        [0.03049, 0.96907, 0.00044, 0.00000, 0.00000, 0.00000, 0.00000],
        [0.00457, 0.97401, 0.02142, 0.00000, 0.00000, 0.00000, 0.00000],
        [0.02002, 0.97271, 0.00727, 0.00000, 0.00000, 0.00000, 0.00000],
        [0.05601, 0.93818, 0.00581, 0.00000, 0.00000, 0.00000, 0.00000],
        [0.06098, 0.93444, 0.00458, 0.00000, 0.00000, 0.00000, 0.00000],
        [0.10635, 0.88298, 0.01067, 0.00000, 0.00000, 0.00000, 0.00000],
        [0.09346, 0.90260, 0.00394, 0.00000, 0.00000, 0.00000, 0.00000]])
    uo2_xsdata.set_beta([2.13333E-04, 1.04514E-03, 6.03969E-04, 1.33963E-03, 
        2.29386E-03, 7.05174E-04, 6.00381E-04, 2.07736E-04])
    uo2_xsdata.set_decay_rate([1.247E-02, 2.829E-02, 4.252E-02, 1.330E-01, 
        2.925E-01, 6.665E-01, 1.635E+00, 3.555E+00])

    h2o_xsdata = openmc.XSdata('LWTR', groups)
    h2o_xsdata.order = 0
    h2o_xsdata.num_delayed_groups = 8
    h2o_xsdata.set_total([0.15920605, 0.412969593, 0.59030986, 0.58435,
        0.718, 1.2544497, 2.650379])
    h2o_xsdata.set_absorption([6.0105E-04, 1.5793E-05, 3.3716E-04,
        1.9406E-03, 5.7416E-03, 1.5001E-02,
        3.7239E-02])

    scatter_matrix = np.array(
        [[[0.0444777, 0.1134000, 0.0007235, 0.0000037, 0.0000001, 0.0000000, 0.0000000],
        [0.0000000, 0.2823340, 0.1299400, 0.0006234, 0.0000480, 0.0000074, 0.0000010],
        [0.0000000, 0.0000000, 0.3452560, 0.2245700, 0.0169990, 0.0026443, 0.0005034],
        [0.0000000, 0.0000000, 0.0000000, 0.0910284, 0.4155100, 0.0637320, 0.0121390],
        [0.0000000, 0.0000000, 0.0000000, 0.0000714, 0.1391380, 0.5118200, 0.0612290],
        [0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0022157, 0.6999130, 0.5373200],
        [0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.1324400, 2.4807000]]])
    scatter_matrix = np.rollaxis(scatter_matrix, 0, 3)
    h2o_xsdata.set_scatter_matrix(scatter_matrix)
    h2o_xsdata.set_inverse_velocity([1./2.23466E+09, 1./5.07347E+08, 1./3.86595E+07,
        1./5.13931E+06, 1./1.67734E+06, 1./7.28603E+05, 1./2.92902E+05])

    mg_cross_sections_file = openmc.MGXSLibrary(groups, 8)
    mg_cross_sections_file.add_xsdatas([uo2_xsdata, h2o_xsdata])
    mg_cross_sections_file.export_to_hdf5()

    uo2_data = openmc.Macroscopic('UO2')
    h2o_data = openmc.Macroscopic('LWTR')

    uo2 = openmc.Material(name='UO2 fuel')
    uo2.set_density('macro', 1.0)
    uo2.add_macroscopic(uo2_data)

    water = openmc.Material(name='Water')
    water.set_density('macro', 1.0)
    water.add_macroscopic(h2o_data)

    model.materials = openmc.Materials([uo2, water])
    model.materials.cross_sections = "mgxs.h5"

    fuel_or = openmc.ZCylinder(r=0.54, name='Fuel OR')

    pitch = 1.26
    box = openmc.rectangular_prism(pitch, pitch, boundary_type='reflective')
    bottom = openmc.ZPlane(z0=-10, boundary_type='vacuum')
    top = openmc.ZPlane(z0= 10, boundary_type='vacuum')

    fuel = openmc.Cell(fill=uo2, region=-fuel_or & +bottom & -top, name='fuel')
    moderator = openmc.Cell(fill=water, region=+fuel_or & box & +bottom & -top, name='moderator')

    model.geometry = openmc.Geometry([fuel, moderator])
    
    model.settings.batches = 100
    model.settings.inactive = 10
    model.settings.particles = 100
    model.settings.output = {'tallies': False}
    model.settings.energy_mode = 'multi-group'

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
    full_pin_cell_mesh.lower_left = [-0.63,-0.63,-10]
    full_pin_cell_mesh.width =[0.63*2,0.63*2,1.0]

    energy_groups = openmc.mgxs.EnergyGroups()
    energy_groups.group_edges = [0, 0.0635, 10.0, 1.0e2, 1.0e3, 0.5e6, 1.0e6, 20.e6]

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

    # Initialize solver object 
    solver = openmc.kinetics.Solver()
    solver.num_delayed_groups = 8
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
    solver.mgxs_lib        = mg_cross_sections_file
    solver.multi_group     = True
    solver.clock           = clock
    solver.run_kwargs      = {'threads': 2, 'mpi_args': None}
    solver.min_outer_iters = 1

    harness = TransientTestHarness(solver, 'log_file.h5', model=model)
    harness.main()
