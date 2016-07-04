from collections import OrderedDict
from xml.etree import ElementTree as ET

import openmc
import openmc.kinetics
from openmc.clean_xml import *
from openmc.checkvalue import check_type
from openmc.kinetics.clock import TIME_POINTS
import numpy as np


class Solver(object):
    """Solver to propagate the neutron flux and power forward in time.

    Attributes
    ----------
    mesh : openmc.mesh.Mesh
        Mesh which specifies the dimensions of coarse mesh.

    geometry : openmc.geometry.Geometry
        Geometry which describes the problem being solved.

    settings_file : openmc.settings.SettingsFile
        Settings file describing the general settings for each simulation.

    materials_file : openmc.materials.MaterialsFile
        Materials file containing the materials info for each simulation.

    executor : openmc.executor.Executor
       Executor object for executing OpenMC simulation.

    clock : openmc.kinetics.Clock
        Clock object.

    energy_groups : openmc.mgxs.groups.EnergyGroups
        EnergyGroups which specifies the energy groups structure.

    A : np.matrix
        Numpy matrix used for storing the destruction terms.

    M : np.matrix
        Numpy matrix used for storing the production terms.

    AM : np.matrix
        Numpy matrix used for storing the combined production/destruction terms.

    flux : np.array
        Numpy array used to store the flux.

    amplitude : np.array
        Numpy array used to store the amplitude.

    shape : np.array
        Numpy array used to store the shape.

    source : np.array
        Numpy array used to store the source.

    power : np.array
        Numpy array used to store the power.

    precursor_conc : np.array
        Numpy array used to store the precursor concentrations.

    sigma_a : OrderedDict of openmc.MGXS.AbsorptionXS
        MGXS absorption multigroup cross-sections.

    nu_sigma_f : OrderedDict of openmc.MGXS.NuFissionXS
        MGXS nu-fission multigroup cross-sections.

    kappa_sigma_f : OrderedDict of openmc.MGXS.NuFissionXS
        MGXS nu-fission multigroup cross-sections.

    dif_coef : OrderedDict of openmc.MGXS.DiffusionCoefficientXS
        MGXS multigroup diffusion coefficients.

    beta : OrderedDict of openmc.MGXS.delayed.Beta
        MGXS multigroup delayed neutron fractions.

    chi_prompt : OrderedDict of openmc.MGXS.delayed.ChiPrompt
        MGXS multigroup prompt neutron spectrums.

    chi_delayed : OrderedDict of openmc.MGXS.delayed.ChiDelayed
        MGXS multigroup delayed neutron spectrums.

    velocity : OrderedDict of openmc.MGXS.Velocity
        MGXS multigroup velocities.

    nu_sigma_s : OrderedDict of openmc.MGXS.NuScatterMatrixXS
        MGXS multigroup nu-scatter matrix.

    flux_xs : OrderedDict openmc.MGXS.Flux
        MGXS multigroup flux.

    k_eff_0 : float
        The initial eigenvalue.

    Methods
    -------
    - initialize_xs()
    take_outer_step()
    take_inner_step()
    solve()
    - extract_xs()
    2 normalize_flux()
    broadcast_to_all()
    broadcast_to_one()
    - compute_shape()
    integrate_precursor_conc()
    3 compute_initial_precursor_conc()
    1 compute_power()
    construct_A()
    construct_M()
    construct_AM()
    interpolate_xs()

    To Do
    -----
    1) Create getters and setters for all attributes
    2) Create method to generate initialize xs
    3) Create method to compute flux
    4) Create method to compute initial precursor concentrations
    5) Create method to compute the initial power

    """

    def __init__(self):

        # Initialize Solver class attributes
        self._mesh = None
        self._geometry = None
        self._settings_file = None
        self._materials_file = None
        self._executor = openmc.Executor()
        self._statepoint = None
        self._summary = None
        self._clock = None
        self._energy_groups = None
        self._A = None
        self._M = None
        self._AM = None
        self._flux = None
        self._amplitude = None
        self._shape = None
        self._source = None
        self._power = None
        self._precursor_conc = None
        self._sigma_a = None
        self._nu_sigma_f = None
        self._kappa_sigma_f = None
        self._dif_coef = None
        self._beta = None
        self._chi_prompt = None
        self._chi_delayed = None
        self._velocity = None
        self._nu_sigma_s = None
        self._flux_xs = None
        self._decay_constants = None
        self._k_eff_0 = None

    @property
    def mesh(self):
        return self._mesh

    @property
    def geometry(self):
        return self._geometry

    @property
    def settings_file(self):
        return self._settings_file

    @property
    def materials_file(self):
        return self._materials_file

    @property
    def executor(self):
        return self._executor

    @property
    def statepoint(self):
        return self._statepoint

    @property
    def summary(self):
        return self._summary

    @property
    def clock(self):
        return self._clock

    @property
    def energy_groups(self):
        return self._energy_groups

    @property
    def A(self):
        return self._A

    @property
    def M(self):
        return self._M

    @property
    def AM(self):
        return self._AM

    @property
    def flux(self):
        return self._flux

    @property
    def amplitude(self):
        return self._amplitude

    @property
    def shape(self):
        return self._shape

    @property
    def source(self):
        return self._source

    @property
    def power(self):
        return self._power

    @property
    def precursor_conc(self):
        return self._precursor_conc

    @property
    def sigma_a(self):
        return self._sigma_a

    @property
    def nu_sigma_f(self):
        return self._nu_sigma_f

    @property
    def kappa_sigma_f(self):
        return self._kappa_sigma_f

    @property
    def dif_coef(self):
        return self._dif_coef

    @property
    def beta(self):
        return self._beta

    @property
    def chi_prompt(self):
        return self._chi_prompt

    @property
    def chi_delayed(self):
        return self._chi_delayed

    @property
    def velocity(self):
        return self._velocity

    @property
    def nu_sigma_s(self):
        return self._nu_sigma_s

    @property
    def flux_xs(self):
        return self._flux_xs

    @property
    def decay_constants(self):
        return self._decay_constants

    @property
    def k_eff_0(self):
        return self._k_eff_0

    @mesh.setter
    def mesh(self, mesh):
        self._mesh = mesh

    @geometry.setter
    def geometry(self, geometry):
        self._geometry = geometry

    @settings_file.setter
    def settings_file(self, settings_file):
        self._settings_file = settings_file

    @materials_file.setter
    def materials_file(self, materials_file):
        self._materials_file = materials_file

    @executor.setter
    def executor(self, exectuor):
        self._executor = executor

    @statepoint.setter
    def statepoint(self, statepoint):
        self._statepoint = statepoint

    @summary.setter
    def summary(self, summary):
        self._summary = summary

    @clock.setter
    def clock(self, clock):
        self._clock = clock

    @energy_groups.setter
    def energy_groups(self, energy_groups):
        self._energy_groups = energy_groups

        # Initialize the arrays
        ng = energy_groups.num_groups
        self._flux = np.zeros(ng)
        self._amplitude = np.zeros(ng)
        self._shape = np.zeros(ng)
        self._source = np.zeros(ng)
        self._power = np.zeros(ng)

    @A.setter
    def A(self):
        self._A = A

    @M.setter
    def M(self, M):
        self._M = M

    @AM.setter
    def AM(self, AM):
        self._AM = AM

    @flux.setter
    def flux(self, flux):
        self._flux = flux

    @amplitude.setter
    def amplitude(self, amplitude):
        self._amplitude = amplitude

    @shape.setter
    def shape(self, shape):
        self._shape = shape

    @source.setter
    def source(self, source):
        self._source = source

    @power.setter
    def power(self, power):
        self._power = power

    @precursor_conc.setter
    def precursor_conc(self, precursor_conc):
        self._precursor_conc = precursor_conc

    @sigma_a.setter
    def sigma_a(self, sigma_a):
        self._sigma_a = sigma_a

    @nu_sigma_f.setter
    def nu_sigma_f(self, nu_sigma_f):
        self._nu_sigma_f = nu_sigma_f

    @kappa_sigma_f.setter
    def kappa_sigma_f(self, kappa_sigma_f):
        self._kappa_sigma_f = kappa_sigma_f

    @dif_coef.setter
    def dif_coef(self, dif_coef):
        self._dif_coef = dif_coef

    @beta.setter
    def beta(self, beta):
        self._beta = beta

    @chi_prompt.setter
    def chi_prompt(self, chi_prompt):
        self._chi_prompt = chi_prompt

    @chi_delayed.setter
    def chi_delayed(self, chi_delayed):
        self._chi_delayed

    @velocity.setter
    def velocity(self, velocity):
        self._velocity = velocity

    @nu_sigma_s.setter
    def nu_sigma_s(self, nu_sigma_s):
        self._nu_sigma_s = nu_sigma_s

    @flux_xs.setter
    def flux_xs(self, flux_xs):
        self._flux_xs = flux_xs

    @decay_constants.setter
    def decay_constants(self, decay_constants):
        self._decay_constants = decay_constants

    @k_eff_0.setter
    def k_eff_0(self, k_eff_0):
        self._k_eff_0 = k_eff_0

    def initialize_xs(self):
        """Initialize all the tallies for the problem.

        """

        self._sigma_a     = {}
        self._nu_sigma_f  = {}
        self._kappa_sigma_f  = {}
        self._dif_coef    = {}
        self._beta        = {}
        self._chi_prompt  = {}
        self._chi_delayed = {}
        self._velocity    = {}
        self._nu_sigma_s  = {}
        self._flux_xs     = {}
        self._precursor_conc = {}

        self._decay_constants = openmc.Tally(name='decay constants')
        self._decay_constants._derived = True
        self._decay_constants.num_score_bins = 1
        self._decay_constants.add_score('None')
        self._decay_constants._mean = np.array([0.012467, 0.028292, 0.042524,\
                                                0.133042, 0.292467, 0.666488,\
                                                1.634781, 3.554600])
        self._decay_constants._mean = np.reshape(self._decay_constants._mean, (8,1,1))
        self._decay_constants._std_dev = np.array([0 for i in range(8)])
        self._decay_constants._std_dev = np.reshape(self._decay_constants._std_dev, (8,1,1))
        self._decay_constants.estimator = 'analog'
        self._decay_constants.add_filter(openmc.Filter('delayedgroup', range(1,9)))
        self._decay_constants._nuclides = ['total']

        # FIXME: replace domain with mesh
        # Get the cell in the geometry
        cells = self.geometry.root_universe.get_all_cells()
        cell = cells.values()[0]

        global TIME_POINTS
        for t in TIME_POINTS:
            print 'computing tallies for time: ' + t
            self._sigma_a[t]    = openmc.mgxs.AbsorptionXS(name='sigma a',
                domain=cell, domain_type='cell', groups=self.energy_groups)
            self._nu_sigma_f[t] = openmc.mgxs.NuFissionXS(name='nu sigma f',
                domain=cell, domain_type='cell', groups=self.energy_groups)
            self._kappa_sigma_f[t] = openmc.mgxs.KappaFissionXS(name='kappa fission',
                domain=cell, domain_type='cell', groups=self.energy_groups)
            self._dif_coef[t] = openmc.mgxs.DiffusionCoefficient(name='dif coef',
                domain=cell, domain_type='cell', groups=self.energy_groups)
            self._beta[t] = openmc.mgxs.Beta(name='beta',
                domain=cell, domain_type='cell', groups=self.energy_groups)
            self._chi_prompt[t] = openmc.mgxs.ChiPrompt(name='chi prompt',
                domain=cell, domain_type='cell', groups=self.energy_groups)
            self._chi_delayed[t] = openmc.mgxs.ChiDelayed(name='chi delayed',
                domain=cell, domain_type='cell', groups=self.energy_groups)
            self._velocity[t] = openmc.mgxs.Velocity(name='velocity',
                domain=cell, domain_type='cell', groups=self.energy_groups)
            self._nu_sigma_s[t] = openmc.mgxs.NuScatterMatrixXS(name='nu scatter',
                domain=cell, domain_type='cell', groups=self.energy_groups)
            self._flux_xs[t] = openmc.mgxs.Flux(name='flux',
                domain=cell, domain_type='cell', groups=self.energy_groups)


    def generate_tallies_file(self, time):
        """Initialize the tallies file.

        """

        tallies_file = openmc.TalliesFile()

        # Add absorption tallies to the tallies file
        for tally in self._sigma_a[time].tallies.values():
            tallies_file.add_tally(tally, merge=True)

        # Add nu-sigma-f tallies to the tallies file
        for tally in self._nu_sigma_f[time].tallies.values():
            tallies_file.add_tally(tally, merge=True)

        # Add kappa-sigma-f tallies to the tallies file
        for tally in self._kappa_sigma_f[time].tallies.values():
            tallies_file.add_tally(tally, merge=True)

        # Add dif-coef tallies to the tallies file
        for tally in self._dif_coef[time].tallies.values():
            tallies_file.add_tally(tally, merge=True)

        # Add beta tallies to the tallies file
        for tally in self._beta[time].tallies.values():
            tallies_file.add_tally(tally, merge=True)

        # Add chi prompt tallies to the tallies file
        for tally in self._chi_prompt[time].tallies.values():
            tallies_file.add_tally(tally, merge=True)

        # Add chi delayed tallies to the tallies file
        for tally in self._chi_delayed[time].tallies.values():
            tallies_file.add_tally(tally, merge=False)

        # Add velocity tallies to the tallies file
        for tally in self._velocity[time].tallies.values():
            tallies_file.add_tally(tally, merge=True)

        # Add nu-sigma-s tallies to the tallies file
        for tally in self._nu_sigma_s[time].tallies.values():
            tallies_file.add_tally(tally, merge=True)

        # Add flux tallies to the tallies file
        for tally in self._flux_xs[time].tallies.values():
            tallies_file.add_tally(tally, merge=True)

        # Export to "tallies.xml"
        tallies_file.export_to_xml()

    def extract_xs(self, time):

        filename = 'statepoint.' + str(self.settings_file.batches) + '.h5'
        self.statepoint = openmc.StatePoint(filename)
        self.summary = openmc.Summary('summary.h5')
        self.statepoint.link_with_summary(self.summary)

        # load xs from statepoint
        self._sigma_a[time].load_from_statepoint(self.statepoint)
        self._nu_sigma_f[time].load_from_statepoint(self.statepoint)
        self._kappa_sigma_f[time].load_from_statepoint(self.statepoint)
        self._dif_coef[time].load_from_statepoint(self.statepoint)
        self._beta[time].load_from_statepoint(self.statepoint)
        self._chi_prompt[time].load_from_statepoint(self.statepoint)
        self._chi_delayed[time].load_from_statepoint(self.statepoint)
        self._velocity[time].load_from_statepoint(self.statepoint)
        self._nu_sigma_s[time].load_from_statepoint(self.statepoint)
        self._flux_xs[time].load_from_statepoint(self.statepoint)
        self.k_eff_0 = self.statepoint.k_combined[0]

        # compute the xs
        self._sigma_a[time].compute_xs()
        self._nu_sigma_f[time].compute_xs()
        self._kappa_sigma_f[time].compute_xs()
        self._dif_coef[time].compute_xs()
        self._beta[time].compute_xs()
        self._chi_prompt[time].compute_xs()
        self._chi_delayed[time].compute_xs()
        self._velocity[time].compute_xs()
        self._nu_sigma_s[time].compute_xs()
        self._flux_xs[time].compute_xs()

        # extract the flux
        #for g in range(self._energy_groups.num_groups):
        #    self._flux[g] = self._flux_xs[time].xs_tally.mean[g][0][0]

    def print_xs(self, time):

        # print the xs
        self._sigma_a[time].print_xs()
        self._nu_sigma_f[time].print_xs()
        self._kappa_sigma_f[time].print_xs()
        self._dif_coef[time].print_xs()
        self._beta[time].print_xs()
        self._chi_prompt[time].print_xs()
        self._chi_delayed[time].print_xs()
        self._velocity[time].print_xs()
        self._nu_sigma_s[time].print_xs()
        self._flux_xs[time].print_xs()

    def compute_shape(self, time):

        geometry_file = openmc.GeometryFile()
        geometry_file.geometry = self.geometry

        # Create the xml files
        self._materials_file.export_to_xml()
        geometry_file.export_to_xml()
        self._settings_file.export_to_xml()
        self.generate_tallies_file(time)

        # Run OpenMC
        self.executor.run_simulation(mpi_procs=4)

    def compute_power(self, time):

        self.power = self._kappa_sigma_f[time].xs_tally * \
                self._flux_xs[time].xs_tally

    def compute_initial_precursor_conc(self, time):

        self.precursor_conc[time] = self.nu_sigma_f[time].xs_tally \
                                    * self.flux_xs[time].xs_tally
        self.precursor_conc[time] = self.precursor_conc[time].\
                                    summation(filter_type='energy', remove_filter=True)
        beta = self.beta[time].xs_tally.\
               summation(filter_type='energy', remove_filter=True)
        self.precursor_conc[time] = beta \
                                    * self.precursor_conc[time]

        self.precursor_conc[time] = self.precursor_conc[time] / self.k_eff_0

        print self.precursor_conc[time]
        print self._decay_constants

        self.precursor_conc[time] = self.precursor_conc[time] / self._decay_constants
        print self.precursor_conc[time]


    def compute_forward_flux(self, time, time_next):

        dt_v = self.clock.dt_outer * self.velocity

        fission_rate = self.nu_sigma_f[time].xs_tally \
                                    * self.flux_xs[time].xs_tally

        fission_rate = fission_rate.summation(filter_type='energy', remove_filter=True)

        self.flux_xs[time_next] = [self.flux_xs[time] + (1 - self.beta) \
                                   / self.k_eff_0 * fission_rate + \

        self.precursor_conc[time] = self.nu_sigma_f[time].xs_tally \
                                    * self.flux_xs[time].xs_tally
        self.precursor_conc[time] = self.precursor_conc[time].\
                                    summation(filter_type='energy', remove_filter=True)
        beta = self.beta[time].xs_tally.\
               summation(filter_type='energy', remove_filter=True)
        self.precursor_conc[time] = beta \
                                    * self.precursor_conc[time]

        self.precursor_conc[time] = self.precursor_conc[time] / self.k_eff_0

        print self.precursor_conc[time]
        print self._decay_constants

        inv_decay_constants = 1.0 / self._decay_constants
        inv_decay_constants.name = 'inverse decay constants'
        self.precursor_conc[time] = inv_decay_constants * self.precursor_conc[time]
        print self.precursor_conc[time]
