from collections import OrderedDict
import copy
import os
import pathlib
from shutil import copyfile, move
import subprocess
import sys
import time

import h5py
import numpy as np
from scipy.sparse.linalg import spsolve

import openmc
import openmc.checkvalue as cv
import openmc.kinetics
from openmc.kinetics.clock import TimePoints
import openmc.mgxs


class Solver:
    """Solver to propagate the neutron flux and power forward in time.

    Parameters
    ----------
    directory : str
        A directory to save the transient simulation data.

    Attributes
    ----------
    directory : pathlib.Path
        A path to the directory where transient simulation data is saved.
    shape_mesh : openmc.RegularMesh
        Mesh by which shape is computed on.
    unity_mesh : openmc.RegularMesh
        Mesh with one cell convering the entire geometry.
    amplitude_mesh : openmc.RegularMesh
        Mesh by which amplitude is computed on.
    tally_mesh : openmc.RegularMesh
        Mesh by which to tally currents
    geometry : openmc.geometry.Geometry
        Geometry which describes the problem being solved.
    settings : openmc.settings.Settings
        Settings file describing the general settings for each simulation.
    materials : openmc.materials.Materials
        Materials file containing the materials info for each simulation.
    transient : OrderedDict()
        Ordered dictionary describing the material changes during the transient.
    mgxs_lib : openmc.materials.MGXSLibrary
        MGXS Library file containing the multi-group xs for mg Monte Carlo.
    clock : openmc.kinetics.Clock
        Clock object.
    one_group : openmc.mgxs.groups.EnergyGroups
        EnergyGroups which specifies the a one-energy-group structure.
    energy_groups : openmc.mgxs.groups.EnergyGroups
        EnergyGroups which specifies the energy groups structure.
    fine_groups : openmc.mgxs.groups.EnergyGroups
        EnergyGroups used to tally the transport cross section that will be
        condensed to get the diffusion coefficients in the coarse group
        structure.
    tally_groups : openmc.mgxs.groups.EnergyGroups
        EnergyGroups used in all tallies except the diffusion coefficients
        (see fine_groups).
    initial_power : float
        The initial core power in [W].
    k_crit : float
        The initial eigenvalue.
    run_kwargs : dict
        Keyword arguments passed to :func:`openmc.run`.
    chi_delayed_by_delayed_group : bool
        Whether to use delayed groups in representing chi-delayed.
    chi_delayed_by_mesh : bool
        Whether to use a mesh in representing chi-delayed.
    use_agd : bool
        Whether to use artificial grid diffusion.
    use_pcmfd : bool
        Whether to use p-CMFD.
    num_delayed_groups : int
        The number of delayed neutron precursor groups.
    states : OrderedDict of openmc.kinetics.State
        States of the problem.
    use_pregenerated_sps : bool
        Whether to use pregenerated statepoint files.
    core_volume : float
        The core volume in [cm^3] used to normalize the initial power.
    log_file_name : str
        Log file name (excluding directory prefix).
    outer_tolerance : float
        Tolerance on the residual when converging outer time steps.
    method : string
        Approximation made for time derivatives in the transient
        scheme. 'adiabatic' allows the use of instantaneous eigenstates
        to approximate the transient.
    min_outer_iters : int
        Minimum number of outer iterations to take.

    """

    def __init__(self, directory='.'):

        # Initialize Solver class attributes
        self.directory = directory
        self._shape_mesh = None
        self._amplitude_mesh = None
        self._unity_mesh = None
        self._tally_mesh = None
        self._geometry = None
        self._settings = None
        self._materials = None
        self._transient = None
        self._mgxs_lib = None
        self._clock = None
        self._one_group = None
        self._energy_groups = None
        self._tally_groups = None
        self._fine_groups = None
        self._initial_power = 1.
        self._k_crit = 1.0
        self._run_kwargs = None
        self._chi_delayed_by_delayed_group = False
        self._chi_delayed_by_mesh = False
        self._num_delayed_groups = 6
        self._states = OrderedDict()
        self._use_pregenerated_sps = False
        self._core_volume = 1.
        self._log_file_name = 'log_file.h5'
        self._outer_tolerance = 1.e-6
        self._method = 'adiabatic'
        self._use_agd = False
        self._use_pcmfd = False
        self._min_outer_iters = 2

    @property
    def directory(self):
        return pathlib.Path(self._directory)

    @property
    def amplitude_mesh(self):
        return self._amplitude_mesh

    @property
    def shape_mesh(self):
        return self._shape_mesh

    @property
    def unity_mesh(self):
        return self._unity_mesh

    @property
    def tally_mesh(self):
        return self._tally_mesh

    @property
    def geometry(self):
        return self._geometry

    @property
    def settings(self):
        return self._settings

    @property
    def materials(self):
        return self._materials

    @property
    def transient(self):
        return self._transient

    @property
    def mgxs_lib(self):
        return self._mgxs_lib

    @property
    def clock(self):
        return self._clock

    @property
    def one_group(self):
        return self._one_group

    @property
    def energy_groups(self):
        return self._energy_groups

    @property
    def fine_groups(self):
        return self._fine_groups

    @property
    def tally_groups(self):
        return self._tally_groups

    @property
    def initial_power(self):
        return self._initial_power

    @property
    def k_crit(self):
        return self._k_crit

    @property
    def run_kwargs(self):
        return self._run_kwargs

    @property
    def chi_delayed_by_delayed_group(self):
        return self._chi_delayed_by_delayed_group

    @property
    def use_agd(self):
        return self._use_agd

    @property
    def use_pcmfd(self):
        return self._use_pcmfd

    @property
    def chi_delayed_by_mesh(self):
        return self._chi_delayed_by_mesh

    @property
    def num_delayed_groups(self):
        return self._num_delayed_groups

    @property
    def states(self):
        return self._states

    @property
    def use_pregenerated_sps(self):
        return self._use_pregenerated_sps

    @property
    def core_volume(self):
        return self._core_volume

    @property
    def log_file_name(self):
        return self._log_file_name

    @property
    def outer_tolerance(self):
        return self._outer_tolerance

    @property
    def method(self):
        return self._method

    @property
    def min_outer_iters(self):
        return self._min_outer_iters

    @directory.setter
    def directory(self, directory):
        cv.check_type('directory', directory, str)
        self._directory = directory

    @amplitude_mesh.setter
    def amplitude_mesh(self, mesh):

        cv.check_type('amplitude mesh', mesh, openmc.RegularMesh)
        self._amplitude_mesh = mesh
        unity_mesh = openmc.RegularMesh()
        unity_mesh.dimension = (1, 1, 1)
        unity_mesh.lower_left  = mesh.lower_left
        if mesh.width is not None:
            unity_mesh.width = [i*j for i,j in zip(mesh.dimension, mesh.width)]
        else:
            unity_mesh.width = [i-j for i,j in zip(mesh.upper_right, mesh.lower_left)]
        self._unity_mesh = unity_mesh

        # Set the power mesh to the shape mesh if it has not be set
        if self.shape_mesh is None:
            self.shape_mesh = mesh

        if self.tally_mesh is None:
            self.tally_mesh = mesh

    @shape_mesh.setter
    def shape_mesh(self, mesh):
        cv.check_type('shape mesh', mesh, openmc.RegularMesh)
        self._shape_mesh = mesh

    @tally_mesh.setter
    def tally_mesh(self, mesh):
        cv.check_type('tally mesh', mesh, openmc.RegularMesh)
        self._tally_mesh = mesh

    @geometry.setter
    def geometry(self, geometry):
        cv.check_type('geometry', geometry, openmc.geometry.Geometry)
        self._geometry = geometry

    @settings.setter
    def settings(self, settings):
        cv.check_type('settings', settings, openmc.settings.Settings)
        self._settings = settings

    @materials.setter
    def materials(self, materials):
        cv.check_type('materials', materials, openmc.Materials)
        self._materials = materials

    @transient.setter
    def transient(self, transient):
        cv.check_type('transients', transient, dict)
        self._transient = transient

    @mgxs_lib.setter
    def mgxs_lib(self, mgxs_lib):
        cv.check_type('mgxs library', mgxs_lib, openmc.mgxs_library.MGXSLibrary)
        self._mgxs_lib = mgxs_lib

    @clock.setter
    def clock(self, clock):
        cv.check_type('clock', clock, openmc.kinetics.Clock)
        self._clock = copy.deepcopy(clock)

    @one_group.setter
    def one_group(self, one_group):
        cv.check_type('one-group structure', one_group, openmc.mgxs.groups.EnergyGroups)
        self._one_group = one_group

    @energy_groups.setter
    def energy_groups(self, energy_groups):
        cv.check_type('energy-group structure', energy_groups, openmc.mgxs.groups.EnergyGroups)
        self._energy_groups = energy_groups

    @fine_groups.setter
    def fine_groups(self, fine_groups):
        cv.check_type('fine-group structure', fine_groups, openmc.mgxs.groups.EnergyGroups)
        self._fine_groups = fine_groups

    @tally_groups.setter
    def tally_groups(self, tally_groups):
        cv.check_type('tally group structure', tally_groups, openmc.mgxs.groups.EnergyGroups)
        self._tally_groups = tally_groups

    @initial_power.setter
    def initial_power(self, initial_power):
        self._initial_power = initial_power

    @k_crit.setter
    def k_crit(self, k_crit):
        self._k_crit = k_crit

    @run_kwargs.setter
    def run_kwargs(self, run_kwargs):
        self._run_kwargs = run_kwargs

    @chi_delayed_by_delayed_group.setter
    def chi_delayed_by_delayed_group(self, chi_delayed_by_delayed_group):
        cv.check_type('chi delayed by group boolean', chi_delayed_by_delayed_group, bool)
        self._chi_delayed_by_delayed_group = chi_delayed_by_delayed_group

    @chi_delayed_by_mesh.setter
    def chi_delayed_by_mesh(self, chi_delayed_by_mesh):
        cv.check_type('chi delayed by mesh boolean', chi_delayed_by_mesh, bool)
        self._chi_delayed_by_mesh = chi_delayed_by_mesh

    @use_agd.setter
    def use_agd(self, use_agd):
        cv.check_type('agd boolean', use_agd, bool)
        self._use_agd = use_agd

    @use_pcmfd.setter
    def use_pcmfd(self, use_pcmfd):
        cv.check_type('pcmfd boolean', use_pcmfd, bool)
        self._use_pcmfd = use_pcmfd

    @num_delayed_groups.setter
    def num_delayed_groups(self, num_delayed_groups):
        cv.check_type('number of delayed groups', num_delayed_groups, int)
        self._num_delayed_groups = num_delayed_groups

    @states.setter
    def states(self, states):
        self._states = states

    @use_pregenerated_sps.setter
    def use_pregenerated_sps(self, use_pregenerated_sps):
        cv.check_type('pregenerated sps boolean', use_pregenerated_sps, bool)
        self._use_pregenerated_sps = use_pregenerated_sps

    @core_volume.setter
    def core_volume(self, core_volume):
        cv.check_type('core volume', core_volume, float)
        self._core_volume = core_volume

    @log_file_name.setter
    def log_file_name(self, name):
        cv.check_type('log file name', name, str)
        self._log_file_name = name

    @outer_tolerance.setter
    def outer_tolerance(self, tolerance):
        cv.check_type('outer tolerance', tolerance, float)
        self._outer_tolerance = tolerance

    @method.setter
    def method(self, method):
        # The omega method solves equations (2.29) and (2.30) in Shaner's thesis.
        # The adiabatic method sets the frequencies = 0 which simplifies those same
        # equations to the instantaneous eigenstate. 
        cv.check_value('method', method, ['adiabatic', 'omega'])
        self._method = method

    @min_outer_iters.setter
    def min_outer_iters(self, iters):
        cv.check_type('minimum outer iterations', iters, int)
        self._min_outer_iters = iters

    @property
    def ng(self):
        return self.energy_groups.num_groups

    @property
    def log_file(self):
        log_file = self.log_file_name.replace(' ', '-')
        log_file = self.directory / log_file
        return log_file

    def _create_log_file(self):
        """Write a log file that stores information about the transient.

        """

        with h5py.File(self.log_file, 'w') as f:

            f.require_group('shape')
            f['shape'].attrs['id'] = self.shape_mesh.id
            f['shape'].attrs['name'] = self.shape_mesh.name
            f['shape'].attrs['dimension'] = self.shape_mesh.dimension
            f['shape'].attrs['lower_left'] = self.shape_mesh.lower_left
            if self.shape_mesh.width:
                f['shape'].attrs['width'] = self.shape_mesh.width
            else:
                f['shape'].attrs['width'] = [(i-j)/k for i,j,k in zip(
                    self.shape_mesh.upper_right, 
                    self.shape_mesh.lower_left, 
                    self.shape_mesh.dimension)]

            f.require_group('amplitude')
            f['amplitude'].attrs['id'] = self.amplitude_mesh.id
            f['amplitude'].attrs['name'] = self.amplitude_mesh.name
            f['amplitude'].attrs['dimension'] = self.amplitude_mesh.dimension
            f['amplitude'].attrs['lower_left'] = self.amplitude_mesh.lower_left
            if self.amplitude_mesh.width:
                f['amplitude'].attrs['width'] = self.amplitude_mesh.width
            else:
                f['amplitude'].attrs['width'] = [(i-j)/k for i,j,k in zip(
                    self.amplitude_mesh.upper_right, 
                    self.amplitude_mesh.lower_left, 
                    self.amplitude_mesh.dimension)]

            for groups,name in \
                zip([self.one_group, self.energy_groups, self.fine_groups, self.tally_groups],
                    ['one_group', 'energy_groups', 'fine_groups', 'tally_groups']):
                f.require_group(name)
                f[name].attrs['group_edges'] = groups.group_edges
                f[name].attrs['num_groups'] = groups.num_groups

            f.attrs['num_delayed_groups'] = self.num_delayed_groups
            f.attrs['chi_delayed_by_delayed_group'] \
                = self.chi_delayed_by_delayed_group
            f.attrs['chi_delayed_by_mesh'] = self.chi_delayed_by_mesh
            f.attrs['use_agd'] = self.use_agd
            f.attrs['use_pcmfd'] = self.use_pcmfd
            f.attrs['num_delayed_groups'] = self.num_delayed_groups
            f.attrs['num_outer_time_steps'] = self.num_outer_time_steps
            f.attrs['use_pregenerated_sps'] = self.use_pregenerated_sps
            f.attrs['core_volume'] = self.core_volume
            f.attrs['k_crit'] = self.k_crit
            f.attrs['method'] = self.method
            f.attrs['pnl'] = self.states['START'].pnl
            f.require_group('clock')
            f['clock'].attrs['dt_outer'] = self.clock.dt_outer
            f['clock'].attrs['dt_inner'] = self.clock.dt_inner
            f.require_group('OUTER_STEPS')
            f.require_group('INNER_STEPS')
            f.create_dataset('beta', data=self.states['START'].beta())
            f.create_dataset('beta_shape', data=self.states['START'].beta_shape())
            f.create_dataset('beta_eff', data=self.states['START'].beta_eff)
            f.create_dataset('power', data=self.states['START'].power)
            f.create_dataset('decay_rate', data=self.states['START'].decay_rate)
            f.create_dataset('inverse_velocity', data=self.states['START'].inverse_velocity)
            total = self.states['START'].absorption + self.states['START'].outscatter
            f.create_dataset('total', data=total)

    def _run_openmc(self, time_point, method='adiabatic'):
        """Run a simulation at a given time point.

        Parameters
        ----------
        time_point : str
            Refers to the present time step of the simulation.
        method : str
            The approximation to the time-derivatives in the Monte Carlo simulation.
            'adiabatic' approximates the time-derivatives as 0.

        """

        self._setup_openmc(time_point, method)

        # Names of the statepoint and summary files
        sp_old_name = self.directory / 'statepoint.{}.h5'.format(self.settings.batches)
        sp_new_name = self.directory / 'statepoint_{:.6f}_sec.{}.h5'\
                      .format(self.clock.times[time_point], self.settings.batches)
        sum_old_name = self.directory / 'summary.h5'
        sum_new_name = self.directory / 'summary_{:.6f}_sec.h5'.format(self.clock.times[time_point])
        old_source = self.directory / 'source.{}.h5'.format(self.settings.batches)
        new_source = self.directory / 'source.h5'

        # Run OpenMC
        if not self.use_pregenerated_sps:
            openmc.run(cwd=self.directory, **self.run_kwargs)

            # Rename the statepoint and summary files
            copyfile(sp_old_name, sp_new_name)
            copyfile(sum_old_name, sum_new_name)
            move(old_source, new_source)

        # Load the summary and statepoint files
        summary_file = openmc.Summary(str(sum_new_name))
        statepoint_file = openmc.StatePoint(str(sp_new_name), False)
        statepoint_file.link_with_summary(summary_file)

        # Load mgxs library
        for mgxs in self.states[time_point].mgxs_lib.values():
            mgxs.load_from_statepoint(statepoint_file)

        self.states[time_point]._load_mgxs()

    def _create_state(self, time_point):
        """Modify state as needed to reflect current values in solver.

        Parameters
        ----------
        time_point : str
            Refers to the present time step of the simulation.

        """

        if time_point in ('START', 'END', 'FORWARD_OUTER', 'PREVIOUS_OUTER'):
            state = openmc.kinetics.OuterState(self.states)
            if time_point in ('FORWARD_OUTER', 'PREVIOUS_OUTER'):
                state.method = self.method
        else:
            state = openmc.kinetics.InnerState(self.states)
            state.method = self.method

        state.solver = self
        state.time_point = time_point

        if self.tally_groups is not None:
            state.tally_groups = self.tally_groups
        else:
            state.tally_groups = self.energy_groups

        self.states[time_point] = state

    def _compute_initial_flux(self):
        """Run a simulation at the starting state.

        """

        # Run OpenMC to obtain XS
        self._run_openmc('START')

        # Get the START state
        state = self.states['START']
        state.shape = state.flux_tallied

        # Compute the initial eigenvalue
        state.amplitude, self.k_crit = openmc.kinetics.compute_eigenvalue(state.destruction_matrix,
                                                                          state.production_matrix,
                                                                          np.ones(state.amplitude_nxyz * self.ng),
                                                                          tolerance = 1.e-6,
                                                                          max_eig_iterations = 10000)

        # Compute the initial adjoint eigenvalue
        state.adjoint_flux, k_adjoint = openmc.kinetics.compute_eigenvalue\
            (state.destruction_matrix.transpose(), state.production_matrix.transpose(),
             np.ones(state.amplitude_nxyz * self.ng), tolerance = 1.e-3, max_eig_iterations= 10000)

        # Compute the power and normalize the amplitude
        state.amplitude    /= np.average(state.amplitude.flatten())
        state.adjoint_flux /= np.average(state.adjoint_flux.flatten())
        norm_factor        = self.initial_power / state.core_power_density
        state.shape       *= norm_factor

        # Compute the initial precursor concentration
        state._compute_initial_precursor_concentration()

        # Copy data to all other states
        for time_point in TimePoints:
            if time_point.name != 'START':
                self._create_state(time_point.name)
                self._copy_states('START', time_point.name)

        # Create hdf5 log file
        self._create_log_file()
        self.states['START'].dump_to_log_file

        # Reset the source to be the last source bank written to file
        self.settings.source   = openmc.Source(filename='source.h5')
        self.settings.batches  = self.settings.batches - self.settings.inactive + 10
        self.settings.inactive = 10
        self.settings.sourcepoint['batches'] = [self.settings.batches]

    def _copy_states(self, time_from, time_to):
        """Use old state as starting point for new state.

        Parameters
        ----------
        time_from : str
            Time point being copied from.
        time_to : str
            Time point being copied into.

        """

        state_from            = self.states[time_from]
        state_to              = self.states[time_to]
        state_to.amplitude    = state_from.amplitude
        state_to.adjoint_flux = state_from.adjoint_flux
        state_to.precursors   = state_from.precursors

        if time_to != 'END':
            self.clock.times[time_to] = self.clock.times[time_from]

        if time_to in ['START', 'END', 'PREVIOUS_OUTER', 'FORWARD_OUTER'] \
                and time_from in ['START', 'END', 'PREVIOUS_OUTER', 'FORWARD_OUTER']:
            state_to.shape    = state_from.shape
            state_to.mgxs_lib = state_from.mgxs_lib
            state_to._load_mgxs()

    def _take_inner_step(self, i):
        """Perform an inner time step using CMFD.

        Parameters
        ----------
        i : int
            Counter for the number of inner time steps.

        """

        # Increment clock
        times     = self.clock.times
        state_pre = self.states['PREVIOUS_INNER']
        state_fwd = self.states['FORWARD_INNER']

        # Update the values for the time step
        self._copy_states('FORWARD_INNER', 'PREVIOUS_INNER')

        # Increment forward in time
        times['FORWARD_INNER'] += self.clock.dt_inner

        # Form the source
        source = state_fwd.transient_source
        matrix = state_fwd.transient_matrix

        # Compute the amplitude at the FORWARD_IN time step
        state_fwd.amplitude = spsolve(matrix, source)

        # Propagate the precursors
        state_fwd.propagate_precursors

        # Dump data at FORWARD_INNER state to log file
        state_fwd.dump_to_log_file

        # Save the core power at FORWARD_IN
        print('t: {0:1.5f} s, P: {1:1.3e} W/cm^3, rho: {2:+1.3f} pcm'
              ', beta_eff: {3:1.5e}, pnl: {4:1.3e} s'.\
                  format(times['FORWARD_INNER'], state_fwd.core_power_density,
                         state_fwd.reactivity * 1.e5,
                         state_fwd.beta_eff.sum(), state_fwd.pnl))

    def _take_outer_step(self, outer_step):
        """Perform a complete outer time step, first by running the next outer 
        time point Monte Carlo simulation, then calling _take_inner_time_step()
        for all the inner time steps. This is done in a loop until the residual
        is below tolerance.

        Parameters
        ----------
        outer_step : int
            Counter for the number of outer time steps.

        """

        # Increment clock
        times     = self.clock.times
        state_fwd = self.states['FORWARD_OUTER']
        state_pre = self.states['PREVIOUS_OUTER']
        iteration = 0

        # Save the old power
        power_old = state_fwd.power

        while True:

            if iteration == 0:
                self.clock.times['FORWARD_OUTER'] = self.clock.times['FORWARD_OUTER'] +\
                    self.clock.dt_outer
                residual = np.inf
            else:

                # Take inner steps
                for i in range(self.num_inner_time_steps):
                    self._take_inner_step(i)

                # Copy the shape, amp, and and precursors to FORWARD_OUTER
                self._copy_states('FORWARD_INNER', 'FORWARD_OUTER')

                new_power = state_fwd.power
                residual_array = (power_old - new_power) / new_power
                residual_array = openmc.kinetics.nan_inf_to_zero(residual_array)
                num_fissile_regions = np.sum(power_old > 0.)
                residual = np.sqrt((residual_array**2).sum() / num_fissile_regions)
                power_old = new_power

            if residual < self.outer_tolerance and iteration >= self.min_outer_iters:
                print('  CONVERGED OUTER residual {}'.format(residual))
                break
            else:

                print('UNCONVERGED OUTER residual {}'.format(residual))

                # Increment outer iteration count
                iteration += 1

                # Run OpenMC on forward out state
                if iteration == 1:
                    self._run_openmc('FORWARD_OUTER')
                else:
                    self._run_openmc('FORWARD_OUTER', self.method)

                state_fwd._extract_shape()

                # Reset the inner states
                self._copy_states('PREVIOUS_OUTER', 'FORWARD_INNER')
                self._copy_states('PREVIOUS_OUTER', 'PREVIOUS_INNER')

        # Dump data at FORWARD_OUT state to log file
        state_fwd.dump_to_log_file

        # Copy the converged forward state to the previous state
        self._copy_states('FORWARD_OUTER', 'PREVIOUS_OUTER')

        # Increment the clock time step
        self.clock.outer_step += 1

    def _generate_tallies_file(self, time_point):
        """Write a file that stores information about the tallies at a given time point.
        Note these files will overwrite each other throughout the transient. 

        Parameters
        ----------
        time_point : str
            Refers to the present time step of the simulation.

        """

        # Generate a new tallies file
        tallies_file = openmc.Tallies()

        # Get the MGXS library
        mgxs_lib = self.states[time_point].mgxs_lib

        # Add the tallies to the file
        for mgxs in mgxs_lib.values():
            tallies = mgxs.tallies.values()
            for tally in tallies:
                tallies_file.append(tally, True)

        # Export the tallies file to xml
        tallies_file.export_to_xml(self.directory / 'tallies.xml')

    def solve(self):
        """Run a transient simulation based on input created in an OpenMC input file.

        """

        # Create run directory
        self.directory.mkdir(exist_ok=True)

        # Create states
        self._create_state('START')

        # Compute the initial steady state flux
        self._compute_initial_flux()

        # Solve the transient
        for i in range(self.clock.num_outer_steps):
            self._take_outer_step(i)

    def _setup_openmc(self, time_point, method):
        """Prepare a simulation at a given time point.

        Parameters
        ----------
        time_point : str
            Refers to the present time step of the simulation.
        method : str
            The approximation to the time-derivatives in the Monte Carlo simulation.
            'adiabatic' approximates the time-derivatives as 0.

        """

        # Get a fresh copy of the settings file
        settings = copy.deepcopy(self.settings)

        if self.mgxs_lib:
            self.materials.cross_sections = './mgxs.h5'
            self.mgxs_lib.export_to_hdf5(str(self.directory / 'mgxs.h5'))

        # Create MGXS
        state = self.states[time_point]

        state._initialize_mgxs()

        # Create the xml files
        self.geometry.time = self.clock.times[time_point]

        materials_list = []
        for material in self.materials:
            time = round(self.geometry.time,4)
            if settings.energy_mode == 'multi-group':
                material.set_density('macro',self.transient[material.name][time]['density'])
            else:
                material.set_density('g/cm3',self.transient[material.name][time]['density'])
            material.temperature = self.transient[material.name][time]['temperature']
            materials_list.append(material)
        self.materials = openmc.Materials(materials_list)
        
        self.geometry.export_to_xml(self.directory / 'geometry.xml')
        self.materials.export_to_xml(self.directory / 'materials.xml')
        settings.export_to_xml(self.directory / 'settings.xml')
        self._generate_tallies_file(time_point)

    @property
    def num_outer_time_steps(self):
        return int(round((self.clock.times['END'] - self.clock.times['START']) \
                             / self.clock.dt_outer))

    @property
    def num_inner_time_steps(self):
        return int(round(self.clock.dt_outer / self.clock.dt_inner))
