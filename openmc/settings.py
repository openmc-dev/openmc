from collections.abc import Iterable, Mapping, MutableSequence
from enum import Enum
import itertools
from math import ceil
from numbers import Integral, Real
from pathlib import Path
import typing  # required to prevent typing.Union namespace overwriting Union
from typing import Optional

import lxml.etree as ET

import openmc.checkvalue as cv
from openmc.stats.multivariate import MeshSpatial
from . import (RegularMesh, SourceBase, MeshSource, IndependentSource,
               VolumeCalculation, WeightWindows, WeightWindowGenerator)
from ._xml import clean_indentation, get_text, reorder_attributes
from openmc.checkvalue import PathLike
from .mesh import _read_meshes


class RunMode(Enum):
    EIGENVALUE = 'eigenvalue'
    FIXED_SOURCE = 'fixed source'
    PLOT = 'plot'
    VOLUME = 'volume'
    PARTICLE_RESTART = 'particle restart'


_RES_SCAT_METHODS = ['dbrc', 'rvs']


class Settings:
    """Settings used for an OpenMC simulation.

    Parameters
    ----------
    **kwargs : dict, optional
        Any keyword arguments are used to set attributes on the instance.

    Attributes
    ----------
    batches : int
        Number of batches to simulate
    confidence_intervals : bool
        If True, uncertainties on tally results will be reported as the
        half-width of the 95% two-sided confidence interval. If False,
        uncertainties on tally results will be reported as the sample standard
        deviation.
    create_fission_neutrons : bool
        Indicate whether fission neutrons should be created or not.
    cutoff : dict
        Dictionary defining weight cutoff, energy cutoff and time cutoff. The
        dictionary may have ten keys, 'weight', 'weight_avg', 'energy_neutron',
        'energy_photon', 'energy_electron', 'energy_positron', 'time_neutron',
        'time_photon', 'time_electron', and 'time_positron'. Value for 'weight'
        should be a float indicating weight cutoff below which particle undergo
        Russian roulette. Value for 'weight_avg' should be a float indicating
        weight assigned to particles that are not killed after Russian roulette.
        Value of energy should be a float indicating energy in eV below which
        particle type will be killed. Value of time should be a float in
        seconds. Particles will be killed exactly at the specified time.
    delayed_photon_scaling : bool
        Indicate whether to scale the fission photon yield by (EGP + EGD)/EGP
        where EGP is the energy release of prompt photons and EGD is the energy
        release of delayed photons.

        .. versionadded:: 0.12
    electron_treatment : {'led', 'ttb'}
        Whether to deposit all energy from electrons locally ('led') or create
        secondary bremsstrahlung photons ('ttb').
    energy_mode : {'continuous-energy', 'multi-group'}
        Set whether the calculation should be continuous-energy or multi-group.
    entropy_mesh : openmc.RegularMesh
        Mesh to be used to calculate Shannon entropy. If the mesh dimensions are
        not specified, OpenMC assigns a mesh such that 20 source sites per mesh
        cell are to be expected on average.
    event_based : bool
        Indicate whether to use event-based parallelism instead of the default
        history-based parallelism.

        .. versionadded:: 0.12
    generations_per_batch : int
        Number of generations per batch
    max_lost_particles : int
        Maximum number of lost particles

        .. versionadded:: 0.12
    rel_max_lost_particles : float
        Maximum number of lost particles, relative to the total number of
        particles

        .. versionadded:: 0.12
    inactive : int
        Number of inactive batches
    keff_trigger : dict
        Dictionary defining a trigger on eigenvalue. The dictionary must have
        two keys, 'type' and 'threshold'. Acceptable values corresponding to
        type are 'variance', 'std_dev', and 'rel_err'. The threshold value
        should be a float indicating the variance, standard deviation, or
        relative error used.
    log_grid_bins : int
        Number of bins for logarithmic energy grid search
    material_cell_offsets : bool
        Generate an "offset table" for material cells by default. These tables
        are necessary when a particular instance of a cell needs to be tallied.

        .. versionadded:: 0.12
    max_particles_in_flight : int
        Number of neutrons to run concurrently when using event-based
        parallelism.

        .. versionadded:: 0.12
    max_particle_events : int
        Maximum number of allowed particle events per source particle.

        .. versionadded:: 0.14.1
    max_order : None or int
        Maximum scattering order to apply globally when in multi-group mode.
    max_history_splits : int
        Maximum number of times a particle can split during a history

        .. versionadded:: 0.13
    max_tracks : int
        Maximum number of tracks written to a track file (per MPI process).

        .. versionadded:: 0.13.1
    max_write_lost_particles : int
        Maximum number of particle restart files (per MPI process) to write for
        lost particles.

        .. versionadded:: 0.14.0
    no_reduce : bool
        Indicate that all user-defined and global tallies should not be reduced
        across processes in a parallel calculation.
    output : dict
        Dictionary indicating what files to output. Acceptable keys are:

        :path: String indicating a directory where output files should be
               written
        :summary: Whether the 'summary.h5' file should be written (bool)
        :tallies: Whether the 'tallies.out' file should be written (bool)
    particles : int
        Number of particles per generation
    photon_transport : bool
        Whether to use photon transport.
    plot_seed : int
       Initial seed for randomly generated plot colors.
    ptables : bool
        Determine whether probability tables are used.
    random_ray : dict
        Options for configuring the random ray solver. Acceptable keys are:

        :distance_inactive:
            Indicates the total active distance in [cm] a ray should travel
        :distance_active:
            Indicates the total active distance in [cm] a ray should travel
        :ray_source:
            Starting ray distribution (must be uniform in space and angle) as
            specified by a :class:`openmc.SourceBase` object.

        .. versionadded:: 0.14.1
    resonance_scattering : dict
        Settings for resonance elastic scattering. Accepted keys are 'enable'
        (bool), 'method' (str), 'energy_min' (float), 'energy_max' (float), and
        'nuclides' (list). The 'method' can be set to 'dbrc' (Doppler broadening
        rejection correction) or 'rvs' (relative velocity sampling). If not
        specified, 'rvs' is the default method. The 'energy_min' and
        'energy_max' values indicate the minimum and maximum energies above and
        below which the resonance elastic scattering method is to be applied.
        The 'nuclides' list indicates what nuclides the method should be applied
        to. In its absence, the method will be applied to all nuclides with 0 K
        elastic scattering data present.
    run_mode : {'eigenvalue', 'fixed source', 'plot', 'volume', 'particle restart'}
        The type of calculation to perform (default is 'eigenvalue')
    seed : int
        Seed for the linear congruential pseudorandom number generator
    source : Iterable of openmc.SourceBase
        Distribution of source sites in space, angle, and energy
    sourcepoint : dict
        Options for writing source points. Acceptable keys are:

        :batches: list of batches at which to write source
        :overwrite: bool indicating whether to overwrite
        :separate: bool indicating whether the source should be written as a
                   separate file
        :write: bool indicating whether or not to write the source
        :mcpl: bool indicating whether to write the source as an MCPL file
    statepoint : dict
        Options for writing state points. Acceptable keys are:

        :batches: list of batches at which to write statepoint files
    surf_source_read : dict
        Options for reading surface source points. Acceptable keys are:

        :path: Path to surface source file (str).
    surf_source_write : dict
        Options for writing surface source points. Acceptable keys are:

        :surface_ids: List of surface ids at which crossing particles are to be
                   banked (int)
        :max_particles: Maximum number of particles to be banked on surfaces per
                   process (int)
        :mcpl: Output in the form of an MCPL-file (bool)
    survival_biasing : bool
        Indicate whether survival biasing is to be used
    tabular_legendre : dict
        Determines if a multi-group scattering moment kernel expanded via
        Legendre polynomials is to be converted to a tabular distribution or
        not. Accepted keys are 'enable' and 'num_points'. The value for 'enable'
        is a bool stating whether the conversion to tabular is performed; the
        value for 'num_points' sets the number of points to use in the tabular
        distribution, should 'enable' be True.
    temperature : dict
        Defines a default temperature and method for treating intermediate
        temperatures at which nuclear data doesn't exist. Accepted keys are
        'default', 'method', 'range', 'tolerance', and 'multipole'. The value
        for 'default' should be a float representing the default temperature in
        Kelvin. The value for 'method' should be 'nearest' or 'interpolation'.
        If the method is 'nearest', 'tolerance' indicates a range of temperature
        within which cross sections may be used. If the method is
        'interpolation', 'tolerance' indicates the range of temperatures outside
        of the available cross section temperatures where cross sections will
        evaluate to the nearer bound. The value for 'range' should be a pair of
        minimum and maximum temperatures which are used to indicate that cross
        sections be loaded at all temperatures within the range. 'multipole' is
        a boolean indicating whether or not the windowed multipole method should
        be used to evaluate resolved resonance cross sections.
    trace : tuple or list
        Show detailed information about a single particle, indicated by three
        integers: the batch number, generation number, and particle number
    track : tuple or list
        Specify particles for which track files should be written. Each particle
        is identified by a tuple with the batch number, generation number, and
        particle number.
    trigger_active : bool
        Indicate whether tally triggers are used
    trigger_batch_interval : int
        Number of batches in between convergence checks
    trigger_max_batches : int
        Maximum number of batches simulated. If this is set, the number of
        batches specified via ``batches`` is interpreted as the minimum number
        of batches
    ufs_mesh : openmc.RegularMesh
        Mesh to be used for redistributing source sites via the uniform fission
        site (UFS) method.
    verbosity : int
        Verbosity during simulation between 1 and 10. Verbosity levels are
        described in :ref:`verbosity`.
    volume_calculations : VolumeCalculation or iterable of VolumeCalculation
        Stochastic volume calculation specifications
    weight_windows : WeightWindows or iterable of WeightWindows
        Weight windows to use for variance reduction

        .. versionadded:: 0.13
    weight_window_checkpoints : dict
        Indicates the checkpoints for weight window split/roulettes. Valid keys
        include "collision" and "surface". Values must be of type bool.

        .. versionadded:: 0.14.0
    weight_window_generators : WeightWindowGenerator or iterable of WeightWindowGenerator
        Weight windows generation parameters to apply during simulation

        .. versionadded:: 0.14.0

    create_delayed_neutrons : bool
        Whether delayed neutrons are created in fission.

        .. versionadded:: 0.13.3
    weight_windows_on : bool
        Whether weight windows are enabled

        .. versionadded:: 0.13

    weight_windows_file: Pathlike
        Path to a weight window file to load during simulation initialization

        .. versionadded::0.14.0
    write_initial_source : bool
        Indicate whether to write the initial source distribution to file
    """

    def __init__(self, **kwargs):
        self._run_mode = RunMode.EIGENVALUE
        self._batches = None
        self._generations_per_batch = None
        self._inactive = None
        self._max_lost_particles = None
        self._rel_max_lost_particles = None
        self._max_write_lost_particles = None
        self._particles = None
        self._keff_trigger = None

        # Energy mode subelement
        self._energy_mode = None
        self._max_order = None

        # Source subelement
        self._source = cv.CheckedList(SourceBase, 'source distributions')

        self._confidence_intervals = None
        self._electron_treatment = None
        self._photon_transport = None
        self._plot_seed = None
        self._ptables = None
        self._seed = None
        self._survival_biasing = None

        # Shannon entropy mesh
        self._entropy_mesh = None

        # Trigger subelement
        self._trigger_active = None
        self._trigger_max_batches = None
        self._trigger_batch_interval = None

        self._output = None

        # Output options
        self._statepoint = {}
        self._sourcepoint = {}

        self._surf_source_read = {}
        self._surf_source_write = {}

        self._no_reduce = None

        self._verbosity = None

        self._trace = None
        self._track = None

        self._tabular_legendre = {}

        self._temperature = {}

        # Cutoff subelement
        self._cutoff = None

        # Uniform fission source subelement
        self._ufs_mesh = None

        self._resonance_scattering = {}
        self._volume_calculations = cv.CheckedList(
            VolumeCalculation, 'volume calculations')

        self._create_fission_neutrons = None
        self._create_delayed_neutrons = None
        self._delayed_photon_scaling = None
        self._material_cell_offsets = None
        self._log_grid_bins = None

        self._event_based = None
        self._max_particles_in_flight = None
        self._max_particle_events = None
        self._write_initial_source = None
        self._weight_windows = cv.CheckedList(WeightWindows, 'weight windows')
        self._weight_window_generators = cv.CheckedList(WeightWindowGenerator, 'weight window generators')
        self._weight_windows_on = None
        self._weight_windows_file = None
        self._weight_window_checkpoints = {}
        self._max_history_splits = None
        self._max_tracks = None

        self._random_ray = {}

        for key, value in kwargs.items():
            setattr(self, key, value)

    @property
    def run_mode(self) -> str:
        return self._run_mode.value

    @run_mode.setter
    def run_mode(self, run_mode: str):
        cv.check_value('run mode', run_mode, {x.value for x in RunMode})
        for mode in RunMode:
            if mode.value == run_mode:
                self._run_mode = mode

    @property
    def batches(self) -> int:
        return self._batches

    @batches.setter
    def batches(self, batches: int):
        cv.check_type('batches', batches, Integral)
        cv.check_greater_than('batches', batches, 0)
        self._batches = batches

    @property
    def generations_per_batch(self) -> int:
        return self._generations_per_batch

    @generations_per_batch.setter
    def generations_per_batch(self, generations_per_batch: int):
        cv.check_type('generations per patch', generations_per_batch, Integral)
        cv.check_greater_than('generations per batch', generations_per_batch, 0)
        self._generations_per_batch = generations_per_batch

    @property
    def inactive(self) -> int:
        return self._inactive

    @inactive.setter
    def inactive(self, inactive: int):
        cv.check_type('inactive batches', inactive, Integral)
        cv.check_greater_than('inactive batches', inactive, 0, True)
        self._inactive = inactive

    @property
    def max_lost_particles(self) -> int:
        return self._max_lost_particles

    @max_lost_particles.setter
    def max_lost_particles(self, max_lost_particles: int):
        cv.check_type('max_lost_particles', max_lost_particles, Integral)
        cv.check_greater_than('max_lost_particles', max_lost_particles, 0)
        self._max_lost_particles = max_lost_particles

    @property
    def rel_max_lost_particles(self) -> float:
        return self._rel_max_lost_particles

    @rel_max_lost_particles.setter
    def rel_max_lost_particles(self, rel_max_lost_particles: float):
        cv.check_type('rel_max_lost_particles', rel_max_lost_particles, Real)
        cv.check_greater_than('rel_max_lost_particles', rel_max_lost_particles, 0)
        cv.check_less_than('rel_max_lost_particles', rel_max_lost_particles, 1)
        self._rel_max_lost_particles = rel_max_lost_particles

    @property
    def max_write_lost_particles(self) -> int:
        return self._max_write_lost_particles

    @max_write_lost_particles.setter
    def max_write_lost_particles(self, max_write_lost_particles: int):
        cv.check_type('max_write_lost_particles', max_write_lost_particles, Integral)
        cv.check_greater_than('max_write_lost_particles', max_write_lost_particles, 0)
        self._max_write_lost_particles = max_write_lost_particles

    @property
    def particles(self) -> int:
        return self._particles

    @particles.setter
    def particles(self, particles: int):
        cv.check_type('particles', particles, Integral)
        cv.check_greater_than('particles', particles, 0)
        self._particles = particles

    @property
    def keff_trigger(self) -> dict:
        return self._keff_trigger

    @keff_trigger.setter
    def keff_trigger(self, keff_trigger: dict):
        if not isinstance(keff_trigger, dict):
            msg = f'Unable to set a trigger on keff from "{keff_trigger}" ' \
                  'which is not a Python dictionary'
            raise ValueError(msg)

        elif 'type' not in keff_trigger:
            msg = f'Unable to set a trigger on keff from "{keff_trigger}" ' \
                  'which does not have a "type" key'
            raise ValueError(msg)

        elif keff_trigger['type'] not in ['variance', 'std_dev', 'rel_err']:
            msg = 'Unable to set a trigger on keff with ' \
                  'type "{0}"'.format(keff_trigger['type'])
            raise ValueError(msg)

        elif 'threshold' not in keff_trigger:
            msg = f'Unable to set a trigger on keff from "{keff_trigger}" ' \
                  'which does not have a "threshold" key'
            raise ValueError(msg)

        elif not isinstance(keff_trigger['threshold'], Real):
            msg = 'Unable to set a trigger on keff with ' \
                  'threshold "{0}"'.format(keff_trigger['threshold'])
            raise ValueError(msg)

        self._keff_trigger = keff_trigger

    @property
    def energy_mode(self) -> str:
        return self._energy_mode

    @energy_mode.setter
    def energy_mode(self, energy_mode: str):
        cv.check_value('energy mode', energy_mode,
                    ['continuous-energy', 'multi-group'])
        self._energy_mode = energy_mode

    @property
    def max_order(self) -> int:
        return self._max_order

    @max_order.setter
    def max_order(self, max_order: Optional[int]):
        if max_order is not None:
            cv.check_type('maximum scattering order', max_order, Integral)
            cv.check_greater_than('maximum scattering order', max_order, 0,
                                  True)
        self._max_order = max_order

    @property
    def source(self) -> typing.List[SourceBase]:
        return self._source

    @source.setter
    def source(self, source: typing.Union[SourceBase, typing.Iterable[SourceBase]]):
        if not isinstance(source, MutableSequence):
            source = [source]
        self._source = cv.CheckedList(SourceBase, 'source distributions', source)

    @property
    def confidence_intervals(self) -> bool:
        return self._confidence_intervals

    @confidence_intervals.setter
    def confidence_intervals(self, confidence_intervals: bool):
        cv.check_type('confidence interval', confidence_intervals, bool)
        self._confidence_intervals = confidence_intervals

    @property
    def electron_treatment(self) -> str:
        return self._electron_treatment

    @electron_treatment.setter
    def electron_treatment(self, electron_treatment: str):
        cv.check_value('electron treatment', electron_treatment, ['led', 'ttb'])
        self._electron_treatment = electron_treatment

    @property
    def ptables(self) -> bool:
        return self._ptables

    @ptables.setter
    def ptables(self, ptables: bool):
        cv.check_type('probability tables', ptables, bool)
        self._ptables = ptables

    @property
    def photon_transport(self) -> bool:
        return self._photon_transport

    @photon_transport.setter
    def photon_transport(self, photon_transport: bool):
        cv.check_type('photon transport', photon_transport, bool)
        self._photon_transport = photon_transport

    @property
    def plot_seed(self):
        return self._plot_seed

    @plot_seed.setter
    def plot_seed(self, seed):
        cv.check_type('random plot color seed', seed, Integral)
        cv.check_greater_than('random plot color seed', seed, 0)
        self._plot_seed = seed

    @property
    def seed(self) -> int:
        return self._seed

    @seed.setter
    def seed(self, seed: int):
        cv.check_type('random number generator seed', seed, Integral)
        cv.check_greater_than('random number generator seed', seed, 0)
        self._seed = seed

    @property
    def survival_biasing(self) -> bool:
        return self._survival_biasing

    @survival_biasing.setter
    def survival_biasing(self, survival_biasing: bool):
        cv.check_type('survival biasing', survival_biasing, bool)
        self._survival_biasing = survival_biasing

    @property
    def entropy_mesh(self) -> RegularMesh:
        return self._entropy_mesh

    @entropy_mesh.setter
    def entropy_mesh(self, entropy: RegularMesh):
        cv.check_type('entropy mesh', entropy, RegularMesh)
        self._entropy_mesh = entropy

    @property
    def trigger_active(self) -> bool:
        return self._trigger_active

    @trigger_active.setter
    def trigger_active(self, trigger_active: bool):
        cv.check_type('trigger active', trigger_active, bool)
        self._trigger_active = trigger_active

    @property
    def trigger_max_batches(self) -> int:
        return self._trigger_max_batches

    @trigger_max_batches.setter
    def trigger_max_batches(self, trigger_max_batches: int):
        cv.check_type('trigger maximum batches', trigger_max_batches, Integral)
        cv.check_greater_than('trigger maximum batches', trigger_max_batches, 0)
        self._trigger_max_batches = trigger_max_batches

    @property
    def trigger_batch_interval(self) -> int:
        return self._trigger_batch_interval

    @trigger_batch_interval.setter
    def trigger_batch_interval(self, trigger_batch_interval: int):
        cv.check_type('trigger batch interval', trigger_batch_interval, Integral)
        cv.check_greater_than('trigger batch interval', trigger_batch_interval, 0)
        self._trigger_batch_interval = trigger_batch_interval

    @property
    def output(self) -> dict:
        return self._output

    @output.setter
    def output(self, output: dict):
        cv.check_type('output', output, Mapping)
        for key, value in output.items():
            cv.check_value('output key', key, ('summary', 'tallies', 'path'))
            if key in ('summary', 'tallies'):
                cv.check_type(f"output['{key}']", value, bool)
            else:
                cv.check_type("output['path']", value, str)
        self._output = output

    @property
    def sourcepoint(self) -> dict:
        return self._sourcepoint

    @sourcepoint.setter
    def sourcepoint(self, sourcepoint: dict):
        cv.check_type('sourcepoint options', sourcepoint, Mapping)
        for key, value in sourcepoint.items():
            if key == 'batches':
                cv.check_type('sourcepoint batches', value, Iterable, Integral)
                for batch in value:
                    cv.check_greater_than('sourcepoint batch', batch, 0)
            elif key == 'separate':
                cv.check_type('sourcepoint separate', value, bool)
            elif key == 'write':
                cv.check_type('sourcepoint write', value, bool)
            elif key == 'overwrite':
                cv.check_type('sourcepoint overwrite', value, bool)
            elif key == 'mcpl':
                cv.check_type('sourcepoint mcpl', value, bool)
            else:
                raise ValueError(f"Unknown key '{key}' encountered when "
                                 "setting sourcepoint options.")
        self._sourcepoint = sourcepoint

    @property
    def statepoint(self) -> dict:
        return self._statepoint

    @statepoint.setter
    def statepoint(self, statepoint: dict):
        cv.check_type('statepoint options', statepoint, Mapping)
        for key, value in statepoint.items():
            if key == 'batches':
                cv.check_type('statepoint batches', value, Iterable, Integral)
                for batch in value:
                    cv.check_greater_than('statepoint batch', batch, 0)
            else:
                raise ValueError(f"Unknown key '{key}' encountered when "
                                 "setting statepoint options.")
        self._statepoint = statepoint

    @property
    def surf_source_read(self) -> dict:
        return self._surf_source_read

    @surf_source_read.setter
    def surf_source_read(self, surf_source_read: dict):
        cv.check_type('surface source reading options', surf_source_read, Mapping)
        for key, value in surf_source_read.items():
            cv.check_value('surface source reading key', key,
                           ('path'))
            if key == 'path':
                cv.check_type('path to surface source file', value, str)
        self._surf_source_read = surf_source_read

    @property
    def surf_source_write(self) -> dict:
        return self._surf_source_write

    @surf_source_write.setter
    def surf_source_write(self, surf_source_write: dict):
        cv.check_type('surface source writing options', surf_source_write, Mapping)
        for key, value in surf_source_write.items():
            cv.check_value('surface source writing key', key,
                           ('surface_ids', 'max_particles', 'mcpl'))
            if key == 'surface_ids':
                cv.check_type('surface ids for source banking', value,
                              Iterable, Integral)
                for surf_id in value:
                    cv.check_greater_than('surface id for source banking',
                                          surf_id, 0)
            elif key == 'max_particles':
                cv.check_type('maximum particle banks on surfaces per process',
                              value, Integral)
                cv.check_greater_than('maximum particle banks on surfaces per process',
                                      value, 0)
            elif key == 'mcpl':
                cv.check_type('write to an MCPL-format file', value, bool)

        self._surf_source_write = surf_source_write

    @property
    def no_reduce(self) -> bool:
        return self._no_reduce

    @no_reduce.setter
    def no_reduce(self, no_reduce: bool):
        cv.check_type('no reduction option', no_reduce, bool)
        self._no_reduce = no_reduce

    @property
    def verbosity(self) -> int:
        return self._verbosity

    @verbosity.setter
    def verbosity(self, verbosity: int):
        cv.check_type('verbosity', verbosity, Integral)
        cv.check_greater_than('verbosity', verbosity, 1, True)
        cv.check_less_than('verbosity', verbosity, 10, True)
        self._verbosity = verbosity

    @property
    def tabular_legendre(self) -> dict:
        return self._tabular_legendre

    @tabular_legendre.setter
    def tabular_legendre(self, tabular_legendre: dict):
        cv.check_type('tabular_legendre settings', tabular_legendre, Mapping)
        for key, value in tabular_legendre.items():
            cv.check_value('tabular_legendre key', key,
                           ['enable', 'num_points'])
            if key == 'enable':
                cv.check_type('enable tabular_legendre', value, bool)
            elif key == 'num_points':
                cv.check_type('num_points tabular_legendre', value, Integral)
                cv.check_greater_than('num_points tabular_legendre', value, 0)
        self._tabular_legendre = tabular_legendre

    @property
    def temperature(self) -> dict:
        return self._temperature

    @temperature.setter
    def temperature(self, temperature: dict):

        cv.check_type('temperature settings', temperature, Mapping)
        for key, value in temperature.items():
            cv.check_value('temperature key', key,
                           ['default', 'method', 'tolerance', 'multipole',
                            'range'])
            if key == 'default':
                cv.check_type('default temperature', value, Real)
            elif key == 'method':
                cv.check_value('temperature method', value,
                               ['nearest', 'interpolation'])
            elif key == 'tolerance':
                cv.check_type('temperature tolerance', value, Real)
            elif key == 'multipole':
                cv.check_type('temperature multipole', value, bool)
            elif key == 'range':
                cv.check_length('temperature range', value, 2)
                for T in value:
                    cv.check_type('temperature', T, Real)

        self._temperature = temperature

    @property
    def trace(self) -> typing.Iterable:
        return self._trace

    @trace.setter
    def trace(self, trace: Iterable):
        cv.check_type('trace', trace, Iterable, Integral)
        cv.check_length('trace', trace, 3)
        cv.check_greater_than('trace batch', trace[0], 0)
        cv.check_greater_than('trace generation', trace[1], 0)
        cv.check_greater_than('trace particle', trace[2], 0)
        self._trace = trace

    @property
    def track(self) -> typing.Iterable[typing.Iterable[int]]:
        return self._track

    @track.setter
    def track(self, track: typing.Iterable[typing.Iterable[int]]):
        cv.check_type('track', track, Iterable)
        for t in track:
            if len(t) != 3:
                msg = f'Unable to set the track to "{t}" since its length is not 3'
                raise ValueError(msg)
            cv.check_greater_than('track batch', t[0], 0)
            cv.check_greater_than('track generation', t[1], 0)
            cv.check_greater_than('track particle', t[2], 0)
            cv.check_type('track batch', t[0], Integral)
            cv.check_type('track generation', t[1], Integral)
            cv.check_type('track particle', t[2], Integral)
        self._track = track

    @property
    def cutoff(self) -> dict:
        return self._cutoff

    @cutoff.setter
    def cutoff(self, cutoff: dict):
        if not isinstance(cutoff, Mapping):
            msg = f'Unable to set cutoff from "{cutoff}" which is not a '\
                  'Python dictionary'
            raise ValueError(msg)
        for key in cutoff:
            if key == 'weight':
                cv.check_type('weight cutoff', cutoff[key], Real)
                cv.check_greater_than('weight cutoff', cutoff[key], 0.0)
            elif key == 'weight_avg':
                cv.check_type('average survival weight', cutoff[key], Real)
                cv.check_greater_than('average survival weight',
                                      cutoff[key], 0.0)
            elif key in ['energy_neutron', 'energy_photon', 'energy_electron',
                         'energy_positron']:
                cv.check_type('energy cutoff', cutoff[key], Real)
                cv.check_greater_than('energy cutoff', cutoff[key], 0.0)
            else:
                msg = f'Unable to set cutoff to "{key}" which is unsupported ' \
                      'by OpenMC'

        self._cutoff = cutoff

    @property
    def ufs_mesh(self) -> RegularMesh:
        return self._ufs_mesh

    @ufs_mesh.setter
    def ufs_mesh(self, ufs_mesh: RegularMesh):
        cv.check_type('UFS mesh', ufs_mesh, RegularMesh)
        cv.check_length('UFS mesh dimension', ufs_mesh.dimension, 3)
        cv.check_length('UFS mesh lower-left corner', ufs_mesh.lower_left, 3)
        cv.check_length('UFS mesh upper-right corner', ufs_mesh.upper_right, 3)
        self._ufs_mesh = ufs_mesh

    @property
    def resonance_scattering(self) -> dict:
        return self._resonance_scattering

    @resonance_scattering.setter
    def resonance_scattering(self, res: dict):
        cv.check_type('resonance scattering settings', res, Mapping)
        keys = ('enable', 'method', 'energy_min', 'energy_max', 'nuclides')
        for key, value in res.items():
            cv.check_value('resonance scattering dictionary key', key, keys)
            if key == 'enable':
                cv.check_type('resonance scattering enable', value, bool)
            elif key == 'method':
                cv.check_value('resonance scattering method', value,
                               _RES_SCAT_METHODS)
            elif key == 'energy_min':
                name = 'resonance scattering minimum energy'
                cv.check_type(name, value, Real)
                cv.check_greater_than(name, value, 0)
            elif key == 'energy_max':
                name = 'resonance scattering minimum energy'
                cv.check_type(name, value, Real)
                cv.check_greater_than(name, value, 0)
            elif key == 'nuclides':
                cv.check_type('resonance scattering nuclides', value,
                              Iterable, str)
        self._resonance_scattering = res

    @property
    def volume_calculations(self) -> typing.List[VolumeCalculation]:
        return self._volume_calculations

    @volume_calculations.setter
    def volume_calculations(
        self, vol_calcs: typing.Union[VolumeCalculation, typing.Iterable[VolumeCalculation]]
    ):
        if not isinstance(vol_calcs, MutableSequence):
            vol_calcs = [vol_calcs]
        self._volume_calculations = cv.CheckedList(
            VolumeCalculation, 'stochastic volume calculations', vol_calcs)

    @property
    def create_fission_neutrons(self) -> bool:
        return self._create_fission_neutrons

    @create_fission_neutrons.setter
    def create_fission_neutrons(self, create_fission_neutrons: bool):
        cv.check_type('Whether create fission neutrons',
                      create_fission_neutrons, bool)
        self._create_fission_neutrons = create_fission_neutrons

    @property
    def create_delayed_neutrons(self) -> bool:
        return self._create_delayed_neutrons

    @create_delayed_neutrons.setter
    def create_delayed_neutrons(self, create_delayed_neutrons: bool):
        cv.check_type('Whether create only prompt neutrons',
                      create_delayed_neutrons, bool)
        self._create_delayed_neutrons = create_delayed_neutrons

    @property
    def delayed_photon_scaling(self) -> bool:
        return self._delayed_photon_scaling

    @delayed_photon_scaling.setter
    def delayed_photon_scaling(self, value: bool):
        cv.check_type('delayed photon scaling', value, bool)
        self._delayed_photon_scaling = value

    @property
    def material_cell_offsets(self) -> bool:
        return self._material_cell_offsets

    @material_cell_offsets.setter
    def material_cell_offsets(self, value: bool):
        cv.check_type('material cell offsets', value, bool)
        self._material_cell_offsets = value

    @property
    def log_grid_bins(self) -> int:
        return self._log_grid_bins

    @log_grid_bins.setter
    def log_grid_bins(self, log_grid_bins: int):
        cv.check_type('log grid bins', log_grid_bins, Real)
        cv.check_greater_than('log grid bins', log_grid_bins, 0)
        self._log_grid_bins = log_grid_bins

    @property
    def event_based(self) -> bool:
        return self._event_based

    @event_based.setter
    def event_based(self, value: bool):
        cv.check_type('event based', value, bool)
        self._event_based = value

    @property
    def max_particles_in_flight(self) -> int:
        return self._max_particles_in_flight

    @max_particles_in_flight.setter
    def max_particles_in_flight(self, value: int):
        cv.check_type('max particles in flight', value, Integral)
        cv.check_greater_than('max particles in flight', value, 0)
        self._max_particles_in_flight = value

    @property
    def max_particle_events(self) -> int:
        return self._max_particle_events

    @max_particle_events.setter
    def max_particle_events(self, value: int):
        cv.check_type('max particle events', value, Integral)
        cv.check_greater_than('max particle events', value, 0)
        self._max_particle_events = value

    @property
    def write_initial_source(self) -> bool:
        return self._write_initial_source

    @write_initial_source.setter
    def write_initial_source(self, value: bool):
        cv.check_type('write initial source', value, bool)
        self._write_initial_source = value

    @property
    def weight_windows(self) -> typing.List[WeightWindows]:
        return self._weight_windows

    @weight_windows.setter
    def weight_windows(self, value: typing.Union[WeightWindows, typing.Iterable[WeightWindows]]):
        if not isinstance(value, MutableSequence):
            value = [value]
        self._weight_windows = cv.CheckedList(WeightWindows, 'weight windows', value)

    @property
    def weight_windows_on(self) -> bool:
        return self._weight_windows_on

    @weight_windows_on.setter
    def weight_windows_on(self, value: bool):
        cv.check_type('weight windows on', value, bool)
        self._weight_windows_on = value

    @property
    def weight_window_checkpoints(self) -> dict:
        return self._weight_window_checkpoints

    @weight_window_checkpoints.setter
    def weight_window_checkpoints(self, weight_window_checkpoints: dict):
        for key in weight_window_checkpoints.keys():
            cv.check_value('weight_window_checkpoints', key, ('collision', 'surface'))
        self._weight_window_checkpoints = weight_window_checkpoints

    @property
    def max_splits(self):
        raise AttributeError('max_splits has been deprecated. Please use max_history_splits instead')

    @property
    def max_history_splits(self) -> int:
        return self._max_history_splits

    @max_history_splits.setter
    def max_history_splits(self, value: int):
        cv.check_type('maximum particle splits', value, Integral)
        cv.check_greater_than('max particle splits', value, 0)
        self._max_history_splits = value

    @property
    def max_tracks(self) -> int:
        return self._max_tracks

    @max_tracks.setter
    def max_tracks(self, value: int):
        cv.check_type('maximum particle tracks', value, Integral)
        cv.check_greater_than('maximum particle tracks', value, 0, True)
        self._max_tracks = value

    @property
    def weight_windows_file(self) -> Optional[PathLike]:
        return self._weight_windows_file

    @weight_windows_file.setter
    def weight_windows_file(self, value: PathLike):
        cv.check_type('weight windows file', value, (str, Path))
        self._weight_windows_file = value

    @property
    def weight_window_generators(self) -> typing.List[WeightWindowGenerator]:
        return self._weight_window_generators

    @weight_window_generators.setter
    def weight_window_generators(self, wwgs):
        if not isinstance(wwgs, MutableSequence):
            wwgs = [wwgs]
        self._weight_window_generators = cv.CheckedList(WeightWindowGenerator, 'weight window generators', wwgs)

    @property
    def random_ray(self) -> dict:
        return self._random_ray

    @random_ray.setter
    def random_ray(self, random_ray: dict):
        if not isinstance(random_ray, Mapping):
            raise ValueError(f'Unable to set random_ray from "{random_ray}" '
                             'which is not a dict.')
        for key in random_ray:
            if key == 'distance_active':
                cv.check_type('active ray length', random_ray[key], Real)
                cv.check_greater_than('active ray length', random_ray[key], 0.0)
            elif key == 'distance_inactive':
                cv.check_type('inactive ray length', random_ray[key], Real)
                cv.check_greater_than('inactive ray length',
                                      random_ray[key], 0.0, True)
            elif key == 'ray_source':
                cv.check_type('random ray source', random_ray[key], SourceBase)
            else:
                raise ValueError(f'Unable to set random ray to "{key}" which is '
                                 'unsupported by OpenMC')

        self._random_ray = random_ray

    def _create_run_mode_subelement(self, root):
        elem = ET.SubElement(root, "run_mode")
        elem.text = self._run_mode.value

    def _create_batches_subelement(self, root):
        if self._batches is not None:
            element = ET.SubElement(root, "batches")
            element.text = str(self._batches)

    def _create_generations_per_batch_subelement(self, root):
        if self._generations_per_batch is not None:
            element = ET.SubElement(root, "generations_per_batch")
            element.text = str(self._generations_per_batch)

    def _create_inactive_subelement(self, root):
        if self._inactive is not None:
            element = ET.SubElement(root, "inactive")
            element.text = str(self._inactive)

    def _create_max_lost_particles_subelement(self, root):
        if self._max_lost_particles is not None:
            element = ET.SubElement(root, "max_lost_particles")
            element.text = str(self._max_lost_particles)

    def _create_rel_max_lost_particles_subelement(self, root):
        if self._rel_max_lost_particles is not None:
            element = ET.SubElement(root, "rel_max_lost_particles")
            element.text = str(self._rel_max_lost_particles)

    def _create_max_write_lost_particles_subelement(self, root):
        if self._max_write_lost_particles is not None:
            element = ET.SubElement(root, "max_write_lost_particles")
            element.text = str(self._max_write_lost_particles)

    def _create_particles_subelement(self, root):
        if self._particles is not None:
            element = ET.SubElement(root, "particles")
            element.text = str(self._particles)

    def _create_keff_trigger_subelement(self, root):
        if self._keff_trigger is not None:
            element = ET.SubElement(root, "keff_trigger")
            for key, value in sorted(self._keff_trigger.items()):
                subelement = ET.SubElement(element, key)
                subelement.text = str(value).lower()

    def _create_energy_mode_subelement(self, root):
        if self._energy_mode is not None:
            element = ET.SubElement(root, "energy_mode")
            element.text = str(self._energy_mode)

    def _create_max_order_subelement(self, root):
        if self._max_order is not None:
            element = ET.SubElement(root, "max_order")
            element.text = str(self._max_order)

    def _create_source_subelement(self, root, mesh_memo=None):
        for source in self.source:
            root.append(source.to_xml_element())
            if isinstance(source, IndependentSource) and isinstance(source.space, MeshSpatial):
                path = f"./mesh[@id='{source.space.mesh.id}']"
                if root.find(path) is None:
                    root.append(source.space.mesh.to_xml_element())
            if isinstance(source, MeshSource):
                path = f"./mesh[@id='{source.mesh.id}']"
                if root.find(path) is None:
                    root.append(source.mesh.to_xml_element())
                    if mesh_memo is not None:
                        mesh_memo.add(source.mesh.id)

    def _create_volume_calcs_subelement(self, root):
        for calc in self.volume_calculations:
            root.append(calc.to_xml_element())

    def _create_output_subelement(self, root):
        if self._output is not None:
            element = ET.SubElement(root, "output")
            for key, value in sorted(self._output.items()):
                subelement = ET.SubElement(element, key)
                if key in ('summary', 'tallies'):
                    subelement.text = str(value).lower()
                else:
                    subelement.text = value

    def _create_verbosity_subelement(self, root):
        if self._verbosity is not None:
            element = ET.SubElement(root, "verbosity")
            element.text = str(self._verbosity)

    def _create_statepoint_subelement(self, root):
        if self._statepoint:
            element = ET.SubElement(root, "state_point")
            if 'batches' in self._statepoint:
                subelement = ET.SubElement(element, "batches")
                subelement.text = ' '.join(
                    str(x) for x in self._statepoint['batches'])

    def _create_sourcepoint_subelement(self, root):
        if self._sourcepoint:
            element = ET.SubElement(root, "source_point")

            if 'batches' in self._sourcepoint:
                subelement = ET.SubElement(element, "batches")
                subelement.text = ' '.join(
                    str(x) for x in self._sourcepoint['batches'])

            if 'separate' in self._sourcepoint:
                subelement = ET.SubElement(element, "separate")
                subelement.text = str(self._sourcepoint['separate']).lower()

            if 'write' in self._sourcepoint:
                subelement = ET.SubElement(element, "write")
                subelement.text = str(self._sourcepoint['write']).lower()

            # Overwrite latest subelement
            if 'overwrite' in self._sourcepoint:
                subelement = ET.SubElement(element, "overwrite_latest")
                subelement.text = str(self._sourcepoint['overwrite']).lower()

            if 'mcpl' in self._sourcepoint:
                subelement = ET.SubElement(element, "mcpl")
                subelement.text = str(self._sourcepoint['mcpl']).lower()

    def _create_surf_source_read_subelement(self, root):
        if self._surf_source_read:
            element = ET.SubElement(root, "surf_source_read")
            if 'path' in self._surf_source_read:
                subelement = ET.SubElement(element, "path")
                subelement.text = self._surf_source_read['path']

    def _create_surf_source_write_subelement(self, root):
        if self._surf_source_write:
            element = ET.SubElement(root, "surf_source_write")
            if 'surface_ids' in self._surf_source_write:
                subelement = ET.SubElement(element, "surface_ids")
                subelement.text = ' '.join(
                    str(x) for x in self._surf_source_write['surface_ids'])
            if 'max_particles' in self._surf_source_write:
                subelement = ET.SubElement(element, "max_particles")
                subelement.text = str(self._surf_source_write['max_particles'])
            if 'mcpl' in self._surf_source_write:
                subelement = ET.SubElement(element, "mcpl")
                subelement.text = str(self._surf_source_write['mcpl']).lower()

    def _create_confidence_intervals(self, root):
        if self._confidence_intervals is not None:
            element = ET.SubElement(root, "confidence_intervals")
            element.text = str(self._confidence_intervals).lower()

    def _create_electron_treatment_subelement(self, root):
        if self._electron_treatment is not None:
            element = ET.SubElement(root, "electron_treatment")
            element.text = str(self._electron_treatment)

    def _create_photon_transport_subelement(self, root):
        if self._photon_transport is not None:
            element = ET.SubElement(root, "photon_transport")
            element.text = str(self._photon_transport).lower()

    def _create_plot_seed_subelement(self, root):
        if self._plot_seed is not None:
            element = ET.SubElement(root, "plot_seed")
            element.text = str(self._plot_seed)

    def _create_ptables_subelement(self, root):
        if self._ptables is not None:
            element = ET.SubElement(root, "ptables")
            element.text = str(self._ptables).lower()

    def _create_seed_subelement(self, root):
        if self._seed is not None:
            element = ET.SubElement(root, "seed")
            element.text = str(self._seed)

    def _create_survival_biasing_subelement(self, root):
        if self._survival_biasing is not None:
            element = ET.SubElement(root, "survival_biasing")
            element.text = str(self._survival_biasing).lower()

    def _create_cutoff_subelement(self, root):
        if self._cutoff is not None:
            element = ET.SubElement(root, "cutoff")
            for key, value in self._cutoff.items():
                subelement = ET.SubElement(element, key)
                subelement.text = str(value)

    def _create_entropy_mesh_subelement(self, root, mesh_memo=None):
        if self.entropy_mesh is None:
            return

        # use default heuristic for entropy mesh if not set by user
        if self.entropy_mesh.dimension is None:
            if self.particles is None:
                raise RuntimeError("Number of particles must be set in order to " \
                    "use entropy mesh dimension heuristic")
            else:
                n = ceil((self.particles / 20.0)**(1.0 / 3.0))
                d = len(self.entropy_mesh.lower_left)
                self.entropy_mesh.dimension = (n,)*d

        # add mesh ID to this element
        subelement = ET.SubElement(root, "entropy_mesh")
        subelement.text = str(self.entropy_mesh.id)

        # If this mesh has already been written outside the
        # settings element, skip writing it again
        if mesh_memo and self.entropy_mesh.id in mesh_memo:
            return

        # See if a <mesh> element already exists -- if not, add it
        path = f"./mesh[@id='{self.entropy_mesh.id}']"
        if root.find(path) is None:
            root.append(self.entropy_mesh.to_xml_element())
            if mesh_memo is not None:
                mesh_memo.add(self.entropy_mesh.id)

    def _create_trigger_subelement(self, root):
        if self._trigger_active is not None:
            trigger_element = ET.SubElement(root, "trigger")
            element = ET.SubElement(trigger_element, "active")
            element.text = str(self._trigger_active).lower()

            if self._trigger_max_batches is not None:
                element = ET.SubElement(trigger_element, "max_batches")
                element.text = str(self._trigger_max_batches)

            if self._trigger_batch_interval is not None:
                element = ET.SubElement(trigger_element, "batch_interval")
                element.text = str(self._trigger_batch_interval)

    def _create_no_reduce_subelement(self, root):
        if self._no_reduce is not None:
            element = ET.SubElement(root, "no_reduce")
            element.text = str(self._no_reduce).lower()

    def _create_tabular_legendre_subelements(self, root):
        if self.tabular_legendre:
            element = ET.SubElement(root, "tabular_legendre")
            subelement = ET.SubElement(element, "enable")
            subelement.text = str(self._tabular_legendre['enable']).lower()
            if 'num_points' in self._tabular_legendre:
                subelement = ET.SubElement(element, "num_points")
                subelement.text = str(self._tabular_legendre['num_points'])

    def _create_temperature_subelements(self, root):
        if self.temperature:
            for key, value in sorted(self.temperature.items()):
                element = ET.SubElement(root, f"temperature_{key}")
                if isinstance(value, bool):
                    element.text = str(value).lower()
                elif key == 'range':
                    element.text = ' '.join(str(T) for T in value)
                else:
                    element.text = str(value)

    def _create_trace_subelement(self, root):
        if self._trace is not None:
            element = ET.SubElement(root, "trace")
            element.text = ' '.join(map(str, self._trace))

    def _create_track_subelement(self, root):
        if self._track is not None:
            element = ET.SubElement(root, "track")
            element.text = ' '.join(map(str, itertools.chain(*self._track)))

    def _create_ufs_mesh_subelement(self, root, mesh_memo=None):
        if self.ufs_mesh is None:
            return

        subelement = ET.SubElement(root, "ufs_mesh")
        subelement.text = str(self.ufs_mesh.id)

        if mesh_memo and self.ufs_mesh.id in mesh_memo:
            return

        # See if a <mesh> element already exists -- if not, add it
        path = f"./mesh[@id='{self.ufs_mesh.id}']"
        if root.find(path) is None:
            root.append(self.ufs_mesh.to_xml_element())
            if mesh_memo is not None: mesh_memo.add(self.ufs_mesh.id)

    def _create_resonance_scattering_subelement(self, root):
        res = self.resonance_scattering
        if res:
            elem = ET.SubElement(root, 'resonance_scattering')
            if 'enable' in res:
                subelem = ET.SubElement(elem, 'enable')
                subelem.text = str(res['enable']).lower()
            if 'method' in res:
                subelem = ET.SubElement(elem, 'method')
                subelem.text = res['method']
            if 'energy_min' in res:
                subelem = ET.SubElement(elem, 'energy_min')
                subelem.text = str(res['energy_min'])
            if 'energy_max' in res:
                subelem = ET.SubElement(elem, 'energy_max')
                subelem.text = str(res['energy_max'])
            if 'nuclides' in res:
                subelem = ET.SubElement(elem, 'nuclides')
                subelem.text = ' '.join(res['nuclides'])

    def _create_create_fission_neutrons_subelement(self, root):
        if self._create_fission_neutrons is not None:
            elem = ET.SubElement(root, "create_fission_neutrons")
            elem.text = str(self._create_fission_neutrons).lower()

    def _create_create_delayed_neutrons_subelement(self, root):
       if self._create_delayed_neutrons is not None:
           elem = ET.SubElement(root, "create_delayed_neutrons")
           elem.text = str(self._create_delayed_neutrons).lower()

    def _create_delayed_photon_scaling_subelement(self, root):
        if self._delayed_photon_scaling is not None:
            elem = ET.SubElement(root, "delayed_photon_scaling")
            elem.text = str(self._delayed_photon_scaling).lower()

    def _create_event_based_subelement(self, root):
        if self._event_based is not None:
            elem = ET.SubElement(root, "event_based")
            elem.text = str(self._event_based).lower()

    def _create_max_particles_in_flight_subelement(self, root):
        if self._max_particles_in_flight is not None:
            elem = ET.SubElement(root, "max_particles_in_flight")
            elem.text = str(self._max_particles_in_flight).lower()

    def _create_max_events_subelement(self, root):
        if self._max_particle_events is not None:
            elem = ET.SubElement(root, "max_particle_events")
            elem.text = str(self._max_particle_events).lower()

    def _create_material_cell_offsets_subelement(self, root):
        if self._material_cell_offsets is not None:
            elem = ET.SubElement(root, "material_cell_offsets")
            elem.text = str(self._material_cell_offsets).lower()

    def _create_log_grid_bins_subelement(self, root):
        if self._log_grid_bins is not None:
            elem = ET.SubElement(root, "log_grid_bins")
            elem.text = str(self._log_grid_bins)

    def _create_write_initial_source_subelement(self, root):
        if self._write_initial_source is not None:
            elem = ET.SubElement(root, "write_initial_source")
            elem.text = str(self._write_initial_source).lower()

    def _create_weight_windows_subelement(self, root, mesh_memo=None):
        for ww in self._weight_windows:
            # Add weight window information
            root.append(ww.to_xml_element())

            # if this mesh has already been written,
            # skip writing the mesh element
            if mesh_memo and ww.mesh.id in mesh_memo:
                continue

            # See if a <mesh> element already exists -- if not, add it
            path = f"./mesh[@id='{ww.mesh.id}']"
            if root.find(path) is None:
                root.append(ww.mesh.to_xml_element())
                if mesh_memo is not None:
                    mesh_memo.add(ww.mesh.id)

        if self._weight_windows_on is not None:
            elem = ET.SubElement(root, "weight_windows_on")
            elem.text = str(self._weight_windows_on).lower()

    def _create_weight_window_generators_subelement(self, root, mesh_memo=None):
        if not self.weight_window_generators:
            return
        elem = ET.SubElement(root, 'weight_window_generators')
        for wwg in self.weight_window_generators:
            elem.append(wwg.to_xml_element())

        # ensure that mesh elements are created if needed
        for wwg in self.weight_window_generators:
            if mesh_memo is not None and wwg.mesh.id in mesh_memo:
                continue

            root.append(wwg.mesh.to_xml_element())
            if mesh_memo is not None:
                mesh_memo.add(wwg.mesh)

    def _create_weight_windows_file_element(self, root):
        if self.weight_windows_file is not None:
            element = ET.Element("weight_windows_file")
            element.text = self.weight_windows_file
            root.append(element)

    def _create_weight_window_checkpoints_subelement(self, root):
        if not self._weight_window_checkpoints:
            return
        element = ET.SubElement(root, "weight_window_checkpoints")

        if 'collision' in self._weight_window_checkpoints:
            subelement = ET.SubElement(element, "collision")
            subelement.text = str(self._weight_window_checkpoints['collision']).lower()

        if 'surface' in self._weight_window_checkpoints:
            subelement = ET.SubElement(element, "surface")
            subelement.text = str(self._weight_window_checkpoints['surface']).lower()

    def _create_max_history_splits_subelement(self, root):
        if self._max_history_splits is not None:
            elem = ET.SubElement(root, "max_history_splits")
            elem.text = str(self._max_history_splits)

    def _create_max_tracks_subelement(self, root):
        if self._max_tracks is not None:
            elem = ET.SubElement(root, "max_tracks")
            elem.text = str(self._max_tracks)

    def _create_random_ray_subelement(self, root):
        if self._random_ray:
            element = ET.SubElement(root, "random_ray")
            for key, value in self._random_ray.items():
                if key == 'ray_source' and isinstance(value, SourceBase):
                    source_element = value.to_xml_element()
                    element.append(source_element)
                else:
                    subelement = ET.SubElement(element, key)
                    subelement.text = str(value)

    def _eigenvalue_from_xml_element(self, root):
        elem = root.find('eigenvalue')
        if elem is not None:
            self._run_mode_from_xml_element(elem)
            self._particles_from_xml_element(elem)
            self._batches_from_xml_element(elem)
            self._inactive_from_xml_element(elem)
            self._max_lost_particles_from_xml_element(elem)
            self._rel_max_lost_particles_from_xml_element(elem)
            self._max_write_lost_particles_from_xml_element(elem)
            self._generations_per_batch_from_xml_element(elem)

    def _run_mode_from_xml_element(self, root):
        text = get_text(root, 'run_mode')
        if text is not None:
            self.run_mode = text

    def _particles_from_xml_element(self, root):
        text = get_text(root, 'particles')
        if text is not None:
            self.particles = int(text)

    def _batches_from_xml_element(self, root):
        text = get_text(root, 'batches')
        if text is not None:
            self.batches = int(text)

    def _inactive_from_xml_element(self, root):
        text = get_text(root, 'inactive')
        if text is not None:
            self.inactive = int(text)

    def _max_lost_particles_from_xml_element(self, root):
        text = get_text(root, 'max_lost_particles')
        if text is not None:
            self.max_lost_particles = int(text)

    def _rel_max_lost_particles_from_xml_element(self, root):
        text = get_text(root, 'rel_max_lost_particles')
        if text is not None:
            self.rel_max_lost_particles = float(text)

    def _max_write_lost_particles_from_xml_element(self, root):
        text = get_text(root, 'max_write_lost_particles')
        if text is not None:
            self.max_write_lost_particles = int(text)

    def _generations_per_batch_from_xml_element(self, root):
        text = get_text(root, 'generations_per_batch')
        if text is not None:
            self.generations_per_batch = int(text)

    def _keff_trigger_from_xml_element(self, root):
        elem = root.find('keff_trigger')
        if elem is not None:
            trigger = get_text(elem, 'type')
            threshold = float(get_text(elem, 'threshold'))
            self.keff_trigger = {'type': trigger, 'threshold': threshold}

    def _source_from_xml_element(self, root, meshes=None):
        for elem in root.findall('source'):
            src = SourceBase.from_xml_element(elem, meshes)
            # add newly constructed source object to the list
            self.source.append(src)

    def _volume_calcs_from_xml_element(self, root):
        volume_elems = root.findall("volume_calc")
        if volume_elems:
            self.volume_calculations = [VolumeCalculation.from_xml_element(elem)
                                        for elem in volume_elems]

    def _output_from_xml_element(self, root):
        elem = root.find('output')
        if elem is not None:
            self.output = {}
            for key in ('summary', 'tallies', 'path'):
                value = get_text(elem, key)
                if value is not None:
                    if key in ('summary', 'tallies'):
                        value = value in ('true', '1')
                    self.output[key] = value

    def _statepoint_from_xml_element(self, root):
        elem = root.find('state_point')
        if elem is not None:
            text = get_text(elem, 'batches')
            if text is not None:
                self.statepoint['batches'] = [int(x) for x in text.split()]

    def _sourcepoint_from_xml_element(self, root):
        elem = root.find('source_point')
        if elem is not None:
            for key in ('separate', 'write', 'overwrite_latest', 'batches', 'mcpl'):
                value = get_text(elem, key)
                if value is not None:
                    if key in ('separate', 'write', 'mcpl'):
                        value = value in ('true', '1')
                    elif key == 'overwrite_latest':
                        value = value in ('true', '1')
                        key = 'overwrite'
                    else:
                        value = [int(x) for x in value.split()]
                    self.sourcepoint[key] = value

    def _surf_source_read_from_xml_element(self, root):
        elem = root.find('surf_source_read')
        if elem is not None:
            value = get_text(elem, 'path')
            if value is not None:
                self.surf_source_read['path'] = value

    def _surf_source_write_from_xml_element(self, root):
        elem = root.find('surf_source_write')
        if elem is not None:
            for key in ('surface_ids', 'max_particles','mcpl'):
                value = get_text(elem, key)
                if value is not None:
                    if key == 'surface_ids':
                        value = [int(x) for x in value.split()]
                    elif key in ('max_particles'):
                        value = int(value)
                    elif key == 'mcpl':
                        value = value in ('true', '1')
                    self.surf_source_write[key] = value

    def _confidence_intervals_from_xml_element(self, root):
        text = get_text(root, 'confidence_intervals')
        if text is not None:
            self.confidence_intervals = text in ('true', '1')

    def _electron_treatment_from_xml_element(self, root):
        text = get_text(root, 'electron_treatment')
        if text is not None:
            self.electron_treatment = text

    def _energy_mode_from_xml_element(self, root):
        text = get_text(root, 'energy_mode')
        if text is not None:
            self.energy_mode = text

    def _max_order_from_xml_element(self, root):
        text = get_text(root, 'max_order')
        if text is not None:
            self.max_order = int(text)

    def _photon_transport_from_xml_element(self, root):
        text = get_text(root, 'photon_transport')
        if text is not None:
            self.photon_transport = text in ('true', '1')

    def _plot_seed_from_xml_element(self, root):
        text = get_text(root, 'plot_seed')
        if text is not None:
            self.plot_seed = int(text)

    def _ptables_from_xml_element(self, root):
        text = get_text(root, 'ptables')
        if text is not None:
            self.ptables = text in ('true', '1')

    def _seed_from_xml_element(self, root):
        text = get_text(root, 'seed')
        if text is not None:
            self.seed = int(text)

    def _survival_biasing_from_xml_element(self, root):
        text = get_text(root, 'survival_biasing')
        if text is not None:
            self.survival_biasing = text in ('true', '1')

    def _cutoff_from_xml_element(self, root):
        elem = root.find('cutoff')
        if elem is not None:
            self.cutoff = {}
            for key in ('energy_neutron', 'energy_photon', 'energy_electron',
                        'energy_positron', 'weight', 'weight_avg', 'time_neutron',
                        'time_photon', 'time_electron', 'time_positron'):
                value = get_text(elem, key)
                if value is not None:
                    self.cutoff[key] = float(value)

    def _entropy_mesh_from_xml_element(self, root, meshes):
        text = get_text(root, 'entropy_mesh')
        if text is None:
            return
        mesh_id = int(text)
        if mesh_id not in meshes:
            raise ValueError(f'Could not locate mesh with ID "{mesh_id}"')
        self.entropy_mesh = meshes[mesh_id]

    def _trigger_from_xml_element(self, root):
        elem = root.find('trigger')
        if elem is not None:
            self.trigger_active = get_text(elem, 'active') in ('true', '1')
            text = get_text(elem, 'max_batches')
            if text is not None:
                self.trigger_max_batches = int(text)
            text = get_text(elem, 'batch_interval')
            if text is not None:
                self.trigger_batch_interval = int(text)

    def _no_reduce_from_xml_element(self, root):
        text = get_text(root, 'no_reduce')
        if text is not None:
            self.no_reduce = text in ('true', '1')

    def _verbosity_from_xml_element(self, root):
        text = get_text(root, 'verbosity')
        if text is not None:
            self.verbosity = int(text)

    def _tabular_legendre_from_xml_element(self, root):
        elem = root.find('tabular_legendre')
        if elem is not None:
            text = get_text(elem, 'enable')
            self.tabular_legendre['enable'] = text in ('true', '1')
            text = get_text(elem, 'num_points')
            if text is not None:
                self.tabular_legendre['num_points'] = int(text)

    def _temperature_from_xml_element(self, root):
        text = get_text(root, 'temperature_default')
        if text is not None:
            self.temperature['default'] = float(text)
        text = get_text(root, 'temperature_tolerance')
        if text is not None:
            self.temperature['tolerance'] = float(text)
        text = get_text(root, 'temperature_method')
        if text is not None:
            self.temperature['method'] = text
        text = get_text(root, 'temperature_range')
        if text is not None:
            self.temperature['range'] = [float(x) for x in text.split()]
        text = get_text(root, 'temperature_multipole')
        if text is not None:
            self.temperature['multipole'] = text in ('true', '1')

    def _trace_from_xml_element(self, root):
        text = get_text(root, 'trace')
        if text is not None:
            self.trace = [int(x) for x in text.split()]

    def _track_from_xml_element(self, root):
        text = get_text(root, 'track')
        if text is not None:
            values = [int(x) for x in text.split()]
            self.track = list(zip(values[::3], values[1::3], values[2::3]))

    def _ufs_mesh_from_xml_element(self, root, meshes):
        text = get_text(root, 'ufs_mesh')
        if text is None:
            return
        mesh_id = int(text)
        if mesh_id not in meshes:
            raise ValueError(f'Could not locate mesh with ID "{mesh_id}"')
        self.ufs_mesh = meshes[mesh_id]

    def _resonance_scattering_from_xml_element(self, root):
        elem = root.find('resonance_scattering')
        if elem is not None:
            keys = ('enable', 'method', 'energy_min', 'energy_max', 'nuclides')
            for key in keys:
                value = get_text(elem, key)
                if value is not None:
                    if key == 'enable':
                        value = value in ('true', '1')
                    elif key in ('energy_min', 'energy_max'):
                        value = float(value)
                    elif key == 'nuclides':
                        value = value.split()
                    self.resonance_scattering[key] = value

    def _create_fission_neutrons_from_xml_element(self, root):
        text = get_text(root, 'create_fission_neutrons')
        if text is not None:
            self.create_fission_neutrons = text in ('true', '1')

    def _create_delayed_neutrons_from_xml_element(self, root):
        text = get_text(root, 'create_delayed_neutrons')
        if text is not None:
            self.create_delayed_neutrons = text in ('true', '1')

    def _delayed_photon_scaling_from_xml_element(self, root):
        text = get_text(root, 'delayed_photon_scaling')
        if text is not None:
            self.delayed_photon_scaling = text in ('true', '1')

    def _event_based_from_xml_element(self, root):
        text = get_text(root, 'event_based')
        if text is not None:
            self.event_based = text in ('true', '1')

    def _max_particles_in_flight_from_xml_element(self, root):
        text = get_text(root, 'max_particles_in_flight')
        if text is not None:
            self.max_particles_in_flight = int(text)

    def _max_particle_events_from_xml_element(self, root):
        text = get_text(root, 'max_particle_events')
        if text is not None:
            self.max_particle_events = int(text)

    def _material_cell_offsets_from_xml_element(self, root):
        text = get_text(root, 'material_cell_offsets')
        if text is not None:
            self.material_cell_offsets = text in ('true', '1')

    def _log_grid_bins_from_xml_element(self, root):
        text = get_text(root, 'log_grid_bins')
        if text is not None:
            self.log_grid_bins = int(text)

    def _write_initial_source_from_xml_element(self, root):
        text = get_text(root, 'write_initial_source')
        if text is not None:
            self.write_initial_source = text in ('true', '1')

    def _weight_window_generators_from_xml_element(self, root, meshes=None):
        for elem in root.iter('weight_windows_generator'):
            wwg = WeightWindowGenerator.from_xml_element(elem, meshes)
            self.weight_window_generators.append(wwg)

    def _weight_windows_from_xml_element(self, root, meshes=None):
        for elem in root.findall('weight_windows'):
            ww = WeightWindows.from_xml_element(elem, meshes)
            self.weight_windows.append(ww)

        text = get_text(root, 'weight_windows_on')
        if text is not None:
            self.weight_windows_on = text in ('true', '1')

    def _weight_window_checkpoints_from_xml_element(self, root):
        elem = root.find('weight_window_checkpoints')
        if elem is None:
            return
        for key in ('collision', 'surface'):
            value = get_text(elem, key)
            if value is not None:
                value = value in ('true', '1')
                self.weight_window_checkpoints[key] = value

    def _max_history_splits_from_xml_element(self, root):
        text = get_text(root, 'max_history_splits')
        if text is not None:
            self.max_history_splits = int(text)

    def _max_tracks_from_xml_element(self, root):
        text = get_text(root, 'max_tracks')
        if text is not None:
            self.max_tracks = int(text)

    def _random_ray_from_xml_element(self, root):
        elem = root.find('random_ray')
        if elem is not None:
            self.random_ray = {}
            for child in elem:
                if child.tag in ('distance_inactive', 'distance_active'):
                    self.random_ray[child.tag] = float(child.text)
                elif child.tag == 'source':
                    source = SourceBase.from_xml_element(child)
                    self.random_ray['ray_source'] = source

    def to_xml_element(self, mesh_memo=None):
        """Create a 'settings' element to be written to an XML file.

        Parameters
        ----------
        mesh_memo : set of ints
            A set of mesh IDs to keep track of whether a mesh has already been written.
        """
        # Reset xml element tree
        element = ET.Element("settings")

        self._create_run_mode_subelement(element)
        self._create_particles_subelement(element)
        self._create_batches_subelement(element)
        self._create_inactive_subelement(element)
        self._create_max_lost_particles_subelement(element)
        self._create_rel_max_lost_particles_subelement(element)
        self._create_max_write_lost_particles_subelement(element)
        self._create_generations_per_batch_subelement(element)
        self._create_keff_trigger_subelement(element)
        self._create_source_subelement(element, mesh_memo)
        self._create_output_subelement(element)
        self._create_statepoint_subelement(element)
        self._create_sourcepoint_subelement(element)
        self._create_surf_source_read_subelement(element)
        self._create_surf_source_write_subelement(element)
        self._create_confidence_intervals(element)
        self._create_electron_treatment_subelement(element)
        self._create_energy_mode_subelement(element)
        self._create_max_order_subelement(element)
        self._create_photon_transport_subelement(element)
        self._create_plot_seed_subelement(element)
        self._create_ptables_subelement(element)
        self._create_seed_subelement(element)
        self._create_survival_biasing_subelement(element)
        self._create_cutoff_subelement(element)
        self._create_entropy_mesh_subelement(element, mesh_memo)
        self._create_trigger_subelement(element)
        self._create_no_reduce_subelement(element)
        self._create_verbosity_subelement(element)
        self._create_tabular_legendre_subelements(element)
        self._create_temperature_subelements(element)
        self._create_trace_subelement(element)
        self._create_track_subelement(element)
        self._create_ufs_mesh_subelement(element, mesh_memo)
        self._create_resonance_scattering_subelement(element)
        self._create_volume_calcs_subelement(element)
        self._create_create_fission_neutrons_subelement(element)
        self._create_create_delayed_neutrons_subelement(element)
        self._create_delayed_photon_scaling_subelement(element)
        self._create_event_based_subelement(element)
        self._create_max_particles_in_flight_subelement(element)
        self._create_max_events_subelement(element)
        self._create_material_cell_offsets_subelement(element)
        self._create_log_grid_bins_subelement(element)
        self._create_write_initial_source_subelement(element)
        self._create_weight_windows_subelement(element, mesh_memo)
        self._create_weight_window_generators_subelement(element, mesh_memo)
        self._create_weight_windows_file_element(element)
        self._create_weight_window_checkpoints_subelement(element)
        self._create_max_history_splits_subelement(element)
        self._create_max_tracks_subelement(element)
        self._create_random_ray_subelement(element)

        # Clean the indentation in the file to be user-readable
        clean_indentation(element)
        reorder_attributes(element)  # TODO: Remove when support is Python 3.8+

        return element

    def export_to_xml(self, path: PathLike = 'settings.xml'):
        """Export simulation settings to an XML file.

        Parameters
        ----------
        path : str
            Path to file to write. Defaults to 'settings.xml'.

        """
        root_element = self.to_xml_element()

        # Check if path is a directory
        p = Path(path)
        if p.is_dir():
            p /= 'settings.xml'

        # Write the XML Tree to the settings.xml file
        tree = ET.ElementTree(root_element)
        tree.write(str(p), xml_declaration=True, encoding='utf-8')

    @classmethod
    def from_xml_element(cls, elem, meshes=None):
        """Generate settings from XML element

        Parameters
        ----------
        elem : lxml.etree._Element
            XML element
        meshes : dict or None
            A dictionary with mesh IDs as keys and mesh instances as values that
            have already been read from XML. Pre-existing meshes are used
            and new meshes are added to when creating tally objects.

        Returns
        -------
        openmc.Settings
            Settings object

        """
        # read all meshes under the settings node and update
        settings_meshes = _read_meshes(elem)
        meshes = {} if meshes is None else meshes
        meshes.update(settings_meshes)

        settings = cls()
        settings._eigenvalue_from_xml_element(elem)
        settings._run_mode_from_xml_element(elem)
        settings._particles_from_xml_element(elem)
        settings._batches_from_xml_element(elem)
        settings._inactive_from_xml_element(elem)
        settings._max_lost_particles_from_xml_element(elem)
        settings._rel_max_lost_particles_from_xml_element(elem)
        settings._max_write_lost_particles_from_xml_element(elem)
        settings._generations_per_batch_from_xml_element(elem)
        settings._keff_trigger_from_xml_element(elem)
        settings._source_from_xml_element(elem, meshes)
        settings._volume_calcs_from_xml_element(elem)
        settings._output_from_xml_element(elem)
        settings._statepoint_from_xml_element(elem)
        settings._sourcepoint_from_xml_element(elem)
        settings._surf_source_read_from_xml_element(elem)
        settings._surf_source_write_from_xml_element(elem)
        settings._confidence_intervals_from_xml_element(elem)
        settings._electron_treatment_from_xml_element(elem)
        settings._energy_mode_from_xml_element(elem)
        settings._max_order_from_xml_element(elem)
        settings._photon_transport_from_xml_element(elem)
        settings._plot_seed_from_xml_element(elem)
        settings._ptables_from_xml_element(elem)
        settings._seed_from_xml_element(elem)
        settings._survival_biasing_from_xml_element(elem)
        settings._cutoff_from_xml_element(elem)
        settings._entropy_mesh_from_xml_element(elem, meshes)
        settings._trigger_from_xml_element(elem)
        settings._no_reduce_from_xml_element(elem)
        settings._verbosity_from_xml_element(elem)
        settings._tabular_legendre_from_xml_element(elem)
        settings._temperature_from_xml_element(elem)
        settings._trace_from_xml_element(elem)
        settings._track_from_xml_element(elem)
        settings._ufs_mesh_from_xml_element(elem, meshes)
        settings._resonance_scattering_from_xml_element(elem)
        settings._create_fission_neutrons_from_xml_element(elem)
        settings._create_delayed_neutrons_from_xml_element(elem)
        settings._delayed_photon_scaling_from_xml_element(elem)
        settings._event_based_from_xml_element(elem)
        settings._max_particles_in_flight_from_xml_element(elem)
        settings._max_particle_events_from_xml_element(elem)
        settings._material_cell_offsets_from_xml_element(elem)
        settings._log_grid_bins_from_xml_element(elem)
        settings._write_initial_source_from_xml_element(elem)
        settings._weight_windows_from_xml_element(elem, meshes)
        settings._weight_window_generators_from_xml_element(elem, meshes)
        settings._weight_window_checkpoints_from_xml_element(elem)
        settings._max_history_splits_from_xml_element(elem)
        settings._max_tracks_from_xml_element(elem)
        settings._random_ray_from_xml_element(elem)

        # TODO: Get volume calculations
        return settings

    @classmethod
    def from_xml(cls, path: PathLike = 'settings.xml'):
        """Generate settings from XML file

        .. versionadded:: 0.13.0

        Parameters
        ----------
        path : str, optional
            Path to settings XML file

        Returns
        -------
        openmc.Settings
            Settings object

        """
        parser = ET.XMLParser(huge_tree=True)
        tree = ET.parse(path, parser=parser)
        root = tree.getroot()
        meshes = _read_meshes(root)
        return cls.from_xml_element(root, meshes)
