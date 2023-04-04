import os
import typing  # imported separately as py3.8 requires typing.Iterable
from collections.abc import Iterable, Mapping, MutableSequence
from enum import Enum
import itertools
from math import ceil
from numbers import Integral, Real
from pathlib import Path
import typing  # required to prevent typing.Union namespace overwriting Union
from typing import Optional
from xml.etree import ElementTree as ET

import openmc.checkvalue as cv
from openmc.stats.multivariate import MeshSpatial

from . import RegularMesh, Source, VolumeCalculation, WeightWindows
from ._xml import clean_indentation, get_text, reorder_attributes, xmlinator, xml_attribute, xml_element, optional_xml_attribute, optional_xml_element
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
        Dictionary defining weight cutoff and energy cutoff. The dictionary may
        have six keys, 'weight', 'weight_avg', 'energy_neutron', 'energy_photon',
        'energy_electron', and 'energy_positron'. Value for 'weight'
        should be a float indicating weight cutoff below which particle undergo
        Russian roulette. Value for 'weight_avg' should be a float indicating
        weight assigned to particles that are not killed after Russian
        roulette. Value of energy should be a float indicating energy in eV
        below which particle type will be killed.
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
        Maximum number of lost particles, relative to the total number of particles

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
    max_order : None or int
        Maximum scattering order to apply globally when in multi-group mode.
    max_splits : int
        Maximum number of times a particle can split during a history

        .. versionadded:: 0.13
    max_tracks : int
        Maximum number of tracks written to a track file (per MPI process).

        .. versionadded:: 0.13.1
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
    ptables : bool
        Determine whether probability tables are used.
    resonance_scattering : dict
        Settings for resonance elastic scattering. Accepted keys are 'enable'
        (bool), 'method' (str), 'energy_min' (float), 'energy_max' (float), and
        'nuclides' (list). The 'method' can be set to 'dbrc' (Doppler broadening
        rejection correction) or 'rvs' (relative velocity sampling). If not
        specified, 'rvs' is the default method. The 'energy_min' and
        'energy_max' values indicate the minimum and maximum energies above and
        below which the resonance elastic scattering method is to be
        applied. The 'nuclides' list indicates what nuclides the method should
        be applied to. In its absence, the method will be applied to all
        nuclides with 0 K elastic scattering data present.
    run_mode : {'eigenvalue', 'fixed source', 'plot', 'volume', 'particle restart'}
        The type of calculation to perform (default is 'eigenvalue')
    seed : int
        Seed for the linear congruential pseudorandom number generator
    source : Iterable of openmc.Source
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

        :batches: list of batches at which to write source
    surf_source_read : dict
        Options for reading surface source points. Acceptable keys are:

        :path: Path to surface source file (str).
    surf_source_write : dict
        Options for writing surface source points. Acceptable keys are:

        :surface_ids: List of surface ids at which crossing particles are to be
                   banked (int)
        :max_particles: Maximum number of particles to be banked on
                   surfaces per process (int)
        :mcpl: Output in the form of an MCPL-file (bool)
    survival_biasing : bool
        Indicate whether survival biasing is to be used
    tabular_legendre : dict
        Determines if a multi-group scattering moment kernel expanded via
        Legendre polynomials is to be converted to a tabular distribution or
        not. Accepted keys are 'enable' and 'num_points'. The value for
        'enable' is a bool stating whether the conversion to tabular is
        performed; the value for 'num_points' sets the number of points to use
        in the tabular distribution, should 'enable' be True.
    temperature : dict
        Defines a default temperature and method for treating intermediate
        temperatures at which nuclear data doesn't exist. Accepted keys are
        'default', 'method', 'range', 'tolerance', and 'multipole'. The value
        for 'default' should be a float representing the default temperature in
        Kelvin. The value for 'method' should be 'nearest' or 'interpolation'.
        If the method is 'nearest', 'tolerance' indicates a range of
        temperature within which cross sections may be used. If the method is
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
    weight_windows : WeightWindows iterable of WeightWindows
        Weight windows to use for variance reduction

        .. versionadded:: 0.13
    create_delayed_neutrons : bool
        Whether delayed neutrons are created in fission.

        .. versionadded:: 0.13.3
    weight_windows_on : bool
        Whether weight windows are enabled

        .. versionadded:: 0.13
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
        self._particles = None
        self._keff_trigger = None

        # Energy mode subelement
        self._energy_mode = None
        self._max_order = None

        # Source subelement
        self._source = cv.CheckedList(Source, 'source distributions')

        self._confidence_intervals = None
        self._electron_treatment = None
        self._photon_transport = None
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
        self._write_initial_source = None
        self._weight_windows = cv.CheckedList(WeightWindows, 'weight windows')
        self._weight_windows_on = None
        self._max_splits = None
        self._max_tracks = None

        for key, value in kwargs.items():
            setattr(self, key, value)

    @property
    @optional_xml_element
    def run_mode(self) -> str:
        return self._run_mode.value

    @property
    @optional_xml_element
    def batches(self) -> int:
        return self._batches

    @property
    @optional_xml_element
    def generations_per_batch(self) -> int:
        return self._generations_per_batch

    @property
    @optional_xml_element
    def inactive(self) -> int:
        return self._inactive

    @property
    @optional_xml_element
    def max_lost_particles(self) -> int:
        return self._max_lost_particles

    @property
    @optional_xml_element
    def rel_max_lost_particles(self) -> float:
        return self._rel_max_lost_particles

    @property
    @optional_xml_element
    def particles(self) -> int:
        return self._particles

    @property
    @optional_xml_element
    def keff_trigger(self) -> dict:
        return self._keff_trigger

    @property
    @optional_xml_element
    def energy_mode(self) -> str:
        return self._energy_mode

    @property
    @optional_xml_element
    def max_order(self) -> int:
        return self._max_order

    @property
    @optional_xml_element
    def source(self) -> typing.List[Source]:
        return self._source

    @property
    @optional_xml_element
    def confidence_intervals(self) -> bool:
        return self._confidence_intervals

    @property
    @optional_xml_element
    def electron_treatment(self) -> str:
        return self._electron_treatment

    @property
    @optional_xml_element
    def ptables(self) -> bool:
        return self._ptables

    @property
    @optional_xml_element
    def photon_transport(self) -> bool:
        return self._photon_transport

    @property
    @optional_xml_element
    def seed(self) -> int:
        return self._seed

    @property
    @optional_xml_element
    def survival_biasing(self) -> bool:
        return self._survival_biasing

    @property
    @optional_xml_element
    def entropy_mesh(self) -> RegularMesh:
        return self._entropy_mesh

    @property
    @optional_xml_element
    def trigger_active(self) -> bool:
        return self._trigger_active

    @property
    @optional_xml_element
    def trigger_max_batches(self) -> int:
        return self._trigger_max_batches

    @property
    @optional_xml_element
    def trigger_batch_interval(self) -> int:
        return self._trigger_batch_interval

    @property
    @optional_xml_element
    def output(self) -> dict:
        return self._output

    @property
    @optional_xml_element
    def sourcepoint(self) -> dict:
        return self._sourcepoint

    @property
    @optional_xml_element
    def statepoint(self) -> dict:
        return self._statepoint

    @property
    @optional_xml_element
    def surf_source_read(self) -> dict:
        return self._surf_source_read

    @property
    @optional_xml_element
    def surf_source_write(self) -> dict:
        return self._surf_source_write

    @property
    @optional_xml_element
    def no_reduce(self) -> bool:
        return self._no_reduce

    @property
    @optional_xml_element
    def verbosity(self) -> int:
        return self._verbosity

    @property
    @optional_xml_element
    def tabular_legendre(self) -> dict:
        return self._tabular_legendre

    @property
    @optional_xml_element
    def temperature(self) -> dict:
        return self._temperature

    @property
    @optional_xml_element
    def trace(self) -> typing.Iterable:
        return self._trace

    @property
    @optional_xml_element
    def track(self) -> typing.Iterable[typing.Iterable[int]]:
        return self._track

    @property
    @optional_xml_element
    def cutoff(self) -> dict:
        return self._cutoff

    @property
    @optional_xml_element
    def ufs_mesh(self) -> RegularMesh:
        return self._ufs_mesh

    @property
    @optional_xml_element
    def resonance_scattering(self) -> dict:
        return self._resonance_scattering

    @property
    @optional_xml_element
    def volume_calculations(self) -> typing.List[VolumeCalculation]:
        return self._volume_calculations

    @property
    @optional_xml_element
    def create_fission_neutrons(self) -> bool:
        return self._create_fission_neutrons

    @property
    @optional_xml_element
    def create_delayed_neutrons(self) -> bool:
        return self._create_delayed_neutrons

    @property
    @optional_xml_element
    def delayed_photon_scaling(self) -> bool:
        return self._delayed_photon_scaling

    @property
    @optional_xml_element
    def material_cell_offsets(self) -> bool:
        return self._material_cell_offsets

    @property
    @optional_xml_element
    def log_grid_bins(self) -> int:
        return self._log_grid_bins

    @property
    @optional_xml_element
    def event_based(self) -> bool:
        return self._event_based

    @property
    @optional_xml_element
    def max_particles_in_flight(self) -> int:
        return self._max_particles_in_flight

    @property
    @optional_xml_element
    def write_initial_source(self) -> bool:
        return self._write_initial_source

    @property
    @optional_xml_element
    def weight_windows(self) -> typing.List[WeightWindows]:
        return self._weight_windows

    @property
    @optional_xml_element
    def weight_windows_on(self) -> bool:
        return self._weight_windows_on

    @property
    @optional_xml_element
    def max_splits(self) -> int:
        return self._max_splits

    @property
    @optional_xml_element
    def max_tracks(self) -> int:
        return self._max_tracks

    @run_mode.setter
    def run_mode(self, run_mode: str):
        cv.check_value('run mode', run_mode, {x.value for x in RunMode})
        for mode in RunMode:
            if mode.value == run_mode:
                self._run_mode = mode

    @batches.setter
    def batches(self, batches: int):
        cv.check_type('batches', batches, Integral)
        cv.check_greater_than('batches', batches, 0)
        self._batches = batches

    @generations_per_batch.setter
    def generations_per_batch(self, generations_per_batch: int):
        cv.check_type('generations per patch', generations_per_batch, Integral)
        cv.check_greater_than('generations per batch',
                              generations_per_batch, 0)
        self._generations_per_batch = generations_per_batch

    @inactive.setter
    def inactive(self, inactive: int):
        cv.check_type('inactive batches', inactive, Integral)
        cv.check_greater_than('inactive batches', inactive, 0, True)
        self._inactive = inactive

    @max_lost_particles.setter
    def max_lost_particles(self, max_lost_particles: int):
        cv.check_type('max_lost_particles', max_lost_particles, Integral)
        cv.check_greater_than('max_lost_particles', max_lost_particles, 0)
        self._max_lost_particles = max_lost_particles

    @rel_max_lost_particles.setter
    def rel_max_lost_particles(self, rel_max_lost_particles: float):
        cv.check_type('rel_max_lost_particles', rel_max_lost_particles, Real)
        cv.check_greater_than('rel_max_lost_particles',
                              rel_max_lost_particles, 0)
        cv.check_less_than('rel_max_lost_particles', rel_max_lost_particles, 1)
        self._rel_max_lost_particles = rel_max_lost_particles

    @particles.setter
    def particles(self, particles: int):
        cv.check_type('particles', particles, Integral)
        cv.check_greater_than('particles', particles, 0)
        self._particles = particles

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

    @energy_mode.setter
    def energy_mode(self, energy_mode: str):
        cv.check_value('energy mode', energy_mode,
                       ['continuous-energy', 'multi-group'])
        self._energy_mode = energy_mode

    @max_order.setter
    def max_order(self, max_order: Optional[int]):
        if max_order is not None:
            cv.check_type('maximum scattering order', max_order, Integral)
            cv.check_greater_than('maximum scattering order', max_order, 0,
                                  True)
        self._max_order = max_order

    @source.setter
    def source(self, source: typing.Union[Source, typing.Iterable[Source]]):
        if not isinstance(source, MutableSequence):
            source = [source]
        self._source = cv.CheckedList(Source, 'source distributions', source)

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

    @verbosity.setter
    def verbosity(self, verbosity: int):
        cv.check_type('verbosity', verbosity, Integral)
        cv.check_greater_than('verbosity', verbosity, 1, True)
        cv.check_less_than('verbosity', verbosity, 10, True)
        self._verbosity = verbosity

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

    @surf_source_read.setter
    def surf_source_read(self, surf_source_read: dict):
        cv.check_type('surface source reading options',
                      surf_source_read, Mapping)
        for key, value in surf_source_read.items():
            cv.check_value('surface source reading key', key,
                           ('path'))
            if key == 'path':
                cv.check_type('path to surface source file', value, str)
        self._surf_source_read = surf_source_read

    @surf_source_write.setter
    def surf_source_write(self, surf_source_write: dict):
        cv.check_type('surface source writing options',
                      surf_source_write, Mapping)
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

    @confidence_intervals.setter
    def confidence_intervals(self, confidence_intervals: bool):
        cv.check_type('confidence interval', confidence_intervals, bool)
        self._confidence_intervals = confidence_intervals

    @electron_treatment.setter
    def electron_treatment(self, electron_treatment: str):
        cv.check_value('electron treatment',
                       electron_treatment, ['led', 'ttb'])
        self._electron_treatment = electron_treatment

    @photon_transport.setter
    def photon_transport(self, photon_transport: bool):
        cv.check_type('photon transport', photon_transport, bool)
        self._photon_transport = photon_transport

    @ptables.setter
    def ptables(self, ptables: bool):
        cv.check_type('probability tables', ptables, bool)
        self._ptables = ptables

    @seed.setter
    def seed(self, seed: int):
        cv.check_type('random number generator seed', seed, Integral)
        cv.check_greater_than('random number generator seed', seed, 0)
        self._seed = seed

    @survival_biasing.setter
    def survival_biasing(self, survival_biasing: bool):
        cv.check_type('survival biasing', survival_biasing, bool)
        self._survival_biasing = survival_biasing

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

    @entropy_mesh.setter
    def entropy_mesh(self, entropy: RegularMesh):
        cv.check_type('entropy mesh', entropy, RegularMesh)
        self._entropy_mesh = entropy

    @trigger_active.setter
    def trigger_active(self, trigger_active: bool):
        cv.check_type('trigger active', trigger_active, bool)
        self._trigger_active = trigger_active

    @trigger_max_batches.setter
    def trigger_max_batches(self, trigger_max_batches: int):
        cv.check_type('trigger maximum batches', trigger_max_batches, Integral)
        cv.check_greater_than('trigger maximum batches',
                              trigger_max_batches, 0)
        self._trigger_max_batches = trigger_max_batches

    @trigger_batch_interval.setter
    def trigger_batch_interval(self, trigger_batch_interval: int):
        cv.check_type('trigger batch interval',
                      trigger_batch_interval, Integral)
        cv.check_greater_than('trigger batch interval',
                              trigger_batch_interval, 0)
        self._trigger_batch_interval = trigger_batch_interval

    @no_reduce.setter
    def no_reduce(self, no_reduce: bool):
        cv.check_type('no reduction option', no_reduce, bool)
        self._no_reduce = no_reduce

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

    @trace.setter
    def trace(self, trace: Iterable):
        cv.check_type('trace', trace, Iterable, Integral)
        cv.check_length('trace', trace, 3)
        cv.check_greater_than('trace batch', trace[0], 0)
        cv.check_greater_than('trace generation', trace[1], 0)
        cv.check_greater_than('trace particle', trace[2], 0)
        self._trace = trace

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

    @ufs_mesh.setter
    def ufs_mesh(self, ufs_mesh: RegularMesh):
        cv.check_type('UFS mesh', ufs_mesh, RegularMesh)
        cv.check_length('UFS mesh dimension', ufs_mesh.dimension, 3)
        cv.check_length('UFS mesh lower-left corner', ufs_mesh.lower_left, 3)
        cv.check_length('UFS mesh upper-right corner', ufs_mesh.upper_right, 3)
        self._ufs_mesh = ufs_mesh

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

    @volume_calculations.setter
    def volume_calculations(
        self, vol_calcs: typing.Union[VolumeCalculation, typing.Iterable[VolumeCalculation]]
    ):
        if not isinstance(vol_calcs, MutableSequence):
            vol_calcs = [vol_calcs]
        self._volume_calculations = cv.CheckedList(
            VolumeCalculation, 'stochastic volume calculations', vol_calcs)

    @create_fission_neutrons.setter
    def create_fission_neutrons(self, create_fission_neutrons: bool):
        cv.check_type('Whether create fission neutrons',
                      create_fission_neutrons, bool)
        self._create_fission_neutrons = create_fission_neutrons

    @create_delayed_neutrons.setter
    def create_delayed_neutrons(self, create_delayed_neutrons: bool):
        cv.check_type('Whether create only prompt neutrons',
                      create_delayed_neutrons, bool)
        self._create_delayed_neutrons = create_delayed_neutrons

    @delayed_photon_scaling.setter
    def delayed_photon_scaling(self, value: bool):
        cv.check_type('delayed photon scaling', value, bool)
        self._delayed_photon_scaling = value

    @event_based.setter
    def event_based(self, value: bool):
        cv.check_type('event based', value, bool)
        self._event_based = value

    @max_particles_in_flight.setter
    def max_particles_in_flight(self, value: int):
        cv.check_type('max particles in flight', value, Integral)
        cv.check_greater_than('max particles in flight', value, 0)
        self._max_particles_in_flight = value

    @material_cell_offsets.setter
    def material_cell_offsets(self, value: bool):
        cv.check_type('material cell offsets', value, bool)
        self._material_cell_offsets = value

    @log_grid_bins.setter
    def log_grid_bins(self, log_grid_bins: int):
        cv.check_type('log grid bins', log_grid_bins, Real)
        cv.check_greater_than('log grid bins', log_grid_bins, 0)
        self._log_grid_bins = log_grid_bins

    @write_initial_source.setter
    def write_initial_source(self, value: bool):
        cv.check_type('write initial source', value, bool)
        self._write_initial_source = value

    @weight_windows.setter
    def weight_windows(self, value: typing.Union[WeightWindows, typing.Iterable[WeightWindows]]):
        if not isinstance(value, MutableSequence):
            value = [value]
        self._weight_windows = cv.CheckedList(
            WeightWindows, 'weight windows', value)

    @weight_windows_on.setter
    def weight_windows_on(self, value):
        cv.check_type('weight windows on', value, bool)
        self._weight_windows_on = value

    @max_splits.setter
    def max_splits(self, value: int):
        cv.check_type('maximum particle splits', value, Integral)
        cv.check_greater_than('max particle splits', value, 0)
        self._max_splits = value

    @max_tracks.setter
    def max_tracks(self, value: int):
        cv.check_type('maximum particle tracks', value, Integral)
        cv.check_greater_than('maximum particle tracks', value, 0, True)
        self._max_tracks = value

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
        elem : xml.etree.ElementTree.Element
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
        settings = cls()
        settings._eigenvalue_from_xml_element(elem)
        settings._run_mode_from_xml_element(elem)
        settings._particles_from_xml_element(elem)
        settings._batches_from_xml_element(elem)
        settings._inactive_from_xml_element(elem)
        settings._max_lost_particles_from_xml_element(elem)
        settings._rel_max_lost_particles_from_xml_element(elem)
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
        settings._material_cell_offsets_from_xml_element(elem)
        settings._log_grid_bins_from_xml_element(elem)
        settings._write_initial_source_from_xml_element(elem)
        settings._weight_windows_from_xml_element(elem, meshes)
        settings._max_splits_from_xml_element(elem)
        settings._max_tracks_from_xml_element(elem)

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
        tree = ET.parse(path)
        root = tree.getroot()
        meshes = _read_meshes(root)
        return cls.from_xml_element(root, meshes)
