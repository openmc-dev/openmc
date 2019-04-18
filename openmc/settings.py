from collections.abc import Iterable, MutableSequence, Mapping
from pathlib import Path
from numbers import Real, Integral
import warnings
from xml.etree import ElementTree as ET
import sys

import numpy as np

from openmc._xml import clean_indentation
import openmc.checkvalue as cv
from openmc import VolumeCalculation, Source, Mesh

_RUN_MODES = ['eigenvalue', 'fixed source', 'plot', 'volume', 'particle restart']
_RES_SCAT_METHODS = ['dbrc', 'rvs']


class Settings(object):
    """Settings used for an OpenMC simulation.

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
    dagmc : bool
        Indicate that a CAD-based DAGMC geometry will be used.
    electron_treatment : {'led', 'ttb'}
        Whether to deposit all energy from electrons locally ('led') or create
        secondary bremsstrahlung photons ('ttb').
    energy_mode : {'continuous-energy', 'multi-group'}
        Set whether the calculation should be continuous-energy or multi-group.
    entropy_mesh : openmc.Mesh
        Mesh to be used to calculate Shannon entropy. If the mesh dimensions are
        not specified. OpenMC assigns a mesh such that 20 source sites per mesh
        cell are to be expected on average.
    generations_per_batch : int
        Number of generations per batch
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
    max_order : None or int
        Maximum scattering order to apply globally when in multi-group mode.
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
    statepoint : dict
        Options for writing state points. Acceptable keys are:

        :batches: list of batches at which to write source
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
        If the method is 'nearest', 'tolerance' indicates a range of temperature
        within which cross sections may be used. The value for 'range' should be
        a pair a minimum and maximum temperatures which are used to indicate
        that cross sections be loaded at all temperatures within the
        range. 'multipole' is a boolean indicating whether or not the windowed
        multipole method should be used to evaluate resolved resonance cross
        sections.
    trace : tuple or list
        Show detailed information about a single particle, indicated by three
        integers: the batch number, generation number, and particle number
    track : tuple or list
        Specify particles for which track files should be written. Each particle
        is identified by a triplet with the batch number, generation number, and
        particle number.
    trigger_active : bool
        Indicate whether tally triggers are used
    trigger_batch_interval : int
        Number of batches in between convergence checks
    trigger_max_batches : int
        Maximum number of batches simulated. If this is set, the number of
        batches specified via ``batches`` is interpreted as the minimum number
        of batches
    ufs_mesh : openmc.Mesh
        Mesh to be used for redistributing source sites via the uniform fision
        site (UFS) method.
    verbosity : int
        Verbosity during simulation between 1 and 10. Verbosity levels are
        described in :ref:`verbosity`.
    volume_calculations : VolumeCalculation or iterable of VolumeCalculation
        Stochastic volume calculation specifications

    """

    def __init__(self):

        # Run mode subelement (default is 'eigenvalue')
        self._run_mode = 'eigenvalue'
        self._batches = None
        self._generations_per_batch = None
        self._inactive = None
        self._particles = None
        self._keff_trigger = None

        # Energy mode subelement
        self._energy_mode = None
        self._max_order = None

        # Source subelement
        self._source = cv.CheckedList(Source, 'source distributions')

        self._confidence_intervals = None
        self._cross_sections = None
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
        self._log_grid_bins = None

        self._dagmc = False

    @property
    def run_mode(self):
        return self._run_mode

    @property
    def batches(self):
        return self._batches

    @property
    def generations_per_batch(self):
        return self._generations_per_batch

    @property
    def inactive(self):
        return self._inactive

    @property
    def particles(self):
        return self._particles

    @property
    def keff_trigger(self):
        return self._keff_trigger

    @property
    def energy_mode(self):
        return self._energy_mode

    @property
    def max_order(self):
        return self._max_order

    @property
    def source(self):
        return self._source

    @property
    def confidence_intervals(self):
        return self._confidence_intervals

    @property
    def electron_treatment(self):
        return self._electron_treatment

    @property
    def ptables(self):
        return self._ptables

    @property
    def photon_transport(self):
        return self._photon_transport

    @property
    def seed(self):
        return self._seed

    @property
    def survival_biasing(self):
        return self._survival_biasing

    @property
    def entropy_mesh(self):
        return self._entropy_mesh

    @property
    def trigger_active(self):
        return self._trigger_active

    @property
    def trigger_max_batches(self):
        return self._trigger_max_batches

    @property
    def trigger_batch_interval(self):
        return self._trigger_batch_interval

    @property
    def output(self):
        return self._output

    @property
    def sourcepoint(self):
        return self._sourcepoint

    @property
    def statepoint(self):
        return self._statepoint

    @property
    def no_reduce(self):
        return self._no_reduce

    @property
    def verbosity(self):
        return self._verbosity

    @property
    def tabular_legendre(self):
        return self._tabular_legendre

    @property
    def temperature(self):
        return self._temperature

    @property
    def trace(self):
        return self._trace

    @property
    def track(self):
        return self._track

    @property
    def cutoff(self):
        return self._cutoff

    @property
    def ufs_mesh(self):
        return self._ufs_mesh

    @property
    def resonance_scattering(self):
        return self._resonance_scattering

    @property
    def volume_calculations(self):
        return self._volume_calculations

    @property
    def create_fission_neutrons(self):
        return self._create_fission_neutrons

    @property
    def log_grid_bins(self):
        return self._log_grid_bins

    @property
    def dagmc(self):
        return self._dagmc

    @run_mode.setter
    def run_mode(self, run_mode):
        cv.check_value('run mode', run_mode, _RUN_MODES)
        self._run_mode = run_mode

    @batches.setter
    def batches(self, batches):
        cv.check_type('batches', batches, Integral)
        cv.check_greater_than('batches', batches, 0)
        self._batches = batches

    @generations_per_batch.setter
    def generations_per_batch(self, generations_per_batch):
        cv.check_type('generations per patch', generations_per_batch, Integral)
        cv.check_greater_than('generations per batch', generations_per_batch, 0)
        self._generations_per_batch = generations_per_batch

    @inactive.setter
    def inactive(self, inactive):
        cv.check_type('inactive batches', inactive, Integral)
        cv.check_greater_than('inactive batches', inactive, 0, True)
        self._inactive = inactive

    @particles.setter
    def particles(self, particles):
        cv.check_type('particles', particles, Integral)
        cv.check_greater_than('particles', particles, 0)
        self._particles = particles

    @keff_trigger.setter
    def keff_trigger(self, keff_trigger):
        if not isinstance(keff_trigger, dict):
            msg = 'Unable to set a trigger on keff from "{0}" which ' \
                  'is not a Python dictionary'.format(keff_trigger)
            raise ValueError(msg)

        elif 'type' not in keff_trigger:
            msg = 'Unable to set a trigger on keff from "{0}" which ' \
                  'does not have a "type" key'.format(keff_trigger)
            raise ValueError(msg)

        elif keff_trigger['type'] not in ['variance', 'std_dev', 'rel_err']:
            msg = 'Unable to set a trigger on keff with ' \
                  'type "{0}"'.format(keff_trigger['type'])
            raise ValueError(msg)

        elif 'threshold' not in keff_trigger:
            msg = 'Unable to set a trigger on keff from "{0}" which ' \
                  'does not have a "threshold" key'.format(keff_trigger)
            raise ValueError(msg)

        elif not isinstance(keff_trigger['threshold'], Real):
            msg = 'Unable to set a trigger on keff with ' \
                  'threshold "{0}"'.format(keff_trigger['threshold'])
            raise ValueError(msg)

        self._keff_trigger = keff_trigger

    @energy_mode.setter
    def energy_mode(self, energy_mode):
        cv.check_value('energy mode', energy_mode,
                    ['continuous-energy', 'multi-group'])
        self._energy_mode = energy_mode

    @max_order.setter
    def max_order(self, max_order):
        if max_order is not None:
            cv.check_type('maximum scattering order', max_order, Integral)
            cv.check_greater_than('maximum scattering order', max_order, 0,
                                  True)
        self._max_order = max_order

    @source.setter
    def source(self, source):
        if not isinstance(source, MutableSequence):
            source = [source]
        self._source = cv.CheckedList(Source, 'source distributions', source)

    @output.setter
    def output(self, output):
        cv.check_type('output', output, Mapping)
        for key, value in output.items():
            cv.check_value('output key', key, ('summary', 'tallies', 'path'))
            if key in ('summary', 'tallies'):
                cv.check_type("output['{}']".format(key), value, bool)
            else:
                cv.check_type("output['path']", value, str)
        self._output = output

    @verbosity.setter
    def verbosity(self, verbosity):
        cv.check_type('verbosity', verbosity, Integral)
        cv.check_greater_than('verbosity', verbosity, 1, True)
        cv.check_less_than('verbosity', verbosity, 10, True)
        self._verbosity = verbosity

    @sourcepoint.setter
    def sourcepoint(self, sourcepoint):
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
            else:
                raise ValueError("Unknown key '{}' encountered when setting "
                                 "sourcepoint options.".format(key))
        self._sourcepoint = sourcepoint

    @statepoint.setter
    def statepoint(self, statepoint):
        cv.check_type('statepoint options', statepoint, Mapping)
        for key, value in statepoint.items():
            if key == 'batches':
                cv.check_type('statepoint batches', value, Iterable, Integral)
                for batch in value:
                    cv.check_greater_than('statepoint batch', batch, 0)
            else:
                raise ValueError("Unknown key '{}' encountered when setting "
                                 "statepoint options.".format(key))
        self._statepoint = statepoint

    @confidence_intervals.setter
    def confidence_intervals(self, confidence_intervals):
        cv.check_type('confidence interval', confidence_intervals, bool)
        self._confidence_intervals = confidence_intervals

    @electron_treatment.setter
    def electron_treatment(self, electron_treatment):
        cv.check_value('electron treatment', electron_treatment, ['led', 'ttb'])
        self._electron_treatment = electron_treatment

    @photon_transport.setter
    def photon_transport(self, photon_transport):
        cv.check_type('photon transport', photon_transport, bool)
        self._photon_transport = photon_transport

    @dagmc.setter
    def dagmc(self, dagmc):
        cv.check_type('dagmc geometry', dagmc, bool)
        self._dagmc = dagmc

    @ptables.setter
    def ptables(self, ptables):
        cv.check_type('probability tables', ptables, bool)
        self._ptables = ptables

    @seed.setter
    def seed(self, seed):
        cv.check_type('random number generator seed', seed, Integral)
        cv.check_greater_than('random number generator seed', seed, 0)
        self._seed = seed

    @survival_biasing.setter
    def survival_biasing(self, survival_biasing):
        cv.check_type('survival biasing', survival_biasing, bool)
        self._survival_biasing = survival_biasing

    @cutoff.setter
    def cutoff(self, cutoff):
        if not isinstance(cutoff, Mapping):
            msg = 'Unable to set cutoff from "{0}" which is not a '\
                  ' Python dictionary'.format(cutoff)
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
                msg = 'Unable to set cutoff to "{0}" which is unsupported by '\
                      'OpenMC'.format(key)

        self._cutoff = cutoff

    @entropy_mesh.setter
    def entropy_mesh(self, entropy):
        cv.check_type('entropy mesh', entropy, Mesh)
        cv.check_length('entropy mesh dimension', entropy.dimension, 3)
        cv.check_length('entropy mesh lower-left corner', entropy.lower_left, 3)
        cv.check_length('entropy mesh upper-right corner', entropy.upper_right, 3)
        self._entropy_mesh = entropy

    @trigger_active.setter
    def trigger_active(self, trigger_active):
        cv.check_type('trigger active', trigger_active, bool)
        self._trigger_active = trigger_active

    @trigger_max_batches.setter
    def trigger_max_batches(self, trigger_max_batches):
        cv.check_type('trigger maximum batches', trigger_max_batches, Integral)
        cv.check_greater_than('trigger maximum batches', trigger_max_batches, 0)
        self._trigger_max_batches = trigger_max_batches

    @trigger_batch_interval.setter
    def trigger_batch_interval(self, trigger_batch_interval):
        cv.check_type('trigger batch interval', trigger_batch_interval, Integral)
        cv.check_greater_than('trigger batch interval', trigger_batch_interval, 0)
        self._trigger_batch_interval = trigger_batch_interval

    @no_reduce.setter
    def no_reduce(self, no_reduce):
        cv.check_type('no reduction option', no_reduce, bool)
        self._no_reduce = no_reduce

    @tabular_legendre.setter
    def tabular_legendre(self, tabular_legendre):
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
    def temperature(self, temperature):

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
    def trace(self, trace):
        cv.check_type('trace', trace, Iterable, Integral)
        cv.check_length('trace', trace, 3)
        cv.check_greater_than('trace batch', trace[0], 0)
        cv.check_greater_than('trace generation', trace[1], 0)
        cv.check_greater_than('trace particle', trace[2], 0)
        self._trace = trace

    @track.setter
    def track(self, track):
        cv.check_type('track', track, Iterable, Integral)
        if len(track) % 3 != 0:
            msg = 'Unable to set the track to "{0}" since its length is ' \
                  'not a multiple of 3'.format(track)
            raise ValueError(msg)
        for t in zip(track[::3], track[1::3], track[2::3]):
            cv.check_greater_than('track batch', t[0], 0)
            cv.check_greater_than('track generation', t[0], 0)
            cv.check_greater_than('track particle', t[0], 0)
        self._track = track

    @ufs_mesh.setter
    def ufs_mesh(self, ufs_mesh):
        cv.check_type('UFS mesh', ufs_mesh, Mesh)
        cv.check_length('UFS mesh dimension', ufs_mesh.dimension, 3)
        cv.check_length('UFS mesh lower-left corner', ufs_mesh.lower_left, 3)
        cv.check_length('UFS mesh upper-right corner', ufs_mesh.upper_right, 3)
        self._ufs_mesh = ufs_mesh

    @resonance_scattering.setter
    def resonance_scattering(self, res):
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
    def volume_calculations(self, vol_calcs):
        if not isinstance(vol_calcs, MutableSequence):
            vol_calcs = [vol_calcs]
        self._volume_calculations = cv.CheckedList(
            VolumeCalculation, 'stochastic volume calculations', vol_calcs)

    @create_fission_neutrons.setter
    def create_fission_neutrons(self, create_fission_neutrons):
        cv.check_type('Whether create fission neutrons',
                      create_fission_neutrons, bool)
        self._create_fission_neutrons = create_fission_neutrons

    @log_grid_bins.setter
    def log_grid_bins(self, log_grid_bins):
        cv.check_type('log grid bins', log_grid_bins, Real)
        cv.check_greater_than('log grid bins', log_grid_bins, 0)
        self._log_grid_bins = log_grid_bins

    def _create_run_mode_subelement(self, root):
        elem = ET.SubElement(root, "run_mode")
        elem.text = self._run_mode

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

    def _create_particles_subelement(self, root):
        if self._particles is not None:
            element = ET.SubElement(root, "particles")
            element.text = str(self._particles)

    def _create_keff_trigger_subelement(self, root):
        if self._keff_trigger is not None:
            element = ET.SubElement(root, "keff_trigger")

            for key in self._keff_trigger:
                subelement = ET.SubElement(element, key)
                subelement.text = str(self._keff_trigger[key]).lower()

    def _create_energy_mode_subelement(self, root):
        if self._energy_mode is not None:
            element = ET.SubElement(root, "energy_mode")
            element.text = str(self._energy_mode)

    def _create_max_order_subelement(self, root):
        if self._max_order is not None:
            element = ET.SubElement(root, "max_order")
            element.text = str(self._max_order)

    def _create_source_subelement(self, root):
        for source in self.source:
            root.append(source.to_xml_element())

    def _create_volume_calcs_subelement(self, root):
        for calc in self.volume_calculations:
            root.append(calc.to_xml_element())

    def _create_output_subelement(self, root):
        if self._output is not None:
            element = ET.SubElement(root, "output")

            for key, value in self._output.items():
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

    def _create_entropy_mesh_subelement(self, root):
        if self.entropy_mesh is not None:
            # See if a <mesh> element already exists -- if not, add it
            path = "./mesh[@id='{}']".format(self.entropy_mesh.id)
            if root.find(path) is None:
                root.append(self.entropy_mesh.to_xml_element())

            subelement = ET.SubElement(root, "entropy_mesh")
            subelement.text = str(self.entropy_mesh.id)

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
                element = ET.SubElement(root,
                                        "temperature_{}".format(key))
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
            element.text = ' '.join(map(str, self._track))

    def _create_ufs_mesh_subelement(self, root):
        if self.ufs_mesh is not None:
            # See if a <mesh> element already exists -- if not, add it
            path = "./mesh[@id='{}']".format(self.ufs_mesh.id)
            if root.find(path) is None:
                root.append(self.ufs_mesh.to_xml_element())

            subelement = ET.SubElement(root, "ufs_mesh")
            subelement.text = str(self.ufs_mesh.id)

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

    def _create_log_grid_bins_subelement(self, root):
        if self._log_grid_bins is not None:
            elem = ET.SubElement(root, "log_grid_bins")
            elem.text = str(self._log_grid_bins)

    def _create_dagmc_subelement(self, root):
        if self._dagmc:
            elem = ET.SubElement(root, "dagmc")
            elem.text = str(self._dagmc).lower()

    def export_to_xml(self, path='settings.xml'):
        """Export simulation settings to an XML file.

        Parameters
        ----------
        path : str
            Path to file to write. Defaults to 'settings.xml'.

        """

        # Reset xml element tree
        root_element = ET.Element("settings")

        self._create_run_mode_subelement(root_element)
        self._create_particles_subelement(root_element)
        self._create_batches_subelement(root_element)
        self._create_inactive_subelement(root_element)
        self._create_generations_per_batch_subelement(root_element)
        self._create_keff_trigger_subelement(root_element)
        self._create_source_subelement(root_element)
        self._create_output_subelement(root_element)
        self._create_statepoint_subelement(root_element)
        self._create_sourcepoint_subelement(root_element)
        self._create_confidence_intervals(root_element)
        self._create_electron_treatment_subelement(root_element)
        self._create_energy_mode_subelement(root_element)
        self._create_max_order_subelement(root_element)
        self._create_photon_transport_subelement(root_element)
        self._create_ptables_subelement(root_element)
        self._create_seed_subelement(root_element)
        self._create_survival_biasing_subelement(root_element)
        self._create_cutoff_subelement(root_element)
        self._create_entropy_mesh_subelement(root_element)
        self._create_trigger_subelement(root_element)
        self._create_no_reduce_subelement(root_element)
        self._create_verbosity_subelement(root_element)
        self._create_tabular_legendre_subelements(root_element)
        self._create_temperature_subelements(root_element)
        self._create_trace_subelement(root_element)
        self._create_track_subelement(root_element)
        self._create_ufs_mesh_subelement(root_element)
        self._create_resonance_scattering_subelement(root_element)
        self._create_volume_calcs_subelement(root_element)
        self._create_create_fission_neutrons_subelement(root_element)
        self._create_log_grid_bins_subelement(root_element)
        self._create_dagmc_subelement(root_element)

        # Clean the indentation in the file to be user-readable
        clean_indentation(root_element)

        # Check if path is a directory
        p = Path(path)
        if p.is_dir():
            p /= 'settings.xml'

        # Write the XML Tree to the settings.xml file
        tree = ET.ElementTree(root_element)
        tree.write(str(p), xml_declaration=True, encoding='utf-8')

    @classmethod
    def from_xml(cls, path='settings.xml'):
        """Generate settings from XML file

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

        settings = cls()

        # Get the run mode
        elem = root.find('run_mode')
        if elem is not None:
            settings.run_mode = elem.text

        # Get number of particles
        elem = root.find('particles')
        if elem is not None:
            settings.particles = int(elem.text)

        # Get number of batches
        elem = root.find('batches')
        if elem is not None:
            settings.batches = int(elem.text)

        # Get number of inactive batches
        elem = root.find('inactive')
        if elem is not None:
            settings.inactive = int(elem.text)

        # Get number of generations per batch
        elem = root.find('generations_per_batch')
        if elem is not None:
            settings.generations_per_batch = int(elem.text)

        # Get keff trigger
        elem = root.find('keff_trigger')
        if elem is not None:
            trigger = elem.findtext('type')
            threshold = float(elem.findtext('threshold'))
            settings.keff_trigger = {'type': trigger, 'threshold': threshold}

        # Get the source
        for elem in root.findall('source'):
            settings.source.append(Source.from_xml_element(elem))

        # Get the output
        elem = root.find('output')
        if elem is not None:
            settings.output = {}
            for entry in elem:
                key = entry.tag
                if key in ('summary', 'tallies'):
                    value = entry.text == 'true'
                else:
                    value = entry.text
                settings.output[key] = value

        # Get the statepoint
        elem = root.find('state_point')
        if elem is not None:
            batches = elem.findtext('batches')
            if batches is not None:
                settings.statepoint['batches'] = [int(x) for x in batches.split()]

        # Get the sourcepoint
        elem = root.find('source_point')
        if elem is not None:
            for entry in elem:
                key = entry.tag
                if key in ('separate', 'write', 'overwrite'):
                    value = entry.text == 'true'
                else:
                    value = [int(x) for x in entry.text.split()]
                settings.sourcepoint[key] = value

        # Get confidence intervals
        elem = root.find('confidence_intervals')
        if elem is not None:
            settings.confidence_intervals = elem.text == 'true'

        # Get electron treatment
        elem = root.find('electron_treatment')
        if elem is not None:
            settings.electron_treatment = elem.text

        # Get energy mode
        elem = root.find('energy_mode')
        if elem is not None:
            settings.energy_mode = elem.text

        # Get max order
        elem = root.find('max_order')
        if elem is not None:
            settings.max_order = int(elem.text)

        # Get photon transport
        elem = root.find('photon_transport')
        if elem is not None:
            settings.photon_transport = elem.text == 'true'

        # Get probability tables
        elem = root.find('ptables')
        if elem is not None:
            settings.ptables = elem.text == 'true'

        # Get seed
        elem = root.find('seed')
        if elem is not None:
            settings.seed = int(elem.text)

        # Get survival biasing
        elem = root.find('survival_biasing')
        if elem is not None:
            settings.survival_biasing = elem.text == 'true'

        # Get cutoff
        elem = root.find('cutoff')
        if elem is not None:
            settings.cutoff = {x.tag: float(x.text) for x in elem}

        # Get entropy mesh
        elem = root.find('entropy_mesh')
        if elem is not None:
            settings.entropy_mesh = Mesh.from_xml_element(elem)

        # Get trigger
        elem = root.find('trigger')
        if elem is not None:
            active = elem.find('active')
            settings.trigger_active = active.text == 'true'
            max_batches = elem.find('max_batches')
            if max_batches is not None:
                settings.trigger_max_batches = int(max_batches.text)
            batch_interval = elem.find('batch_interval')
            if batch_interval is not None:
                settings.trigger_batch_interval = int(batch_interval.text)

        # Get no reduce
        elem = root.find('no_reduce')
        if elem is not None:
            settings.no_reduce = elem.text == 'true'

        # Get verbosity
        elem = root.find('verbosity')
        if elem is not None:
            settings.verbosity = int(elem.text)

        # Get tabular legendre
        elem = root.find('tabular_legendre')
        if elem is not None:
            enable = elem.findtext('eneable')
            settings.tabular_legendre['enable'] = enable == 'true'
            num_points = elem.findtext('num_points')
            if num_points is not None:
                settings.tabular_legendre['num_points'] = int(num_points)

        # Get temperature
        elem = root.findtext('temperature_default')
        if elem is not None:
            settings.temperature['default'] = float(elem)
        elem = root.findtext('temperature_tolerance')
        if elem is not None:
            settings.temperature['tolerance'] = float(elem)
        elem = root.findtext('temperature_method')
        if elem is not None:
            settings.temperature['method'] = elem
        elem = root.findtext('temperature_range')
        if elem is not None:
            settings.temperature['range'] = [float(x) for x in elem.split()]
        elem = root.findtext('temperature_multipole')
        if elem is not None:
            settings.temperature['multipole'] = elem == 'true'

        # Get trace
        elem = root.find('trace')
        if elem is not None:
            settings.trace = [int(x) for x in elem.text.split()]

        # Get track
        elem = root.find('track')
        if elem is not None:
            settings.track = [int(x) for x in elem.text.split()]

        # Get UFS mesh
        elem = root.find('ufs_mesh')
        if elem is not None:
            settings.ufs_mesh = Mesh.from_xml_element(elem)

        # Get resonance scattering
        elem = root.find('resonance_scattering')
        if elem is not None:
            for entry in elem:
                key = entry.tag
                if key == 'enable':
                    value = entry.text == 'true'
                elif key == 'method':
                    value = entry.text
                elif key == 'energy_min':
                    value = float(entry.text)
                elif key == 'energy_max':
                    value = float(entry.text)
                elif key == 'nuclides':
                    value = entry.text.split()
                settings.resonance_scattering[key] = value

        # TODO: Get volume calculations

        # Get fission neutrons
        elem = root.find('create_fission_neutrons')
        if elem is not None:
            settings.create_fission_neutrons = elem.text == 'true'

        # Get log grid bins
        elem = root.find('log_grid_bins')
        if elem is not None:
            settings.log_grid_bins = int(elem.text)

        # Get dagmc
        elem = root.find('dagmc')
        if elem is not None:
            settings.dagmc = elem.text == 'true'

        return settings
