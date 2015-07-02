from collections import Iterable
from numbers import Real, Integral
import warnings
from xml.etree import ElementTree as ET
import sys

import numpy as np

from openmc.clean_xml import *
from openmc.checkvalue import (check_type, check_length, check_value,
                               check_greater_than, check_less_than)

if sys.version_info[0] >= 3:
    basestring = str


class SettingsFile(object):
    """Settings file used for an OpenMC simulation. Corresponds directly to the
    settings.xml input file.

    Attributes
    ----------
    batches : int
        Number of batches to simulate
    generations_per_batch : int
        Number of generations per batch
    inactive : int
        Number of inactive batches
    particles : int
        Number of particles per generation
    keff_trigger : dict
        Dictionary defining a trigger on eigenvalue. The dictionary must have
        two keys, 'type' and 'threshold'. Acceptable values corresponding to
        type are 'variance', 'std_dev', and 'rel_err'. The threshold value
        should be a float indicating the variance, standard deviation, or
        relative error used.
    source_file : str
        Path to a source file
    output : dict
        Dictionary indicating what files to output. Valid keys are 'summary',
        'cross_sections', 'tallies', and 'distribmats'. Values corresponding to
        each key should be given as a boolean value.
    output_path : str
        Path to write output to
    verbosity : int
        Verbosity during simulation between 1 and 10
    statepoint_batches : Iterable of int
        List of batches at which to write statepoint files
    statepoint_interval : int
        Number of batches after which a new statepoint file should be written
    sourcepoint_batches : Iterable of int
        List of batches at which to write source files
    sourcepoint_interval : int
        Number of batches after which a new source file should be written
    sourcepoint_separate : bool
        Indicate whether the souce should be written as part of the statepoint
        file or on its own
    sourcepoint_write : bool
        Indicate whether the source should be written at all
    sourcepoint_overwrite : bool
        Indicate whether to
    confidence_intervals : bool
        If True, uncertainties on tally results will be reported as the
        half-width of the 95% two-sided confidence interval. If False,
        uncertainties on tally results will be reported as the sample standard
        deviation.
    cross_sections : str
        Indicates the path to an XML cross section listing file (usually named
        cross_sections.xml). If it is not set, the :envvar:`CROSS_SECTIONS`
        environment variable will be used to find the path to the XML cross
        section listing.
    energy_grid : str
        Set the method used to search energy grids. Acceptable values are
        'nuclide', 'logarithm', and 'material-union'.
    ptables : bool
        Determine whether probability tables are used.
    run_cmfd : bool
        Indicate if coarse mesh finite difference acceleration is to be used
    seed : int
        Seed for the linear congruential pseudorandom number generator
    survival_biasing : bool
        Indicate whether survival biasing is to be used
    weight : float
        Weight cutoff below which particle undergo Russian roulette
    weight_avg : float
        Weight assigned to particles that are not killed after Russian roulette
    entropy_dimension : tuple or list
        Number of Shannon entropy mesh cells in the x, y, and z directions,
        respectively
    entropy_lower_left : tuple or list
        Coordinates of the lower-left point of the Shannon entropy mesh
    entropy_upper_right : tuple or list
        Coordinates of the upper-right point of the Shannon entropy mesh
    trigger_active : bool
        Indicate whether tally triggers are used
    trigger_max_batches : int
        Maximum number of batches simulated. If this is set, the number of
        batches specified via ``batches`` is interpreted as the minimum number
        of batches
    trigger_batch_interval : int
        Number of batches in between convergence checks
    no_reduce : bool
        Indicate that all user-defined and global tallies should not be reduced
        across processes in a parallel calculation.
    threads : int
        Number of OpenMP threads
    trace : tuple or list
        Show detailed information about a single particle, indicated by three
        integers: the batch number, generation number, and particle number
    track : tuple or list
        Specify particles for which track files should be written. Each particle
        is identified by a triplet with the batch number, generation number, and
        particle number.
    ufs_dimension : tuple or list
        Number of uniform fission site (UFS) mesh cells in the x, y, and z
        directions, respectively
    ufs_lower_left : tuple or list
        Coordinates of the lower-left point of the UFS mesh
    ufs_upper_right : tuple or list
        Coordinates of the upper-right point of the UFS mesh

    """

    def __init__(self):
        # Eigenvalue subelement
        self._batches = None
        self._generations_per_batch = None
        self._inactive = None
        self._particles = None
        self._keff_trigger = None

        # Source subelement
        self._source_subelement = None
        self._source_file = None
        self._source_space_type = None
        self._source_space_params = None
        self._source_angle_type = None
        self._source_angle_params = None
        self._source_energy_type = None
        self._source_energy_params = None

        self._confidence_intervals = None
        self._cross_sections = None
        self._energy_grid = None
        self._ptables = None
        self._run_cmfd = None
        self._seed = None
        self._survival_biasing = None

        # Entropy subelement
        self._entropy_dimension = None
        self._entropy_lower_left = None
        self._entropy_upper_right = None

        # Trigger subelement
        self._trigger_subelement = None
        self._trigger_active = None
        self._trigger_max_batches = None
        self._trigger_batch_interval = None

        self._output = None
        self._output_path = None

        # Statepoint subelement
        self._statepoint_batches = None
        self._statepoint_interval = None
        self._sourcepoint_batches = None
        self._sourcepoint_interval = None
        self._sourcepoint_separate = None
        self._sourcepoint_write = None
        self._sourcepoint_overwrite = None

        self._threads = None
        self._no_reduce = None

        self._verbosity = None

        self._trace = None
        self._track = None

        # Cutoff subelement
        self._weight = None
        self._weight_avg = None

        # Uniform fission source subelement
        self._ufs_dimension = 1
        self._ufs_lower_left = None
        self._ufs_upper_right = None

        # Domain decomposition subelement
        self._dd_mesh_dimension = None
        self._dd_mesh_lower_left = None
        self._dd_mesh_upper_right = None
        self._dd_nodemap = None
        self._dd_allow_leakage = False
        self._dd_count_interactions = False

        self._settings_file = ET.Element("settings")
        self._eigenvalue_subelement = None
        self._source_element = None

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
    def source_file(self):
        return self._source_file

    @property
    def source_space_type(self):
        return self._source_space_type

    @property
    def source_space_params(self):
        return self._source_space_params

    @property
    def source_angle_type(self):
        return self._source_angle_type

    @property
    def source_angle_params(self):
        return self._source_angle_params

    @property
    def source_energy_type(self):
        return self._source_energy_type

    @property
    def source_energy_params(self):
        return self._source_energy_params

    @property
    def confidence_intervals(self):
        return self._confidence_intervals

    @property
    def cross_sections(self):
        return self._cross_sections

    @property
    def energy_grid(self):
        return self._energy_grid

    @property
    def ptables(self):
        return self._ptables

    @property
    def run_cmfd(self):
        return self._run_cmfd

    @property
    def seed(self):
        return self._seed

    @property
    def survival_biasing(self):
        return self._survival_biasing

    @property
    def entropy_dimension(self):
        return self._entropy_dimension

    @property
    def entropy_lower_left(self):
        return self._entropy_lower_left

    @property
    def entropy_upper_right(self):
        return self._entropy_upper_right

    @property
    def trigger_active(self):
        return self._trigger_active

    @property
    def trigger_max_batches(self):
        return self._trigger_max_batches

    @property
    def trigger_batch_interval(self):
        return self._batch_interval

    @property
    def output(self):
        return self._output

    @property
    def output_path(self):
        return self._output_path

    @property
    def statepoint_batches(self):
        return self._statepoint_batches

    @property
    def statepoint_interval(self):
        return self._statepoint_interval

    @property
    def sourcepoint_batches(self):
        return self._sourcepoint_interval

    @property
    def sourcepoint_interval(self):
        return self._sourcepoint_interval

    @property
    def sourcepoint_separate(self):
        return self._sourcepoint_separate

    @property
    def sourcepoint_write(self):
        return self._sourcepoint_write

    @property
    def sourcepoint_overwrite(self):
        return self._sourcepoint_overwrite

    @property
    def threads(self):
        return self._threads

    @property
    def no_reduce(self):
        return self._no_reduce

    @property
    def verbosity(self):
        return self._verbosity

    @property
    def trace(self):
        return self._trace

    @property
    def track(self):
        return self._track

    @property
    def weight(self):
        return self._weight

    @property
    def weight_avg(self):
        return self._weight_avg

    @property
    def ufs_dimension(self):
        return self._ufs_dimension

    @property
    def ufs_lower_left(self):
        return self._ufs_lower_left

    @property
    def ufs_upper_right(self):
        return self._ufs_upper_right

    @property
    def dd_mesh_dimension(self):
        return self._dd_mesh_dimension

    @property
    def dd_mesh_lower_left(self):
        return self._dd_mesh_lower_left

    @property
    def dd_mesh_upper_right(self):
        return self._dd_mesh_upper_right

    @property
    def dd_nodemap(self):
        return self._dd_nodemap

    @property
    def dd_allow_leakage(self):
        return self._dd_allow_leakage

    @property
    def dd_count_interactions(self):
        return self._dd_count_interactions

    @batches.setter
    def batches(self, batches):
        check_type('batches', batches, Integral)
        check_greater_than('batches', batches, 0)
        self._batches = batches

    @generations_per_batch.setter
    def generations_per_batch(self, generations_per_batch):
        check_type('generations per patch', generations_per_batch, Integral)
        check_greater_than('generations per batch', generations_per_batch, 0)
        self._generations_per_batch = generations_per_batch

    @inactive.setter
    def inactive(self, inactive):
        check_type('inactive batches', inactive, Integral)
        check_greater_than('inactive batches', inactive, 0, True)
        self._inactive = inactive

    @particles.setter
    def particles(self, particles):
        check_type('particles', particles, Integral)
        check_greater_than('particles', particles, 0)
        self._particles = particles

    @keff_trigger.setter
    def keff_trigger(self, keff_trigger):
        if not isinstance(keff_trigger, dict):
            msg = 'Unable to set a trigger on keff from {0} which ' \
                  'is not a Python dictionary'.format(keff_trigger)
            raise ValueError(msg)

        elif 'type' not in keff_trigger:
            msg = 'Unable to set a trigger on keff from {0} which ' \
                  'does not have a "type" key'.format(keff_trigger)
            raise ValueError(msg)

        elif keff_trigger['type'] not in ['variance', 'std_dev', 'rel_err']:
            msg = 'Unable to set a trigger on keff with ' \
                  'type {0}'.format(keff_trigger['type'])
            raise ValueError(msg)

        elif 'threshold' not in keff_trigger:
            msg = 'Unable to set a trigger on keff from {0} which ' \
                  'does not have a "threshold" key'.format(keff_trigger)
            raise ValueError(msg)

        elif not isinstance(keff_trigger['threshold'], Real):
            msg = 'Unable to set a trigger on keff with ' \
                  'threshold {0}'.format(keff_trigger['threshold'])
            raise ValueError(msg)

        self._keff_trigger = keff_trigger

    @source_file.setter
    def source_file(self, source_file):
        check_type('source file', source_file, basestring)
        self._source_file = source_file

    def set_source_space(self, stype, params):
        """Defined the spatial bounds of the external/starting source.

        Parameters
        ----------
        stype : str
            The type of spatial distribution. Valid options are "box",
            "fission", and "point". A "box" spatial distribution has coordinates
            sampled uniformly in a parallelepiped. A "fission" spatial
            distribution samples locations from a "box" distribution but only
            locations in fissionable materials are accepted. A "point" spatial
            distribution has coordinates specified by a triplet.
        params : Iterable of float
            For a "box" or "fission" spatial distribution, ``params`` should be
            given as six real numbers, the first three of which specify the
            lower-left corner of a parallelepiped and the last three of which
            specify the upper-right corner. Source sites are sampled uniformly
            through that parallelepiped.

            For a "point" spatial distribution, ``params`` should be given as
            three real numbers which specify the (x,y,z) location of an
            isotropic point source

        """

        check_type('source space type', stype, basestring)
        check_value('source space type', stype, ['box', 'fission', 'point'])
        check_type('source space parameters', params, Iterable, Real)
        if stype in ['box', 'fission']:
            check_length('source space parameters for a '
                         'box/fission distribution', params, 6)
        elif stype == 'point':
            check_length('source space parameters for a point source',
                         params, 3)

        self._source_space_type = stype
        self._source_space_params = params

    def set_source_angle(self, stype, params=[]):
        """Defined the angular distribution of the external/starting source.

        Parameters
        ----------
        stype : str
            The type of angular distribution. Valid options are "isotropic" and
            "monodirectional". The angle of the particle emitted from a source
            site is isotropic if the "isotropic" option is given. The angle of
            the particle emitted from a source site is the direction specified
            in ``params`` if the "monodirectional" option is given.
        params : Iterable of float
            For an "isotropic" angular distribution, ``params`` should not
            be specified.

            For a "monodirectional" angular distribution, ``params`` should
            be given as three floats which specify the angular cosines
            with respect to each axis.

        """

        check_type('source angle type', stype, basestring)
        check_value('source angle type', stype,
                    ['isotropic', 'monodirectional'])
        check_type('source angle parameters', params, Iterable, Real)
        if stype == 'isotropic' and params is not None:
            msg = 'Unable to set source angle parameters since they are not ' \
                  'it is not supported for isotropic type sources'
            raise ValueError(msg)
        elif stype == 'monodirectional':
            check_length('source angle parameters for a monodirectional '
                         'source', params, 3)

        self._source_angle_type = stype
        self._source_angle_params = params

    def set_source_energy(self, stype, params=[]):
        """Defined the energy distribution of the external/starting source.

        Parameters
        ----------
        stype : str
            The type of energy distribution. Valid options are "monoenergetic",
            "watt", and "maxwell". The "monoenergetic" option produces source
            sites at a single energy. The "watt" option produces source sites
            whose energy is sampled from a Watt fission spectrum. The "maxwell"
            option produce source sites whose energy is sampled from a Maxwell
            fission spectrum.
        params : Iterable of float
            For a "monoenergetic" energy distribution, ``params`` should be
            given as the energy in MeV of the source sites.

            For a "watt" energy distribution, ``params`` should be given as two
            real numbers :math:`a` and :math:`b` that parameterize the
            distribution :math:`p(E) dE = c e^{-E/a} \sinh \sqrt{b \, E} dE`.

            For a "maxwell" energy distribution, ``params`` should be given as
            one real number :math:`a` that parameterizes the distribution
            :math:`p(E) dE = c E e^{-E/a} dE`.

        """

        check_type('source energy type', stype, basestring)
        check_value('source energy type', stype,
                    ['monoenergetic', 'watt', 'maxwell'])
        check_type('source energy parameters', params, Iterable, Real)
        if stype in ['monoenergetic', 'maxwell']:
            check_length('source energy parameters for a monoenergetic '
                         'or Maxwell source', params, 1)
        elif stype == 'watt':
            check_length('source energy parameters for a Watt source',
                         params, 2)

        self._source_energy_type = stype
        self._source_energy_params = params

    @output.setter
    def output(self, output):
        if not isinstance(output, dict):
            msg = 'Unable to set output to {0} which is not a Python ' \
                  'dictionary of string keys and boolean values'.format(output)
            raise ValueError(msg)

        for element in output:
            keys = ['summary', 'cross_sections', 'tallies', 'distribmats']
            if element not in keys:
                msg = 'Unable to set output to {0} which is unsupported by ' \
                      'OpenMC'.format(element)
                raise ValueError(msg)

            if not isinstance(output[element], (bool, np.bool)):
                msg = 'Unable to set output for {0} to a non-boolean ' \
                      'value {1}'.format(element, output[element])
                raise ValueError(msg)

        self._output = output

    @output_path.setter
    def output_path(self, output_path):
        check_type('output path', output_path, basestring)
        self._output_path = output_path

    @verbosity.setter
    def verbosity(self, verbosity):
        check_type('verbosity', verbosity, Integral)
        check_greater_than('verbosity', verbosity, 1, True)
        check_less_than('verbosity', verbosity, 10, True)
        self._verbosity = verbosity

    @statepoint_batches.setter
    def statepoint_batches(self, batches):
        check_type('statepoint batches', batches, Iterable, Integral)
        for batch in batches:
            check_greater_than('statepoint batch', batch, 0)
        self._statepoint_batches = batches

    @statepoint_interval.setter
    def statepoint_interval(self, interval):
        check_type('statepoint interval', interval, Integral)
        self._statepoint_interval = interval

    @sourcepoint_batches.setter
    def sourcepoint_batches(self, batches):
        check_type('sourcepoint batches', batches, Iterable, Integral)
        for batch in batches:
            check_greater_than('sourcepoint batch', batch, 0)
        self._sourcepoint_batches = batches

    @sourcepoint_interval.setter
    def sourcepoint_interval(self, interval):
        check_type('sourcepoint interval', interval, Integral)
        self._sourcepoint_interval = interval

    @sourcepoint_separate.setter
    def sourcepoint_separate(self, source_separate):
        check_type('sourcepoint separate', source_separate, bool)
        self._sourcepoint_separate = source_separate

    @sourcepoint_write.setter
    def sourcepoint_write(self, source_write):
        check_type('sourcepoint write', source_write, bool)
        self._sourcepoint_write = source_write

    @sourcepoint_overwrite.setter
    def sourcepoint_overwrite(self, source_overwrite):
        check_type('sourcepoint overwrite', source_overwrite, bool)
        self._sourcepoint_overwrite = source_overwrite

    @confidence_intervals.setter
    def confidence_intervals(self, confidence_intervals):
        check_type('confidence interval', confidence_intervals, bool)
        self._confidence_intervals = confidence_intervals

    @cross_sections.setter
    def cross_sections(self, cross_sections):
        check_type('cross sections', cross_sections, basestring)
        self._cross_sections = cross_sections

    @energy_grid.setter
    def energy_grid(self, energy_grid):
        check_value('energy grid', energy_grid,
                    ['nuclide', 'logarithm', 'material-union'])
        self._energy_grid = energy_grid

    @ptables.setter
    def ptables(self, ptables):
        check_type('probability tables', ptables, bool)
        self._ptables = ptables

    @run_cmfd.setter
    def run_cmfd(self, run_cmfd):
        check_type('run_cmfd', run_cmfd, bool)
        self._run_cmfd = run_cmfd

    @seed.setter
    def seed(self, seed):
        check_type('random number generator seed', seed, Integral)
        check_greater_than('random number generator seed', seed, 0)
        self._seed = seed

    @survival_biasing.setter
    def survival_biasing(self, survival_biasing):
        check_type('survival biasing', survival_biasing, bool)
        self._survival_biasing = survival_biasing

    @weight.setter
    def weight(self, weight):
        check_type('weight cutoff', weight, Real)
        check_greater_than('weight cutoff', weight, 0.0)
        self._weight = weight

    @weight_avg.setter
    def weight_avg(self, weight_avg):
        check_type('average survival weight', weight_avg, Real)
        check_greater_than('average survival weight', weight_avg, 0.0)
        self._weight_avg = weight_avg

    @entropy_dimension.setter
    def entropy_dimension(self, dimension):
        check_type('entropy mesh dimension', dimension, Iterable, Integral)
        check_length('entropy mesh dimension', dimension, 3)
        self._entropy_dimension = dimension

    @entropy_lower_left.setter
    def entropy_lower_left(self, lower_left):
        check_type('entropy mesh lower left corner', lower_left,
                   Iterable, Real)
        check_length('entropy mesh lower left corner', lower_left, 3)
        self._entropy_lower_left = lower_left

    @entropy_upper_right.setter
    def entropy_upper_right(self, upper_right):
        check_type('entropy mesh upper right corner', upper_right,
                   Iterable, Real)
        check_length('entropy mesh upper right corner', upper_right, 3)
        self._entropy_upper_right = upper_right

    @trigger_active.setter
    def trigger_active(self, trigger_active):
        check_type('trigger active', trigger_active, bool)
        self._trigger_active = trigger_active

    @trigger_max_batches.setter
    def trigger_max_batches(self, trigger_max_batches):
        check_type('trigger maximum batches', trigger_max_batches, Integral)
        check_greater_than('trigger maximum batches', trigger_max_batches, 0)
        self._trigger_max_batches = trigger_max_batches

    @trigger_batch_interval.setter
    def trigger_batch_interval(self, trigger_batch_interval):
        check_type('trigger batch interval', trigger_batch_interval, Integral)
        check_greater_than('trigger batch interval', trigger_batch_interval, 0)
        self._trigger_batch_interval = trigger_batch_interval

    @no_reduce.setter
    def no_reduce(self, no_reduce):
        check_type('no reduction option', no_reduce, bool)
        self._no_reduce = no_reduce

    @threads.setter
    def threads(self, threads):
        check_type('number of threads', threads, Integral)
        check_greater_than('number of threads', threads, 0)
        self._threads = threads

    @trace.setter
    def trace(self, trace):
        check_type('trace', trace, Iterable, Integral)
        check_length('trace', trace, 3)
        check_greater_than('trace batch', trace[0], 0)
        check_greater_than('trace generation', trace[1], 0)
        check_greater_than('trace particle', trace[2], 0)
        self._trace = trace

    @track.setter
    def track(self, track):
        check_type('track', track, Iterable, Integral)
        if len(track) % 3 != 0:
            msg = 'Unable to set the track to {0} since its length is ' \
                  'not a multiple of 3'.format(track)
            raise ValueError(msg)
        for t in zip(track[::3], track[1::3], track[2::3]):
            check_greater_than('track batch', t[0], 0)
            check_greater_than('track generation', t[0], 0)
            check_greater_than('track particle', t[0], 0)
        self._track = track

    @ufs_dimension.setter
    def ufs_dimension(self, dimension):
        check_type('UFS mesh dimension', dimension, Iterable, Integral)
        check_length('UFS mesh dimension', dimension, 3)
        for dim in dimension:
            check_greater_than('UFS mesh dimension', dim, 1, True)
        self._ufs_dimension = dimension

    @ufs_lower_left.setter
    def ufs_lower_left(self, lower_left):
        check_type('UFS mesh lower left corner', lower_left, Iterable, Real)
        check_length('UFS mesh lower left corner', lower_left, 3)
        self._ufs_lower_left = lower_left

    @ufs_upper_right.setter
    def ufs_upper_right(self, upper_right):
        check_type('UFS mesh upper right corner', upper_right, Iterable, Real)
        check_length('UFS mesh upper right corner', upper_right, 3)
        self._ufs_upper_right = upper_right

    @dd_mesh_dimension.setter
    def dd_mesh_dimension(self, dimension):
        # TODO: remove this when domain decomposition is merged
        warnings.warn('This feature is not yet implemented in a release '
                      'version of openmc')

        check_type('DD mesh dimension', dimension, Iterable, Integral)
        check_length('DD mesh dimension', dimension, 3)

        self._dd_mesh_dimension = dimension

    @dd_mesh_lower_left.setter
    def dd_mesh_lower_left(self, lower_left):
        # TODO: remove this when domain decomposition is merged
        warnings.warn('This feature is not yet implemented in a release '
                      'version of openmc')

        check_type('DD mesh lower left corner', lower_left, Iterable, Real)
        check_length('DD mesh lower left corner', lower_left, 3)

        self._dd_mesh_lower_left = lower_left

    @dd_mesh_upper_right.setter
    def dd_mesh_upper_right(self, upper_right):
        # TODO: remove this when domain decomposition is merged
        warnings.warn('This feature is not yet implemented in a release '
                      'version of openmc')

        check_type('DD mesh upper right corner', upper_right, Iterable, Real)
        check_length('DD mesh upper right corner', upper_right, 3)

        self._dd_mesh_upper_right = upper_right

    @dd_nodemap.setter
    def dd_nodemap(self, nodemap):
        # TODO: remove this when domain decomposition is merged
        warnings.warn('This feature is not yet implemented in a release '
                      'version of openmc')

        check_type('DD nodemap', nodemap, Iterable)

        nodemap = np.array(nodemap).flatten()

        if self._dd_mesh_dimension is None:
            msg = 'Must set DD mesh dimension before setting the nodemap'
            raise ValueError(msg)
        else:
            len_nodemap = np.prod(self._dd_mesh_dimension)

        if len(nodemap) < len_nodemap or len(nodemap) > len_nodemap:
            msg = 'Unable to set DD nodemap with length {0} which ' \
                  'does not have the same dimensionality as the domain ' \
                  'mesh'.format(len(nodemap))
            raise ValueError(msg)

        self._dd_nodemap = nodemap

    @dd_allow_leakage.setter
    def dd_allow_leakage(self, allow):

        # TODO: remove this when domain decomposition is merged
        warnings.warn('This feature is not yet implemented in a release '
                      'version of openmc')

        check_type('DD allow leakage', allow, bool)

        self._dd_allow_leakage = allow

    @dd_count_interactions.setter
    def dd_count_interactions(self, interactions):

        # TODO: remove this when domain decomposition is merged
        warnings.warn('This feature is not yet implemented in a release '
                      'version of openmc')

        check_type('DD count interactions', interactions, bool)

        self._dd_count_interactions = interactions

    def _create_eigenvalue_subelement(self):
        self._create_particles_subelement()
        self._create_batches_subelement()
        self._create_inactive_subelement()
        self._create_generations_per_batch_subelement()
        self._create_keff_trigger_subelement()

    def _create_batches_subelement(self):
        if self._batches is not None:
            if self._eigenvalue_subelement is None:
                self._eigenvalue_subelement = ET.SubElement(self._settings_file,
                                                         "eigenvalue")

            element = ET.SubElement(self._eigenvalue_subelement, "batches")
            element.text = str(self._batches)

    def _create_generations_per_batch_subelement(self):
        if self._generations_per_batch is not None:
            if self._eigenvalue_subelement is None:
                self._eigenvalue_subelement = ET.SubElement(self._settings_file,
                                                         "eigenvalue")

            element = ET.SubElement(self._eigenvalue_subelement,
                                    "generations_per_batch")
            element.text = str(self._generations_per_batch)

    def _create_inactive_subelement(self):
        if self._inactive is not None:
            if self._eigenvalue_subelement is None:
                self._eigenvalue_subelement = ET.SubElement(self._settings_file,
                                                         "eigenvalue")

            element = ET.SubElement(self._eigenvalue_subelement, "inactive")
            element.text = str(self._inactive)

    def _create_particles_subelement(self):
        if self._particles is not None:
            if self._eigenvalue_subelement is None:
                self._eigenvalue_subelement = ET.SubElement(self._settings_file,
                                                         "eigenvalue")

            element = ET.SubElement(self._eigenvalue_subelement, "particles")
            element.text = str(self._particles)

    def _create_keff_trigger_subelement(self):
        if self._keff_trigger is not None:
            if self._eigenvalue_subelement is None:
                self._eigenvalue_subelement = ET.SubElement(self._settings_file,
                                                         "eigenvalue")

            element = ET.SubElement(self._eigenvalue_subelement, "keff_trigger")

            for key in self._keff_trigger:
                subelement = ET.SubElement(element, key)
                subelement.text = str(self._keff_trigger[key]).lower()

    def _create_source_subelement(self):
        self._create_source_space_subelement()
        self._create_source_energy_subelement()
        self._create_source_angle_subelement()

    def _create_source_space_subelement(self):
        if self._source_space_params is not None:
            if self._source_subelement is None:
                self._source_subelement = ET.SubElement(self._settings_file,
                                                        "source")

            element = ET.SubElement(self._source_subelement, "space")
            element.set("type", self._source_space_type)

            subelement = ET.SubElement(element, "parameters")
            subelement.text = ' '.join(map(str, self._source_space_params))

    def _create_source_angle_subelement(self):
        if self._source_angle_params is not None:
            if self._source_subelement is None:
                self._source_subelement = ET.SubElement(self._settings_file,
                                                        "source")

            element = ET.SubElement(self._source_subelement, "angle")
            element.set("type", self._source_angle_type)

            subelement = ET.SubElement(element, "parameters")
            subelement.text = ' '.join(map(str, self._source_angle_params))

    def _create_source_energy_subelement(self):
        if self._source_energy_params is not None:
            if self._source_subelement is None:
                self._source_subelement = ET.SubElement(self._settings_file,
                                                        "source")

            element = ET.SubElement(self._source_subelement, "energy")
            element.set("type", self._source_energy_type)

            subelement = ET.SubElement(element, "parameters")
            subelement.text = ' '.join(map(str, self._source_energy_params))

    def _create_output_subelement(self):
        if self._output is not None:
            element = ET.SubElement(self._settings_file, "output")

            for key in self._output:
                subelement = ET.SubElement(element, key)
                subelement.text = str(self._output[key]).lower()

            self._create_output_path_subelement()

    def _create_output_path_subelement(self):
        if self._output_path is not None:
            element = ET.SubElement(self._settings_file, "output_path")
            element.text = self._output_path

    def _create_verbosity_subelement(self):
        if self._verbosity is not None:
            element = ET.SubElement(self._settings_file, "verbosity")
            element.text = str(self._verbosity)

    def _create_statepoint_subelement(self):
        # Batches subelement
        if self._statepoint_batches is not None:
            element = ET.SubElement(self._settings_file, "state_point")
            subelement = ET.SubElement(element, "batches")
            subelement.text = ' '.join(map(str, self._statepoint_batches))

        # Interval subelement
        elif self._statepoint_interval is not None:
            element = ET.SubElement(self._settings_file, "state_point")
            subelement = ET.SubElement(element, "interval")
            subelement.text = str(self._statepoint_interval)

    def _create_sourcepoint_subelement(self):
        # Batches subelement
        if self._sourcepoint_batches is not None:
            element = ET.SubElement(self._settings_file, "source_point")
            subelement = ET.SubElement(element, "batches")
            subelement.text = ' '.join(map(str, self._sourcepoint_batches))

        # Interval subelement
        elif self._sourcepoint_interval is not None:
            element = ET.SubElement(self._settings_file, "source_point")
            subelement = ET.SubElement(element, "interval")
            subelement.text = str(self._sourcepoint_interval)

        # Separate subelement
        if self._sourcepoint_separate is not None:
            subelement = ET.SubElement(element, "separate")
            subelement.text = str(self._sourcepoint_separate).lower()

        # Write subelement
        if self._sourcepoint_write is not None:
            subelement = ET.SubElement(element, "write")
            subelement.text = str(self._sourcepoint_write).lower()

        # Overwrite latest subelement
        if self._sourcepoint_overwrite is not None:
            subelement = ET.SubElement(element, "overwrite_latest")
            subelement.text = str(self._sourcepoint_overwrite).lower()

    def _create_confidence_intervals(self):
        if self._confidence_intervals is not None:
            element = ET.SubElement(self._settings_file, "confidence_intervals")
            element.text = str(self._confidence_intervals).lower()

    def _create_cross_sections_subelement(self):
        if self._cross_sections is not None:
            element = ET.SubElement(self._settings_file, "cross_sections")
            element.text = str(self._cross_sections)

    def _create_energy_grid_subelement(self):
        if self._energy_grid is not None:
            element = ET.SubElement(self._settings_file, "energy_grid")
            element.text = str(self._energy_grid)

    def _create_ptables_subelement(self):
        if self._ptables is not None:
            element = ET.SubElement(self._settings_file, "ptables")
            element.text = str(self._ptables).lower()

    def _create_run_cmfd_subelement(self):
        if self._run_cmfd is not None:
            element = ET.SubElement(self._settings_file, "run_cmfd")
            element.text = str(self._run_cmfd).lower()

    def _create_seed_subelement(self):
        if self._seed is not None:
            element = ET.SubElement(self._settings_file, "seed")
            element.text = str(self._seed)

    def _create_survival_biasing_subelement(self):
        if self._survival_biasing is not None:
            element = ET.SubElement(self._settings_file, "survival_biasing")
            element.text = str(self._survival_biasing).lower()

    def _create_cutoff_subelement(self):
        if self._weight is not None:
            element = ET.SubElement(self._settings_file, "cutoff")

            subelement = ET.SubElement(element, "weight")
            subelement.text = str(self._weight)

            subelement = ET.SubElement(element, "weight_avg")
            subelement.text = str(self._weight_avg)

    def _create_entropy_subelement(self):
        if self._entropy_lower_left is not None and \
           self._entropy_upper_right is not None:

            element = ET.SubElement(self._settings_file, "entropy")

            subelement = ET.SubElement(element, "dimension")
            subelement.text = ' '.join(map(str, self._entropy_dimension))

            subelement = ET.SubElement(element, "lower_left")
            subelement.text = ' '.join(map(str, self._entropy_lower_left))

            subelement = ET.SubElement(element, "upper_right")
            subelement.text = ' '.join(map(str, self._entropy_upper_right))

    def _create_trigger_subelement(self):
        self._create_trigger_active_subelement()
        self._create_trigger_max_batches_subelement()
        self._create_trigger_batch_interval_subelement()

    def _create_trigger_active_subelement(self):
        if self._trigger_active is not None:
            if self._trigger_subelement is None:
                self._trigger_subelement = ET.SubElement(self._settings_file,
                                                      "trigger")

            element = ET.SubElement(self._trigger_subelement, "active")
            element.text = str(self._trigger_active).lower()

    def _create_trigger_max_batches_subelement(self):
        if self._trigger_max_batches is not None:
            if self._trigger_subelement is None:
                self._trigger_subelement = ET.SubElement(self._settings_file,
                                                      "trigger")

            element = ET.SubElement(self._trigger_subelement, "max_batches")
            element.text = str(self._trigger_max_batches)

    def _create_trigger_batch_interval_subelement(self):
        if self._trigger_batch_interval is not None:
            if self._trigger_subelement is None:
                self._trigger_subelement = ET.SubElement(self._settings_file,
                                                      "trigger")

            element = ET.SubElement(self._trigger_subelement, "batch_interval")
            element.text = str(self._trigger_batch_interval)

    def _create_no_reduce_subelement(self):
        if self._no_reduce is not None:
            element = ET.SubElement(self._settings_file, "no_reduce")
            element.text = str(self._no_reduce).lower()

    def _create_threads_subelement(self):
        if self._threads is not None:
            element = ET.SubElement(self._settings_file, "threads")
            element.text = str(self._threads)

    def _create_trace_subelement(self):
        if self._trace is not None:
            element = ET.SubElement(self._settings_file, "trace")
            element.text = ' '.join(map(str, self._trace))

    def _create_track_subelement(self):
        if self._track is not None:
            element = ET.SubElement(self._settings_file, "track")
            element.text = ' '.join(map(str, self._track))

    def _create_ufs_subelement(self):
        if self._ufs_lower_left is not None and \
           self._ufs_upper_right is not None:

            element = ET.SubElement(self._settings_file, "uniform_fs")

            subelement = ET.SubElement(element, "dimension")
            subelement.text = str(self._ufs_dimension)

            subelement = ET.SubElement(element, "lower_left")
            subelement.text = ' '.join(map(str, self._ufs_lower_left))

            subelement = ET.SubElement(element, "upper_right")
            subelement.text = ' '.join(map(str, self._ufs_upper_right))

    def _create_dd_subelement(self):
        if self._dd_mesh_lower_left is not None and \
           self._dd_mesh_upper_right is not None and \
           self._dd_mesh_dimension is not None:

            element = ET.SubElement(self._settings_file, "domain_decomposition")

            subelement = ET.SubElement(element, "mesh")
            subsubelement = ET.SubElement(subelement, "dimension")
            subsubelement.text = ' '.join(map(str, self._dd_mesh_dimension))

            subsubelement = ET.SubElement(subelement, "lower_left")
            subsubelement.text = ' '.join(map(str, self._dd_mesh_lower_left))

            subsubelement = ET.SubElement(subelement, "upper_right")
            subsubelement.text = ' '.join(map(str, self._dd_mesh_upper_right))

            if self._dd_nodemap is not None:
                subelement = ET.SubElement(element, "nodemap")
                subelement.text = ' '.join(map(str, self._dd_nodemap))

            subelement = ET.SubElement(element, "allow_leakage")
            subelement.text = str(self._dd_allow_leakage).lower()

            subelement = ET.SubElement(element, "count_interactions")
            subelement.text = str(self._dd_count_interactions).lower()

    def export_to_xml(self):
        """Create a settings.xml file that can be used for a simulation.

        """

        self._create_eigenvalue_subelement()
        self._create_source_subelement()
        self._create_output_subelement()
        self._create_statepoint_subelement()
        self._create_sourcepoint_subelement()
        self._create_confidence_intervals()
        self._create_cross_sections_subelement()
        self._create_energy_grid_subelement()
        self._create_ptables_subelement()
        self._create_run_cmfd_subelement()
        self._create_seed_subelement()
        self._create_survival_biasing_subelement()
        self._create_cutoff_subelement()
        self._create_entropy_subelement()
        self._create_trigger_subelement()
        self._create_no_reduce_subelement()
        self._create_threads_subelement()
        self._create_verbosity_subelement()
        self._create_trace_subelement()
        self._create_track_subelement()
        self._create_ufs_subelement()
        self._create_dd_subelement()

        # Clean the indentation in the file to be user-readable
        clean_xml_indentation(self._settings_file)

        # Write the XML Tree to the settings.xml file
        tree = ET.ElementTree(self._settings_file)
        tree.write("settings.xml", xml_declaration=True,
                             encoding='utf-8', method="xml")
