import collections
import warnings
from xml.etree import ElementTree as ET

import numpy as np

from openmc.checkvalue import *
from openmc.clean_xml import *


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
    statepoint_batches : tuple or list or ndarray
        List of batches at which to write statepoint files
    statepoint_interval : int
        Number of batches after which a new statepoint file should be written
    sourcepoint_batches : tuple or list or ndarray
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
        if not is_integer(batches):
            msg = 'Unable to set batches to a non-integer ' \
                  'value {0}'.format(batches)
            raise ValueError(msg)

        if batches <= 0:
            msg = 'Unable to set batches to a negative ' \
                  'value {0}'.format(batches)
            raise ValueError(msg)

        self._batches = batches

    @generations_per_batch.setter
    def generations_per_batch(self, generations_per_batch):
        if not is_integer(generations_per_batch):
            msg = 'Unable to set generations per batch to a non-integer ' \
                  'value {0}'.format(generations_per_batch)
            raise ValueError(msg)

        if generations_per_batch <= 0:
            msg = 'Unable to set generations per batch to a negative ' \
                  'value {0}'.format(generations_per_batch)
            raise ValueError(msg)

        self._generations_per_batch = generations_per_batch

    @inactive.setter
    def inactive(self, inactive):
        if not is_integer(inactive):
            msg = 'Unable to set inactive batches to a non-integer ' \
                  'value {0}'.format(inactive)
            raise ValueError(msg)

        if inactive <= 0:
            msg = 'Unable to set inactive batches to a negative ' \
                  'value {0}'.format(inactive)
            raise ValueError(msg)

        self._inactive = inactive

    @particles.setter
    def particles(self, particles):
        if not is_integer(particles):
            msg = 'Unable to set particles to a non-integer ' \
                  'value {0}'.format(particles)
            raise ValueError(msg)

        if particles <= 0:
            msg = 'Unable to set particles to a negative ' \
                  'value {0}'.format(particles)
            raise ValueError(msg)

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

        elif not is_float(keff_trigger['threshold']):
            msg = 'Unable to set a trigger on keff with ' \
                  'threshold {0}'.format(keff_trigger['threshold'])
            raise ValueError(msg)

        self._keff_trigger = keff_trigger

    @source_file.setter
    def source_file(self, source_file):
        if not is_string(source_file):
            msg = 'Unable to set source file to a non-string ' \
                  'value {0}'.format(source_file)
            raise ValueError(msg)

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
        params : tuple or list or ndarray
            For a "box" or "fission" spatial distribution, ``params`` should be
            given as six real numbers, the first three of which specify the
            lower-left corner of a parallelepiped and the last three of which
            specify the upper-right corner. Source sites are sampled uniformly
            through that parallelepiped.

            For a "point" spatial distribution, ``params`` should be given as
            three real numbers which specify the (x,y,z) location of an
            isotropic point source

        """

        if not is_string(stype):
            msg = 'Unable to set source space type to a non-string ' \
                  'value {0}'.format(stype)
            raise ValueError(msg)

        elif stype not in ['box', 'fission', 'point']:
            msg = 'Unable to set source space type to {0} since it is not ' \
                  'box or point'.format(stype)
            raise ValueError(msg)

        elif not isinstance(params, (tuple, list, np.ndarray)):
            msg = 'Unable to set source space parameters to {0} since it is ' \
                  'not a Python tuple, list or NumPy array'.format(params)
            raise ValueError(msg)

        elif stype in ['box', 'fission'] and len(params) != 6:
            msg = 'Unable to set source space parameters for a box/fission ' \
                  'distribution to {0} since it does not contain 6 values'\
                      .format(params)
            raise ValueError(msg)

        elif stype == 'point' and len(params) != 3:
            msg = 'Unable to set source space parameters for a point to {0} ' \
                  'since it does not contain 3 values'.format(params)
            raise ValueError(msg)

        for param in params:
            if not is_integer(param) and not is_float(param):
                msg = 'Unable to set source space parameters to {0} since it ' \
                      'is not an integer or floating point value'.format(param)
                raise ValueError(msg)

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
        params : tuple or list or ndarray
            For an "isotropic" angular distribution, ``params`` should not
            be specified.

            For a "monodirectional" angular distribution, ``params`` should
            be given as three floats which specify the angular cosines
            with respect to each axis.

        """

        if not is_string(stype):
            msg = 'Unable to set source angle type to a non-string ' \
                  'value {0}'.format(stype)
            raise ValueError(msg)

        elif stype not in ['isotropic', 'monodirectional']:
            msg = 'Unable to set source angle type to {0} since it is not ' \
                  'isotropic or monodirectional'.format(stype)
            raise ValueError(msg)

        elif not isinstance(params, (tuple, list, np.ndarray)):
            msg = 'Unable to set source angle parameters to {0} since it is ' \
                  'not a Python list/tuple or NumPy array'.format(params)
            raise ValueError(msg)

        elif stype == 'isotropic' and params is not None:
            msg = 'Unable to set source angle parameters since they are not ' \
                  'it is not supported for isotropic type sources'
            raise ValueError(msg)

        elif stype == 'monodirectional' and len(params) != 3:
            msg = 'Unable to set source angle parameters to {0} ' \
                  'since 3 parameters are required for monodirectional ' \
                  'sources'.format(params)
            raise ValueError(msg)

        for param in params:
            if not is_integer(param) and not is_float(param):
                msg = 'Unable to set source angle parameters to {0} since it ' \
                      'is not an integer or floating point value'.format(param)
                raise ValueError(msg)

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
        params : tuple or list or ndarray
            For a "monoenergetic" energy distribution, ``params`` should be
            given as the energy in MeV of the source sites.

            For a "watt" energy distribution, ``params`` should be given as two
            real numbers :math:`a` and :math:`b` that parameterize the
            distribution :math:`p(E) dE = c e^{-E/a} \sinh \sqrt{b \, E} dE`.

            For a "maxwell" energy distribution, ``params`` should be given as
            one real number :math:`a` that parameterizes the distribution
            :math:`p(E) dE = c E e^{-E/a} dE`.

        """

        if not is_string(stype):
            msg = 'Unable to set source energy type to a non-string ' \
                   'value {0}'.format(stype)
            raise ValueError(msg)

        elif stype not in ['monoenergetic', 'watt', 'maxwell']:
            msg = 'Unable to set source energy type to {0} since it is not ' \
                  'monoenergetic, watt or maxwell'.format(stype)
            raise ValueError(msg)

        elif not isinstance(params, (tuple, list, np.ndarray)):
            msg = 'Unable to set source energy params to {0} since it ' \
                  'is not a Python list/tuple or NumPy array'.format(params)
            raise ValueError(msg)

        elif stype == 'monoenergetic' and not len(params) != 1:
            msg = 'Unable to set source energy params to {0} ' \
                  'since 1 paramater is required for monenergetic ' \
                  'sources'.format(params)
            raise ValueError(msg)

        elif stype == 'watt' and len(params) != 2:
            msg = 'Unable to set source energy params to {0} ' \
                  'since 2 params are required for monoenergetic ' \
                  'sources'.format(params)
            raise ValueError(msg)

        elif stype == 'maxwell' and len(params) != 2:
            msg = 'Unable to set source energy params to {0} since 1 ' \
                  'parameter is required for maxwell sources'.format(params)
            raise ValueError(msg)

        for param in params:
            if not is_integer(param) and not is_float(param):
                msg = 'Unable to set source energy params to {0} ' \
                      'since it is not an integer or floating point ' \
                      'value'.format(param)
                raise ValueError(msg)

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
        if not is_string(output_path):
            msg = 'Unable to set output path to non-string ' \
                  'value {0}'.format(output_path)
            raise ValueError(msg)

        self._output_path = output_path

    @verbosity.setter
    def verbosity(self, verbosity):
        if not is_integer(verbosity):
            msg = 'Unable to set verbosity to non-integer ' \
                  'value {0}'.format(verbosity)
            raise ValueError(msg)

        if verbosity < 1 or verbosity > 10:
            msg = 'Unable to set verbosity to {0} which is not between ' \
                  '1 and 10'.format(verbosity)
            raise ValueError(msg)

        self._verbosity = verbosity

    @statepoint_batches.setter
    def statepoint_batches(self, batches):
        if not isinstance(batches, (tuple, list, np.ndarray)):
            msg = 'Unable to set statepoint batches to {0} which is not a ' \
                  'Python tuple/list or NumPy array'.format(batches)
            raise ValueError(msg)

        for batch in batches:
            if not is_integer(batch):
                msg = 'Unable to set statepoint batches with non-integer ' \
                      'value {0}'.format(batch)
                raise ValueError(msg)

            if batch <= 0:
                msg = 'Unable to set statepoint batches with {0} which is ' \
                      'less than or equal to zero'.format(batch)
                raise ValueError(msg)

        self._statepoint_batches = batches

    @statepoint_interval.setter
    def statepoint_interval(self, interval):
        if not is_integer(interval):
            msg = 'Unable to set statepoint interval to non-integer ' \
                        'value {0}'.format(interval)
            raise ValueError(msg)

        self._statepoint_interval = interval

    @sourcepoint_batches.setter
    def sourcepoint_batches(self, batches):
        if not isinstance(batches, (tuple, list, np.ndarray)):
            msg = 'Unable to set sourcepoint batches to {0} which is ' \
                  'not a Python tuple/list or NumPy array'.format(batches)
            raise ValueError(msg)

        for batch in batches:
            if not is_integer(batch):
                msg = 'Unable to set sourcepoint batches with non-integer ' \
                      'value {0}'.format(batch)
                raise ValueError(msg)

            if batch <= 0:
                msg = 'Unable to set sourcepoint batches with {0} which is ' \
                      'less than or equal to zero'.format(batch)
                raise ValueError(msg)

        self._sourcepoint_batches = batches

    @sourcepoint_interval.setter
    def sourcepoint_interval(self, interval):
        if not is_integer(interval):
            msg = 'Unable to set sourcepoint interval to non-integer ' \
                  'value {0}'.format(interval)
            raise ValueError(msg)

        self._sourcepoint_interval = interval

    @sourcepoint_separate.setter
    def sourcepoint_separate(self, source_separate):
        if not isinstance(source_separate, (bool, np.bool)):
            msg = 'Unable to set sourcepoint separate to non-boolean ' \
                  'value {0}'.format(source_separate)
            raise ValueError(msg)

        self._sourcepoint_separate = source_separate

    @sourcepoint_write.setter
    def sourcepoint_write(self, source_write):
        if not isinstance(source_write, (bool, np.bool)):
            msg = 'Unable to set sourcepoint write to non-boolean ' \
                  'value {0}'.format(source_write)
            raise ValueError(msg)

        self._sourcepoint_write = source_write

    @sourcepoint_overwrite.setter
    def sourcepoint_overwrite(self, source_overwrite):
        if not isinstance(source_overwrite, (bool, np.bool)):
            msg = 'Unable to set sourcepoint overwrite to non-boolean ' \
                  'value {0}'.format(source_overwrite)
            raise ValueError(msg)

        self._sourcepoint_overwrite = source_overwrite

    @confidence_intervals.setter
    def confidence_intervals(self, confidence_intervals):
        if not isinstance(confidence_intervals, (bool, np.bool)):
            msg = 'Unable to set confidence interval to non-boolean ' \
                  'value {0}'.format(confidence_intervals)
            raise ValueError(msg)

        self._confidence_intervals = confidence_intervals

    @cross_sections.setter
    def cross_sections(self, cross_sections):
        if not is_string(cross_sections):
            msg = 'Unable to set cross sections to non-string ' \
                  'value {0}'.format(cross_sections)
            raise ValueError(msg)

        self._cross_sections = cross_sections

    @energy_grid.setter
    def energy_grid(self, energy_grid):
        if energy_grid not in ['nuclide', 'logarithm', 'material-union']:
            msg = 'Unable to set energy grid to {0} which is neither ' \
                  'nuclide, logarithm, nor material-union'.format(energy_grid)
            raise ValueError(msg)

        self._energy_grid = energy_grid

    @ptables.setter
    def ptables(self, ptables):
        if not isinstance(ptables, (bool, np.bool)):
            msg = 'Unable to set ptables to non-boolean ' \
                  'value {0}'.format(ptables)
            raise ValueError(msg)

        self._ptables = ptables

    @run_cmfd.setter
    def run_cmfd(self, run_cmfd):
        if not isinstance(run_cmfd, (bool, np.bool)):
            msg = 'Unable to set run_cmfd to non-boolean ' \
                  'value {0}'.format(run_cmfd)
            raise ValueError(msg)

        self._run_cmfd = run_cmfd

    @seed.setter
    def seed(self, seed):
        if not is_integer(seed):
            msg = 'Unable to set seed to non-integer value {0}'.format(seed)
            raise ValueError(msg)

        elif seed <= 0:
            msg = 'Unable to set seed to non-positive integer {0}'.format(seed)
            raise ValueError(msg)

        self._seed = seed

    @survival_biasing.setter
    def survival_biasing(self, survival_biasing):
        if not isinstance(survival_biasing, (bool, np.bool)):
            msg = 'Unable to set survival biasing to non-boolean ' \
                  'value {0}'.format(survival_biasing)
            raise ValueError(msg)

        self._survival_biasing = survival_biasing

    @weight.setter
    def weight(self, weight):
        if not is_float(weight):
            msg = 'Unable to set weight cutoff to non-floating point ' \
                  'value {0}'.format(weight)
            raise ValueError(msg)

        elif weight < 0.0:
            msg = 'Unable to set weight cutoff to negative ' \
                  'value {0}'.format(weight)
            raise ValueError(msg)

        self._weight = weight

    @weight_avg.setter
    def weight_avg(self, weight_avg):
        if not is_float(weight_avg):
            msg = 'Unable to set weight avg. to non-floating point ' \
                  'value {0}'.format(weight_avg)
            raise ValueError(msg)
        elif weight_avg < 0.0:
            msg = 'Unable to set weight avg. to negative ' \
                  'value {0}'.format(weight_avg)
            raise ValueError(msg)

        self._weight_avg = weight_avg

    @entropy_dimension.setter
    def entropy_dimension(self, dimension):
        if not isinstance(dimension, (tuple, list)):
            msg = 'Unable to set entropy mesh dimension to {0} which is ' \
                  'not a Python tuple or list'.format(dimension)
            raise ValueError(msg)

        elif len(dimension) != 3:
            msg = 'Unable to set entropy mesh dimension to {0} which is ' \
                  'not a set of 3 integer dimensions'.format(dimension)
            raise ValueError(msg)

        for dim in dimension:
            if not is_integer(dim) and not is_float(dim):
                msg = 'Unable to set entropy mesh dimension to a ' \
                      'non-integer or floating point value {0}'.format(dim)
                raise ValueError(msg)

        self._entropy_dimension = dimension

    @entropy_lower_left.setter
    def entropy_lower_left(self, lower_left):
        if not isinstance(lower_left, (tuple, list)):
            msg = 'Unable to set entropy mesh lower left corner to {0} which ' \
                  'is not a Python tuple or list'.format(lower_left)
            raise ValueError(msg)

        elif len(lower_left) != 3:
            msg = 'Unable to set entropy mesh lower left corner to {0} which ' \
                  'is not a 3D point'.format(lower_left)
            raise ValueError(msg)

        for coord in lower_left:
            if not is_integer(coord) and not is_float(coord):
                msg = 'Unable to set entropy mesh lower left corner to a ' \
                      'non-integer or floating point value {0}'.format(coord)
                raise ValueError(msg)

        self._entropy_lower_left = lower_left

    @entropy_upper_right.setter
    def entropy_upper_right(self, upper_right):
        if not isinstance(upper_right, (tuple, list)):
            msg = 'Unable to set entropy mesh upper right corner to {0} ' \
                  'which is not a Python tuple or list'.format(upper_right)
            raise ValueError(msg)

        elif len(upper_right) < 3 or len(upper_right) > 3:
            msg = 'Unable to set entropy mesh upper right corner to {0} ' \
                  'which is not a 3D point'.format(upper_right)
            raise ValueError(msg)

        for coord in upper_right:
            if not is_integer(coord) and not is_float(coord):
                msg = 'Unable to set entropy mesh upper right corner to a ' \
                      'non-integer or floating point value {0}'.format(coord)
                raise ValueError(msg)

        self._entropy_upper_right = upper_right

    @trigger_active.setter
    def trigger_active(self, trigger_active):
        if not isinstance(trigger_active, bool):
            msg = 'Unable to set trigger active to a ' \
                  'non-boolean value {0}'.format(trigger_active)
            raise ValueError(msg)

        self._trigger_active = trigger_active

    @trigger_max_batches.setter
    def trigger_max_batches(self, trigger_max_batches):
        if not is_integer(trigger_max_batches):
            msg = 'Unable to set trigger max batches to a non-integer ' \
                  'value {0}'.format(trigger_max_batches)
            raise ValueError(msg)

        elif trigger_max_batches <= 0:
            msg = 'Unable to set trigger max batches to a non-positive ' \
                  'value {0}'.format(trigger_max_batches)
            raise ValueError(msg)

        self._trigger_max_batches = trigger_max_batches

    @trigger_batch_interval.setter
    def trigger_batch_interval(self, trigger_batch_interval):
        if not is_integer(trigger_batch_interval):
            msg = 'Unable to set trigger batch interval to a non-integer ' \
                  'value {0}'.format(trigger_batch_interval)
            raise ValueError(msg)

        elif trigger_batch_interval <= 0:
            msg = 'Unable to set trigger batch interval to a non-positive ' \
                  'value {0}'.format(trigger_batch_interval)
            raise ValueError(msg)

        self._trigger_batch_interval = trigger_batch_interval

    @no_reduce.setter
    def no_reduce(self, no_reduce):
        if not isinstance(no_reduce, (bool, np.bool)):
            msg = 'Unable to set the no_reduce to a non-boolean ' \
                  'value {0}'.format(no_reduce)
            raise ValueError(msg)

        self._no_reduce = no_reduce

    @threads.setter
    def threads(self, threads):
        if not is_integer(threads):
            msg = 'Unable to set the threads to a non-integer ' \
                  'value {0}'.format(threads)
            raise ValueError(msg)

        elif threads <= 0:
            msg = 'Unable to set the threads to a negative ' \
                  'value {0}'.format(threads)
            raise ValueError(msg)

        self._threads = threads

    @trace.setter
    def trace(self, trace):
        if not isinstance(trace, (list, tuple)):
            msg = 'Unable to set the trace to {0} which is not a Python ' \
                  'tuple or list'.format(trace)
            raise ValueError(msg)

        elif len(trace) != 3:
            msg = 'Unable to set the trace to {0} since it does not contain ' \
                  '3 elements - batch, generation, and particle'.format(trace)
            raise ValueError(msg)

        elif trace[0] < 1:
            msg = 'Unable to set the trace batch to {0} since it must be ' \
                  'greater than or equal to 1'.format(trace[0])
            raise ValueError(msg)

        elif trace[1] < 1:
            msg = 'Unable to set the trace generation to {0} since it ' \
                  'must be greater than or equal to 1'.format(trace[1])
            raise ValueError(msg)

        elif trace[2] < 1:
            msg = 'Unable to set the trace particle to {0} since it ' \
                  'must be greater than or equal to 1'.format(trace[2])
            raise ValueError(msg)

        self._trace = trace

    @track.setter
    def track(self, track):
        if not isinstance(track, (list, tuple)):
            msg = 'Unable to set the track to {0} which is not a Python ' \
                  'tuple or list'.format(track)
            raise ValueError(msg)

        elif len(track) != 3:
            msg = 'Unable to set the track to {0} since it does not contain ' \
                  '3 elements - batch, generation, and particle'.format(track)
            raise ValueError(msg)

        elif track[0] < 1:
            msg = 'Unable to set the track batch to {0} since it must be ' \
                  'greater than or equal to 1'.format(track[0])
            raise ValueError(msg)

        elif track[1] < 1:
            msg = 'Unable to set the track generation to {0} since it must ' \
                  'be greater than or equal to 1'.format(track[1])
            raise ValueError(msg)

        elif track[2] < 1:
            msg = 'Unable to set the track particle to {0} since it must ' \
                  'be greater than or equal to 1'.format(track[2])
            raise ValueError(msg)

        self._track = track

    @ufs_dimension.setter
    def ufs_dimension(self, dimension):
        if not isinstance(dimension, (tuple, list)):
            msg = 'Unable to set UFS mesh dimension to {0} which is ' \
                  'not a Python tuple or list'.format(dimension)
            raise ValueError(msg)

        elif len(dimension) != 3:
            msg = 'Unable to set UFS mesh dimension to {0} which is ' \
                  'not a set of 3 integer dimensions'.format(dimension)
            raise ValueError(msg)

        for dim in dimension:
            if not is_integer(dim):
                msg = 'Unable to set entropy mesh dimension to a ' \
                      'non-integer {0}'.format(dim)
                raise ValueError(msg)
            elif dim < 1:
                msg = 'Unable to set UFS dimension to value {0} which is ' \
                      'less than one'.format(dimension)
                raise ValueError(msg)

        self._ufs_dimension = dimension

    @ufs_lower_left.setter
    def ufs_lower_left(self, lower_left):
        if not isinstance(lower_left, (tuple, list, np.ndarray)):
            msg = 'Unable to set UFS mesh lower left corner to {0} which is ' \
                  'not a Python tuple or list'.format(lower_left)
            raise ValueError(msg)

        elif len(lower_left) != 3:
            msg = 'Unable to set UFS mesh lower left corner to {0} which ' \
                  'is not a 3D point'.format(lower_left)
            raise ValueError(msg)

        self._ufs_lower_left = lower_left

    @ufs_upper_right.setter
    def ufs_upper_right(self, upper_right):
        if not isinstance(upper_right, (tuple, list)):
            msg = 'Unable to set UFs mesh upper right corner to {0} which is ' \
                  'not a Python tuple or list'.format(upper_right)
            raise ValueError(msg)

        if len(upper_right) != 3:
            msg = 'Unable to set UFS mesh upper right corner to {0} which ' \
                  'is not a 3D point'.format(upper_right)
            raise ValueError(msg)

        self._ufs_upper_right = upper_right

    @dd_mesh_dimension.setter
    def dd_mesh_dimension(self, dimension):
        # TODO: remove this when domain decomposition is merged
        warnings.warn('This feature is not yet implemented in a release '
                      'version of openmc')

        if not isinstance(dimension, (tuple, list)):
            msg = 'Unable to set DD mesh upper right corner to {0} which is ' \
                  'not a Python tuple or list'.format(dimension)
            raise ValueError(msg)

        if len(dimension) != 3:
            msg = 'Unable to set DD mesh upper right corner to {0} which ' \
                  'is not a 3D point'.format(dimension)
            raise ValueError(msg)

        self._dd_mesh_dimension = dimension

    @dd_mesh_lower_left.setter
    def dd_mesh_lower_left(self, lower_left):
        # TODO: remove this when domain decomposition is merged
        warnings.warn('This feature is not yet implemented in a release '
                      'version of openmc')

        if not isinstance(lower_left, (tuple, list, np.ndarray)):
            msg = 'Unable to set DD mesh lower left corner to {0} which is ' \
                  'not a Python tuple or list'.format(lower_left)
            raise ValueError(msg)

        elif len(lower_left) < 3 or len(lower_left) > 3:
            msg = 'Unable to set DD mesh lower left corner to {0} which ' \
                  'is not a 3D point'.format(lower_left)
            raise ValueError(msg)

        self._dd_mesh_lower_left = lower_left

    @dd_mesh_upper_right.setter
    def dd_mesh_upper_right(self, upper_right):
        # TODO: remove this when domain decomposition is merged
        warnings.warn('This feature is not yet implemented in a release '
                      'version of openmc')

        if not isinstance(upper_right, tuple) and \
          not isinstance(upper_right, list):
            msg = 'Unable to set DD mesh upper right corner to {0} which is ' \
                  'not a Python tuple or list'.format(upper_right)
            raise ValueError(msg)

        if len(upper_right) < 3 or len(upper_right) > 3:
            msg = 'Unable to set DD mesh upper right corner to {0} which ' \
                  'is not a 3D point'.format(upper_right)
            raise ValueError(msg)

        self._dd_mesh_upper_right = upper_right

    @dd_nodemap.setter
    def dd_nodemap(self, nodemap):
        # TODO: remove this when domain decomposition is merged
        warnings.warn('This feature is not yet implemented in a release '
                      'version of openmc')

        if not isinstance(nodemap, (tuple, list)):
            msg = 'Unable to set DD nodemap {0} which is ' \
                  'not a Python tuple or list'.format(nodemap)
            raise ValueError(msg)

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

        if not isinstance(allow, bool):
            msg = 'Unable to set DD allow_leakage {0} which is ' \
                  'not a Python bool'.format(allow)
            raise ValueError(msg)

        self._dd_allow_leakage = allow

    @dd_count_interactions.setter
    def dd_count_interactions(self, interactions):

        # TODO: remove this when domain decomposition is merged
        warnings.warn('This feature is not yet implemented in a release '
                      'version of openmc')

        if not isinstance(interactions, bool):
            msg = 'Unable to set DD count_interactions {0} which is ' \
                  'not a Python bool'.format(interactions)
            raise ValueError(msg)

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
