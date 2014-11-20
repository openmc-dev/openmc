#!/usr/bin/env python

import warnings

from openmc.checkvalue import *
from openmc.clean_xml import *
from xml.etree import ElementTree as ET
import numpy as np
import multiprocessing


class SettingsFile(object):

    def __init__(self):

        # Eigenvalue subelement
        self._batches = None
        self._generations_per_batch = None
        self._inactive = None
        self._particles = None

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

        self._settings_file = ET.Element("settings")
        self._eigenvalue_element = None
        self._source_element = None


    def set_batches(self, batches):

        if not is_integer(batches):
            msg = 'Unable to set batches to a non-integer ' \
                  'value {0}'.format(batches)
            raise ValueError(msg)

        if batches <= 0:
            msg = 'Unable to set batches to a negative ' \
                  'value {0}'.format(batches)
            raise ValueError(msg)

        self._batches = batches


    def set_generations_per_batch(self, generations_per_batch):

        if not is_integer(generations_per_batch):
            msg = 'Unable to set generations per batch to a non-integer ' \
                  'value {0}'.format(generations_per_batch)
            raise ValueError(msg)

        if generations_per_batch <= 0:
            msg = 'Unable to set generations per batch to a negative ' \
                  'value {0}'.format(generations_per_batch)
            raise ValueError(msg)

        self._generations_per_batch = generations_per_batch


    def set_inactive(self, inactive):

        if not is_integer(inactive):
            msg = 'Unable to set inactive batches to a non-integer ' \
                  'value {0}'.format(inactive)
            raise ValueError(msg)

        if inactive <= 0:
            msg = 'Unable to set inactive batches to a negative ' \
                  'value {0}'.format(inactive)
            raise ValueError(msg)

        self._inactive = inactive


    def set_particles(self, particles):

        if not is_integer(particles):
            msg = 'Unable to set particles to a non-integer ' \
                  'value {0}'.format(particles)
            raise ValueError(msg)

        if particles <= 0:
            msg = 'Unable to set particles to a negative ' \
                  'value {0}'.format(particles)
            raise ValueError(msg)

        self._particles = particles


    def set_source_file(self, source_file):

        if not is_string(source_file):
            msg = 'Unable to set source file to a non-string ' \
                  'value {0}'.format(source_file)
            raise ValueError(msg)

        self._source_file = source_file


    def set_source_space(self, type, params):

        if not is_string(type):
            msg = 'Unable to set source space type to a non-string ' \
                  'value {0}'.format(type)
            raise ValueError(msg)

        elif not type in ['box', 'point']:
            msg = 'Unable to set source space type to {0} since it is not ' \
                  'box or point'.format(type)
            raise ValueError(msg)

        elif not isinstance(params, (tuple, list, np.ndarray)):
            msg = 'Unable to set source space parameters to {0} since it is ' \
                  'not a Python tuple, list or NumPy array'.format(params)
            raise ValueError(msg)

        elif len(params) != 6:
            msg = 'Unable to set source space parameters to {0} since it ' \
                  'does not contain 6 values'.format(params)
            raise ValueError(msg)

        for param in params:

            if not is_integer(param) and not is_float(param):
                msg = 'Unable to set source space parameters to {0} since it ' \
                      'is not an integer or floating point value'.format(param)
                raise ValueError(msg)

        self._source_space_type = type
        self._source_space_params = params


    def set_source_angle(self, type, params=[]):

        if not is_string(type):
            msg = 'Unable to set source angle type to a non-string ' \
                  'value {0}'.format(type)
            raise ValueError(msg)

        elif not type in ['isotropic', 'monodirectional']:
            msg = 'Unable to set source angle type to {0} since it is not ' \
                  'isotropic or monodirectional'.format(type)
            raise ValueError(msg)

        elif not isinstance(params, (tuple, list, np.ndarray)):
            msg = 'Unable to set source angle parameters to {0} since it is ' \
                  'not a Python list/tuple or NumPy array'.format(params)
            raise ValueError(msg)

        elif type is 'isotropic' and not params is None:
            msg = 'Unable to set source angle parameters since they are not ' \
                  'it is not supported for isotropic type sources'
            raise ValueError(msg)

        elif type is 'monodirectional' and len(params) != 3:
            msg = 'Unable to set source angle parameters to {0} ' \
                  'since 3 parameters are required for monodirectional ' \
                  'sources'.format(params)
            raise ValueError(msg)

        for param in params:

            if not is_integer(param) and not is_float(param):
                msg = 'Unable to set source angle parameters to {0} since it ' \
                      'is not an integer or floating point value'.format(param)
                raise ValueError(msg)

        self._source_angle_type = type
        self._source_angle_params = params


    def set_source_energy(self, type, params=[]):

        if not is_string(type):
            msg = 'Unable to set source energy type to a non-string ' \
                   'value {0}'.format(type)
            raise ValueError(msg)

        elif not type in ['monoenergetic', 'watt', 'maxwell']:
            msg = 'Unable to set source energy type to {0} since it is not ' \
                  'monoenergetic, watt or maxwell'.format(type)
            raise ValueError(msg)

        elif not isinstance(params, (tuple, list, np.ndarray)):
            msg = 'Unable to set source energy parameters to {0} since it ' \
                  'is not a Python list/tuple or NumPy array'.format(params)
            raise ValueError(msg)

        elif type is 'monoenergetic' and not len(params) != 1:
            msg = 'Unable to set source energy parameters to {0} ' \
                  'since 1 paramater is required for monenergetic ' \
                  'sources'.format(params)
            raise ValueError(msg)

        elif type is 'watt' and len(params) != 2:
            msg = 'Unable to set source energy parameters to {0} ' \
                  'since 2 parameters are required for monoenergetic ' \
                  'sources'.format(params)
            raise ValueError(msg)

        elif type is 'maxwell' and len(params) != 2:
            msg = 'Unable to set source energy parameters to {0} since 1 ' \
                  'parameter is required for maxwell sources'.format(params)
            raise ValueError(msg)

        for param in params:

            if not is_integer(param) and not is_float(param):
                msg = 'Unable to set source energy parameters to {0} ' \
                      'since it is not an integer or floating point ' \
                      'value'.format(param)
                raise ValueError(msg)

        self._source_energy_type = type
        self._source_energy_params = params


    def set_output(self, output):

        if not isinstance(output, dict):
            msg = 'Unable to set output to {0} which is not a Python ' \
                  'dictionary of string keys and boolean values'.format(output)
            raise ValueError(msg)


        for element in output:

            keys = ['summary', 'cross_sections', 'tallies', 'distribmats']
            if not element in keys:
                msg = 'Unable to set output to {0} which is unsupported by ' \
                      'OpenMC'.format(element)
                raise ValueError(msg)

            if not isinstance(output[element], (bool, np.bool)):
                msg = 'Unable to set output for {0} to a non-boolean ' \
                      'value {1}'.format(element, output[element])
                raise ValueError(msg)

        self._output = output


    def set_output_path(self, output_path):

        if not is_string(output_path):
            msg = 'Unable to set output path to non-string ' \
                  'value {0}'.format(output_path)
            raise ValueError(msg)

        self._output_path = output_path


    def set_verbosity(self, verbosity):

        if not is_integer(verbosity):
            msg = 'Unable to set verbosity to non-integer ' \
                  'value {0}'.format(verbosity)
            raise ValueError(msg)

        if verbosity < 1 or verbosity > 10:
            msg = 'Unable to set verbosity to {0} which is not between ' \
                  '1 and 10'.format(verbosity)
            raise ValueError(msg)

        self._verbosity = verbosity


    def set_statepoint_batches(self, batches):

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


    def set_statepoint_interval(self, interval):

        if not is_integer(interval):
            msg = 'Unable to set statepoint interval to non-integer ' \
                        'value {0}'.format(interval)
            raise ValueError(msg)

        self._statepoint_interval = interval


    def set_sourcepoint_batches(self, batches):

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


    def set_sourcepoint_interval(self, interval):

        if not is_integer(interval):
            msg = 'Unable to set sourcepoint interval to non-integer ' \
                  'value {0}'.format(interval)
            raise ValueError(msg)

        self._sourcepoint_interval = interval


    def set_sourcepoint_separate(self, source_separate):

        if not isinstance(source_separate, (bool, np.bool)):
            msg = 'Unable to set sourcepoint separate to non-boolean ' \
                  'value {0}'.format(source_separate)
            raise ValueError(msg)

        self._sourcepoint_separate = source_separate


    def set_sourcepoint_write(self, source_write):

        if not isinstance(source_write, (bool, np.bool)):
            msg = 'Unable to set sourcepoint write to non-boolean ' \
                  'value {0}'.format(source_write)
            raise ValueError(msg)

        self._sourcepoint_write = source_write


    def set_sourcepoint_overwrite(self, source_overwrite):

        if not isinstance(source_overwrite, (bool, np.bool)):
            msg = 'Unable to set sourcepoint overwrite to non-boolean ' \
                  'value {0}'.format(source_overwrite)
            raise ValueError(msg)

        self._sourcepoint_overwrite = source_overwrite


    def set_confidence_intervals(self, confidence_intervals):

        if not isinstance(confidence_intervals, (bool, np.bool)):
            msg = 'Unable to set confidence interval to non-boolean ' \
                  'value {0}'.format(confidence_intervals)
            raise ValueError(msg)

        self._confidence_intervals = confidence_intervals


    def set_cross_sections(self, cross_sections):

        if not is_string(cross_sections):
            msg = 'Unable to set cross sections to non-string ' \
                  'value {0}'.format(cross_sections)
            raise ValueError(msg)

        self._cross_sections = cross_sections


    def set_energy_grid(self, energy_grid):

        if not energy_grid in ['union', 'nuclide']:
            msg = 'Unable to set energy grid to {0} which is neither ' \
                  'union nor nuclide'.format(energy_grid)
            raise ValueError(msg)

        self._energy_grid = energy_grid


    def set_ptables(self, ptables):

        if not isinstance(ptables, (bool, np.bool)):
            msg = 'Unable to set ptables to non-boolean ' \
                  'value {0}'.format(ptables)
            raise ValueError(msg)

        self._ptables = ptables


    def set_run_cmfd(self, run_cmfd):

        if not isinstance(run_cmfd, (bool, np.bool)):
            msg = 'Unable to set run_cmfd to non-boolean ' \
                  'value {0}'.format(run_cmfd)
            raise ValueError(msg)

        self._run_cmfd = run_cmfd


    def set_seed(self, seed):

        if not is_integer(seed):
            msg = 'Unable to set seed to non-integer value {0}'.format(seed)
            raise ValueError(msg)

        elif seed <= 0:
            msg = 'Unable to set seed to non-positive integer {0}'.format(seed)
            raise ValueError(msg)

        self._seed = seed


    def set_survival_biasing(self, survival_biasing):

        if not isinstance(survival_biasing, (bool, np.bool)):
            msg = 'Unable to set survival biasing to non-boolean ' \
                  'value {0}'.format(survival_biasing)
            raise ValueError(msg)

        self._survival_biasing = survival_biasing


    def set_weight(self, weight, weight_avg):

        if not is_float(weight):
            msg = 'Unable to set weight cutoff to non-floating point ' \
                  'value {0}'.format(weight)
            raise ValueError(msg)

        elif not is_float(weight_avg):
            msg = 'Unable to set weight avg. to non-floating point ' \
                  'value {0}'.format(weight_avg)
            raise ValueError(msg)

        elif weight < 0.0:
            msg = 'Unable to set weight cutoff to negative ' \
                  'value {0}'.format(weight)
            raise ValueError(msg)

        elif weight_avg < 0.0:
            msg = 'Unable to set weight avg. to negative ' \
                  'value {0}'.format(weight_avg)
            raise ValueError(msg)

        self._weight = weight
        self._weight_avg = weight_avg


    def set_entropy_dimension(self, dimension):

        if not isinstance(dimension, (tuple, list)):
            msg = 'Unable to set entropy mesh dimension to {0} which is ' \
                  'not a Python tuple or list'.format(dimension)
            raise ValueError(msg)

        elif len(dimension) < 3 or len(dimension) > 3:
            msg = 'Unable to set entropy mesh dimension to {0} which is ' \
                  'not a set of 3 integer dimensions'.format(dimension)
            raise ValueError(msg)

        for dim in dimension:

            if not is_integer(dim) and not is_float(dim):
                msg = 'Unable to set entropy mesh dimension to a ' \
                      'non-integer or floating point value {0}'.format(dim)
                raise ValueError(msg)

        self._entropy_dimension = dimension


    def set_entropy_lower_left(self, lower_left):

        if not isinstance(lower_left, (tuple, list)):
            msg = 'Unable to set entropy mesh lower left corner to {0} which ' \
                  'is not a Python tuple or list'.format(lower_left)
            raise ValueError(msg)

        elif len(lower_left) < 3 or len(lower_left) > 3:
            msg = 'Unable to set entropy mesh lower left corner to {0} which ' \
                  'is not a 3D point'.format(lower_left)
            raise ValueError(msg)

        for coord in lower_left:

            if not is_integer(coord) and not is_float(coord):

                msg = 'Unable to set entropy mesh lower left corner to a ' \
                      'non-integer or floating point value {0}'.format(coord)
                raise ValueError(msg)

        self._entropy_lower_left = lower_left


    def set_entropy_upper_right(self, upper_right):

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


    def set_no_reduce(self, no_reduce):

        if not isinstance(no_reduce, (bool, np.bool)):
            msg = 'Unable to set the no_reduce to a non-boolean ' \
                  'value {0}'.format(no_reduce)
            raise ValueError(msg)

        self._no_reduce = no_reduce


    def set_threads(self, threads):

        if not is_integer(threads):
            msg = 'Unable to set the threads to a non-integer ' \
                  'value {0}'.format(threads)
            raise ValueError(msg)

        elif threads <=0:
            msg = 'Unable to set the threads to a negative ' \
                  'value {0}'.format(threads)
            raise ValueError(msg)

        self._threads = threads


    def set_trace(self, trace):

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


    def set_track(self, track):

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


    def set_ufs_dimension(self, dimension):

        if not is_integer(dimension) and not is_float(dimension):
            msg = 'Unable to set UFS dimension to non-integer or ' \
                  'non-floating point value {0}'.format(dimension)
            raise ValueError(msg)

        elif dimension < 1:
            msg = 'Unable to set UFS dimension to value {0} which is ' \
                  'less than one'.format(dimension)
            raise ValueError(msg)

        self._ufs_dimension = dimension


    def set_ufs_lower_left(self, lower_left):

        if not isinstance(lower_left, (tuple, list, np.ndarray)):
            msg = 'Unable to set UFS mesh lower left corner to {0} which is ' \
                  'not a Python tuple or list'.format(lower_left)
            raise ValueError(msg)

        elif len(lower_left) < 3 or len(lower_left) > 3:
            msg = 'Unable to set UFS mesh lower left corner to {0} which ' \
                  'is not a 3D point'.format(lower_left)
            raise ValueError(msg)

        self._ufs_lower_left = lower_left


    def set_ufs_upper_right(self, upper_right):

        if not isinstance(upper_right, tuple) and \
          not isinstance(upper_right, list):
            msg = 'Unable to set UFs mesh upper right corner to {0} which is ' \
                  'not a Python tuple or list'.format(upper_right)
            raise ValueError(msg)

        if len(upper_right) < 3 or len(upper_right) > 3:
            msg = 'Unable to set UFS mesh upper right corner to {0} which ' \
                  'is not a 3D point'.format(upper_right)
            raise ValueError(msg)

        self._ufs_upper_right = upper_right


    def set_dd_mesh_dimension(self, dimension):

        # TODO: remove this when domain decomposition is merged
        warnings.warn('This feature is not yet implemented in a release ' \
                      'version of openmc')

        if not isinstance(dimension, tuple) and \
          not isinstance(dimension, list):
            msg = 'Unable to set DD mesh upper right corner to {0} which is ' \
                  'not a Python tuple or list'.format(dimension)
            raise ValueError(msg)

        if len(dimension) < 3 or len(dimension) > 3:
            msg = 'Unable to set DD mesh upper right corner to {0} which ' \
                  'is not a 3D point'.format(dimension)
            raise ValueError(msg)

        self._dd_mesh_dimension = dimension


    def set_dd_mesh_lower_left(self, lower_left):

        # TODO: remove this when domain decomposition is merged
        warnings.warn('This feature is not yet implemented in a release ' \
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


    def set_dd_mesh_upper_right(self, upper_right):

        # TODO: remove this when domain decomposition is merged
        warnings.warn('This feature is not yet implemented in a release ' \
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


    def set_dd_nodemap(self, nodemap):

        # TODO: remove this when domain decomposition is merged
        warnings.warn('This feature is not yet implemented in a release ' \
                      'version of openmc')

        if not isinstance(nodemap, tuple) and \
          not isinstance(nodemap, list):
            msg = 'Unable to set DD nodemap {0} which is ' \
                  'not a Python tuple or list'.format(dimension)
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

        self._dd_mesh_dimension = dimension

    def set_dd_allow_leakage(self, allow):

        # TODO: remove this when domain decomposition is merged
        warnings.warn('This feature is not yet implemented in a release ' \
                      'version of openmc')

        if not type(allow) == bool:
            msg = 'Unable to set DD allow_leakage {0} which is ' \
                  'not a Python bool'.format(dimension)
            raise ValueError(msg)

        self._dd_allow_leakage = allow

    def create_eigenvalue_subelement(self):

        self.create_particles_subelement()
        self.create_batches_subelement()
        self.create_inactive_subelement()
        self.create_generations_per_batch_subelement()


    def create_batches_subelement(self):

        if not self._batches is None:

            if self._eigenvalue_element is None:
                self._eigenvalue_element = ET.SubElement(self._settings_file,
                                                         "eigenvalue")

            element = ET.SubElement(self._eigenvalue_element, "batches")
            element.text = '{0}'.format(self._batches)


    def create_generations_per_batch_subelement(self):

        if not self._generations_per_batch is None:

            if self._eigenvalue_element is None:
                self._eigenvalue_element = ET.SubElement(self._settings_file,
                                                         "eigenvalue")

            element = ET.SubElement(self._eigenvalue_element,
                                    "generations_per_batch")
            element.text = '{0}'.format(self._generations_per_batch)


    def create_inactive_subelement(self):

        if not self._inactive is None:

            if self._eigenvalue_element is None:
                self._eigenvalue_element = ET.SubElement(self._settings_file,
                                                         "eigenvalue")

            element = ET.SubElement(self._eigenvalue_element, "inactive")
            element.text = '{0}'.format(self._inactive)


    def create_particles_subelement(self):

        if not self._particles is None:

            if self._eigenvalue_element is None:
                self._eigenvalue_element = ET.SubElement(self._settings_file,
                                                         "eigenvalue")

            element = ET.SubElement(self._eigenvalue_element, "particles")
            element.text = '{0}'.format(self._particles)


    def create_source_subelement(self):

        self.create_source_space_subelement()
        self.create_source_energy_subelement()
        self.create_source_angle_subelement()


    def create_source_space_subelement(self):


        if not self._source_space_params is None:

            if self._source_subelement is None:
                self._source_subelement = ET.SubElement(self._settings_file,
                                                        "source")

            element = ET.SubElement(self._source_subelement, "space")
            element.set("type", self._source_space_type)

            subelement = ET.SubElement(element, "parameters")

            text = ''
            for param in self._source_space_params:
                text += '{0} '.format(param)
            subelement.text = text.rstrip(' ')


    def create_source_angle_subelement(self):

        if not self._source_angle_params is None:

            if self._source_subelement is None:
                self._source_subelement = ET.SubElement(self._settings_file,
                                                        "source")

            element = ET.SubElement(self._source_subelement, "angle")
            element.set("type", self._source_angle_type)

            subelement = ET.SubElement(element, "parameters")

            text = ''
            for param in self._source_angle_params:
                text += '{0} '.format(param)
            subelement.text = text.rstrip(' ')


    def create_source_energy_subelement(self):

        if not self._source_energy_params is None:

            if self._source_subelement is None:
                self._source_subelement = ET.SubElement(self._settings_file,
                                                        "source")

            element = ET.SubElement(self._source_subelement, "energy")
            element.set("type", self._source_energy_type)

            subelement = ET.SubElement(element, "parameters")

            text = ''
            for param in self._source_energy_params:
                text += '{0} '.format(param)
            subelement.text = text.rstrip(' ')


    def create_output_subelement(self):

        if not self._output is None:
            element = ET.SubElement(self._settings_file, "output")

            for key in self._output:
                subelement = ET.SubElement(element, key)
                subelement.text = str(self._output[key]).lower()

            self.create_output_path_subelement()


    def create_output_path_subelement(self):

        if not self._output_path is None:
            element = ET.SubElement(self._settings_file, "output_path")
            element.text = self._output_path


    def create_verbosity_subelement(self):

        if not self._verbosity is None:
            element = ET.SubElement(self._settings_file, "verbosity")
            element.text = '{0}'.format(self._verbosity)


    def create_statepoint_subelement(self):

        # Batches subelement
        if not self._statepoint_batches is None:
            element = ET.SubElement(self._settings_file, "state_point")
            subelement = ET.SubElement(element, "batches")
            text = ''
            for batch in self._statepoint_batches:
                text += '{0} '.format(batch)
            subelement.text = text.rstrip(' ')

        # Interval subelement
        elif not self._statepoint_interval is None:
            element = ET.SubElement(self._settings_file, "state_point")
            subelement = ET.SubElement(element, "interval")
            subelement.text = '{0}'.format(self._statepoint_interval)


    def create_sourcepoint_subelement(self):

        # Batches subelement
        if not self._sourcepoint_batches is None:
            element = ET.SubElement(self._settings_file, "source_point")
            subelement = ET.SubElement(element, "batches")
            text = ''
            for batch in self._sourcepoint_batches:
                text += '{0} '.format(batch)
            subelement.text = text.rstrip(' ')

        # Interval subelement
        elif not self._sourcepoint_interval is None:
            element = ET.SubElement(self._settings_file, "source_point")
            subelement = ET.SubElement(element, "interval")
            subelement.text = '{0}'.format(self._sourcepoint_interval)

        # Separate subelement
        if not self._sourcepoint_separate is None:
            subelement = ET.SubElement(element, "separate")
            subelement.text = '{0}'.format(str(self._sourcepoint_separate).lower())

        # Write subelement
        if not self._sourcepoint_write is None:
            subelement = ET.SubElement(element, "write")
            subelement.text = '{0}'.format(str(self._sourcepoint_write).lower())

        # Overwrite latest subelement
        if not self._sourcepoint_overwrite is None:
            subelement = ET.SubElement(element, "overwrite_latest")
            subelement.text = '{0}'.format(str(self._sourcepoint_overwrite).lower())


    def create_confidence_intervals(self):

        if not self._confidence_intervals is None:
            element = ET.SubElement(self._settings_file, "confidence_intervals")
            element.text = '{0}'.format(str(self._confidence_intervals).lower())


    def create_cross_sections_subelement(self):

        if not self._cross_sections is None:
            element = ET.SubElement(self._settings_file, "cross_sections")
            element.text = '{0}'.format(self._cross_sections)


    def create_energy_grid_subelement(self):

        if not self._energy_grid is None:
            element = ET.SubElement(self._settings_file, "energy_grid")
            element.text = '{0}'.format(self._energy_grid)


    def create_ptables_subelement(self):

        if not self._ptables is None:
            element = ET.SubElement(self._settings_file, "ptables")
            element.text = '{0}'.format(str(self._ptables).lower())


    def create_run_cmfd_subelement(self):

        if not self._run_cmfd is None:
            element = ET.SubElement(self._settings_file, "run_cmfd")
            element.text = '{0}'.format(str(self._run_cmfd).lower())


    def create_seed_subelement(self):

        if not self._seed is None:
            element = ET.SubElement(self._settings_file, "seed")
            element.text = '{0}'.format(self._seed)


    def create_survival_biasing_subelement(self):

        if not self._survival_biasing is None:
            element = ET.SubElement(self._settings_file, "survival_biasing")
            element.text = '{0}'.format(str(self._survival_biasing).lower())


    def create_cutoff_subelement(self):

        if not self._weight is None:
            element = ET.SubElement(self._settings_file, "cutoff")

            subelement = ET.SubElement(element, "weight")
            subelement.text = '{0}'.format(self._weight)

            subelement = ET.SubElement(element, "weight_avg")
            subelement.text = '{0}'.format(self._weight_avg)


    def create_entropy_subelement(self):

        if not self._entropy_lower_left is None and \
            not self._entropy_upper_right is None:

            element = ET.SubElement(self._settings_file, "entropy")

            subelement = ET.SubElement(element, "dimension")
            subelement.text = '{0} {1} {2}'.format(self._entropy_dimension[0],
                                                   self._entropy_dimension[1],
                                                   self._entropy_dimension[2])

            subelement = ET.SubElement(element, "lower_left")
            subelement.text = '{0} {1} {2}'.format(self._entropy_lower_left[0],
                                                   self._entropy_lower_left[1],
                                                   self._entropy_lower_left[2])

            subelement = ET.SubElement(element, "upper_right")
            subelement.text = '{0} {1} {2}'.format(self._entropy_upper_right[0],
                                                   self._entropy_upper_right[1],
                                                   self._entropy_upper_right[2])


    def create_no_reduce_subelement(self):

        if not self._no_reduce is None:
            element = ET.SubElement(self._settings_file, "no_reduce")
            element.text = '{0}'.format(str(self._no_reduce).lower())


    def create_threads_subelement(self):

        if not self._threads is None:
            element = ET.SubElement(self._settings_file, "threads")
            element.text = '{0}'.format(self._threads)


    def create_trace_subelement(self):

        if not self._trace is None:
            element = ET.SubElement(self._settings_file, "trace")

            text = ''
            for item in self._trace:
                text += '{0} '.format(item)
            element.text = text.rstrip(' ')


    def create_track_subelement(self):

        if not self._track is None:
            element = ET.SubElement(self._settings_file, "track")

            text = ''
            for item in self._track:
                text += '{0} '.format(item)
            element.text = text.rstrip(' ')


    def create_ufs_subelement(self):

        if not self._ufs_lower_left is None and \
            not self._ufs_upper_right is None:

            element = ET.SubElement(self._settings_file, "uniform_fs")

            subelement = ET.SubElement(element, "dimension")
            subelement.text = '{0}'.format(self._ufs_dimension)

            subelement = ET.SubElement(element, "lower_left")
            subelement.text = '{0} {1} {2}'.format(self._ufs_lower_left[0],
                                                   self._ufs_lower_left[1],
                                                   self._ufs_lower_left[2])

            subelement = ET.SubElement(element, "upper_right")
            subelement.text = '{0} {1} {2}'.format(self._ufs_upper_right[0],
                                                   self._ufs_upper_right[1],
                                                   self._ufs_upper_right[2])


    def create_dd_subelement(self):

        if not self._dd_mesh_lower_left is None and \
            not self._dd_mesh_upper_right is None and \
            not self._dd_mesh_dimension is None:

            element = ET.SubElement(self._settings_file, "domain_decomposition")

            subelement = ET.SubElement(element, "mesh")
            subsubelement = ET.SubElement(subelement, "dimension")
            subsubelement.text = '{0} {1} {2}'.format(
                    self._dd_mesh_dimension[0],
                    self._dd_mesh_dimension[1],
                    self._dd_mesh_dimension[2])

            subsubelement = ET.SubElement(subelement, "lower_left")
            subsubelement.text = '{0} {1} {2}'.format(
                    self._dd_mesh_lower_left[0],
                    self._dd_mesh_lower_left[1],
                    self._dd_mesh_lower_left[2])

            subsubelement = ET.SubElement(subelement, "upper_right")
            subsubelement.text = '{0} {1} {2}'.format(
                    self._dd_mesh_upper_right[0],
                    self._dd_mesh_upper_right[1],
                    self._dd_mesh_upper_right[2])

            if not self._dd_nodemap is None:
                subelement = ET.SubElement(element, "nodemap")
                subsubelement.text = ' '.join([str(n) for n in self._dd_nodemap])

            subelement = ET.SubElement(element, "allow_leakage")
            if self._dd_allow_leakage:
                subelement.text = 'true'
            else:
                subelement.text = 'false'
                

    def export_to_xml(self):

        self.create_eigenvalue_subelement()
        self.create_source_subelement()
        self.create_output_subelement()
        self.create_statepoint_subelement()
        self.create_sourcepoint_subelement()
        self.create_confidence_intervals()
        self.create_cross_sections_subelement()
        self.create_energy_grid_subelement()
        self.create_ptables_subelement()
        self.create_run_cmfd_subelement()
        self.create_seed_subelement()
        self.create_survival_biasing_subelement()
        self.create_cutoff_subelement()
        self.create_entropy_subelement()
        self.create_no_reduce_subelement()
        self.create_threads_subelement()
        self.create_verbosity_subelement()
        self.create_trace_subelement()
        self.create_track_subelement()
        self.create_ufs_subelement()
        self.create_dd_subelement()

        # Clean the indentation in the file to be user-readable
        clean_xml_indentation(self._settings_file)

        # Write the XML Tree to the settings.xml file
        tree = ET.ElementTree(self._settings_file)
        tree.write("settings.xml", xml_declaration=True,
                             encoding='utf-8', method="xml")
