#!/usr/bin/env python

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

        self._settings_file = ET.Element("settings")
        self._eigenvalue_element = None
        self._source_element = None


    def setBatches(self, batches):

        if not is_integer(batches):
            msg = 'Unable to set batches to a non-integer ' \
                  'value {0}'.format(batches)
            raise ValueError(msg)

        if batches <= 0:
            msg = 'Unable to set batches to a negative ' \
                  'value {0}'.format(batches)
            raise ValueError(msg)

        self._batches = batches


    def setGenerationsPerBatch(self, generations_per_batch):

        if not is_integer(generations_per_batch):
            msg = 'Unable to set generations per batch to a non-integer ' \
                  'value {0}'.format(generations_per_batch)
            raise ValueError(msg)

        if generations_per_batch <= 0:
            msg = 'Unable to set generations per batch to a negative ' \
                  'value {0}'.format(generations_per_batch)
            raise ValueError(msg)

        self._generations_per_batch = generations_per_batch


    def setInactive(self, inactive):

        if not is_integer(inactive):
            msg = 'Unable to set inactive batches to a non-integer ' \
                  'value {0}'.format(inactive)
            raise ValueError(msg)

        if inactive <= 0:
            msg = 'Unable to set inactive batches to a negative ' \
                  'value {0}'.format(inactive)
            raise ValueError(msg)

        self._inactive = inactive


    def setParticles(self, particles):

        if not is_integer(particles):
            msg = 'Unable to set particles to a non-integer ' \
                  'value {0}'.format(particles)
            raise ValueError(msg)

        if particles <= 0:
            msg = 'Unable to set particles to a negative ' \
                  'value {0}'.format(particles)
            raise ValueError(msg)

        self._particles = particles


    def setSourceFile(self, source_file):

        if not is_string(source_file):
            msg = 'Unable to set source file to a non-string ' \
                  'value {0}'.format(source_file)
            raise ValueError(msg)

        self._source_file = source_file


    def setSourceSpace(self, type, params):

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


    def setSourceAngle(self, type, params=[]):

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


    def setSourceEnergy(self, type, params=[]):

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


    def setOutput(self, output):

        if not isinstance(output, dict):
            msg = 'Unable to set output to {0} which is not a Python ' \
                  'dictionary of string keys and boolean values'.format(output)
            raise ValueError(msg)


        for element in output:

            if not element in ['summary', 'cross_sections', 'tallies']:
                msg = 'Unable to set output to {0} which is unsupported by ' \
                      'OpenMC'.format(element)
                raise ValueError(msg)

            if not isinstance(output[element], (bool, np.bool)):
                msg = 'Unable to set output for {0} to a non-boolean ' \
                      'value {1}'.format(element, output[element])
                raise ValueError(msg)

        self._output = output


    def setOutputPath(self, output_path):

        if not is_string(output_path):
            msg = 'Unable to set output path to non-string ' \
                  'value {0}'.format(output_path)
            raise ValueError(msg)

        self._output_path = output_path


    def setVerbosity(self, verbosity):

        if not is_integer(verbosity):
            msg = 'Unable to set verbosity to non-integer ' \
                  'value {0}'.format(verbosity)
            raise ValueError(msg)

        if verbosity < 1 or verbosity > 10:
            msg = 'Unable to set verbosity to {0} which is not between ' \
                  '1 and 10'.format(verbosity)
            raise ValueError(msg)

        self._verbosity = verbosity


    def setStatepointBatches(self, batches):

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


    def setStatepointInterval(self, interval):

        if not is_integer(interval):
            msg = 'Unable to set statepoint interval to non-integer ' \
                        'value {0}'.format(interval)
            raise ValueError(msg)

        self._statepoint_interval = interval


    def setSourcepointBatches(self, batches):

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


    def setSourcepointInterval(self, interval):

        if not is_integer(interval):
            msg = 'Unable to set sourcepoint interval to non-integer ' \
                  'value {0}'.format(interval)
            raise ValueError(msg)

        self._sourcepoint_interval = interval


    def setSourcepointSeparate(self, source_separate):

        if not isinstance(source_separate, (bool, np.bool)):
            msg = 'Unable to set sourcepoint separate to non-boolean ' \
                  'value {0}'.format(source_separate)
            raise ValueError(msg)

        self._sourcepoint_separate = source_separate


    def setSourcepointWrite(self, source_write):

        if not isinstance(source_write, (bool, np.bool)):
            msg = 'Unable to set sourcepoint write to non-boolean ' \
                  'value {0}'.format(source_write)
            raise ValueError(msg)

        self._sourcepoint_write = source_write


    def setSourcepointOverwrite(self, source_overwrite):

        if not isinstance(source_overwrite, (bool, np.bool)):
            msg = 'Unable to set sourcepoint overwrite to non-boolean ' \
                  'value {0}'.format(source_overwrite)
            raise ValueError(msg)

        self._sourcepoint_overwrite = source_overwrite


    def setConfidenceIntervals(self, confidence_intervals):

        if not isinstance(confidence_intervals, (bool, np.bool)):
            msg = 'Unable to set confidence interval to non-boolean ' \
                  'value {0}'.format(confidence_intervals)
            raise ValueError(msg)

        self._confidence_intervals = confidence_intervals


    def setCrossSections(self, cross_sections):

        if not is_string(cross_sections):
            msg = 'Unable to set cross sections to non-string ' \
                  'value {0}'.format(cross_sections)
            raise ValueError(msg)

        self._cross_sections = cross_sections


    def setEnergyGrid(self, energy_grid):

        if not energy_grid in ['union', 'nuclide']:
            msg = 'Unable to set energy grid to {0} which is neither ' \
                  'union nor nuclide'.format(energy_grid)
            raise ValueError(msg)

        self._energy_grid = energy_grid


    def setPTables(self, ptables):

        if not isinstance(ptables, (bool, np.bool)):
            msg = 'Unable to set ptables to non-boolean ' \
                  'value {0}'.format(ptables)
            raise ValueError(msg)

        self._ptables = ptables


    def setRunCMFD(self, run_cmfd):

        if not isinstance(run_cmfd, (bool, np.bool)):
            msg = 'Unable to set run_cmfd to non-boolean ' \
                  'value {0}'.format(run_cmfd)
            raise ValueError(msg)

        self._run_cmfd = run_cmfd


    def setSeed(self, seed):

        if not is_integer(seed):
            msg = 'Unable to set seed to non-integer value {0}'.format(seed)
            raise ValueError(msg)

        elif seed <= 0:
            msg = 'Unable to set seed to non-positive integer {0}'.format(seed)
            raise ValueError(msg)

        self._seed = seed


    def setSurvivalBiasing(self, survival_biasing):

        if not isinstance(survival_biasing, (bool, np.bool)):
            msg = 'Unable to set survival biasing to non-boolean ' \
                  'value {0}'.format(survival_biasing)
            raise ValueError(msg)

        self._survival_biasing = survival_biasing


    def setWeight(self, weight, weight_avg):

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


    def setEntropyDimension(self, dimension):

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


    def setEntropyLowerLeft(self, lower_left):

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


    def setEntropyUpperRight(self, upper_right):

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


    def setNoReduce(self, no_reduce):

        if not isinstance(no_reduce, (bool, np.bool)):
            msg = 'Unable to set the no_reduce to a non-boolean ' \
                  'value {0}'.format(no_reduce)
            raise ValueError(msg)

        self._no_reduce = no_reduce


    def setThreads(self, threads):

        if not is_integer(threads):
            msg = 'Unable to set the threads to a non-integer ' \
                  'value {0}'.format(threads)
            raise ValueError(msg)

        elif threads <=0:
            msg = 'Unable to set the threads to a negative ' \
                  'value {0}'.format(threads)
            raise ValueError(msg)

        self._threads = threads


    def setTrace(self, trace):

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


    def setTrack(self, track):

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


    def setUFSDimension(self, dimension):

        if not is_integer(dimension) and not is_float(dimension):
            msg = 'Unable to set UFS dimension to non-integer or ' \
                  'non-floating point value {0}'.format(dimension)
            raise ValueError(msg)

        elif dimension < 1:
            msg = 'Unable to set UFS dimension to value {0} which is ' \
                  'less than one'.format(dimension)
            raise ValueError(msg)

        self._ufs_dimension = dimension


    def setUFSLowerLeft(self, lower_left):

        if not isinstance(lower_left, (tuple, list, np.ndarray)):
            msg = 'Unable to set UFS mesh lower left corner to {0} which is ' \
                  'not a Python tuple or list'.format(lower_left)
            raise ValueError(msg)

        elif len(lower_left) < 3 or len(lower_left) > 3:
            msg = 'Unable to set UFS mesh lower left corner to {0} which ' \
                  'is not a 3D point'.format(lower_left)
            raise ValueError(msg)

        self._ufs_lower_left = lower_left


    def setUFSUpperRight(self, upper_right):

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


    def createEigenvalueSubelement(self):

        self.createParticlesSubelement()
        self.createBatchesSubelement()
        self.createInactiveSubelement()
        self.createGenerationsPerBatchSubelement()


    def createBatchesSubelement(self):

        if not self._batches is None:

            if self._eigenvalue_element is None:
                self._eigenvalue_element = ET.SubElement(self._settings_file,
                                                         "eigenvalue")

            element = ET.SubElement(self._eigenvalue_element, "batches")
            element.text = '{0}'.format(self._batches)


    def createGenerationsPerBatchSubelement(self):

        if not self._generations_per_batch is None:

            if self._eigenvalue_element is None:
                self._eigenvalue_element = ET.SubElement(self._settings_file,
                                                         "eigenvalue")

            element = ET.SubElement(self._eigenvalue_element,
                                    "generations_per_batch")
            element.text = '{0}'.format(self._generations_per_batch)


    def createInactiveSubelement(self):

        if not self._inactive is None:

            if self._eigenvalue_element is None:
                self._eigenvalue_element = ET.SubElement(self._settings_file,
                                                         "eigenvalue")

            element = ET.SubElement(self._eigenvalue_element, "inactive")
            element.text = '{0}'.format(self._inactive)


    def createParticlesSubelement(self):

        if not self._particles is None:

            if self._eigenvalue_element is None:
                self._eigenvalue_element = ET.SubElement(self._settings_file,
                                                         "eigenvalue")

            element = ET.SubElement(self._eigenvalue_element, "particles")
            element.text = '{0}'.format(self._particles)


    def createSourceSubelement(self):

        self.createSourceSpaceSubelement()
        self.createSourceEnergySubelement()
        self.createSourceAngleSubelement()


    def createSourceSpaceSubelement(self):


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


    def createSourceAngleSubelement(self):

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


    def createSourceEnergySubelement(self):

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


    def createOutputSubelement(self):

        if not self._output is None:
            element = ET.SubElement(self._settings_file, "output")

            for key in self._output:
                subelement = ET.SubElement(element, key)
                subelement.text = str(self._output[key]).lower()

            self.createOuputPathSubelement()


    def createOuputPathSubelement(self):

        if not self._output_path is None:
            element = ET.SubElement(self._settings_file, "output_path")
            element.text = self._output_path


    def createVerbositySubelement(self):

        if not self._verbosity is None:
            element = ET.SubElement(self._settings_file, "verbosity")
            element.text = '{0}'.format(self._verbosity)


    def createStatepointSubelement(self):

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


    def createSourcepointSubelement(self):

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


    def createConfidenceIntervalsSubelement(self):

        if not self._confidence_intervals is None:
            element = ET.SubElement(self._settings_file, "confidence_intervals")
            element.text = '{0}'.format(str(self._confidence_intervals).lower())


    def createCrossSectionsSubelement(self):

        if not self._cross_sections is None:
            element = ET.SubElement(self._settings_file, "cross_sections")
            element.text = '{0}'.format(self._cross_sections)


    def createEnergyGridSubelement(self):

        if not self._energy_grid is None:
            element = ET.SubElement(self._settings_file, "energy_grid")
            element.text = '{0}'.format(self._energy_grid)


    def createPTablesSubelement(self):

        if not self._ptables is None:
            element = ET.SubElement(self._settings_file, "ptables")
            element.text = '{0}'.format(str(self._ptables).lower())


    def createRunCMFDSubelement(self):

        if not self._run_cmfd is None:
            element = ET.SubElement(self._settings_file, "run_cmfd")
            element.text = '{0}'.format(str(self._run_cmfd).lower())


    def createSeedSubelement(self):

        if not self._seed is None:
            element = ET.SubElement(self._settings_file, "seed")
            element.text = '{0}'.format(self._seed)


    def createSurvivalBiasingSubelement(self):

        if not self._survival_biasing is None:
            element = ET.SubElement(self._settings_file, "survival_biasing")
            element.text = '{0}'.format(str(self._survival_biasing).lower())


    def createCutoffSubelement(self):

        if not self._weight is None:
            element = ET.SubElement(self._settings_file, "cutoff")

            subelement = ET.SubElement(element, "weight")
            subelement.text = '{0}'.format(self._weight)

            subelement = ET.SubElement(element, "weight_avg")
            subelement.text = '{0}'.format(self._weight_avg)


    def createEntropySubelement(self):

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


    def createNoReduceSubelement(self):

        if not self._no_reduce is None:
            element = ET.SubElement(self._settings_file, "no_reduce")
            element.text = '{0}'.format(str(self._no_reduce).lower())


    def createThreadsSubelement(self):

        if not self._threads is None:
            element = ET.SubElement(self._settings_file, "threads")
            element.text = '{0}'.format(self._threads)


    def createTraceSubelement(self):

        if not self._trace is None:
            element = ET.SubElement(self._settings_file, "trace")

            text = ''
            for item in self._trace:
                text += '{0} '.format(item)
            element.text = text.rstrip(' ')


    def createTrackSubelement(self):

        if not self._track is None:
            element = ET.SubElement(self._settings_file, "track")

            text = ''
            for item in self._track:
                text += '{0} '.format(item)
            element.text = text.rstrip(' ')


    def createUFSSubelement(self):

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


    def exportToXML(self):

        self.createEigenvalueSubelement()
        self.createSourceSubelement()
        self.createOutputSubelement()
        self.createStatepointSubelement()
        self.createSourcepointSubelement()
        self.createConfidenceIntervalsSubelement()
        self.createCrossSectionsSubelement()
        self.createEnergyGridSubelement()
        self.createPTablesSubelement()
        self.createRunCMFDSubelement()
        self.createSeedSubelement()
        self.createSurvivalBiasingSubelement()
        self.createCutoffSubelement()
        self.createEntropySubelement()
        self.createNoReduceSubelement()
        self.createThreadsSubelement()
        self.createVerbositySubelement()
        self.createTraceSubelement()
        self.createTrackSubelement()
        self.createUFSSubelement()

        # Clean the indentation in the file to be user-readable
        clean_xml_indentation(self._settings_file)

        # Write the XML Tree to the settings.xml file
        tree = ET.ElementTree(self._settings_file)
        tree.write("settings.xml", xml_declaration=True,
                             encoding='utf-8', method="xml")
