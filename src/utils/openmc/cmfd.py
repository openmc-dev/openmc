from xml.etree import ElementTree as ET

import numpy as np

from openmc.checkvalue import *
from openmc.clean_xml import *


class CMFDMesh(object):

    def __init__(self):

        self._lower_left = None
        self._upper_right = None
        self._dimension = None
        self._width = None
        self._energy = None
        self._albedo = None
        self._map = None


    @property
    def lower_left(self):
        return self._lower_left


    @property
    def upper_right(self):
        return self._upper_right


    @property
    def dimension(self):
        return self._dimension


    @property
    def width(self):
        return self._width


    @property
    def energy(self):
        return self._energy


    @property
    def albedo(self):
        return self._albedo


    @property
    def map(self):
        return self._mape


    @lower_left.setter
    def lower_left(self, lower_left):

        if not isinstance(lower_left, (tuple, list, np.ndarray)):
            msg = 'Unable to set CMFD Mesh with lower_left {0} which is ' \
                  'not a Python list, tuple or NumPy array'.format(lower_left)
            raise ValueError(msg)

        elif len(lower_left) != 2 and len(lower_left) != 3:
            msg = 'Unable to set CMFD Mesh with lower_left {0} since it ' \
                   'must include 2 or 3 dimensions'.format(lower_left)
            raise ValueError(msg)

        for coord in lower_left:

            if not is_integer(coord) and not is_float(coord):
                msg = 'Unable to set CMFD Mesh with lower_left {0} which is ' \
                      'not an integer or a floating point value'.format(coord)
                raise ValueError(msg)

        self._lower_left = lower_left


    @upper_right.setter
    def upper_right(self, upper_right):

        if not isinstance(upper_right, (tuple, list, np.ndarray)):
            msg = 'Unable to set CMFD Mesh with upper_right {0} which is ' \
                  'not a Python list, tuple or NumPy array'.format(upper_right)
            raise ValueError(msg)

        if len(upper_right) != 2 and len(upper_right) != 3:
            msg = 'Unable to set CMFD Mesh with upper_right {0} since it ' \
                  'must include 2 or 3 dimensions'.format(upper_right)
            raise ValueError(msg)

        for coord in upper_right:

            if not is_integer(coord) and not is_float(coord):
                msg = 'Unable to set CMFD Mesh with upper_right {0} which ' \
                      'is not an integer or floating point value'.format(coord)
                raise ValueError(msg)

        self._upper_right = upper_right


    @dimension.setter
    def dimension(self, dimension):

        if not isinstance(dimension, (tuple, list, np.ndarray)):
            msg = 'Unable to set CMFD Mesh with dimension {0} which is ' \
                  'not a Python list, tuple or NumPy array'.format(dimension)
            raise ValueError(msg)

        elif len(dimension) != 2 and len(dimension) != 3:
            msg = 'Unable to set CMFD Mesh with dimension {0} since it ' \
                  'must include 2 or 3 dimensions'.format(dimension)
            raise ValueError(msg)

        for dim in dimension:

            if not is_integer(dim):
                msg = 'Unable to set CMFD Mesh with dimension {0} which ' \
                      'is a non-integer'.format(dim)
                raise ValueError(msg)

        self._dimension = dimension


    @width.setter
    def width(self, width):

        if not width is None:

            if not isinstance(width, (tuple, list, np.ndarray)):
                msg = 'Unable to set CMFD Mesh with width {0} which ' \
                      'is not a Python list, tuple or NumPy array'.format(width)
                raise ValueError(msg)

        if len(width) != 2 and len(width) != 3:
            msg = 'Unable to set CMFD Mesh with width {0} since it must ' \
                  'include 2 or 3 dimensions'.format(width)
            raise ValueError(msg)

        for dim in width:

            if not is_integer(dim) and not is_float(dim):
                msg = 'Unable to set CMFD Mesh with width {0} which is ' \
                      'not an integer or floating point value'.format(width)
                raise ValueError(msg)

        self._width = width


    @energy.setter
    def energy(self, energy):

        if not isinstance(energy, (tuple, list, np.ndarray)):
            msg = 'Unable to set CMFD Mesh energy to {0} which is not ' \
                  'a Python tuple/list or NumPy array'.format(energy)
            raise ValueError(msg)

        for e in energy:

            if not is_integer(e) and not is_float(e):
                msg = 'Unable to set CMFD Mesh energy to {0} which is not ' \
                      'an integer or floating point value'.format(e)
                raise ValueError(msg)

            elif e < 0:
                msg = 'Unable to set CMFD Mesh energy to {0} which is ' \
                      'is a negative integer'.format(e)
                raise ValueError(msg)

        self._energy = energy


    @albedo.setter
    def albedo(self, albedo):

        if not isinstance(albedo, (tuple, list, np.ndarray)):
            msg = 'Unable to set CMFD Mesh albedo to {0} which is not ' \
                  'a Python tuple/list or NumPy array'.format(albedo)
            raise ValueError(msg)

        if not len(albedo) == 6:
            msg = 'Unable to set CMFD Mesh albedo to {0} which is not ' \
                  'length 6 for +/-x,y,z'.format(albedo)
            raise ValueError(msg)

        for a in albedo:

            if not is_integer(a) and not is_float(a):
                msg = 'Unable to set CMFD Mesh albedo to {0} which is not ' \
                      'an integer or floating point value'.format(a)
                raise ValueError(msg)

            elif a < 0 or a > 1:
                msg = 'Unable to set CMFD Mesh albedo to {0} which is ' \
                      'is not in [0,1]'.format(a)
                raise ValueError(msg)

        self._albedo = albedo


    @map.setter
    def map(self, map):

        if not isinstance(map, (tuple, list, np.ndarray)):
            msg = 'Unable to set CMFD Mesh map to {0} which is not ' \
                  'a Python tuple/list or NumPy array'.format(map)
            raise ValueError(msg)

        for m in map:

            if m != 1 and m != 2:
                msg = 'Unable to set CMFD Mesh map to {0} which is ' \
                      'is not 1 or 2'.format(m)
                raise ValueError(msg)

        self._map = map


    def get_mesh_xml(self):

        element = ET.Element("mesh")

        if len(self._lower_left) == 2:
            subelement = ET.SubElement(element, "lower_left")
            subelement.text = '{0} {1}'.format(self._lower_left[0],
                                               self._lower_left[1])
        else:
            subelement = ET.SubElement(element, "lower_left")
            subelement.text = '{0} {1} {2}'.format(self._lower_left[0],
                                                   self._lower_left[1],
                                                   self._lower_left[2])

        if not self._upper_right is None:
            if len(self._upper_right) == 2:
                subelement = ET.SubElement(element, "upper_right")
                subelement.text = '{0} {1}'.format(self._upper_right[0],
                                                   self._upper_right[1])
            else:
                subelement = ET.SubElement(element, "upper_right")
                subelement.text = '{0} {1} {2}'.format(self._upper_right[0],
                                                       self._upper_right[1],
                                                       self._upper_right[2])

        if len(self._dimension) == 2:
            subelement = ET.SubElement(element, "dimension")
            subelement.text = '{0} {1}'.format(self._dimension[0],
                                               self._dimension[1])
        else:
            subelement = ET.SubElement(element, "dimension")
            subelement.text = '{0} {1} {2}'.format(self._dimension[0],
                                                   self._dimension[1],
                                                   self._dimension[2])

        if not self._width is None:
            if len(self._width) == 2:
                subelement = ET.SubElement(element, "width")
                subelement.text = '{0} {1}'.format(self._width[0],
                                                   self._width[1])
            else:
                subelement = ET.SubElement(element, "width")
                subelement.text = '{0} {1} {2}'.format(self._width[0],
                                                       self._width[1],
                                                       self._width[2])

        if not self._energy is None:

            subelement = ET.SubElement(element, "energy")

            energy = ''
            for e in self._energy:
                energy += '{0} '.format(e)

            subelement.set("energy", energy.rstrip(' '))

        if not self._albedo is None:

            subelement = ET.SubElement(element, "albedo")

            albedo = ''
            for a in self._albedo:
                albedo += '{0} '.format(a)

            subelement.set("albedo", albedo.rstrip(' '))

        if not self._map is None:

            subelement = ET.SubElement(element, "map")

            map = ''
            for m in self._map:
                map += '{0} '.format(m)

            subelement.set("map", map.rstrip(' '))

        return element


class CMFDFile(object):

    def __init__(self):

        self._active_flush = None
        self._begin = None
        self._display = None
        self._feedback = None
        self._inactive = None
        self._inactive_flush = None
        self._ksp_monitor = None
        self._cmfd_mesh = None
        self._norm = None
        self._num_flushes = None
        self._power_monitor = None
        self._run_adjoint = None
        self._snes_monitor = None
        self._solver = None
        self._write_matrices = None

        self._cmfd_file = ET.Element("cmfd")
        self._cmfd_mesh_element = None


    @property
    def active_flush(self):
        return self._active_flush


    @property
    def begin(self):
        return self._begin


    @property
    def display(self):
        return self._display


    @property
    def feedback(self):
        return self._feedback


    @property
    def inactive(self):
        return self._inactive


    @property
    def inactive_flush(self):
        return self._inactive_flush


    @property
    def ksp_monitor(self):
        return self._ksp_monitor


    @property
    def cmfd_mesh(self):
        return self._cmfd_mesh


    @property
    def norm(self):
        return self._norm


    @property
    def num_flushes(self):
        return self._num_flushes


    @property
    def power_monitor(self):
        return self._power_monitor


    @property
    def run_adjoint(self):
        return self._run_adjoint


    @property
    def snes_monitor(self):
        return self._snes_monitor


    @property
    def solver(self):
        return self._solver


    @property
    def write_matrices(self):
        return self._write_matrices


    @active_flush.setter
    def active_flush(self, active_flush):

        if not is_integer(active_flush):
            msg = 'Unable to set CMFD active flush batch to a non-integer ' \
                  'value {0}'.format(active_flush)
            raise ValueError(msg)

        if active_flush < 0:
            msg = 'Unable to set CMFD active flush batch to a negative ' \
                  'value {0}'.format(active_flush)
            raise ValueError(msg)

        self._active_flush = active_flush


    @begin.setter
    def begin(self, begin):

        if not is_integer(begin):
            msg = 'Unable to set CMFD begin batch to a non-integer ' \
                  'value {0}'.format(begin)
            raise ValueError(msg)

        if begin <= 0:
            msg = 'Unable to set CMFD begin batch batch to a negative ' \
                  'value {0}'.format(begin)
            raise ValueError(msg)

        self._begin = begin


    @display.setter
    def display(self, display):

        if not is_string(display):
            msg = 'Unable to set CMFD display to a non-string ' \
                  'value'.format(display)
            raise ValueError(msg)

        if display not in ['balance', 'dominance', 'entropy', 'source']:
            msg = 'Unable to set CMFD display to {0} which is ' \
                  'not an accepted value'.format(display)
            raise ValueError(msg)

        self._display = display


    @feedback.setter
    def feedback(self, feedback):

        if not isinstance(feedback, bool):
            msg = 'Unable to set CMFD feedback to {0} which is ' \
                  'a non-boolean value'.format(feedback)
            raise ValueError(msg)

        self._feedback = feedback


    @inactive.setter
    def inactive(self, inactive):

        if not isinstance(inactive, bool):
            msg = 'Unable to set CMFD inactive batch to {0} which is ' \
                  ' a non-boolean value'.format(inactive)
            raise ValueError(msg)

        self._inactive = inactive


    @inactive_flush.setter
    def inactive_flush(self, inactive_flush):

        if not is_integer(inactive_flush):
            msg = 'Unable to set CMFD inactive flush batch to {0} which is ' \
                  'a non-integer value'.format(inactive_flush)
            raise ValueError(msg)

        if inactive_flush <= 0:
            msg = 'Unable to set CMFD inactive flush batch to {0} which is ' \
                  'a negative value {0}'.format(inactive_flush)
            raise ValueError(msg)

        self._inactive_flush = inactive_flush


    @ksp_monitor.setter
    def ksp_monitor(self, ksp_monitor):

        if not isinstance(ksp_monitor, bool):
            msg = 'Unable to set CMFD ksp monitor to {0} which is a ' \
                  'non-boolean value'.format(ksp_monitor)
            raise ValueError(msg)

        self._ksp_monitor = ksp_monitor


    @cmfd_mesh.setter
    def cmfd_mesh(self, mesh):

        if not isinstance(mesh, CMFDMesh):
            msg = 'Unable to set CMFD mesh to {0} which is not a ' \
                  'CMFDMesh object'.format(mesh)
            raise ValueError(msg)

        self._mesh = mesh


    @norm.setter
    def norm(self, norm):

        if not is_integer(norm) and not is_float(norm):
            msg = 'Unable to set the CMFD norm to {0} which is not ' \
                  'an integer or floating point value'.format(norm)
            raise ValueError(msg)

        self._norm = norm


    @num_flushes.setter
    def num_flushes(self, num_flushes):

        if not is_integer(num_flushes):
            msg = 'Unable to set the CMFD number of flushes to {0} ' \
                  'which is not an integer value'.format(num_flushes)
            raise ValueError(msg)

        if num_flushes < 0:
            msg = 'Unable to set CMFD number of flushes to a negative ' \
                  'value {0}'.format(num_flushes)
            raise ValueError(msg)

        self._num_flushes = num_flushes


    @power_monitor.setter
    def power_monitor(self, power_monitor):

        if not isinstance(power_monitor, bool):
            msg = 'Unable to set CMFD power monitor to {0} which is a ' \
                  'non-boolean value'.format(power_monitor)
            raise ValueError(msg)

        self._power_monitor = power_monitor


    @run_adjoint.setter
    def run_adjoint(self, run_adjoint):

        if not isinstance(run_adjoint, bool):
            msg = 'Unable to set CMFD run adjoint to {0} which is a ' \
                  'non-boolean value'.format(run_adjoint)
            raise ValueError(msg)

        self._run_adjoint = run_adjoint


    @snes_monitor.setter
    def snes_monitor(self, snes_monitor):

        if not isinstance(snes_monitor, bool):
            msg = 'Unable to set CMFD snes monitor to {0} which is a ' \
                  'non-boolean value'.format(snes_monitor)
            raise ValueError(msg)

        self._snes_monitor = snes_monitor


    @solver.setter
    def solver(self, solver):

        if not solver in ['power', 'jfnk']:
            msg = 'Unable to set CMFD solver to {0} which is not ' \
                  '"power" or "jfnk"'.format(solver)
            raise ValueError(msg)

        self._solver = solver


    @write_matrices.setter
    def write_matrices(self, write_matrices):

        if not isinstance(write_matrices, bool):
            msg = 'Unable to set CMFD write matrices to {0} which is a ' \
                  'non-boolean value'.format(write_matrices)
            raise ValueError(msg)

        self._write_matrices = write_matrices


    def get_active_flush_subelement(self):

        if not self._active_flush is None:
            element = ET.SubElement(self._cmfd_file, "active_flush")
            element.text = '{0}'.format(str(self._active_flush))


    def get_begin_subelement(self):

        if not self._begin is None:
            element = ET.SubElement(self._cmfd_file, "begin")
            element.text = '{0}'.format(str(self._begin))


    def create_display_subelement(self):

        if not self._display is None:
            element = ET.SubElement(self._cmfd_file, "display")
            element.text = '{0}'.format(str(self._display))


    def create_feedback_subelement(self):

        if not self._feedback is None:
            element = ET.SubElement(self._cmfd_file, "feeback")
            element.text = '{0}'.format(str(self._feedback).lower())


    def create_inactive_subelement(self):

        if not self._inactive is None:
            element = ET.SubElement(self._cmfd_file, "inactive")
            element.text = '{0}'.format(str(self._inactive).lower())


    def create_inactive_flush_subelement(self):

        if not self._inactive_flush is None:
            element = ET.SubElement(self._cmfd_file, "inactive_flush")
            element.text = '{0}'.format(str(self._inactive_flush))


    def create_ksp_monitor_subelement(self):

        if not self._ksp_monitor is None:
            element = ET.SubElement(self._cmfd_file, "ksp_monitor")
            element.text = '{0}'.format(str(self._ksp_monitor).lower())


    def create_mesh_subelement(self):

        if not self._mesh is None:
            xml_element = self._mesh.get_mesh_xml()
            self._cmfd_file.append(xml_element)


    def create_norm_subelement(self):

        if not self._num_flushes is None:
            element = ET.SubElement(self._cmfd_file, "norm")
            element.text = '{0}'.format(str(self._norm))


    def create_num_flushes_subelement(self):

        if not self._num_flushes is None:
            element = ET.SubElement(self._cmfd_file, "num_flushes")
            element.text = '{0}'.format(str(self._num_flushes))


    def create_power_monitor_subelement(self):

        if not self._power_monitor is None:
            element = ET.SubElement(self._cmfd_file, "power_monitor")
            element.text = '{0}'.format(str(self._power_monitor).lower())


    def create_run_adjoint_subelement(self):

        if not self._run_adjoint is None:
            element = ET.SubElement(self._cmfd_file, "run_adjoint")
            element.text = '{0}'.format(str(self._run_adjoint).lower())


    def create_snes_monitor_subelement(self):

        if not self._snes_monitor is None:
            element = ET.SubElement(self._cmfd_file, "snes_monitor")
            element.text = '{0}'.format(str(self._snes_monitor).lower())


    def create_solver_subelement(self):

        if not self._solver is None:
            element = ET.SubElement(self._cmfd_file, "solver")
            element.text = '{0}'.format(str(self._solver))


    def create_write_matrices_subelement(self):

        if not self._write_matrices is None:
            element = ET.SubElement(self._cmfd_file, "write_matrices")
            element.text = '{0}'.format(str(self._write_matrices).lower())


    def export_to_xml(self):

        self.create_active_flush_subelement()
        self.create_begin_subelement()
        self.create_display_subelement()
        self.create_feedback_subelement()
        self.create_inactive_subelement()
        self.create_inactive_flush_subelement()
        self.create_ksp_monitor_subelement()
        self.create_mesh_subelement()
        self.create_norm_subelement()
        self.create_num_flushes_subelement()
        self.create_power_monitor_subelement()
        self.create_run_adjoint_subelement()
        self.create_snes_monitor_subelement()
        self.create_solver_subelement()
        self.create_write_matrices_subelement()

        # Clean the indentation in the file to be user-readable
        clean_xml_indentation(self._cmfd_file)

        # Write the XML Tree to the cmfd.xml file
        tree = ET.ElementTree(self._cmfd_file)
        tree.write("cmfd.xml", xml_declaration=True,
                   encoding='utf-8', method="xml")
