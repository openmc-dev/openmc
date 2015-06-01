from xml.etree import ElementTree as ET

import numpy as np

from openmc.checkvalue import *
from openmc.clean_xml import *


class CMFDMesh(object):
    """A structured Cartesian mesh used for coarse mesh finite difference
    acceleration.

    Attributes
    ----------
    lower_left : tuple or list or ndarray
        The lower-left corner of the structured mesh. If only two coordinate are
        given, it is assumed that the mesh is an x-y mesh.
    upper_right : tuple or list or ndarray
        The upper-right corner of the structrued mesh. If only two coordinate
        are given, it is assumed that the mesh is an x-y mesh.
    dimension : tuple or list or ndarray
        The number of mesh cells in each direction.
    width : tuple or list or ndarray
        The width of mesh cells in each direction.
    energy : tuple or list or ndarray
        Energy bins in MeV, listed in ascending order (e.g. [0.0, 0.625e-7,
        20.0]) for CMFD tallies and acceleration. If no energy bins are listed,
        OpenMC automatically assumes a one energy group calculation over the
        entire energy range.
    albedo : tuple or list or ndarray
        Surface ratio of incoming to outgoing partial currents on global
        boundary conditions. They are listed in the following order: -x +x -y +y
        -z +z.
    map : tuple or list or ndarray

        An optional acceleration map can be specified to overlay on the coarse
        mesh spatial grid. If this option is used, a ``1`` is used for a
        non-accelerated region and a ``2`` is used for an accelerated region.
        For a simple 4x4 coarse mesh with a 2x2 fuel lattice surrounded by
        reflector, the map is:

        ::

            [1, 1, 1, 1,
             1, 2, 2, 1,
             1, 2, 2, 1,
             1, 1, 1, 1]

        Therefore a 2x2 system of equations is solved rather than a 4x4. This is
        extremely important to use in reflectors as neutrons will not contribute
        to any tallies far away from fission source neutron regions.  A ``2``
        must be used to identify any fission source region.

    """

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

    def _get_xml_element(self):
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
    """Parameters that control the use of coarse-mesh finite difference acceleration
    in OpenMC. This corresponds directly to the cmfd.xml input file.

    Attributes
    ----------
    begin : int
        Batch number at which CMFD calculations should begin
    display : {'balance', 'dominance', 'entropy', 'source'}
        Set one additional CMFD output column. Options are:

          * "balance" - prints the RMS [%] of the resdiual from the neutron balance
             equation on CMFD tallies.
          * "dominance" - prints the estimated dominance ratio from the CMFD
            iterations.
          * "entropy" - prints the *entropy* of the CMFD predicted fission source.
          * "source" - prints the RMS [%] between the OpenMC fission source and
            CMFD fission source.
    feedback : bool
        Indicate or not the CMFD diffusion result is used to adjust the weight
        of fission source neutrons on the next OpenMC batch. Defaults to False.
    cmfd_mesh : CMFDMesh
        Structured mesh to be used for acceleration
    norm : float
        Normalization factor applied to the CMFD fission source distribution
    power_monitor : bool
        View convergence of power iteration during CMFD acceleration
    run_adjoint : bool
        Perform adjoint calculation on the last batch
    write_matrices : bool
        Write sparse matrices that are used during CMFD acceleration (loss,
        production) to file

    """

    def __init__(self):

        self._begin = None
        self._display = None
        self._feedback = None
        self._cmfd_mesh = None
        self._norm = None
        self._power_monitor = None
        self._run_adjoint = None
        self._write_matrices = None

        self._cmfd_file = ET.Element("cmfd")
        self._cmfd_mesh_element = None

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
    def cmfd_mesh(self):
        return self._cmfd_mesh

    @property
    def norm(self):
        return self._norm

    @property
    def power_monitor(self):
        return self._power_monitor

    @property
    def run_adjoint(self):
        return self._run_adjoint

    @property
    def solver(self):
        return self._solver

    @property
    def write_matrices(self):
        return self._write_matrices

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

    @write_matrices.setter
    def write_matrices(self, write_matrices):

        if not isinstance(write_matrices, bool):
            msg = 'Unable to set CMFD write matrices to {0} which is a ' \
                  'non-boolean value'.format(write_matrices)
            raise ValueError(msg)

        self._write_matrices = write_matrices

    def _create_begin_subelement(self):

        if not self._begin is None:
            element = ET.SubElement(self._cmfd_file, "begin")
            element.text = '{0}'.format(str(self._begin))

    def _create_display_subelement(self):

        if not self._display is None:
            element = ET.SubElement(self._cmfd_file, "display")
            element.text = '{0}'.format(str(self._display))

    def _create_feedback_subelement(self):

        if not self._feedback is None:
            element = ET.SubElement(self._cmfd_file, "feeback")
            element.text = '{0}'.format(str(self._feedback).lower())

    def _create_mesh_subelement(self):

        if not self._mesh is None:
            xml_element = self._mesh._get_xml_element()
            self._cmfd_file.append(xml_element)

    def _create_norm_subelement(self):

        if not self._norm is None:
            element = ET.SubElement(self._cmfd_file, "norm")
            element.text = '{0}'.format(str(self._norm))

    def _create_power_monitor_subelement(self):

        if not self._power_monitor is None:
            element = ET.SubElement(self._cmfd_file, "power_monitor")
            element.text = '{0}'.format(str(self._power_monitor).lower())

    def _create_run_adjoint_subelement(self):

        if not self._run_adjoint is None:
            element = ET.SubElement(self._cmfd_file, "run_adjoint")
            element.text = '{0}'.format(str(self._run_adjoint).lower())

    def _create_write_matrices_subelement(self):

        if not self._write_matrices is None:
            element = ET.SubElement(self._cmfd_file, "write_matrices")
            element.text = '{0}'.format(str(self._write_matrices).lower())

    def export_to_xml(self):
        """Create a cmfd.xml file using the class data that can be used for an OpenMC
        simulation.

        """

        self._create_begin_subelement()
        self._create_display_subelement()
        self._create_feedback_subelement()
        self._create_mesh_subelement()
        self._create_norm_subelement()
        self._create_power_monitor_subelement()
        self._create_run_adjoint_subelement()
        self._create_write_matrices_subelement()

        # Clean the indentation in the file to be user-readable
        clean_xml_indentation(self._cmfd_file)

        # Write the XML Tree to the cmfd.xml file
        tree = ET.ElementTree(self._cmfd_file)
        tree.write("cmfd.xml", xml_declaration=True,
                   encoding='utf-8', method="xml")
