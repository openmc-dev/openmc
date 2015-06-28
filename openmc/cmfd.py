"""This module can be used to specify parameters used for coarse mesh finite
difference (CMFD) acceleration in OpenMC. CMFD was first proposed by [Smith]_
and is widely used in accelerating neutron transport problems.

References
----------

.. [Smith] K. Smith, "Nodal method storage reduction by non-linear
   iteration", *Trans. Am. Nucl. Soc.*, **44**, 265 (1983).

"""

from collections import Sequence
from numbers import Real, Integral
from xml.etree import ElementTree as ET
import sys

import numpy as np

from openmc.clean_xml import *

if sys.version_info[0] >= 3:
    basestring = str


class CMFDMesh(object):
    """A structured Cartesian mesh used for Coarse Mesh Finite Difference (CMFD)
    acceleration.

    Attributes
    ----------
    lower_left : Sequence of float
        The lower-left corner of the structured mesh. If only two coordinates are
        given, it is assumed that the mesh is an x-y mesh.
    upper_right : Sequence of float
        The upper-right corner of the structrued mesh. If only two coordinates
        are given, it is assumed that the mesh is an x-y mesh.
    dimension : Sequence of int
        The number of mesh cells in each direction.
    width : Sequence of float
        The width of mesh cells in each direction.
    energy : Sequence of float
        Energy bins in MeV, listed in ascending order (e.g. [0.0, 0.625e-7,
        20.0]) for CMFD tallies and acceleration. If no energy bins are listed,
        OpenMC automatically assumes a one energy group calculation over the
        entire energy range.
    albedo : Sequence of float
        Surface ratio of incoming to outgoing partial currents on global
        boundary conditions. They are listed in the following order: -x +x -y +y
        -z +z.
    map : Sequence of int
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
        return self._map

    @lower_left.setter
    def lower_left(self, lower_left):
        if not isinstance(lower_left, Sequence):
            msg = 'Unable to set CMFD Mesh with lower_left {0} which is ' \
                  'not a sequence'.format(lower_left)
            raise ValueError(msg)

        elif len(lower_left) != 2 and len(lower_left) != 3:
            msg = 'Unable to set CMFD Mesh with lower_left {0} since it ' \
                  'must include 2 or 3 dimensions'.format(lower_left)
            raise ValueError(msg)

        for coord in lower_left:
            if not isinstance(coord, Real):
                msg = 'Unable to set CMFD Mesh with lower_left {0} which is ' \
                      'not a real number'.format(coord)
                raise ValueError(msg)

        self._lower_left = lower_left

    @upper_right.setter
    def upper_right(self, upper_right):

        if not isinstance(upper_right, Sequence):
            msg = 'Unable to set CMFD Mesh with upper_right {0} which is ' \
                  'not a sequence'.format(upper_right)
            raise ValueError(msg)

        if len(upper_right) != 2 and len(upper_right) != 3:
            msg = 'Unable to set CMFD Mesh with upper_right {0} since it ' \
                  'must include 2 or 3 dimensions'.format(upper_right)
            raise ValueError(msg)

        for coord in upper_right:
            if not isinstance(coord, Real):
                msg = 'Unable to set CMFD Mesh with upper_right {0} which ' \
                      'is not a real number'.format(coord)
                raise ValueError(msg)

        self._upper_right = upper_right

    @dimension.setter
    def dimension(self, dimension):
        if not isinstance(dimension, Sequence):
            msg = 'Unable to set CMFD Mesh with dimension {0} which is ' \
                  'not a sequence'.format(dimension)
            raise ValueError(msg)

        elif len(dimension) != 2 and len(dimension) != 3:
            msg = 'Unable to set CMFD Mesh with dimension {0} since it ' \
                  'must include 2 or 3 dimensions'.format(dimension)
            raise ValueError(msg)

        for dim in dimension:
            if not isinstance(dim, Integral):
                msg = 'Unable to set CMFD Mesh with dimension {0} which ' \
                      'is a non-integer'.format(dim)
                raise ValueError(msg)

        self._dimension = dimension

    @width.setter
    def width(self, width):
        if width is not None:
            if not isinstance(width, Sequence):
                msg = 'Unable to set CMFD Mesh with width {0} which ' \
                      'is not a sequence'.format(width)
                raise ValueError(msg)

        if len(width) != 2 and len(width) != 3:
            msg = 'Unable to set CMFD Mesh with width {0} since it must ' \
                  'include 2 or 3 dimensions'.format(width)
            raise ValueError(msg)

        for dim in width:
            if not isinstance(dim, Real):
                msg = 'Unable to set CMFD Mesh with width {0} which is ' \
                      'not a real number'.format(width)
                raise ValueError(msg)

        self._width = width

    @energy.setter
    def energy(self, energy):
        if not isinstance(energy, Sequence):
            msg = 'Unable to set CMFD Mesh energy to {0} which is not ' \
                  'a sequence'.format(energy)
            raise ValueError(msg)

        for e in energy:
            if not isinstance(e, Real):
                msg = 'Unable to set CMFD Mesh energy to {0} which is not ' \
                      'a real number'.format(e)
                raise ValueError(msg)
            elif e < 0:
                msg = 'Unable to set CMFD Mesh energy to {0} which is ' \
                      'is a negative integer'.format(e)
                raise ValueError(msg)

        self._energy = energy

    @albedo.setter
    def albedo(self, albedo):
        if not isinstance(albedo, Sequence):
            msg = 'Unable to set CMFD Mesh albedo to {0} which is not ' \
                  'a sequence'.format(albedo)
            raise ValueError(msg)

        if not len(albedo) == 6:
            msg = 'Unable to set CMFD Mesh albedo to {0} which is not ' \
                  'length 6 for +/-x,y,z'.format(albedo)
            raise ValueError(msg)

        for a in albedo:
            if not isinstance(a, Real):
                msg = 'Unable to set CMFD Mesh albedo to {0} which is not ' \
                      'a real number'.format(a)
                raise ValueError(msg)
            elif a < 0 or a > 1:
                msg = 'Unable to set CMFD Mesh albedo to {0} which is ' \
                      'is not in [0,1]'.format(a)
                raise ValueError(msg)

        self._albedo = albedo

    @map.setter
    def map(self, map):

        if not isinstance(map, Sequence):
            msg = 'Unable to set CMFD Mesh map to {0} which is not ' \
                  'a sequence'.format(map)
            raise ValueError(msg)

        for m in map:
            if m != 1 and m != 2:
                msg = 'Unable to set CMFD Mesh map to {0} which is ' \
                      'is not 1 or 2'.format(m)
                raise ValueError(msg)

        self._map = map

    def _get_xml_element(self):
        element = ET.Element("mesh")

        subelement = ET.SubElement(element, "lower_left")
        subelement.text = ' '.join(map(str, self._lower_left))

        if self.upper_right is not None:
            subelement = ET.SubElement(element, "upper_right")
            subelement.text = ' '.join(map(str, self.upper_right))

        subelement = ET.SubElement(element, "dimension")
        subelement.text = ' '.join(map(str, self.dimension))

        if self.width is not None:
            subelement = ET.SubElement(element, "width")
            subelement.text = ' '.join(map(str, self.width))

        if self.energy is not None:
            subelement = ET.SubElement(element, "energy")
            subelement.text = ' '.join(map(str, self.energy))

        if self.albedo is not None:
            subelement = ET.SubElement(element, "albedo")
            subelement.text = ' '.join(map(str, self.albedo))

        if self.map is not None:
            subelement = ET.SubElement(element, "map")
            subelement.text = ' '.join(map(str, self.map))

        return element


class CMFDFile(object):
    """Parameters that control the use of coarse-mesh finite difference acceleration
    in OpenMC. This corresponds directly to the cmfd.xml input file.

    Attributes
    ----------
    begin : int
        Batch number at which CMFD calculations should begin
    dhat_reset : bool
        Indicate whether :math:`\widehat{D}` nonlinear CMFD parameters should be
        reset to zero before solving CMFD eigenproblem.
    display : {'balance', 'dominance', 'entropy', 'source'}
        Set one additional CMFD output column. Options are:

          * "balance" - prints the RMS [%] of the resdiual from the neutron balance
             equation on CMFD tallies.
          * "dominance" - prints the estimated dominance ratio from the CMFD
            iterations.
          * "entropy" - prints the *entropy* of the CMFD predicted fission source.
          * "source" - prints the RMS [%] between the OpenMC fission source and
            CMFD fission source.
    downscatter : bool
        Indicate whether an effective downscatter cross section should be used
        when using 2-group CMFD.
    feedback : bool
        Indicate or not the CMFD diffusion result is used to adjust the weight
        of fission source neutrons on the next OpenMC batch. Defaults to False.
    gauss_seidel_tolerance : Sequence of float
        Two parameters specifying the absolute inner tolerance and the relative
        inner tolerance for Gauss-Seidel iterations when performing CMFD.
    ktol : float
        Tolerance on the eigenvalue when performing CMFD power iteration
    cmfd_mesh : CMFDMesh
        Structured mesh to be used for acceleration
    norm : float
        Normalization factor applied to the CMFD fission source distribution
    power_monitor : bool
        View convergence of power iteration during CMFD acceleration
    run_adjoint : bool
        Perform adjoint calculation on the last batch
    shift : float
        Optional Wielandt shift parameter for accelerating power iterations. By
        default, it is very large so there is effectively no impact.
    spectral : float
        Optional spectral radius that can be used to accelerate the convergence
        of Gauss-Seidel iterations during CMFD power iteration.
    stol : float
        Tolerance on the fission source when performing CMFD power iteration
    tally_reset : list of int
        List of batch numbers at which CMFD tallies should be reset
    write_matrices : bool
        Write sparse matrices that are used during CMFD acceleration (loss,
        production) to file

    """

    def __init__(self):
        self._begin = None
        self._dhat_reset = None
        self._display = None
        self._downscatter = None
        self._feedback = None
        self._gauss_seidel_tolerance = None
        self._ktol = None
        self._cmfd_mesh = None
        self._norm = None
        self._power_monitor = None
        self._run_adjoint = None
        self._shift = None
        self._spectral = None
        self._stol = None
        self._tally_reset = None
        self._write_matrices = None

        self._cmfd_file = ET.Element("cmfd")
        self._cmfd_mesh_element = None

    @property
    def begin(self):
        return self._begin

    @property
    def dhat_reset(self):
        return self._dhat_reset

    @property
    def display(self):
        return self._display

    @property
    def downscatter(self):
        return self._downscatter

    @property
    def feedback(self):
        return self._feedback

    @property
    def gauss_seidel_tolerance(self):
        return self._gauss_seidel_tolerance

    @property
    def ktol(self):
        return self._ktol

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
    def shift(self):
        return self._shift

    @property
    def spectral(self):
        return self._spectral

    @property
    def stol(self):
        return self._stol

    @property
    def tally_reset(self):
        return self._tally_reset

    @property
    def write_matrices(self):
        return self._write_matrices

    @begin.setter
    def begin(self, begin):
        if not isinstance(begin, Integral):
            msg = 'Unable to set CMFD begin batch to a non-integer ' \
                  'value {0}'.format(begin)
            raise ValueError(msg)

        if begin <= 0:
            msg = 'Unable to set CMFD begin batch batch to a negative ' \
                  'value {0}'.format(begin)
            raise ValueError(msg)

        self._begin = begin

    @dhat_reset.setter
    def dhat_reset(self, dhat_reset):
        if not isinstance(dhat_reset, bool):
            msg = 'Unable to set Dhat reset to {0} which is ' \
                  'a non-boolean value'.format(dhat_reset)
            raise ValueError(msg)

        self._dhat_reset = dhat_reset

    @display.setter
    def display(self, display):
        if not isinstance(basestring):
            msg = 'Unable to set CMFD display to a non-string ' \
                  'value'.format(display)
            raise ValueError(msg)

        if display not in ['balance', 'dominance', 'entropy', 'source']:
            msg = 'Unable to set CMFD display to {0} which is ' \
                  'not an accepted value'.format(display)
            raise ValueError(msg)

        self._display = display

    @downscatter.setter
    def downscatter(self, downscatter):
        if not isinstance(downscatter, bool):
            msg = 'Unable to set downscatter to {0} which is ' \
                  'a non-boolean value'.format(downscatter)
            raise ValueError(msg)

        self._downscatter = downscatter

    @feedback.setter
    def feedback(self, feedback):
        if not isinstance(feedback, bool):
            msg = 'Unable to set CMFD feedback to {0} which is ' \
                  'a non-boolean value'.format(feedback)
            raise ValueError(msg)

        self._feedback = feedback

    @gauss_seidel_tolerance.setter
    def gauss_seidel_tolerance(self, gauss_seidel_tolerance):
        if not isinstance(gauss_seidel_tolerance, Sequence):
            msg = 'Unable to set Gauss-Seidel tolerance to {0} which is ' \
                  'not a sequence'.format(gauss_seidel_tolerance)
            raise ValueError(msg)

        if len(gauss_seidel_tolerance) != 2:
            msg = 'Unable to set Gauss-Seidel tolerance with {0} since ' \
                  'it must be of length 2'.format(width)
            raise ValueError(msg)

        for t in gauss_seidel_tolerance:
            if not isinstance(t, Real):
                msg = 'Unable to set Gauss-Seidel tolerance with {0} which ' \
                      'is not a real number'.format(t)
                raise ValueError(msg)

        self._gauss_seidel_tolerance = gauss_seidel_tolerance

    @ktol.setter
    def ktol(self, ktol):
        if not isinstance(ktol, Real):
            msg = 'Unable to set the eigenvalue tolerance to {0} which is ' \
                  'not a real number'.format(ktol)
            raise ValueError(msg)

        self._ktol = ktol

    @cmfd_mesh.setter
    def cmfd_mesh(self, mesh):
        if not isinstance(mesh, CMFDMesh):
            msg = 'Unable to set CMFD mesh to {0} which is not a ' \
                  'CMFDMesh object'.format(mesh)
            raise ValueError(msg)

        self._mesh = mesh

    @norm.setter
    def norm(self, norm):
        if not isinstance(norm, Real):
            msg = 'Unable to set the CMFD norm to {0} which is not ' \
                  'a real number'.format(norm)
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

    @shift.setter
    def shift(self, shift):
        if not isinstance(shift, Real):
            msg = 'Unable to set the Wielandt shift to {0} which is ' \
                  'not a real number'.format(shift)
            raise ValueError(msg)

        self._shift = shift

    @spectral.setter
    def spectral(self, spectral):
        if not isinstance(spectral, Real):
            msg = 'Unable to set the spectral radius to {0} which is ' \
                  'not a real number'.format(spectral)
            raise ValueError(msg)

        self._spectral = spectral

    @stol.setter
    def stol(self, stol):
        if not isinstance(stol, Real):
            msg = 'Unable to set the fission source tolerance to {0} which ' \
                  'is not a real number'.format(stol)
            raise ValueError(msg)

        self._stol = stol

    @tally_reset.setter
    def tally_reset(self, tally_reset):
        if not isinstance(tally_reset, Sequence):
            msg = 'Unable to set tally reset batches to {0} which is ' \
                  'not a sequence'.format(tally_reset)
            raise ValueError(msg)

        for t in tally_reset:
            if not isinstance(t, Integral):
                msg = 'Unable to set tally reset batch to {0} which ' \
                      'is not an integer'.format(t)
                raise ValueError(msg)

        self._tally_reset = tally_reset

    @write_matrices.setter
    def write_matrices(self, write_matrices):
        if not isinstance(write_matrices, bool):
            msg = 'Unable to set CMFD write matrices to {0} which is a ' \
                  'non-boolean value'.format(write_matrices)
            raise ValueError(msg)

        self._write_matrices = write_matrices

    def _create_begin_subelement(self):
        if self._begin is not None:
            element = ET.SubElement(self._cmfd_file, "begin")
            element.text = str(self._begin)

    def _create_dhat_reset_subelement(self):
        if self._dhat_reset is not None:
            element = ET.SubElement(self._cmfd_file, "dhat_reset")
            element.text = str(self._dhat_reset).lower()

    def _create_display_subelement(self):
        if self._display is not None:
            element = ET.SubElement(self._cmfd_file, "display")
            element.text = str(self._display)

    def _create_downscatter_subelement(self):
        if self._downscatter is not None:
            element = ET.SubElement(self._cmfd_file, "downscatter")
            element.text = str(self._downscatter).lower()

    def _create_feedback_subelement(self):
        if self._feedback is not None:
            element = ET.SubElement(self._cmfd_file, "feeback")
            element.text = str(self._feedback).lower()

    def _create_gauss_seidel_tolerance_subelement(self):
        if self._gauss_seidel_tolerance is not None:
            element = ET.SubElement(self._cmfd_file, "gauss_seidel_tolerance")
            element.text = ' '.join(map(str, self._gauss_seidel_tolerance))

    def _create_ktol_subelement(self):
        if self._ktol is not None:
            element = ET.SubElement(self._ktol, "ktol")
            element.text = str(self._ktol)

    def _create_mesh_subelement(self):
        if self._mesh is not None:
            xml_element = self._mesh._get_xml_element()
            self._cmfd_file.append(xml_element)

    def _create_norm_subelement(self):
        if self._norm is not None:
            element = ET.SubElement(self._cmfd_file, "norm")
            element.text = str(self._norm)

    def _create_power_monitor_subelement(self):
        if self._power_monitor is not None:
            element = ET.SubElement(self._cmfd_file, "power_monitor")
            element.text = str(self._power_monitor).lower()

    def _create_run_adjoint_subelement(self):
        if self._run_adjoint is not None:
            element = ET.SubElement(self._cmfd_file, "run_adjoint")
            element.text = str(self._run_adjoint).lower()

    def _create_shift_subelement(self):
        if self._shift is not None:
            element = ET.SubElement(self._shift, "shift")
            element.text = str(self._shift)

    def _create_spectral_subelement(self):
        if self._spectral is not None:
            element = ET.SubElement(self._spectral, "spectral")
            element.text = str(self._spectral)

    def _create_stol_subelement(self):
        if self._stol is not None:
            element = ET.SubElement(self._stol, "stol")
            element.text = str(self._stol)

    def _create_tally_reset_subelement(self):
        if self._tally_reset is not None:
            element = ET.SubElement(self._tally_reset, "tally_reset")
            element.text = ' '.join(map(str, self._tally_reset))

    def _create_write_matrices_subelement(self):
        if self._write_matrices is not None:
            element = ET.SubElement(self._cmfd_file, "write_matrices")
            element.text = str(self._write_matrices).lower()

    def export_to_xml(self):
        """Create a cmfd.xml file using the class data that can be used for an OpenMC
        simulation.

        """

        self._create_begin_subelement()
        self._create_dhat_reset_subelement()
        self._create_display_subelement()
        self._create_downscatter_subelement()
        self._create_feedback_subelement()
        self._create_gauss_seidel_tolerance_subelement()
        self._create_ktol_subelement()
        self._create_mesh_subelement()
        self._create_norm_subelement()
        self._create_power_monitor_subelement()
        self._create_run_adjoint_subelement()
        self._create_shift_subelement()
        self._create_spectral_subelement()
        self._create_stol_subelement()
        self._create_tally_reset_subelement()
        self._create_write_matrices_subelement()

        # Clean the indentation in the file to be user-readable
        clean_xml_indentation(self._cmfd_file)

        # Write the XML Tree to the cmfd.xml file
        tree = ET.ElementTree(self._cmfd_file)
        tree.write("cmfd.xml", xml_declaration=True,
                   encoding='utf-8', method="xml")
