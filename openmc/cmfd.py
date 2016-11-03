"""This module can be used to specify parameters used for coarse mesh finite
difference (CMFD) acceleration in OpenMC. CMFD was first proposed by [Smith]_
and is widely used in accelerating neutron transport problems.

References
----------

.. [Smith] K. Smith, "Nodal method storage reduction by non-linear
   iteration", *Trans. Am. Nucl. Soc.*, **44**, 265 (1983).

"""

from collections import Iterable
from numbers import Real, Integral
from xml.etree import ElementTree as ET
import sys

from six import string_types

from openmc.clean_xml import clean_xml_indentation
from openmc.checkvalue import (check_type, check_length, check_value,
                               check_greater_than, check_less_than)


class CMFDMesh(object):
    """A structured Cartesian mesh used for Coarse Mesh Finite Difference (CMFD)
    acceleration.

    Attributes
    ----------
    lower_left : Iterable of float
        The lower-left corner of the structured mesh. If only two coordinates are
        given, it is assumed that the mesh is an x-y mesh.
    upper_right : Iterable of float
        The upper-right corner of the structrued mesh. If only two coordinates
        are given, it is assumed that the mesh is an x-y mesh.
    dimension : Iterable of int
        The number of mesh cells in each direction.
    width : Iterable of float
        The width of mesh cells in each direction.
    energy : Iterable of float
        Energy bins in eV, listed in ascending order (e.g. [0.0, 0.625e-1,
        20.0e6]) for CMFD tallies and acceleration. If no energy bins are listed,
        OpenMC automatically assumes a one energy group calculation over the
        entire energy range.
    albedo : Iterable of float
        Surface ratio of incoming to outgoing partial currents on global
        boundary conditions. They are listed in the following order: -x +x -y +y
        -z +z.
    map : Iterable of int
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
        check_type('CMFD mesh lower_left', lower_left, Iterable, Real)
        check_length('CMFD mesh lower_left', lower_left, 2, 3)
        self._lower_left = lower_left

    @upper_right.setter
    def upper_right(self, upper_right):
        check_type('CMFD mesh upper_right', upper_right, Iterable, Real)
        check_length('CMFD mesh upper_right', upper_right, 2, 3)
        self._upper_right = upper_right

    @dimension.setter
    def dimension(self, dimension):
        check_type('CMFD mesh dimension', dimension, Iterable, Integral)
        check_length('CMFD mesh dimension', dimension, 2, 3)
        self._dimension = dimension

    @width.setter
    def width(self, width):
        check_type('CMFD mesh width', width, Iterable, Real)
        check_length('CMFD mesh width', width, 2, 3)
        self._width = width

    @energy.setter
    def energy(self, energy):
        check_type('CMFD mesh energy', energy, Iterable, Real)
        for e in energy:
            check_greater_than('CMFD mesh energy', e, 0, True)
        self._energy = energy

    @albedo.setter
    def albedo(self, albedo):
        check_type('CMFD mesh albedo', albedo, Iterable, Real)
        check_length('CMFD mesh albedo', albedo, 6)
        for a in albedo:
            check_greater_than('CMFD mesh albedo', a, 0, True)
            check_less_than('CMFD mesh albedo', a, 1, True)
        self._albedo = albedo

    @map.setter
    def map(self, meshmap):
        check_type('CMFD mesh map', meshmap, Iterable, Integral)
        for m in meshmap:
            check_value('CMFD mesh map', m, [1, 2])
        self._map = meshmap

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


class CMFD(object):
    r"""Parameters that control the use of coarse-mesh finite difference acceleration
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
    gauss_seidel_tolerance : Iterable of float
        Two parameters specifying the absolute inner tolerance and the relative
        inner tolerance for Gauss-Seidel iterations when performing CMFD.
    ktol : float
        Tolerance on the eigenvalue when performing CMFD power iteration
    cmfd_mesh : openmc.CMFDMesh
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
        check_type('CMFD begin batch', begin, Integral)
        check_greater_than('CMFD begin batch', begin, 0)
        self._begin = begin

    @dhat_reset.setter
    def dhat_reset(self, dhat_reset):
        check_type('CMFD Dhat reset', dhat_reset, bool)
        self._dhat_reset = dhat_reset

    @display.setter
    def display(self, display):
        check_type('CMFD display', display, string_types)
        check_value('CMFD display', display,
                    ['balance', 'dominance', 'entropy', 'source'])
        self._display = display

    @downscatter.setter
    def downscatter(self, downscatter):
        check_type('CMFD downscatter', downscatter, bool)
        self._downscatter = downscatter

    @feedback.setter
    def feedback(self, feedback):
        check_type('CMFD feedback', feedback, bool)
        self._feedback = feedback

    @gauss_seidel_tolerance.setter
    def gauss_seidel_tolerance(self, gauss_seidel_tolerance):
        check_type('CMFD Gauss-Seidel tolerance', gauss_seidel_tolerance,
                   Iterable, Real)
        check_length('Gauss-Seidel tolerance', gauss_seidel_tolerance, 2)
        self._gauss_seidel_tolerance = gauss_seidel_tolerance

    @ktol.setter
    def ktol(self, ktol):
        check_type('CMFD eigenvalue tolerance', ktol, Real)
        self._ktol = ktol

    @cmfd_mesh.setter
    def cmfd_mesh(self, mesh):
        check_type('CMFD mesh', mesh, CMFDMesh)
        self._mesh = mesh

    @norm.setter
    def norm(self, norm):
        check_type('CMFD norm', norm, Real)
        self._norm = norm

    @power_monitor.setter
    def power_monitor(self, power_monitor):
        check_type('CMFD power monitor', power_monitor, bool)
        self._power_monitor = power_monitor

    @run_adjoint.setter
    def run_adjoint(self, run_adjoint):
        check_type('CMFD run adjoint', run_adjoint, bool)
        self._run_adjoint = run_adjoint

    @shift.setter
    def shift(self, shift):
        check_type('CMFD Wielandt shift', shift, Real)
        self._shift = shift

    @spectral.setter
    def spectral(self, spectral):
        check_type('CMFD spectral radius', spectral, Real)
        self._spectral = spectral

    @stol.setter
    def stol(self, stol):
        check_type('CMFD fission source tolerance', stol, Real)
        self._stol = stol

    @tally_reset.setter
    def tally_reset(self, tally_reset):
        check_type('tally reset batches', tally_reset, Iterable, Integral)
        self._tally_reset = tally_reset

    @write_matrices.setter
    def write_matrices(self, write_matrices):
        check_type('CMFD write matrices', write_matrices, bool)
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
