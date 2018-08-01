"""This module can be used to specify parameters used for coarse mesh finite
difference (CMFD) acceleration in OpenMC. CMFD was first proposed by [Smith]_
and is widely used in accelerating neutron transport problems.

References
----------

.. [Smith] K. Smith, "Nodal method storage reduction by non-linear
   iteration", *Trans. Am. Nucl. Soc.*, **44**, 265 (1983).

"""

from collections.abc import Iterable
from numbers import Real, Integral
from xml.etree import ElementTree as ET
import sys
import numpy as np
import openmc.capi
from openmc.exceptions import AllocationError

from openmc.clean_xml import clean_xml_indentation
from openmc.checkvalue import (check_type, check_length, check_value,
                               check_greater_than, check_less_than)

# Maximum/minimum neutron energies, from src/api.F90
ENERGY_MAX_NEUTRON = np.inf
ENERGY_MIN_NEUTRON = 0.

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
        for d in dimension:
            check_greater_than('CMFD mesh dimension', d, 0)
        self._dimension = dimension

    @width.setter
    def width(self, width):
        check_type('CMFD mesh width', width, Iterable, Real)
        check_length('CMFD mesh width', width, 2, 3)
        for w in width:
            check_greater_than('CMFD mesh width', w, 0)
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
        check_type('CMFD display', display, str)
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
        self._cmfd_mesh = mesh

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
        if self._cmfd_mesh is not None:
            xml_element = self._cmfd_mesh._get_xml_element()
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

class CMFDRun(object):
    r"""Class to run openmc with CMFD acceleration through the C API. Running
    openmc through this manner obviates the need of defining CMFD parameters
    through a cmfd.xml file. Instead, all input parameters should be passed through
    the CMFDRun initializer.

    Self notes:
    Have it so that only required parameters part of __init__, all other parameters
        need to be set through setter functions
    Make sure all input variables are used / have coveragage somewhere in code


    Input parameters: Fill out, specifying which ones are required for CMFDRun
        Required: cmfd_mesh, cmfd_mesh.dimension, cmfd_mesh.lower_left

    Ignore:
        1) MG mode - read_cmfd_xml
        2) configure_cmfd
            call time_cmfd % reset()
            call time_cmfdbuild % reset()
            call time_cmfdsolve % reset()
        3) cmfd_adjoint_type input parameter, either set to 'physical' or 'math' - see cmfd_solver.F90



    Attributes
    ----------
    To add: cmfd_coremap: Flag for active core map
            indices: Stores spatial and group dimensions as [nx, ny, nz, ng]
            egrid: energy grid used for CMFD acceleration
            albedo: Albedo for global boundary conditions, taken from CMFD mesh. Set to [1,1,1,1,1,1] if not specified by user
            cmfd_coremap: Optional acceleration map to overlay on coarse mesh spatial grid, taken from CMFD mesh
            n_cmfd_resets: Number of elements in tally_reset, list that stores batches where CMFD tallies should be reset
            cmfd_atoli: Absolute GS tolerance, set by gauss_seidel_tolerance
            cmfd_rtoli: Relative GS tolerance, set by gauss_seidel_tolerance
            cmfd_mesh_id: Mesh id of openmc.capi.Mesh object that corresponds to the CMFD mesh
            energy_filters: Boolean that stores whether energy filters should be created or not.
                Set to true if user specifies energy grid in CMFDMesh, false otherwise

    TODO: Put descriptions for all methods in CMFDRun

    """

    def __init__(self):

        # Set CMFD default parameters based on cmfd_header.F90
        # Input parameters that user can define
        self._cmfd_begin = 1
        self._dhat_reset = False
        self._cmfd_display = 'balance'
        self._cmfd_downscatter = False
        self._cmfd_feedback = False
        self._gauss_seidel_tolerance = [1.e-10, 1.e-5]
        self._cmfd_ktol = 1.e-8
        self._cmfd_mesh = None
        self._norm = 1.
        self._cmfd_power_monitor = False
        self._cmfd_run_adjoint = False
        self._cmfd_shift = 1.e-6
        self._cmfd_spectral = 0.
        self._cmfd_stol = 1.e-8
        self._cmfd_reset = None
        self._cmfd_write_matrices = False

        # External variables used during runtime but users don't have control over
        self._cmfd_coremap = False
        self._indices = np.zeros(4, dtype=int)
        self._egrid = None
        self._albedo = None
        self._coremap = None
        self._n_cmfd_resets = 0
        self._cmfd_atoli = None
        self._cmfd_rtoli = None
        self._cmfd_mesh_id = None
        self._energy_filters = None

        # TODO All timing variables



    @property
    def cmfd_begin(self):
        return self._cmfd_begin

    @property
    def dhat_reset(self):
        return self._dhat_reset

    @property
    def display(self):
        return self._display

    @property
    def cmfd_downscatter(self):
        return self._cmfd_downscatter

    @property
    def cmfd_feedback(self):
        return self._cmfd_feedback

    @property
    def gauss_seidel_tolerance(self):
        return self._gauss_seidel_tolerance

    @property
    def cmfd_ktol(self):
        return self._cmfd_ktol

    @property
    def cmfd_mesh(self):
        return self._cmfd_mesh

    @property
    def norm(self):
        return self._norm

    @property
    def cmfd_power_monitor(self):
        return self._cmfd_power_monitor

    @property
    def cmfd_run_adjoint(self):
        return self._cmfd_run_adjoint

    @property
    def cmfd_shift(self):
        return self._cmfd_shift

    @property
    def cmfd_spectral(self):
        return self._cmfd_spectral

    @property
    def cmfd_stol(self):
        return self._cmfd_stol

    @property
    def cmfd_reset(self):
        return self._cmfd_reset

    @property
    def cmfd_write_matrices(self):
        return self._cmfd_write_matrices

    @cmfd_begin.setter
    def cmfd_begin(self, cmfd_begin):
        check_type('CMFD begin batch', cmfd_begin, Integral)
        check_greater_than('CMFD begin batch', cmfd_begin, 0)
        self._cmfd_begin = begin

    @dhat_reset.setter
    def dhat_reset(self, dhat_reset):
        check_type('CMFD Dhat reset', dhat_reset, bool)
        self._dhat_reset = dhat_reset

    @display.setter
    def display(self, display):
        check_type('CMFD display', display, str)
        check_value('CMFD display', display,
                    ['balance', 'dominance', 'entropy', 'source'])
        self._display = display

    @cmfd_downscatter.setter
    def cmfd_downscatter(self, cmfd_downscatter):
        check_type('CMFD downscatter', cmfd_downscatter, bool)
        self._cmfd_downscatter = cmfd_downscatter

    @cmfd_feedback.setter
    def cmfd_feedback(self, cmfd_feedback):
        check_type('CMFD feedback', cmfd_feedback, bool)
        self._cmfd_feedback = cmfd_feedback

    @gauss_seidel_tolerance.setter
    def gauss_seidel_tolerance(self, gauss_seidel_tolerance):
        check_type('CMFD Gauss-Seidel tolerance', gauss_seidel_tolerance,
                   Iterable, Real)
        check_length('Gauss-Seidel tolerance', gauss_seidel_tolerance, 2)
        self._gauss_seidel_tolerance = gauss_seidel_tolerance

    @cmfd_ktol.setter
    def cmfd_ktol(self, cmfd_ktol):
        check_type('CMFD eigenvalue tolerance', cmfd_ktol, Real)
        self._cmfd_ktol = cmfd_ktol

    @cmfd_mesh.setter
    def cmfd_mesh(self, mesh):
        check_type('CMFD mesh', mesh, CMFDMesh)

        # Check dimension defined
        if mesh.dimension is None:
            raise AllocationError('CMFD mesh requires spatial '
                                  'dimensions to be specified')

        # Check lower left defined
        if mesh.lower_left is None:
            raise AllocationError('CMFD mesh requires lower left coordinates '
                                  'to be specified')

        # Check that both upper right and width both not defined
        if mesh.upper_right is not None and mesh.width is not None:
            raise AllocationError('Both upper right coordinates and width '
                                  'cannot be specified for CMFD mesh')

        # Check that at least one of width or upper right is defined
        if mesh.upper_right is None and mesh.width is None:
            raise AllocationError('CMFD mesh requires either upper right '
                                  'coordinates or width to be specified')

        # Check width and lower length are same dimension and define upper_right
        if mesh.width is not None:
            check_length('CMFD mesh width', mesh.width, len(mesh.lower_left))
            mesh.upper_right = np.array(mesh.lower_left) + \
                               np.array(mesh.width) * np.array(mesh.dimension)

        # Check upper_right and lower length are same dimension and define width
        elif mesh.upper_right is not None:
            check_length('CMFD mesh upper right', mesh.upper_right, \
                         len(mesh.lower_left))
            # Check upper right coordinates are greater than lower left
            if np.any(np.array(mesh.upper_right) <= np.array(mesh.lower_left)):
                raise AllocationError('CMFD mesh requires upper right '
                                      'coordinates to be greater than lower '
                                      'left coordinates')
            mesh.width = np.true_divide(
                         (np.array(mesh.upper_right) - np.array(mesh.lower_left)), \
                         np.array(mesh.dimension))
        self._cmfd_mesh = mesh

    @norm.setter
    def norm(self, norm):
        check_type('CMFD norm', norm, Real)
        self._norm = norm

    @cmfd_power_monitor.setter
    def cmfd_power_monitor(self, cmfd_power_monitor):
        check_type('CMFD power monitor', cmfd_power_monitor, bool)
        self._cmfd_power_monitor = cmfd_power_monitor

    @cmfd_run_adjoint.setter
    def cmfd_run_adjoint(self, cmfd_run_adjoint):
        check_type('CMFD run adjoint', cmfd_run_adjoint, bool)
        self._cmfd_run_adjoint = cmfd_run_adjoint

    @cmfd_shift.setter
    def cmfd_shift(self, cmfd_shift):
        check_type('CMFD Wielandt shift', cmfd_shift, Real)
        self._cmfd_shift = cmfd_shift

    @cmfd_spectral.setter
    def cmfd_spectral(self, cmfd_spectral):
        check_type('CMFD spectral radius', cmfd_spectral, Real)
        self._cmfd_spectral = cmfd_spectral

    @cmfd_stol.setter
    def cmfd_stol(self, cmfd_stol):
        check_type('CMFD fission source tolerance', cmfd_stol, Real)
        self._cmfd_stol = cmfd_stol

    @cmfd_reset.setter
    def cmfd_reset(self, cmfd_reset):
        check_type('tally reset batches', cmfd_reset, Iterable, Integral)
        self._cmfd_reset = cmfd_reset

    @cmfd_write_matrices.setter
    def cmfd_write_matrices(self, cmfd_write_matrices):
        check_type('CMFD write matrices', cmfd_write_matrices, bool)
        self._cmfd_write_matrices = cmfd_write_matrices

    def run(self):
        openmc.capi.init()
        self._configure_cmfd()
        openmc.capi.simulation_init()
        openmc.capi.run()
        openmc.capi.simulation_finalize()
        openmc.capi.finalize()

    def _configure_cmfd(self):
        # TODO Tell user
        # print('Configuring CMFD parameters for simulation')

        # Check if CMFD mesh is defined
        if self._cmfd_mesh is None:
            raise AllocationError('No CMFD mesh has been specified for '
                                  'simulation')

        # Set spatial dimensions of CMFD object
        # Iterate through each element of self._cmfd_mesh.dimension as could be
        # length 2 or 3
        for i, n in enumerate(self._cmfd_mesh.dimension):
            self._indices[i] = n

        # Set number of energy groups
        if self._cmfd_mesh.energy is not None:
            ng = len(self._cmfd_mesh.energy)
            self._egrid = np.array(self._cmfd_mesh.energy)
            self._indices[3] = ng - 1
            self._energy_filters = True
            # TODO: MG mode check
        else:
            self._egrid = np.array([ENERGY_MIN_NEUTRON, ENERGY_MAX_NEUTRON])
            self._indices[3] = 1
            self._energy_filters = False

        # Set global albedo
        if self._cmfd_mesh.albedo is not None:
            self._albedo = np.array(self._cmfd_mesh.albedo)
        else:
            self._albedo = np.array([1.,1.,1.,1.,1.,1.])

        # Get acceleration map
        if self._cmfd_mesh.map is not None:
            check_length('CMFD coremap', self._cmfd_mesh.map,
                         np.product(self._indices[0:3]))
            self._coremap = np.array(self._cmfd_mesh.map).reshape(( \
                            self._indices[0], self._indices[1], \
                            self._indices[2]))
            self._cmfd_coremap = True

        # Set number of batches where cmfd tallies should be reset
        if self._cmfd_reset is not None:
            self._n_cmfd_resets = len(self._cmfd_reset)

        # Set Gauss Sidel tolerances
        self._cmfd_atoli = self._gauss_seidel_tolerance[0]
        self._cmfd_rtoli = self._gauss_seidel_tolerance[1]

        # Create tally objects
        self._create_cmfd_tally()

    def _create_cmfd_tally(self):
        # Create Mesh object based on CMFDMesh, stored internally
        cmfd_mesh = openmc.capi.Mesh()
        # Store id of Mesh object
        self._cmfd_mesh_id = cmfd_mesh.id
        # Set dimension and parameters of Mesh object
        cmfd_mesh.dimension = self._cmfd_mesh.dimension
        cmfd_mesh.set_parameters(lower_left=self._cmfd_mesh.lower_left,
                                 upper_right=self._cmfd_mesh.upper_right,
                                 width=self._cmfd_mesh.width)

        # Create Mesh Filter object, stored internally
        mesh_filter = openmc.capi.MeshFilter()
        # Set mesh for Mesh Filter
        mesh_filter.mesh = cmfd_mesh

        # Set up energy filters, if applicable
        if self._energy_filters:
            # Create Energy Filter object, stored internally
            energy_filter = openmc.capi.EnergyFilter()
            # Set bins for Energy Filter
            energy_filter.bins = self._egrid

            # Create Energy Out Filter object, stored internally
            energyout_filter = openmc.capi.EnergyoutFilter()
            # Set bins for Energy Filter
            energyout_filter.bins = self._egrid

        # Create Mesh Surface Filter object, stored internally
        meshsurface_filter = openmc.capi.MeshSurfaceFilter()
        # Set mesh for Mesh Surface Filter
        meshsurface_filter.mesh = cmfd_mesh

        # Create Legendre Filter object, stored internally
        legendre_filter = openmc.capi.LegendreFilter()
        # Set order for Legendre Filter
        legendre_filter.order = 1

        # Create CMFD tallies, stored internally
        n_tallies = 4
        for i in range(n_tallies):
            tally = openmc.capi.Tally()
            # Set nuclide bins

            # Set attributes of CMFD flux, total tally
            if i == 0:
                # Set filters for tally
                if self._energy_filters:
                    tally.filters = [mesh_filter, energy_filter]
                else:
                    tally.filters = [mesh_filter]
                # Set scores for tally
                tally.scores = ['flux', 'total']

            # Set attributes of CMFD neutron production tally
            if i == 1:
                # Set filters for tally
                if self._energy_filters:
                    tally.filters = [mesh_filter, energy_filter, energyout_filter]
                else:
                    tally.filters = [mesh_filter]
                # Set scores for tally
                tally.scores = ['nu-scatter', 'nu-fission']

            # Set attributes of CMFD surface current tally
            if i == 2:
                # Set filters for tally
                if self._energy_filters:
                    tally.filters = [meshsurface_filter, energy_filter]
                else:
                    tally.filters = [meshsurface_filter]
                # Set scores for tally
                tally.scores = ['current']

            # Set attributes of CMFD P1 scatter tally
            if i == 3:
                # Set filters for tally
                if self._energy_filters:
                    tally.filters = [mesh_filter, legendre_filter, energy_filter]
                else:
                    tally.filters = [mesh_filter, legendre_filter]
                # Set scores for tally
                tally.scores = ['scatter']

            # Set all tallies to be active from beginning
            tally.active = True