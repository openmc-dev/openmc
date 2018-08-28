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
from xml.etree import ElementTree as ET # TODO Remove
import sys   # TODO Remove
import numpy as np
# TODO Include comment for this
np.seterr(divide='ignore', invalid='ignore')
from scipy import sparse
import time
from mpi4py import MPI

import openmc.capi
from openmc.clean_xml import clean_xml_indentation  # TODO Remove
from openmc.checkvalue import (check_type, check_length, check_value,
                               check_greater_than, check_less_than)
from openmc.exceptions import OpenMCError



"""
--------------
CMFD CONSTANTS
--------------
"""
# Maximum/minimum neutron energies
_ENERGY_MAX_NEUTRON = np.inf
_ENERGY_MIN_NEUTRON = 0.

# Tolerance for detecting zero flux values
_TINY_BIT = 1.e-8

# For non-accelerated regions on coarse mesh overlay
_CMFD_NOACCEL = 99999

# Constant to represent a zero flux "albedo"
_ZERO_FLUX = 999.0

# Map that returns index of current direction in current matrix
_CURRENTS = {
    'out_left'  : 0, 'in_left'  : 1, 'out_right':  2, 'in_right': 3,
    'out_back'  : 4, 'in_back'  : 5, 'out_front':  6, 'in_front': 7,
    'out_bottom': 8, 'in_bottom': 9, 'out_top'  : 10, 'in_top'  : 11
}

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
        mesh spatial grid. If this option is used, a ``0`` is used for a
        non-accelerated region and a ``1`` is used for an accelerated region.
        For a simple 4x4 coarse mesh with a 2x2 fuel lattice surrounded by
        reflector, the map is:

        ::

            [0, 0, 0, 0,
             0, 1, 1, 0,
             0, 1, 1, 0,
             0, 0, 0, 0]

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
            check_value('CMFD mesh map', m, [0, 1])
        self._map = meshmap

    # REMOVE
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
            subelement.text = ' '.join(map(str, [self.map[i]+1 for i in range(len(self.map))]))

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

        # Check dimension defined
        if mesh.dimension is None:
            raise ValueError('CMFD mesh requires spatial '
                             'dimensions to be specified')

        # Check lower left defined
        if mesh.lower_left is None:
            raise ValueError('CMFD mesh requires lower left coordinates '
                             'to be specified')

        # Check that both upper right and width both not defined
        if mesh.upper_right is not None and mesh.width is not None:
            raise ValueError('Both upper right coordinates and width '
                             'cannot be specified for CMFD mesh')

        # Check that at least one of width or upper right is defined
        if mesh.upper_right is None and mesh.width is None:
            raise ValueError('CMFD mesh requires either upper right '
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
                raise ValueError('CMFD mesh requires upper right '
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
    openmc through this manner obviates the need for defining CMFD parameters
    through a cmfd.xml file. Instead, all input parameters should be passed through
    the CMFDRun initializer.


    Attributes
    ----------
    To add: indices: Stores spatial and group dimensions as [nx, ny, nz, ng]
            egrid: energy grid used for CMFD acceleration
            albedo: Albedo for global boundary conditions, taken from CMFD mesh. Set to [1,1,1,1,1,1] if not specified by user
            n_cmfd_resets: Number of elements in tally_reset, list that stores batches where CMFD tallies should be reset
            cmfd_mesh_id: Mesh id of openmc.capi.Mesh object that corresponds to the CMFD mesh
            cmfd_tally_ids: list of ids corresponding to CMFD tallies (details:)
            energy_filters: Boolean that stores whether energy filters should be created or not.
                Set to true if user specifies energy grid in CMFDMesh, false otherwise
            cmfd_on: Boolean to tell if cmfd solver should be initiated, based on whether current batch has reached variable
                cmfd_begin
            # Look at cmfd_header.F90 for description
            flux
            totalxs
            p1scattxs
            scattxs
            nfissxs
            diffcof
            dtilde
            dhat
            hxyz: mesh width
            current
            cmfd_src
            openmc_src
            sourcecounts
            weightfactors
            entropy
            balance
            src_cmp
            dom
            k_cmfd
            keff_bal
            mat_dim

    TODO: Put descriptions for all methods in CMFDRun
    TODO Get rid of CMFD constants
    TODO Get rid of unused variables defined in init
    TODO Make sure all self variables defined in init
    TODO Check to make sure no compatibility issues with numpy arrays for input variables
    TODO Create write_vector function in cmfd_solver.F90

    """

    def __init__(self):

        # Set CMFD default parameters based on cmfd_header.F90
        # Input parameters that user can define
        self._cmfd_begin = 1
        self._dhat_reset = False
        self._cmfd_display = 'balance'
        self._cmfd_downscatter = False
        self._cmfd_feedback = False
        self._cmfd_ktol = 1.e-8
        self._cmfd_mesh = None
        self._norm = 1.
        self._cmfd_power_monitor = False
        self._cmfd_run_adjoint = False
        self._cmfd_shift = 1.e6
        self._cmfd_stol = 1.e-8
        self._cmfd_reset = []
        self._cmfd_write_matrices = False

        # External variables used during runtime but users don't have control over
        self._indices = np.zeros(4, dtype=int)
        self._egrid = None
        self._albedo = None
        self._coremap = None
        self._n_cmfd_resets = 0
        self._cmfd_mesh_id = None
        self._cmfd_filter_ids = None
        self._cmfd_tally_ids = None
        self._energy_filters = None
        self._cmfd_on = False
        self._mat_dim = _CMFD_NOACCEL
        self._keff_bal = None
        self._cmfd_adjoint_type = "physical"
        self._keff = None
        self._adj_keff = None
        self._phi = None
        self._adj_phi = None

        # Numpy arrays used to build CMFD matrices
        self._flux = None
        self._totalxs = None
        self._p1scattxs = None
        self._scattxs = None
        self._nfissxs = None
        self._diffcof = None
        self._dtilde = None
        self._dhat = None
        self._hxyz = None
        self._current = None

        self._cmfd_src = None
        self._openmc_src = None
        self._sourcecounts = None
        self._weightfactors = None
        self._entropy = None
        self._balance = None
        self._src_cmp = None
        self._dom = None
        self._k_cmfd = None
        self._resnb = None

        self._time_cmfd = None
        self._time_cmfdbuild = None
        self._time_cmfdsolve = None

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
    def cmfd_ktol(self):
        return self._cmfd_ktol

    @property
    def cmfd_mesh(self):
        return self._cmfd_mesh

    @property
    def norm(self):
        return self._norm

    @property
    def cmfd_adjoint_type(self):
        return self._cmfd_adjoint_type

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
        self._cmfd_begin = cmfd_begin

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

    @cmfd_ktol.setter
    def cmfd_ktol(self, cmfd_ktol):
        check_type('CMFD eigenvalue tolerance', cmfd_ktol, Real)
        self._cmfd_ktol = cmfd_ktol

    @cmfd_mesh.setter
    def cmfd_mesh(self, mesh):
        check_type('CMFD mesh', mesh, CMFDMesh)

        # Check dimension defined
        if mesh.dimension is None:
            raise ValueError('CMFD mesh requires spatial '
                             'dimensions to be specified')

        # Check lower left defined
        if mesh.lower_left is None:
            raise ValueError('CMFD mesh requires lower left coordinates '
                             'to be specified')

        # Check that both upper right and width both not defined
        if mesh.upper_right is not None and mesh.width is not None:
            raise ValueError('Both upper right coordinates and width '
                             'cannot be specified for CMFD mesh')

        # Check that at least one of width or upper right is defined
        if mesh.upper_right is None and mesh.width is None:
            raise ValueError('CMFD mesh requires either upper right '
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
                raise ValueError('CMFD mesh requires upper right '
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

    @cmfd_adjoint_type.setter
    def cmfd_adjoint_type(self, adjoint_type):
        check_type('CMFD adjoint type', adjoint_type, str)
        check_value('CMFD adjoint type', adjoint_type,
                    ['math', 'phyical'])
        self._cmfd_adjoint_type = adjoint_type

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

    def run(self, mpi_procs=None, omp_num_threads=None):
        # Check number of OpenMP threads is valid input and initialize C API
        if omp_num_threads is not None:
            check_type('OpenMP num threads', omp_num_threads, Integral)
            openmc.capi.init(args=['-s',str(omp_num_threads)])
        else:
            comm = MPI.COMM_WORLD
            openmc.capi.init(intracomm=comm)

        # Configure cmfd parameters and tallies
        self._configure_cmfd()

        # Initialize simulation
        openmc.capi.simulation_init()

        while(True):
            # Run everything in next batch before initializing cmfd
            openmc.capi.next_batch_before_cmfd_init()

            # Initialize CMFD batch
            self._cmfd_init_batch()

            # Run everything in next batch in between initializing and
            # executing CMFD
            status = openmc.capi.next_batch_between_cmfd_init_execute()

            # Status determines whether batch should continue with a
            # CMFD update or skip it entirely if it is a restart run
            if status != 0:

                # Perform CMFD calculation if on
                if self._cmfd_on:
                    self._execute_cmfd()

                # Run everything in next batch after executing CMFD. Status
                # now determines whether another batch should be run or
                # simulation should be terminated.
                status = openmc.capi.next_batch_after_cmfd_execute()
            if status != 0:
                break

        # Finalize simuation
        openmc.capi.simulation_finalize()

        # Finalize and free memory
        openmc.capi.finalize()

    def _configure_cmfd(self):
        # Read in cmfd input defined in Python
        self._read_cmfd_input()

        # Initialize timers
        self._time_cmfd = 0.0
        self._time_cmfdbuild = 0.0
        self._time_cmfdsolve = 0.0

        # Initialize all numpy arrays used for cmfd solver
        self._allocate_cmfd()

    def _read_cmfd_input(self):
        # Print message to user
        if openmc.capi.settings.verbosity >= 7 and openmc.capi.settings.master:
            print(' Configuring CMFD parameters for simulation')

        # Check if CMFD mesh is defined
        if self._cmfd_mesh is None:
            raise ValueError('No CMFD mesh has been specified for '
                             'simulation')

        # Set spatial dimensions of CMFD object
        # Iterate through each element of self._cmfd_mesh.dimension
        # as could be length 2 or 3
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
            self._egrid = np.array([_ENERGY_MIN_NEUTRON, _ENERGY_MAX_NEUTRON])
            self._indices[3] = 1
            self._energy_filters = False

        # Set global albedo
        if self._cmfd_mesh.albedo is not None:
            self._albedo = np.array(self._cmfd_mesh.albedo)
        else:
            self._albedo = np.array([1.,1.,1.,1.,1.,1.])

        # Get acceleration map, otherwise set all regions to be accelerated
        if self._cmfd_mesh.map is not None:
            check_length('CMFD coremap', self._cmfd_mesh.map,
                         np.product(self._indices[0:3]))
            self._coremap = np.array(self._cmfd_mesh.map)
        else:
            self._coremap = np.ones((np.product(self._indices[0:3])))

        # Set number of batches where cmfd tallies should be reset
        if self._cmfd_reset is not None:
            self._n_cmfd_resets = len(self._cmfd_reset)

        # Create tally objects
        self._create_cmfd_tally()

    def _allocate_cmfd(self):
        # Extract spatial and energy indices
        nx = self._indices[0]
        ny = self._indices[1]
        nz = self._indices[2]
        ng = self._indices[3]

        # Allocate flux, cross sections and diffusion coefficient
        self._flux = np.zeros((nx, ny, nz, ng))
        self._totalxs = np.zeros((nx, ny, nz, ng))
        self._p1scattxs = np.zeros((nx, ny, nz, ng))
        self._scattxs = np.zeros((nx, ny, nz, ng, ng))  # Incoming, outgoing
        self._nfissxs = np.zeros((nx, ny, nz, ng, ng))  # Incoming, outgoing
        self._diffcof = np.zeros((nx, ny, nz, ng))

        # Allocate dtilde and dhat
        self._dtilde = np.zeros((nx, ny, nz, ng, 6))
        self._dhat = np.zeros((nx, ny, nz, ng, 6))

        # Allocate dimensions for each box (assume fixed mesh dimensions)
        # TODO Update this
        self._hxyz = np.zeros((3))

        # Allocate surface currents
        self._current = np.zeros((nx, ny, nz, ng, 12))

        # Allocate source distributions
        self._cmfd_src = np.zeros((nx, ny, nz, ng))
        self._openmc_src = np.zeros((nx, ny, nz, ng))

        # Allocate source weight modification variables
        self._sourcecounts = np.zeros((nx*ny*nz, ng))
        self._weightfactors = np.ones((nx, ny, nz, ng))

        # Allocate batchwise parameters
        self._entropy = []
        self._balance = []
        self._src_cmp = []
        self._dom = []
        self._k_cmfd = []

    def _cmfd_init_batch(self):
        # Get simulation parameters through C API
        current_batch = openmc.capi.settings.current_batch
        restart_run = openmc.capi.settings.restart_run
        restart_batch = openmc.capi.settings.restart_batch

        # Check to activate CMFD diffusion and possible feedback
        if self._cmfd_begin == current_batch:
            self._cmfd_on = True

        # TODO: Test restart_batch
        # If this is a restart run we are just replaying batches so don't
        # execute anything
        if restart_run and current_batch <= restart_batch:
            return

        # Check to reset tallies
        if (self._n_cmfd_resets > 0
                and current_batch in self._cmfd_reset):
            self._cmfd_tally_reset()

    def _execute_cmfd(self):
        # Run CMFD on single processor on master
        if openmc.capi.settings.master:
            #! Start cmfd timer
            time_start_cmfd = time.time()

            # Create cmfd data from OpenMC tallies
            self._set_up_cmfd()

            # Call solver
            self._cmfd_solver_execute()

            # Save k-effective
            self._k_cmfd.append(self._keff)

            # Check to perform adjoint on last batch
            if (openmc.capi.settings.current_batch == \
                    openmc.capi.settings.batches and self._cmfd_run_adjoint):
                self._cmfd_solver_execute(adjoint=True)

            # Calculate fission source
            self._calc_fission_source()

        # Calculate weight factors
        self._cmfd_reweight(True)

        # Stop cmfd timer
        if openmc.capi.settings.master:
            time_stop_cmfd = time.time()
            self._time_cmfd += time_stop_cmfd - time_start_cmfd

    def _cmfd_tally_reset(self):
        # Print message
        if openmc.capi.settings.verbosity >= 6 and openmc.capi.settings.master:
            print(' CMFD tallies reset')

        # Reset CMFD tallies
        tallies = openmc.capi.tallies
        for tally_id in self._cmfd_tally_ids:
            tallies[tally_id].reset()

    def _set_up_cmfd(self):
        # Set up CMFD coremap
        if (self._mat_dim == _CMFD_NOACCEL):
            self._set_coremap()

        # Calculate all cross sections based on reaction rates from last batch
        self._compute_xs()

        # Compute effective downscatter cross section
        if (self._cmfd_downscatter): self._compute_effective_downscatter()

        # Check neutron balance
        self._neutron_balance()

        # Calculate dtilde
        self._compute_dtilde()

        #self._compute_dtilde2()

        # Calculate dhat
        self._compute_dhat()

        #self._compute_dhat2()

    def _cmfd_solver_execute(self, adjoint=False):
        # Check for physical adjoint
        physical_adjoint = adjoint and self._cmfd_adjoint_type == 'physical'

        # Start timer for build
        time_start_buildcmfd = time.time()

        # Build loss and production matrices
        loss, prod = self._build_matrices(physical_adjoint)

        # Check for mathematical adjoint calculation
        if adjoint and self._cmfd_adjoint_type == 'math':
            loss, prod = self._compute_adjoint(loss, prod)

        # Stop timer for build
        time_stop_buildcmfd = time.time()
        self._time_cmfdbuild += time_stop_buildcmfd - time_start_buildcmfd

        # Begin power iteration
        time_start_solvecmfd = time.time()
        phi, keff, dom = self._execute_power_iter(loss, prod)
        time_stop_solvecmfd = time.time()
        self._time_cmfdsolve += time_stop_solvecmfd - time_start_solvecmfd

        # Save results, normalizing phi to sum to 1
        if adjoint:
            self._adj_keff = keff
            self._adj_phi = phi/np.sqrt(np.sum(phi*phi))
        else:
            self._keff = keff
            self._phi = phi/np.sqrt(np.sum(phi*phi))

        self._dom.append(dom)

        # TODO Write out flux vector
        '''
        if (cmfd_write_matrices) then
          if (adjoint_calc) then
            filename = 'adj_fluxvec.dat'
          else
            filename = 'fluxvec.dat'
          end if
          ! TODO: call phi_n % write(filename)
        end if
        '''

    def _calc_fission_source(self):
        # Extract number of groups and number of accelerated regions
        nx = self._indices[0]
        ny = self._indices[1]
        nz = self._indices[2]
        ng = self._indices[3]
        n = self._mat_dim

        # Reset cmfd source to 0
        self._cmfd_src.fill(0.)

        # Calculate volume
        vol = np.product(self._hxyz)

        # Reshape phi by number of groups
        phi = self._phi.reshape((n, ng))

        # Extract indices of coremap that are accelerated
        idx = np.where(self._coremap != _CMFD_NOACCEL)

        # Initialize CMFD flux map that maps phi to actualy spatial and group
        # indices of problem
        cmfd_flux = np.zeros((nx, ny, nz, ng))

        # Loop over all groups and set CMFD flux based on indices of coremap
        # and values of phi
        for g in range(ng):
            cmfd_flux[idx + (np.full((n,),g),)]  = phi[:,g]

        # Compute fission source
        self._cmfd_src = np.sum(self._nfissxs[:,:,:,:,:] * \
                         cmfd_flux[:,:,:,:,np.newaxis], axis=3) * vol

    def _cmfd_reweight(self, new_weights):
        # Compute new weight factors
        if new_weights:

            # Set weight factors to default 1.0
            self._weightfactors.fill(1.0)

            # Count bank site in mesh and reverse due to egrid structured
            outside = self._count_bank_sites()

            if openmc.capi.settings.master and outside:
                raise OpenMCError('Source sites outside of the CMFD mesh!')

            if openmc.capi.settings.master:
                norm = np.sum(self._sourcecounts) / np.sum(self._cmfd_src)
                sourcecounts = np.flip(
                        self._sourcecounts.reshape(self._cmfd_src.shape), \
                        axis=3)
                divide_condition = np.logical_and(sourcecounts > 0,
                                                  self._cmfd_src > 0)
                self._weightfactors = np.divide(self._cmfd_src * norm, \
                        sourcecounts, where=divide_condition, \
                        out=np.ones_like(self._cmfd_src))

            if not self._cmfd_feedback:
                return

            comm = MPI.COMM_WORLD
            self._weightfactors = comm.bcast(self._weightfactors, root=0)

            m = openmc.capi.meshes[self._cmfd_mesh_id]
            bank = openmc.capi.source_bank()
            energy = self._egrid
            ng = self._indices[3]

            for i in range(bank.size):
                # Get xyz location of source point
                xyz = bank[i][1]
                # Get energy of source point
                source_energy = bank[i][3]

                # Determine which mesh location of source point
                ijk = np.floor((xyz - m.lower_left) / m.width).astype(int)

                # Determine energy bin of source point
                e_bin = np.where(source_energy < energy[0], 0,
                        np.where(source_energy > energy[-1], ng - 1, \
                        np.digitize(source_energy, energy) - 1))

                # Issue warnings if energy below or above defined grid
                if openmc.capi.settings.master and source_energy < energy[0]:
                    print(' WARNING: Source pt below energy grid')
                if openmc.capi.settings.master and source_energy > energy[-1]:
                    print(' WARNING: Source pt above energy grid')

                # Reverse index of bin due to egrid structure
                e_bin = ng - e_bin - 1

                # Reweight particle
                bank[i][0] *= self._weightfactors[ijk[0], ijk[1], ijk[2], e_bin]

    def _count_bank_sites(self):
        comm = MPI.COMM_WORLD

        m = openmc.capi.meshes[self._cmfd_mesh_id]
        bank = openmc.capi.source_bank()
        energy = self._egrid
        sites_outside = np.zeros(1, dtype=bool)
        ng = self._indices[3]

        # Initiate variables
        outside = np.zeros(1, dtype=bool)
        count = np.zeros(self._sourcecounts.shape)

        source_xyz = np.array([bank[i][1] for i in range(bank.size)])
        source_energies = np.array([bank[i][3] for i in range(bank.size)])

        mesh_locations = np.floor((source_xyz - m.lower_left) / m.width)
        mesh_bins = mesh_locations[:,2] * m.dimension[2] + \
                    mesh_locations[:,1] * m.dimension[1] + mesh_locations[:,0]

        if np.any(mesh_bins < 0) or np.any(mesh_bins >= np.prod(m.dimension)):
            outside[0] = True

        energy_bins = np.where(source_energies < energy[0], 0,
                np.where(source_energies > energy[-1], ng - 1, \
                         np.digitize(source_energies, energy) - 1))

        idx, counts = np.unique(np.array([mesh_bins, energy_bins]), axis=1,
                                return_counts=True)

        count[idx[0].astype(int), idx[1].astype(int)] = counts

        comm.Reduce(count, self._sourcecounts, MPI.SUM, 0)
        comm.Reduce(outside, sites_outside, MPI.LOR, 0)

        # TODO ifdef mpi?
        # how to define comm?
        # issue with output out of order?
        return sites_outside[0]




    def _build_matrices(self, adjoint):
        # Build loss matrix
        loss = self._build_loss_matrix(adjoint)
        self._build_loss_matrix2(adjoint)

        # Build production matrix
        prod = self._build_prod_matrix(adjoint)
        self._build_prod_matrix2(adjoint)

        # TODO Write out matrices
        #if self._cmfd_write_matrices:
        #    self._write_matrix(loss, 'loss.dat')
        #    self._write_matrix(prod, 'prod.dat')

        return loss, prod

    def _compute_adjoint(self, loss, prod):
        # Transpose matrices
        loss = np.transpose(loss)
        prod = np.transpose(prod)

        # TODO Write out matrices
        #if self._cmfd_write_matrices:
        #    self._write_matrix(loss, 'adj_loss.dat')
        #    self._write_matrix(prod, 'adj_prod.dat')

        return loss, prod

    def _build_loss_matrix(self, adjoint):
        # Extract spatial and energy indices and define matrix dimension
        nx = self._indices[0]
        ny = self._indices[1]
        nz = self._indices[2]
        ng = self._indices[3]
        n = self._mat_dim*ng

        # Allocate matrix
        loss = np.zeros((n, n))

        # Create single vector of these indices for boundary calculation
        nxyz = np.array([[0,nx-1], [0,ny-1], [0,nz-1]])

        # Allocate leakage coefficients in front of cell flux
        jo = np.zeros((6,))

        for irow in range(n):
            # Get indices for row in matrix
            i,j,k,g = self._matrix_to_indices(irow, nx, ny, nz, ng)

            # Retrieve cell data
            totxs = self._totalxs[i,j,k,g]
            scattxsgg = self._scattxs[i,j,k,g,g]
            dtilde = self._dtilde[i,j,k,g,:]
            hxyz = self._hxyz
            dhat = self._dhat[i,j,k,g,:]

            # Create boundary vector
            bound = np.repeat([i,j,k], 2)

            # Begin loop over leakages
            for l in range(6):
                # Define (x,y,z) and (-,+) indices
                xyz_idx = int(l/2)  # x=0, y=1, z=2
                dir_idx = l % 2     # -=0, +=1

                # Calculate spatial indices of neighbor
                neig_idx = [i,j,k]  # Begin with i,j,k
                shift_idx = 2*(l % 2) - 1 # shift neig by -1 or +1
                neig_idx[xyz_idx] += shift_idx

                # Check for global boundary
                if bound[l] != nxyz[xyz_idx, dir_idx]:

                    # Check that neighbor is not reflector
                    if self._coremap[tuple(neig_idx)] != _CMFD_NOACCEL:
                        # Compute leakage coefficient for neighbor
                        jn = -1.0 * dtilde[l] + shift_idx*dhat[l]

                        # Get neighbor matrix index
                        neig_mat_idx = self._indices_to_matrix(neig_idx[0], \
                                neig_idx[1], neig_idx[2], g, ng)
                        # Compute value and record to bank
                        val = jn/hxyz[xyz_idx]
                        loss[irow, neig_mat_idx] = val

                # Compute leakage coefficient for target
                jo[l] = shift_idx*dtilde[l] + dhat[l]

                # Calculate net leakage coefficient for target
                jnet = (jo[1] - jo[0])/hxyz[0] + (jo[3] - jo[2])/hxyz[1] + \
                    (jo[5] - jo[4])/hxyz[2]

                # Calculate loss of neutrons
                val = jnet + totxs - scattxsgg
                loss[irow, irow] = val

                # Begin loop over off diagonal in-scattering
                for h in range(ng):
                    # Cycle though if h=g, value already banked in removal xs
                    if h == g:
                        continue

                    # Get neighbor matrix index
                    scatt_mat_idx = self._indices_to_matrix(i,j,k, h, ng)

                    # Check for adjoint
                    if adjoint:
                        # Get scattering macro xs, transposed!
                        scattxshg = self._scattxs[i, j, k, g, h]
                    else:
                        # Get scattering macro xs
                        scattxshg = self._scattxs[i, j, k, h, g]

                    # Negate the scattering xs
                    val = -1.0*scattxshg

                    # Record value in matrix
                    loss[irow, scatt_mat_idx] = val

        '''
        print("Loss matrix:")
        for i in range(loss.shape[0]):
            print_str = ""
            for j in range(loss.shape[1]):
                if loss[i,j] != 0.0:
                    print_str += "%.3g\t" % loss[i, j]
                else:
                    print_str += "%.3f\t" % loss[i, j]
            print(print_str)
        '''

        return loss

    def _build_loss_matrix2(self, adjoint):
        # TODO build as csr matrix
        # Extract spatial and energy indices and define matrix dimension
        ng = self._indices[3]
        n = self._mat_dim*ng

        # Allocate matrix
        loss2 = np.zeros((n, n))

        # Define net leakage coefficient for each matrix element
        jnet = ((1.0 * self._dtilde[:,:,:,:,1] + self._dhat[:,:,:,:,1]) - \
               (-1.0 * self._dtilde[:,:,:,:,0] + self._dhat[:,:,:,:,0])) / \
               self._hxyz[0] + \
               ((1.0 * self._dtilde[:,:,:,:,3] + self._dhat[:,:,:,:,3]) - \
               (-1.0 * self._dtilde[:,:,:,:,2] + self._dhat[:,:,:,:,2])) / \
               self._hxyz[1] + \
               ((1.0 * self._dtilde[:,:,:,:,5] + self._dhat[:,:,:,:,5]) - \
               (-1.0 * self._dtilde[:,:,:,:,4] + self._dhat[:,:,:,:,4])) / \
               self._hxyz[2]

        # Shift coremap in all directions to determine whether leakage term
        # should be defined for particular cell in matrix
        coremap_minusx = np.pad(self._coremap, ((1,0),(0,0),(0,0)), mode='constant',
                            constant_values=_CMFD_NOACCEL)[:-1,:,:]

        coremap_plusx = np.pad(self._coremap, ((0,1),(0,0),(0,0)), mode='constant',
                            constant_values=_CMFD_NOACCEL)[1:,:,:]

        coremap_minusy = np.pad(self._coremap, ((0,0),(1,0),(0,0)), mode='constant',
                            constant_values=_CMFD_NOACCEL)[:,:-1,:]

        coremap_plusy = np.pad(self._coremap, ((0,0),(0,1),(0,0)), mode='constant',
                            constant_values=_CMFD_NOACCEL)[:,1:,:]

        coremap_minusz = np.pad(self._coremap, ((0,0),(0,0),(1,0)), mode='constant',
                            constant_values=_CMFD_NOACCEL)[:,:,:-1]

        coremap_plusz = np.pad(self._coremap, ((0,0),(0,0),(0,1)), mode='constant',
                            constant_values=_CMFD_NOACCEL)[:,:,1:]

        condition = np.logical_and(self._coremap != _CMFD_NOACCEL,
                                   coremap_minusy != _CMFD_NOACCEL)

        for g in range(ng):
            # Leakage terms
            condition = np.logical_and(self._coremap != _CMFD_NOACCEL,
                                       coremap_minusx != _CMFD_NOACCEL)
            idx_x = ng * (self._coremap[condition]) + g
            idx_y = ng * (coremap_minusx[condition]) + g
            vals = (-1.0 * self._dtilde[:,:,:,g,0] -
                   self._dhat[:,:,:,g,0])[condition] / self._hxyz[0]
            loss2[idx_x, idx_y] = vals

            condition = np.logical_and(self._coremap != _CMFD_NOACCEL,
                                       coremap_plusx != _CMFD_NOACCEL)
            idx_x = ng * (self._coremap[condition]) + g
            idx_y = ng * (coremap_plusx[condition]) + g
            vals = (-1.0 * self._dtilde[:,:,:,g,1] +
                   self._dhat[:,:,:,g,1])[condition] / self._hxyz[0]
            loss2[idx_x, idx_y] = vals

            condition = np.logical_and(self._coremap != _CMFD_NOACCEL,
                                       coremap_minusy != _CMFD_NOACCEL)
            idx_x = ng * (self._coremap[condition]) + g
            idx_y = ng * (coremap_minusy[condition]) + g
            vals = (-1.0 * self._dtilde[:,:,:,g,2] -
                   self._dhat[:,:,:,g,2])[condition] / self._hxyz[1]
            loss2[idx_x, idx_y] = vals

            condition = np.logical_and(self._coremap != _CMFD_NOACCEL,
                                       coremap_plusy != _CMFD_NOACCEL)
            idx_x = ng * (self._coremap[condition]) + g
            idx_y = ng * (coremap_plusy[condition]) + g
            vals = (-1.0 * self._dtilde[:,:,:,g,3] +
                   self._dhat[:,:,:,g,3])[condition] / self._hxyz[1]
            loss2[idx_x, idx_y] = vals

            condition = np.logical_and(self._coremap != _CMFD_NOACCEL,
                                       coremap_minusz != _CMFD_NOACCEL)
            idx_x = ng * (self._coremap[condition]) + g
            idx_y = ng * (coremap_minusz[condition]) + g
            vals = (-1.0 * self._dtilde[:,:,:,g,4] -
                   self._dhat[:,:,:,g,4])[condition] / self._hxyz[1]
            loss2[idx_x, idx_y] = vals

            condition = np.logical_and(self._coremap != _CMFD_NOACCEL,
                                       coremap_plusz != _CMFD_NOACCEL)
            idx_x = ng * (self._coremap[condition]) + g
            idx_y = ng * (coremap_plusz[condition]) + g
            vals = (-1.0 * self._dtilde[:,:,:,g,5] +
                   self._dhat[:,:,:,g,5])[condition] / self._hxyz[1]
            loss2[idx_x, idx_y] = vals

            # Loss of neutrons
            condition = self._coremap != _CMFD_NOACCEL
            idx_x = ng * (self._coremap[condition]) + g
            idx_y = idx_x
            vals = (jnet[:,:,:,g] + self._totalxs[:,:,:,g] - \
                   self._scattxs[:,:,:,g,g])[condition]
            loss2[idx_x, idx_y] = vals

            # Off diagonal in-scattering
            for h in range(ng):
                #TODO check for adjoint (see cmfd_loss_operator.F90)
                if h != g:
                    condition = self._coremap != _CMFD_NOACCEL
                    idx_x = ng * (self._coremap[condition]) + g
                    idx_y = ng * (self._coremap[condition]) + h
                    vals = (-1.0 * self._scattxs[:, :, :, h, g])[condition]
                    loss2[idx_x, idx_y] = vals


    def _build_prod_matrix(self, adjoint):
        # Extract spatial and energy indices and define matrix dimension
        nx = self._indices[0]
        ny = self._indices[1]
        nz = self._indices[2]
        ng = self._indices[3]
        n = self._mat_dim*ng

        # Allocate matrix
        prod = np.zeros((n, n))

        for irow in range(n):
            # Get indices for row in matrix
            i,j,k,g = self._matrix_to_indices(irow, nx, ny, nz, ng)

            # Check if at a reflector
            if self._coremap[i,j,k] == _CMFD_NOACCEL:
                continue

            # Loop around all other groups
            for h in range(ng):
                # Get matrix column location
                hmat_idx = self._indices_to_matrix(i,j,k, h, ng)
                # Check for adjoint and bank val
                if adjoint:
                    # Get nu-fission cross section from cell, transposed!
                    nfissxs = self._nfissxs[i, j, k, g, h]
                else:
                    # Get nu-fission cross section from cell
                    nfissxs = self._nfissxs[i, j, k, h, g]

                # Set as value to be recorded
                val = nfissxs

                # record value in matrix
                prod[irow, hmat_idx] = val

        '''
        print("Prod matrix:")
        for i in range(prod.shape[0]):
            print_str = ""
            for j in range(prod.shape[1]):
                if prod[i,j] != 0.0:
                    print_str += "%.3g\t" % prod[i, j]
                else:
                    print_str += "%.3f\t" % prod[i, j]
            print(print_str)
        '''
        return prod

    def _build_prod_matrix2(self, adjoint):
        # Extract spatial and energy indices and define matrix dimension
        ng = self._indices[3]
        n = self._mat_dim*ng

        # Allocate matrix
        prod2 = np.zeros((n, n))

        for g in range(ng):
            for h in range(ng):
                #TODO check for adjoint (see cmfd_prod_operator.F90)
                condition = self._coremap != _CMFD_NOACCEL
                idx_x = ng * (self._coremap[condition]) + g
                idx_y = ng * (self._coremap[condition]) + h
                vals = (self._nfissxs[:, :, :, h, g])[condition]
                prod2[idx_x, idx_y] = vals

        prod = prod2

    def _matrix_to_indices(self, irow, nx, ny, nz, ng):
        # Get indices from coremap
        g = irow % ng
        spatial_idx = np.where(self._coremap == int(irow/ng))
        i = spatial_idx[0][0]
        j = spatial_idx[1][0]
        k = spatial_idx[2][0]

        return i, j, k, g

    def _indices_to_matrix(self, i, j, k, g, ng):
        # Get matrix index from coremap
        matidx = ng*(self._coremap[i,j,k]) + g
        return matidx

    def _execute_power_iter(self, loss, prod):
        # Get problem size
        n = loss.shape[0]

        # Set up flux vectors, intital guess set to 1
        phi_n = np.ones((n,))
        phi_o = np.ones((n,))

        # Set up source vectors
        s_n = np.zeros((n,))
        s_o = np.zeros((n,))
        serr_v = np.zeros((n,))

        # Set initial guess
        k_n = openmc.capi.keff()[0]
        k_o = k_n
        dw = self._cmfd_shift
        k_s = k_o + dw
        k_ln = 1.0/(1.0/k_n - 1.0/k_s)
        k_lo = k_ln

        # Set norms to 0
        norm_n = 0.0
        norm_o = 0.0

        # Maximum number of power iterations
        maxits = 10000

        # Perform Wielandt shift
        loss -= 1.0/k_s*prod

        # Convert matrices to csr matrix in order to use scipy sparse solver
        prod = sparse.csr_matrix(prod)
        loss = sparse.csr_matrix(loss)

        # Begin power iteration
        for i in range(maxits):
            # Check if reach max number of iterations
            if i == maxits - 1:
                raise OpenMCError('Reached maximum iterations in CMFD power '
                                  'iteration solver.')

            # Compute source vector
            s_o = prod.dot(phi_o)

            # Normalize source vector
            s_o /= k_lo

            # Compute new flux vector with scipy sparse solver
            phi_n = sparse.linalg.spsolve(loss, s_o)

            # Compute new source vector
            s_n = prod.dot(phi_n)

            # Compute new shifted eigenvalue
            k_ln = np.sum(s_n) / np.sum(s_o)

            # Compute new eigenvalue
            k_n = 1.0/(1.0/k_ln + 1.0/k_s)

            # Renormalize the old source
            s_o *= k_lo

            # Check convergence
            iconv, norm_n = self._check_convergence(s_n, s_o, k_n, k_o)

            # If converged, calculate dominance ratio and break from loop
            if iconv:
                dom = norm_n / norm_o
                return phi_n, k_n, dom

            # Record old values if not converged
            phi_o = phi_n
            k_o = k_n
            k_lo = k_ln
            norm_o = norm_n

    def _check_convergence(self, s_n, s_o, k_n, k_o):
        # Calculate error in keff
        kerr = abs(k_o - k_n)/k_n

        # Calculate max error in source
        serr = np.sqrt(np.sum(np.where(s_n>0, ((s_n-s_o)/s_n)**2,0))/len(s_n))

        # Check for convergence
        iconv = kerr < self._cmfd_ktol and serr < self._cmfd_stol

        # TODO Print out to user
        #if (cmfd_power_monitor .and. master) then
        #  write(OUTPUT_UNIT,FMT='(I0,":",T10,"k-eff: ",F0.8,T30,"k-error: ", &
        #       &1PE12.5,T55, "src-error: ",1PE12.5,T80,I0)') iter, k_n, kerr, &
        #        serr, innerits

        # Return source error and convergence logical back to solver
        return iconv, serr

    def _set_coremap(self):
        self._mat_dim = np.sum(self._coremap)

        # Define coremap as cumulative sum over accelerated regions,
        # otherwise set value to _CMFD_NOACCEL
        self._coremap = np.where(self._coremap==0, _CMFD_NOACCEL,
                        np.cumsum(self._coremap)-1)

    def _compute_xs(self):
        # Extract energy indices
        ng = self._indices[3]

        # Set flux object and source distribution all to zeros
        self._flux.fill(0.)
        self._openmc_src.fill(0.)

        # Set mesh widths
        self._hxyz = openmc.capi.meshes[self._cmfd_mesh_id].width
        # self._hxyz[:,:,:,] = openmc.capi.meshes[self._cmfd_mesh_id].width

        # Reset keff_bal to zero
        self._keff_bal = 0.

        # Get tallies in-memory
        tallies = openmc.capi.tallies

        # Set conditional numpy array as boolean vector based on coremap
        # Repeat each value for number of groups in problem
        is_cmfd_accel = np.repeat(self._coremap != _CMFD_NOACCEL, ng)

        # Get flux from CMFD tally 0
        tally_id = self._cmfd_tally_ids[0]
        tally_results = tallies[tally_id].results[:,0,1]
        flux = np.where(is_cmfd_accel, tally_results, 0.)

        # Detect zero flux, abort if located
        if np.any(flux[is_cmfd_accel] < _TINY_BIT):
            # Get index of zero flux in flux array
            idx = np.argmax(np.where(is_cmfd_accel, flux, 1) < _TINY_BIT)

            # Convert scalar idx to index in flux matrix
            mat_idx = np.unravel_index(idx, self._flux.shape)

            # Throw error message (one-based indexing)
            # Index of group is flipped
            err_message = 'Detected zero flux without coremap overlay' + \
                          ' at mesh: (' + \
                          ', '.join(str(i+1) for i in mat_idx[:-1]) + \
                          ') in group ' + str(ng-mat_idx[-1])
            raise OpenMCError(err_message)

        # Store flux and reshape
        # Flux is flipped in energy axis as tally results are given in reverse
        # order of energy group
        self._flux = np.flip(flux.reshape(self._flux.shape), axis=3)

        # Get total rr and convert to total xs from CMFD tally 0
        tally_results = tallies[tally_id].results[:,1,1]
        totalxs = np.divide(tally_results, flux, \
                            where=flux>0, out=np.zeros_like(tally_results))

        # Store total xs and reshape
        # Total xs is flipped in energy axis as tally results are given in
        # reverse order of energy group
        self._totalxs = np.flip(totalxs.reshape(self._totalxs.shape), axis=3)

        # Get scattering xs from CMFD tally 1
        # flux is repeated to account for extra dimensionality of scattering xs
        tally_id = self._cmfd_tally_ids[1]
        tally_results = tallies[tally_id].results[:,0,1]
        scattxs = np.divide(tally_results, \
                            np.repeat(flux, ng), \
                            where=np.repeat(flux>0, ng), \
                            out=np.zeros_like(tally_results))

        # Store scattxs and reshape
        # Scattering xs is flipped in both incoming and outgoing energy axes
        # as tally results are given in reverse order of energy group
        self._scattxs = np.flip(scattxs.reshape(self._scattxs.shape), axis=3)
        self._scattxs = np.flip(self._scattxs, axis=4)

        # Get nu-fission xs from CMFD tally 1
        # flux is repeated to account for extra dimensionality of nu-fission xs
        tally_results = tallies[tally_id].results[:,1,1]
        num_realizations = tallies[tally_id].num_realizations
        nfissxs = np.divide(tally_results, \
                            np.repeat(flux, ng), \
                            where=np.repeat(flux>0, ng), \
                            out=np.zeros_like(tally_results))

        # Store nfissxs and reshape
        # Nu-fission xs is flipped in both incoming and outgoing energy axes
        # as tally results are given in reverse order of energy group
        self._nfissxs = np.flip(nfissxs.reshape(self._nfissxs.shape), axis=3)
        self._nfissxs = np.flip(self._nfissxs, axis=4)

        # Openmc source distribution is sum of nu-fission rr in incoming energies
        self._openmc_src = np.sum(self._nfissxs*self._flux[:,:,:,:,np.newaxis],
                                  axis=3)

        # Compute k_eff from source distribution
        self._keff_bal = np.sum(self._openmc_src) / num_realizations

        # Normalize openmc source distribution
        self._openmc_src /= np.sum(self._openmc_src) * self._norm

        # Get surface currents from CMFD tally 2
        tally_id = self._cmfd_tally_ids[2]
        tally_results = tallies[tally_id].results[:,0,1]

        # Filter tally results to include only accelerated regions
        tally_results = np.where(np.repeat(flux>0, 12), tally_results, 0.)

        # Reshape and store current
        # Current is flipped in energy axis as tally results are given in
        # reverse order of energy group
        self._current = np.flip(tally_results.reshape(self._current.shape), \
                                axis=3)

        # Get p1 scatter xs from CMFD tally 3
        tally_id = self._cmfd_tally_ids[3]
        tally_results = tallies[tally_id].results[:,0,1]

        # Reshape and extract only p1 data from tally results (no need for p0 data)
        p1scattrr = tally_results.reshape(self._p1scattxs.shape+(2,))[:,:,:,1]

        # Store p1 scatter xs
        # p1 scatter xs is flipped in energy axis as tally results are given in
        # reverse order of energy group
        self._p1scattxs = np.divide(np.flip(p1scattrr, axis=3), self._flux, \
                                    where=self._flux>0, \
                                    out=np.zeros_like(p1scattrr))

        # Calculate and store diffusion coefficient
        self._diffcof = np.where(self._flux>0, 1.0 / (3.0 * \
                                 (self._totalxs - self._p1scattxs)), 0.)

        # Reshape coremap to three dimensional array as all cross section data
        # has been reshaped
        self._coremap = self._coremap.reshape(self._indices[0:3])

    def _compute_effective_downscatter(self):
        # Extract energy index
        ng = self._indices[3]

        # Return if not two groups
        if ng != 2:
            return

        # Extract cross sections and flux for each group
        flux1 = self._flux[:,:,:,0]
        flux2 = self._flux[:,:,:,1]
        sigt1 = self._totalxs[:,:,:,0]
        sigt2 = self._totalxs[:,:,:,1]
        # First energy index is incoming, second is outgoing
        sigs11 = self._scattxs[:,:,:,0,0]
        sigs21 = self._scattxs[:,:,:,1,0]
        sigs12 = self._scattxs[:,:,:,0,1]
        sigs22 = self._scattxs[:,:,:,1,1]

        # Compute absorption xs
        siga1 = sigt1 - sigs11 - sigs12
        siga2 = sigt2 - sigs22 - sigs21

        # Compute effective downscatter XS
        sigs12_eff = sigs12 - sigs21 * np.divide(flux2, flux1, where=flux1>0,
                                                 out=np.zeros_like(flux2))

        # Recompute total cross sections and record
        self._totalxs[:,:,:,0] = siga1 + sigs11 + sigs12_eff
        self._totalxs[:,:,:,1] = siga2 + sigs22

        # Record effective dowmscatter xs
        self._scattxs[:,:,:,0,1] = sigs12_eff

        # Zero out upscatter cross section
        self._scattxs[:,:,:,1,0] = 0.0

    def _neutron_balance(self):
        # Extract energy indices
        ng = self._indices[3]

        # Get openmc k-effective
        keff = openmc.capi.keff()[0]

        # Define leakage in each mesh cell and energy group
        leakage = ((self._current[:,:,:,:,_CURRENTS['out_right']] - \
            self._current[:,:,:,:,_CURRENTS['in_right']]) - \
            (self._current[:,:,:,:,_CURRENTS['in_left']] - \
            self._current[:,:,:,:,_CURRENTS['out_left']])) + \
            ((self._current[:,:,:,:,_CURRENTS['out_front']] - \
            self._current[:,:,:,:,_CURRENTS['in_front']]) - \
            (self._current[:,:,:,:,_CURRENTS['in_back']] - \
            self._current[:,:,:,:,_CURRENTS['out_back']])) + \
            ((self._current[:,:,:,:,_CURRENTS['out_top']] - \
            self._current[:,:,:,:,_CURRENTS['in_top']]) - \
            (self._current[:,:,:,:,_CURRENTS['in_bottom']] - \
            self._current[:,:,:,:,_CURRENTS['out_bottom']]))

        # Compute total rr
        interactions = self._totalxs * self._flux

        # Compute scattering rr by broadcasting flux in outgoing energy and
        # summing over incoming energy
        scattering = np.sum(self._scattxs * self._flux[:,:,:,:,np.newaxis],
                            axis=3)

        # Compute fission rr by broadcasting flux in outgoing energy and
        # summing over incoming energy
        fission = np.sum(self._nfissxs * self._flux[:,:,:,:,np.newaxis],
                         axis=3)

        # Compute residual
        res = leakage + interactions - scattering - (1.0 / keff) * fission

        # Normalize res by flux and bank res
        self._resnb = np.divide(res, self._flux, where=self._flux>0)

        # Calculate RMS and record for this batch
        self._balance.append(np.sqrt(
            np.sum(np.multiply(self._resnb, self._resnb)) / \
            np.count_nonzero(self._resnb)))

    def _compute_dtilde2(self):
        #TODO add coments for this method
        dtilde2 = np.zeros(self._dtilde.shape)

        is_accel = self._coremap != _CMFD_NOACCEL
        is_zero_flux_alb = abs(self._albedo - _ZERO_FLUX) < _TINY_BIT

        dtilde2[0,:,:,:,0] = np.where(is_accel[0,:,:,np.newaxis],
            np.where(is_zero_flux_alb[0], 2.0 * self._diffcof[0,:,:,:] / \
                     self._hxyz[0],
                     (2.0 * self._diffcof[0,:,:,:] * \
                     (1.0 - self._albedo[0])) / \
                     (4.0 * self._diffcof[0,:,:,:] * \
                     (1.0 + self._albedo[0]) + \
                     (1.0 - self._albedo[0]) * self._hxyz[0])), 0)

        dtilde2[-1,:,:,:,1] = np.where(is_accel[-1,:,:,np.newaxis],
            np.where(is_zero_flux_alb[1], 2.0 * self._diffcof[-1,:,:,:] / \
                     self._hxyz[0],
                     (2.0 * self._diffcof[-1,:,:,:] * \
                     (1.0 - self._albedo[1])) / \
                     (4.0 * self._diffcof[-1,:,:,:] * \
                     (1.0 + self._albedo[1]) + \
                     (1.0 - self._albedo[1]) * self._hxyz[0])), 0)

        dtilde2[:,0,:,:,2] = np.where(is_accel[:,0,:,np.newaxis],
            np.where(is_zero_flux_alb[2], 2.0 * self._diffcof[:,0,:,:] / \
                     self._hxyz[1],
                     (2.0 * self._diffcof[:,0,:,:] * \
                     (1.0 - self._albedo[2])) / \
                     (4.0 * self._diffcof[:,0,:,:] * \
                     (1.0 + self._albedo[2]) + \
                     (1.0 - self._albedo[2]) * self._hxyz[1])), 0)

        dtilde2[:,-1,:,:,3] = np.where(is_accel[:,-1,:,np.newaxis],
            np.where(is_zero_flux_alb[3], 2.0 * self._diffcof[:,-1,:,:] / \
                     self._hxyz[1],
                     (2.0 * self._diffcof[:,-1,:,:] * \
                     (1.0 - self._albedo[3])) / \
                     (4.0 * self._diffcof[:,-1,:,:] * \
                     (1.0 + self._albedo[3]) + \
                     (1.0 - self._albedo[3]) * self._hxyz[1])), 0)

        dtilde2[:,:,0,:,4] = np.where(is_accel[:,:,0,np.newaxis],
            np.where(is_zero_flux_alb[4], 2.0 * self._diffcof[:,:,0,:] / \
                     self._hxyz[2],
                     (2.0 * self._diffcof[:,:,0,:] * \
                     (1.0 - self._albedo[4])) / \
                     (4.0 * self._diffcof[:,:,0,:] * \
                     (1.0 + self._albedo[4]) + \
                     (1.0 - self._albedo[4]) * self._hxyz[2])), 0)

        dtilde2[:,:,-1,:,5] = np.where(is_accel[:,:,-1,np.newaxis],
            np.where(is_zero_flux_alb[5], 2.0 * self._diffcof[:,:,-1,:] / \
                     self._hxyz[2],
                     (2.0 * self._diffcof[:,:,-1,:] * \
                     (1.0 - self._albedo[5])) / \
                     (4.0 * self._diffcof[:,:,-1,:] * \
                     (1.0 + self._albedo[5]) + \
                     (1.0 - self._albedo[5]) * self._hxyz[2])), 0)

        ref_albedo = np.divide(self._current[:,:,:,:,_CURRENTS['in_left']],
                self._current[:,:,:,:,_CURRENTS['out_left']],
                where=self._current[:,:,:,:,_CURRENTS['out_left']] > 1.0e-10,
                out=np.ones_like(self._current[:,:,:,:,_CURRENTS['out_left']]))
        adj_reflector = np.roll(self._coremap, 1, axis=0) == _CMFD_NOACCEL
        neig_dc = np.roll(self._diffcof, 1, axis=0)
        # Define neg_hxyz

        dtilde2[1:,:,:,:,0] = np.where(is_accel[1:,:,:,np.newaxis], \
            np.where(adj_reflector[1:,:,:,np.newaxis],
                     (2.0 * self._diffcof[1:,:,:,:] * \
                     (1.0 - ref_albedo[1:,:,:,:])) / \
                     (4.0 * self._diffcof[1:,:,:,:] * \
                     (1.0 + ref_albedo[1:,:,:,:]) + \
                     (1.0 - ref_albedo[1:,:,:,:]) * self._hxyz[0]),
                     (2.0 * self._diffcof[1:,:,:,:] * neig_dc[1:,:,:,:]) / \
                     (self._hxyz[0] * self._diffcof[1:,:,:,:] + \
                     self._hxyz[0] * neig_dc[1:,:,:,:])), 0.0)

        ref_albedo = np.divide(self._current[:,:,:,:,_CURRENTS['in_right']],
                self._current[:,:,:,:,_CURRENTS['out_right']],
                where=self._current[:,:,:,:,_CURRENTS['out_right']] > 1.0e-10,
                out=np.ones_like(self._current[:,:,:,:,_CURRENTS['out_right']]))
        adj_reflector = np.roll(self._coremap, -1, axis=0) == _CMFD_NOACCEL
        neig_dc = np.roll(self._diffcof, -1, axis=0)
        # Define neg_hxyz

        dtilde2[:-1,:,:,:,1] = np.where(is_accel[:-1,:,:,np.newaxis], \
            np.where(adj_reflector[:-1,:,:,np.newaxis],
                     (2.0 * self._diffcof[:-1,:,:,:] * \
                     (1.0 - ref_albedo[:-1,:,:,:])) / \
                     (4.0 * self._diffcof[:-1,:,:,:] * \
                     (1.0 + ref_albedo[:-1,:,:,:]) + \
                     (1.0 - ref_albedo[:-1,:,:,:]) * self._hxyz[0]),
                     (2.0 * self._diffcof[:-1,:,:,:] * neig_dc[:-1,:,:,:]) / \
                     (self._hxyz[0] * self._diffcof[:-1,:,:,:] + \
                     self._hxyz[0] * neig_dc[:-1,:,:,:])), 0.0)

        ref_albedo = np.divide(self._current[:,:,:,:,_CURRENTS['in_back']],
                self._current[:,:,:,:,_CURRENTS['out_back']],
                where=self._current[:,:,:,:,_CURRENTS['out_back']] > 1.0e-10,
                out=np.ones_like(self._current[:,:,:,:,_CURRENTS['out_back']]))
        adj_reflector = np.roll(self._coremap, 1, axis=1) == _CMFD_NOACCEL
        neig_dc = np.roll(self._diffcof, 1, axis=1)
        # Define neg_hxyz

        dtilde2[:,1:,:,:,2] = np.where(is_accel[:,1:,:,np.newaxis], \
            np.where(adj_reflector[:,1:,:,np.newaxis],
                     (2.0 * self._diffcof[:,1:,:,:] * \
                     (1.0 - ref_albedo[:,1:,:,:])) / \
                     (4.0 * self._diffcof[:,1:,:,:] * \
                     (1.0 + ref_albedo[:,1:,:,:]) + \
                     (1.0 - ref_albedo[:,1:,:,:]) * self._hxyz[1]),
                     (2.0 * self._diffcof[:,1:,:,:] * neig_dc[:,1:,:,:]) / \
                     (self._hxyz[1] * self._diffcof[:,1:,:,:] + \
                     self._hxyz[1] * neig_dc[:,1:,:,:])), 0.0)

        ref_albedo = np.divide(self._current[:,:,:,:,_CURRENTS['in_front']],
                self._current[:,:,:,:,_CURRENTS['out_front']],
                where=self._current[:,:,:,:,_CURRENTS['out_front']] > 1.0e-10,
                out=np.ones_like(self._current[:,:,:,:,_CURRENTS['out_front']]))
        adj_reflector = np.roll(self._coremap, -1, axis=1) == _CMFD_NOACCEL
        neig_dc = np.roll(self._diffcof, -1, axis=1)
        # Define neg_hxyz

        dtilde2[:,:-1,:,:,3] = np.where(is_accel[:,:-1,:,np.newaxis], \
            np.where(adj_reflector[:,:-1,:,np.newaxis],
                     (2.0 * self._diffcof[:,:-1,:,:] * \
                     (1.0 - ref_albedo[:,:-1,:,:])) / \
                     (4.0 * self._diffcof[:,:-1,:,:] * \
                     (1.0 + ref_albedo[:,:-1,:,:]) + \
                     (1.0 - ref_albedo[:,:-1,:,:]) * self._hxyz[1]),
                     (2.0 * self._diffcof[:,:-1,:,:] * neig_dc[:,:-1,:,:]) / \
                     (self._hxyz[1] * self._diffcof[:,:-1,:,:] + \
                     self._hxyz[1] * neig_dc[:,:-1,:,:])), 0.0)

        ref_albedo = np.divide(self._current[:,:,:,:,_CURRENTS['in_bottom']],
                self._current[:,:,:,:,_CURRENTS['out_bottom']],
                where=self._current[:,:,:,:,_CURRENTS['out_bottom']] > 1.0e-10,
                out=np.ones_like(self._current[:,:,:,:,_CURRENTS['out_bottom']]))
        adj_reflector = np.roll(self._coremap, 1, axis=2) == _CMFD_NOACCEL
        neig_dc = np.roll(self._diffcof, 1, axis=2)
        # Define neg_hxyz

        dtilde2[:,:,1:,:,4] = np.where(is_accel[:,:,1:,np.newaxis], \
            np.where(adj_reflector[:,:,1:,np.newaxis],
                     (2.0 * self._diffcof[:,:,1:,:] * \
                     (1.0 - ref_albedo[:,:,1:,:])) / \
                     (4.0 * self._diffcof[:,:,1:,:] * \
                     (1.0 + ref_albedo[:,:,1:,:]) + \
                     (1.0 - ref_albedo[:,:,1:,:]) * self._hxyz[2]),
                     (2.0 * self._diffcof[:,:,1:,:] * neig_dc[:,:,1:,:]) / \
                     (self._hxyz[2] * self._diffcof[:,:,1:,:] + \
                     self._hxyz[2] * neig_dc[:,:,1:,:])), 0.0)

        ref_albedo = np.divide(self._current[:,:,:,:,_CURRENTS['in_top']],
                self._current[:,:,:,:,_CURRENTS['out_top']],
                where=self._current[:,:,:,:,_CURRENTS['out_top']] > 1.0e-10,
                out=np.ones_like(self._current[:,:,:,:,_CURRENTS['out_top']]))
        adj_reflector = np.roll(self._coremap, -1, axis=2) == _CMFD_NOACCEL
        neig_dc = np.roll(self._diffcof, -1, axis=2)
        # Define neg_hxyz

        dtilde2[:,:,:-1,:,5] = np.where(is_accel[:,:,:-1,np.newaxis], \
            np.where(adj_reflector[:,:,:-1,np.newaxis],
                     (2.0 * self._diffcof[:,:,:-1,:] * \
                     (1.0 - ref_albedo[:,:,:-1,:])) / \
                     (4.0 * self._diffcof[:,:,:-1,:] * \
                     (1.0 + ref_albedo[:,:,:-1,:]) + \
                     (1.0 - ref_albedo[:,:,:-1,:]) * self._hxyz[2]),
                     (2.0 * self._diffcof[:,:,:-1,:] * neig_dc[:,:,:-1,:]) / \
                     (self._hxyz[2] * self._diffcof[:,:,:-1,:] + \
                     self._hxyz[2] * neig_dc[:,:,:-1,:])), 0.0)

    def _compute_dtilde(self):
        # Get maximum of spatial and group indices
        nx = self._indices[0]
        ny = self._indices[1]
        nz = self._indices[2]
        ng = self._indices[3]

        # Create single vector of these indices for boundary calculation
        nxyz = np.array([[0,nx-1], [0,ny-1], [0,nz-1]])

        # Get boundary condition information
        albedo = self._albedo

        # Loop over group and spatial indices
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    for g in range(ng):
                        # Cycle if non-accelration region
                        if self._coremap[i,j,k] == _CMFD_NOACCEL:
                            continue

                        # Get cell data
                        cell_dc = self._diffcof[i,j,k,g]
                        cell_hxyz = self._hxyz

                        # Setup of vector to identify boundary conditions
                        bound = np.repeat([i,j,k], 2)

                        # Begin loop around sides of cell for leakage
                        for l in range(6):
                            xyz_idx = int(l/2)  # x=0, y=1, z=2
                            dir_idx = l % 2     # -=0, +=1

                            # Check if at boundary
                            if bound[l] == nxyz[xyz_idx, dir_idx]:
                                # Compute dtilde with albedo boundary condition
                                dtilde = (2*cell_dc*(1-albedo[l]))/ \
                                    (4*cell_dc*(1+albedo[l]) + \
                                    (1-albedo[l])*cell_hxyz[xyz_idx])

                                # Check for zero flux albedo
                                if abs(albedo[l] - _ZERO_FLUX) < _TINY_BIT:
                                    dtilde = 2*cell_dc / cell_hxyz[xyz_idx]

                            else:  # Not at a boundary
                                shift_idx = 2*(l % 2) - 1 # shift neig by -1 or +1

                                # Compute neighboring cell indices
                                neig_idx = [i,j,k]  # Begin with i,j,k
                                neig_idx[xyz_idx] += shift_idx

                                # Get neighbor cell data
                                neig_dc = self._diffcof[tuple(neig_idx) + (g,)]

                                # Check for fuel-reflector interface
                                if (self._coremap[tuple(neig_idx)] ==
                                         _CMFD_NOACCEL):
                                    # Get albedo
                                    ref_albedo = self._get_reflector_albedo(l,g,i,j,k)
                                    dtilde = (2*cell_dc*(1-ref_albedo))/(4*cell_dc*(1+ \
                                             ref_albedo)+(1-ref_albedo)*cell_hxyz[xyz_idx])

                                else:  # Not next to a reflector
                                    # Compute dtilde to neighbor cell
                                    dtilde = (2*cell_dc*neig_dc)/(cell_hxyz[xyz_idx]*cell_dc + \
                                             cell_hxyz[xyz_idx]*neig_dc)

                            # Record dtilde
                            self._dtilde[i, j, k, g, l] = dtilde

    def _compute_dhat2(self):
        # TODO Write comments for this function
        dhat2 = np.zeros(self._dhat.shape)

        net_current_minusx = ((self._current[:,:,:,:,_CURRENTS['in_left']] - \
                             self._current[:,:,:,:,_CURRENTS['out_left']]) / \
                             np.prod(self._hxyz)*self._hxyz[0])
        net_current_plusx = ((self._current[:,:,:,:,_CURRENTS['out_right']] - \
                             self._current[:,:,:,:,_CURRENTS['in_right']]) / \
                             np.prod(self._hxyz)*self._hxyz[0])
        net_current_minusy = ((self._current[:,:,:,:,_CURRENTS['in_back']] - \
                             self._current[:,:,:,:,_CURRENTS['out_back']]) / \
                             np.prod(self._hxyz)*self._hxyz[1])
        net_current_plusy = ((self._current[:,:,:,:,_CURRENTS['out_front']] - \
                             self._current[:,:,:,:,_CURRENTS['in_front']]) / \
                             np.prod(self._hxyz)*self._hxyz[1])
        net_current_minusz = ((self._current[:,:,:,:,_CURRENTS['in_bottom']] - \
                             self._current[:,:,:,:,_CURRENTS['out_bottom']]) / \
                             np.prod(self._hxyz)*self._hxyz[2])
        net_current_plusz = ((self._current[:,:,:,:,_CURRENTS['out_top']] - \
                             self._current[:,:,:,:,_CURRENTS['in_top']]) / \
                             np.prod(self._hxyz)*self._hxyz[2])

        cell_flux = self._flux / np.prod(self._hxyz)
        is_accel = self._coremap != _CMFD_NOACCEL

        dhat2[0,:,:,:,0] = np.where(is_accel[0,:,:,np.newaxis],
            (net_current_minusx[0,:,:,:] + self._dtilde[0,:,:,:,0] * \
            cell_flux[0,:,:,:]) / cell_flux[0,:,:,:], 0)
        dhat2[-1,:,:,:,1] = np.where(is_accel[-1,:,:,np.newaxis],
            (net_current_plusx[-1,:,:,:] - self._dtilde[-1,:,:,:,1] * \
            cell_flux[-1,:,:,:]) / cell_flux[-1,:,:,:], 0)
        dhat2[:,0,:,:,2] = np.where(is_accel[:,0,:,np.newaxis],
            (net_current_minusy[:,0,:,:] + self._dtilde[:,0,:,:,2] * \
            cell_flux[:,0,:,:]) / cell_flux[:,0,:,:], 0)
        dhat2[:,-1,:,:,3] = np.where(is_accel[:,-1,:,np.newaxis],
            (net_current_plusy[:,-1,:,:] + self._dtilde[:,-1,:,:,3] * \
            cell_flux[:,-1,:,:]) / cell_flux[:,-1,:,:], 0)
        dhat2[:,:,0,:,4] = np.where(is_accel[:,:,0,np.newaxis],
            (net_current_minusz[:,:,0,:] + self._dtilde[:,:,0,:,4] * \
            cell_flux[:,:,0,:]) / cell_flux[:,:,0,:], 0)
        dhat2[:,:,-1,:,5] = np.where(is_accel[:,:,-1,np.newaxis],
            (net_current_minusz[:,:,-1,:] + self._dtilde[:,:,-1,:,5] * \
            cell_flux[:,:,-1,:]) / cell_flux[:,:,-1,:], 0)

        # Minus x direction
        adj_reflector = np.roll(self._coremap, 1, axis=0) == _CMFD_NOACCEL
        neig_flux = np.roll(self._flux, 1, axis=0) / np.prod(self._hxyz)
        dhat2[1:,:,:,:,0] = np.where(is_accel[1:,:,:,np.newaxis], \
            np.where(adj_reflector[1:,:,:,np.newaxis],
                     (net_current_minusx[1:,:,:,:] + self._dtilde[1:,:,:,:,0] * \
                     cell_flux[1:,:,:,:]) / cell_flux[1:,:,:,:],
                     (net_current_minusx[1:,:,:,:] - self._dtilde[1:,:,:,:,0] * \
                     (neig_flux[1:,:,:,:] - cell_flux[1:,:,:,:])) / \
                     (neig_flux[1:,:,:,:] + cell_flux[1:,:,:,:])), 0.0)

        # Plus x direction
        adj_reflector = np.roll(self._coremap, -1, axis=0) == _CMFD_NOACCEL
        neig_flux = np.roll(self._flux, -1, axis=0) / np.prod(self._hxyz)
        dhat2[:-1,:,:,:,1] = np.where(is_accel[:-1,:,:,np.newaxis], \
            np.where(adj_reflector[:-1,:,:,np.newaxis],
                     (net_current_plusx[:-1,:,:,:] - self._dtilde[:-1,:,:,:,1] * \
                     cell_flux[:-1,:,:,:]) / cell_flux[:-1,:,:,:],
                     (net_current_plusx[:-1,:,:,:] + self._dtilde[:-1,:,:,:,1] * \
                     (neig_flux[:-1,:,:,:] - cell_flux[:-1,:,:,:])) / \
                     (neig_flux[:-1,:,:,:] + cell_flux[:-1,:,:,:])), 0.0)

        # Minus y direction
        adj_reflector = np.roll(self._coremap, 1, axis=1) == _CMFD_NOACCEL
        neig_flux = np.roll(self._flux, 1, axis=1) / np.prod(self._hxyz)
        dhat2[:,1:,:,:,2] = np.where(is_accel[:,1:,:,np.newaxis], \
            np.where(adj_reflector[:,1:,:,np.newaxis],
                     (net_current_minusy[:,1:,:,:] + self._dtilde[:,1:,:,:,2] * \
                     cell_flux[:,1:,:,:]) / cell_flux[:,1:,:,:],
                     (net_current_minusy[:,1:,:,:] - self._dtilde[:,1:,:,:,2] * \
                     (neig_flux[:,1:,:,:] - cell_flux[:,1:,:,:])) / \
                     (neig_flux[:,1:,:,:] + cell_flux[:,1:,:,:])), 0.0)

        # Plus y direction
        adj_reflector = np.roll(self._coremap, -1, axis=1) == _CMFD_NOACCEL
        neig_flux = np.roll(self._flux, -1, axis=1) / np.prod(self._hxyz)
        dhat2[:,:-1,:,:,3] = np.where(is_accel[:,:-1,:,np.newaxis], \
            np.where(adj_reflector[:,:-1,:,np.newaxis],
                     (net_current_plusy[:,:-1,:,:] - self._dtilde[:,:-1,:,:,3] * \
                     cell_flux[:,:-1,:,:]) / cell_flux[:,:-1,:,:],
                     (net_current_plusy[:,:-1,:,:] + self._dtilde[:,:-1,:,:,3] * \
                     (neig_flux[:,:-1,:,:] - cell_flux[:,:-1,:,:])) / \
                     (neig_flux[:,:-1,:,:] + cell_flux[:,:-1,:,:])), 0.0)

        # Minus z direction
        adj_reflector = np.roll(self._coremap, 1, axis=2) == _CMFD_NOACCEL
        neig_flux = np.roll(self._flux, 1, axis=2) / np.prod(self._hxyz)
        dhat2[:,:,1:,:,4] = np.where(is_accel[:,:,1:,np.newaxis], \
            np.where(adj_reflector[:,:,1:,np.newaxis],
                     (net_current_minusz[:,:,1:,:] + self._dtilde[:,:,1:,:,4] * \
                     cell_flux[:,:,1:,:]) / cell_flux[:,:,1:,:],
                     (net_current_minusz[:,:,1:,:] - self._dtilde[:,:,1:,:,4] * \
                     (neig_flux[:,:,1:,:] - cell_flux[:,:,1:,:])) / \
                     (neig_flux[:,:,1:,:] + cell_flux[:,:,1:,:])), 0.0)

        # Plus z direction
        adj_reflector = np.roll(self._coremap, -1, axis=2) == _CMFD_NOACCEL
        neig_flux = np.roll(self._flux, -1, axis=2) / np.prod(self._hxyz)
        dhat2[:,:,:-1,:,5] = np.where(is_accel[:,:,:-1,np.newaxis], \
            np.where(adj_reflector[:,:,:-1,np.newaxis],
                     (net_current_plusz[:,:,:-1,:] - self._dtilde[:,:,:-1,:,5] * \
                     cell_flux[:,:,:-1,:]) / cell_flux[:,:,:-1,:],
                     (net_current_plusz[:,:,:-1,:] + self._dtilde[:,:,:-1,:,5] * \
                     (neig_flux[:,:,:-1,:] - cell_flux[:,:,:-1,:])) / \
                     (neig_flux[:,:,:-1,:] + cell_flux[:,:,:-1,:])), 0.0)

    def _compute_dhat(self):
        #TODO compute dhat and dtilde for general case with hxyz (just define as repeated but use in formulas)
        # Get maximum of spatial and group indices
        nx = self._indices[0]
        ny = self._indices[1]
        nz = self._indices[2]
        ng = self._indices[3]

        # Create single vector of these indices for boundary calculation
        nxyz = np.array([[0,nx-1], [0,ny-1], [0,nz-1]])

        # Loop over group and spatial indices
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    for g in range(ng):
                        # Cycle if non-accelration region
                        if self._coremap[i,j,k] == _CMFD_NOACCEL:
                            continue

                        # Get cell data
                        cell_dtilde = self._dtilde[i,j,k,g,:]
                        cell_flux = self._flux[i,j,k,g]/np.product(self._hxyz)
                        current = self._current[i,j,k,g,:]

                        # Setup of vector to identify boundary conditions
                        bound = np.repeat([i,j,k], 2)

                        # Begin loop around sides of cell for leakage
                        for l in range(6):
                            xyz_idx = int(l/2)  # x=0, y=1, z=2
                            dir_idx = l % 2     # -=0, +=1
                            shift_idx = 2*(l % 2) - 1 # shift neig by -1 or +1

                            # Calculate net current on l face (divided by surf area)
                            net_current = shift_idx*(current[2*l] - current[2*l+1]) / \
                                np.product(self._hxyz) * self._hxyz[xyz_idx]

                            # Check if at boundary
                            if bound[l] == nxyz[xyz_idx, dir_idx]:
                                # Compute dhat
                                dhat = (net_current - shift_idx*cell_dtilde[l]*cell_flux) / \
                                     cell_flux
                                if l == 1:
                                    print(dhat, i, j, k, g, l)
                                    print(net_current, cell_dtilde, cell_flux)
                                    print("yo")

                            else:  # Not at a boundary
                                # Compute neighboring cell indices
                                neig_idx = [i,j,k]  # Begin with i,j,k
                                neig_idx[xyz_idx] += shift_idx

                                # Get neigbor flux
                                neig_flux = self._flux[tuple(neig_idx)+(g,)] / \
                                            np.product(self._hxyz)

                                # Check for fuel-reflector interface
                                if (self._coremap[tuple(neig_idx)] ==
                                         _CMFD_NOACCEL):
                                    # Compute dhat
                                    dhat = (net_current - shift_idx*cell_dtilde[l]*cell_flux) / \
                                         cell_flux
                                    #if l==0:
                                    #    print("hit")
                                    #    print(dhat, net_current, cell_dtilde[l], cell_flux)
                                else:  # not a fuel-reflector interface
                                    # Compute dhat
                                    dhat = (net_current + shift_idx*cell_dtilde[l]* \
                                        (neig_flux - cell_flux))/(neig_flux + cell_flux)

                            # Record dhat
                            self._dhat[i, j, k, g, l] = dhat

                            # check for dhat reset
                            if self._dhat_reset:
                                self._dhat[i, j, k, g, l] = 0.0

        # Write that dhats are zero
        if self._dhat_reset and openmc.capi.settings.verbosity >= 8 and \
                openmc.capi.settings.master:
            print(' Dhats reset to zero')

    def _get_reflector_albedo(self, l, g, i, j, k):
        # Get partial currents from object
        current = self._current[i,j,k,g,:]

        # Calculate albedo
        if current[2*l] < 1.0e-10:
            return 1.0
        else:
          return current[2*l+1]/current[2*l]

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
        self._cmfd_tally_ids = []
        for i in range(n_tallies):
            tally = openmc.capi.Tally()
            # Set nuclide bins
            tally.nuclides = ['total']
            self._cmfd_tally_ids.append(tally.id)

            # Set attributes of CMFD flux, total tally
            if i == 0:
                # Set filters for tally
                if self._energy_filters:
                    tally.filters = [mesh_filter, energy_filter]
                else:
                    tally.filters = [mesh_filter]
                # Set scores, type, and estimator for tally
                tally.scores = ['flux', 'total']
                tally.type = 'volume'
                tally.estimator = 'analog'

            # Set attributes of CMFD neutron production tally
            elif i == 1:
                # Set filters for tally
                if self._energy_filters:
                    tally.filters = [mesh_filter, energy_filter, energyout_filter]
                else:
                    tally.filters = [mesh_filter]
                # Set scores, type, and estimator for tally
                tally.scores = ['nu-scatter', 'nu-fission']
                tally.type = 'volume'
                tally.estimator = 'analog'

            # Set attributes of CMFD surface current tally
            elif i == 2:
                # Set filters for tally
                if self._energy_filters:
                    tally.filters = [meshsurface_filter, energy_filter]
                else:
                    tally.filters = [meshsurface_filter]
                # Set scores, type, and estimator for tally
                tally.scores = ['current']
                tally.type = 'mesh-surface'
                tally.estimator = 'analog'

            # Set attributes of CMFD P1 scatter tally
            elif i == 3:
                # Set filters for tally
                if self._energy_filters:
                    tally.filters = [mesh_filter, legendre_filter, energy_filter]
                else:
                    tally.filters = [mesh_filter, legendre_filter]
                # Set scores for tally
                tally.scores = ['scatter']
                tally.type = 'volume'
                tally.estimator = 'analog'

            # Set all tallies to be active from beginning
            tally.active = True