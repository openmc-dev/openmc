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
import sys
import numpy as np
# Line below is added to suppress warnings when using numpy.divide to
# divide by numpy arrays that contain entries with zero
np.seterr(divide='ignore', invalid='ignore')
from scipy import sparse
import time
# See if mpi4py module can be imported, define have_mpi global variable
try:
    from mpi4py import MPI
    have_mpi = True
except ImportError:
    have_mpi = False

import openmc.capi
from openmc._xml import clean_indentation  #TODO Remove
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

# Map that returns index of current direction in numpy current matrix
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
        to any tallies far away from fission source neutron regions.  A ``1``
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

    # REMOVE this method
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


# REMOVE this entire class
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
        clean_indentation(self._cmfd_file)

        # Write the XML Tree to the cmfd.xml file
        tree = ET.ElementTree(self._cmfd_file)
        tree.write("cmfd.xml", xml_declaration=True,
                   encoding='utf-8', method="xml")


class CMFDRun(object):
    r"""Class to run openmc with CMFD acceleration through the C API. Running
    openmc in this manner obviates the need for defining CMFD parameters
    through a cmfd.xml file. Instead, all input parameters should be passed
    through the CMFDRun instance.

    Attributes
    ----------
    cmfd_begin : int
        Batch number at which CMFD calculations should begin
    dhat_reset : bool
        Indicate whether :math:`\widehat{D}` nonlinear CMFD parameters should be
        reset to zero before solving CMFD eigenproblem.
    cmfd_display : {'balance', 'dominance', 'entropy', 'source'}
        Set one additional CMFD output column. Options are:

          * "balance" - prints the RMS [%] of the resdiual from the neutron balance
             equation on CMFD tallies.
          * "dominance" - prints the estimated dominance ratio from the CMFD
            iterations.
          * "entropy" - prints the *entropy* of the CMFD predicted fission source.
          * "source" - prints the RMS [%] between the OpenMC fission source and
            CMFD fission source.
    cmfd_downscatter : bool
        Indicate whether an effective downscatter cross section should be used
        when using 2-group CMFD.
    cmfd_feedback : bool
        Indicate or not the CMFD diffusion result is used to adjust the weight
        of fission source neutrons on the next OpenMC batch. Defaults to False.
    cmfd_ktol : float
        Tolerance on the eigenvalue when performing CMFD power iteration
    cmfd_mesh : openmc.CMFDMesh
        Structured mesh to be used for acceleration
    norm : float
        Normalization factor applied to the CMFD fission source distribution
    cmfd_power_monitor : bool
        View convergence of power iteration during CMFD acceleration
    cmfd_run_adjoint : bool
        Perform adjoint calculation on the last batch
    cmfd_shift : float
        Optional Wielandt shift parameter for accelerating power iterations. By
        default, it is very large so there is effectively no impact.
    cmfd_stol : float
        Tolerance on the fission source when performing CMFD power iteration
    cmfd_reset : list of int
        List of batch numbers at which CMFD tallies should be reset
    cmfd_write_matrices : bool
        # TODO update this to read "resultant normalized flux vector"
        Write sparse matrices that are used during CMFD acceleration (loss,
        production) and resultant flux vector phi to file
    indices : numpy.ndarray
        Stores spatial and group dimensions as [nx, ny, nz, ng]
    egrid : numpy.ndarray
        Energy grid used for CMFD acceleration
    albedo : numpy.ndarray
        Albedo for global boundary conditions, taken from CMFD mesh. It is
        listed in the following order: [-x +x -y +y -z +z]. Set to
        [1, 1, 1, 1, 1, 1] if not specified by user.
    coremap : numpy.ndarray
        Coremap for coarse mesh overlay, defining acceleration regions in CMFD
        problem, created from map variable in CMFD class. Each accelerated
        region is given a unique index to map each spatial mesh to a
        corresponding row within the CMFD matrices. Non-accelerated regions are
        set with value _CMFD_NOACCEL. If a user does not specify a CMFD map,
        coremap is set to treat each spatial mesh as an accelerated region.
        For non-accelerated regions on coarse mesh overlay
    n_cmfd_resets : int
        Number of elements in cmfd_reset to store number of times cmfd tallies
        will be reset.
    cmfd_mesh_id : int
        Mesh id of CMFD mesh, stored to access CMFD mesh object in memory
    cmfd_tally_ids : list of ints
        List that stores the id's of all CMFD tallies. The tallies that each
        id corresponds to are:

          * cmfd_tally_id[0] - CMFD flux, total tally
          * cmfd_tally_id[1] - CMFD neutron production tally
          * cmfd_tally_id[2] - CMFD surface current tally
          * cmfd_tally_id[3] - CMFD P1 scatter tally
    energy_filters : bool
        Stores whether energy filters should be created or not, based on whether
        a user specifies an energy grid
    cmfd_on : bool
        Stores whether cmfd solver should be invoked, based on whether current
        batch has reached variable ``cmfd_begin``
    mat_dim : int
        Number of accelerated regions exist in problem. Value of ``mat_dim`` is
        tied to size of CMFD matrices
    keff_bal : float
        Balance k-effective computed from OpenMC source
    cmfd_adjoint_type : {'physical', 'math'}
        Stores type of adjoint calculation that should be performed.
        ``cmfd_run_adjoint`` must be true for an adjoint calculation to be
        perfomed. Options are:

          * "physical" - Create adjoint matrices from physical parameters of
             CMFD problem
          * "math" - Create adjoint matrices mathematically as the transpose of
             loss and production CMFD matrices
    keff : float
        K-effective from solving CMFD matrix equations
    adj_keff : float
        K-effective from solving adjoint CMFD matrix equations
    phi : numpy.ndarray
        Final flux vector from solving CMFD matrix equations
    adj_phi : numpy.ndarray
        Final flux vector from solving adjoint CMFD matrix equations
    flux : numpy.ndarray
        Flux computed from tally data
    totalxs : numpy.ndarray
        Total cross section computed from tally data
    p1scattxs : numpy.ndarray
        P1 scattering cross section computed from tally data
    scattxs : numpy.ndarray
        Scattering cross section computed from tally data
    nfissxs : numpy.ndarray
        Nu-fission cross section computed from tally data
    diffcof : numpy.ndarray
        Diffusion coefficients computed from tally data
    dtilde : numpy.ndarray
        Array of diffusion coupling coefficients
    dhat : numpy.ndarray
        Array of nonlinear coupling coefficients
    hxyz : numpy.ndarray
        Dimensions of mesh cells, stored as (xloc,yloc,zloc,[hu,hv,hw])
    current : numpy.ndarray
        Surface currents for each mesh cell in each incoming and outgoing
        direction on each mesh surface, computed from tally data
    cmfd_src : numpy.ndarray
        CMFD source distribution calculated from solving CMFD equations
    openmc_src : numpy.ndarray
        OpenMC source distribution computed from tally data
    sourcecounts : numpy.ndarray
        Number of neutrons in each spatial and energy group, used to calculate
        normalizing factors and weight factors when reweighting each neutron
    weightfactors : numpy.ndarray
        Weight factors for each spatial and energy group, used to reweight
        neutrons at the end of each batch
    entropy : list of floats
        "Shannon entropy" from cmfd fission source, stored for each generation
        that CMFD is invoked
    balance : list of floats
        RMS of neutron balance equations, stored for each generation that CMFD is
        invoked
    src_cmp : list of floats
        RMS deviation of OpenMC and CMFD normalized source, stored for each
        generation that CMFD is invoked
    dom : list of floats
        Dominance ratio from solving CMFD matrix equations, stored for each
        generation that CMFD is invoked
    k_cmfd : list of floats
        List of CMFD k-effectives, stored for each generation that CMFD is
        invoked
    resnb : numpy.ndarray
        Residual from solving neutron balance equations
    time_cmfd : float
        Time for entire CMFD calculation, in seconds
    time_cmfdbuild : float
        Time for building CMFD matrices, in seconds
    time_cmfdsolve : float
        Time for solving CMFD matrix equations, in seconds
    intracomm : mpi4py.MPI.Intracomm or None
        MPI intercommunicator for running MPI commands

    TODO Remove CMFD constants, timing variables in Fortran
    TODO Remove cmfd_data.F90, cmfd_execute.F90, cmfd_header.F90,
         cmfd_input.F90, cmfd_loss_operator.F90, cmfd_prod_operator.F90,
         cmfd_solver.F90
    TODO Remove instances of cmfd in sourcepoint.F90, output.F90,
         simulation.F90, api.F90, settings.F90, input_xml.F90, plots.F90

    """

    def __init__(self):
        """Constructor for CMFDRun class. Default values for instance variables
        set in this method.

        """
        # Variables that users can modify
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

        # External variables used during runtime but users cannot control
        self._indices = np.zeros(4, dtype=int)
        self._egrid = None
        self._albedo = None
        self._coremap = None
        self._n_cmfd_resets = 0
        self._cmfd_mesh_id = None
        self._cmfd_tally_ids = None
        self._energy_filters = None
        self._cmfd_on = False
        self._mat_dim = _CMFD_NOACCEL
        self._keff_bal = None
        self._cmfd_adjoint_type = 'physical'
        self._keff = None
        self._adj_keff = None
        self._phi = None
        self._adj_phi = None
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

        # Add all index-related variables, for numpy vectorization
        self._last_x_accel = None
        self._first_y_accel = None
        self._last_y_accel = None
        self._first_z_accel = None
        self._last_z_accel = None
        self._notfirst_x_accel = None
        self._notlast_x_accel = None
        self._notfirst_y_accel = None
        self._notlast_y_accel = None
        self._notfirst_z_accel = None
        self._notlast_z_accel = None
        self._adj_reflector_left = None
        self._adj_reflector_right = None
        self._adj_reflector_back = None
        self._adj_reflector_front = None
        self._adj_reflector_bottom = None
        self._adj_reflector_top = None
        self._accel_idxs = None
        self._accel_neig_left_idxs = None
        self._accel_neig_right_idxs = None
        self._accel_neig_back_idxs = None
        self._accel_neig_front_idxs = None
        self._accel_neig_bottom_idxs = None
        self._accel_neig_top_idxs = None
        self._loss_row = None
        self._loss_col = None
        self._intracomm = None


    @property
    def cmfd_begin(self):
        return self._cmfd_begin

    @property
    def dhat_reset(self):
        return self._dhat_reset

    @property
    def cmfd_display(self):
        return self._cmfd_display

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

    @cmfd_display.setter
    def cmfd_display(self, display):
        check_type('CMFD display', display, str)
        check_value('CMFD display', display,
                    ['balance', 'dominance', 'entropy', 'source'])
        self._cmfd_display = display

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

    # TEMP
    def run(self, omp_num_threads=None, intracomm=None, vectorized=True, updated=False):
        """Public method to run OpenMC with CMFD

        This method is called by user to run CMFD once instance variables of
        CMFDRun class are set

        Parameters
        ----------
        omp_num_threads : int
            Number of OpenMP threads to use for OpenMC simulation
        intracomm : mpi4py.MPI.Intracomm or None
            MPI intercommunicator to pass through C API. Set to MPI.COMM_WORLD
            by default

        """
        # Store intracomm for part of CMFD routine where MPI reduce and
        # broadcast calls are made
        if intracomm is not None:
            self._intracomm = intracomm
        elif intracomm is None and have_mpi:
            self._intracomm = MPI.COMM_WORLD

        # TEMP
        self._updated = updated

        # Check number of OpenMP threads is valid input and initialize C API
        if omp_num_threads is not None:
            check_type('OpenMP num threads', omp_num_threads, Integral)
            openmc.capi.init(args=['-s',str(omp_num_threads)])
        else:
            openmc.capi.init()

        # Configure CMFD parameters and tallies
        self._configure_cmfd()

        # Set up CMFD coremap
        self._set_coremap()

        # Compute and store array indices used to build cross section arrays
        self._precompute_array_indices()

        # Compute and store row and column indices used to build CMFD matrices
        self._precompute_matrix_indices()

        # Initialize simulation
        openmc.capi.simulation_init()

        while(True):
            # Initialize CMFD batch
            self._cmfd_init_batch()

            # Run next batch
            status = openmc.capi.next_batch()

            # Perform CMFD calculation if on
            if self._cmfd_on:
                self._execute_cmfd(vectorized)

                # Write CMFD output if CMFD on for current batch
                if openmc.capi.master():
                    self._write_cmfd_output()

            if status != 0:
                break

        # Finalize simuation
        openmc.capi.simulation_finalize()

        # Print out CMFD timing statistics
        self._write_cmfd_timing_stats()

        # Finalize and free memory
        openmc.capi.finalize()

    def _write_cmfd_output(self):
        """Write CMFD output to buffer at the end of each batch"""
        # Display CMFD k-effective
        str1 = 'CMFD k:   {:0.5f}'.format(self._k_cmfd[-1])
        # Display value of additional field based on value of cmfd_display
        if self._cmfd_display == 'dominance':
            str2 = 'Dom Rat:  {:0.5f}'.format(self._dom[-1])
        elif self._cmfd_display == 'entropy':
            str2 = 'CMFD Ent: {:0.5f}'.format(self._entropy[-1])
        elif self._cmfd_display == 'source':
            str2 = 'RMS Src:  {:0.5f}'.format(self._src_cmp[-1])
        else:
            str2 = 'RMS Bal:  {:0.5f}'.format(self._balance[-1])

        print('{0:>76s}\n{1:>76s}'.format(str1, str2))
        sys.stdout.flush()

    def _write_cmfd_timing_stats(self):
        """Write CMFD timing statistics to buffer after finalizing simulation"""
        if openmc.capi.master():
          print(
"""=====================>     CMFD TIMING STATISTICS     <====================

   Time in CMFD                    =  {0:.5E} seconds
     Building matrices             =  {1:.5E} seconds
     Solving matrices              =  {2:.5E} seconds
""".format(self._time_cmfd, self._time_cmfdbuild, self._time_cmfdsolve))
          sys.stdout.flush()

    def _configure_cmfd(self):
        """Initialize CMFD parameters and set CMFD input variables"""
        # Read in cmfd input defined in Python
        self._read_cmfd_input()

        # Initialize timers
        self._time_cmfd = 0.0
        self._time_cmfdbuild = 0.0
        self._time_cmfdsolve = 0.0

        # Initialize all numpy arrays used for CMFD solver
        self._allocate_cmfd()

    def _read_cmfd_input(self):
        """Sets values of additional instance variables based on user input"""
        # Print message to user and flush output to stdout
        if openmc.capi.settings.verbosity >= 7 and openmc.capi.master():
            print(' Configuring CMFD parameters for simulation')
            sys.stdout.flush()

        # Check if CMFD mesh is defined
        if self._cmfd_mesh is None:
            raise ValueError('No CMFD mesh has been specified for '
                             'simulation')

        # Set spatial dimensions of CMFD object
        for i, n in enumerate(self._cmfd_mesh.dimension):
            self._indices[i] = n

        # Check if in continuous energy mode
        if not openmc.capi.settings.run_CE:
            raise OpenMCError('CMFD must be run in continuous energy mode')

        # Set number of energy groups
        if self._cmfd_mesh.energy is not None:
            ng = len(self._cmfd_mesh.energy)
            self._egrid = np.array(self._cmfd_mesh.energy)
            self._indices[3] = ng - 1
            self._energy_filters = True
        else:
            self._egrid = np.array([_ENERGY_MIN_NEUTRON, _ENERGY_MAX_NEUTRON])
            self._indices[3] = 1
            self._energy_filters = False

        # Set global albedo
        if self._cmfd_mesh.albedo is not None:
            self._albedo = np.array(self._cmfd_mesh.albedo)
        else:
            self._albedo = np.array([1., 1., 1., 1., 1., 1.])

        # Get acceleration map, otherwise set all regions to be accelerated
        if self._cmfd_mesh.map is not None:
            check_length('CMFD coremap', self._cmfd_mesh.map,
                         np.product(self._indices[0:3]))
            self._coremap = np.array(self._cmfd_mesh.map)
        else:
            self._coremap = np.ones((np.product(self._indices[0:3])), dtype=int)

        # Set number of batches where cmfd tallies should be reset
        if self._cmfd_reset is not None:
            self._n_cmfd_resets = len(self._cmfd_reset)

        # Create tally objects
        self._create_cmfd_tally()

    def _allocate_cmfd(self):
        """Allocates all numpy arrays and lists used in CMFD algorithm"""
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

        # Allocate dimensions for each mesh cell
        self._hxyz = np.zeros((nx, ny, nz, 3))

        # Allocate surface currents
        self._current = np.zeros((nx, ny, nz, 12, ng))

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
        """Handles CMFD options at the beginning of each batch"""
        # Get current batch through C API
        # Add 1 as next_batch has not been called yet
        current_batch = openmc.capi.current_batch() + 1

        # Check to activate CMFD diffusion and possible feedback
        if self._cmfd_begin == current_batch:
            self._cmfd_on = True

        # Check to reset tallies
        if self._n_cmfd_resets > 0 and current_batch in self._cmfd_reset:
            self._cmfd_tally_reset()

    def _execute_cmfd(self, vectorized):
        """Runs CMFD calculation on master node"""
        # Run CMFD on single processor on master
        if openmc.capi.master():
            # Start CMFD timer
            time_start_cmfd = time.time()

            # Create CMFD data from OpenMC tallies
            self._set_up_cmfd(vectorized=vectorized)

            # Call solver
            self._cmfd_solver_execute(vectorized=vectorized)

            # Store k-effective
            self._k_cmfd.append(self._keff)

            # Check to perform adjoint on last batch
            if (openmc.capi.current_batch() == openmc.capi.settings.batches
                and self._cmfd_run_adjoint):
                self._cmfd_solver_execute(adjoint=True)

            # Calculate fission source
            self._calc_fission_source(vectorized=vectorized)

        # Calculate weight factors
        self._cmfd_reweight(True)

        # Stop CMFD timer
        if openmc.capi.master():
            time_stop_cmfd = time.time()
            self._time_cmfd += time_stop_cmfd - time_start_cmfd

    def _cmfd_tally_reset(self):
        """Resets all CMFD tallies in memory"""
        # Print message
        if openmc.capi.settings.verbosity >= 6 and openmc.capi.master():
            print(' CMFD tallies reset')
            sys.stdout.flush()

        # Reset CMFD tallies
        tallies = openmc.capi.tallies
        for tally_id in self._cmfd_tally_ids:
            tallies[tally_id].reset()

    def _set_up_cmfd(self, vectorized=True):
        """Configures CMFD object for a CMFD eigenvalue calculation

        Parameters
        ----------
        vectorized : bool
            Whether to compute dhat and dtilde using a vectorized numpy approach
            or with traditional for loops

        """
        # Calculate all cross sections based on reaction rates from last batch
        self._compute_xs()

        # Compute effective downscatter cross section
        if (self._cmfd_downscatter): self._compute_effective_downscatter()

        # Check neutron balance
        self._neutron_balance()

        # Calculate dtilde
        if vectorized:
            # TEMP
            if self._updated:
                self._compute_dtilde_vectorized_updated()
            else:
                self._compute_dtilde_vectorized()
        else:
            self._compute_dtilde()

        # Calculate dhat
        if vectorized:
            if self._updated:
                self._compute_dhat_vectorized_updated()
            else:
                self._compute_dhat_vectorized()
        else:
            self._compute_dhat()

    def _cmfd_solver_execute(self, adjoint=False, vectorized=True):
        """Sets up and runs power iteration solver for CMFD

        Parameters
        ----------
        adjoint : bool
            Whether or not to run an adjoint calculation
        vectorized : bool
            Whether to build CMFD matrices using a vectorized numpy approach
            or with traditional for loops

        """
        # Check for physical adjoint
        physical_adjoint = adjoint and self._cmfd_adjoint_type == 'physical'

        # Start timer for build
        time_start_buildcmfd = time.time()

        # Build loss and production matrices
        loss, prod = self._build_matrices(physical_adjoint, vectorized)

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

        # Write out flux vector
        if self._cmfd_write_matrices:
            if adjoint:
                # TODO Change to self._adj_phi
                self._write_vector(phi, 'adj_fluxvec')
            else:
                # TODO Change to self._phi
                self._write_vector(phi, 'fluxvec')

    def _write_vector(self, vector, base_filename):
        """Write a 1-D numpy array to file and also save it in .npy format. This
        particular format allows users to load the variable directly in a Python
        session with np.load()

        Parameters
        ----------
        vector : numpy.ndarray
            Vector that will be saved
        base_filename : string
            Filename to save vector as, without any file extension at the end.
            Vector will be saved to file [base_filename].dat and in numpy format
            as [base_filename].npy

        """
        # Write each element in vector to file
        with open(base_filename+'.dat', 'w') as fh:
            for val in vector:
                fh.write('{:0.8f}\n'.format(val))

        # Save as numpy format
        np.save(base_filename, vector)

    def _write_matrix(self, matrix, base_filename):
        """Write a numpy matrix to file and also save it in .npz format. This
        particular format allows users to load the variable directly in a Python
        session with scipy.sparse.load_npz()

        Parameters
        ----------
        matrix : scipy.sparse.spmatrix
            Sparse matrix that will be saved
        base_filename : string
            Filename to save matrix entries, without any file extension at the
            end. Matrix entries will be saved to file [base_filename].dat and in
            scipy format as [base_filename].npz

        """
        # Write row, col, and data of each entry in sparse matrix. This ignores
        # all zero-entries, and indices are written with zero-based indexing
        with open(base_filename+'.dat', 'w') as fh:
            for row in range(matrix.shape[0]):
                # Get all cols for particular row in matrix
                cols = matrix.indices[matrix.indptr[row]:matrix.indptr[row+1]]
                # Get all data entries for particular row in matrix
                data = matrix.data[matrix.indptr[row]:matrix.indptr[row+1]]
                for i in range(len(cols)):
                    fh.write('({0:3d}, {1:3d}): {2:0.8f}\n'.format(
                        row, cols[i], data[i]))

        # Save matrix in scipy format
        sparse.save_npz(base_filename, matrix)

    def _calc_fission_source(self, vectorized=True):
        """Calculates CMFD fission source from CMFD flux. If a coremap is
        defined, there will be a discrepancy between the spatial indices in the
        variables ``phi`` and ``nfissxs``, so ``phi`` needs to be mapped to the
        spatial indices of the cross sections. This can be done in a vectorized
        numpy manner or with for loops

        Parameters
        ----------
        vectorized : bool
            Whether to compute CMFD fission source in a vectorized manner or by
            looping over all spatial regions and energy groups. This distinction
            is made because ``phi`` does not correspond to the same spatial
            domain as ``nfissxs``

        Returns
        -------
        bool
            Whether or not the other filter is a subset of this filter

        """
        # Extract number of groups and number of accelerated regions
        nx = self._indices[0]
        ny = self._indices[1]
        nz = self._indices[2]
        ng = self._indices[3]
        n = self._mat_dim

        # Reset CMFD source to 0
        self._cmfd_src.fill(0.)

        # Compute cmfd_src in a vecotorized manner by phi to the spatial indices
        # of the actual problem so that cmfd_flux can be multiplied by nfissxs
        if vectorized:
            # Calculate volume
            vol = np.product(self._hxyz, axis=3)

            # Reshape phi by number of groups
            phi = self._phi.reshape((n, ng))

            # Extract indices of coremap that are accelerated
            idx = self._accel_idxs

            # Initialize CMFD flux map that maps phi to actual spatial and
            # group indices of problem
            cmfd_flux = np.zeros((nx, ny, nz, ng))

            # Loop over all groups and set CMFD flux based on indices of
            # coremap and values of phi
            for g in range(ng):
                phi_g = phi[:,g]
                # TEMP
                if self._updated:
                    cmfd_flux[idx + (g,)]  = phi_g[self._coremap[idx]]
                else:
                    cmfd_flux[idx + (np.full((n,),g),)]  = phi_g[self._coremap[idx]]

            # Compute fission source
            cmfd_src = np.sum(self._nfissxs[:,:,:,:,:] * \
                       cmfd_flux[:,:,:,:,np.newaxis], axis=3) * \
                       vol[:,:,:,np.newaxis]

        # Otherwise compute cmfd_src by looping over all groups and spatial
        # regions
        else:
            cmfd_src = np.zeros((nx, ny, nz, ng))
            for k in range(nz):
                for j in range(ny):
                    for i in range(nx):
                        for g in range(ng):
                            # Cycle through if non-accelerated region
                            if self._coremap[i,j,k] == _CMFD_NOACCEL:
                                continue

                            # Calculate volume
                            vol = np.product(self._hxyz[i,j,k])

                            # Get index in matrix
                            idx = self._indices_to_matrix(i, j, k, 0, ng)

                            # Compute fission source
                            cmfd_src[i,j,k,g] = np.sum(
                                self._nfissxs[i,j,k,:,g] \
                                * self._phi[idx:idx+ng]) * vol

        # Normalize source such that it sums to 1.0
        self._cmfd_src = cmfd_src / np.sum(cmfd_src)

        # Compute entropy
        if openmc.capi.settings.entropy_on:
            # Compute source times log_2(source)
            source = self._cmfd_src[self._cmfd_src > 0] \
                * np.log(self._cmfd_src[self._cmfd_src > 0])/np.log(2)

            # Sum source and store
            self._entropy.append(-1.0 * np.sum(source))

        # Normalize source so average is 1.0
        self._cmfd_src = self._cmfd_src/np.sum(self._cmfd_src) * self._norm

        # Calculate differences between normalized sources
        self._src_cmp.append(np.sqrt(1.0 / self._norm \
            * np.sum((self._cmfd_src - self._openmc_src)**2)))

    def _cmfd_reweight(self, new_weights):
        """Performs weighting of particles in source bank

        Parameters
        ----------
        new_weights : bool
            Whether to reweight particles or not

        """
        # Compute new weight factors
        if new_weights:

            # Get spatial dimensions and energy groups
            nx = self._indices[0]
            ny = self._indices[1]
            nz = self._indices[2]
            ng = self._indices[3]

            # Set weight factors to default 1.0
            self._weightfactors.fill(1.0)

            # Count bank site in mesh and reverse due to egrid structured
            outside = self._count_bank_sites()

            # Check and raise error if source sites exist outside of CMFD mesh
            if openmc.capi.master() and outside:
                raise OpenMCError('Source sites outside of the CMFD mesh')

            # Have master compute weight factors, ignore any zeros in
            # sourcecounts or cmfd_src
            if openmc.capi.master():
                # Compute normalization factor
                norm = np.sum(self._sourcecounts) / np.sum(self._cmfd_src)

                # Define target reshape dimensions for sourcecounts. This defines
                # how self._sourcecounts is ordered by dimension
                target_shape = [nz, ny, nx, ng]

                # Reshape sourcecounts to target shape. Swap x and z axes so that
                # shape is now [nx, ny, nz, ng]
                sourcecounts = np.swapaxes(
                        self._sourcecounts.reshape(target_shape), 0, 2)

                # Flip index of energy dimension
                sourcecounts = np.flip(sourcecounts, axis=3)

                # Compute weight factors
                divide_condition = np.logical_and(sourcecounts > 0,
                                                  self._cmfd_src > 0)
                self._weightfactors = np.divide(self._cmfd_src * norm, \
                        sourcecounts, where=divide_condition, \
                        out=np.ones_like(self._cmfd_src))

            if not self._cmfd_feedback:
                return

            # Broadcast weight factors to all procs
            if have_mpi:
                self._weightfactors = self._intracomm.bcast(
                                      self._weightfactors, root=0)

            m = openmc.capi.meshes[self._cmfd_mesh_id]
            energy = self._egrid
            ng = self._indices[3]

            # Get xyz locations and energies of all particles in source bank
            source_xyz = openmc.capi.source_bank()['xyz']
            source_energies = openmc.capi.source_bank()['E']

            # Convert xyz location to the CMFD mesh index
            mesh_ijk = np.floor((source_xyz - m.lower_left) / m.width).astype(int)

            # TEMP
            if self._updated:
                # Determine which energy bin each particle's energy belongs to
                # Separate into cases bases on where source energies lies on egrid
                energy_bins = np.zeros(len(source_energies), dtype=int)
                idx = np.where(source_energies < energy[0])
                energy_bins[idx] = ng - 1
                idx = np.where(source_energies > energy[-1])
                energy_bins[idx] = 0
                idx = np.where((source_energies >= energy[0]) & (source_energies <= energy[-1]))
                energy_bins[idx] = ng - np.digitize(source_energies, energy)
            else:
                # Determine which energy bin each particle's energy belongs to
                energy_bins = np.where(source_energies < energy[0], ng - 1,
                        np.where(source_energies > energy[-1], 0, \
                        ng - np.digitize(source_energies, energy)))

            # Determine weight factor of each particle based on its mesh index
            # and energy bin and updates its weight
            openmc.capi.source_bank()['wgt'] *= self._weightfactors[ \
                    mesh_ijk[:,0], mesh_ijk[:,1], mesh_ijk[:,2], energy_bins]

            if openmc.capi.master() and np.any(source_energies < energy[0]):
                print(' WARNING: Source pt below energy grid')
                sys.stdout.flush()
            if openmc.capi.master() and np.any(source_energies > energy[-1]):
                print(' WARNING: Source pt above energy grid')
                sys.stdout.flush()

    def _count_bank_sites(self):
        """Determines the number of fission bank sites in each cell of a given
        mesh and energy group structure.
        Returns
        -------
        bool
            Wheter any source sites outside of CMFD mesh were found

        """
        # Initialize variables
        m = openmc.capi.meshes[self._cmfd_mesh_id]
        bank = openmc.capi.source_bank()
        energy = self._egrid
        sites_outside = np.zeros(1, dtype=bool)
        ng = self._indices[3]

        outside = np.zeros(1, dtype=bool)
        count = np.zeros(self._sourcecounts.shape)

        # Get xyz locations and energies of each particle in source bank
        source_xyz = openmc.capi.source_bank()['xyz']
        source_energies = openmc.capi.source_bank()['E']

        # Convert xyz location to mesh index and ravel index to scalar
        mesh_locations = np.floor((source_xyz - m.lower_left) / m.width)
        mesh_bins = mesh_locations[:,2] * m.dimension[1] * m.dimension[0] + \
                    mesh_locations[:,1] * m.dimension[0] + mesh_locations[:,0]

        # Check if any source locations lie outside of defined CMFD mesh
        if np.any(mesh_bins < 0) or np.any(mesh_bins >= np.prod(m.dimension)):
            outside[0] = True

        # TEMP
        if self._updated:
            # Determine which energy bin each particle's energy belongs to
            # Separate into cases bases on where source energies lies on egrid
            energy_bins = np.zeros(len(source_energies), dtype=int)
            idx = np.where(source_energies < energy[0])
            energy_bins[idx] = 0
            idx = np.where(source_energies > energy[-1])
            energy_bins[idx] = ng - 1
            idx = np.where((source_energies >= energy[0]) & (source_energies <= energy[-1]))
            energy_bins[idx] = np.digitize(source_energies, energy) - 1
        else:
            # Determine which energy bin each particle's energy corresponds to
            energy_bins = np.where(source_energies < energy[0], 0,
                    np.where(source_energies > energy[-1], ng - 1, \
                             np.digitize(source_energies, energy) - 1))

        # Determine all unique combinations of mesh bin and energy bin, and
        # count number of particles that belong to these combinations
        idx, counts = np.unique(np.array([mesh_bins, energy_bins]), axis=1,
                                return_counts=True)

        # Store counts to appropriate mesh-energy combination
        count[idx[0].astype(int), idx[1].astype(int)] = counts

        if have_mpi:
            # Collect values of count from all processors
            self._intracomm.Reduce(count, self._sourcecounts, MPI.SUM, root=0)
            # Check if there were sites outside the mesh for any processor
            self._intracomm.Reduce(outside, sites_outside, MPI.LOR, root=0)
        # Deal with case if MPI not defined (only one proc)
        else:
            sites_outside = outside
            self._sourcecounts = count

        return sites_outside[0]

    def _build_matrices(self, adjoint, vectorized):
        """Build loss and production matrices and write these matrices

        Parameters
        ----------
        adjoint : bool
            Whether or not to run an adjoint calculation
        vectorized : bool
            Whether to build CMFD matrices using a vectorized numpy approach
            or with traditional for loops

        Returns
        -------
        loss : scipy.sparse.spmatrix
            Sparse matrix storing elements of CMFD loss matrix
        prod : scipy.sparse.spmatrix
            Sparse matrix storing elements of CMFD production matrix

        """
        # Build loss and production matrices
        if vectorized:
            if self._updated:
                loss = self._build_loss_matrix_vectorized_updated(adjoint)
            else:
                loss = self._build_loss_matrix_vectorized(adjoint)
            if self._updated:
                prod = self._build_prod_matrix_vectorized_updated(adjoint)
            else:
                prod = self._build_prod_matrix_vectorized(adjoint)
        else:
            loss = self._build_loss_matrix(adjoint)
            prod = self._build_prod_matrix(adjoint)

        # Write out matrices
        if self._cmfd_write_matrices:
            if adjoint:
                self._write_matrix(loss, 'adj_loss')
                self._write_matrix(prod, 'adj_prod')
            else:
                self._write_matrix(loss, 'loss')
                self._write_matrix(prod, 'prod')

        return loss, prod

    def _compute_adjoint(self, loss, prod):
        """Computes a mathematical adjoint of the CMFD problem by transposing
        production and loss matrices passed in as arguments

        Parameters
        ----------
        loss : scipy.sparse.spmatrix
            Sparse matrix storing elements of CMFD loss matrix
        prod : scipy.sparse.spmatrix
            Sparse matrix storing elements of CMFD production matrix

        Returns
        -------
        loss : scipy.sparse.spmatrix
            Sparse matrix storing elements of adjoint CMFD loss matrix
        prod : scipy.sparse.spmatrix
            Sparse matrix storing elements of adjoint CMFD production matrix

        """
        # Transpose matrices
        loss = np.transpose(loss)
        prod = np.transpose(prod)

        # Write out matrices
        if self._cmfd_write_matrices:
            self._write_matrix(loss, 'adj_loss')
            self._write_matrix(prod, 'adj_prod')

        return loss, prod

    def _build_loss_matrix_vectorized_updated(self, adjoint):
        # Extract spatial and energy indices and define matrix dimension
        ng = self._indices[3]
        n = self._mat_dim*ng

        # Define data entries used to build csr matrix
        data = np.array([])

        dtilde_left = self._dtilde[:,:,:,:,0]
        dtilde_right = self._dtilde[:,:,:,:,1]
        dtilde_back = self._dtilde[:,:,:,:,2]
        dtilde_front = self._dtilde[:,:,:,:,3]
        dtilde_bottom = self._dtilde[:,:,:,:,4]
        dtilde_top = self._dtilde[:,:,:,:,5]
        dhat_left = self._dhat[:,:,:,:,0]
        dhat_right = self._dhat[:,:,:,:,1]
        dhat_back = self._dhat[:,:,:,:,2]
        dhat_front = self._dhat[:,:,:,:,3]
        dhat_bottom = self._dhat[:,:,:,:,4]
        dhat_top = self._dhat[:,:,:,:,5]

        dx = self._hxyz[:,:,:,np.newaxis,0]
        dy = self._hxyz[:,:,:,np.newaxis,1]
        dz = self._hxyz[:,:,:,np.newaxis,2]

        # Define net leakage coefficient for each surface in each matrix element
        jnet = (((dtilde_right + dhat_right) - (-1.0 * dtilde_left + dhat_left)) / dx
             +  ((dtilde_front + dhat_front) - (-1.0 * dtilde_back + dhat_back)) / dy
             +  ((dtilde_top + dhat_top) - (-1.0 * dtilde_bottom + dhat_bottom)) / dz)

        for g in range(ng):
            # Define leakage terms that relate terms to their neighbors to the
            # left
            dtilde = self._dtilde[:,:,:,g,0][self._accel_neig_left_idxs]
            dhat = self._dhat[:,:,:,g,0][self._accel_neig_left_idxs]
            dx = self._hxyz[:,:,:,0][self._accel_neig_left_idxs]
            vals = (-1.0 * dtilde - dhat) / dx
            # Store data to add to CSR matrix
            data = np.append(data, vals)

            # Define leakage terms that relate terms to their neighbors to the
            # right
            dtilde = self._dtilde[:,:,:,g,1][self._accel_neig_right_idxs]
            dhat = self._dhat[:,:,:,g,1][self._accel_neig_right_idxs]
            dx = self._hxyz[:,:,:,0][self._accel_neig_right_idxs]
            vals = (-1.0 * dtilde + dhat) / dx
            # Store data to add to CSR matrix
            data = np.append(data, vals)

            # Define leakage terms that relate terms to their neighbors in the
            # back
            dtilde = self._dtilde[:,:,:,g,2][self._accel_neig_back_idxs]
            dhat = self._dhat[:,:,:,g,2][self._accel_neig_back_idxs]
            dy = self._hxyz[:,:,:,1][self._accel_neig_back_idxs]
            vals = (-1.0 * dtilde - dhat) / dy
            # Store data to add to CSR matrix
            data = np.append(data, vals)

            # Define leakage terms that relate terms to their neighbors in the
            # front
            dtilde = self._dtilde[:,:,:,g,3][self._accel_neig_front_idxs]
            dhat = self._dhat[:,:,:,g,3][self._accel_neig_front_idxs]
            dy = self._hxyz[:,:,:,1][self._accel_neig_front_idxs]
            vals = (-1.0 * dtilde + dhat) / dy
            # Store data to add to CSR matrix
            data = np.append(data, vals)

            # Define leakage terms that relate terms to their neighbors to the
            # bottom
            dtilde = self._dtilde[:,:,:,g,4][self._accel_neig_bottom_idxs]
            dhat = self._dhat[:,:,:,g,4][self._accel_neig_bottom_idxs]
            dz = self._hxyz[:,:,:,2][self._accel_neig_bottom_idxs]
            vals = (-1.0 * dtilde - dhat) / dz
            # Store data to add to CSR matrix
            data = np.append(data, vals)

            # Define leakage terms that relate terms to their neighbors to the
            # top
            dtilde = self._dtilde[:,:,:,g,5][self._accel_neig_top_idxs]
            dhat = self._dhat[:,:,:,g,5][self._accel_neig_top_idxs]
            dz = self._hxyz[:,:,:,2][self._accel_neig_top_idxs]
            vals = (-1.0 * dtilde + dhat) / dz
            # Store data to add to CSR matrix
            data = np.append(data, vals)

            # Define terms that relate to loss of neutrons in a cell. These
            # correspond to all the diagonal entries of the loss matrix
            jnet_g = jnet[:,:,:,g][self._accel_idxs]
            total_xs = self._totalxs[:,:,:,g][self._accel_idxs]
            scatt_xs = self._scattxs[:,:,:,g,g][self._accel_idxs]
            vals = jnet_g + total_xs - scatt_xs
            # Store data to add to CSR matrix
            data = np.append(data, vals)

            # Define terms that relate to in-scattering from group to group.
            # These terms relate a mesh index to all mesh indices with the same
            # spatial dimensions but belong to a different energy group
            for h in range(ng):
                if h != g:
                    # Get scattering macro xs, transposed
                    if adjoint:
                        scatt_xs = self._scattxs[:,:,:,g,h][self._accel_idxs]
                    # Get scattering macro xs
                    else:
                        scatt_xs = self._scattxs[:,:,:,h,g][self._accel_idxs]
                    vals = -1.0 * scatt_xs
                    # Store data to add to CSR matrix
                    data = np.append(data, vals)

        # Create csr matrix
        loss = sparse.csr_matrix((data, (self._loss_row, self._loss_col)), shape=(n, n))
        return loss

    def _build_prod_matrix_vectorized_updated(self, adjoint):
        # Extract spatial and energy indices and define matrix dimension
        ng = self._indices[3]
        n = self._mat_dim*ng

        # Define rows, columns, and data used to build csr matrix
        data = np.array([])

        # Define terms that relate to fission production from group to group.
        for g in range(ng):
            for h in range(ng):
                # Get nu-fission macro xs, transposed
                if adjoint:
                    vals = (self._nfissxs[:, :, :, g, h])[self._accel_idxs]
                # Get nu-fission macro xs
                else:
                    vals = (self._nfissxs[:, :, :, h, g])[self._accel_idxs]
                # Store rows, cols, and data to add to CSR matrix
                data = np.append(data, vals)

        # Create csr matrix
        prod = sparse.csr_matrix((data, (self._prod_row, self._prod_col)), shape=(n, n))
        return prod

    def _execute_power_iter(self, loss, prod):
        """Main power iteration routine for the CMFD calculation

        Parameters
        ----------
        loss : scipy.sparse.spmatrix
            Sparse matrix storing elements of CMFD loss matrix
        prod : scipy.sparse.spmatrix
            Sparse matrix storing elements of CMFD production matrix

        Returns
        -------
        phi_n : numpy.ndarray
            Flux vector of CMFD problem
        k_n : float
            Eigenvalue of CMFD problem
        dom : float
            Dominance ratio of CMFD problem

        """
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
        k_n = openmc.capi.keff_temp()[0]
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
            iconv, norm_n = self._check_convergence(s_n, s_o, k_n, k_o, i+1)

            # If converged, calculate dominance ratio and break from loop
            if iconv:
                dom = norm_n / norm_o
                return phi_n, k_n, dom

            # Record old values if not converged
            phi_o = phi_n
            k_o = k_n
            k_lo = k_ln
            norm_o = norm_n

    def _check_convergence(self, s_n, s_o, k_n, k_o, iter):
        """Checks the convergence of the CMFD problem

        Parameters
        ----------
        s_n : numpy.ndarray
            Source vector from current iteration
        s_o : numpy.ndarray
            Source vector from previous iteration
        k_n : float
            K-effective  from current iteration
        k_o : float
            K-effective from previous iteration
        iter: int
            Iteration number

        Returns
        -------
        iconv : bool
            Whether the power iteration has reached convergence
        serr : float
            Error in source from previous iteration to current iteration, used
            for dominance ratio calculations

        """
        # Calculate error in keff
        kerr = abs(k_o - k_n)/k_n

        # Calculate max error in source
        serr = np.sqrt(np.sum(np.where(s_n>0, ((s_n-s_o)/s_n)**2,0))/len(s_n))

        # Check for convergence
        iconv = kerr < self._cmfd_ktol and serr < self._cmfd_stol

        # Print out to user
        if self._cmfd_power_monitor and openmc.capi.master():
            str1 = ' {:d}:'.format(iter)
            str2 = 'k-eff: {:0.8f}'.format(k_n)
            str3 = 'k-error:  {0:.5E}'.format(kerr)
            str4 = 'src-error:  {0:.5E}'.format(serr)
            print('{0:8s}{1:20s}{2:25s}{3:s}'.format(str1, str2, str3, str4))
            sys.stdout.flush()

        return iconv, serr

    def _set_coremap(self):
        """Sets the core mapping information. All regions marked with zero
        are set to CMFD_NO_ACCEL, while all regions marked with 1 are set to a
        unique index that maps each fuel region to a row number when building
        CMFD matrices

        """
        # Set number of accelerated regions in problem. This will be related to
        # the dimension of CMFD matrices
        self._mat_dim = np.sum(self._coremap)

        # Define coremap as cumulative sum over accelerated regions,
        # otherwise set value to _CMFD_NOACCEL
        self._coremap = np.where(self._coremap==0, _CMFD_NOACCEL,
                        np.cumsum(self._coremap)-1)

        # Reshape coremap to three dimensional array
        # Indices of coremap in user input switched in x and z axes
        nx = self._indices[0]
        ny = self._indices[1]
        nz = self._indices[2]
        self._coremap = self._coremap.reshape(nz, ny, nx)
        self._coremap = np.swapaxes(self._coremap, 0, 2)

    def _compute_xs(self):
        """Takes CMFD tallies from OpenMC and computes macroscopic cross
        sections, flux, and diffusion coefficients for each mesh cell

        """
        # Extract spatial and energy indices
        nx = self._indices[0]
        ny = self._indices[1]
        nz = self._indices[2]
        ng = self._indices[3]

        # Set flux object and source distribution all to zeros
        self._flux.fill(0.)
        self._openmc_src.fill(0.)

        # Set mesh widths
        self._hxyz[:,:,:,:] = openmc.capi.meshes[self._cmfd_mesh_id].width

        # Reset keff_bal to zero
        self._keff_bal = 0.

        # Get tallies in-memory
        tallies = openmc.capi.tallies

        # Ravel coremap as 1d array similar to how tally data is arranged
        coremap = np.ravel(self._coremap.swapaxes(0, 2))

        # Set conditional numpy array as boolean vector based on coremap
        # Repeat each value for number of groups in problem
        is_cmfd_accel = np.repeat(coremap != _CMFD_NOACCEL, ng)

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

        # Define target tally reshape dimensions. This defines how openmc
        # tallies are ordered by dimension
        target_tally_shape = [nz, ny, nx, ng]

        # Reshape flux array to target shape. Swap x and z axes so that
        # flux shape is now [nx, ny, nz, ng]
        reshape_flux = np.swapaxes(flux.reshape(target_tally_shape), 0, 2)

        # Flip energy axis as tally results are given in reverse order of
        # energy group
        self._flux = np.flip(reshape_flux, axis=3)

        # Get total rr and convert to total xs from CMFD tally 0
        tally_results = tallies[tally_id].results[:,1,1]
        totalxs = np.divide(tally_results, flux, \
                            where=flux>0, out=np.zeros_like(tally_results))

        # Reshape totalxs array to target shape. Swap x and z axes so that
        # shape is now [nx, ny, nz, ng]
        reshape_totalxs = np.swapaxes(totalxs.reshape(target_tally_shape), 0, 2)

        # Total xs is flipped in energy axis as tally results are given in
        # reverse order of energy group
        self._totalxs = np.flip(reshape_totalxs, axis=3)

        # Get scattering xs from CMFD tally 1
        # flux is repeated to account for extra dimensionality of scattering xs
        tally_id = self._cmfd_tally_ids[1]
        tally_results = tallies[tally_id].results[:,0,1]
        scattxs = np.divide(tally_results, \
                            np.repeat(flux, ng), \
                            where=np.repeat(flux>0, ng), \
                            out=np.zeros_like(tally_results))

        # Define target tally reshape dimensions for xs with incoming
        # and outgoing energies
        target_tally_shape = [nz, ny, nx, ng, ng]

        # Reshape scattxs array to target shape. Swap x and z axes so that
        # shape is now [nx, ny, nz, ng, ng]
        reshape_scattxs = np.swapaxes(scattxs.reshape(target_tally_shape), 0, 2)

        # Scattering xs is flipped in both incoming and outgoing energy axes
        # as tally results are given in reverse order of energy group
        self._scattxs = np.flip(reshape_scattxs, axis=3)
        self._scattxs = np.flip(self._scattxs, axis=4)

        # Get nu-fission xs from CMFD tally 1
        # flux is repeated to account for extra dimensionality of nu-fission xs
        tally_results = tallies[tally_id].results[:,1,1]
        num_realizations = tallies[tally_id].num_realizations
        nfissxs = np.divide(tally_results, \
                            np.repeat(flux, ng), \
                            where=np.repeat(flux>0, ng), \
                            out=np.zeros_like(tally_results))

        # Reshape nfissxs array to target shape. Swap x and z axes so that
        # shape is now [nx, ny, nz, ng, ng]
        reshape_nfissxs = np.swapaxes(nfissxs.reshape(target_tally_shape), 0, 2)

        # Nu-fission xs is flipped in both incoming and outgoing energy axes
        # as tally results are given in reverse order of energy group
        self._nfissxs = np.flip(reshape_nfissxs, axis=3)
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
        current = np.where(np.repeat(flux>0, 12), tally_results, 0.)

        # Define target tally reshape dimensions for current
        target_tally_shape = [nz, ny, nx, 12, ng]

        # Reshape current array to target shape. Swap x and z axes so that
        # shape is now [nx, ny, nz, ng, 12]
        reshape_current = np.swapaxes(current.reshape(target_tally_shape), 0, 2)

        # Current is flipped in energy axis as tally results are given in
        # reverse order of energy group
        self._current = np.flip(reshape_current, axis=4)

        # Get p1 scatter xs from CMFD tally 3
        tally_id = self._cmfd_tally_ids[3]
        tally_results = tallies[tally_id].results[:,0,1]

        # Define target tally reshape dimensions for p1 scatter tally
        target_tally_shape = [nz, ny, nx, 2, ng]

        # Reshape and extract only p1 data from tally results (no need for p0 data)
        p1scattrr = np.swapaxes(tally_results.reshape(target_tally_shape),
                                0, 2)[:,:,:,1,:]

        # Store p1 scatter xs
        # p1 scatter xs is flipped in energy axis as tally results are given in
        # reverse order of energy group
        self._p1scattxs = np.divide(np.flip(p1scattrr, axis=3), self._flux, \
                                    where=self._flux>0, \
                                    out=np.zeros_like(p1scattrr))

        # Calculate and store diffusion coefficient
        self._diffcof = np.where(self._flux>0, 1.0 / (3.0 * \
                                 (self._totalxs - self._p1scattxs)), 0.)

    def _compute_effective_downscatter(self):
        """Changes downscatter rate for zero upscatter"""
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

        # First energy index is incoming energy, second is outgoing energy
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
        """Computes the RMS neutron balance over the CMFD mesh"""
        # Extract energy indices
        ng = self._indices[3]

        # Get number of accelerated regions
        num_accel = self._mat_dim

        # Get openmc k-effective
        keff = openmc.capi.keff_temp()[0]

        # Define leakage in each mesh cell and energy group
        leakage = ((self._current[:,:,:,_CURRENTS['out_right'],:] - \
            self._current[:,:,:,_CURRENTS['in_right'],:]) - \
            (self._current[:,:,:,_CURRENTS['in_left'],:] - \
            self._current[:,:,:,_CURRENTS['out_left'],:])) + \
            ((self._current[:,:,:,_CURRENTS['out_front'],:] - \
            self._current[:,:,:,_CURRENTS['in_front'],:]) - \
            (self._current[:,:,:,_CURRENTS['in_back'],:] - \
            self._current[:,:,:,_CURRENTS['out_back'],:])) + \
            ((self._current[:,:,:,_CURRENTS['out_top'],:] - \
            self._current[:,:,:,_CURRENTS['in_top'],:]) - \
            (self._current[:,:,:,_CURRENTS['in_bottom'],:] - \
            self._current[:,:,:,_CURRENTS['out_bottom'],:]))

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
            (ng * num_accel)))

    def _precompute_array_indices(self):
        """Computes the indices used to populate certain cross section arrays.
        These indices are used in _compute_dtilde and _compute_dhat

        """
        # Extract spatial indices
        nx = self._indices[0]
        ny = self._indices[1]
        nz = self._indices[2]


        # Logical for determining whether region of interest is accelerated region
        is_accel = self._coremap != _CMFD_NOACCEL
        # Logical for determining whether a zero flux "albedo" b.c. should be
        # applied
        is_zero_flux_alb = abs(self._albedo - _ZERO_FLUX) < _TINY_BIT
        x_inds, y_inds, z_inds = np.indices((nx, ny, nz))

        # Define slice equivalent to _accel[0,:,:]
        slice_x = x_inds[:1,:,:]
        slice_y = y_inds[:1,:,:]
        slice_z = z_inds[:1,:,:]
        bndry_accel = is_accel[(slice_x, slice_y, slice_z)]
        self._first_x_accel = (slice_x[bndry_accel], slice_y[bndry_accel], slice_z[bndry_accel])

        # Define slice equivalent to is_accel[-1,:,:]
        slice_x = x_inds[-1:,:,:]
        slice_y = y_inds[-1:,:,:]
        slice_z = z_inds[-1:,:,:]
        bndry_accel = is_accel[(slice_x, slice_y, slice_z)]
        self._last_x_accel = (slice_x[bndry_accel], slice_y[bndry_accel], slice_z[bndry_accel])

        # Define slice equivalent to is_accel[:,0,:]
        slice_x = x_inds[:,:1,:]
        slice_y = y_inds[:,:1,:]
        slice_z = z_inds[:,:1,:]
        bndry_accel = is_accel[(slice_x, slice_y, slice_z)]
        self._first_y_accel = (slice_x[bndry_accel], slice_y[bndry_accel], slice_z[bndry_accel])

        # Define slice equivalent to is_accel[:,-1,:]
        slice_x = x_inds[:,-1:,:]
        slice_y = y_inds[:,-1:,:]
        slice_z = z_inds[:,-1:,:]
        bndry_accel = is_accel[(slice_x, slice_y, slice_z)]
        self._last_y_accel = (slice_x[bndry_accel], slice_y[bndry_accel], slice_z[bndry_accel])

        # Define slice equivalent to is_accel[:,:,0]
        slice_x = x_inds[:,:,:1]
        slice_y = y_inds[:,:,:1]
        slice_z = z_inds[:,:,:1]
        bndry_accel = is_accel[(slice_x, slice_y, slice_z)]
        self._first_z_accel = (slice_x[bndry_accel], slice_y[bndry_accel], slice_z[bndry_accel])

        # Define slice equivalent to is_accel[:,:,-1]
        slice_x = x_inds[:,:,-1:]
        slice_y = y_inds[:,:,-1:]
        slice_z = z_inds[:,:,-1:]
        bndry_accel = is_accel[(slice_x, slice_y, slice_z)]
        self._last_z_accel = (slice_x[bndry_accel], slice_y[bndry_accel], slice_z[bndry_accel])

        # Define slice equivalent to is_accel[1:,:,:]
        slice_x = x_inds[1:,:,:]
        slice_y = y_inds[1:,:,:]
        slice_z = z_inds[1:,:,:]
        bndry_accel = is_accel[(slice_x, slice_y, slice_z)]
        self._notfirst_x_accel = (slice_x[bndry_accel], slice_y[bndry_accel],
            slice_z[bndry_accel])

        # Define slice equivalent to is_accel[:-1,:,:]
        slice_x = x_inds[:-1,:,:]
        slice_y = y_inds[:-1,:,:]
        slice_z = z_inds[:-1,:,:]
        bndry_accel = is_accel[(slice_x, slice_y, slice_z)]
        self._notlast_x_accel = (slice_x[bndry_accel], slice_y[bndry_accel],
             slice_z[bndry_accel])

        # Define slice equivalent to is_accel[:,1:,:]
        slice_x = x_inds[:,1:,:]
        slice_y = y_inds[:,1:,:]
        slice_z = z_inds[:,1:,:]
        bndry_accel = is_accel[(slice_x, slice_y, slice_z)]
        self._notfirst_y_accel = (slice_x[bndry_accel], slice_y[bndry_accel],
             slice_z[bndry_accel])

        # Define slice equivalent to is_accel[:,:-1,:]
        slice_x = x_inds[:,:-1,:]
        slice_y = y_inds[:,:-1,:]
        slice_z = z_inds[:,:-1,:]
        bndry_accel = is_accel[(slice_x, slice_y, slice_z)]
        self._notlast_y_accel = (slice_x[bndry_accel], slice_y[bndry_accel],
             slice_z[bndry_accel])

        # Define slice equivalent to is_accel[:,:,1:]
        slice_x = x_inds[:,:,1:]
        slice_y = y_inds[:,:,1:]
        slice_z = z_inds[:,:,1:]
        bndry_accel = is_accel[(slice_x, slice_y, slice_z)]
        self._notfirst_z_accel = (slice_x[bndry_accel], slice_y[bndry_accel],
             slice_z[bndry_accel])

        # Define slice equivalent to is_accel[:,:,:-1]
        slice_x = x_inds[:,:,:-1]
        slice_y = y_inds[:,:,:-1]
        slice_z = z_inds[:,:,:-1]
        bndry_accel = is_accel[(slice_x, slice_y, slice_z)]
        self._notlast_z_accel = (slice_x[bndry_accel], slice_y[bndry_accel],
             slice_z[bndry_accel])

        # Store logical for whether neighboring cell is reflector region
        # in all directions
        self._adj_reflector_left = np.roll(self._coremap, 1, axis=0) == _CMFD_NOACCEL
        self._adj_reflector_right = np.roll(self._coremap, -1, axis=0) == _CMFD_NOACCEL
        self._adj_reflector_back = np.roll(self._coremap, 1, axis=1) == _CMFD_NOACCEL
        self._adj_reflector_front = np.roll(self._coremap, -1, axis=1) == _CMFD_NOACCEL
        self._adj_reflector_bottom = np.roll(self._coremap, 1, axis=2) == _CMFD_NOACCEL
        self._adj_reflector_top = np.roll(self._coremap, -1, axis=2) == _CMFD_NOACCEL

    def _precompute_matrix_indices(self):
        """Computes the indices and row/column data used to populate CMFD CSR matrices.
        These indices are used in _build_loss_matrix and _build_prod_matrix

        """
        # Extract energy group indices
        ng = self._indices[3]

        # Shift coremap in all directions to determine whether leakage term
        # should be defined for particular cell in matrix
        coremap_shift_left = np.pad(self._coremap, ((1,0),(0,0),(0,0)),
                mode='constant', constant_values=_CMFD_NOACCEL)[:-1,:,:]

        coremap_shift_right = np.pad(self._coremap, ((0,1),(0,0),(0,0)),
                mode='constant', constant_values=_CMFD_NOACCEL)[1:,:,:]

        coremap_shift_back = np.pad(self._coremap, ((0,0),(1,0),(0,0)),
                mode='constant', constant_values=_CMFD_NOACCEL)[:,:-1,:]

        coremap_shift_front = np.pad(self._coremap, ((0,0),(0,1),(0,0)),
                mode='constant', constant_values=_CMFD_NOACCEL)[:,1:,:]

        coremap_shift_bottom = np.pad(self._coremap, ((0,0),(0,0),(1,0)),
                mode='constant', constant_values=_CMFD_NOACCEL)[:,:,:-1]

        coremap_shift_top = np.pad(self._coremap, ((0,0),(0,0),(0,1)),
                mode='constant', constant_values=_CMFD_NOACCEL)[:,:,1:]

        # Create empty row and column vectors to store for loss matrix
        row = np.array([], dtype=int)
        col = np.array([], dtype=int)

        # Store all indices used to populate production and loss matrix
        self._accel_idxs = np.where(self._coremap != _CMFD_NOACCEL)
        self._accel_neig_left_idxs = np.where((self._coremap != _CMFD_NOACCEL) &
             (coremap_shift_left != _CMFD_NOACCEL))
        self._accel_neig_right_idxs = np.where((self._coremap != _CMFD_NOACCEL) &
             (coremap_shift_right != _CMFD_NOACCEL))
        self._accel_neig_back_idxs = np.where((self._coremap != _CMFD_NOACCEL) &
             (coremap_shift_back != _CMFD_NOACCEL))
        self._accel_neig_front_idxs = np.where((self._coremap != _CMFD_NOACCEL) &
             (coremap_shift_front != _CMFD_NOACCEL))
        self._accel_neig_bottom_idxs = np.where((self._coremap != _CMFD_NOACCEL) &
             (coremap_shift_bottom != _CMFD_NOACCEL))
        self._accel_neig_top_idxs = np.where((self._coremap != _CMFD_NOACCEL) &
             (coremap_shift_top != _CMFD_NOACCEL))

        for g in range(ng):
            # Extract row and column data of regions where a cell and its
            # neighbor to the left are both fuel regions
            idx_x = ng * (self._coremap[self._accel_neig_left_idxs]) + g
            idx_y = ng * (coremap_shift_left[self._accel_neig_left_idxs]) + g
            row = np.append(row, idx_x)
            col = np.append(col, idx_y)

            # Extract row and column data of regions where a cell and its
            # neighbor to the right are both fuel regions
            idx_x = ng * (self._coremap[self._accel_neig_right_idxs]) + g
            idx_y = ng * (coremap_shift_right[self._accel_neig_right_idxs]) + g
            row = np.append(row, idx_x)
            col = np.append(col, idx_y)

            # Extract row and column data of regions where a cell and its
            # neighbor to the back are both fuel regions
            idx_x = ng * (self._coremap[self._accel_neig_back_idxs]) + g
            idx_y = ng * (coremap_shift_back[self._accel_neig_back_idxs]) + g
            row = np.append(row, idx_x)
            col = np.append(col, idx_y)

            # Extract row and column data of regions where a cell and its
            # neighbor to the front are both fuel regions
            idx_x = ng * (self._coremap[self._accel_neig_front_idxs]) + g
            idx_y = ng * (coremap_shift_front[self._accel_neig_front_idxs]) + g
            row = np.append(row, idx_x)
            col = np.append(col, idx_y)

            # Extract row and column data of regions where a cell and its
            # neighbor to the bottom are both fuel regions
            idx_x = ng * (self._coremap[self._accel_neig_bottom_idxs]) + g
            idx_y = ng * (coremap_shift_bottom[self._accel_neig_bottom_idxs]) + g
            row = np.append(row, idx_x)
            col = np.append(col, idx_y)

            # Extract row and column data of regions where a cell and its
            # neighbor to the top are both fuel regions
            idx_x = ng * (self._coremap[self._accel_neig_top_idxs]) + g
            idx_y = ng * (coremap_shift_top[self._accel_neig_top_idxs]) + g
            row = np.append(row, idx_x)
            col = np.append(col, idx_y)

            # Extract all regions where a cell is a fuel region
            idx_x = ng * (self._coremap[self._accel_idxs]) + g
            idx_y = idx_x
            row = np.append(row, idx_x)
            col = np.append(col, idx_y)

            for h in range(ng):
                if h != g:
                    # Extract all regions where a cell is a fuel region
                    idx_x = ng * (self._coremap[self._accel_idxs]) + g
                    idx_y = ng * (self._coremap[self._accel_idxs]) + h
                    row = np.append(row, idx_x)
                    col = np.append(col, idx_y)

        # Store row and col as rows and columns of production matrix
        self._loss_row = row
        self._loss_col = col

        # Create empty row and column vectors to store for production matrix
        row = np.array([], dtype=int)
        col = np.array([], dtype=int)

        for g in range(ng):
            for h in range(ng):
                # Extract all regions where a cell is a fuel region
                idx_x = ng * (self._coremap[self._accel_idxs]) + g
                idx_y = ng * (self._coremap[self._accel_idxs]) + h
                # Store rows, cols, and data to add to CSR matrix
                row = np.append(row, idx_x)
                col = np.append(col, idx_y)

        # Store row and col as rows and columns of production matrix
        self._prod_row = row
        self._prod_col = col

    def _compute_dtilde_vectorized_updated(self):
        # TODO: Update comment

        # Logical for determining whether a zero flux "albedo" b.c. should be
        # applied
        is_zero_flux_alb = abs(self._albedo - _ZERO_FLUX) < _TINY_BIT

        # Define dtilde at left surface for all mesh cells on left boundary
        # Separate between zero flux b.c. and alebdo b.c.
        boundary = self._first_x_accel
        boundary_grps = boundary + (slice(None),)
        D = self._diffcof[boundary_grps]
        dx = self._hxyz[boundary + (np.newaxis, 0)]
        if is_zero_flux_alb[0]:
            self._dtilde[boundary_grps + (0,)] = 2.0 * D / dx
        else:
            alb = self._albedo[0]
            self._dtilde[boundary_grps + (0,)] = ((2.0 * D * (1.0 - alb))
                 / (4.0 * D * (1.0 + alb) + (1.0 - alb) * dx))

        # Define dtilde at right surface for all mesh cells on right boundary
        # Separate between zero flux b.c. and alebdo b.c.
        boundary = self._last_x_accel
        boundary_grps = boundary + (slice(None),)
        D = self._diffcof[boundary_grps]
        dx = self._hxyz[boundary + (np.newaxis, 0)]
        if is_zero_flux_alb[1]:
            self._dtilde[boundary_grps + (1,)] = 2.0 * D / dx
        else:
            alb = self._albedo[1]
            self._dtilde[boundary_grps + (1,)] = ((2.0 * D * (1.0 - alb))
                 / (4.0 * D * (1.0 + alb) + (1.0 - alb) * dx))

        # Define dtilde at back surface for all mesh cells on back boundary
        # Separate between zero flux b.c. and alebdo b.c.
        boundary = self._first_y_accel
        boundary_grps = boundary + (slice(None),)
        D = self._diffcof[boundary_grps]
        dy = self._hxyz[boundary + (np.newaxis, 1)]
        if is_zero_flux_alb[2]:
            self._dtilde[boundary_grps + (2,)] = 2.0 * D / dy
        else:
            alb = self._albedo[2]
            self._dtilde[boundary_grps + (2,)] = ((2.0 * D * (1.0 - alb))
                 / (4.0 * D * (1.0 + alb) + (1.0 - alb) * dy))

        # Define dtilde at front surface for all mesh cells on front boundary
        # Separate between zero flux b.c. and alebdo b.c.
        boundary = self._last_y_accel
        boundary_grps = boundary + (slice(None),)
        D = self._diffcof[boundary_grps]
        dy = self._hxyz[boundary + (np.newaxis, 1)]
        if is_zero_flux_alb[3]:
            self._dtilde[boundary_grps + (3,)] = 2.0 * D / dy
        else:
            alb = self._albedo[3]
            self._dtilde[boundary_grps + (3,)] = ((2.0 * D * (1.0 - alb))
                 / (4.0 * D * (1.0 + alb) + (1.0 - alb) * dy))

        # Define dtilde at bottom surface for all mesh cells on bottom boundary
        # Separate between zero flux b.c. and alebdo b.c.
        boundary = self._first_z_accel
        boundary_grps = boundary + (slice(None),)
        D = self._diffcof[boundary_grps]
        dz = self._hxyz[boundary + (np.newaxis, 2)]
        if is_zero_flux_alb[4]:
            self._dtilde[boundary_grps + (4,)] = 2.0 * D / dz
        else:
            alb = self._albedo[4]
            self._dtilde[boundary_grps + (4,)] = ((2.0 * D * (1.0 - alb))
                 / (4.0 * D * (1.0 + alb) + (1.0 - alb) * dz))

        # Define dtilde at top surface for all mesh cells on top boundary
        # Separate between zero flux b.c. and alebdo b.c.
        boundary = self._last_z_accel
        boundary_grps = boundary + (slice(None),)

        D = self._diffcof[boundary_grps]
        dz = self._hxyz[boundary + (np.newaxis, 2)]
        if is_zero_flux_alb[5]:
            self._dtilde[boundary_grps + (5,)] = 2.0 * D / dz
        else:
            alb = self._albedo[5]
            self._dtilde[boundary_grps + (5,)] = ((2.0 * D * (1 - alb))
                 / (4.0 * D * (1.0 + alb) + (1.0 - alb) * dz))

        # Define reflector albedo for all cells on the left surface, in case
        # a cell borders a reflector region on the left
        current_in_left = self._current[:,:,:,_CURRENTS['in_left'],:]
        current_out_left = self._current[:,:,:,_CURRENTS['out_left'],:]
        ref_albedo = np.divide(current_in_left, current_out_left,
                               where=current_out_left > 1.0e-10,
                               out=np.ones_like(current_out_left))

        # Diffusion coefficient of neighbor to left
        neig_dc = np.roll(self._diffcof, 1, axis=0)
        # Cell dimensions of neighbor to left
        neig_hxyz = np.roll(self._hxyz, 1, axis=0)

        # Define dtilde at left surface for all mesh cells not on left boundary
        # Dtilde is defined differently for regions that do and don't neighbor
        # reflector regions
        boundary = self._notfirst_x_accel
        boundary_grps = boundary + (slice(None),)
        D = self._diffcof[boundary_grps]
        dx = self._hxyz[boundary + (np.newaxis, 0)]
        neig_D = neig_dc[boundary_grps]
        neig_dx = neig_hxyz[boundary + (np.newaxis, 0)]
        alb = ref_albedo[boundary_grps]
        is_adj_ref_left = self._adj_reflector_left[boundary + (np.newaxis,)]
        self._dtilde[boundary_grps + (0,)] = np.where(is_adj_ref_left,
             (2.0 * D * (1.0 - alb)) / (4.0 * D * (1.0 + alb) + (1.0 - alb) * dx),
             (2.0 * D * neig_D) / (neig_dx * D + dx * neig_D))

        # Define reflector albedo for all cells on the right surface, in case
        # a cell borders a reflector region on the right
        current_in_right = self._current[:,:,:,_CURRENTS['in_right'],:]
        current_out_right = self._current[:,:,:,_CURRENTS['out_right'],:]
        ref_albedo = np.divide(current_in_right, current_out_right,
                               where=current_out_right > 1.0e-10,
                               out=np.ones_like(current_out_right))

        # Diffusion coefficient of neighbor to right
        neig_dc = np.roll(self._diffcof, -1, axis=0)
        # Cell dimensions of neighbor to right
        neig_hxyz = np.roll(self._hxyz, -1, axis=0)

        # Define dtilde at right surface for all mesh cells not on right boundary
        # Dtilde is defined differently for regions that do and don't neighbor
        # reflector regions
        boundary = self._notlast_x_accel
        boundary_grps = boundary + (slice(None),)
        D = self._diffcof[boundary_grps]
        dx = self._hxyz[boundary + (np.newaxis, 0)]
        neig_D = neig_dc[boundary_grps]
        neig_dx = neig_hxyz[boundary + (np.newaxis, 0)]
        alb = ref_albedo[boundary_grps]
        is_adj_ref_right = self._adj_reflector_right[boundary + (np.newaxis,)]
        self._dtilde[boundary_grps + (1,)] = np.where(is_adj_ref_right,
             (2.0 * D * (1.0 - alb)) / (4.0 * D * (1.0 + alb) + (1.0 - alb) * dx),
             (2.0 * D * neig_D) / (neig_dx * D + dx * neig_D))

        # Define reflector albedo for all cells on the back surface, in case
        # a cell borders a reflector region on the back
        current_in_back = self._current[:,:,:,_CURRENTS['in_back'],:]
        current_out_back = self._current[:,:,:,_CURRENTS['out_back'],:]
        ref_albedo = np.divide(current_in_back, current_out_back,
                               where=current_out_back > 1.0e-10,
                               out=np.ones_like(current_out_back))

        # Diffusion coefficient of neighbor to back
        neig_dc = np.roll(self._diffcof, 1, axis=1)
        # Cell dimensions of neighbor to back
        neig_hxyz = np.roll(self._hxyz, 1, axis=1)

        # Define dtilde at back surface for all mesh cells not on back boundary
        # Dtilde is defined differently for regions that do and don't neighbor
        # reflector regions
        boundary = self._notfirst_y_accel
        boundary_grps = boundary + (slice(None),)
        D = self._diffcof[boundary_grps]
        dy = self._hxyz[boundary + (np.newaxis, 1)]
        neig_D = neig_dc[boundary_grps]
        neig_dy = neig_hxyz[boundary + (np.newaxis, 1)]
        alb = ref_albedo[boundary_grps]
        is_adj_ref_back = self._adj_reflector_back[boundary + (np.newaxis,)]
        self._dtilde[boundary_grps + (2,)] = np.where(is_adj_ref_back,
             (2.0 * D * (1.0 - alb)) / (4.0 * D * (1.0 + alb) + (1.0 - alb) * dy),
             (2.0 * D * neig_D) / (neig_dy * D + dy * neig_D))

        # Define reflector albedo for all cells on the front surface, in case
        # a cell borders a reflector region in the front
        current_in_front = self._current[:,:,:,_CURRENTS['in_front'],:]
        current_out_front = self._current[:,:,:,_CURRENTS['out_front'],:]
        ref_albedo = np.divide(current_in_front, current_out_front,
                               where=current_out_front > 1.0e-10,
                               out=np.ones_like(current_out_front))

        # Diffusion coefficient of neighbor to front
        neig_dc = np.roll(self._diffcof, -1, axis=1)
        # Cell dimensions of neighbor to front
        neig_hxyz = np.roll(self._hxyz, -1, axis=1)

        # Define dtilde at front surface for all mesh cells not on front boundary
        # Dtilde is defined differently for regions that do and don't neighbor
        # reflector regions
        boundary = self._notlast_y_accel
        boundary_grps = boundary + (slice(None),)
        D = self._diffcof[boundary_grps]
        dy = self._hxyz[boundary + (np.newaxis, 1)]
        neig_D = neig_dc[boundary_grps]
        neig_dy = neig_hxyz[boundary + (np.newaxis, 1)]
        alb = ref_albedo[boundary_grps]
        is_adj_ref_front = self._adj_reflector_front[boundary + (np.newaxis,)]
        self._dtilde[boundary_grps + (3,)] = np.where(is_adj_ref_front,
             (2.0 * D * (1.0 - alb)) / (4.0 * D * (1.0 + alb) + (1.0 - alb) * dy),
             (2.0 * D * neig_D) / (neig_dy * D + dy * neig_D))

        # Define reflector albedo for all cells on the bottom surface, in case
        # a cell borders a reflector region on the bottom
        current_in_bottom = self._current[:,:,:,_CURRENTS['in_bottom'],:]
        current_out_bottom = self._current[:,:,:,_CURRENTS['out_bottom'],:]
        ref_albedo = np.divide(current_in_bottom, current_out_bottom,
                               where=current_out_bottom > 1.0e-10,
                               out=np.ones_like(current_out_bottom))

        # Diffusion coefficient of neighbor to bottom
        neig_dc = np.roll(self._diffcof, 1, axis=2)
        # Cell dimensions of neighbor to bottom
        neig_hxyz = np.roll(self._hxyz, 1, axis=2)

        # Define dtilde at bottom surface for all mesh cells not on bottom boundary
        # Dtilde is defined differently for regions that do and don't neighbor
        # reflector regions
        boundary = self._notfirst_z_accel
        boundary_grps = boundary + (slice(None),)
        D = self._diffcof[boundary_grps]
        dz = self._hxyz[boundary + (np.newaxis, 2)]
        neig_D = neig_dc[boundary_grps]
        neig_dz = neig_hxyz[boundary + (np.newaxis, 2)]
        alb = ref_albedo[boundary_grps]
        is_adj_ref_bottom = self._adj_reflector_bottom[boundary + (np.newaxis,)]
        self._dtilde[boundary_grps + (4,)] = np.where(is_adj_ref_bottom,
             (2.0 * D * (1.0 - alb)) / (4.0 * D * (1.0 + alb) + (1.0 - alb) * dz),
             (2.0 * D * neig_D) / (neig_dz * D + dz * neig_D))

        # Define reflector albedo for all cells on the top surface, in case
        # a cell borders a reflector region on the top
        current_in_top = self._current[:,:,:,_CURRENTS['in_top'],:]
        current_out_top = self._current[:,:,:,_CURRENTS['out_top'],:]
        ref_albedo = np.divide(current_in_top, current_out_top,
                               where=current_out_top > 1.0e-10,
                               out=np.ones_like(current_out_top))

        # Diffusion coefficient of neighbor to top
        neig_dc = np.roll(self._diffcof, -1, axis=2)
        # Cell dimensions of neighbor to top
        neig_hxyz = np.roll(self._hxyz, -1, axis=2)

        # Define dtilde at top surface for all mesh cells not on top boundary
        # Dtilde is defined differently for regions that do and don't neighbor
        # reflector regions
        boundary = self._notlast_z_accel
        boundary_grps = boundary + (slice(None),)
        D = self._diffcof[boundary_grps]
        dz = self._hxyz[boundary + (np.newaxis, 2)]
        neig_D = neig_dc[boundary_grps]
        neig_dz = neig_hxyz[boundary + (np.newaxis, 2)]
        alb = ref_albedo[boundary_grps]
        is_adj_ref_top = self._adj_reflector_top[boundary + (np.newaxis,)]
        self._dtilde[boundary_grps + (5,)] = np.where(is_adj_ref_top,
             (2.0 * D * (1.0 - alb)) / (4.0 * D * (1.0 + alb) + (1.0 - alb) * dz),
             (2.0 * D * neig_D) / (neig_dz * D + dz * neig_D))

    def _compute_dhat_vectorized_updated(self):
        #TODO update comment
        # Define current in each direction
        current_in_left = self._current[:,:,:,_CURRENTS['in_left'],:]
        current_out_left = self._current[:,:,:,_CURRENTS['out_left'],:]
        current_in_right = self._current[:,:,:,_CURRENTS['in_right'],:]
        current_out_right = self._current[:,:,:,_CURRENTS['out_right'],:]
        current_in_back = self._current[:,:,:,_CURRENTS['in_back'],:]
        current_out_back = self._current[:,:,:,_CURRENTS['out_back'],:]
        current_in_front = self._current[:,:,:,_CURRENTS['in_front'],:]
        current_out_front = self._current[:,:,:,_CURRENTS['out_front'],:]
        current_in_bottom = self._current[:,:,:,_CURRENTS['in_bottom'],:]
        current_out_bottom = self._current[:,:,:,_CURRENTS['out_bottom'],:]
        current_in_top = self._current[:,:,:,_CURRENTS['in_top'],:]
        current_out_top = self._current[:,:,:,_CURRENTS['out_top'],:]

        dx = self._hxyz[:,:,:,np.newaxis,0]
        dy = self._hxyz[:,:,:,np.newaxis,1]
        dz = self._hxyz[:,:,:,np.newaxis,2]
        dxdydz = np.prod(self._hxyz, axis=3)[:,:,:,np.newaxis]

        # Define net current on each face
        net_current_left = (current_in_left - current_out_left) / dxdydz * dx
        net_current_right = (current_out_right - current_in_right) / dxdydz * dx
        net_current_back = (current_in_back - current_out_back) / dxdydz * dy
        net_current_front = (current_out_front - current_in_front) / dxdydz * dy
        net_current_bottom = (current_in_bottom - current_out_bottom) / dxdydz * dz
        net_current_top = (current_out_top - current_in_top) / dxdydz * dz

        # Define flux in each cell
        cell_flux = self._flux / dxdydz
        # Extract indices of coremap that are accelerated
        is_accel = self._coremap != _CMFD_NOACCEL

        # Define dhat at left surface for all mesh cells on left boundary
        boundary = self._first_x_accel
        boundary_grps = boundary + (slice(None),)
        net_current = net_current_left[boundary_grps]
        dtilde = self._dtilde[boundary + (slice(None), 0)]
        flux = cell_flux[boundary_grps]
        self._dhat[boundary_grps + (0,)] = (net_current + dtilde * flux) / flux

        # Define dhat at right surface for all mesh cells on right boundary
        boundary = self._last_x_accel
        boundary_grps = boundary + (slice(None),)
        net_current = net_current_right[boundary_grps]
        dtilde = self._dtilde[boundary + (slice(None), 1)]
        flux = cell_flux[boundary_grps]
        self._dhat[boundary_grps + (1,)] = (net_current - dtilde * flux) / flux

        # Define dhat at back surface for all mesh cells on back boundary
        boundary = self._first_y_accel
        boundary_grps = boundary + (slice(None),)
        net_current = net_current_back[boundary_grps]
        dtilde = self._dtilde[boundary + (slice(None), 2)]
        flux = cell_flux[boundary_grps]
        self._dhat[boundary_grps + (2,)] = (net_current + dtilde * flux) / flux

        # Define dhat at front surface for all mesh cells on front boundary
        boundary = self._last_y_accel
        boundary_grps = boundary + (slice(None),)
        net_current = net_current_front[boundary_grps]
        dtilde = self._dtilde[boundary + (slice(None), 3)]
        flux = cell_flux[boundary_grps]
        self._dhat[boundary_grps + (3,)] = (net_current - dtilde * flux) / flux

        # Define dhat at bottom surface for all mesh cells on bottom boundary
        boundary = self._first_z_accel
        boundary_grps = boundary + (slice(None),)
        net_current = net_current_bottom[boundary_grps]
        dtilde = self._dtilde[boundary + (slice(None), 4)]
        flux = cell_flux[boundary_grps]
        self._dhat[boundary_grps + (4,)] = (net_current + dtilde * flux) / flux

        # Define dhat at top surface for all mesh cells on top boundary
        boundary = self._last_z_accel
        boundary_grps = boundary + (slice(None),)
        net_current = net_current_top[boundary_grps]
        dtilde = self._dtilde[boundary + (slice(None), 5)]
        flux = cell_flux[boundary_grps]
        self._dhat[boundary_grps + (5,)] = (net_current - dtilde * flux) / flux

        # Cell flux of neighbor to left
        neig_flux = np.roll(self._flux, 1, axis=0) / dxdydz

        # Define dhat at left surface for all mesh cells not on left boundary
        # Dhat is defined differently for regions that do and don't neighbor
        # reflector regions
        boundary = self._notfirst_x_accel
        boundary_grps = boundary + (slice(None),)
        net_current = net_current_left[boundary_grps]
        dtilde = self._dtilde[boundary_grps + (0,)]
        flux = cell_flux[boundary_grps]
        flux_left = neig_flux[boundary_grps]
        is_adj_ref_left = self._adj_reflector_left[boundary + (np.newaxis,)]
        self._dhat[boundary_grps + (0,)] = np.where(is_adj_ref_left,
             (net_current + dtilde * flux) / flux,
             (net_current - dtilde * (flux_left - flux)) / (flux_left + flux))

        # Cell flux of neighbor to right
        neig_flux = np.roll(self._flux, -1, axis=0) / dxdydz

        # Define dhat at right surface for all mesh cells not on right boundary
        # Dhat is defined differently for regions that do and don't neighbor
        # reflector regions
        boundary = self._notlast_x_accel
        boundary_grps = boundary + (slice(None),)
        net_current = net_current_right[boundary_grps]
        dtilde = self._dtilde[boundary_grps + (1,)]
        flux = cell_flux[boundary_grps]
        flux_right = neig_flux[boundary_grps]
        is_adj_ref_right = self._adj_reflector_right[boundary + (np.newaxis,)]
        self._dhat[boundary_grps + (1,)] = np.where(is_adj_ref_right,
             (net_current - dtilde * flux) / flux,
             (net_current + dtilde * (flux_right - flux)) / (flux_right + flux))

        # Cell flux of neighbor to back
        neig_flux = np.roll(self._flux, 1, axis=1) / dxdydz

        # Define dhat at back surface for all mesh cells not on back boundary
        # Dhat is defined differently for regions that do and don't neighbor
        # reflector regions
        boundary = self._notfirst_y_accel
        boundary_grps = boundary + (slice(None),)
        net_current = net_current_back[boundary_grps]
        dtilde = self._dtilde[boundary_grps + (2,)]
        flux = cell_flux[boundary_grps]
        flux_back = neig_flux[boundary_grps]
        is_adj_ref_back = self._adj_reflector_back[boundary + (np.newaxis,)]
        self._dhat[boundary_grps + (2,)] = np.where(is_adj_ref_back,
             (net_current + dtilde * flux) / flux,
             (net_current - dtilde * (flux_back - flux)) / (flux_back + flux))

        # Cell flux of neighbor to front
        neig_flux = np.roll(self._flux, -1, axis=1) / dxdydz

        # Define dhat at front surface for all mesh cells not on front boundary
        # Dhat is defined differently for regions that do and don't neighbor
        # reflector regions
        boundary = self._notlast_y_accel
        boundary_grps = boundary + (slice(None),)
        net_current = net_current_front[boundary_grps]
        dtilde = self._dtilde[boundary_grps + (3,)]
        flux = cell_flux[boundary_grps]
        flux_front = neig_flux[boundary_grps]
        is_adj_ref_front = self._adj_reflector_front[boundary + (np.newaxis,)]
        self._dhat[boundary_grps + (3,)] = np.where(is_adj_ref_front,
             (net_current - dtilde * flux) / flux,
             (net_current + dtilde * (flux_front - flux)) / (flux_front + flux))

        # Cell flux of neighbor to bottom
        neig_flux = np.roll(self._flux, 1, axis=2) / dxdydz

        # Define dhat at bottom surface for all mesh cells not on bottom boundary
        # Dhat is defined differently for regions that do and don't neighbor
        # reflector regions
        boundary = self._notfirst_z_accel
        boundary_grps = boundary + (slice(None),)
        net_current = net_current_bottom[boundary_grps]
        dtilde = self._dtilde[boundary_grps + (4,)]
        flux = cell_flux[boundary_grps]
        flux_bottom = neig_flux[boundary_grps]
        is_adj_ref_bottom = self._adj_reflector_bottom[boundary + (np.newaxis,)]
        self._dhat[boundary_grps + (4,)] = np.where(is_adj_ref_bottom,
             (net_current + dtilde * flux) / flux,
             (net_current - dtilde * (flux_bottom - flux)) / (flux_bottom + flux))

        # Cell flux of neighbor to top
        neig_flux = np.roll(self._flux, -1, axis=2) / dxdydz

        # Define dhat at top surface for all mesh cells not on top boundary
        # Dhat is defined differently for regions that do and don't neighbor
        # reflector regions
        boundary = self._notlast_z_accel
        boundary_grps = boundary + (slice(None),)
        net_current = net_current_top[boundary_grps]
        dtilde = self._dtilde[boundary_grps + (5,)]
        flux = cell_flux[boundary_grps]
        flux_top = neig_flux[boundary_grps]
        is_adj_ref_top = self._adj_reflector_top[boundary + (np.newaxis,)]
        self._dhat[boundary_grps + (5,)] = np.where(is_adj_ref_top,
             (net_current - dtilde * flux) / flux,
             (net_current + dtilde * (flux_top - flux)) / (flux_top + flux))

    def _create_cmfd_tally(self):
        """Creates all tallies in-memory that are used to solve CMFD problem"""
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

    def _build_loss_matrix(self, adjoint):
        """Creates matrix representing loss of neutrons. This method uses for
        loops to loop over all spatial indices and energy groups for easier
        readability. Matrix is returned in CSR format by populating all entries
        with a numpy array and later converting to CSR format

        Parameters
        ----------
        adjoint : bool
            Whether or not to run an adjoint calculation

        Returns
        -------
        loss : scipy.sparse.spmatrix
            Sparse matrix storing elements of CMFD loss matrix

        """
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
            hxyz = self._hxyz[i,j,k,:]
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

        # Convert matrix to csr matrix in order to use scipy sparse solver
        loss = sparse.csr_matrix(loss)
        return loss

    def _build_prod_matrix(self, adjoint):
        """Creates matrix representing production of neutrons. This method uses for
        loops to loop over all spatial indices and energy groups for easier
        readability. Matrix is returned in CSR format by populating all entries
        with a numpy array and later converting to CSR format

        Parameters
        ----------
        adjoint : bool
            Whether or not to run an adjoint calculation

        Returns
        -------
        loss : scipy.sparse.spmatrix
            Sparse matrix storing elements of CMFD production matrix

        """
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

        # Convert matrix to csr matrix in order to use scipy sparse solver
        prod = sparse.csr_matrix(prod)

        return prod

    def _build_loss_matrix_vectorized(self, adjoint):
        """Creates matrix representing loss of neutrons. Since the end shape
        of the loss matrix is known a priori, this method uses a
        vectorized numpy approach to define all rows, columns, and data, which
        is then used to initialize the loss matrix in CSR format.

        Parameters
        ----------
        adjoint : bool
            Whether or not to run an adjoint calculation

        Returns
        -------
        loss : scipy.sparse.spmatrix
            Sparse matrix storing elements of CMFD loss matrix

        """
        # Extract spatial and energy indices and define matrix dimension
        ng = self._indices[3]
        n = self._mat_dim*ng

        # Define rows, columns, and data used to build csr matrix
        row = np.array([], dtype=int)
        col = np.array([], dtype=int)
        data = np.array([])

        # Define net leakage coefficient for each surface in each matrix element
        jnet = ((1.0 * self._dtilde[:,:,:,:,1] + self._dhat[:,:,:,:,1]) - \
               (-1.0 * self._dtilde[:,:,:,:,0] + self._dhat[:,:,:,:,0])) / \
               self._hxyz[:,:,:,np.newaxis,0] + \
               ((1.0 * self._dtilde[:,:,:,:,3] + self._dhat[:,:,:,:,3]) - \
               (-1.0 * self._dtilde[:,:,:,:,2] + self._dhat[:,:,:,:,2])) / \
               self._hxyz[:,:,:,np.newaxis,1] + \
               ((1.0 * self._dtilde[:,:,:,:,5] + self._dhat[:,:,:,:,5]) - \
               (-1.0 * self._dtilde[:,:,:,:,4] + self._dhat[:,:,:,:,4])) / \
               self._hxyz[:,:,:,np.newaxis,2]

        # Shift coremap in all directions to determine whether leakage term
        # should be defined for particular cell in matrix
        coremap_shift_left = np.pad(self._coremap, ((1,0),(0,0),(0,0)),
                mode='constant', constant_values=_CMFD_NOACCEL)[:-1,:,:]

        coremap_shift_right = np.pad(self._coremap, ((0,1),(0,0),(0,0)),
                mode='constant', constant_values=_CMFD_NOACCEL)[1:,:,:]

        coremap_shift_back = np.pad(self._coremap, ((0,0),(1,0),(0,0)),
                mode='constant', constant_values=_CMFD_NOACCEL)[:,:-1,:]

        coremap_shift_front = np.pad(self._coremap, ((0,0),(0,1),(0,0)),
                mode='constant', constant_values=_CMFD_NOACCEL)[:,1:,:]

        coremap_shift_bottom = np.pad(self._coremap, ((0,0),(0,0),(1,0)),
                mode='constant', constant_values=_CMFD_NOACCEL)[:,:,:-1]

        coremap_shift_top = np.pad(self._coremap, ((0,0),(0,0),(0,1)),
                mode='constant', constant_values=_CMFD_NOACCEL)[:,:,1:]

        for g in range(ng):
            # Define leakage terms that relate terms to their neighbors to the
            # left

            # Extract all regions where a cell and its neighbor to the left
            # are both fuel regions
            condition = np.logical_and(self._coremap != _CMFD_NOACCEL,
                                       coremap_shift_left != _CMFD_NOACCEL)
            idx_x = ng * (self._coremap[condition]) + g
            idx_y = ng * (coremap_shift_left[condition]) + g
            # Compute leakage term associated with these regions
            vals = (-1.0 * self._dtilde[:,:,:,g,0] -
                   self._dhat[:,:,:,g,0])[condition] / \
                   self._hxyz[:,:,:,0][condition]
            # Store rows, cols, and data to add to CSR matrix
            row = np.append(row, idx_x)
            col = np.append(col, idx_y)
            data = np.append(data, vals)

            # Define leakage terms that relate terms to their neighbors to the
            # right

            # Extract all regions where a cell and its neighbor to the right
            # are both fuel regions
            condition = np.logical_and(self._coremap != _CMFD_NOACCEL,
                                       coremap_shift_right != _CMFD_NOACCEL)
            idx_x = ng * (self._coremap[condition]) + g
            idx_y = ng * (coremap_shift_right[condition]) + g
            # Compute leakage term associated with these regions
            vals = (-1.0 * self._dtilde[:,:,:,g,1] +
                   self._dhat[:,:,:,g,1])[condition] / \
                   self._hxyz[:,:,:,0][condition]
            # Store rows, cols, and data to add to CSR matrix
            row = np.append(row, idx_x)
            col = np.append(col, idx_y)
            data = np.append(data, vals)


            # Define leakage terms that relate terms to their neighbors in the
            # back

            # Extract all regions where a cell and its neighbor in the back
            # are both fuel regions
            condition = np.logical_and(self._coremap != _CMFD_NOACCEL,
                                       coremap_shift_back != _CMFD_NOACCEL)
            idx_x = ng * (self._coremap[condition]) + g
            idx_y = ng * (coremap_shift_back[condition]) + g
            # Compute leakage term associated with these regions
            vals = (-1.0 * self._dtilde[:,:,:,g,2] -
                   self._dhat[:,:,:,g,2])[condition] / \
                   self._hxyz[:,:,:,1][condition]
            # Store rows, cols, and data to add to CSR matrix
            row = np.append(row, idx_x)
            col = np.append(col, idx_y)
            data = np.append(data, vals)

            # Define leakage terms that relate terms to their neighbors in the
            # front

            # Extract all regions where a cell and its neighbor in the front
            # are both fuel regions
            condition = np.logical_and(self._coremap != _CMFD_NOACCEL,
                                       coremap_shift_front != _CMFD_NOACCEL)
            idx_x = ng * (self._coremap[condition]) + g
            idx_y = ng * (coremap_shift_front[condition]) + g
            # Compute leakage term associated with these regions
            vals = (-1.0 * self._dtilde[:,:,:,g,3] +
                   self._dhat[:,:,:,g,3])[condition] / \
                   self._hxyz[:,:,:,1][condition]
            # Store rows, cols, and data to add to CSR matrix
            row = np.append(row, idx_x)
            col = np.append(col, idx_y)
            data = np.append(data, vals)

            # Define leakage terms that relate terms to their neighbors to the
            # bottom

            # Extract all regions where a cell and its neighbor to the bottom
            # are both fuel regions
            condition = np.logical_and(self._coremap != _CMFD_NOACCEL,
                                       coremap_shift_bottom != _CMFD_NOACCEL)
            idx_x = ng * (self._coremap[condition]) + g
            idx_y = ng * (coremap_shift_bottom[condition]) + g
            # Compute leakage term associated with these regions
            vals = (-1.0 * self._dtilde[:,:,:,g,4] -
                   self._dhat[:,:,:,g,4])[condition] / \
                   self._hxyz[:,:,:,2][condition]
            # Store rows, cols, and data to add to CSR matrix
            row = np.append(row, idx_x)
            col = np.append(col, idx_y)
            data = np.append(data, vals)

            # Define leakage terms that relate terms to their neighbors to the
            # top

            # Extract all regions where a cell and its neighbor to the
            # top are both fuel regions
            condition = np.logical_and(self._coremap != _CMFD_NOACCEL,
                                       coremap_shift_top != _CMFD_NOACCEL)
            idx_x = ng * (self._coremap[condition]) + g
            idx_y = ng * (coremap_shift_top[condition]) + g
            # Compute leakage term associated with these regions
            vals = (-1.0 * self._dtilde[:,:,:,g,5] +
                   self._dhat[:,:,:,g,5])[condition] / \
                   self._hxyz[:,:,:,2][condition]
            # Store rows, cols, and data to add to CSR matrix
            row = np.append(row, idx_x)
            col = np.append(col, idx_y)
            data = np.append(data, vals)

            # Define terms that relate to loss of neutrons in a cell. These
            # correspond to all the diagonal entries of the loss matrix

            # Extract all regions where a cell is a fuel region
            condition = self._coremap != _CMFD_NOACCEL
            idx_x = ng * (self._coremap[condition]) + g
            idx_y = idx_x
            vals = (jnet[:,:,:,g] + self._totalxs[:,:,:,g] - \
                   self._scattxs[:,:,:,g,g])[condition]
            # Store rows, cols, and data to add to CSR matrix
            row = np.append(row, idx_x)
            col = np.append(col, idx_y)
            data = np.append(data, vals)

            # Define terms that relate to in-scattering from group to group.
            # These terms relate a mesh index to all mesh indices with the same
            # spatial dimensions but belong to a different energy group
            for h in range(ng):
                if h != g:
                    # Extract all regions where a cell is a fuel region
                    condition = self._coremap != _CMFD_NOACCEL
                    idx_x = ng * (self._coremap[condition]) + g
                    idx_y = ng * (self._coremap[condition]) + h
                    # Get scattering macro xs, transposed
                    if adjoint:
                        vals = (-1.0 * self._scattxs[:, :, :, g, h])[condition]
                    # Get scattering macro xs
                    else:
                        vals = (-1.0 * self._scattxs[:, :, :, h, g])[condition]
                    # Store rows, cols, and data to add to CSR matrix
                    row = np.append(row, idx_x)
                    col = np.append(col, idx_y)
                    data = np.append(data, vals)

        # Create csr matrix
        loss = sparse.csr_matrix((data, (row, col)), shape=(n, n))
        return loss

    def _build_prod_matrix_vectorized(self, adjoint):
        """Creates matrix representing production of neutrons. Since the end shape
        of the production matrix is known a priori, this method uses a
        vectorized numpy approach to define all rows, columns, and data, which
        is then used to initialize the loss matrix in CSR format.

        Parameters
        ----------
        adjoint : bool
            Whether or not to run an adjoint calculation

        Returns
        -------
        prod : scipy.sparse.spmatrix
            Sparse matrix storing elements of CMFD production matrix

        """
        # Extract spatial and energy indices and define matrix dimension
        ng = self._indices[3]
        n = self._mat_dim*ng

        # Define rows, columns, and data used to build csr matrix
        row = np.array([], dtype=int)
        col = np.array([], dtype=int)
        data = np.array([])

        # Define terms that relate to fission production from group to group.
        for g in range(ng):
            for h in range(ng):
                # Extract all regions where a cell is a fuel region
                condition = self._coremap != _CMFD_NOACCEL
                idx_x = ng * (self._coremap[condition]) + g
                idx_y = ng * (self._coremap[condition]) + h
                # Get nu-fission macro xs, transposed
                if adjoint:
                    vals = (self._nfissxs[:, :, :, g, h])[condition]
                # Get nu-fission macro xs
                else:
                    vals = (self._nfissxs[:, :, :, h, g])[condition]
                # Store rows, cols, and data to add to CSR matrix
                row = np.append(row, idx_x)
                col = np.append(col, idx_y)
                data = np.append(data, vals)

        # Create csr matrix
        prod = sparse.csr_matrix((data, (row, col)), shape=(n, n))
        return prod

    def _matrix_to_indices(self, irow, nx, ny, nz, ng):
        """Converts matrix index in CMFD matrices to spatial and group indices
        of actual problem, based on values from coremap

        Parameters
        ----------
        irow : int
            Row in CMFD matrix
        nx : int
            Total number of mesh cells in problem in x direction
        ny : int
            Total number of mesh cells in problem in y direction
        nz : int
            Total number of mesh cells in problem in z direction
        ng : int
            Total number of energy groups in problem

        Returns
        -------
        i : int
            Corresponding x-index in CMFD problem
        j : int
            Corresponding y-index in CMFD problem
        k : int
            Corresponding z-index in CMFD problem
        g : int
            Corresponding group in CMFD problem

        """
        g = irow % ng
        # Get indices from coremap
        spatial_idx = np.where(self._coremap == int(irow/ng))
        i = spatial_idx[0][0]
        j = spatial_idx[1][0]
        k = spatial_idx[2][0]

        return i, j, k, g

    def _indices_to_matrix(self, i, j, k, g, ng):
        """Takes (i,j,k,g) indices and computes location in CMFD matrix

        Parameters
        ----------
        i : int
            Corresponding x-index in CMFD problem
        j : int
            Corresponding y-index in CMFD problem
        k : int
            Corresponding z-index in CMFD problem
        g : int
            Corresponding group in CMFD problem
        ng : int
            Total number of energy groups in problem

        Returns
        -------
        matidx : int
            Row in CMFD matrix

        """
        # Get matrix index from coremap
        matidx = ng*(self._coremap[i,j,k]) + g
        return matidx

    def _compute_dtilde_vectorized(self):
        """Computes the diffusion coupling coefficient using a vectorized numpy
        approach. Aggregate values for the dtilde multidimensional array are
        populated by first defining values on the problem boundary, and then for
        all other regions. For indices not lying on a boundary, dtilde values
        are distinguished between regions that neighbor a reflector region and
        regions that don't neighbor a reflector

        """
        # Logical for determining whether region of interest is accelerated region
        is_accel = self._coremap != _CMFD_NOACCEL
        # Logical for determining whether a zero flux "albedo" b.c. should be
        # applied
        is_zero_flux_alb = abs(self._albedo - _ZERO_FLUX) < _TINY_BIT

        # Define dtilde at left surface for all mesh cells on left boundary
        self._dtilde[0,:,:,:,0] = np.where(is_accel[0,:,:,np.newaxis],
            np.where(is_zero_flux_alb[0], 2.0 * self._diffcof[0,:,:,:] / \
                     self._hxyz[0,:,:,np.newaxis,0],
                     (2.0 * self._diffcof[0,:,:,:] * \
                     (1.0 - self._albedo[0])) / \
                     (4.0 * self._diffcof[0,:,:,:] * \
                     (1.0 + self._albedo[0]) + \
                     (1.0 - self._albedo[0]) * \
                     self._hxyz[0,:,:,np.newaxis,0])), 0)

        # Define dtilde at right surface for all mesh cells on right boundary
        self._dtilde[-1,:,:,:,1] = np.where(is_accel[-1,:,:,np.newaxis],
            np.where(is_zero_flux_alb[1], 2.0 * self._diffcof[-1,:,:,:] / \
                     self._hxyz[-1,:,:,np.newaxis,0],
                     (2.0 * self._diffcof[-1,:,:,:] * \
                     (1.0 - self._albedo[1])) / \
                     (4.0 * self._diffcof[-1,:,:,:] * \
                     (1.0 + self._albedo[1]) + \
                     (1.0 - self._albedo[1]) * \
                     self._hxyz[-1,:,:,np.newaxis,0])), 0)

        # Define dtilde at back surface for all mesh cells on back boundary
        self._dtilde[:,0,:,:,2] = np.where(is_accel[:,0,:,np.newaxis],
            np.where(is_zero_flux_alb[2], 2.0 * self._diffcof[:,0,:,:] / \
                     self._hxyz[:,0,:,np.newaxis,1],
                     (2.0 * self._diffcof[:,0,:,:] * \
                     (1.0 - self._albedo[2])) / \
                     (4.0 * self._diffcof[:,0,:,:] * \
                     (1.0 + self._albedo[2]) + \
                     (1.0 - self._albedo[2]) * \
                     self._hxyz[:,0,:,np.newaxis,1])), 0)

        # Define dtilde at front surface for all mesh cells on front boundary
        self._dtilde[:,-1,:,:,3] = np.where(is_accel[:,-1,:,np.newaxis],
            np.where(is_zero_flux_alb[3], 2.0 * self._diffcof[:,-1,:,:] / \
                     self._hxyz[:,-1,:,np.newaxis,1],
                     (2.0 * self._diffcof[:,-1,:,:] * \
                     (1.0 - self._albedo[3])) / \
                     (4.0 * self._diffcof[:,-1,:,:] * \
                     (1.0 + self._albedo[3]) + \
                     (1.0 - self._albedo[3]) * \
                     self._hxyz[:,-1,:,np.newaxis,1])), 0)

        # Define dtilde at bottom surface for all mesh cells on bottom boundary
        self._dtilde[:,:,0,:,4] = np.where(is_accel[:,:,0,np.newaxis],
            np.where(is_zero_flux_alb[4], 2.0 * self._diffcof[:,:,0,:] / \
                     self._hxyz[:,:,0,np.newaxis,2],
                     (2.0 * self._diffcof[:,:,0,:] * \
                     (1.0 - self._albedo[4])) / \
                     (4.0 * self._diffcof[:,:,0,:] * \
                     (1.0 + self._albedo[4]) + \
                     (1.0 - self._albedo[4]) * \
                     self._hxyz[:,:,0,np.newaxis,2])), 0)

        # Define dtilde at top surface for all mesh cells on top boundary
        self._dtilde[:,:,-1,:,5] = np.where(is_accel[:,:,-1,np.newaxis],
            np.where(is_zero_flux_alb[5], 2.0 * self._diffcof[:,:,-1,:] / \
                     self._hxyz[:,:,-1,np.newaxis,2],
                     (2.0 * self._diffcof[:,:,-1,:] * \
                     (1.0 - self._albedo[5])) / \
                     (4.0 * self._diffcof[:,:,-1,:] * \
                     (1.0 + self._albedo[5]) + \
                     (1.0 - self._albedo[5]) * \
                     self._hxyz[:,:,-1,np.newaxis,2])), 0)

        # Define reflector albedo for all cells on the left surface, in case
        # a cell borders a reflector region on the left
        ref_albedo = np.divide(self._current[:,:,:,_CURRENTS['in_left'],:],
                self._current[:,:,:,_CURRENTS['out_left'],:],
                where=self._current[:,:,:,_CURRENTS['out_left'],:] > 1.0e-10,
                out=np.ones_like(self._current[:,:,:,_CURRENTS['out_left'],:]))
        # Logical for whether neighboring cell to the left is reflector region
        adj_reflector = np.roll(self._coremap, 1, axis=0) == _CMFD_NOACCEL
        # Diffusion coefficient of neighbor to left
        neig_dc = np.roll(self._diffcof, 1, axis=0)
        # Cell dimensions of neighbor to left
        neig_hxyz = np.roll(self._hxyz, 1, axis=0)

        # Define dtilde at left surface for all mesh cells not on left boundary
        self._dtilde[1:,:,:,:,0] = np.where(is_accel[1:,:,:,np.newaxis], \
            np.where(adj_reflector[1:,:,:,np.newaxis],
                     (2.0 * self._diffcof[1:,:,:,:] * \
                     (1.0 - ref_albedo[1:,:,:,:])) / \
                     (4.0 * self._diffcof[1:,:,:,:] * \
                     (1.0 + ref_albedo[1:,:,:,:]) + \
                     (1.0 - ref_albedo[1:,:,:,:]) * \
                     self._hxyz[1:,:,:,np.newaxis,0]), \
                     (2.0 * self._diffcof[1:,:,:,:] * neig_dc[1:,:,:,:]) / \
                     (neig_hxyz[1:,:,:,np.newaxis,0] * self._diffcof[1:,:,:,:] + \
                     self._hxyz[1:,:,:,np.newaxis,0] * neig_dc[1:,:,:,:])), 0.0)

        # Define reflector albedo for all cells on the right surface, in case
        # a cell borders a reflector region on the right
        ref_albedo = np.divide(self._current[:,:,:,_CURRENTS['in_right'],:],
                self._current[:,:,:,_CURRENTS['out_right'],:],
                where=self._current[:,:,:,_CURRENTS['out_right'],:] > 1.0e-10,
                out=np.ones_like(self._current[:,:,:,_CURRENTS['out_right'],:]))
        # Logical for whether neighboring cell to the right is reflector region
        adj_reflector = np.roll(self._coremap, -1, axis=0) == _CMFD_NOACCEL
        # Diffusion coefficient of neighbor to right
        neig_dc = np.roll(self._diffcof, -1, axis=0)
        # Cell dimensions of neighbor to right
        neig_hxyz = np.roll(self._hxyz, -1, axis=0)

        # Define dtilde at right surface for all mesh cells not on right boundary
        self._dtilde[:-1,:,:,:,1] = np.where(is_accel[:-1,:,:,np.newaxis], \
            np.where(adj_reflector[:-1,:,:,np.newaxis],
                     (2.0 * self._diffcof[:-1,:,:,:] * \
                     (1.0 - ref_albedo[:-1,:,:,:])) / \
                     (4.0 * self._diffcof[:-1,:,:,:] * \
                     (1.0 + ref_albedo[:-1,:,:,:]) + \
                     (1.0 - ref_albedo[:-1,:,:,:]) * \
                     self._hxyz[:-1,:,:,np.newaxis,0]),
                     (2.0 * self._diffcof[:-1,:,:,:] * neig_dc[:-1,:,:,:]) / \
                     (neig_hxyz[:-1,:,:,np.newaxis,0] * self._diffcof[:-1,:,:,:] + \
                     self._hxyz[:-1,:,:,np.newaxis,0] * neig_dc[:-1,:,:,:])), 0.0)

        # Define reflector albedo for all cells on the back surface, in case
        # a cell borders a reflector region on the back
        ref_albedo = np.divide(self._current[:,:,:,_CURRENTS['in_back'],:],
                self._current[:,:,:,_CURRENTS['out_back'],:],
                where=self._current[:,:,:,_CURRENTS['out_back'],:] > 1.0e-10,
                out=np.ones_like(self._current[:,:,:,_CURRENTS['out_back'],:]))
        # Logical for whether neighboring cell to the back is reflector region
        adj_reflector = np.roll(self._coremap, 1, axis=1) == _CMFD_NOACCEL
        # Diffusion coefficient of neighbor to back
        neig_dc = np.roll(self._diffcof, 1, axis=1)
        # Cell dimensions of neighbor to back
        neig_hxyz = np.roll(self._hxyz, 1, axis=1)

        # Define dtilde at back surface for all mesh cells not on back boundary
        self._dtilde[:,1:,:,:,2] = np.where(is_accel[:,1:,:,np.newaxis], \
            np.where(adj_reflector[:,1:,:,np.newaxis],
                     (2.0 * self._diffcof[:,1:,:,:] * \
                     (1.0 - ref_albedo[:,1:,:,:])) / \
                     (4.0 * self._diffcof[:,1:,:,:] * \
                     (1.0 + ref_albedo[:,1:,:,:]) + \
                     (1.0 - ref_albedo[:,1:,:,:]) * \
                     self._hxyz[:,1:,:,np.newaxis,1]),
                     (2.0 * self._diffcof[:,1:,:,:] * neig_dc[:,1:,:,:]) / \
                     (neig_hxyz[:,1:,:,np.newaxis,1] * self._diffcof[:,1:,:,:] + \
                     self._hxyz[:,1:,:,np.newaxis,1] * neig_dc[:,1:,:,:])), 0.0)

        # Define reflector albedo for all cells on the front surface, in case
        # a cell borders a reflector region in the front
        ref_albedo = np.divide(self._current[:,:,:,_CURRENTS['in_front'],:],
                self._current[:,:,:,_CURRENTS['out_front'],:],
                where=self._current[:,:,:,_CURRENTS['out_front'],:] > 1.0e-10,
                out=np.ones_like(self._current[:,:,:,_CURRENTS['out_front'],:]))
        # Logical for whether neighboring cell to the front is reflector region
        adj_reflector = np.roll(self._coremap, -1, axis=1) == _CMFD_NOACCEL
        # Diffusion coefficient of neighbor to front
        neig_dc = np.roll(self._diffcof, -1, axis=1)
        # Cell dimensions of neighbor to front
        neig_hxyz = np.roll(self._hxyz, -1, axis=1)

        # Define dtilde at front surface for all mesh cells not on front boundary
        self._dtilde[:,:-1,:,:,3] = np.where(is_accel[:,:-1,:,np.newaxis], \
            np.where(adj_reflector[:,:-1,:,np.newaxis],
                     (2.0 * self._diffcof[:,:-1,:,:] * \
                     (1.0 - ref_albedo[:,:-1,:,:])) / \
                     (4.0 * self._diffcof[:,:-1,:,:] * \
                     (1.0 + ref_albedo[:,:-1,:,:]) + \
                     (1.0 - ref_albedo[:,:-1,:,:]) * \
                     self._hxyz[:,:-1,:,np.newaxis,1]),
                     (2.0 * self._diffcof[:,:-1,:,:] * neig_dc[:,:-1,:,:]) / \
                     (neig_hxyz[:,:-1,:,np.newaxis,1] * self._diffcof[:,:-1,:,:] + \
                     self._hxyz[:,:-1,:,np.newaxis,1] * neig_dc[:,:-1,:,:])), 0.0)

        # Define reflector albedo for all cells on the bottom surface, in case
        # a cell borders a reflector region on the bottom
        ref_albedo = np.divide(self._current[:,:,:,_CURRENTS['in_bottom'],:],
                self._current[:,:,:,_CURRENTS['out_bottom'],:],
                where=self._current[:,:,:,_CURRENTS['out_bottom'],:] > 1.0e-10,
                out=np.ones_like(self._current[:,:,:,_CURRENTS['out_bottom'],:]))
        # Logical for whether neighboring cell to the bottom is reflector region
        adj_reflector = np.roll(self._coremap, 1, axis=2) == _CMFD_NOACCEL
        # Diffusion coefficient of neighbor to bottom
        neig_dc = np.roll(self._diffcof, 1, axis=2)
        # Cell dimensions of neighbor to bottom
        neig_hxyz = np.roll(self._hxyz, 1, axis=2)

        # Define dtilde at bottom surface for all mesh cells not on bottom boundary
        self._dtilde[:,:,1:,:,4] = np.where(is_accel[:,:,1:,np.newaxis], \
            np.where(adj_reflector[:,:,1:,np.newaxis],
                     (2.0 * self._diffcof[:,:,1:,:] * \
                     (1.0 - ref_albedo[:,:,1:,:])) / \
                     (4.0 * self._diffcof[:,:,1:,:] * \
                     (1.0 + ref_albedo[:,:,1:,:]) + \
                     (1.0 - ref_albedo[:,:,1:,:]) * \
                     self._hxyz[:,:,1:,np.newaxis,2]),
                     (2.0 * self._diffcof[:,:,1:,:] * neig_dc[:,:,1:,:]) / \
                     (neig_hxyz[:,:,1:,np.newaxis,2] * self._diffcof[:,:,1:,:] + \
                     self._hxyz[:,:,1:,np.newaxis,2] * neig_dc[:,:,1:,:])), 0.0)

        # Define reflector albedo for all cells on the top surface, in case
        # a cell borders a reflector region on the top
        ref_albedo = np.divide(self._current[:,:,:,_CURRENTS['in_top'],:],
                self._current[:,:,:,_CURRENTS['out_top'],:],
                where=self._current[:,:,:,_CURRENTS['out_top'],:] > 1.0e-10,
                out=np.ones_like(self._current[:,:,:,_CURRENTS['out_top'],:]))
        # Logical for whether neighboring cell to the top is reflector region
        adj_reflector = np.roll(self._coremap, -1, axis=2) == _CMFD_NOACCEL
        # Diffusion coefficient of neighbor to top
        neig_dc = np.roll(self._diffcof, -1, axis=2)
        # Cell dimensions of neighbor to top
        neig_hxyz = np.roll(self._hxyz, -1, axis=2)

        # Define dtilde at top surface for all mesh cells not on top boundary
        self._dtilde[:,:,:-1,:,5] = np.where(is_accel[:,:,:-1,np.newaxis], \
            np.where(adj_reflector[:,:,:-1,np.newaxis],
                     (2.0 * self._diffcof[:,:,:-1,:] * \
                     (1.0 - ref_albedo[:,:,:-1,:])) / \
                     (4.0 * self._diffcof[:,:,:-1,:] * \
                     (1.0 + ref_albedo[:,:,:-1,:]) + \
                     (1.0 - ref_albedo[:,:,:-1,:]) * \
                     self._hxyz[:,:,:-1,np.newaxis,2]),
                     (2.0 * self._diffcof[:,:,:-1,:] * neig_dc[:,:,:-1,:]) / \
                     (neig_hxyz[:,:,:-1,np.newaxis,2] * self._diffcof[:,:,:-1,:] + \
                     self._hxyz[:,:,:-1,np.newaxis,2] * neig_dc[:,:,:-1,:])), 0.0)

    def _compute_dtilde(self):
        """Computes the diffusion coupling coefficient by looping over all
        spatial regions and energy groups

        """
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
                        cell_hxyz = self._hxyz[i,j,k,:]

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
                                neig_hxyz = self._hxyz[tuple(neig_idx)]

                                # Check for fuel-reflector interface
                                if (self._coremap[tuple(neig_idx)] ==
                                         _CMFD_NOACCEL):
                                    # Get albedo
                                    ref_albedo = self._get_reflector_albedo(l,g,i,j,k)
                                    dtilde = (2*cell_dc*(1-ref_albedo))/(4*cell_dc*(1+ \
                                             ref_albedo)+(1-ref_albedo)*cell_hxyz[xyz_idx])

                                else:  # Not next to a reflector
                                    # Compute dtilde to neighbor cell
                                    dtilde = (2*cell_dc*neig_dc)/(neig_hxyz[xyz_idx]*cell_dc + \
                                             cell_hxyz[xyz_idx]*neig_dc)

                            # Record dtilde
                            self._dtilde[i, j, k, g, l] = dtilde

    def _compute_dhat_vectorized(self):
        """Computes the nonlinear coupling coefficient using a vectorized numpy
        approach. Aggregate values for the dhat multidimensional array are
        populated by first defining values on the problem boundary, and then for
        all other regions. For indices not lying by a boundary, dhat values
        are distinguished between regions that neighbor a reflector region and
        regions that don't neighbor a reflector

        """
        # Define net current on each face, divided by surface area
        net_current_left = ((self._current[:,:,:,_CURRENTS['in_left'],:] - \
                             self._current[:,:,:,_CURRENTS['out_left'],:]) / \
                             np.prod(self._hxyz, axis=3)[:,:,:,np.newaxis] * \
                             self._hxyz[:,:,:,np.newaxis,0])
        net_current_right = ((self._current[:,:,:,_CURRENTS['out_right'],:] - \
                             self._current[:,:,:,_CURRENTS['in_right'],:]) / \
                             np.prod(self._hxyz, axis=3)[:,:,:,np.newaxis] * \
                             self._hxyz[:,:,:,np.newaxis,0])
        net_current_back = ((self._current[:,:,:,_CURRENTS['in_back'],:] - \
                             self._current[:,:,:,_CURRENTS['out_back'],:]) / \
                             np.prod(self._hxyz, axis=3)[:,:,:,np.newaxis] * \
                             self._hxyz[:,:,:,np.newaxis,1])
        net_current_front = ((self._current[:,:,:,_CURRENTS['out_front'],:] - \
                             self._current[:,:,:,_CURRENTS['in_front'],:]) / \
                             np.prod(self._hxyz, axis=3)[:,:,:,np.newaxis] * \
                             self._hxyz[:,:,:,np.newaxis,1])
        net_current_bottom = ((self._current[:,:,:,_CURRENTS['in_bottom'],:] - \
                             self._current[:,:,:,_CURRENTS['out_bottom'],:]) / \
                             np.prod(self._hxyz, axis=3)[:,:,:,np.newaxis] * \
                             self._hxyz[:,:,:,np.newaxis,2])
        net_current_top = ((self._current[:,:,:,_CURRENTS['out_top'],:] - \
                             self._current[:,:,:,_CURRENTS['in_top'],:]) / \
                             np.prod(self._hxyz, axis=3)[:,:,:,np.newaxis] * \
                             self._hxyz[:,:,:,np.newaxis,2])

        # Define flux in each cell
        cell_flux = self._flux / np.prod(self._hxyz, axis=3)[:,:,:,np.newaxis]
        # Extract indices of coremap that are accelerated
        is_accel = self._coremap != _CMFD_NOACCEL

        # Define dhat at left surface for all mesh cells on left boundary
        self._dhat[0,:,:,:,0] = np.where(is_accel[0,:,:,np.newaxis],
            (net_current_left[0,:,:,:] + self._dtilde[0,:,:,:,0] * \
            cell_flux[0,:,:,:]) / cell_flux[0,:,:,:], 0)
        # Define dhat at right surface for all mesh cells on right boundary
        self._dhat[-1,:,:,:,1] = np.where(is_accel[-1,:,:,np.newaxis],
            (net_current_right[-1,:,:,:] - self._dtilde[-1,:,:,:,1] * \
            cell_flux[-1,:,:,:]) / cell_flux[-1,:,:,:], 0)
        # Define dhat at back surface for all mesh cells on back boundary
        self._dhat[:,0,:,:,2] = np.where(is_accel[:,0,:,np.newaxis],
            (net_current_back[:,0,:,:] + self._dtilde[:,0,:,:,2] * \
            cell_flux[:,0,:,:]) / cell_flux[:,0,:,:], 0)
        # Define dhat at front surface for all mesh cells on front boundary
        self._dhat[:,-1,:,:,3] = np.where(is_accel[:,-1,:,np.newaxis],
            (net_current_front[:,-1,:,:] - self._dtilde[:,-1,:,:,3] * \
            cell_flux[:,-1,:,:]) / cell_flux[:,-1,:,:], 0)
        # Define dhat at bottom surface for all mesh cells on bottom boundary
        self._dhat[:,:,0,:,4] = np.where(is_accel[:,:,0,np.newaxis],
            (net_current_bottom[:,:,0,:] + self._dtilde[:,:,0,:,4] * \
            cell_flux[:,:,0,:]) / cell_flux[:,:,0,:], 0)
        # Define dhat at top surface for all mesh cells on top boundary
        self._dhat[:,:,-1,:,5] = np.where(is_accel[:,:,-1,np.newaxis],
            (net_current_top[:,:,-1,:] - self._dtilde[:,:,-1,:,5] * \
            cell_flux[:,:,-1,:]) / cell_flux[:,:,-1,:], 0)

        # Logical for whether neighboring cell to the left is reflector region
        adj_reflector = np.roll(self._coremap, 1, axis=0) == _CMFD_NOACCEL
        # Cell flux of neighbor to left
        neig_flux = np.roll(self._flux, 1, axis=0) / \
                    np.prod(self._hxyz, axis=3)[:,:,:,np.newaxis]
        # Define dhat at left surface for all mesh cells not on left boundary
        self._dhat[1:,:,:,:,0] = np.where(is_accel[1:,:,:,np.newaxis], \
            np.where(adj_reflector[1:,:,:,np.newaxis],
                     (net_current_left[1:,:,:,:] + self._dtilde[1:,:,:,:,0] * \
                     cell_flux[1:,:,:,:]) / cell_flux[1:,:,:,:],
                     (net_current_left[1:,:,:,:] - self._dtilde[1:,:,:,:,0] * \
                     (neig_flux[1:,:,:,:] - cell_flux[1:,:,:,:])) / \
                     (neig_flux[1:,:,:,:] + cell_flux[1:,:,:,:])), 0.0)

        # Logical for whether neighboring cell to the right is reflector region
        adj_reflector = np.roll(self._coremap, -1, axis=0) == _CMFD_NOACCEL
        # Cell flux of neighbor to right
        neig_flux = np.roll(self._flux, -1, axis=0) / \
                    np.prod(self._hxyz, axis=3)[:,:,:,np.newaxis]
        # Define dhat at right surface for all mesh cells not on right boundary
        self._dhat[:-1,:,:,:,1] = np.where(is_accel[:-1,:,:,np.newaxis], \
            np.where(adj_reflector[:-1,:,:,np.newaxis],
                     (net_current_right[:-1,:,:,:] - self._dtilde[:-1,:,:,:,1] * \
                     cell_flux[:-1,:,:,:]) / cell_flux[:-1,:,:,:],
                     (net_current_right[:-1,:,:,:] + self._dtilde[:-1,:,:,:,1] * \
                     (neig_flux[:-1,:,:,:] - cell_flux[:-1,:,:,:])) / \
                     (neig_flux[:-1,:,:,:] + cell_flux[:-1,:,:,:])), 0.0)

        # Logical for whether neighboring cell to the back is reflector region
        adj_reflector = np.roll(self._coremap, 1, axis=1) == _CMFD_NOACCEL
        # Cell flux of neighbor to back
        neig_flux = np.roll(self._flux, 1, axis=1) / \
                    np.prod(self._hxyz, axis=3)[:,:,:,np.newaxis]
        # Define dhat at back surface for all mesh cells not on back boundary
        self._dhat[:,1:,:,:,2] = np.where(is_accel[:,1:,:,np.newaxis], \
            np.where(adj_reflector[:,1:,:,np.newaxis],
                     (net_current_back[:,1:,:,:] + self._dtilde[:,1:,:,:,2] * \
                     cell_flux[:,1:,:,:]) / cell_flux[:,1:,:,:],
                     (net_current_back[:,1:,:,:] - self._dtilde[:,1:,:,:,2] * \
                     (neig_flux[:,1:,:,:] - cell_flux[:,1:,:,:])) / \
                     (neig_flux[:,1:,:,:] + cell_flux[:,1:,:,:])), 0.0)

        # Logical for whether neighboring cell to the front is reflector region
        adj_reflector = np.roll(self._coremap, -1, axis=1) == _CMFD_NOACCEL
        # Cell flux of neighbor to front
        neig_flux = np.roll(self._flux, -1, axis=1) / \
                    np.prod(self._hxyz, axis=3)[:,:,:,np.newaxis]
        # Define dhat at front surface for all mesh cells not on front boundary
        self._dhat[:,:-1,:,:,3] = np.where(is_accel[:,:-1,:,np.newaxis], \
            np.where(adj_reflector[:,:-1,:,np.newaxis],
                     (net_current_front[:,:-1,:,:] - self._dtilde[:,:-1,:,:,3] * \
                     cell_flux[:,:-1,:,:]) / cell_flux[:,:-1,:,:],
                     (net_current_front[:,:-1,:,:] + self._dtilde[:,:-1,:,:,3] * \
                     (neig_flux[:,:-1,:,:] - cell_flux[:,:-1,:,:])) / \
                     (neig_flux[:,:-1,:,:] + cell_flux[:,:-1,:,:])), 0.0)

        # Logical for whether neighboring cell to the bottom is reflector region
        adj_reflector = np.roll(self._coremap, 1, axis=2) == _CMFD_NOACCEL
        # Cell flux of neighbor to bottom
        neig_flux = np.roll(self._flux, 1, axis=2) / \
                    np.prod(self._hxyz, axis=3)[:,:,:,np.newaxis]
        # Define dhat at bottom surface for all mesh cells not on bottom boundary
        self._dhat[:,:,1:,:,4] = np.where(is_accel[:,:,1:,np.newaxis], \
            np.where(adj_reflector[:,:,1:,np.newaxis],
                     (net_current_bottom[:,:,1:,:] + self._dtilde[:,:,1:,:,4] * \
                     cell_flux[:,:,1:,:]) / cell_flux[:,:,1:,:],
                     (net_current_bottom[:,:,1:,:] - self._dtilde[:,:,1:,:,4] * \
                     (neig_flux[:,:,1:,:] - cell_flux[:,:,1:,:])) / \
                     (neig_flux[:,:,1:,:] + cell_flux[:,:,1:,:])), 0.0)

        # Logical for whether neighboring cell to the top is reflector region
        adj_reflector = np.roll(self._coremap, -1, axis=2) == _CMFD_NOACCEL
        # Cell flux of neighbor to top
        neig_flux = np.roll(self._flux, -1, axis=2) / \
                    np.prod(self._hxyz, axis=3)[:,:,:,np.newaxis]
        # Define dhat at top surface for all mesh cells not on top boundary
        self._dhat[:,:,:-1,:,5] = np.where(is_accel[:,:,:-1,np.newaxis], \
            np.where(adj_reflector[:,:,:-1,np.newaxis],
                     (net_current_top[:,:,:-1,:] - self._dtilde[:,:,:-1,:,5] * \
                     cell_flux[:,:,:-1,:]) / cell_flux[:,:,:-1,:],
                     (net_current_top[:,:,:-1,:] + self._dtilde[:,:,:-1,:,5] * \
                     (neig_flux[:,:,:-1,:] - cell_flux[:,:,:-1,:])) / \
                     (neig_flux[:,:,:-1,:] + cell_flux[:,:,:-1,:])), 0.0)

    def _compute_dhat(self):
        """Computes the nonlinear coupling coefficient by looping over all
        spatial regions and energy groups

        """
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
                        cell_flux = self._flux[i,j,k,g]/np.product(self._hxyz[i,j,k,:])
                        current = self._current[i,j,k,:,g]

                        # Setup of vector to identify boundary conditions
                        bound = np.repeat([i,j,k], 2)

                        # Begin loop around sides of cell for leakage
                        for l in range(6):
                            xyz_idx = int(l/2)  # x=0, y=1, z=2
                            dir_idx = l % 2     # -=0, +=1
                            shift_idx = 2*(l % 2) - 1 # shift neig by -1 or +1

                            # Calculate net current on l face (divided by surf area)
                            net_current = shift_idx*(current[2*l] - current[2*l+1]) / \
                                np.product(self._hxyz[i,j,k,:]) * self._hxyz[i,j,k,xyz_idx]

                            # Check if at boundary
                            if bound[l] == nxyz[xyz_idx, dir_idx]:
                                # Compute dhat
                                dhat = (net_current - shift_idx*cell_dtilde[l]*cell_flux) / \
                                     cell_flux

                            else:  # Not at a boundary
                                # Compute neighboring cell indices
                                neig_idx = [i,j,k]  # Begin with i,j,k
                                neig_idx[xyz_idx] += shift_idx

                                # Get neigbor flux
                                neig_flux = self._flux[tuple(neig_idx)+(g,)] / \
                                            np.product(self._hxyz[tuple(neig_idx)])

                                # Check for fuel-reflector interface
                                if (self._coremap[tuple(neig_idx)] ==
                                         _CMFD_NOACCEL):
                                    # Compute dhat
                                    dhat = (net_current - shift_idx*cell_dtilde[l]*cell_flux) / \
                                         cell_flux

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
                openmc.capi.master():
            print(' Dhats reset to zero')
            sys.stdout.flush()

    def _get_reflector_albedo(self, l, g, i, j, k):
        """Calculates the albedo to the reflector by returning ratio of
        incoming to outgoing ratio

        Parameters
        ----------
        l : int
            leakage index, see _CURRENTS for map between leakage index and leakage
            direction
        g : int
            index of energy group
        i : int
            index of mesh location in x-direction
        j : int
            index of mesh location in y-direction
        k : int
            index of mesh location in z-direction

        Returns
        -------
        float
            reflector albedo

        """
        # Get partial currents from object
        current = self._current[i,j,k,:,g]

        # Calculate albedo
        if current[2*l] < 1.0e-10:
            return 1.0
        else:
          return current[2*l+1]/current[2*l]


