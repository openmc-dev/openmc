module constants

  implicit none

  ! ============================================================================
  ! VERSIONING NUMBERS

  integer, parameter :: VERSION_MAJOR   = 0
  integer, parameter :: VERSION_MINOR   = 3
  integer, parameter :: VERSION_RELEASE = 3

  ! ============================================================================
  ! ADJUSTABLE PARAMETERS 

  ! NOTE: This is the only section of the constants module that should ever be
  ! adjusted. Modifying constants in other sections may cause the code to fail.

  ! Monoatomic ideal-gas scattering treatment threshold
  real(8), parameter :: FREE_GAS_THRESHOLD = 400.0

  ! Used for surface current tallies
  real(8), parameter :: TINY_BIT = 1e-8_8

  ! User for precision in geometry
  real(8), parameter :: FP_PRECISION = 1e-7_8

  ! Maximum number of collisions/crossings
  integer, parameter :: MAX_EVENTS = 10000
  integer, parameter :: MAX_SAMPLE = 100000

  ! Maximum number of words in a single line, length of line, and length of
  ! single word
  integer, parameter :: MAX_WORDS    = 500
  integer, parameter :: MAX_LINE_LEN = 250
  integer, parameter :: MAX_WORD_LEN = 150
  integer, parameter :: MAX_FILE_LEN = 255

  ! ============================================================================
  ! PHYSICAL CONSTANTS

  real(8), parameter ::            &
       PI           = 3.1415926535898_8, & ! pi
       MASS_NEUTRON = 1.0086649156,      & ! mass of a neutron
       MASS_PROTON  = 1.00727646677,     & ! mass of a proton
       AMU          = 1.66053873e-27,    & ! 1 amu in kg
       N_AVOGADRO   = 0.602214179,       & ! Avogadro's number in 10^24/mol
       K_BOLTZMANN  = 8.617342e-11,      & ! Boltzmann constant in MeV/K
       INFINITY     = huge(0.0_8),       & ! positive infinity
       ZERO         = 0.0_8,             &
       ONE          = 1.0_8,             &
       TWO          = 2.0_8

  ! ============================================================================
  ! GEOMETRY-RELATED CONSTANTS

  ! Boundary conditions
  integer, parameter ::    &
       & BC_TRANSMIT = 0,  & ! Transmission boundary condition (default)
       & BC_VACUUM   = 1,  & ! Vacuum boundary condition
       & BC_REFLECT  = 2,  & ! Reflecting boundary condition
       & BC_PERIODIC = 3     ! Periodic boundary condition

  ! Logical operators for cell definitions
  integer, parameter ::                &
       & OP_LEFT_PAREN  = huge(0),     & ! Left parentheses
       & OP_RIGHT_PAREN = huge(0) - 1, & ! Right parentheses
       & OP_UNION       = huge(0) - 2, & ! Union operator
       & OP_DIFFERENCE  = huge(0) - 3    ! Difference operator

  ! Cell types
  integer, parameter ::    &
       & CELL_NORMAL  = 1, & ! Cell with a specified material
       & CELL_FILL    = 2, & ! Cell filled by a separate universe
       & CELL_LATTICE = 3, & ! Cell filled with a lattice
       & CELL_VOID    = -1

  ! Lattice types
  integer, parameter :: &
       & LATTICE_RECT = 1, &
       & LATTICE_HEX  = 2

  ! Surface types
  integer, parameter ::    &
       & SURF_PX     =  1, & ! Plane parallel to x-plane 
       & SURF_PY     =  2, & ! Plane parallel to y-plane 
       & SURF_PZ     =  3, & ! Plane parallel to z-plane 
       & SURF_PLANE  =  4, & ! Arbitrary plane
       & SURF_CYL_X  =  5, & ! Cylinder along x-axis
       & SURF_CYL_Y  =  6, & ! Cylinder along y-axis
       & SURF_CYL_Z  =  7, & ! Cylinder along z-axis
       & SURF_SPHERE =  8, & ! Sphere
       & SURF_BOX_X  =  9, & ! Box extending infinitely in x-direction
       & SURF_BOX_Y  = 10, & ! Box extending infinitely in y-direction
       & SURF_BOX_Z  = 11, & ! Box extending infinitely in z-direction
       & SURF_BOX    = 12, & ! Rectangular prism
       & SURF_GQ     = 13    ! General quadratic surface

  ! Surface senses
  integer, parameter ::      &
       & SENSE_POSITIVE = 1, &
       & SENSE_NEGATIVE = -1

  ! ============================================================================
  ! CROSS SECTION RELATED CONSTANTS

  ! Interpolation flag
  integer, parameter ::     &
       HISTOGRAM     = 1, & ! y is constant in x
       LINEAR_LINEAR = 2, & ! y is linear in x
       LINEAR_LOG    = 3, & ! y is linear in ln(x)
       LOG_LINEAR    = 4, & ! ln(y) is linear in x
       LOG_LOG       = 5    ! ln(y) is linear in ln(x)

  ! Particle type
  integer, parameter :: &
       NEUTRON  = 1, &
       PHOTON   = 2, &
       ELECTRON = 3

  ! Angular distribution type
  integer, parameter :: & 
       & ANGLE_ISOTROPIC = 1, & ! Isotropic angular distribution
       & ANGLE_32_EQUI   = 2, & ! 32 equiprobable bins
       & ANGLE_TABULAR   = 3    ! Tabular angular distribution

  ! Secondary energy mode for S(a,b) inelastic scattering
  integer, parameter :: &
       & SAB_SECONDARY_EQUAL  = 0, & ! Equally-likely outgoing energy bins
       & SAB_SECONDARY_SKEWED = 1    ! Skewed outgoing energy bins

  ! Elastic mode for S(a,b) elastic scattering
  integer, parameter :: &
       & SAB_ELASTIC_DISCRETE = 3, & ! Sample from discrete cosines
       & SAB_ELASTIC_EXACT    = 4    ! Exact treatment for coherent elastic

  ! Reaction types
  integer, parameter :: &
       TOTAL_XS    = 1, &
       ELASTIC     = 2, &
       N_LEVEL     = 4, &
       MISC        = 5, &
       N_2ND       = 11, &
       N_2N        = 16, &
       N_3N        = 17, &
       N_FISSION   = 18, &
       N_F         = 19, &
       N_NF        = 20, &
       N_2NF       = 21, &
       N_NA        = 22, &
       N_N3A       = 23, &
       N_2NA       = 24, &
       N_3NA       = 25, &
       N_NP        = 28, &
       N_N2A       = 29, &
       N_2N2A      = 30, &
       N_ND        = 32, &
       N_NT        = 33, &
       N_N3HE      = 34, &
       N_ND2A      = 35, &
       N_NT2A      = 36, &
       N_4N        = 37, &
       N_3NF       = 38, &
       N_2NP       = 41, &
       N_3NP       = 42, &
       N_N2P       = 44, &
       N_NPA       = 45, &
       N_N1        = 51, &
       N_N40       = 90, &
       N_NC        = 91, &
       N_DISAPPEAR = 101, &
       N_GAMMA     = 102, &
       N_P         = 103, &
       N_D         = 104, &
       N_T         = 105, &
       N_3HE       = 106, &
       N_A         = 107, &
       N_2A        = 108, &
       N_3A        = 109, &
       N_2P        = 111, &
       N_PA        = 112, &
       N_T2A       = 113, &
       N_D2A       = 114, &
       N_PD        = 115, &
       N_PT        = 116, &
       N_DA        = 117

  ! ACE table types
  integer, parameter :: &
       ACE_NEUTRON   = 1, & ! continuous-energy neutron
       ACE_THERMAL   = 2, & ! thermal S(a,b) scattering data
       ACE_DOSIMETRY = 3    ! dosimetry cross sections

  ! Fission neutron emission (nu) type
  integer, parameter ::   &
       NU_NONE       = 0, & ! No nu values (non-fissionable)
       NU_POLYNOMIAL = 1, & ! Nu values given by polynomial
       NU_TABULAR    = 2    ! Nu values given by tabular distribution

  ! Cross section filetypes
  integer, parameter :: &
       ASCII  = 1, & ! ASCII cross section file
       BINARY = 2    ! Binary cross section file

  ! Maximum number of partial fission reactions
  integer, parameter :: PARTIAL_FISSION_MAX = 4

  ! ============================================================================
  ! TALLY-RELATED CONSTANTS

  ! Tally macro reactions
  integer, parameter :: N_MACRO_TYPES = 15
  integer, parameter :: &
       MACRO_FLUX       = -1,  & ! flux
       MACRO_TOTAL      = -2,  & ! total reaction rate
       MACRO_SCATTER    = -3,  & ! scattering rate
       MACRO_NU_SCATTER = -4,  & ! scattering production rate
       MACRO_SCATTER_1  = -5,  & ! first scattering moment
       MACRO_SCATTER_2  = -6,  & ! second scattering moment
       MACRO_SCATTER_3  = -7,  & ! third scattering moment
       MACRO_N_1N       = -8,  & ! (n,1n) rate
       MACRO_N_2N       = -9,  & ! (n,2n) rate
       MACRO_N_3N       = -10, & ! (n,3n) rate
       MACRO_N_4N       = -11, & ! (n,4n) rate
       MACRO_ABSORPTION = -12, & ! absorption rate
       MACRO_FISSION    = -13, & ! fission rate
       MACRO_NU_FISSION = -14, & ! neutron production rate
       MACRO_CURRENT    = -15    ! partial current

  ! Tally map bin finding
  integer, parameter :: NO_BIN_FOUND = -1

  ! Tally filter and map types
  integer, parameter :: TALLY_TYPES = 8
  integer, parameter :: &
       T_UNIVERSE  = 1, &
       T_MATERIAL  = 2, &
       T_CELL      = 3, &
       T_CELLBORN  = 4, &
       T_SURFACE   = 5, &
       T_MESH      = 6, &
       T_ENERGYIN  = 7, &
       T_ENERGYOUT = 8

  ! Filter types for surface current tallies
  integer, parameter :: &
       TS_MESH_X   = 1, &
       TS_MESH_Y   = 2, &
       TS_MESH_Z   = 3, &
       TS_ENERGYIN = 4, &
       TS_SURFACE  = 5

  ! Tally surface current directions
  integer, parameter :: &
       IN_RIGHT  = 1,   &
       OUT_RIGHT = 2,   &
       IN_FRONT  = 3,   &
       OUT_FRONT = 4,   &
       IN_TOP    = 5,   &
       OUT_TOP   = 6

  ! ============================================================================
  ! MISCELLANEOUS CONSTANTS

  ! indicates that an array index hasn't been set
  integer, parameter :: NONE = 0

  ! Codes for read errors -- better hope these numbers are never used in an
  ! input file!
  integer, parameter :: ERROR_INT  = -huge(0)
  real(8), parameter :: ERROR_REAL = -huge(0.0_8) * 0.917826354_8

  ! Source types
  integer, parameter ::   &
       SRC_BOX     = 1, & ! Source in a rectangular prism
       SRC_CELL    = 2, & ! Source in a cell
       SRC_SURFACE = 3    ! Source on a surface

  integer, parameter ::        &
       PROB_SOURCE      = 1, & ! External source problem
       PROB_CRITICALITY = 2    ! Criticality problem

  ! Unit numbers
  integer, parameter :: UNIT_LOG   = 11 ! unit # for writing log file
  integer, parameter :: UNIT_TALLY = 12 ! unit # for writing tally file
  integer, parameter :: UNIT_PLOT  = 13 ! unit # for writing plot file

end module constants

