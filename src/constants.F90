module constants

  use, intrinsic :: ISO_C_BINDING

  implicit none

  ! ============================================================================
  ! VERSIONING NUMBERS

  ! OpenMC major, minor, and release numbers
  integer, parameter :: VERSION_MAJOR   = 0
  integer, parameter :: VERSION_MINOR   = 10
  integer, parameter :: VERSION_RELEASE = 0
  integer, parameter :: &
       VERSION(3) = [VERSION_MAJOR, VERSION_MINOR, VERSION_RELEASE]

  ! HDF5 data format
  integer, parameter :: HDF5_VERSION(2) = [1, 0]

  ! Version numbers for binary files
  integer, parameter :: VERSION_STATEPOINT(2)       = [17, 0]
  integer, parameter :: VERSION_PARTICLE_RESTART(2) = [2, 0]
  integer, parameter :: VERSION_TRACK(2)            = [2, 0]
  integer, parameter :: VERSION_SUMMARY(2)          = [5, 0]
  integer, parameter :: VERSION_VOLUME(2)           = [1, 0]
  integer, parameter :: VERSION_VOXEL(2)            = [1, 0]
  integer, parameter :: VERSION_MGXS_LIBRARY(2)     = [1, 0]
  character(10), parameter :: VERSION_MULTIPOLE     = "v0.2"

  ! ============================================================================
  ! ADJUSTABLE PARAMETERS

  ! NOTE: This is the only section of the constants module that should ever be
  ! adjusted. Modifying constants in other sections may cause the code to fail.

  ! Monoatomic ideal-gas scattering treatment threshold
  real(8), parameter :: FREE_GAS_THRESHOLD = 400.0_8

  ! Significance level for confidence intervals
  real(8), parameter :: CONFIDENCE_LEVEL = 0.95_8

  ! Used for surface current tallies
  real(8), parameter :: TINY_BIT = 1e-8_8

  ! User for precision in geometry
  real(C_DOUBLE), bind(C, name='FP_PRECISION') :: FP_PRECISION = 1e-14_8
  real(C_DOUBLE), bind(C, name='FP_REL_PRECISION') :: FP_REL_PRECISION = 1e-5_8
  real(C_DOUBLE), bind(C, name='FP_COINCIDENT') :: FP_COINCIDENT = 1e-12_8

  ! Maximum number of collisions/crossings
  integer, parameter :: MAX_EVENTS = 1000000
  integer, parameter :: MAX_SAMPLE = 100000

  ! Maximum number of secondary particles created
  integer, parameter :: MAX_SECONDARY = 1000

  ! Maximum number of words in a single line, length of line, and length of
  ! single word
  integer, parameter :: MAX_WORDS       = 500
  integer, parameter :: MAX_LINE_LEN    = 250
  integer, parameter :: MAX_WORD_LEN    = 150
  integer, parameter :: MAX_FILE_LEN    = 255

  ! Maximum number of external source spatial resamples to encounter before an
  ! error is thrown.
  integer, parameter :: EXTSRC_REJECT_THRESHOLD = 10000
  real(8), parameter :: EXTSRC_REJECT_FRACTION = 0.05

  ! ============================================================================
  ! PHYSICAL CONSTANTS

  ! Values here are from the Committee on Data for Science and Technology
  ! (CODATA) 2014 recommendation (doi:10.1103/RevModPhys.88.035009).

  real(8), parameter ::                      &
       PI               = 3.1415926535898_8, & ! pi
       SQRT_PI          = 1.7724538509055_8, & ! square root of pi
       MASS_NEUTRON     = 1.00866491588_8,   & ! mass of a neutron in amu
       MASS_NEUTRON_EV  = 939.5654133e6_8,   & ! mass of a neutron in eV/c^2
       MASS_PROTON      = 1.007276466879_8,  & ! mass of a proton in amu
       AMU              = 1.660539040e-27_8, & ! 1 amu in kg
       C_LIGHT          = 2.99792458e8_8,    & ! speed of light in m/s
       N_AVOGADRO       = 0.6022140857_8,    & ! Avogadro's number in 10^24/mol
       K_BOLTZMANN      = 8.6173303e-5_8,    & ! Boltzmann constant in eV/K
       INFINITY         = huge(0.0_8),       & ! positive infinity
       ZERO             = 0.0_8,             &
       HALF             = 0.5_8,             &
       ONE              = 1.0_8,             &
       TWO              = 2.0_8,             &
       THREE            = 3.0_8,             &
       FOUR             = 4.0_8
  complex(8), parameter :: ONEI = (ZERO, ONE)

  ! ============================================================================
  ! GEOMETRY-RELATED CONSTANTS

  ! Boundary conditions
  integer(C_INT), bind(C, name="BC_TRANSMIT") :: BC_TRANSMIT
  integer(C_INT), bind(C, name="BC_VACUUM") :: BC_VACUUM
  integer(C_INT), bind(C, name="BC_REFLECT") :: BC_REFLECT
  integer(C_INT), bind(C, name="BC_PERIODIC") :: BC_PERIODIC

  ! Logical operators for cell definitions
  integer, parameter ::              &
       OP_LEFT_PAREN   = huge(0),     & ! Left parentheses
       OP_RIGHT_PAREN  = huge(0) - 1, & ! Right parentheses
       OP_COMPLEMENT   = huge(0) - 2, & ! Complement operator (~)
       OP_INTERSECTION = huge(0) - 3, & ! Intersection operator
       OP_UNION        = huge(0) - 4    ! Union operator (^)

  ! Cell fill types
  integer, parameter ::  &
       FILL_MATERIAL = 1, & ! Cell with a specified material
       FILL_UNIVERSE = 2, & ! Cell filled by a separate universe
       FILL_LATTICE  = 3    ! Cell filled with a lattice

  ! Void material
  integer, parameter :: MATERIAL_VOID = -1

  ! Lattice types
  integer, parameter ::  &
       LATTICE_RECT = 1, & ! Rectangular lattice
       LATTICE_HEX  = 2    ! Hexagonal lattice

  ! Lattice boundary crossings
  integer, parameter ::    &
       LATTICE_LEFT   = 1, & ! Flag for crossing left (x) lattice boundary
       LATTICE_RIGHT  = 2, & ! Flag for crossing right (x) lattice boundary
       LATTICE_BACK   = 3, & ! Flag for crossing back (y) lattice boundary
       LATTICE_FRONT  = 4, & ! Flag for crossing front (y) lattice boundary
       LATTICE_BOTTOM = 5, & ! Flag for crossing bottom (z) lattice boundary
       LATTICE_TOP    = 6    ! Flag for crossing top (z) lattice boundary

  ! Surface types
  integer, parameter ::  &
       SURF_PX     =  1, & ! Plane parallel to x-plane
       SURF_PY     =  2, & ! Plane parallel to y-plane
       SURF_PZ     =  3, & ! Plane parallel to z-plane
       SURF_PLANE  =  4, & ! Arbitrary plane
       SURF_CYL_X  =  5, & ! Cylinder along x-axis
       SURF_CYL_Y  =  6, & ! Cylinder along y-axis
       SURF_CYL_Z  =  7, & ! Cylinder along z-axis
       SURF_SPHERE =  8, & ! Sphere
       SURF_CONE_X =  9, & ! Cone parallel to x-axis
       SURF_CONE_Y = 10, & ! Cone parallel to y-axis
       SURF_CONE_Z = 11    ! Cone parallel to z-axis

  ! Flag to say that the outside of a lattice is not defined
  integer, parameter :: NO_OUTER_UNIVERSE = -1

  ! Maximum number of lost particles
  integer, parameter :: MAX_LOST_PARTICLES = 10

  ! Maximum number of lost particles, relative to the total number of particles
  real(8), parameter :: REL_MAX_LOST_PARTICLES = 1e-6_8

  ! ============================================================================
  ! CROSS SECTION RELATED CONSTANTS

  ! Interpolation flag
  integer, parameter ::   &
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
       ANGLE_ISOTROPIC = 1, & ! Isotropic angular distribution (CE)
       ANGLE_32_EQUI   = 2, & ! 32 equiprobable bins (CE)
       ANGLE_TABULAR   = 3, & ! Tabular angular distribution (CE or MG)
       ANGLE_LEGENDRE  = 4, & ! Legendre angular distribution (MG)
       ANGLE_HISTOGRAM = 5    ! Histogram angular distribution (MG)

  ! Number of mu bins to use when converting Legendres to tabular type
  integer, parameter :: DEFAULT_NMU = 33

  ! Secondary energy mode for S(a,b) inelastic scattering
  integer, parameter :: &
       SAB_SECONDARY_EQUAL  = 0, & ! Equally-likely outgoing energy bins
       SAB_SECONDARY_SKEWED = 1, & ! Skewed outgoing energy bins
       SAB_SECONDARY_CONT   = 2    ! Continuous, linear-linear interpolation

  ! Elastic mode for S(a,b) elastic scattering
  integer, parameter :: &
       SAB_ELASTIC_DISCRETE = 3, & ! Sample from discrete cosines
       SAB_ELASTIC_EXACT    = 4    ! Exact treatment for coherent elastic

  ! Reaction types
  integer, parameter :: &
       TOTAL_XS = 1,  ELASTIC = 2,   N_NONELASTIC = 3, N_LEVEL = 4, MISC  = 5,   &
       N_2ND   = 11,  N_2N    = 16,  N_3N   = 17,  N_FISSION = 18, N_F    = 19,  &
       N_NF    = 20,  N_2NF   = 21,  N_NA   = 22,  N_N3A   = 23,  N_2NA   = 24,  &
       N_3NA   = 25,  N_NP    = 28,  N_N2A  = 29,  N_2N2A  = 30,  N_ND    = 32,  &
       N_NT    = 33,  N_N3HE  = 34,  N_ND2A = 35,  N_NT2A  = 36,  N_4N    = 37,  &
       N_3NF   = 38,  N_2NP   = 41,  N_3NP  = 42,  N_N2P   = 44,  N_NPA   = 45,  &
       N_N1    = 51,  N_N40   = 90,  N_NC   = 91,  N_DISAPPEAR = 101, N_GAMMA = 102, &
       N_P     = 103, N_D     = 104, N_T    = 105, N_3HE   = 106, N_A     = 107, &
       N_2A    = 108, N_3A    = 109, N_2P   = 111, N_PA    = 112, N_T2A   = 113, &
       N_D2A   = 114, N_PD    = 115, N_PT   = 116, N_DA    = 117, N_5N    = 152, &
       N_6N    = 153, N_2NT   = 154, N_TA   = 155, N_4NP   = 156, N_3ND   = 157, &
       N_NDA   = 158, N_2NPA  = 159, N_7N   = 160, N_8N    = 161, N_5NP   = 162, &
       N_6NP   = 163, N_7NP   = 164, N_4NA  = 165, N_5NA   = 166, N_6NA   = 167, &
       N_7NA   = 168, N_4ND   = 169, N_5ND  = 170, N_6ND   = 171, N_3NT   = 172, &
       N_4NT   = 173, N_5NT   = 174, N_6NT  = 175, N_2N3HE = 176, N_3N3HE = 177, &
       N_4N3HE = 178, N_3N2P  = 179, N_3N3A = 180, N_3NPA  = 181, N_DT    = 182, &
       N_NPD   = 183, N_NPT   = 184, N_NDT  = 185, N_NP3HE = 186, N_ND3HE = 187, &
       N_NT3HE = 188, N_NTA   = 189, N_2N2P = 190, N_P3HE  = 191, N_D3HE  = 192, &
       N_3HEA  = 193, N_4N2P  = 194, N_4N2A = 195, N_4NPA  = 196, N_3P    = 197, &
       N_N3P   = 198, N_3N2PA = 199, N_5N2P = 200, N_P0    = 600, N_PC    = 649, &
       N_D0    = 650, N_DC    = 699, N_T0   = 700, N_TC    = 749, N_3HE0  = 750, &
       N_3HEC  = 799, N_A0    = 800, N_AC   = 849, N_2N0   = 875, N_2NC   = 891

  ! Depletion reactions
  integer, parameter :: DEPLETION_RX(6) = [N_2N, N_3N, N_4N, N_GAMMA, N_P, N_A]

  ! ACE table types
  integer, parameter :: &
       ACE_NEUTRON   = 1, & ! continuous-energy neutron
       ACE_THERMAL   = 2, & ! thermal S(a,b) scattering data
       ACE_DOSIMETRY = 3    ! dosimetry cross sections

  ! MGXS Table Types
  integer, parameter :: &
       MGXS_ISOTROPIC   = 1, & ! Isotropically Weighted Data
       MGXS_ANGLE       = 2    ! Data by Angular Bins

  ! Fission neutron emission (nu) type
  integer, parameter ::   &
       NU_NONE       = 0, & ! No nu values (non-fissionable)
       NU_POLYNOMIAL = 1, & ! Nu values given by polynomial
       NU_TABULAR    = 2    ! Nu values given by tabular distribution

  ! Secondary particle emission type
  integer, parameter :: &
       EMISSION_PROMPT = 1,  & ! Prompt emission of secondary particle
       EMISSION_DELAYED = 2, & ! Delayed emission of secondary particle
       EMISSION_TOTAL = 3      ! Yield represents total emission (prompt + delayed)

  ! Cross section filetypes
  integer, parameter :: &
       ASCII  = 1, & ! ASCII cross section file
       BINARY = 2    ! Binary cross section file

  ! Library types
  integer, parameter :: &
       LIBRARY_NEUTRON = 1, &
       LIBRARY_THERMAL = 2, &
       LIBRARY_PHOTON = 3, &
       LIBRARY_MULTIGROUP = 4

  ! Probability table parameters
  integer, parameter :: &
       URR_CUM_PROB = 1, &
       URR_TOTAL    = 2, &
       URR_ELASTIC  = 3, &
       URR_FISSION  = 4, &
       URR_N_GAMMA  = 5, &
       URR_HEATING  = 6

  ! Maximum number of partial fission reactions
  integer, parameter :: PARTIAL_FISSION_MAX = 4

  ! Temperature treatment method
  integer, parameter :: &
       TEMPERATURE_NEAREST = 1, &
       TEMPERATURE_INTERPOLATION = 2

  ! Resonance elastic scattering methods
  integer, parameter :: &
       RES_SCAT_ARES = 1, &
       RES_SCAT_DBRC = 2, &
       RES_SCAT_WCM = 3, &
       RES_SCAT_CXS = 4

  ! ============================================================================
  ! TALLY-RELATED CONSTANTS

  ! Tally result entries
  integer, parameter :: &
       RESULT_VALUE  = 1, &
       RESULT_SUM    = 2, &
       RESULT_SUM_SQ = 3

  ! Tally type
  integer, parameter :: &
       TALLY_VOLUME          = 1, &
       TALLY_MESH_CURRENT    = 2, &
       TALLY_SURFACE         = 3

  ! Tally estimator types
  integer, parameter :: &
       ESTIMATOR_ANALOG      = 1, &
       ESTIMATOR_TRACKLENGTH = 2, &
       ESTIMATOR_COLLISION   = 3

  ! Event types for tallies
  integer, parameter :: &
       EVENT_SURFACE = -2, &
       EVENT_LATTICE = -1, &
       EVENT_SCATTER =  1, &
       EVENT_ABSORB  =  2

  ! Tally score type -- if you change these, make sure you also update the
  ! _SCORES dictionary in openmc/capi/tally.py
  integer, parameter :: N_SCORE_TYPES = 24
  integer, parameter :: &
       SCORE_FLUX               = -1,  & ! flux
       SCORE_TOTAL              = -2,  & ! total reaction rate
       SCORE_SCATTER            = -3,  & ! scattering rate
       SCORE_NU_SCATTER         = -4,  & ! scattering production rate
       SCORE_SCATTER_N          = -5,  & ! arbitrary scattering moment
       SCORE_SCATTER_PN         = -6,  & ! system for scoring 0th through nth moment
       SCORE_NU_SCATTER_N       = -7,  & ! arbitrary nu-scattering moment
       SCORE_NU_SCATTER_PN      = -8,  & ! system for scoring 0th through nth nu-scatter moment
       SCORE_ABSORPTION         = -9,  & ! absorption rate
       SCORE_FISSION            = -10, & ! fission rate
       SCORE_NU_FISSION         = -11, & ! neutron production rate
       SCORE_KAPPA_FISSION      = -12, & ! fission energy production rate
       SCORE_CURRENT            = -13, & ! current
       SCORE_FLUX_YN            = -14, & ! angular moment of flux
       SCORE_TOTAL_YN           = -15, & ! angular moment of total reaction rate
       SCORE_SCATTER_YN         = -16, & ! angular flux-weighted scattering moment (0:N)
       SCORE_NU_SCATTER_YN      = -17, & ! angular flux-weighted nu-scattering moment (0:N)
       SCORE_EVENTS             = -18, & ! number of events
       SCORE_DELAYED_NU_FISSION = -19, & ! delayed neutron production rate
       SCORE_PROMPT_NU_FISSION  = -20, & ! prompt neutron production rate
       SCORE_INVERSE_VELOCITY   = -21, & ! flux-weighted inverse velocity
       SCORE_FISS_Q_PROMPT      = -22, & ! prompt fission Q-value
       SCORE_FISS_Q_RECOV       = -23, & ! recoverable fission Q-value
       SCORE_DECAY_RATE         = -24    ! delayed neutron precursor decay rate

  ! Maximum scattering order supported
  integer, parameter :: MAX_ANG_ORDER = 10

  ! Names of *-PN & *-YN scores (MOMENT_STRS) and *-N moment scores
  character(*), parameter :: &
       MOMENT_STRS(6)    = (/ "scatter-p   ",   &
                              "nu-scatter-p",   &
                              "flux-y      ",   &
                              "total-y     ",   &
                              "scatter-y   ",   &
                              "nu-scatter-y"/), &
       MOMENT_N_STRS(2)  = (/ "scatter-    ",   &
                              "nu-scatter- "/)

  ! Location in MOMENT_STRS where the YN data begins
  integer, parameter :: YN_LOC = 3

  ! Tally map bin finding
  integer, parameter :: NO_BIN_FOUND = -1

  ! Tally filter and map types
  integer, parameter :: N_FILTER_TYPES = 15
  integer, parameter :: &
       FILTER_UNIVERSE       = 1,  &
       FILTER_MATERIAL       = 2,  &
       FILTER_CELL           = 3,  &
       FILTER_CELLBORN       = 4,  &
       FILTER_SURFACE        = 5,  &
       FILTER_MESH           = 6,  &
       FILTER_ENERGYIN       = 7,  &
       FILTER_ENERGYOUT      = 8,  &
       FILTER_DISTRIBCELL    = 9,  &
       FILTER_MU             = 10, &
       FILTER_POLAR          = 11, &
       FILTER_AZIMUTHAL      = 12, &
       FILTER_DELAYEDGROUP   = 13, &
       FILTER_ENERGYFUNCTION = 14, &
       FILTER_CELLFROM       = 15

  ! Mesh types
  integer, parameter :: &
       MESH_REGULAR = 1

  ! Tally surface current directions
  integer, parameter :: &
       OUT_LEFT   = 1,   &   ! x min
       IN_LEFT    = 2,   &   ! x min
       OUT_RIGHT  = 3,   &   ! x max
       IN_RIGHT   = 4,   &   ! x max
       OUT_BACK   = 5,   &   ! y min
       IN_BACK    = 6,   &   ! y min
       OUT_FRONT  = 7,   &   ! y max
       IN_FRONT   = 8,   &   ! y max
       OUT_BOTTOM = 9,   &   ! z min
       IN_BOTTOM  = 10,  &   ! z min
       OUT_TOP    = 11,  &   ! z max
       IN_TOP     = 12       ! z max

  ! Tally trigger types and threshold
  integer, parameter :: &
       VARIANCE           = 1, &
       RELATIVE_ERROR     = 2, &
       STANDARD_DEVIATION = 3

  ! Global tally parameters
  integer, parameter :: N_GLOBAL_TALLIES = 4
  integer, parameter :: &
       K_COLLISION   = 1, &
       K_ABSORPTION  = 2, &
       K_TRACKLENGTH = 3, &
       LEAKAGE       = 4

  ! Differential tally independent variables
  integer, parameter :: &
       DIFF_DENSITY = 1, &
       DIFF_NUCLIDE_DENSITY = 2, &
       DIFF_TEMPERATURE = 3

  ! ============================================================================
  ! RANDOM NUMBER STREAM CONSTANTS

  integer(C_INT), bind(C, name='N_STREAMS') :: N_STREAMS
  integer(C_INT), bind(C, name='STREAM_TRACKING') :: STREAM_TRACKING
  integer(C_INT), bind(C, name='STREAM_TALLIES') :: STREAM_TALLIES
  integer(C_INT), bind(C, name='STREAM_SOURCE') :: STREAM_SOURCE
  integer(C_INT), bind(C, name='STREAM_URR_PTABLE') :: STREAM_URR_PTABLE
  integer(C_INT), bind(C, name='STREAM_VOLUME') :: STREAM_VOLUME
  integer(C_INT64_T), parameter :: DEFAULT_SEED = 1_8

  ! ============================================================================
  ! MISCELLANEOUS CONSTANTS

  ! indicates that an array index hasn't been set
  integer, parameter :: NONE = 0

  ! Codes for read errors -- better hope these numbers are never used in an
  ! input file!
  integer, parameter :: ERROR_INT  = -huge(0)
  real(8), parameter :: ERROR_REAL = -huge(0.0_8) * 0.917826354_8

  ! Running modes
  integer, parameter ::        &
       MODE_FIXEDSOURCE = 1, & ! Fixed source mode
       MODE_EIGENVALUE  = 2, & ! K eigenvalue mode
       MODE_PLOTTING    = 3, & ! Plotting mode
       MODE_PARTICLE    = 4, & ! Particle restart mode
       MODE_VOLUME      = 5    ! Volume calculation mode

  !=============================================================================
  ! CMFD CONSTANTS

  ! for non-accelerated regions on coarse mesh overlay
  integer, parameter :: CMFD_NOACCEL = 99999

  ! constant to represent a zero flux "albedo"
  real(8), parameter :: ZERO_FLUX = 999.0_8

  ! constant for writing out no residual
  real(8), parameter :: CMFD_NORES = 99999.0_8

  !=============================================================================
  ! DELAYED NEUTRON PRECURSOR CONSTANTS

  ! Since cross section libraries come with different numbers of delayed groups
  ! (e.g. ENDF/B-VII.1 has 6 and JEFF 3.1.1 has 8 delayed groups) and we don't
  ! yet know what cross section library is being used when the tallies.xml file
  ! is read in, we want to have an upper bound on the size of the array we
  ! use to store the bins for delayed group tallies.
  integer, parameter :: MAX_DELAYED_GROUPS = 8

end module constants
