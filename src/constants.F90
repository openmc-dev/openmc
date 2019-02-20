module constants

  use, intrinsic :: ISO_C_BINDING

  implicit none

  ! ============================================================================
  ! VERSIONING NUMBERS

  ! Version numbers for binary files
  integer, parameter :: VERSION_TRACK(2) = [2, 0]

  ! ============================================================================
  ! ADJUSTABLE PARAMETERS

  ! NOTE: This is the only section of the constants module that should ever be
  ! adjusted. Modifying constants in other sections may cause the code to fail.

  ! Maximum number of secondary particles created
  integer, parameter :: MAX_SECONDARY = 1000

  ! Maximum number of words in a single line, length of line, and length of
  ! single word
  integer, parameter :: MAX_LINE_LEN    = 250
  integer, parameter :: MAX_WORD_LEN    = 150
  integer, parameter :: MAX_FILE_LEN    = 255

  ! ============================================================================
  ! PHYSICAL CONSTANTS

  ! Values here are from the Committee on Data for Science and Technology
  ! (CODATA) 2014 recommendation (doi:10.1103/RevModPhys.88.035009).

  real(8), parameter ::                      &
       PI               = 3.1415926535898_8, & ! pi
       INFINITY         = huge(0.0_8),       & ! positive infinity
       ZERO             = 0.0_8,             &
       HALF             = 0.5_8,             &
       ONE              = 1.0_8,             &
       TWO              = 2.0_8

  ! ============================================================================
  ! GEOMETRY-RELATED CONSTANTS

  ! Boundary conditions
  integer(C_INT), bind(C, name="BC_TRANSMIT") :: BC_TRANSMIT
  integer(C_INT), bind(C, name="BC_VACUUM") :: BC_VACUUM
  integer(C_INT), bind(C, name="BC_REFLECT") :: BC_REFLECT
  integer(C_INT), bind(C, name="BC_PERIODIC") :: BC_PERIODIC

  ! Cell fill types
  integer, parameter ::  &
       FILL_MATERIAL = 1, & ! Cell with a specified material
       FILL_UNIVERSE = 2, & ! Cell filled by a separate universe
       FILL_LATTICE  = 3    ! Cell filled with a lattice
  integer(C_INT), bind(C, name='FILL_MATERIAL') :: FILL_MATERIAL_C = FILL_MATERIAL
  integer(C_INT), bind(C, name='FILL_UNIVERSE') :: FILL_UNIVERSE_C = FILL_UNIVERSE
  integer(C_INT), bind(C, name='FILL_LATTICE')  :: FILL_LATTICE_C  = FILL_LATTICE

  ! Void material
  integer, parameter :: MATERIAL_VOID = -1

  ! ============================================================================
  ! CROSS SECTION RELATED CONSTANTS

  ! Particle type
  integer, parameter :: &
       NEUTRON  = 0, &
       PHOTON   = 1, &
       ELECTRON = 2, &
       POSITRON = 3

  ! MGXS Table Types
  integer, parameter :: &
       MGXS_ISOTROPIC   = 1, & ! Isotropically Weighted Data
       MGXS_ANGLE       = 2    ! Data by Angular Bins

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
       TALLY_MESH_SURFACE    = 2, &
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
  integer, parameter :: N_SCORE_TYPES = 16
  integer, parameter :: &
       SCORE_FLUX               = -1,  & ! flux
       SCORE_TOTAL              = -2,  & ! total reaction rate
       SCORE_SCATTER            = -3,  & ! scattering rate
       SCORE_NU_SCATTER         = -4,  & ! scattering production rate
       SCORE_ABSORPTION         = -5,  & ! absorption rate
       SCORE_FISSION            = -6,  & ! fission rate
       SCORE_NU_FISSION         = -7,  & ! neutron production rate
       SCORE_KAPPA_FISSION      = -8,  & ! fission energy production rate
       SCORE_CURRENT            = -9,  & ! current
       SCORE_EVENTS             = -10, & ! number of events
       SCORE_DELAYED_NU_FISSION = -11, & ! delayed neutron production rate
       SCORE_PROMPT_NU_FISSION  = -12, & ! prompt neutron production rate
       SCORE_INVERSE_VELOCITY   = -13, & ! flux-weighted inverse velocity
       SCORE_FISS_Q_PROMPT      = -14, & ! prompt fission Q-value
       SCORE_FISS_Q_RECOV       = -15, & ! recoverable fission Q-value
       SCORE_DECAY_RATE         = -16    ! delayed neutron precursor decay rate

  ! Tally map bin finding
  integer, parameter :: NO_BIN_FOUND = -1

  ! Tally filter and map types
  integer, parameter :: &
       FILTER_UNIVERSE       = 1,  &
       FILTER_MATERIAL       = 2,  &
       FILTER_CELL           = 3

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
  ! MISCELLANEOUS CONSTANTS

  ! indicates that an array index hasn't been set
  integer, parameter :: NONE = 0
  integer, parameter :: C_NONE = -1

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
  ! DELAYED NEUTRON PRECURSOR CONSTANTS

  ! Since cross section libraries come with different numbers of delayed groups
  ! (e.g. ENDF/B-VII.1 has 6 and JEFF 3.1.1 has 8 delayed groups) and we don't
  ! yet know what cross section library is being used when the tallies.xml file
  ! is read in, we want to have an upper bound on the size of the array we
  ! use to store the bins for delayed group tallies.
  integer, parameter :: MAX_DELAYED_GROUPS = 8

end module constants
