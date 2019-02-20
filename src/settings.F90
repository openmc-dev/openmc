module settings

  use, intrinsic :: ISO_C_BINDING

  use constants

  implicit none

  ! ============================================================================
  ! CONTINUOUS-ENERGY CROSS SECTION RELATED VARIABLES

  logical(C_BOOL), bind(C) :: photon_transport

  ! ============================================================================
  ! SIMULATION VARIABLES

  ! Assume all tallies are spatially distinct
  logical(C_BOOL), bind(C) :: assume_separate

  integer(C_INT64_T), bind(C) :: n_particles   ! # of particles per generation
  integer(C_INT32_T), bind(C) :: n_inactive    ! # of inactive batches

  logical(C_BOOL), bind(C) :: trigger_on      ! flag for turning triggers on/off

  ! Mode to run in (fixed source, eigenvalue, plotting, etc)
  integer(C_INT), bind(C) :: run_mode

  ! flag for use of DAGMC geometry
  logical(C_BOOL), bind(C) :: dagmc

  ! The verbosity controls how much information will be printed to the screen
  ! and in logs
  integer(C_INT), bind(C) :: verbosity

  character(MAX_FILE_LEN) :: path_input               ! Path to input file
  character(MAX_FILE_LEN) :: path_cross_sections = '' ! Path to cross_sections.xml

  ! No reduction at end of batch
  logical(C_BOOL), bind(C) :: reduce_tallies

end module settings
