module settings

  use, intrinsic :: ISO_C_BINDING

  use constants
  use set_header,    only: SetInt

  implicit none

  ! ============================================================================
  ! ENERGY TREATMENT RELATED VARIABLES
  logical :: run_CE = .true.  ! Run in CE mode?

  ! ============================================================================
  ! CONTINUOUS-ENERGY CROSS SECTION RELATED VARIABLES

  ! Unreoslved resonance probablity tables
  logical :: urr_ptables_on = .true.

  ! Default temperature and method for choosing temperatures
  integer :: temperature_method = TEMPERATURE_NEAREST
  logical :: temperature_multipole = .false.
  real(8) :: temperature_tolerance = 10.0_8
  real(8) :: temperature_default = 293.6_8
  real(8) :: temperature_range(2) = [ZERO, ZERO]

  integer :: n_log_bins  ! number of bins for logarithmic grid

  ! ============================================================================
  ! MULTI-GROUP CROSS SECTION RELATED VARIABLES

  ! Maximum Data Order
  integer :: max_order

  ! Whether or not to convert Legendres to tabulars
  logical :: legendre_to_tabular = .true.

  ! Number of points to use in the Legendre to tabular conversion
  integer :: legendre_to_tabular_points = 33

  ! Assume all tallies are spatially distinct
  logical :: assume_separate = .false.

  ! Use confidence intervals for results instead of standard deviations
  logical :: confidence_intervals = .false.

  ! ============================================================================
  ! SIMULATION VARIABLES

  integer(C_INT64_T), bind(C) :: n_particles = 0   ! # of particles per generation
  integer(C_INT32_T), bind(C) :: n_batches         ! # of batches
  integer(C_INT32_T), bind(C) :: n_inactive        ! # of inactive batches
  integer(C_INT32_T), bind(C) :: gen_per_batch = 1 ! # of generations per batch

  integer :: n_max_batches             ! max # of batches
  integer :: n_batch_interval = 1      ! batch interval for triggers
  logical :: pred_batches = .false.    ! predict batches for triggers
  logical :: trigger_on = .false.      ! flag for turning triggers on/off

  logical :: entropy_on = .false.
  integer :: index_entropy_mesh = -1

  logical :: ufs = .false.
  integer :: index_ufs_mesh = -1

  ! Write source at end of simulation
  logical :: source_separate = .false.
  logical :: source_write = .true.
  logical :: source_latest = .false.

  ! Variance reduction settins
  logical :: survival_biasing = .false.
  real(8) :: weight_cutoff = 0.25_8
  real(8) :: energy_cutoff = ZERO
  real(8) :: weight_survive = ONE

  ! Mode to run in (fixed source, eigenvalue, plotting, etc)
  integer(C_INT), bind(C) :: run_mode = NONE

  ! Restart run
  logical :: restart_run = .false.

  ! The verbosity controls how much information will be printed to the screen
  ! and in logs
  integer(C_INT), bind(C) :: verbosity = 7

  logical :: check_overlaps = .false.

  ! Trace for single particle
  integer    :: trace_batch
  integer    :: trace_gen
  integer(8) :: trace_particle

  ! Particle tracks
  logical :: write_all_tracks = .false.
  integer, allocatable :: track_identifiers(:,:)

  ! Particle restart run
  logical :: particle_restart_run = .false.

  ! Write out initial source
  logical :: write_initial_source = .false.

  ! Whether create fission neutrons or not. Only applied for MODE_FIXEDSOURCE
  logical :: create_fission_neutrons = .true.

  ! Information about state points to be written
  integer :: n_state_points = 0
  type(SetInt) :: statepoint_batch

  ! Information about source points to be written
  integer :: n_source_points = 0
  type(SetInt) :: sourcepoint_batch

  character(MAX_FILE_LEN) :: path_input               ! Path to input file
  character(MAX_FILE_LEN) :: path_cross_sections = '' ! Path to cross_sections.xml
  character(MAX_FILE_LEN) :: path_multipole           ! Path to wmp library
  character(MAX_FILE_LEN) :: path_source = ''         ! Path to binary source
  character(MAX_FILE_LEN) :: path_state_point         ! Path to binary state point
  character(MAX_FILE_LEN) :: path_source_point        ! Path to binary source point
  character(MAX_FILE_LEN) :: path_particle_restart    ! Path to particle restart
  character(MAX_FILE_LEN) :: path_output = ''         ! Path to output directory

  ! Various output options
  logical :: output_summary = .true.
  logical :: output_tallies = .true.

  ! Resonance scattering settings
  logical :: res_scat_on = .false. ! is resonance scattering treated?
  integer :: res_scat_method = RES_SCAT_ARES  ! resonance scattering method
  real(8) :: res_scat_energy_min = 0.01_8
  real(8) :: res_scat_energy_max = 1000.0_8
  character(10), allocatable :: res_scat_nuclides(:)

  ! Is CMFD active
  logical :: cmfd_run = .false.

  ! No reduction at end of batch
  logical :: reduce_tallies = .true.

contains

!===============================================================================
! FREE_MEMORY_SETTINGS deallocates global arrays defined in this module
!===============================================================================

  subroutine free_memory_settings()
    if (allocated(res_scat_nuclides)) deallocate(res_scat_nuclides)
    if (allocated(track_identifiers)) deallocate(track_identifiers)

    call statepoint_batch % clear()
    call sourcepoint_batch % clear()
  end subroutine free_memory_settings

end module settings
