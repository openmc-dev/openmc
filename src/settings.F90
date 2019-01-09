module settings

  use, intrinsic :: ISO_C_BINDING

  use constants

  implicit none

  ! ============================================================================
  ! ENERGY TREATMENT RELATED VARIABLES
  logical(C_BOOL), bind(C, name='run_CE') :: run_CE  ! Run in CE mode?

  ! ============================================================================
  ! CONTINUOUS-ENERGY CROSS SECTION RELATED VARIABLES

  ! Default temperature and method for choosing temperatures
  integer(C_INT), bind(C) :: temperature_method
  logical(C_BOOL), bind(C) :: temperature_multipole
  real(C_DOUBLE), bind(C) :: temperature_tolerance
  real(C_DOUBLE), bind(C) :: temperature_default
  real(C_DOUBLE), bind(C) :: temperature_range(2)

  integer(C_INT), bind(C) :: n_log_bins  ! number of bins for logarithmic grid

  logical(C_BOOL), bind(C) :: photon_transport
  integer(C_INT), bind(C) :: electron_treatment

  ! ============================================================================
  ! MULTI-GROUP CROSS SECTION RELATED VARIABLES

  ! Maximum Data Order
  integer(C_INT), bind(C) :: max_order

  ! Whether or not to convert Legendres to tabulars
  logical(C_BOOL), bind(C) :: legendre_to_tabular

  ! Number of points to use in the Legendre to tabular conversion
  integer(C_INT), bind(C) :: legendre_to_tabular_points

  ! ============================================================================
  ! SIMULATION VARIABLES

  ! Assume all tallies are spatially distinct
  logical(C_BOOL), bind(C) :: assume_separate

  ! Use confidence intervals for results instead of standard deviations
  logical(C_BOOL), bind(C) :: confidence_intervals

  integer(C_INT64_T), bind(C) :: n_particles   ! # of particles per generation
  integer(C_INT32_T), bind(C) :: n_batches     ! # of batches
  integer(C_INT32_T), bind(C) :: n_inactive    ! # of inactive batches
  integer(C_INT32_T), bind(C) :: gen_per_batch ! # of generations per batch

  integer(C_INT), bind(C) :: n_max_batches             ! max # of batches
  integer(C_INT), bind(C, name='trigger_batch_interval') :: n_batch_interval      ! batch interval for triggers
  logical(C_BOOL), bind(C, name='trigger_predict') :: pred_batches   ! predict batches for triggers
  logical(C_BOOL), bind(C) :: trigger_on      ! flag for turning triggers on/off

  logical(C_BOOL), bind(C) :: entropy_on

  ! Write source at end of simulation

  ! Variance reduction settins
  logical(C_BOOL), bind(C) :: survival_biasing
  real(C_DOUBLE), bind(C) :: energy_cutoff(4)

  ! Mode to run in (fixed source, eigenvalue, plotting, etc)
  integer(C_INT), bind(C) :: run_mode

  ! flag for use of DAGMC geometry
  logical(C_BOOL), bind(C) :: dagmc

  ! Restart run
  logical(C_BOOL), bind(C) :: restart_run

  ! The verbosity controls how much information will be printed to the screen
  ! and in logs
  integer(C_INT), bind(C) :: verbosity

  logical(C_BOOL), bind(C) :: check_overlaps

  ! Trace for single particle
  integer(C_INT), bind(C)     :: trace_batch
  integer(C_INT), bind(C)     :: trace_gen
  integer(C_INT64_T), bind(C) :: trace_particle

  ! Particle tracks
  logical(C_BOOL), bind(C) :: write_all_tracks

  character(MAX_FILE_LEN) :: path_input               ! Path to input file
  character(MAX_FILE_LEN) :: path_cross_sections = '' ! Path to cross_sections.xml
  character(MAX_FILE_LEN) :: path_state_point         ! Path to binary state point
  character(MAX_FILE_LEN) :: path_source_point        ! Path to binary source point
  character(MAX_FILE_LEN) :: path_particle_restart    ! Path to particle restart
  character(MAX_FILE_LEN) :: path_output = ''         ! Path to output directory

  ! Various output options
  logical(C_BOOL), bind(C) :: output_summary

  ! No reduction at end of batch
  logical(C_BOOL), bind(C) :: reduce_tallies

end module settings
