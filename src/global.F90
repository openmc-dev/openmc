module global

  use, intrinsic :: ISO_C_BINDING

#ifdef MPIF08
  use mpi_f08
#endif

  use bank_header,      only: Bank
  use cmfd_header
  use constants
  use set_header,       only: SetInt
  use stl_vector,       only: VectorInt
  use source_header,    only: SourceDistribution
  use trigger_header,   only: KTrigger
  use timer_header,     only: Timer

  ! Inherit module variables from other modules
  use geometry_header
  use material_header
  use mesh_header
  use mgxs_header
  use nuclide_header
  use plot_header
  use sab_header
  use surface_header
  use tally_filter_header
  use tally_header
  use volume_header

  implicit none

  ! ============================================================================
  ! GEOMETRY-RELATED VARIABLES

  ! Number of lost particles
  integer :: n_lost_particles

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
  real(8) :: log_spacing ! spacing on logarithmic grid

  ! ============================================================================
  ! MULTI-GROUP CROSS SECTION RELATED VARIABLES

  ! Maximum Data Order
  integer :: max_order

  ! Whether or not to convert Legendres to tabulars
  logical :: legendre_to_tabular = .true.

  ! Number of points to use in the Legendre to tabular conversion
  integer :: legendre_to_tabular_points = 33

  ! ============================================================================
  ! TALLY-RELATED VARIABLES

  ! Pointers for different tallies
  type(TallyObject), pointer :: user_tallies(:) => null()
  type(TallyObject), pointer :: cmfd_tallies(:) => null()

  ! Starting index (minus 1) in tallies for each tally group
  integer :: i_user_tallies = -1
  integer :: i_cmfd_tallies = -1

  ! Active tally lists
  type(VectorInt) :: active_analog_tallies
  type(VectorInt) :: active_tracklength_tallies
  type(VectorInt) :: active_current_tallies
  type(VectorInt) :: active_collision_tallies
  type(VectorInt) :: active_tallies
  type(VectorInt) :: active_surface_tallies

  integer :: n_user_meshes  = 0 ! # of structured user meshes
  integer :: n_user_filters = 0 ! # of user filters
  integer :: n_user_tallies = 0 ! # of user tallies

  ! Normalization for statistics
  integer :: n_realizations = 0 ! # of independent realizations
  real(8) :: total_weight       ! total starting particle weight in realization

  ! Flag for turning tallies on
  logical :: tallies_on = .false.
  logical :: active_batches = .false.

  ! Assume all tallies are spatially distinct
  logical :: assume_separate = .false.

  ! Use confidence intervals for results instead of standard deviations
  logical :: confidence_intervals = .false.

  ! ============================================================================
  ! EIGENVALUE SIMULATION VARIABLES

  integer(8) :: n_particles = 0   ! # of particles per generation
  integer    :: n_batches         ! # of batches
  integer    :: n_inactive        ! # of inactive batches
  integer    :: n_active          ! # of active batches
  integer    :: gen_per_batch = 1 ! # of generations per batch
  integer    :: current_batch     ! current batch
  integer    :: current_gen       ! current generation within a batch
  integer    :: total_gen     = 0 ! total number of generations simulated

  ! ============================================================================
  ! TALLY PRECISION TRIGGER VARIABLES

  integer        :: n_max_batches             ! max # of batches
  integer        :: n_batch_interval = 1      ! batch interval for triggers
  logical        :: pred_batches = .false.    ! predict batches for triggers
  logical        :: trigger_on = .false.      ! flag for turning triggers on/off
  type(KTrigger) :: keff_trigger              ! trigger for k-effective
  logical :: satisfy_triggers = .false.       ! whether triggers are satisfied

  ! External source
  type(SourceDistribution), allocatable :: external_source(:)

  ! Source and fission bank
  type(Bank), allocatable, target :: source_bank(:)
  type(Bank), allocatable, target :: fission_bank(:)
#ifdef _OPENMP
  type(Bank), allocatable, target :: master_fission_bank(:)
#endif
  integer(8) :: n_bank       ! # of sites in fission bank
  integer(8) :: work         ! number of particles per processor
  integer(8), allocatable :: work_index(:) ! starting index in source bank for each process
  integer(8) :: current_work ! index in source bank of current history simulated

  ! Temporary k-effective values
  real(8), allocatable :: k_generation(:) ! single-generation estimates of k
  real(8) :: keff = ONE       ! average k over active batches
  real(8) :: keff_std         ! standard deviation of average k
  real(8) :: k_col_abs = ZERO ! sum over batches of k_collision * k_absorption
  real(8) :: k_col_tra = ZERO ! sum over batches of k_collision * k_tracklength
  real(8) :: k_abs_tra = ZERO ! sum over batches of k_absorption * k_tracklength

  ! Shannon entropy
  logical :: entropy_on = .false.
  real(8), allocatable :: entropy(:)         ! shannon entropy at each generation
  real(8), allocatable :: entropy_p(:,:,:,:) ! % of source sites in each cell
  type(RegularMesh), pointer :: entropy_mesh

  ! Uniform fission source weighting
  logical :: ufs = .false.
  type(RegularMesh), pointer :: ufs_mesh => null()
  real(8), allocatable :: source_frac(:,:,:,:)

  ! Write source at end of simulation
  logical :: source_separate = .false.
  logical :: source_write = .true.
  logical :: source_latest = .false.

  ! ============================================================================
  ! PARALLEL PROCESSING VARIABLES

#ifdef _OPENMP
  integer :: n_threads = NONE      ! number of OpenMP threads
  integer :: thread_id             ! ID of a given thread
#endif

  ! No reduction at end of batch
  logical :: reduce_tallies = .true.

  ! ============================================================================
  ! TIMING VARIABLES

  type(Timer) :: time_total         ! timer for total run
  type(Timer) :: time_initialize    ! timer for initialization
  type(Timer) :: time_read_xs       ! timer for reading cross sections
  type(Timer) :: time_unionize      ! timer for material xs-energy grid union
  type(Timer) :: time_bank          ! timer for fission bank synchronization
  type(Timer) :: time_bank_sample   ! timer for fission bank sampling
  type(Timer) :: time_bank_sendrecv ! timer for fission bank SEND/RECV
  type(Timer) :: time_tallies       ! timer for accumulate tallies
  type(Timer) :: time_inactive      ! timer for inactive batches
  type(Timer) :: time_active        ! timer for active batches
  type(Timer) :: time_transport     ! timer for transport only
  type(Timer) :: time_finalize      ! timer for finalization

  ! ===========================================================================
  ! VARIANCE REDUCTION VARIABLES

  logical :: survival_biasing = .false.
  real(8) :: weight_cutoff = 0.25_8
  real(8) :: energy_cutoff = ZERO
  real(8) :: weight_survive = ONE

  ! ============================================================================
  ! MISCELLANEOUS VARIABLES

  ! Mode to run in (fixed source, eigenvalue, plotting, etc)
  integer :: run_mode = NONE

  ! Restart run
  logical :: restart_run = .false.
  integer :: restart_batch

  character(MAX_FILE_LEN) :: path_input               ! Path to input file
  character(MAX_FILE_LEN) :: path_cross_sections = '' ! Path to cross_sections.xml
  character(MAX_FILE_LEN) :: path_multipole           ! Path to wmp library
  character(MAX_FILE_LEN) :: path_source = ''         ! Path to binary source
  character(MAX_FILE_LEN) :: path_state_point         ! Path to binary state point
  character(MAX_FILE_LEN) :: path_source_point        ! Path to binary source point
  character(MAX_FILE_LEN) :: path_particle_restart    ! Path to particle restart
  character(MAX_FILE_LEN) :: path_output = ''         ! Path to output directory

  ! The verbosity controls how much information will be printed to the
  ! screen and in logs
  integer :: verbosity = 7

  ! Flag for enabling cell overlap checking during transport
  logical                  :: check_overlaps = .false.
  integer(8), allocatable  :: overlap_check_cnt(:)

  ! Trace for single particle
  logical    :: trace
  integer    :: trace_batch
  integer    :: trace_gen
  integer(8) :: trace_particle

  ! Particle tracks
  logical :: write_all_tracks = .false.
  integer, allocatable :: track_identifiers(:,:)

  ! Particle restart run
  logical :: particle_restart_run = .false.

  ! Number of distribcell maps
  integer :: n_maps

  ! Write out initial source
  logical :: write_initial_source = .false.

  ! Whether create fission neutrons or not. Only applied for MODE_FIXEDSOURCE
  logical :: create_fission_neutrons = .true.

  ! ============================================================================
  ! CMFD VARIABLES

  ! Main object
  type(cmfd_type) :: cmfd

  ! Is CMFD active
  logical :: cmfd_run = .false.

  ! Timing objects
  type(Timer) :: time_cmfd      ! timer for whole cmfd calculation
  type(Timer) :: time_cmfdbuild ! timer for matrix build
  type(Timer) :: time_cmfdsolve ! timer for solver

  ! Flag for active core map
  logical :: cmfd_coremap = .false.

  ! Flag to reset dhats to zero
  logical :: dhat_reset = .false.

  ! Flag to activate neutronic feedback via source weights
  logical :: cmfd_feedback = .false.

  ! User-defined tally information
  integer :: n_cmfd_meshes  = 1 ! # of structured meshes
  integer :: n_cmfd_filters = 0 ! # of filters
  integer :: n_cmfd_tallies = 3 ! # of user-defined tallies

  ! Adjoint method type
  character(len=10) :: cmfd_adjoint_type = 'physical'

  ! Number of incomplete ilu factorization levels
  integer :: cmfd_ilu_levels = 1

  ! Batch to begin cmfd
  integer :: cmfd_begin = 1

  ! Tally reset list
  integer :: n_cmfd_resets
  type(SetInt) :: cmfd_reset

  ! Compute effective downscatter cross section
  logical :: cmfd_downscatter = .false.

  ! Convergence monitoring
  logical :: cmfd_power_monitor = .false.

  ! Cmfd output
  logical :: cmfd_write_matrices = .false.

  ! Run an adjoint calculation (last batch only)
  logical :: cmfd_run_adjoint = .false.

  ! CMFD run logicals
  logical :: cmfd_on             = .false.

  ! CMFD display info
  character(len=25) :: cmfd_display = 'balance'

  ! Estimate of spectral radius of CMFD matrices and tolerances
  real(8) :: cmfd_spectral = ZERO
  real(8) :: cmfd_shift = 1.e6
  real(8) :: cmfd_ktol = 1.e-8_8
  real(8) :: cmfd_stol = 1.e-8_8
  real(8) :: cmfd_atoli = 1.e-10_8
  real(8) :: cmfd_rtoli = 1.e-5_8

  ! Information about state points to be written
  integer :: n_state_points = 0
  type(SetInt) :: statepoint_batch

  ! Information about source points to be written
  integer :: n_source_points = 0
  type(SetInt) :: sourcepoint_batch

  ! Various output options
  logical :: output_summary = .true.
  logical :: output_tallies = .true.

  ! ============================================================================
  ! RESONANCE SCATTERING VARIABLES

  logical :: res_scat_on = .false. ! is resonance scattering treated?
  integer :: res_scat_method = RES_SCAT_ARES  ! resonance scattering method
  real(8) :: res_scat_energy_min = 0.01_8
  real(8) :: res_scat_energy_max = 1000.0_8
  character(10), allocatable :: res_scat_nuclides(:)

!$omp threadprivate(micro_xs, material_xs, fission_bank, n_bank, &
!$omp&              trace, thread_id, current_work)

contains

!===============================================================================
! FREE_MEMORY deallocates and clears  all global allocatable arrays in the
! program
!===============================================================================

  subroutine free_memory()

    integer :: i ! Loop Index

    ! Deallocate cells, surfaces, materials
    if (allocated(cells)) deallocate(cells)
    if (allocated(universes)) deallocate(universes)
    if (allocated(lattices)) deallocate(lattices)
    if (allocated(surfaces)) deallocate(surfaces)
    if (allocated(materials)) deallocate(materials)
    if (allocated(plots)) deallocate(plots)
    if (allocated(volume_calcs)) deallocate(volume_calcs)

    ! Deallocate geometry debugging information
    if (allocated(overlap_check_cnt)) deallocate(overlap_check_cnt)

    ! Deallocate cross section data, listings, and cache
    if (allocated(nuclides)) then
    ! First call the clear routines
      do i = 1, size(nuclides)
        call nuclides(i) % clear()
      end do
      deallocate(nuclides)
    end if
    if (allocated(libraries)) deallocate(libraries)

    if (allocated(res_scat_nuclides)) deallocate(res_scat_nuclides)

    if (allocated(nuclides_MG)) deallocate(nuclides_MG)

    if (allocated(macro_xs)) deallocate(macro_xs)

    if (allocated(sab_tables)) deallocate(sab_tables)

    ! Deallocate external source
    if (allocated(external_source)) deallocate(external_source)

    ! Deallocate k and entropy
    if (allocated(k_generation)) deallocate(k_generation)
    if (allocated(entropy)) deallocate(entropy)
    if (allocated(entropy_p)) deallocate(entropy_p)

    ! Deallocate tally-related arrays
    if (allocated(global_tallies)) deallocate(global_tallies)
    if (allocated(meshes)) deallocate(meshes)
    if (allocated(filters)) deallocate(filters)
    if (allocated(tallies)) deallocate(tallies)

    ! Deallocate fission and source bank and entropy
!$omp parallel
    if (allocated(fission_bank)) deallocate(fission_bank)
    if (allocated(tally_derivs)) deallocate(tally_derivs)
!$omp end parallel
#ifdef _OPENMP
    if (allocated(master_fission_bank)) deallocate(master_fission_bank)
#endif
    if (allocated(source_bank)) deallocate(source_bank)

    ! Deallocate array of work indices
    if (allocated(work_index)) deallocate(work_index)

    if (allocated(energy_bins)) deallocate(energy_bins)
    if (allocated(energy_bin_avg)) deallocate(energy_bin_avg)

    ! Deallocate CMFD
    call deallocate_cmfd(cmfd)

    ! Deallocate tally node lists
    call active_analog_tallies % clear()
    call active_tracklength_tallies % clear()
    call active_current_tallies % clear()
    call active_collision_tallies % clear()
    call active_surface_tallies % clear()
    call active_tallies % clear()

    ! Deallocate track_identifiers
    if (allocated(track_identifiers)) deallocate(track_identifiers)

    ! Deallocate dictionaries
    call cell_dict % clear()
    call universe_dict % clear()
    call lattice_dict % clear()
    call surface_dict % clear()
    call material_dict % clear()
    call mesh_dict % clear()
    call filter_dict % clear()
    call tally_dict % clear()
    call plot_dict % clear()
    call nuclide_dict % clear()
    call sab_dict % clear()
    call library_dict % clear()

    ! Clear statepoint and sourcepoint batch set
    call statepoint_batch % clear()
    call sourcepoint_batch % clear()

    ! Deallocate entropy mesh
    if (associated(entropy_mesh)) then
      if (allocated(entropy_mesh % lower_left)) &
           deallocate(entropy_mesh % lower_left)
      if (allocated(entropy_mesh % upper_right)) &
           deallocate(entropy_mesh % upper_right)
      if (allocated(entropy_mesh % width)) deallocate(entropy_mesh % width)
      deallocate(entropy_mesh)
    end if

    ! Deallocate ufs
    if (allocated(source_frac)) deallocate(source_frac)
    if (associated(ufs_mesh)) then
        if (allocated(ufs_mesh % lower_left)) deallocate(ufs_mesh % lower_left)
        if (allocated(ufs_mesh % upper_right)) &
             deallocate(ufs_mesh % upper_right)
        if (allocated(ufs_mesh % width)) deallocate(ufs_mesh % width)
        deallocate(ufs_mesh)
    end if

  end subroutine free_memory

!===============================================================================
! OVERALL_GENERATION determines the overall generation number
!===============================================================================

  pure function overall_generation() result(gen)
    integer :: gen
    gen = gen_per_batch*(current_batch - 1) + current_gen
  end function overall_generation

end module global
