module global

  use ace_header,       only: Nuclide, SAlphaBeta, xsListing, NuclideMicroXS, &
                              MaterialMacroXS
  use bank_header,      only: Bank
  use cmfd_header
  use constants
  use dict_header,      only: DictCharInt, DictIntInt
  use geometry_header,  only: Cell, Universe, Lattice, Surface
  use material_header,  only: Material
  use mesh_header,      only: StructuredMesh
  use particle_header,  only: Particle
  use plot_header,      only: ObjectPlot
  use set_header,       only: SetInt
  use source_header,    only: ExtSource
  use tally_header,     only: TallyObject, TallyMap, TallyResult
  use timer_header,     only: Timer

#ifdef MPI
  use mpi
#endif

#ifdef HDF5
  use hdf5_interface,  only: HID_T
#endif

  implicit none
  save

  ! ============================================================================
  ! THE PARTICLE

  type(Particle), pointer :: p => null()

  ! ============================================================================
  ! GEOMETRY-RELATED VARIABLES

  ! Main arrays
  type(Cell),      allocatable, target :: cells(:)
  type(Universe),  allocatable, target :: universes(:)
  type(Lattice),   allocatable, target :: lattices(:)
  type(Surface),   allocatable, target :: surfaces(:)
  type(Material),  allocatable, target :: materials(:)
  type(ObjectPlot),allocatable, target :: plots(:)

  ! Size of main arrays
  integer :: n_cells     ! # of cells
  integer :: n_universes ! # of universes
  integer :: n_lattices  ! # of lattices
  integer :: n_surfaces  ! # of surfaces
  integer :: n_materials ! # of materials
  integer :: n_plots     ! # of plots

  ! These dictionaries provide a fast lookup mechanism -- the key is the
  ! user-specified identifier and the value is the index in the corresponding
  ! array
  type(DictIntInt) :: cell_dict
  type(DictIntInt) :: universe_dict
  type(DictIntInt) :: lattice_dict
  type(DictIntInt) :: surface_dict
  type(DictIntInt) :: material_dict
  type(DictIntInt) :: mesh_dict
  type(DictIntInt) :: tally_dict
  type(DictIntInt) :: plot_dict

  ! Number of lost particles
  integer :: n_lost_particles

  ! ============================================================================
  ! CROSS SECTION RELATED VARIABLES

  ! Cross section arrays
  type(Nuclide),    allocatable, target :: nuclides(:)    ! Nuclide cross-sections
  type(SAlphaBeta), allocatable, target :: sab_tables(:)  ! S(a,b) tables
  type(XsListing),  allocatable, target :: xs_listings(:) ! cross_sections.xml listings 

  ! Cross section caches
  type(NuclideMicroXS), allocatable :: micro_xs(:)  ! Cache for each nuclide
  type(MaterialMacroXS)             :: material_xs  ! Cache for current material

  integer :: n_nuclides_total ! Number of nuclide cross section tables
  integer :: n_sab_tables     ! Number of S(a,b) thermal scattering tables
  integer :: n_listings       ! Number of listings in cross_sections.xml

  ! Dictionaries to look up cross sections and listings
  type(DictCharInt) :: nuclide_dict
  type(DictCharInt) :: sab_dict
  type(DictCharInt) :: xs_listing_dict

  ! Unionized energy grid
  integer :: grid_method ! how to treat the energy grid
  integer :: n_grid      ! number of points on unionized grid
  real(8), allocatable :: e_grid(:) ! energies on unionized grid

  ! Unreoslved resonance probablity tables
  logical :: urr_ptables_on = .true.

  ! Default xs identifier (e.g. 70c)
  character(3):: default_xs

  ! ============================================================================
  ! TALLY-RELATED VARIABLES

  type(StructuredMesh), allocatable, target :: meshes(:)
  type(TallyObject),    allocatable, target :: tallies(:)

  ! Pointers for different tallies
  type(TallyObject), pointer :: user_tallies(:) => null()
  type(TallyObject), pointer :: cmfd_tallies(:) => null()

  ! Starting index (minus 1) in tallies for each tally group
  integer :: i_user_tallies = -1
  integer :: i_cmfd_tallies = -1

  ! Active tally lists
  type(SetInt) :: active_analog_tallies
  type(SetInt) :: active_tracklength_tallies
  type(SetInt) :: active_current_tallies
  type(SetInt) :: active_tallies

  ! Global tallies
  !   1) collision estimate of k-eff
  !   2) track-length estimate of k-eff
  !   3) leakage fraction

  type(TallyResult), target :: global_tallies(N_GLOBAL_TALLIES)

  ! Tally map structure
  type(TallyMap), allocatable :: tally_maps(:)

  integer :: n_meshes       = 0 ! # of structured meshes
  integer :: n_user_meshes  = 0 ! # of structured user meshes
  integer :: n_tallies      = 0 ! # of tallies
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
  integer    :: current_batch = 0 ! current batch
  integer    :: current_gen   = 0 ! current generation within a batch
  integer    :: overall_gen   = 0 ! overall generation in the run

  ! External source
  type(ExtSource), target :: external_source

  ! Source and fission bank
  type(Bank), allocatable, target :: source_bank(:)
  type(Bank), allocatable, target :: fission_bank(:)
  integer(8) :: n_bank       ! # of sites in fission bank
  integer(8) :: bank_first   ! index of first particle in bank
  integer(8) :: bank_last    ! index of last particle in bank
  integer(8) :: work         ! number of particles per processor
  integer(8) :: maxwork      ! maximum number of particles per processor
  integer(8) :: current_work ! index in source bank of current history simulated

  ! Temporary k-effective values
  real(8), allocatable :: k_generation(:) ! single-generation estimates of k
  real(8) :: keff = ONE       ! average k over active batches
  real(8) :: keff_std         ! standard deviation of average k
  real(8) :: k_col_abs = ZERO ! sum over batches of k_collision * k_absorption
  real(8) :: k_col_tra = ZERO ! sum over batches of k_collision * k_tracklength
  real(8) :: k_abs_tra = ZERO ! sum over batches of k_absorption * k_tracklength
  real(8) :: k_combined(2)    ! combined best estimate of k-effective

  ! Shannon entropy
  logical :: entropy_on = .false.
  real(8), allocatable :: entropy(:)         ! shannon entropy at each generation
  real(8), allocatable :: entropy_p(:,:,:,:) ! % of source sites in each cell
  type(StructuredMesh), pointer :: entropy_mesh

  ! Uniform fission source weighting
  logical :: ufs = .false.
  type(StructuredMesh), pointer :: ufs_mesh => null()
  real(8), allocatable :: source_frac(:,:,:,:)

  ! Write source at end of simulation
  logical :: source_separate = .false.
  logical :: source_write = .true.

  ! ============================================================================
  ! PARALLEL PROCESSING VARIABLES

  ! The defaults set here for the number of processors, rank, and master and
  ! mpi_enabled flag are for when MPI is not being used at all, i.e. a serial
  ! run. In this case, these variables are still used at times.

  integer :: n_procs     = 1       ! number of processes
  integer :: rank        = 0       ! rank of process
  logical :: master      = .true.  ! master process?
  logical :: mpi_enabled = .false. ! is MPI in use and initialized?
  integer :: mpi_err               ! MPI error code
  integer :: MPI_BANK              ! MPI datatype for fission bank
  integer :: MPI_TALLYRESULT       ! MPI datatype for TallyResult

  ! No reduction at end of batch
  logical :: reduce_tallies = .true.

  ! ============================================================================
  ! TIMING VARIABLES

  type(Timer) :: time_total         ! timer for total run
  type(Timer) :: time_initialize    ! timer for initialization
  type(Timer) :: time_read_xs       ! timer for reading cross sections
  type(Timer) :: time_unionize      ! timer for unionizing energy grid
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
  real(8) :: weight_cutoff = 0.25
  real(8) :: weight_survive = 1.0

  ! ============================================================================
  ! HDF5 VARIABLES

#ifdef HDF5
  integer(HID_T) :: hdf5_output_file   ! identifier for output file
  integer(HID_T) :: hdf5_tallyresult_t ! Compound type for TallyResult
  integer(HID_T) :: hdf5_bank_t        ! Compound type for Bank
  integer(HID_T) :: hdf5_integer8_t    ! type for integer(8)
#endif

  ! ============================================================================
  ! MISCELLANEOUS VARIABLES

  ! Mode to run in (fixed source, eigenvalue, plotting, etc)
  integer :: run_mode = NONE

  ! Restart run
  logical :: restart_run = .false.
  integer :: restart_batch

  character(MAX_FILE_LEN) :: path_input            ! Path to input file
  character(MAX_FILE_LEN) :: path_cross_sections   ! Path to cross_sections.xml
  character(MAX_FILE_LEN) :: path_source = ''      ! Path to binary source
  character(MAX_FILE_LEN) :: path_state_point      ! Path to binary state point
  character(MAX_FILE_LEN) :: path_particle_restart ! Path to particle restart
  character(MAX_FILE_LEN) :: path_output = ''      ! Path to output directory

  ! Message used in message/warning/fatal_error
  character(MAX_LINE_LEN) :: message

  ! Random number seed
  integer(8) :: seed = 1_8

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

  ! Particle restart run
  logical :: particle_restart_run = .false.

  ! ============================================================================
  ! CMFD VARIABLES 

  ! Main object
  type(cmfd_type) :: cmfd

  ! Is CMFD active
  logical :: cmfd_run = .false.
 
  ! Timing objects
  type(Timer) :: time_cmfd   ! timer for whole cmfd calculation
  type(Timer) :: time_solver ! timer for solver 

  ! Flag for CMFD only
  logical :: cmfd_only = .false.

  ! Flag for coremap accelerator
  logical :: cmfd_coremap = .false.

  ! number of processors for cmfd
  integer :: n_procs_cmfd

  ! reset dhats to zero
  logical :: dhat_reset = .false.

  ! activate neutronic feedback
  logical :: cmfd_feedback = .false.

  ! activate auto-balance of tallies (2grp only)
! logical :: cmfd_balance = .false.

  ! calculate effective downscatter
! logical :: cmfd_downscatter = .false.

  ! user-defined tally information
  integer :: n_cmfd_meshes              = 1 ! # of structured meshes
  integer :: n_cmfd_tallies             = 3 ! # of user-defined tallies

  ! overwrite with 2grp xs
  logical :: cmfd_run_2grp = .false.

  ! hold cmfd weight adjustment factors
  logical :: cmfd_hold_weights = .false.

  ! eigenvalue solver type
  character(len=10) :: cmfd_solver_type = 'power'

  ! adjoint method type
  character(len=10) :: cmfd_adjoint_type = 'physical'

  ! number of incomplete ilu factorization levels
  integer :: cmfd_ilu_levels = 1

  ! batch to begin cmfd
  integer :: cmfd_begin = 1

  ! when and how long to flush cmfd tallies during inactive batches
  integer :: cmfd_inact_flush(2) = (/9999,1/)

  ! batch to last flush before active batches
  integer :: cmfd_act_flush = 0

  ! compute effective downscatter cross section
  logical :: cmfd_downscatter = .false.

  ! convergence monitoring
  logical :: cmfd_snes_monitor  = .false.
  logical :: cmfd_ksp_monitor   = .false.
  logical :: cmfd_power_monitor = .false.

  ! cmfd output
  logical :: cmfd_write_balance  = .false.
  logical :: cmfd_write_matrices = .false.
  logical :: cmfd_write_hdf5     = .false.

  ! run an adjoint calculation (last batch only)
  logical :: cmfd_run_adjoint = .false.

  ! cmfd run logicals
  logical :: cmfd_on             = .false.
  logical :: cmfd_tally_on       = .true. 

  ! tolerance on keff to run cmfd
  real(8) :: cmfd_keff_tol = 0.005_8

  ! Information about state points to be written
  integer :: n_state_points = 0
  type(SetInt) :: statepoint_batch

  ! Various output options
  logical :: output_summary = .false.
  logical :: output_xs      = .false.
  logical :: output_tallies = .true.

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
    if (allocated(sab_tables)) deallocate(sab_tables)
    if (allocated(xs_listings)) deallocate(xs_listings)
    if (allocated(micro_xs)) deallocate(micro_xs)

    ! Deallocate tally-related arrays
    if (allocated(meshes)) deallocate(meshes)
    if (allocated(tallies)) then
    ! First call the clear routines
      do i = 1, size(tallies)
        call tallies(i) % clear()
      end do
      ! Now deallocate the tally array
      deallocate(tallies)
    end if
    if (allocated(tally_maps)) deallocate(tally_maps)

    ! Deallocate energy grid
    if (allocated(e_grid)) deallocate(e_grid)

    ! Deallocate fission and source bank and entropy
    if (allocated(fission_bank)) deallocate(fission_bank)
    if (allocated(source_bank)) deallocate(source_bank)
    if (allocated(entropy_p)) deallocate(entropy_p)

    ! Deallocate cmfd
    call deallocate_cmfd(cmfd)

    ! Deallocate tally node lists
    call active_analog_tallies % clear()
    call active_tracklength_tallies % clear()
    call active_current_tallies % clear()
    call active_tallies % clear()
    
    ! Deallocate dictionaries
    call cell_dict % clear()
    call universe_dict % clear()
    call lattice_dict % clear()
    call surface_dict % clear()
    call material_dict % clear()
    call mesh_dict % clear()
    call tally_dict % clear()
    call plot_dict % clear()
    call nuclide_dict % clear()
    call sab_dict % clear()
    call xs_listing_dict % clear()
    
  end subroutine free_memory

end module global
