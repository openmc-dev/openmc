module global

  use ace_header,       only: Nuclide, SAB_Table, xsListing, NuclideMicroXS, &
                              MaterialMacroXS
  use bank_header,      only: Bank
  use constants
  use datatypes_header, only: DictionaryII, DictionaryCI
  use geometry_header,  only: Cell, Universe, Lattice, Surface
  use material_header,  only: Material
  use mesh_header,      only: StructuredMesh
  use particle_header,  only: Particle
  use plot_header,      only: Plot
  use source_header,    only: ExtSource
  use tally_header,     only: TallyObject, TallyMap, TallyScore
  use timing,           only: Timer

#ifdef MPI
  use mpi
#endif

#ifdef HDF5
  use hdf5
#endif

  implicit none
  save

  ! ============================================================================
  ! THE PARTICLE

  type(Particle), pointer :: p => null()

  ! ============================================================================
  ! GEOMETRY-RELATED VARIABLES

  ! Main arrays
  type(Cell),     allocatable, target :: cells(:)
  type(Universe), allocatable, target :: universes(:)
  type(Lattice),  allocatable, target :: lattices(:)
  type(Surface),  allocatable, target :: surfaces(:)
  type(Material), allocatable, target :: materials(:)
  type(Plot),     allocatable, target :: plots(:)

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
  type(DictionaryII), pointer :: cell_dict     => null()
  type(DictionaryII), pointer :: universe_dict => null()
  type(DictionaryII), pointer :: lattice_dict  => null()
  type(DictionaryII), pointer :: surface_dict  => null()
  type(DictionaryII), pointer :: material_dict => null()
  type(DictionaryII), pointer :: mesh_dict     => null()
  type(DictionaryII), pointer :: tally_dict    => null()

  ! ============================================================================
  ! CROSS SECTION RELATED VARIABLES

  ! Cross section arrays
  type(Nuclide),   allocatable, target :: nuclides(:)    ! Nuclide cross-sections
  type(SAB_Table), allocatable, target :: sab_tables(:)  ! S(a,b) tables
  type(XsListing), allocatable, target :: xs_listings(:) ! cross_sections.xml listings 

  ! Cross section caches
  type(NuclideMicroXS), allocatable :: micro_xs(:)  ! Cache for each nuclide
  type(MaterialMacroXS)             :: material_xs  ! Cache for current material

  integer :: n_nuclides_total ! Number of nuclide cross section tables
  integer :: n_sab_tables     ! Number of S(a,b) thermal scattering tables
  integer :: n_listings       ! Number of listings in cross_sections.xml

  ! Dictionaries to look up cross sections and listings
  type(DictionaryCI), pointer :: nuclide_dict    => null()
  type(DictionaryCI), pointer :: sab_dict        => null()
  type(DictionaryCI), pointer :: xs_listing_dict => null()

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

  ! Pointers for analog, track-length, and surface-current tallies
  integer, allocatable :: analog_tallies(:)
  integer, allocatable :: tracklength_tallies(:)
  integer, allocatable :: current_tallies(:)

  ! Global tallies
  !   1) analog estimate of k-eff
  !   2) collision estimate of k-eff
  !   3) track-length estimate of k-eff
  !   4) leakage fraction

  type(TallyScore) :: global_tallies(N_GLOBAL_TALLIES)

  ! Tally map structure
  type(TallyMap), allocatable :: tally_maps(:)

  integer :: n_meshes                  ! # of structured meshes
  integer :: n_tallies                 ! # of tallies
  integer :: n_analog_tallies      = 0 ! # of analog tallies
  integer :: n_tracklength_tallies = 0 ! # of track-length tallies
  integer :: n_current_tallies     = 0 ! # of surface current tallies

  ! Normalization for statistics
  integer :: n_realizations = 0 ! # of independent realizations
  real(8) :: total_weight       ! total starting particle weight in realization

  ! Flag for turning tallies on
  logical :: tallies_on

  ! Assume all tallies are spatially distinct
  logical :: assume_separate = .false.

  ! ============================================================================
  ! CRITICALITY SIMULATION VARIABLES

  integer(8) :: n_particles = 0   ! # of particles per generation
  integer    :: n_batches         ! # of batches
  integer    :: n_inactive        ! # of inactive batches
  integer    :: n_active          ! # of active batches
  integer    :: gen_per_batch = 1 ! # of generations per batch
  integer    :: current_batch = 0 ! current batch
  integer    :: current_gen   = 0 ! current generation

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

  ! Temporary k-effective values
  real(8), allocatable :: k_batch(:) ! batch estimates of k
  real(8) :: keff = ONE ! average k over active cycles
  real(8) :: keff_std   ! standard deviation of average k

  ! Shannon entropy
  logical :: entropy_on = .false.
  real(8) :: entropy                         ! value of shannon entropy
  real(8), allocatable :: entropy_p(:,:,:,:) ! % of source sites in each cell
  type(StructuredMesh), pointer :: entropy_mesh

  ! Uniform fission source weighting
  logical :: ufs = .false.
  type(StructuredMesh), pointer :: ufs_mesh => null()
  real(8), allocatable :: source_frac(:,:,:,:)

  ! Write source at end of simulation
  logical :: write_source = .false.

  ! ============================================================================
  ! PARALLEL PROCESSING VARIABLES

  integer :: n_procs     ! number of processes
  integer :: rank        ! rank of process
  logical :: master      ! master process?
  logical :: mpi_enabled ! is MPI in use and initialized?
  integer :: mpi_err     ! MPI error code
  integer :: MPI_BANK    ! MPI datatype for fission bank

  ! No reduction at end of batch
  logical :: reduce_tallies = .true.

  ! ============================================================================
  ! TIMING VARIABLES

  type(Timer) :: time_total       ! timer for total run
  type(Timer) :: time_initialize  ! timer for initialization
  type(Timer) :: time_read_xs     ! timer for reading cross sections
  type(Timer) :: time_unionize    ! timer for unionizing energy grid
  type(Timer) :: time_intercycle  ! timer for intercycle synchronization
  type(Timer) :: time_ic_tallies  ! timer for intercycle accumulate tallies
  type(Timer) :: time_ic_sample   ! timer for intercycle sampling
  type(Timer) :: time_ic_sendrecv ! timer for intercycle SEND/RECV
  type(Timer) :: time_inactive    ! timer for inactive cycles
  type(Timer) :: time_active      ! timer for active cycles
  type(Timer) :: time_transport   ! timer for transport only
  type(Timer) :: time_finalize    ! timer for finalization

  ! ===========================================================================
  ! VARIANCE REDUCTION VARIABLES

  logical :: survival_biasing = .false.
  real(8) :: weight_cutoff = 0.25
  real(8) :: weight_survive = 1.0

  ! ============================================================================
  ! HDF5 VARIABLES

#ifdef HDF5
  integer(HID_T) :: hdf5_output_file ! identifier for output file
  integer        :: hdf5_err         ! error flag 
#endif

  ! ============================================================================
  ! MISCELLANEOUS VARIABLES

  ! Mode to run in (fixed source, criticality, plotting, etc)
  integer :: run_mode = MODE_CRITICALITY

  ! Restart run
  logical :: restart_run = .false.
  integer :: restart_batch

  character(MAX_FILE_LEN) :: path_input          ! Path to input file
  character(MAX_FILE_LEN) :: path_cross_sections ! Path to cross_sections.xml
  character(MAX_FILE_LEN) :: path_source = ''    ! Path to binary source
  character(MAX_FILE_LEN) :: path_state_point    ! Path to binary state point

  ! Message used in message/warning/fatal_error
  character(MAX_LINE_LEN) :: message

  ! Random number seed
  integer(8) :: seed = 1_8

  ! The verbosity controls how much information will be printed to the
  ! screen and in logs
  integer :: verbosity = 7

  ! Trace for single particle
  logical    :: trace
  integer    :: trace_batch
  integer    :: trace_gen
  integer(8) :: trace_particle

  ! Information about state points to be written
  integer :: n_state_points = 0
  integer, allocatable :: statepoint_batch(:)

contains

!===============================================================================
! FREE_MEMORY deallocates all global allocatable arrays in the program
!===============================================================================

  subroutine free_memory()

    ! Deallocate cells, surfaces, materials
    if (allocated(cells)) deallocate(cells)
    if (allocated(universes)) deallocate(universes)
    if (allocated(lattices)) deallocate(lattices)
    if (allocated(surfaces)) deallocate(surfaces)
    if (allocated(materials)) deallocate(materials)
    if (allocated(plots)) deallocate(plots)

    ! Deallocate cross section data, listings, and cache
    if (allocated(nuclides)) deallocate(nuclides)
    if (allocated(sab_tables)) deallocate(sab_tables)
    if (allocated(xs_listings)) deallocate(xs_listings)
    if (allocated(micro_xs)) deallocate(micro_xs)

    ! Deallocate tally-related arrays
    if (allocated(meshes)) deallocate(meshes)
    if (allocated(tallies)) deallocate(tallies)
    if (allocated(analog_tallies)) deallocate(analog_tallies)
    if (allocated(tracklength_tallies)) deallocate(tracklength_tallies)
    if (allocated(current_tallies)) deallocate(current_tallies)
    if (allocated(tally_maps)) deallocate(tally_maps)

    ! Deallocate energy grid
    if (allocated(e_grid)) deallocate(e_grid)

    ! Deallocate fission and source bank and entropy
    if (allocated(fission_bank)) deallocate(fission_bank)
    if (allocated(source_bank)) deallocate(source_bank)
    if (allocated(entropy_p)) deallocate(entropy_p)

  end subroutine free_memory

end module global
