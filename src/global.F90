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
  use source_header,    only: ExtSource
  use tally_header,     only: TallyObject, TallyMap
  use timing,           only: Timer

#ifdef MPI
  use mpi
#endif

  implicit none
  save

  ! ============================================================================
  ! GEOMETRY-RELATED VARIABLES

  ! Main arrays
  type(Cell),     allocatable, target :: cells(:)
  type(Universe), allocatable, target :: universes(:)
  type(Lattice),  allocatable, target :: lattices(:)
  type(Surface),  allocatable, target :: surfaces(:)
  type(Material), allocatable, target :: materials(:)

  ! Size of main arrays
  integer :: n_cells     ! # of cells
  integer :: n_universes ! # of universes
  integer :: n_lattices  ! # of lattices
  integer :: n_surfaces  ! # of surfaces
  integer :: n_materials ! # of materials

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
  integer              :: n_grid    ! number of points on unionized grid
  real(8), allocatable :: e_grid(:) ! energies on unionized grid

  ! Unreoslved resonance probablity tables
  logical :: urr_ptables_on = .false.

  ! ============================================================================
  ! TALLY-RELATED VARIABLES

  type(StructuredMesh), allocatable, target :: meshes(:)
  type(TallyObject),    allocatable, target :: tallies(:)

  ! Tally map structure
  type(TallyMap), allocatable :: tally_maps(:)

  integer :: n_meshes         ! # of structured meshes
  integer :: n_tallies        ! # of tallies

  ! Flag for turning tallies on
  logical :: tallies_on

  ! ============================================================================
  ! CRITICALITY SIMULATION VARIABLES

  integer(8) :: n_particles = 10000 ! # of particles per cycle
  integer    :: n_cycles    = 500   ! # of cycles
  integer    :: n_inactive  = 50    ! # of inactive cycles

  ! External source
  type(ExtSource), target :: external_source

  ! Source and fission bank
  type(Particle), allocatable, target :: source_bank(:)
  type(Bank),     allocatable, target :: fission_bank(:)
  integer(8) :: n_bank       ! # of sites in fission bank
  integer(8) :: bank_first   ! index of first particle in bank
  integer(8) :: bank_last    ! index of last particle in bank
  integer(8) :: work         ! number of particles per processor
  integer(8) :: source_index ! index for source particles

  ! cycle keff
  real(8) :: keff = ONE
  real(8) :: keff_std

  ! Shannon entropy
  logical :: entropy_on = .false.
  real(8) :: entropy                ! value of shannon entropy
  real(8) :: entropy_lower_left(3)  ! lower-left corner for entropy box
  real(8) :: entropy_upper_right(3) ! upper-right corner for entropy box
  real(8), allocatable :: entropy_p(:,:,:)

  ! ============================================================================
  ! PARALLEL PROCESSING VARIABLES

  integer :: n_procs     ! number of processes
  integer :: rank        ! rank of process
  logical :: master      ! master process?
  logical :: mpi_enabled ! is MPI in use and initialized?
  integer :: mpi_err     ! MPI error code
  integer :: MPI_BANK    ! MPI datatype for fission bank

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
  type(Timer) :: time_ic_rebuild  ! timer for intercycle source bank rebuild
  type(Timer) :: time_inactive    ! timer for inactive cycles
  type(Timer) :: time_active      ! timer for active cycles
  type(Timer) :: time_compute     ! timer for computation

  ! ===========================================================================
  ! VARIANCE REDUCTION VARIABLES

  logical :: survival_biasing = .false.
  real(8) :: weight_cutoff = 0.25
  real(8) :: weight_survive = 1.0

  ! ============================================================================
  ! PLOTTING VARIABLES

  logical :: plotting = .false.
  real(8) :: plot_origin(3)
  real(8) :: plot_width(2)
  real(8) :: plot_basis(6)
  real(8) :: pixel

  ! ============================================================================
  ! MISCELLANEOUS VARIABLES

  character(MAX_FILE_LEN) :: path_input          ! Path to input file
  character(MAX_FILE_LEN) :: path_cross_sections ! Path to cross_sections.xml

  ! Message used in message/warning/fatal_error
  character(MAX_LINE_LEN) :: message

  ! Problem type
  integer :: problem_type = PROB_CRITICALITY

  ! The verbosity controls how much information will be printed to the
  ! screen and in logs
  integer :: verbosity = 7

  ! Trace for single particle
  logical    :: trace
  integer    :: trace_cycle
  integer(8) :: trace_particle

contains

!===============================================================================
! FREE_MEMORY deallocates all allocatable arrays in the program, namely the
! cells, surfaces, materials, and sources
!===============================================================================

  subroutine free_memory()

    ! Deallocate cells, surfaces, materials
    if (allocated(cells)) deallocate(cells)
    if (allocated(surfaces)) deallocate(surfaces)
    if (allocated(materials)) deallocate(materials)
    if (allocated(lattices)) deallocate(lattices)

    ! Deallocate cross section data and listings
    if (allocated(nuclides)) deallocate(nuclides)
    if (allocated(sab_tables)) deallocate(sab_tables)
    if (allocated(xs_listings)) deallocate(xs_listings)

    ! Deallocate energy grid
    if (allocated(e_grid)) deallocate(e_grid)

    ! Deallocate fission and source bank
    if (allocated(fission_bank)) deallocate(fission_bank)
    if (allocated(source_bank)) deallocate(source_bank)

  end subroutine free_memory

end module global
