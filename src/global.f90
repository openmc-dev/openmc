module global

  use bank_header,          only: Bank
  use constants
  use cross_section_header, only: Nuclide, SAB_Table, xsData
  use datatypes_header,     only: DictionaryII, DictionaryCI
  use geometry_header,      only: Cell, Universe, Lattice, Surface
  use material_header,      only: Material
  use particle_header,      only: Particle
  use source_header,        only: ExtSource
  use tally_header,         only: Tally
  use timing,               only: Timer

#ifdef MPI
  use mpi
#endif

  implicit none

  ! Main arrays for cells, surfaces, materials
  type(Cell),     allocatable, target :: cells(:)
  type(Universe), allocatable, target :: universes(:)
  type(Lattice),  allocatable, target :: lattices(:)
  type(Surface),  allocatable, target :: surfaces(:)
  type(Material), allocatable, target :: materials(:)
  type(xsData),   allocatable, target :: xsdatas(:)
  type(Tally),    allocatable, target :: tallies(:)
  type(Tally),    allocatable, target :: tallies_global(:)
  integer :: n_cells          ! # of cells
  integer :: n_universes      ! # of universes
  integer :: n_lattices       ! # of lattices
  integer :: n_surfaces       ! # of surfaces
  integer :: n_materials      ! # of materials
  integer :: n_tallies        ! # of tallies
  integer :: n_tallies_global ! # of global tallies

  ! These dictionaries provide a fast lookup mechanism
  type(DictionaryII), pointer :: cell_dict
  type(DictionaryII), pointer :: universe_dict
  type(DictionaryII), pointer :: lattice_dict
  type(DictionaryII), pointer :: surface_dict
  type(DictionaryII), pointer :: material_dict
  type(DictionaryII), pointer :: tally_dict
  type(DictionaryCI), pointer :: xsdata_dict
  type(DictionaryCI), pointer :: nuclide_dict

  ! Cross section arrays
  type(Nuclide),   allocatable, target :: nuclides(:)
  type(SAB_Table), allocatable, target :: sab_tables(:)
  integer :: n_nuclides
  integer :: n_sab_tables

  ! Current cell, surface, material
  type(Cell),     pointer :: cCell
  type(Universe), pointer :: cUniverse
  type(Lattice),  pointer :: cLattice
  type(Surface),  pointer :: cSurface
  type(Material), pointer :: cMaterial

  ! unionized energy grid
  integer              :: n_grid    ! number of points on unionized grid
  real(8), allocatable :: e_grid(:) ! energies on unionized grid

  ! Histories/cycles/etc for both external source and criticality
  integer(8) :: n_particles ! # of particles (per cycle for criticality)
  integer    :: n_cycles    ! # of cycles
  integer    :: n_inactive  ! # of inactive cycles

  ! External source
  type(ExtSource), target :: external_source

  ! Source and fission bank
  type(Particle), allocatable, target :: source_bank(:)
  type(Bank),     allocatable, target :: fission_bank(:)
  integer(8) :: n_bank      ! # of sites in fission bank
  integer(8) :: bank_first  ! index of first particle in bank
  integer(8) :: bank_last   ! index of last particle in bank
  integer(8) :: work        ! number of particles per processor

  ! cycle keff
  real(8) :: keff

  logical :: tallies_on

  ! Parallel processing variables
  integer :: n_procs     ! number of processes
  integer :: rank        ! rank of process
  logical :: master      ! master process?
  logical :: mpi_enabled ! is MPI in use and initialized?

  ! Timing variables
  type(Timer) :: time_total       ! timer for total run
  type(Timer) :: time_init        ! timer for initialization
  type(Timer) :: time_intercycle  ! timer for intercycle synchronization
  type(Timer) :: time_compute     ! timer for computation

  ! Paths to input file, cross section data, etc
  character(MAX_WORD_LEN) :: & 
       & path_input,         &
       & path_xsdata

  ! Problem type
  integer :: problem_type

  ! The verbosity controls how much information will be printed to the
  ! screen and in logs
  integer :: verbosity

contains

!===============================================================================
! SET_DEFAULTS gives default values for many global parameters
!===============================================================================

  subroutine set_defaults()

    ! Default problem type is external source
    problem_type = PROB_SOURCE
    
    ! Default number of particles
    n_particles = 10000

    ! Default verbosity
    verbosity = 5

    ! Defualt multiplication factor
    keff = ONE

  end subroutine set_defaults

!===============================================================================
! FREE_MEMORY deallocates all allocatable arrays in the program, namely the
! cells, surfaces, materials, and sources
!===============================================================================

  subroutine free_memory()

    integer :: ierr

    ! Deallocate cells, surfaces, materials
    if (allocated(cells)) deallocate(cells)
    if (allocated(surfaces)) deallocate(surfaces)
    if (allocated(materials)) deallocate(materials)
    if (allocated(lattices)) deallocate(lattices)

    ! Deallocate cross section data
    if (allocated(xsdatas)) deallocate(xsdatas)
    if (allocated(nuclides)) deallocate(nuclides)
    if (allocated(sab_tables)) deallocate(sab_tables)

    ! Deallocate energy grid
    if (allocated(e_grid)) deallocate(e_grid)

    ! Deallocate fission and source bank
    if (allocated(fission_bank)) deallocate(fission_bank)
    if (allocated(source_bank)) deallocate(source_bank)

#ifdef MPI
    ! If MPI is in use and enabled, terminate it
    call MPI_FINALIZE(ierr)
#endif

    ! End program
    stop
   
  end subroutine free_memory

end module global
