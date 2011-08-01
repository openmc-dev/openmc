module global

  use types

#ifdef MPI
  use mpi
#endif

  implicit none

  ! Maximum number of words in a single line, length of line, and length of
  ! single word
  integer, parameter :: max_words    = 500
  integer, parameter :: max_line_len = 250
  integer, parameter :: max_word_len = 150

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
  type(DictionaryII), pointer :: isotope_dict
  type(DictionaryII), pointer :: tally_dict
  type(DictionaryCI), pointer :: xsdata_dict
  type(DictionaryCI), pointer :: ace_dict

  ! Cross section arrays
  type(AceContinuous), allocatable, target :: xs_continuous(:)
  type(AceThermal),    allocatable, target :: xs_thermal(:)
  integer :: n_continuous
  integer :: n_thermal

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
  type(TimerObj) :: time_total       ! timer for total run
  type(TimerObj) :: time_init        ! timer for initialization
  type(TimerObj) :: time_intercycle  ! timer for intercycle synchronization
  type(TimerObj) :: time_compute     ! timer for computation

  ! Paths to input file, cross section data, etc
  character(max_word_len) :: & 
       & path_input,         &
       & path_xsdata
  integer, parameter :: UNIT_LOG = 9 ! unit # for writing log file

  ! Physical constants
  real(8), parameter ::            &
       & PI           = 3.1415926535898_8, & ! pi
       & MASS_NEUTRON = 1.0086649156,      & ! mass of a neutron
       & MASS_PROTON  = 1.00727646677,     & ! mass of a proton
       & AMU          = 1.66053873e-27,    & ! 1 amu in kg
       & N_AVOGADRO   = 0.602214179,       & ! Avogadro's number in 10^24/mol
       & K_BOLTZMANN  = 8.617342e-5,       & ! Boltzmann constant in eV/K
       & INFINITY     = huge(0.0_8),       & ! positive infinity
       & ZERO         = 0.0_8,             &
       & ONE          = 1.0_8,             &
       & TWO          = 2.0_8

  ! Boundary conditions
  integer, parameter ::    &
       & BC_TRANSMIT = 0,  & ! Transmission boundary condition (default)
       & BC_VACUUM   = 1,  & ! Vacuum boundary condition
       & BC_REFLECT  = 2,  & ! Reflecting boundary condition
       & BC_PERIODIC = 3     ! Periodic boundary condition

  ! Logical operators for cell definitions
  integer, parameter ::                &
       & OP_LEFT_PAREN  = huge(0),     & ! Left parentheses
       & OP_RIGHT_PAREN = huge(0) - 1, & ! Right parentheses
       & OP_UNION       = huge(0) - 2, & ! Union operator
       & OP_DIFFERENCE  = huge(0) - 3    ! Difference operator

  ! Cell types
  integer, parameter ::    &
       & CELL_NORMAL  = 1, & ! Cell with a specified material
       & CELL_FILL    = 2, & ! Cell filled by a separate universe
       & CELL_LATTICE = 3    ! Cell filled with a lattice

  ! Lattice types
  integer, parameter :: &
       & LATTICE_RECT = 1, &
       & LATTICE_HEX  = 2

  ! array index of universe 0
  integer :: BASE_UNIVERSE

  ! Surface types
  integer, parameter ::    &
       & SURF_PX     =  1, & ! Plane parallel to x-plane 
       & SURF_PY     =  2, & ! Plane parallel to y-plane 
       & SURF_PZ     =  3, & ! Plane parallel to z-plane 
       & SURF_PLANE  =  4, & ! Arbitrary plane
       & SURF_CYL_X  =  5, & ! Cylinder along x-axis
       & SURF_CYL_Y  =  6, & ! Cylinder along y-axis
       & SURF_CYL_Z  =  7, & ! Cylinder along z-axis
       & SURF_SPHERE =  8, & ! Sphere
       & SURF_BOX_X  =  9, & ! Box extending infinitely in x-direction
       & SURF_BOX_Y  = 10, & ! Box extending infinitely in y-direction
       & SURF_BOX_Z  = 11, & ! Box extending infinitely in z-direction
       & SURF_BOX    = 12, & ! Rectangular prism
       & SURF_GQ     = 13    ! General quadratic surface

  ! Surface senses
  integer, parameter ::      &
       & SENSE_POSITIVE = 1, &
       & SENSE_NEGATIVE = -1

  ! Source types
  integer, parameter ::   &
       & SRC_BOX     = 1, & ! Source in a rectangular prism
       & SRC_CELL    = 2, & ! Source in a cell
       & SRC_SURFACE = 3    ! Source on a surface

  ! Problem type
  integer :: problem_type
  integer, parameter ::        &
       & PROB_SOURCE      = 1, & ! External source problem
       & PROB_CRITICALITY = 2    ! Criticality problem

  ! Angular distribution type
  integer, parameter :: & 
       & ANGLE_ISOTROPIC = 1, & ! Isotropic angular distribution
       & ANGLE_32_EQUI   = 2, & ! 32 equiprobable bins
       & ANGLE_TABULAR   = 3    ! Tabular angular distribution

  ! Interpolation flag
  integer, parameter ::     &
       & HISTOGRAM     = 1, & ! y is constant in x
       & LINEAR_LINEAR = 2, & ! y is linear in x
       & LINEAR_LOG    = 3, & ! y is linear in ln(x)
       & LOG_LINEAR    = 4, & ! ln(y) is linear in x
       & LOG_LOG       = 5    ! ln(y) is linear in ln(x)

  ! Particle type
  integer, parameter :: &
       & NEUTRON  = 1, &
       & PHOTON   = 2, &
       & ELECTRON = 3

  ! Reaction types
  integer, parameter :: &
       & TOTAL_XS  = 1, &
       & ELASTIC   = 2, &
       & N_LEVEL   = 4, &
       & MISC      = 5, &
       & N_2ND     = 11, &
       & N_2N      = 16, &
       & N_3N      = 17, &
       & FISSION   = 18, &
       & N_F       = 19, &
       & N_NF      = 20, &
       & N_2NF     = 21, &
       & N_NA      = 22, &
       & N_N3A     = 23, &
       & N_2NA     = 24, &
       & N_3NA     = 25, &
       & N_NP      = 28, &
       & N_N2A     = 29, &
       & N_2N2A    = 30, &
       & N_ND      = 32, &
       & N_NT      = 33, &
       & N_N3HE    = 34, &
       & N_ND2A    = 35, &
       & N_NT2A    = 36, &
       & N_4N      = 37, &
       & N_3NF     = 38, &
       & N_2NP     = 41, &
       & N_3NP     = 42, &
       & N_N2P     = 44, &
       & N_NPA     = 45, &
       & N_N1      = 51, &
       & N_N40     = 90, &
       & N_NC      = 91, &
       & N_GAMMA   = 102, &
       & N_P       = 103, &
       & N_D       = 104, &
       & N_T       = 105, &
       & N_3HE     = 106, &
       & N_A       = 107, &
       & N_2A      = 108, &
       & N_3A      = 109, &
       & N_2P      = 111, &
       & N_PA      = 112, &
       & N_T2A     = 113, &
       & N_D2A     = 114, &
       & N_PD      = 115, &
       & N_PT      = 116, &
       & N_DA      = 117

  ! Tally distinguishers
  integer, parameter :: &
       & TALLY_FLUX = -1, & ! tally the flux only
       & TALLY_ALL  =  0, & ! tally all reactions
       & TALLY_BINS =  1, & ! tally individual (reactions, cells)
       & TALLY_SUM  =  2    ! tally sum of specified (reactions, cells)

  ! Fission neutron emission (nu) type
  integer, parameter ::     &
       & NU_NONE       = 0, & ! No nu values (non-fissionable)
       & NU_POLYNOMIAL = 1, & ! Nu values given by polynomial
       & NU_TABULAR    = 2    ! Nu values given by tabular distribution

  ! Integer code for read error -- better hope this number is never
  ! used in an input file!
  integer, parameter :: &
       & ERROR_INT = -huge(0)

  real(8), parameter :: &
       & ERROR_REAL = -huge(0.0_8) * 0.917826354_8

  ! The verbosity controls how much information will be printed to the
  ! screen and in logs
  integer :: verbosity

  ! Versioning numbers
  integer, parameter :: VERSION_MAJOR = 0
  integer, parameter :: VERSION_MINOR = 2
  integer, parameter :: VERSION_RELEASE = 2

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
    if (allocated(xs_continuous)) deallocate(xs_continuous)
    if (allocated(xs_thermal)) deallocate(xs_thermal)

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
