module global

  use types

  implicit none

  ! Main arrays for cells, surfaces, materials
  type(Cell),     allocatable, target :: cells(:)
  type(Universe), allocatable, target :: universes(:)
  type(Lattice),  allocatable, target :: lattices(:)
  type(Surface),  allocatable, target :: surfaces(:)
  type(Material), allocatable, target :: materials(:)
  type(Isotope),  allocatable, target :: isotopes(:)
  type(xsData),   allocatable, target :: xsdatas(:)
  integer :: n_cells     ! # of cells
  integer :: n_universes ! # of universes
  integer :: n_lattices  ! # of lattices
  integer :: n_surfaces  ! # of surfaces
  integer :: n_materials ! # of materials

  ! These dictionaries provide a fast lookup mechanism
  type(DictionaryII), pointer :: cell_dict
  type(DictionaryII), pointer :: universe_dict
  type(DictionaryII), pointer :: lattice_dict
  type(DictionaryII), pointer :: surface_dict
  type(DictionaryII), pointer :: material_dict
  type(DictionaryII), pointer :: isotope_dict
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
  integer :: n_particles ! # of particles (per cycle for criticality)
  integer :: n_cycles    ! # of cycles
  integer :: n_inactive  ! # of inactive cycles

  ! External source
  type(ExtSource), target :: external_source

  ! Source and fission bank
  type(Particle), allocatable, target :: source_bank(:)
  type(Bank),     allocatable, target :: fission_bank(:)
  integer :: n_bank      ! # of sites in fission bank
  integer :: bank_first  ! index of first particle in bank
  integer :: bank_last   ! index of last particle in bank

  ! cycle keff
  real(8) :: keff

  ! Parallel processing variables
  integer :: n_procs     ! number of processes
  integer :: rank        ! rank of process
  logical :: master      ! master process?
  logical :: mpi_enabled ! is MPI in use and initialized?

  ! Paths to input file, cross section data, etc
  character(100) :: & 
       & path_input,    &
       & path_xsdata

  ! Physical constants
  real(8), parameter ::            &
       & PI           = 2.0_8*acos(0.0_8), & ! pi
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
       & BC_REFLECT  = 2     ! Reflecting boundary condition

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
  integer, parameter :: &
       & HISTOGRAM    = 1, & ! Histogram interpolation on angle
       & LINEARLINEAR = 2    ! Linear-linear interpolation on angle

  ! Particle type
  integer, parameter :: &
       & NEUTRON  = 1, &
       & PHOTON   = 2, &
       & ELECTRON = 3

  ! Fission neutron emission (nu) type
  integer, parameter ::     &
       & NU_NONE       = 0, & ! No nu values (non-fissionable)
       & NU_POLYNOMIAL = 1, & ! Nu values given by polynomial
       & NU_TABULAR    = 2    ! Nu values given by tabular distribution

  ! Integer code for read error -- better hope this number is never
  ! used in an input file!
  integer, parameter :: &
       & ERROR_CODE = -huge(0)

  ! The verbosity controls how much information will be printed to the
  ! screen and in logs
  integer :: verbosity

  integer, parameter :: max_words = 500
  integer, parameter :: max_line  = 250

  ! Versioning numbers
  integer, parameter :: VERSION_MAJOR = 0
  integer, parameter :: VERSION_MINOR = 2
  integer, parameter :: VERSION_RELEASE = 1

contains

!=====================================================================
! SET_DEFAULTS gives default values for many global parameters
!=====================================================================

  subroutine set_defaults()

    ! Default problem type is external source
    problem_type = PROB_SOURCE
    
    ! Default number of particles
    n_particles = 10000

    ! Default verbosity
    verbosity = 5

  end subroutine set_defaults

!=====================================================================
! FREE_MEMORY deallocates all allocatable arrays in the program,
! namely the cells, surfaces, materials, and sources
!=====================================================================

  subroutine free_memory()

    integer :: ierr

    ! Deallocate cells, surfaces, materials
    if (allocated(cells)) deallocate(cells)
    if (allocated(surfaces)) deallocate(surfaces)
    if (allocated(materials)) deallocate(materials)

    ! Deallocate xsdata list
    if (allocated(xsdatas)) deallocate(xsdatas)

    ! Deallocate fission and source bank
    if (allocated(fission_bank)) deallocate(fission_bank)
    if (allocated(source_bank)) deallocate(source_bank)

    ! If MPI is in use and enabled, terminate it
    call MPI_FINALIZE(ierr)

    ! End program
    stop
   
  end subroutine free_memory

!=====================================================================
! INT_TO_STR converts an integer to a string. Right now, it is limited
! to integers less than 10 billion.
!=====================================================================

  function int_to_str(num) result(str)

    integer, intent(in) :: num
    character(10) :: str

    write ( str, '(I10)' ) num
    str = adjustl(str)

  end function int_to_str

!=====================================================================
! INT_TO_STR converts a string to an integer. 
!=====================================================================

  function str_to_int(str) result(num)

    character(*), intent(in) :: str
    integer :: num
    
    character(5) :: fmt
    integer      :: w
    integer      :: ioError

    ! Determine width of string
    w = len_trim(str)
    
    ! Create format specifier for reading string
    write(UNIT=fmt, FMT='("(I",I2,")")') w

    ! read string into integer
    read(UNIT=str, FMT=fmt, IOSTAT=ioError) num
    if (ioError > 0) num = ERROR_CODE

  end function str_to_int

end module global
