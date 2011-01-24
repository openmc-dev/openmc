module global

  use types

  implicit none

  ! Main arrays for cells, surfaces, materials
  type(Cell),      allocatable, target :: cells(:)
  type(Surface),   allocatable, target :: surfaces(:)
  type(Material),  allocatable, target :: materials(:)
  integer :: n_cells     ! # of cells
  integer :: n_surfaces  ! # of surfaces
  integer :: n_materials ! # of materials

  ! Histories/cycles/etc for both external source and criticality
  integer :: n_particles ! # of particles (per cycle for criticality)
  integer :: n_cycles    ! # of cycles
  integer :: n_inactive  ! # of inactive cycles

  ! External source
  type(ExtSource), target :: external_source

  ! Source and fission bank
  type(Neutron), allocatable, target :: source_bank(:)
  type(Bank),    allocatable, target :: fission_bank(:)

  ! Physical constants
  real(8), parameter ::           &
       & PI       = 2.*acos(0.0), &  ! pi
       & INFINITY = huge(0.0_8)      ! positive infinity

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

  ! Particle type
  integer, parameter :: &
       & NEUTRON_  = 1, &
       & PHOTON_   = 2, &
       & ELECTRON_ = 3

  character(32) :: inputfile

  integer, parameter :: verbosity = 5
  integer, parameter :: max_words = 100

  ! Versioning numbers
  integer, parameter :: VERSION_MAJOR = 0
  integer, parameter :: VERSION_MINOR = 1
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

  end subroutine set_defaults

!=====================================================================
! FREE_MEMORY deallocates all allocatable arrays in the program,
! namely the cells, surfaces, materials, and sources
!=====================================================================

  subroutine free_memory()

    ! Deallocate cells, surfaces, materials
    if (allocated(cells)) then
       deallocate(cells)
    end if
    if (allocated(surfaces)) then
       deallocate(surfaces)
    end if
    if (allocated(materials)) then
       deallocate(materials)
    end if

    ! Deallocate fission and source bank
    if (allocated(fission_bank)) then
       deallocate(fission_bank)
    end if
    if (allocated(source_bank)) then
       deallocate(source_bank)
    end if

    ! End program
    stop
   
  end subroutine free_memory

!=====================================================================
! INT_TO_STR converts an integer to a string. Right now, it is limited
! to integers less than 10 billion.
!=====================================================================

  function int_to_str( num )

    integer, intent(in) :: num
    character(10) :: int_to_str

    write ( int_to_str, '(I10)' ) num
    int_to_str = adjustl(int_to_str)

  end function int_to_str

end module global
