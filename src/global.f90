module global

  use types

  implicit none

  type(Cell),     allocatable, target :: cells(:)
  type(Surface),  allocatable, target :: surfaces(:)
  type(Material), allocatable, target :: materials(:)

  integer :: ncell    ! Number of cells
  integer :: nsurf    ! Number of surfaces
  integer :: nmat     ! Number of materials

  ! Physical constants
  real(8), parameter :: pi = 2.*acos(0.0)  ! pi

  ! Boundary conditions
  integer, parameter :: &
       & TRANSMIT = 0,  & ! Transmission boundary condition (default)
       & VACUUM   = 1,  & ! Vacuum boundary condition
       & REFLECT  = 2     ! Reflecting boundary condition

  ! Logical operators for cell definitions
  integer, parameter :: &
       & OP_LEFT_PAREN  = huge(0),     & ! Left parentheses
       & OP_RIGHT_PAREN = huge(0) - 1, & ! Right parentheses
       & OP_UNION       = huge(0) - 2, & ! Union operator
       & OP_DIFFERENCE  = huge(0) - 3    ! Difference operator

  ! Surface types
  integer, parameter :: &
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

  integer, parameter :: &
       & SENSE_POSITIVE = 1, &
       & SENSE_NEGATIVE = -1

  real(8), parameter :: &
       & INFINITY = huge(0.0_8)

  character(32) :: inputfile

  integer, parameter :: verbosity = 5
  integer, parameter :: max_words = 100

  ! Versioning numbers
  integer, parameter :: VERSION_MAJOR = 0
  integer, parameter :: VERSION_MINOR = 1
  integer, parameter :: VERSION_RELEASE = 1

contains

!-----------------------------------------------------------------------

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

    ! End program
    stop
   
  end subroutine free_memory

!-----------------------------------------------------------------------

end module global
