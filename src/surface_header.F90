module surface_header

  use, intrinsic :: ISO_C_BINDING

  use constants, only: NONE
  use dict_header, only: DictIntInt

  implicit none

!===============================================================================
! SURFACE type defines a first- or second-order surface that can be used to
! construct closed volumes (cells)
!===============================================================================

  type, abstract :: Surface
    integer :: id                     ! Unique ID
    integer, allocatable :: &
         neighbor_pos(:), &           ! List of cells on positive side
         neighbor_neg(:)              ! List of cells on negative side
    integer :: bc                     ! Boundary condition
    integer :: i_periodic = NONE      ! Index of corresponding periodic surface
    character(len=104) :: name = ""   ! User-defined name
  end type Surface

!===============================================================================
! SURFACECONTAINER allows us to store an array of different types of surfaces
!===============================================================================

  type :: SurfaceContainer
    class(Surface), allocatable :: obj
  end type SurfaceContainer

!===============================================================================
! All the derived types below are extensions of the abstract Surface type. They
! inherent the reflect() type-bound procedure and must implement normal()
!===============================================================================

  type, extends(Surface) :: SurfaceXPlane
  end type SurfaceXPlane

  type, extends(Surface) :: SurfaceYPlane
  end type SurfaceYPlane

  type, extends(Surface) :: SurfaceZPlane
  end type SurfaceZPlane

  type, extends(Surface) :: SurfacePlane
  end type SurfacePlane

  type, extends(Surface) :: SurfaceXCylinder
  end type SurfaceXCylinder

  type, extends(Surface) :: SurfaceYCylinder
  end type SurfaceYCylinder

  type, extends(Surface) :: SurfaceZCylinder
  end type SurfaceZCylinder

  type, extends(Surface) :: SurfaceSphere
  end type SurfaceSphere

  type, extends(Surface) :: SurfaceXCone
  end type SurfaceXCone

  type, extends(Surface) :: SurfaceYCone
  end type SurfaceYCone

  type, extends(Surface) :: SurfaceZCone
  end type SurfaceZCone

  type, extends(Surface) :: SurfaceQuadric
  end type SurfaceQuadric

  integer(C_INT32_T), bind(C) :: n_surfaces  ! # of surfaces

  type(SurfaceContainer), allocatable, target :: surfaces(:)

  ! Dictionary that maps user IDs to indices in 'surfaces'
  type(DictIntInt) :: surface_dict

contains

!===============================================================================
! FREE_MEMORY_SURFACES deallocates global arrays defined in this module
!===============================================================================

  subroutine free_memory_surfaces()
    n_surfaces = 0
    if (allocated(surfaces)) deallocate(surfaces)
    call surface_dict % clear()
  end subroutine free_memory_surfaces

end module surface_header
