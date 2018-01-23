module surface_header

  use, intrinsic :: ISO_C_BINDING

  use constants, only: NONE
  use dict_header, only: DictIntInt

  implicit none

!===============================================================================
! SURFACE type defines a first- or second-order surface that can be used to
! construct closed volumes (cells)
!===============================================================================

  type :: Surface
    integer :: id                     ! Unique ID
    integer, allocatable :: &
         neighbor_pos(:), &           ! List of cells on positive side
         neighbor_neg(:)              ! List of cells on negative side
    integer :: bc                     ! Boundary condition
  end type Surface

  integer(C_INT32_T), bind(C) :: n_surfaces  ! # of surfaces

  type(Surface), allocatable, target :: surfaces(:)

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
