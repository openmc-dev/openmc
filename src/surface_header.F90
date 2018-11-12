module surface_header

  use, intrinsic :: ISO_C_BINDING

  use hdf5_interface

  implicit none

  interface
    pure function surface_pointer(surf_ind) bind(C) result(ptr)
      use ISO_C_BINDING
      implicit none
      integer(C_INT), intent(in), value :: surf_ind
      type(C_PTR)                       :: ptr
    end function surface_pointer

    pure function surface_id_c(surf_ptr) bind(C, name='surface_id') result(id)
      use ISO_C_BINDING
      implicit none
      type(C_PTR), intent(in), value :: surf_ptr
      integer(C_INT)                 :: id
    end function surface_id_c

    pure function surface_bc_c(surf_ptr) bind(C, name='surface_bc') result(bc)
      use ISO_C_BINDING
      implicit none
      type(C_PTR), intent(in), value :: surf_ptr
      integer(C_INT)                 :: bc
    end function surface_bc_c

    pure subroutine surface_reflect_c(surf_ptr, xyz, uvw) &
         bind(C, name='surface_reflect')
      use ISO_C_BINDING
      implicit none
      type(C_PTR),    intent(in), value :: surf_ptr
      real(C_DOUBLE), intent(in)        :: xyz(3);
      real(C_DOUBLE), intent(inout)     :: uvw(3);
    end subroutine surface_reflect_c

    pure function surface_i_periodic_c(surf_ptr) &
         bind(C, name="surface_i_periodic") result(i_periodic)
      use ISO_C_BINDING
      implicit none
      type(C_PTR),    intent(in), value :: surf_ptr
      integer(C_INT)                    :: i_periodic
    end function surface_i_periodic_c

    function surface_periodic_c(surf_ptr, other_ptr, xyz, uvw) &
         bind(C, name="surface_periodic") result(rotational)
      use ISO_C_BINDING
      implicit none
      type(C_PTR),    intent(in), value :: surf_ptr
      type(C_PTR),    intent(in), value :: other_ptr
      real(C_DOUBLE), intent(inout)     :: xyz(3);
      real(C_DOUBLE), intent(inout)     :: uvw(3);
      logical(C_BOOL)                   :: rotational
    end function surface_periodic_c

    subroutine free_memory_surfaces_c() bind(C)
    end subroutine free_memory_surfaces_c
  end interface

!===============================================================================
! SURFACE type defines a first- or second-order surface that can be used to
! construct closed volumes (cells)
!===============================================================================

  type :: Surface
    type(C_PTR) :: ptr

  contains

    procedure :: id => surface_id
    procedure :: bc => surface_bc
    procedure :: reflect => surface_reflect
    procedure :: i_periodic => surface_i_periodic
    procedure :: periodic_translate => surface_periodic

  end type Surface

  integer(C_INT32_T), bind(C) :: n_surfaces  ! # of surfaces

  type(Surface), allocatable, target :: surfaces(:)

contains

  pure function surface_id(this) result(id)
    class(Surface), intent(in) :: this
    integer(C_INT)             :: id
    id = surface_id_c(this % ptr)
  end function surface_id

  pure function surface_bc(this) result(bc)
    class(Surface), intent(in) :: this
    integer(C_INT)             :: bc
    bc = surface_bc_c(this % ptr)
  end function surface_bc

  pure subroutine surface_reflect(this, xyz, uvw)
    class(Surface), intent(in)    :: this
    real(C_DOUBLE), intent(in)    :: xyz(3);
    real(C_DOUBLE), intent(inout) :: uvw(3);
    call surface_reflect_c(this % ptr, xyz, uvw)
  end subroutine surface_reflect

  pure function surface_i_periodic(this) result(i_periodic)
    class(Surface), intent(in)  :: this
    integer(C_INT)              :: i_periodic
    i_periodic = surface_i_periodic_c(this % ptr)
  end function surface_i_periodic

  function surface_periodic(this, other, xyz, uvw) result(rotational)
    class(Surface), intent(in)    :: this
    class(Surface), intent(in)    :: other
    real(C_DOUBLE), intent(inout) :: xyz(3);
    real(C_DOUBLE), intent(inout) :: uvw(3);
    logical(C_BOOL)               :: rotational
    rotational = surface_periodic_c(this % ptr, other % ptr, xyz, uvw)
  end function surface_periodic

!===============================================================================
! FREE_MEMORY_SURFACES deallocates global arrays defined in this module
!===============================================================================

  subroutine free_memory_surfaces()
    if (allocated(surfaces)) deallocate(surfaces)
    call free_memory_surfaces_c()
  end subroutine free_memory_surfaces

end module surface_header
