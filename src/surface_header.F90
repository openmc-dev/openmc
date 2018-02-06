module surface_header

  use, intrinsic :: ISO_C_BINDING
  use hdf5

  use dict_header, only: DictIntInt

  implicit none

  interface
    pure function surface_pointer_c(surf_ind) &
         bind(C, name='surface_pointer') result(ptr)
      use ISO_C_BINDING
      implicit none
      integer(C_INT), intent(in), value :: surf_ind
      type(C_PTR)                       :: ptr
    end function surface_pointer_c

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

    pure function surface_sense_c(surf_ptr, xyz, uvw) &
         bind(C, name='surface_sense') result(sense)
      use ISO_C_BINDING
      implicit none
      type(C_PTR),    intent(in), value :: surf_ptr
      real(C_DOUBLE), intent(in)        :: xyz(3)
      real(C_DOUBLE), intent(in)        :: uvw(3)
      logical(C_BOOL)                   :: sense
    end function surface_sense_c

    pure subroutine surface_reflect_c(surf_ptr, xyz, uvw) &
         bind(C, name='surface_reflect')
      use ISO_C_BINDING
      implicit none
      type(C_PTR),    intent(in), value :: surf_ptr
      real(C_DOUBLE), intent(in)        :: xyz(3);
      real(C_DOUBLE), intent(inout)     :: uvw(3);
    end subroutine surface_reflect_c

    pure function surface_distance_c(surf_ptr, xyz, uvw, coincident) &
         bind(C, name='surface_distance') result(d)
      use ISO_C_BINDING
      implicit none
      type(C_PTR),     intent(in), value :: surf_ptr
      real(C_DOUBLE),  intent(in)        :: xyz(3);
      real(C_DOUBLE),  intent(in)        :: uvw(3);
      logical(C_BOOL), intent(in), value :: coincident;
      real(C_DOUBLE)                     :: d;
    end function surface_distance_c

    pure subroutine surface_normal_c(surf_ptr, xyz, uvw) &
         bind(C, name='surface_normal')
      use ISO_C_BINDING
      implicit none
      type(C_PTR),    intent(in), value :: surf_ptr
      real(C_DOUBLE), intent(in)        :: xyz(3);
      real(C_DOUBLE), intent(out)       :: uvw(3);
    end subroutine surface_normal_c

    subroutine surface_to_hdf5_c(surf_ptr, group) &
         bind(C, name='surface_to_hdf5')
      use ISO_C_BINDING
      use hdf5
      implicit none
      type(C_PTR),    intent(in), value :: surf_ptr
      integer(HID_T), intent(in), value :: group
    end subroutine surface_to_hdf5_c

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
    integer, allocatable :: &
         neighbor_pos(:), &           ! List of cells on positive side
         neighbor_neg(:)              ! List of cells on negative side
    type(C_PTR) :: ptr

  contains

    procedure :: id => surface_id
    procedure :: bc => surface_bc
    procedure :: sense => surface_sense
    procedure :: reflect => surface_reflect
    procedure :: distance => surface_distance
    procedure :: normal => surface_normal
    procedure :: to_hdf5 => surface_to_hdf5
    procedure :: i_periodic => surface_i_periodic
    procedure :: periodic_translate => surface_periodic

  end type Surface

  integer(C_INT32_T), bind(C) :: n_surfaces  ! # of surfaces

  type(Surface), allocatable, target :: surfaces(:)

  ! Dictionary that maps user IDs to indices in 'surfaces'
  type(DictIntInt) :: surface_dict

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

  pure function surface_sense(this, xyz, uvw) result(sense)
    class(Surface), intent(in) :: this
    real(C_DOUBLE), intent(in) :: xyz(3)
    real(C_DOUBLE), intent(in) :: uvw(3)
    logical(C_BOOL)            :: sense
    sense = surface_sense_c(this % ptr, xyz, uvw)
  end function surface_sense

  pure subroutine surface_reflect(this, xyz, uvw)
    class(Surface), intent(in)    :: this
    real(C_DOUBLE), intent(in)    :: xyz(3);
    real(C_DOUBLE), intent(inout) :: uvw(3);
    call surface_reflect_c(this % ptr, xyz, uvw)
  end subroutine surface_reflect

  pure function surface_distance(this, xyz, uvw, coincident) result(d)
    class(Surface),  intent(in) :: this
    real(C_DOUBLE),  intent(in) :: xyz(3);
    real(C_DOUBLE),  intent(in) :: uvw(3);
    logical(C_BOOL), intent(in) :: coincident;
    real(C_DOUBLE)              :: d;
    d = surface_distance_c(this % ptr, xyz, uvw, coincident)
  end function surface_distance

  pure subroutine surface_normal(this, xyz, uvw)
    class(Surface), intent(in)  :: this
    real(C_DOUBLE), intent(in)  :: xyz(3);
    real(C_DOUBLE), intent(out) :: uvw(3);
    call surface_normal_c(this % ptr, xyz, uvw)
  end subroutine surface_normal

  subroutine surface_to_hdf5(this, group)
    class(Surface), intent(in) :: this
    integer(HID_T), intent(in) :: group
    call surface_to_hdf5_c(this % ptr, group)
  end subroutine surface_to_hdf5

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
    call surface_dict % clear()
    call free_memory_surfaces_c()
  end subroutine free_memory_surfaces

end module surface_header
