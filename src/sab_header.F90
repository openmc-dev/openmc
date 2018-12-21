module sab_header

  use, intrinsic :: ISO_C_BINDING

  use dict_header, only: DictCharInt
  use hdf5_interface

  implicit none
  private

  public :: free_memory_sab, sab_from_hdf5, sab_has_nuclide, sab_threshold

  ! S(a,b) tables
  integer, public :: n_sab_tables
  type(DictCharInt), public :: sab_dict

  interface
    subroutine sab_from_hdf5(group_id, temperature, n) bind(C)
      import HID_T, C_DOUBLE, C_INT
      integer(HID_T), value :: group_id
      real(C_DOUBLE), intent(in) :: temperature
      integer(C_INT), value :: n
    end subroutine

    subroutine sab_clear() bind(C)
    end subroutine

    function sab_has_nuclide(i_sab, name) result(val) bind(C)
      import C_INT, C_CHAR, C_BOOL
      integer(C_INT), value :: i_sab
      character(kind=C_CHAR), intent(in) :: name(*)
      logical(C_BOOL) :: val
    end function

    function sab_threshold(i_sab) result(threshold) bind(C)
      import C_INT, C_DOUBLE
      integer(C_INT), value :: i_sab
      real(C_DOUBLE) :: threshold
    end function
  end interface

contains

!===============================================================================
! FREE_MEMORY_SAB deallocates global arrays defined in this module
!===============================================================================

  subroutine free_memory_sab()
    n_sab_tables = 0
    call sab_clear()
    call sab_dict % clear()
  end subroutine free_memory_sab

end module sab_header
