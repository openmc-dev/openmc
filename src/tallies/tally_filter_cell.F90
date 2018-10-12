module tally_filter_cell

  use, intrinsic :: ISO_C_BINDING

  use error
  use tally_filter_header

  implicit none
  private
  public :: openmc_cell_filter_get_bins

!===============================================================================
! CELLFILTER specifies which geometric cells tally events reside in.
!===============================================================================

  type, public, extends(CppTallyFilter) :: CellFilter
  end type CellFilter

contains

!===============================================================================
!                               C API FUNCTIONS
!===============================================================================

  function openmc_cell_filter_get_bins(index, cells, n) result(err) bind(C)
    ! Return the cells associated with a cell filter
    integer(C_INT32_T), value :: index
    type(C_PTR), intent(out) :: cells
    integer(C_INT32_T), intent(out) :: n
    integer(C_INT) :: err

    interface
      subroutine cell_filter_get_bins(filt, cells, n) bind(C)
        import C_PTR, C_INT
        type(C_PTR),    intent(in), value :: filt
        type(C_PTR),    intent(out)       :: cells
        integer(C_INT), intent(out)       :: n
      end subroutine cell_filter_get_bins
    end interface

    err = verify_filter(index)
    if (err == 0) then
       select type (f => filters(index) % obj)
       type is (CellFilter)
          call cell_filter_get_bins(f % ptr, cells, n)
       class default
          err = E_INVALID_TYPE
          call set_errmsg("Tried to get cells from a non-cell filter.")
       end select
    end if
  end function openmc_cell_filter_get_bins

end module tally_filter_cell
