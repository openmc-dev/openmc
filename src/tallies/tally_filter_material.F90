module tally_filter_material

  use, intrinsic :: ISO_C_BINDING

  use error
  use tally_filter_header

  implicit none

!===============================================================================
! MATERIAL specifies which material tally events reside in.
!===============================================================================

  type, extends(CppTallyFilter) :: MaterialFilter
  end type MaterialFilter

contains

!===============================================================================
!                               C API FUNCTIONS
!===============================================================================

  function openmc_material_filter_get_bins(index, bins, n) result(err) bind(C)
    ! Return the bins for a material filter
    integer(C_INT32_T), value :: index
    type(C_PTR), intent(out) :: bins
    integer(C_INT32_T), intent(out) :: n
    integer(C_INT) :: err

    interface
      subroutine material_filter_get_bins(filt, bins, n) bind(C)
        import C_PTR, C_INT32_T
        type(C_PTR), value :: filt
        type(C_PTR)        :: bins
        integer(C_INT32_T) :: n
      end subroutine material_filter_get_bins
    end interface

    err = verify_filter(index)
    if (err == 0) then
      select type (f => filters(index) % obj)
      type is (MaterialFilter)
        call material_filter_get_bins(f % ptr, bins, n)
        err = 0
        class default
        err = E_INVALID_TYPE
        call set_errmsg("Tried to get material filter bins on a &
             &non-material filter.")
      end select
    end if
  end function openmc_material_filter_get_bins


  function openmc_material_filter_set_bins(index, n, bins) result(err) bind(C)
    ! Set the materials for the filter
    integer(C_INT32_T), value, intent(in) :: index
    integer(C_INT32_T), value, intent(in) :: n
    integer(C_INT32_T), intent(in) :: bins(n)
    integer(C_INT) :: err

    interface
      subroutine material_filter_set_bins(filt, n, bins) bind(C)
        import C_PTR, C_INT32_T
        type(C_PTR),        value :: filt
        integer(C_INT32_T), value :: n
        integer(C_INT32_T)        :: bins(n)
      end subroutine material_filter_set_bins
    end interface

    err = verify_filter(index)
    if (err == 0) then
      select type (f => filters(index) % obj)
      type is (MaterialFilter)
        call material_filter_set_bins(f % ptr, n, bins)
        f % n_bins = f % n_bins_cpp()

        class default
        err = E_INVALID_TYPE
        call set_errmsg("Tried to set material filter bins on a &
             &non-material filter.")
      end select
    end if
  end function openmc_material_filter_set_bins

end module tally_filter_material
