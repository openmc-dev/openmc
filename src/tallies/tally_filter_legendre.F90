module tally_filter_legendre

  use, intrinsic :: ISO_C_BINDING

  use error
  use tally_filter_header

  implicit none
  private
  public :: openmc_legendre_filter_get_order
  public :: openmc_legendre_filter_set_order

!===============================================================================
! LEGENDREFILTER gives Legendre moments of the change in scattering angle
!===============================================================================

  type, public, extends(CppTallyFilter) :: LegendreFilter
  end type LegendreFilter

contains

!===============================================================================
!                               C API FUNCTIONS
!===============================================================================

  function openmc_legendre_filter_get_order(index, order) result(err) bind(C)
    ! Get the order of an expansion filter
    integer(C_INT32_T), value       :: index
    integer(C_INT),     intent(out) :: order
    integer(C_INT) :: err

    interface
      function legendre_filter_get_order(filt) result(order) bind(C)
        import C_PTR, C_INT
        type(C_PTR), value :: filt
        integer(C_INT)     :: order
      end function
    end interface

    err = verify_filter(index)
    if (err == 0) then
      select type (f => filters(index) % obj)
      type is (LegendreFilter)
        order = legendre_filter_get_order(f % ptr)
      class default
        err = E_INVALID_TYPE
        call set_errmsg("Tried to get order on a non-expansion filter.")
      end select
    end if
  end function openmc_legendre_filter_get_order


  function openmc_legendre_filter_set_order(index, order) result(err) bind(C)
    ! Set the order of an expansion filter
    integer(C_INT32_T), value :: index
    integer(C_INT),     value :: order
    integer(C_INT) :: err

    interface
      subroutine legendre_filter_set_order(filt, order) bind(C)
        import C_PTR, C_INT
        type(C_PTR),    value :: filt
        integer(C_INT), value :: order
      end subroutine legendre_filter_set_order
    end interface

    err = verify_filter(index)
    if (err == 0) then
      select type (f => filters(index) % obj)
      type is (LegendreFilter)
        call legendre_filter_set_order(f % ptr, order)
        f % n_bins = f % n_bins_cpp()
      class default
        err = E_INVALID_TYPE
        call set_errmsg("Tried to set order on a non-expansion filter.")
      end select
    end if
  end function openmc_legendre_filter_set_order

end module tally_filter_legendre
