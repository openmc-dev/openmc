module tally_filter_legendre

  use, intrinsic :: ISO_C_BINDING

  use tally_filter_header

  implicit none

  interface
    function openmc_legendre_filter_set_order(index, order) result(err) bind(C)
      import C_INT32_T, C_INT
      integer(C_INT32_T), value :: index
      integer(C_INT),     value :: order
      integer(C_INT) :: err
    end function
  end interface

  type, extends(TallyFilter) :: LegendreFilter
  end type LegendreFilter

end module tally_filter_legendre
