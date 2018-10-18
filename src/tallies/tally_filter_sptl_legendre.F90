module tally_filter_sptl_legendre

  use, intrinsic :: ISO_C_BINDING

  use error
  use tally_filter_header

  implicit none
  private
  public :: openmc_spatial_legendre_filter_get_order
  public :: openmc_spatial_legendre_filter_get_params
  public :: openmc_spatial_legendre_filter_set_order
  public :: openmc_spatial_legendre_filter_set_params

  interface
    subroutine sptl_legendre_filter_get_params(filt, order, axis, min, max) &
         bind(C)
      import C_PTR, C_INT, C_DOUBLE
      type(C_PTR), value :: filt
      integer(C_INT)     :: order
      integer(C_INT)     :: axis
      real(C_DOUBLE)     :: min
      real(C_DOUBLE)     :: max
    end subroutine

    subroutine sptl_legendre_filter_set_params(filt, order, axis, min, max) &
         bind(C)
      import C_PTR, C_INT, C_DOUBLE
      type(C_PTR),    value :: filt
      integer(C_INT), value :: order
      integer(C_INT), value :: axis
      real(C_DOUBLE), value :: min
      real(C_DOUBLE), value :: max
    end subroutine
  end interface

!===============================================================================
! SPATIALLEGENDREFILTER gives Legendre moments of the particle's normalized
! position along an axis
!===============================================================================

  type, public, extends(CppTallyFilter) :: SpatialLegendreFilter
  end type SpatialLegendreFilter

contains

!===============================================================================
!                               C API FUNCTIONS
!===============================================================================

  function openmc_spatial_legendre_filter_get_order(index, order) result(err) bind(C)
    ! Get the order of an expansion filter
    integer(C_INT32_T), value       :: index
    integer(C_INT),     intent(out) :: order
    integer(C_INT) :: err

    integer(C_INT) :: axis
    real(C_DOUBLE) :: min, max

    err = verify_filter(index)
    if (err == 0) then
      select type (f => filters(index) % obj)
      type is (SpatialLegendreFilter)
        call sptl_legendre_filter_get_params(f % ptr, order, axis, min, max)
      class default
        err = E_INVALID_TYPE
        call set_errmsg("Not a spatial Legendre filter.")
      end select
    end if
  end function openmc_spatial_legendre_filter_get_order


  function openmc_spatial_legendre_filter_set_order(index, order) result(err) bind(C)
    ! Set the order of an expansion filter
    integer(C_INT32_T), value :: index
    integer(C_INT),     value :: order
    integer(C_INT) :: err

    integer(C_INT) :: axis, order_
    real(C_DOUBLE) :: min, max

    err = verify_filter(index)
    if (err == 0) then
      select type (f => filters(index) % obj)
      type is (SpatialLegendreFilter)
        call sptl_legendre_filter_get_params(f % ptr, order_, axis, min, max)
        call sptl_legendre_filter_set_params(f % ptr, order, axis, min, max)
        f % n_bins = f % n_bins_cpp()
      class default
        err = E_INVALID_TYPE
        call set_errmsg("Not a spatial Legendre filter.")
      end select
    end if
  end function openmc_spatial_legendre_filter_set_order


  function openmc_spatial_legendre_filter_get_params(index, axis, min, max) &
       result(err) bind(C)
    ! Get the parameters for a spatial Legendre filter
    integer(C_INT32_T), value :: index
    integer(C_INT), intent(out) :: axis
    real(C_DOUBLE), intent(out) :: min
    real(C_DOUBLE), intent(out) :: max
    integer(C_INT) :: err

    integer(C_INT) :: order

    err = verify_filter(index)
    if (err == 0) then
      select type(f => filters(index) % obj)
      type is (SpatialLegendreFilter)
        call sptl_legendre_filter_get_params(f % ptr, order, axis, min, max)
      class default
        err = E_INVALID_TYPE
        call set_errmsg("Not a spatial Legendre filter.")
      end select
    end if
  end function openmc_spatial_legendre_filter_get_params


  function openmc_spatial_legendre_filter_set_params(index, axis, min, max) &
       result(err) bind(C)
    ! Set the parameters for a spatial Legendre filter
    integer(C_INT32_T), value :: index
    integer(C_INT), intent(in), optional :: axis
    real(C_DOUBLE), intent(in), optional :: min
    real(C_DOUBLE), intent(in), optional :: max
    integer(C_INT) :: err

    integer(C_INT) :: order, axis_
    real(C_DOUBLE) :: min_, max_

    err = verify_filter(index)
    if (err == 0) then
      select type(f => filters(index) % obj)
      type is (SpatialLegendreFilter)
        call sptl_legendre_filter_get_params(f % ptr, order, axis_, min_, max_)
        if (present(axis)) axis_ = axis
        if (present(min)) min_ = min
        if (present(max)) max_ = max
        call sptl_legendre_filter_set_params(f % ptr, order, axis_, min_, max_)
      class default
        err = E_INVALID_TYPE
        call set_errmsg("Not a spatial Legendre filter.")
      end select
    end if
  end function openmc_spatial_legendre_filter_set_params

end module tally_filter_sptl_legendre
