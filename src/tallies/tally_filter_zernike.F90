module tally_filter_zernike

  use, intrinsic :: ISO_C_BINDING

  use error
  use tally_filter_header

  implicit none
  private
  public :: openmc_zernike_filter_get_order
  public :: openmc_zernike_filter_get_params
  public :: openmc_zernike_filter_set_order
  public :: openmc_zernike_filter_set_params

  interface
    subroutine zernike_filter_get_params(filt, order, x, y, r) bind(C)
      import C_PTR, C_INT, C_DOUBLE
      type(C_PTR), value :: filt
      integer(C_INT)     :: order
      real(C_DOUBLE)     :: x, y, r
    end subroutine

    subroutine zernike_filter_set_params(filt, order, x, y, r) bind(C)
      import C_PTR, C_INT, C_DOUBLE
      type(C_PTR),    value :: filt
      integer(C_INT), value :: order
      real(C_DOUBLE), value :: x, y, r
    end subroutine
  end interface

!===============================================================================
! ZERNIKEFILTER gives Zernike polynomial moments of a particle's position
!===============================================================================

  type, public, extends(CppTallyFilter) :: ZernikeFilter
  end type ZernikeFilter

!===============================================================================
! ZERNIKERADIALFILTER gives even order radial Zernike polynomial moments of a
! particle's position
!===============================================================================

  type, public, extends(ZernikeFilter) :: ZernikeRadialFilter
  end type ZernikeRadialFilter

contains

!===============================================================================
!                               C API FUNCTIONS
!===============================================================================

  function openmc_zernike_filter_get_order(index, order) result(err) bind(C)
    ! Get the order of an expansion filter
    integer(C_INT32_T), value       :: index
    integer(C_INT),     intent(out) :: order
    integer(C_INT) :: err

    real(C_DOUBLE) :: x, y, r

    err = verify_filter(index)
    if (err == 0) then
      select type (f => filters(index) % obj)
      class is (ZernikeFilter)
        call zernike_filter_get_params(f % ptr, order, x, y, r)
      class default
        err = E_INVALID_TYPE
        call set_errmsg("Not a Zernike filter.")
      end select
    end if
  end function openmc_zernike_filter_get_order


  function openmc_zernike_filter_get_params(index, x, y, r) result(err) bind(C)
    ! Get the Zernike filter parameters
    integer(C_INT32_T), value :: index
    real(C_DOUBLE), intent(out) :: x
    real(C_DOUBLE), intent(out) :: y
    real(C_DOUBLE), intent(out) :: r
    integer(C_INT) :: err

    integer(C_INT) :: order

    err = verify_filter(index)
    if (err == 0) then
      select type (f => filters(index) % obj)
      class is (ZernikeFilter)
        call zernike_filter_get_params(f % ptr, order, x, y, r)
      class default
        err = E_INVALID_TYPE
        call set_errmsg("Not a Zernike filter.")
      end select
    end if
  end function openmc_zernike_filter_get_params


  function openmc_zernike_filter_set_order(index, order) result(err) bind(C)
    ! Set the order of an expansion filter
    integer(C_INT32_T), value :: index
    integer(C_INT),     value :: order
    integer(C_INT) :: err

    integer(C_INT) :: order_
    real(C_DOUBLE) :: x, y, r

    err = verify_filter(index)
    if (err == 0) then
      select type (f => filters(index) % obj)
      class is (ZernikeFilter)
        call zernike_filter_get_params(f % ptr, order_, x, y, r)
        call zernike_filter_set_params(f % ptr, order, x, y, r)
      class default
        err = E_INVALID_TYPE
        call set_errmsg("Not a Zernike filter.")
      end select
    end if
  end function openmc_zernike_filter_set_order


  function openmc_zernike_filter_set_params(index, x, y, r) result(err) bind(C)
    ! Set the Zernike filter parameters
    integer(C_INT32_T), value :: index
    real(C_DOUBLE), intent(in), optional :: x
    real(C_DOUBLE), intent(in), optional :: y
    real(C_DOUBLE), intent(in), optional :: r
    integer(C_INT) :: err

    integer(C_INT) :: order
    real(C_DOUBLE) :: x_, y_, r_

    err = verify_filter(index)
    if (err == 0) then
      select type (f => filters(index) % obj)
      class is (ZernikeFilter)
        call zernike_filter_get_params(f % ptr, order, x_, y_, r_)
        if (present(x)) x_ = x
        if (present(y)) y_ = y
        if (present(r)) r_ = r
        call zernike_filter_set_params(f % ptr, order, x_, y_, r_)
      class default
        err = E_INVALID_TYPE
        call set_errmsg("Not a Zernike filter.")
      end select
    end if
  end function openmc_zernike_filter_set_params

end module tally_filter_zernike
