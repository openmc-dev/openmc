module math

  use, intrinsic :: ISO_C_BINDING

  use constants
  use random_lcg, only: prn

  implicit none
  private
  public :: t_percentile
  public :: calc_pn
  public :: calc_rn
  public :: calc_zn
  public :: calc_zn_rad
  public :: evaluate_legendre
  public :: rotate_angle
  public :: maxwell_spectrum
  public :: watt_spectrum
  public :: faddeeva
  public :: w_derivative
  public :: broaden_wmp_polynomials
  public :: spline
  public :: spline_interpolate
  public :: spline_integrate

  interface

    pure function t_percentile(p, df) bind(C, name='t_percentile_c') &
         result(t)
      use ISO_C_BINDING
      implicit none
      real(C_DOUBLE), value, intent(in) :: p
      integer(C_INT), value, intent(in) :: df
      real(C_DOUBLE) :: t
    end function t_percentile

    pure subroutine calc_pn(n, x, pnx) bind(C, name='calc_pn_c')
      use ISO_C_BINDING
      implicit none
      integer(C_INT), value, intent(in) :: n
      real(C_DOUBLE), value, intent(in) :: x
      real(C_DOUBLE), intent(out) :: pnx(n + 1)
    end subroutine calc_pn

    pure function evaluate_legendre_c_intfc(n, data, x) &
         bind(C, name='evaluate_legendre_c') result(val)
      use ISO_C_BINDING
      implicit none
      integer(C_INT), value, intent(in) :: n
      real(C_DOUBLE), intent(in) :: data(n)
      real(C_DOUBLE), value, intent(in) :: x
      real(C_DOUBLE) :: val
    end function evaluate_legendre_c_intfc

    pure subroutine calc_rn(n, uvw, rn) bind(C, name='calc_rn_c')
      use ISO_C_BINDING
      implicit none
      integer(C_INT), value, intent(in) :: n
      real(C_DOUBLE), intent(in)  :: uvw(3)
      real(C_DOUBLE), intent(out) :: rn(2 * n + 1)
    end subroutine calc_rn

    pure subroutine calc_zn(n, rho, phi, zn) bind(C, name='calc_zn_c')
      use ISO_C_BINDING
      implicit none
      integer(C_INT), value, intent(in) :: n
      real(C_DOUBLE), value, intent(in) :: rho
      real(C_DOUBLE), value, intent(in) :: phi
      real(C_DOUBLE), intent(out) :: zn(((n + 1) * (n + 2)) / 2)
    end subroutine calc_zn

    pure subroutine calc_zn_rad(n, rho, zn_rad) bind(C, name='calc_zn_rad_c')
      use ISO_C_BINDING
      implicit none
      integer(C_INT), value, intent(in) :: n
      real(C_DOUBLE), value, intent(in) :: rho
      real(C_DOUBLE), intent(out) :: zn_rad((n / 2) + 1)
    end subroutine calc_zn_rad

    subroutine rotate_angle_c_intfc(uvw, mu, phi) bind(C, name='rotate_angle_c')
      use ISO_C_BINDING
      implicit none
      real(C_DOUBLE), intent(inout) :: uvw(3)
      real(C_DOUBLE), value, intent(in)    :: mu
      real(C_DOUBLE), optional, intent(in) :: phi
    end subroutine rotate_angle_c_intfc

    function maxwell_spectrum(T) bind(C, name='maxwell_spectrum_c') &
         result(E_out)
      use ISO_C_BINDING
      implicit none
      real(C_DOUBLE), value, intent(in) :: T
      real(C_DOUBLE) :: E_out
    end function maxwell_spectrum

    function watt_spectrum(a, b) bind(C, name='watt_spectrum_c') &
         result(E_out)
      use ISO_C_BINDING
      implicit none
      real(C_DOUBLE), value, intent(in) :: a
      real(C_DOUBLE), value, intent(in) :: b
      real(C_DOUBLE) :: E_out
    end function watt_spectrum

    subroutine broaden_wmp_polynomials(E, dopp, n, factors) &
           bind(C, name='broaden_wmp_polynomials_c')
      use ISO_C_BINDING
      implicit none
      real(C_DOUBLE), value, intent(in) :: E
      real(C_DOUBLE), value, intent(in) :: dopp
      integer(C_INT), value, intent(in) :: n
      real(C_DOUBLE), intent(inout) :: factors(n)
    end subroutine broaden_wmp_polynomials

    subroutine spline(n, x, y, z) bind(C, name='spline_c')
      use ISO_C_BINDING
      implicit none
      integer(C_INT), value, intent(in) :: n
      real(C_DOUBLE), intent(in)        :: x(n)
      real(C_DOUBLE), intent(in)        :: y(n)
      real(C_DOUBLE), intent(in)        :: z(n)
    end subroutine spline

    function spline_interpolate(n, x, y, z, xint) &
         bind(C, name='spline_interpolate_c') result(yint)
      use ISO_C_BINDING
      implicit none
      integer(C_INT), value, intent(in) :: n
      real(C_DOUBLE), intent(in)        :: x(n)
      real(C_DOUBLE), intent(in)        :: y(n)
      real(C_DOUBLE), intent(in)        :: z(n)
      real(C_DOUBLE), value, intent(in) :: xint
      real(C_DOUBLE)                    :: yint
    end function spline_interpolate

    function spline_integrate(n, x, y, z, xa, xb) &
         bind(C, name='spline_integrate_c') result(s)
      use ISO_C_BINDING
      implicit none
      integer(C_INT), value, intent(in) :: n
      real(C_DOUBLE), intent(in)        :: x(n)
      real(C_DOUBLE), intent(in)        :: y(n)
      real(C_DOUBLE), intent(in)        :: z(n)
      real(C_DOUBLE), value, intent(in) :: xa
      real(C_DOUBLE), value, intent(in) :: xb
      real(C_DOUBLE)                    :: s
    end function spline_integrate

    function faddeeva_w(z, relerr) bind(C, name='Faddeeva_w') result(w)
      use ISO_C_BINDING
      implicit none
      complex(C_DOUBLE_COMPLEX), value :: z
      real(C_DOUBLE),            value :: relerr
      complex(C_DOUBLE_COMPLEX)        :: w
    end function faddeeva_w
  end interface

contains

!===============================================================================
! EVALUATE_LEGENDRE Find the value of f(x) given a set of Legendre coefficients
! and the value of x
!===============================================================================

  pure function evaluate_legendre(data, x) result(val) bind(C)
    real(C_DOUBLE), intent(in) :: data(:)
    real(C_DOUBLE), intent(in) :: x
    real(C_DOUBLE)             :: val

    val = evaluate_legendre_c_intfc(size(data) - 1, data, x)

  end function evaluate_legendre

!===============================================================================
! ROTATE_ANGLE rotates direction cosines through a polar angle whose cosine is
! mu and through an azimuthal angle sampled uniformly. Note that this is done
! with direct sampling rather than rejection as is done in MCNP and SERPENT.
!===============================================================================

  function rotate_angle(uvw0, mu, phi) result(uvw)
    real(C_DOUBLE), intent(in) :: uvw0(3)       ! directional cosine
    real(C_DOUBLE), intent(in) :: mu            ! cosine of angle in lab or CM
    real(C_DOUBLE), intent(in), optional :: phi ! azimuthal angle

    real(C_DOUBLE) :: uvw(3)  ! rotated directional cosine

    uvw = uvw0
    call rotate_angle_c_intfc(uvw, mu, phi)

  end function rotate_angle

!===============================================================================
! FADDEEVA the Faddeeva function, using Stephen Johnson's implementation
!===============================================================================

  function faddeeva(z) result(wv) bind(C)
    complex(C_DOUBLE_COMPLEX), intent(in) :: z  ! The point to evaluate Z at
    complex(C_DOUBLE_COMPLEX)             :: wv ! The resulting w(z) value
    real(C_DOUBLE) :: relerr ! Target relative error in inner loop of MIT
                             !  Faddeeva

    ! Technically, the value we want is given by the equation:
    ! w(z) = I/Pi * Integrate[Exp[-t^2]/(z-t), {t, -Infinity, Infinity}]
    ! as shown in Equation 63 from Hwang, R. N. "A rigorous pole
    ! representation of multilevel cross sections and its practical
    ! applications." Nuclear Science and Engineering 96.3 (1987): 192-209.
    !
    ! The MIT Faddeeva function evaluates w(z) = exp(-z^2)erfc(-iz). These
    ! two forms of the Faddeeva function are related by a transformation.
    !
    ! If we call the integral form w_int, and the function form w_fun:
    ! For imag(z) > 0, w_int(z) = w_fun(z)
    ! For imag(z) < 0, w_int(z) = -conjg(w_fun(conjg(z)))

    ! Note that faddeeva_w will interpret zero as machine epsilon

    relerr = ZERO
    if (aimag(z) > ZERO) then
      wv = faddeeva_w(z, relerr)
    else
      wv = -conjg(faddeeva_w(conjg(z), relerr))
    end if

  end function faddeeva

  recursive function w_derivative(z, order) result(wv) bind(C)
    complex(C_DOUBLE_COMPLEX), intent(in) :: z ! The point to evaluate Z at
    integer(C_INT),            intent(in) :: order
    complex(C_DOUBLE_COMPLEX)     :: wv     ! The resulting w(z) value

    select case(order)
    case (0)
      wv = faddeeva(z)
    case (1)
      wv = -TWO * z * faddeeva(z) + TWO * ONEI / SQRT_PI
    case default
      wv = -TWO * z * w_derivative(z, order-1) &
           - TWO * (order-1) * w_derivative(z, order-2)
    end select
  end function w_derivative

end module math
