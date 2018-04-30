module math

  use, intrinsic :: ISO_C_BINDING

  use constants
  use random_lcg, only: prn

  implicit none
  private
  public :: normal_percentile
  public :: t_percentile
  public :: calc_pn
  public :: calc_rn
  public :: calc_zn
  public :: evaluate_legendre
  public :: rotate_angle
  public :: maxwell_spectrum
  public :: watt_spectrum
  public :: faddeeva
  public :: w_derivative
  public :: broaden_wmp_polynomials

!===============================================================================
! FADDEEVA_W evaluates the scaled complementary error function.  This
! interfaces with the MIT C library
!===============================================================================

  interface

    pure function t_percentile_cc(p, df) bind(C, name='t_percentile_c') &
         result(t)
      use ISO_C_BINDING
      implicit none
      real(C_DOUBLE), value, intent(in) :: p
      integer(C_INT), value, intent(in) :: df
      real(C_DOUBLE) :: t
    end function t_percentile_cc

    pure function calc_pn_cc(n, x) bind(C, name='calc_pn_c') result(pnx)
      use ISO_C_BINDING
      implicit none
      integer(C_INT), value, intent(in) :: n
      real(C_DOUBLE), value, intent(in) :: x
      real(C_DOUBLE) :: pnx
    end function calc_pn_cc

    subroutine calc_rn_cc(n, uvw, rn) bind(C, name='calc_rn_c')
      use ISO_C_BINDING
      implicit none
      integer(C_INT), value, intent(in) :: n
      real(C_DOUBLE), intent(in) :: uvw(3)
      real(C_DOUBLE), intent(in) :: rn(2 * n + 1)
    end subroutine calc_rn_cc

    pure function evaluate_legendre_cc(n, data, x) &
         bind(C, name='evaluate_legendre_c') result(val)
      use ISO_C_BINDING
      implicit none
      integer(C_INT), value, intent(in) :: n
      real(C_DOUBLE), intent(in) :: data(n)
      real(C_DOUBLE), value, intent(in) :: x
      real(C_DOUBLE) :: val
    end function evaluate_legendre_cc

    subroutine rotate_angle_cc(uvw, mu, phi) bind(C, name='rotate_angle_c')
      use ISO_C_BINDING
      implicit none
      real(C_DOUBLE), intent(inout) :: uvw(3)
      real(C_DOUBLE), value, intent(in)    :: mu
      real(C_DOUBLE), value, intent(in) :: phi
    end subroutine rotate_angle_cc

    function maxwell_spectrum_cc(T) bind(C, name='maxwell_spectrum_c') &
         result(E_out)
      use ISO_C_BINDING
      implicit none
      real(C_DOUBLE), value, intent(in) :: T
      real(C_DOUBLE) :: E_out
    end function maxwell_spectrum_cc

    function watt_spectrum_cc(a, b) bind(C, name='watt_spectrum_c') &
         result(E_out)
      use ISO_C_BINDING
      implicit none
      real(C_DOUBLE), value, intent(in) :: a
      real(C_DOUBLE), value, intent(in) :: b
      real(C_DOUBLE) :: E_out
    end function watt_spectrum_cc

    function faddeeva_cc(z) bind(C, name='faddeeva_c') result(wv)
      use ISO_C_BINDING
      implicit none
      complex(C_DOUBLE_COMPLEX), value, intent(in) :: z
      complex(C_DOUBLE_COMPLEX) :: wv
    end function faddeeva_cc

    function w_derivative_cc(z, order) bind(C, name='w_derivative_c') &
           result(wv)
      use ISO_C_BINDING
      implicit none
      complex(C_DOUBLE_COMPLEX), value, intent(in) :: z
      integer(C_INT), value,            intent(in) :: order
      complex(C_DOUBLE_COMPLEX) :: wv
    end function w_derivative_cc

    subroutine broaden_wmp_polynomials_cc(E, dopp, n, factors) &
           bind(C, name='broaden_wmp_polynomials_c')
      use ISO_C_BINDING
      implicit none
      real(C_DOUBLE), value, intent(in) :: E
      real(C_DOUBLE), value, intent(in) :: dopp
      integer(C_INT), value, intent(in) :: n
      real(C_DOUBLE), intent(inout) :: factors(n)
    end subroutine broaden_wmp_polynomials_cc

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
! NORMAL_PERCENTILE calculates the percentile of the standard normal
! distribution with a specified probability level
!===============================================================================

  pure function normal_percentile(p) result(z) bind(C)

    real(C_DOUBLE), intent(in) :: p ! probability level
    real(C_DOUBLE)             :: z ! corresponding z-value

    real(C_DOUBLE)            :: q
    real(C_DOUBLE)            :: r
    real(C_DOUBLE), parameter :: p_low  = 0.02425_8
    real(C_DOUBLE), parameter :: a(6) = (/ &
         -3.969683028665376e1_8, 2.209460984245205e2_8, -2.759285104469687e2_8, &
         1.383577518672690e2_8, -3.066479806614716e1_8, 2.506628277459239e0_8 /)
    real(C_DOUBLE), parameter :: b(5) = (/ &
         -5.447609879822406e1_8, 1.615858368580409e2_8, -1.556989798598866e2_8, &
         6.680131188771972e1_8, -1.328068155288572e1_8 /)
    real(C_DOUBLE), parameter :: c(6) = (/ &
         -7.784894002430293e-3_8, -3.223964580411365e-1_8, -2.400758277161838_8, &
         -2.549732539343734_8, 4.374664141464968_8, 2.938163982698783_8 /)
    real(C_DOUBLE), parameter :: d(4) = (/ &
         7.784695709041462e-3_8, 3.224671290700398e-1_8, &
         2.445134137142996_8,    3.754408661907416_8 /)

    ! The rational approximation used here is from an unpublished work at
    ! http://home.online.no/~pjacklam/notes/invnorm/

    if (p < p_low) then
      ! Rational approximation for lower region.

      q = sqrt(-TWO*log(p))
      z = (((((c(1)*q + c(2))*q + c(3))*q + c(4))*q + c(5))*q + c(6)) / &
           ((((d(1)*q + d(2))*q + d(3))*q + d(4))*q + ONE)

    elseif (p <= ONE - p_low) then
      ! Rational approximation for central region

      q = p - HALF
      r = q*q
      z = (((((a(1)*r + a(2))*r + a(3))*r + a(4))*r + a(5))*r + a(6))*q / &
           (((((b(1)*r + b(2))*r + b(3))*r + b(4))*r + b(5))*r + ONE)

    else
      ! Rational approximation for upper region

      q = sqrt(-TWO*log(ONE - p))
      z = -(((((c(1)*q + c(2))*q + c(3))*q + c(4))*q + c(5))*q + c(6)) / &
           ((((d(1)*q + d(2))*q + d(3))*q + d(4))*q + ONE)
    endif

    ! Refinement based on Newton's method
#ifndef NO_F2008
    z = z - (HALF * erfc(-z/sqrt(TWO)) - p) * sqrt(TWO*PI) * exp(HALF*z*z)
#endif

  end function normal_percentile

!===============================================================================
! T_PERCENTILE calculates the percentile of the Student's t distribution with a
! specified probability level and number of degrees of freedom
!===============================================================================

  pure function t_percentile(p, df) result(t) bind(C)

    real(C_DOUBLE), intent(in)  :: p  ! probability level
    integer(C_INT), intent(in)  :: df ! degrees of freedom
    real(C_DOUBLE)              :: t  ! corresponding t-value

    real(C_DOUBLE)             :: n  ! degrees of freedom as a real(8)
    real(C_DOUBLE)             :: k  ! n - 2
    real(C_DOUBLE)             :: z  ! percentile of normal distribution
    real(C_DOUBLE)             :: z2 ! z * z

    if (df == 1) then
      ! For one degree of freedom, the t-distribution becomes a Cauchy
      ! distribution whose cdf we can invert directly

      t = tan(PI*(p - HALF))

    elseif (df == 2) then
      ! For two degrees of freedom, the cdf is given by 1/2 + x/(2*sqrt(x^2 +
      ! 2)). This can be directly inverted to yield the solution below

      t = TWO*sqrt(TWO)*(p - HALF)/sqrt(ONE - FOUR*(p - HALF)**2)

    else

      ! This approximation is from E. Olusegun George and Meenakshi Sivaram, "A
      ! modification of the Fisher-Cornish approximation for the student t
      ! percentiles," Communication in Statistics - Simulation and Computation,
      ! 16 (4), pp. 1123-1132 (1987).

      n = real(df,8)
      k = ONE/(n - TWO)
      z = normal_percentile(p)
      z2 = z * z
      t = sqrt(n*k) * (z + (z2 - THREE)*z*k/FOUR + ((5._8*z2 - 56._8)*z2 + &
           75._8)*z*k*k/96._8 + (((z2 - 27._8)*THREE*z2 + 417._8)*z2 - 315._8) &
           *z*k*k*k/384._8)

    end if

  end function t_percentile

!===============================================================================
! CALC_PN calculates the n-th order Legendre polynomial at the value of x.
! Since this function is called repeatedly during the neutron transport process,
! neither n or x is checked to see if they are in the applicable range.
! This is left to the client developer to use where applicable. x is to be in
! the domain of [-1,1], and 0<=n<=5. If x is outside of the range, the return
! value will be outside the expected range; if n is outside the stated range,
! the return value will be 1.0.
!===============================================================================

  pure function calc_pn(n,x) result(pnx) bind(C)

    integer(C_INT), intent(in) :: n   ! Legendre order requested
    real(C_DOUBLE), intent(in) :: x   ! Independent variable the Legendre is to
                                      ! be evaluated at; x must be in the
                                      ! domain [-1,1]
    real(C_DOUBLE)             :: pnx ! The Legendre poly of order n evaluated
                                      ! at x

    select case(n)
    case(1)
      pnx = x
    case(2)
      pnx = 1.5_8 * x * x - HALF
    case(3)
      pnx = 2.5_8 * x * x * x - 1.5_8 * x
    case(4)
      pnx = 4.375_8 * (x ** 4) - 3.75_8 * x * x + 0.375_8
    case(5)
      pnx = 7.875_8 * (x ** 5) - 8.75_8 * x * x * x + 1.875 * x
    case(6)
      pnx = 14.4375_8 * (x ** 6) - 19.6875_8 * (x ** 4) + &
           6.5625_8 * x * x - 0.3125_8
    case(7)
      pnx = 26.8125_8 * (x ** 7) - 43.3125_8 * (x ** 5) + &
           19.6875_8 * x * x * x - 2.1875_8 * x
    case(8)
      pnx = 50.2734375_8 * (x ** 8) - 93.84375_8 * (x ** 6) + &
           54.140625 * (x ** 4) - 9.84375_8 * x * x + 0.2734375_8
    case(9)
      pnx = 94.9609375_8 * (x ** 9) - 201.09375_8 * (x ** 7) + &
           140.765625_8 * (x ** 5) - 36.09375_8 * x * x * x + 2.4609375_8 * x
    case(10)
      pnx = 180.42578125_8 * (x ** 10) - 427.32421875_8 * (x ** 8) + &
           351.9140625_8 * (x ** 6) - 117.3046875_8 * (x ** 4) + &
           13.53515625_8 * x * x - 0.24609375_8
    case default
      pnx = ONE ! correct for case(0), incorrect for the rest
    end select

  end function calc_pn

!===============================================================================
! CALC_RN calculates the n-th order real spherical harmonics for a given angle
! (in terms of (u,v,w)).  All Rn,m values are provided (where -n<=m<=n)
!===============================================================================

  subroutine calc_rn(n, uvw, rn) bind(C)

    integer(C_INT), intent(in) :: n      ! Order requested
    real(C_DOUBLE), intent(in) :: uvw(3) ! Direction of travel;
                                         ! assumed to be on unit sphere
    real(C_DOUBLE)             :: rn(2*n + 1) ! The resultant R_n(uvw)


    real(C_DOUBLE) :: phi, w ! Azimuthal and Cosine of Polar angles (from uvw)
    real(C_DOUBLE) :: w2m1   ! (w^2 - 1), frequently used in these

    w = uvw(3) ! z = cos(polar)
    if (uvw(1) == ZERO) then
      phi = ZERO
    else
      phi = atan2(uvw(2), uvw(1))
    end if

    w2m1 = (ONE - w**2)
    select case(n)
    case (0)
      ! l = 0, m = 0
      rn(1) = ONE
    case (1)
      ! l = 1, m = -1
      rn(1) = -(ONE*sqrt(w2m1) * sin(phi))
      ! l = 1, m = 0
      rn(2) = ONE * w
      ! l = 1, m = 1
      rn(3) = -(ONE*sqrt(w2m1) * cos(phi))
    case (2)
      ! l = 2, m = -2
      rn(1) = 0.288675134594813_8 * (-THREE * w**2 + THREE) * sin(TWO*phi)
      ! l = 2, m = -1
      rn(2) = -(1.73205080756888_8 * w*sqrt(w2m1) * sin(phi))
      ! l = 2, m = 0
      rn(3) = 1.5_8 * w**2 - HALF
      ! l = 2, m = 1
      rn(4) = -(1.73205080756888_8 * w*sqrt(w2m1) * cos(phi))
      ! l = 2, m = 2
      rn(5) = 0.288675134594813_8 * (-THREE * w**2 + THREE) * cos(TWO*phi)
    case (3)
      ! l = 3, m = -3
      rn(1) = -(0.790569415042095_8 * (w2m1)**(THREE/TWO) * sin(THREE * phi))
      ! l = 3, m = -2
      rn(2) = 1.93649167310371_8 * w*(w2m1) * sin(TWO*phi)
      ! l = 3, m = -1
      rn(3) = -(0.408248290463863_8*sqrt(w2m1)*((15.0_8/TWO)*w**2 - THREE/TWO) * &
           sin(phi))
      ! l = 3, m = 0
      rn(4) = 2.5_8 * w**3 - 1.5_8 * w
      ! l = 3, m = 1
      rn(5) = -(0.408248290463863_8*sqrt(w2m1)*((15.0_8/TWO)*w**2 - THREE/TWO) * &
           cos(phi))
      ! l = 3, m = 2
      rn(6) = 1.93649167310371_8 * w*(w2m1) * cos(TWO*phi)
      ! l = 3, m = 3
      rn(7) = -(0.790569415042095_8 * (w2m1)**(THREE/TWO) * cos(THREE* phi))
    case (4)
      ! l = 4, m = -4
      rn(1) = 0.739509972887452_8 * (w2m1)**2 * sin(4.0_8*phi)
      ! l = 4, m = -3
      rn(2) = -(2.09165006633519_8 * w*(w2m1)**(THREE/TWO) * sin(THREE* phi))
      ! l = 4, m = -2
      rn(3) = 0.074535599249993_8 * (w2m1)*((105.0_8/TWO)*w**2 - 15.0_8/TWO) * &
           sin(TWO*phi)
      ! l = 4, m = -1
      rn(4) = -(0.316227766016838_8*sqrt(w2m1)*((35.0_8/TWO)*w**3 - 15.0_8/TWO*w)&
           * sin(phi))
      ! l = 4, m = 0
      rn(5) = 4.375_8 * w**4 - 3.75_8 * w**2 + 0.375_8
      ! l = 4, m = 1
      rn(6) = -(0.316227766016838_8*sqrt(w2m1)*((35.0_8/TWO)*w**3 - 15.0_8/TWO*w)&
           * cos(phi))
      ! l = 4, m = 2
      rn(7) = 0.074535599249993_8 * (w2m1)*((105.0_8/TWO)*w**2 - 15.0_8/TWO) * &
           cos(TWO*phi)
      ! l = 4, m = 3
      rn(8) = -(2.09165006633519_8 * w*(w2m1)**(THREE/TWO) * cos(THREE* phi))
      ! l = 4, m = 4
      rn(9) = 0.739509972887452_8 * (w2m1)**2 * cos(4.0_8*phi)
    case (5)
      ! l = 5, m = -5
      rn(1) = -(0.701560760020114_8 * (w2m1)**(5.0_8/TWO) * sin(5.0_8*phi))
      ! l = 5, m = -4
      rn(2) = 2.21852991866236_8 * w*(w2m1)**2 * sin(4.0_8*phi)
      ! l = 5, m = -3
      rn(3) = -(0.00996023841111995_8 * (w2m1)**(THREE/TWO)* &
           ((945.0_8 /TWO)*w**2 - 105.0_8/TWO) * sin(THREE*phi))
      ! l = 5, m = -2
      rn(4) = 0.0487950036474267_8 * (w2m1) &
           * ((315.0_8/TWO)*w**3 - 105.0_8/TWO*w) * sin(TWO*phi)
      ! l = 5, m = -1
      rn(5) = -(0.258198889747161_8*sqrt(w2m1)* &
           ((315.0_8/8.0_8)*w**4 - 105.0_8/4.0_8 * w**2 + 15.0_8/8.0_8) &
           * sin(phi))
      ! l = 5, m = 0
      rn(6) = 7.875_8 * w**5 - 8.75_8 * w**3 + 1.875_8 * w
      ! l = 5, m = 1
      rn(7) = -(0.258198889747161_8*sqrt(w2m1)* &
           ((315.0_8/8.0_8)*w**4 - 105.0_8/4.0_8 * w**2 + 15.0_8/8.0_8) &
           * cos(phi))
      ! l = 5, m = 2
      rn(8) = 0.0487950036474267_8 * (w2m1)* &
           ((315.0_8/TWO)*w**3 - 105.0_8/TWO*w) * cos(TWO*phi)
      ! l = 5, m = 3
      rn(9) = -(0.00996023841111995_8 * (w2m1)**(THREE/TWO)* &
           ((945.0_8 /TWO)*w**2 - 105.0_8/TWO) * cos(THREE*phi))
      ! l = 5, m = 4
      rn(10) = 2.21852991866236_8 * w*(w2m1)**2 * cos(4.0_8*phi)
      ! l = 5, m = 5
      rn(11) = -(0.701560760020114_8 * (w2m1)**(5.0_8/TWO) * cos(5.0_8* phi))
    case (6)
      ! l = 6, m = -6
      rn(1) = 0.671693289381396_8 * (w2m1)**3 * sin(6.0_8*phi)
      ! l = 6, m = -5
      rn(2) = -(2.32681380862329_8 * w*(w2m1)**(5.0_8/TWO) * sin(5.0_8*phi))
      ! l = 6, m = -4
      rn(3) = 0.00104990131391452_8 * (w2m1)**2 * &
           ((10395.0_8/TWO)*w**2 - 945.0_8/TWO) * sin(4.0_8*phi)
      ! l = 6, m = -3
      rn(4) = -(0.00575054632785295_8 * (w2m1)**(THREE/TWO) * &
           ((3465.0_8/TWO)*w**3 - 945.0_8/TWO*w) * sin(THREE*phi))
      ! l = 6, m = -2
      rn(5) = 0.0345032779671177_8 * (w2m1) * &
           ((3465.0_8/8.0_8)*w**4 - 945.0_8/4.0_8 * w**2 + 105.0_8/8.0_8) &
           * sin(TWO*phi)
      ! l = 6, m = -1
      rn(6) = -(0.218217890235992_8*sqrt(w2m1) * &
           ((693.0_8/8.0_8)*w**5- 315.0_8/4.0_8 * w**3 + (105.0_8/8.0_8)*w) &
           * sin(phi))
      ! l = 6, m = 0
      rn(7) = 14.4375_8 * w**6 - 19.6875_8 * w**4 + 6.5625_8 * w**2 - 0.3125_8
      ! l = 6, m = 1
      rn(8) = -(0.218217890235992_8*sqrt(w2m1) * &
           ((693.0_8/8.0_8)*w**5- 315.0_8/4.0_8 * w**3 + (105.0_8/8.0_8)*w) &
           * cos(phi))
      ! l = 6, m = 2
      rn(9) = 0.0345032779671177_8 * (w2m1) * &
           ((3465.0_8/8.0_8)*w**4 -945.0_8/4.0_8 * w**2 + 105.0_8/8.0_8) &
           * cos(TWO*phi)
      ! l = 6, m = 3
      rn(10) = -(0.00575054632785295_8 * (w2m1)**(THREE/TWO) * &
           ((3465.0_8/TWO)*w**3 - 945.0_8/TWO*w) * cos(THREE*phi))
      ! l = 6, m = 4
      rn(11) = 0.00104990131391452_8 * (w2m1)**2 * &
           ((10395.0_8/TWO)*w**2 - 945.0_8/TWO) * cos(4.0_8*phi)
      ! l = 6, m = 5
      rn(12) = -(2.32681380862329_8 * w*(w2m1)**(5.0_8/TWO) * cos(5.0_8*phi))
      ! l = 6, m = 6
      rn(13) = 0.671693289381396_8 * (w2m1)**3 * cos(6.0_8*phi)
    case (7)
      ! l = 7, m = -7
      rn(1) = -(0.647259849287749_8 * (w2m1)**(7.0_8/TWO) * sin(7.0_8*phi))
      ! l = 7, m = -6
      rn(2) = 2.42182459624969_8 * w*(w2m1)**3 * sin(6.0_8*phi)
      ! l = 7, m = -5
      rn(3) = -(9.13821798555235d-5*(w2m1)**(5.0_8/TWO)* &
                 ((135135.0_8/TWO)*w**2 - 10395.0_8/TWO) * sin(5.0_8*phi))
      ! l = 7, m = -4
      rn(4) = 0.000548293079133141_8 * (w2m1)**2* &
           ((45045.0_8/TWO)*w**3 - 10395.0_8/TWO*w) * sin(4.0_8*phi)
      ! l = 7, m = -3
      rn(5) = -(0.00363696483726654_8 * (w2m1)**(THREE/TWO)* &
                 ((45045.0_8/8.0_8)*w**4 - 10395.0_8/4.0_8 * w**2 + 945.0_8/8.0_8)* &
                 sin(THREE*phi))
      ! l = 7, m = -2
      rn(6) = 0.025717224993682_8 * (w2m1)* &
           ((9009.0_8/8.0_8)*w**5 -3465.0_8/4.0_8 * w**3 + (945.0_8/8.0_8)*w)* &
           sin(TWO*phi)
      ! l = 7, m = -1
      rn(7) = -(0.188982236504614_8*sqrt(w2m1)* &
                 ((3003.0_8/16.0_8)*w**6 - 3465.0_8/16.0_8 * w**4 + &
                 (945.0_8/16.0_8)*w**2 - 35.0_8/16.0_8) * sin(phi))
      ! l = 7, m = 0
      rn(8) = 26.8125_8 * w**7 - 43.3125_8 * w**5 + 19.6875_8 * w**3 -2.1875_8 &
           * w
      ! l = 7, m = 1
      rn(9) = -(0.188982236504614_8*sqrt(w2m1)* &
                 ((3003.0_8/16.0_8)*w**6 - 3465.0_8/16.0_8 * w**4 + &
                  (945.0_8/16.0_8)*w**2 - 35.0_8/16.0_8) * cos(phi))
      ! l = 7, m = 2
      rn(10) = 0.025717224993682_8 * (w2m1)* &
           ((9009.0_8/8.0_8)*w**5 -3465.0_8/4.0_8 * w**3 + (945.0_8/8.0_8)*w)* &
           cos(TWO*phi)
      ! l = 7, m = 3
      rn(11) = -(0.00363696483726654_8 * (w2m1)**(THREE/TWO)* &
                 ((45045.0_8/8.0_8)*w**4 - 10395.0_8/4.0_8 * w**2 + 945.0_8/8.0_8)* &
                 cos(THREE*phi))
      ! l = 7, m = 4
      rn(12) = 0.000548293079133141_8 * (w2m1)**2 * &
           ((45045.0_8/TWO)*w**3 - 10395.0_8/TWO*w) * cos(4.0_8*phi)
      ! l = 7, m = 5
      rn(13) = -(9.13821798555235d-5*(w2m1)**(5.0_8/TWO)* &
                 ((135135.0_8/TWO)*w**2 - 10395.0_8/TWO) * cos(5.0_8*phi))
      ! l = 7, m = 6
      rn(14) = 2.42182459624969_8 * w*(w2m1)**3 * cos(6.0_8*phi)
      ! l = 7, m = 7
      rn(15) = -(0.647259849287749_8 * (w2m1)**(7.0_8/TWO) * cos(7.0_8*phi))
    case (8)
      ! l = 8, m = -8
      rn(1) = 0.626706654240044_8 * (w2m1)**4 * sin(8.0_8*phi)
      ! l = 8, m = -7
      rn(2) = -(2.50682661696018_8 * w*(w2m1)**(7.0_8/TWO) * sin(7.0_8*phi))
      ! l = 8, m = -6
      rn(3) = 6.77369783729086d-6*(w2m1)**3* &
           ((2027025.0_8/TWO)*w**2 - 135135.0_8/TWO) * sin(6.0_8*phi)
      ! l = 8, m = -5
      rn(4) = -(4.38985792528482d-5*(w2m1)**(5.0_8/TWO)* &
                 ((675675.0_8/TWO)*w**3 - 135135.0_8/TWO*w) * sin(5.0_8*phi))
      ! l = 8, m = -4
      rn(5) = 0.000316557156832328_8 * (w2m1)**2* &
           ((675675.0_8/8.0_8)*w**4 - 135135.0_8/4.0_8 * w**2 &
           + 10395.0_8/8.0_8) * sin(4.0_8*phi)
      ! l = 8, m = -3
      rn(6) = -(0.00245204119306875_8 * (w2m1)**(THREE/TWO)* &
                 ((135135.0_8/8.0_8)*w**5 - 45045.0_8/4.0_8 * w**3 &
                 + (10395.0_8/8.0_8)*w) * sin(THREE*phi))
      ! l = 8, m = -2
      rn(7) = 0.0199204768222399_8 * (w2m1)* &
           ((45045.0_8/16.0_8)*w**6- 45045.0_8/16.0_8 * w**4 + &
           (10395.0_8/16.0_8)*w**2 - 315.0_8/16.0_8) * sin(TWO*phi)
      ! l = 8, m = -1
      rn(8) = -(0.166666666666667_8*sqrt(w2m1)* &
                 ((6435.0_8/16.0_8)*w**7 - 9009.0_8/16.0_8 * w**5 + &
                 (3465.0_8/16.0_8)*w**3 - 315.0_8/16.0_8 * w) * sin(phi))
      ! l = 8, m = 0
      rn(9) = 50.2734375_8 * w**8 - 93.84375_8 * w**6 + 54.140625_8 * w**4 -&
           9.84375_8 * w**2 + 0.2734375_8
      ! l = 8, m = 1
      rn(10) = -(0.166666666666667_8*sqrt(w2m1)* &
                 ((6435.0_8/16.0_8)*w**7 - 9009.0_8/16.0_8 * w**5 + &
                 (3465.0_8/16.0_8)*w**3 - 315.0_8/16.0_8 * w) * cos(phi))
      ! l = 8, m = 2
      rn(11) = 0.0199204768222399_8 * (w2m1)*((45045.0_8/16.0_8)*w**6- &
           45045.0_8/16.0_8 * w**4 + (10395.0_8/16.0_8)*w**2 - &
           315.0_8/16.0_8) * cos(TWO*phi)
      ! l = 8, m = 3
      rn(12) = -(0.00245204119306875_8 * (w2m1)**(THREE/TWO)* &
                 ((135135.0_8/8.0_8)*w**5 - 45045.0_8/4.0_8 * w**3 + &
                 (10395.0_8/8.0_8)*w) * cos(THREE*phi))
      ! l = 8, m = 4
      rn(13) = 0.000316557156832328_8 * (w2m1)**2*((675675.0_8/8.0_8)*w**4 - &
           135135.0_8/4.0_8 * w**2 + 10395.0_8/8.0_8) * cos(4.0_8*phi)
      ! l = 8, m = 5
      rn(14) = -(4.38985792528482d-5*(w2m1)**(5.0_8/TWO)*((675675.0_8/TWO)*w**3 -&
                 135135.0_8/TWO*w) * cos(5.0_8*phi))
      ! l = 8, m = 6
      rn(15) = 6.77369783729086d-6*(w2m1)**3*((2027025.0_8/TWO)*w**2 - &
           135135.0_8/TWO) * cos(6.0_8*phi)
      ! l = 8, m = 7
      rn(16) = -(2.50682661696018_8 * w*(w2m1)**(7.0_8/TWO) * cos(7.0_8*phi))
      ! l = 8, m = 8
      rn(17) = 0.626706654240044_8 * (w2m1)**4 * cos(8.0_8*phi)
    case (9)
      ! l = 9, m = -9
      rn(1) = -(0.609049392175524_8 * (w2m1)**(9.0_8/TWO) * sin(9.0_8*phi))
      ! l = 9, m = -8
      rn(2) = 2.58397773170915_8 * w*(w2m1)**4 * sin(8.0_8*phi)
      ! l = 9, m = -7
      rn(3) = -(4.37240315267812d-7*(w2m1)**(7.0_8/TWO)* &
                 ((34459425.0_8/TWO)*w**2 - 2027025.0_8/TWO) * sin(7.0_8*phi))
      ! l = 9, m = -6
      rn(4) = 3.02928976464514d-6*(w2m1)**3* &
           ((11486475.0_8/TWO)*w**3 - 2027025.0_8/TWO*w) * sin(6.0_8*phi)
      ! l = 9, m = -5
      rn(5) = -(2.34647776186144d-5*(w2m1)**(5.0_8/TWO)* &
                 ((11486475.0_8/8.0_8)*w**4 - 2027025.0_8/4.0_8 * w**2 + &
                 135135.0_8/8.0_8) * sin(5.0_8*phi))
      ! l = 9, m = -4
      rn(6) = 0.000196320414650061_8 * (w2m1)**2*((2297295.0_8/8.0_8)*w**5 - &
           675675.0_8/4.0_8 * w**3 + (135135.0_8/8.0_8)*w) * sin(4.0_8*phi)
      ! l = 9, m = -3
      rn(7) = -(0.00173385495536766_8 * (w2m1)**(THREE/TWO)* &
                 ((765765.0_8/16.0_8)*w**6 - 675675.0_8/16.0_8 * w**4 + &
                 (135135.0_8/16.0_8)*w**2 - 3465.0_8/16.0_8) * sin(THREE*phi))
      ! l = 9, m = -2
      rn(8) = 0.0158910431540932_8 * (w2m1)*((109395.0_8/16.0_8)*w**7- &
           135135.0_8/16.0_8 * w**5 + (45045.0_8/16.0_8)*w**3 &
           - 3465.0_8/16.0_8 * w) * sin(TWO*phi)
      ! l = 9, m = -1
      rn(9) = -(0.149071198499986_8*sqrt(w2m1)*((109395.0_8/128.0_8)*w**8 - &
                 45045.0_8/32.0_8 * w**6 + (45045.0_8/64.0_8)*w**4 - 3465.0_8/32.0_8 &
                 * w**2 + 315.0_8/128.0_8) * sin(phi))
      ! l = 9, m = 0
      rn(10) = 94.9609375_8 * w**9 - 201.09375_8 * w**7 + 140.765625_8 * w**5- &
           36.09375_8 * w**3 + 2.4609375_8 * w
      ! l = 9, m = 1
      rn(11) = -(0.149071198499986_8*sqrt(w2m1)*((109395.0_8/128.0_8)*w**8 - &
                 45045.0_8/32.0_8 * w**6 + (45045.0_8/64.0_8)*w**4 -3465.0_8/32.0_8 &
                 * w**2 + 315.0_8/128.0_8) * cos(phi))
      ! l = 9, m = 2
      rn(12) = 0.0158910431540932_8 * (w2m1)*((109395.0_8/16.0_8)*w**7 - &
           135135.0_8/16.0_8 * w**5 + (45045.0_8/16.0_8)*w**3 &
           - 3465.0_8/ 16.0_8 * w) * cos(TWO*phi)
      ! l = 9, m = 3
      rn(13) = -(0.00173385495536766_8 * (w2m1)**(THREE/TWO)*((765765.0_8/16.0_8)&
                 *w**6 - 675675.0_8/16.0_8 * w**4 + (135135.0_8/16.0_8)*w**2 &
                 - 3465.0_8/16.0_8)* cos(THREE*phi))
      ! l = 9, m = 4
      rn(14) = 0.000196320414650061_8 * (w2m1)**2*((2297295.0_8/8.0_8)*w**5 - &
           675675.0_8/4.0_8 * w**3 + (135135.0_8/8.0_8)*w) * cos(4.0_8*phi)
      ! l = 9, m = 5
      rn(15) = -(2.34647776186144d-5*(w2m1)**(5.0_8/TWO)*((11486475.0_8/8.0_8)* &
                 w**4 - 2027025.0_8/4.0_8 * w**2 + 135135.0_8/8.0_8) * cos(5.0_8*phi))
      ! l = 9, m = 6
      rn(16) = 3.02928976464514d-6*(w2m1)**3*((11486475.0_8/TWO)*w**3 - &
           2027025.0_8/TWO*w) * cos(6.0_8*phi)
      ! l = 9, m = 7
      rn(17) = -(4.37240315267812d-7*(w2m1)**(7.0_8/TWO)* &
                 ((34459425.0_8/TWO)*w**2 - 2027025.0_8/TWO) * cos(7.0_8*phi))
      ! l = 9, m = 8
      rn(18) = 2.58397773170915_8 * w*(w2m1)**4 * cos(8.0_8*phi)
      ! l = 9, m = 9
      rn(19) = -(0.609049392175524_8 * (w2m1)**(9.0_8/TWO) * cos(9.0_8*phi))
    case (10)
      ! l = 10, m = -10
      rn(1) = 0.593627917136573_8 * (w2m1)**5 * sin(10.0_8*phi)
      ! l = 10, m = -9
      rn(2) = -(2.65478475211798_8 * w*(w2m1)**(9.0_8/TWO) * sin(9.0_8*phi))
      ! l = 10, m = -8
      rn(3) = 2.49953651452314d-8*(w2m1)**4*((654729075.0_8/TWO)*w**2 - &
           34459425.0_8/TWO) * sin(8.0_8*phi)
      ! l = 10, m = -7
      rn(4) = -(1.83677671621093d-7*(w2m1)**(7.0_8/TWO)* &
                 ((218243025.0_8/TWO)*w**3 - 34459425.0_8/TWO*w) * sin(7.0_8*phi))
      ! l = 10, m = -6
      rn(5) = 1.51464488232257d-6*(w2m1)**3*((218243025.0_8/8.0_8)*w**4 - &
           34459425.0_8/4.0_8 * w**2 + 2027025.0_8/8.0_8) * sin(6.0_8*phi)
      ! l = 10, m = -5
      rn(6) = -(1.35473956745817d-5*(w2m1)**(5.0_8/TWO)* &
                 ((43648605.0_8/8.0_8)*w**5 - 11486475.0_8/4.0_8 * w**3 + &
                 (2027025.0_8/8.0_8)*w) * sin(5.0_8*phi))
      ! l = 10, m = -4
      rn(7) = 0.000128521880085575_8 * (w2m1)**2*((14549535.0_8/16.0_8)*w**6 - &
           11486475.0_8/16.0_8 * w**4 + (2027025.0_8/16.0_8)*w**2 - &
           45045.0_8/16.0_8) * sin(4.0_8*phi)
      ! l = 10, m = -3
      rn(8) = -(0.00127230170115096_8 * (w2m1)**(THREE/TWO)* &
                 ((2078505.0_8/16.0_8)*w**7 - 2297295.0_8/16.0_8 * w**5 + &
                 (675675.0_8/16.0_8)*w**3 - 45045.0_8/16.0_8 * w) * sin(THREE*phi))
      ! l = 10, m = -2
      rn(9) = 0.012974982402692_8 * (w2m1)*((2078505.0_8/128.0_8)*w**8 - &
           765765.0_8/32.0_8 * w**6 + (675675.0_8/64.0_8)*w**4 - &
           45045.0_8/32.0_8 * w**2 + 3465.0_8/128.0_8) * sin(TWO*phi)
      ! l = 10, m = -1
      rn(10) = -(0.134839972492648_8*sqrt(w2m1)*((230945.0_8/128.0_8)*w**9 - &
                 109395.0_8/32.0_8 * w**7 + (135135.0_8/64.0_8)*w**5 - &
                 15015.0_8/32.0_8 * w**3 + (3465.0_8/128.0_8)*w) * sin(phi))
      ! l = 10, m = 0
      rn(11) = 180.42578125_8 * w**10 - 427.32421875_8 * w**8 +351.9140625_8 &
           * w**6 - 117.3046875_8 * w**4 + 13.53515625_8 * w**2 -0.24609375_8
      ! l = 10, m = 1
      rn(12) = -(0.134839972492648_8*sqrt(w2m1)*((230945.0_8/128.0_8)*w**9 - &
                 109395.0_8/32.0_8 * w**7 + (135135.0_8/64.0_8)*w**5 -15015.0_8/ &
                 32.0_8 * w**3 + (3465.0_8/128.0_8)*w) * cos(phi))
      ! l = 10, m = 2
      rn(13) = 0.012974982402692_8 * (w2m1)*((2078505.0_8/128.0_8)*w**8 - &
           765765.0_8/32.0_8 * w**6 + (675675.0_8/64.0_8)*w**4 -&
           45045.0_8/32.0_8 * w**2 + 3465.0_8/128.0_8) * cos(TWO*phi)
      ! l = 10, m = 3
      rn(14) = -(0.00127230170115096_8 * (w2m1)**(THREE/TWO)* &
                 ((2078505.0_8/16.0_8)*w**7 - 2297295.0_8/16.0_8 * w**5 + &
                 (675675.0_8/16.0_8)*w**3 - 45045.0_8/16.0_8 * w) * cos(THREE*phi))
      ! l = 10, m = 4
      rn(15) = 0.000128521880085575_8 * (w2m1)**2*((14549535.0_8/16.0_8)*w**6 -&
           11486475.0_8/16.0_8 * w**4 + (2027025.0_8/16.0_8)*w**2 - &
           45045.0_8/16.0_8) * cos(4.0_8*phi)
      ! l = 10, m = 5
      rn(16) = -(1.35473956745817d-5*(w2m1)**(5.0_8/TWO)* &
                 ((43648605.0_8/8.0_8)*w**5 - 11486475.0_8/4.0_8 * w**3 + &
                 (2027025.0_8/8.0_8)*w) * cos(5.0_8*phi))
      ! l = 10, m = 6
      rn(17) = 1.51464488232257d-6*(w2m1)**3*((218243025.0_8/8.0_8)*w**4 - &
           34459425.0_8/4.0_8 * w**2 + 2027025.0_8/8.0_8) * cos(6.0_8*phi)
      ! l = 10, m = 7
      rn(18) = -(1.83677671621093d-7*(w2m1)**(7.0_8/TWO)* &
                 ((218243025.0_8/TWO)*w**3 - 34459425.0_8/TWO*w) * cos(7.0_8*phi))
      ! l = 10, m = 8
      rn(19) = 2.49953651452314d-8*(w2m1)**4* &
           ((654729075.0_8/TWO)*w**2 - 34459425.0_8/TWO) * cos(8.0_8*phi)
      ! l = 10, m = 9
      rn(20) = -(2.65478475211798_8 * w*(w2m1)**(9.0_8/TWO) * cos(9.0_8*phi))
      ! l = 10, m = 10
      rn(21) = 0.593627917136573_8 * (w2m1)**5 * cos(10.0_8*phi)
    case default
      rn = ONE
    end select

  end subroutine calc_rn

!===============================================================================
! CALC_ZN calculates the n-th order modified Zernike polynomial moment for a
! given angle (rho, theta) location in the unit disk. The normlization of the
! polynomials is such that the integral of Z_pq*Z_pq over the unit disk is
! exactly pi
!===============================================================================

  subroutine calc_zn(n, rho, phi, zn) bind(C)
    ! This procedure uses the modified Kintner's method for calculating Zernike
    ! polynomials as outlined in Chong, C. W., Raveendran, P., & Mukundan,
    ! R. (2003). A comparative analysis of algorithms for fast computation of
    ! Zernike moments. Pattern Recognition, 36(3), 731-742.

    integer(C_INT), intent(in) :: n      ! Maximum order
    real(C_DOUBLE), intent(in) :: rho    ! Radial location in the unit disk
    real(C_DOUBLE), intent(in) :: phi    ! Theta (radians) location in the unit disk
    real(C_DOUBLE), intent(out) :: zn(:) ! The resulting list of coefficients

    real(C_DOUBLE) :: sin_phi, cos_phi   ! Sine and Cosine of phi
    real(C_DOUBLE) :: sin_phi_vec(n+1)   ! Contains sin(n*phi)
    real(C_DOUBLE) :: cos_phi_vec(n+1)   ! Contains cos(n*phi)
    real(C_DOUBLE) :: zn_mat(n+1, n+1)   ! Matrix form of the coefficients which is
                                         ! easier to work with
    real(C_DOUBLE) :: k1, k2, k3, k4     ! Variables for R_m_n calculation
    real(C_DOUBLE) :: sqrt_norm          ! normalization for radial moments
    integer(C_INT) :: i,p,q              ! Loop counters

    real(C_DOUBLE), parameter :: SQRT_N_1(0:10) = [&
         sqrt(1.0_8), sqrt(2.0_8), sqrt(3.0_8), sqrt(4.0_8), &
         sqrt(5.0_8), sqrt(6.0_8), sqrt(7.0_8), sqrt(8.0_8), &
         sqrt(9.0_8), sqrt(10.0_8), sqrt(11.0_8)]
    real(C_DOUBLE), parameter :: SQRT_2N_2(0:10) = SQRT_N_1*sqrt(2.0_8)

    ! n == radial degree
    ! m == azimuthal frequency

    ! ==========================================================================
    ! Determine vector of sin(n*phi) and cos(n*phi). This takes advantage of the
    ! following recurrence relations so that only a single sin/cos have to be
    ! evaluated (http://mathworld.wolfram.com/Multiple-AngleFormulas.html)
    !
    ! sin(nx) = 2 cos(x) sin((n-1)x) - sin((n-2)x)
    ! cos(nx) = 2 cos(x) cos((n-1)x) - cos((n-2)x)

    sin_phi = sin(phi)
    cos_phi = cos(phi)

    sin_phi_vec(1) = 1.0_8
    cos_phi_vec(1) = 1.0_8

    sin_phi_vec(2) = 2.0_8 * cos_phi
    cos_phi_vec(2) = cos_phi

    do i = 3, n+1
      sin_phi_vec(i) = 2.0_8 * cos_phi * sin_phi_vec(i-1) - sin_phi_vec(i-2)
      cos_phi_vec(i) = 2.0_8 * cos_phi * cos_phi_vec(i-1) - cos_phi_vec(i-2)
    end do

    do i = 1, n+1
      sin_phi_vec(i) = sin_phi_vec(i) * sin_phi
    end do

    ! ==========================================================================
    ! Calculate R_pq(rho)

    ! Fill the main diagonal first (Eq. 3.9 in Chong)
    do p = 0, n
      zn_mat(p+1, p+1) = rho**p
    end do

    ! Fill in the second diagonal (Eq. 3.10 in Chong)
    do q = 0, n-2
      zn_mat(q+2+1, q+1) = (q+2) * zn_mat(q+2+1, q+2+1) - (q+1) * zn_mat(q+1, q+1)
    end do

    ! Fill in the rest of the values using the original results (Eq. 3.8 in Chong)
    do p = 4, n
      k2 = 2 * p * (p - 1) * (p - 2)
      do q = p-4, 0, -2
        k1 = (p + q) * (p - q) * (p - 2) / 2
        k3 = -q**2*(p - 1) - p * (p - 1) * (p - 2)
        k4 = -p * (p + q - 2) * (p - q - 2) / 2
        zn_mat(p+1, q+1) = ((k2 * rho**2 + k3) * zn_mat(p-2+1, q+1) + k4 * zn_mat(p-4+1, q+1)) / k1
      end do
    end do

    ! Roll into a single vector for easier computation later
    ! The vector is ordered (0,0), (1,-1), (1,1), (2,-2), (2,0),
    ! (2, 2), ....   in (n,m) indices
    ! Note that the cos and sin vectors are offset by one
    ! sin_phi_vec = [sin(x), sin(2x), sin(3x) ...]
    ! cos_phi_vec = [1.0, cos(x), cos(2x)... ]
    i = 1
    do p = 0, n
      do q = -p, p, 2
        if (q < 0) then
          zn(i) = zn_mat(p+1, abs(q)+1) * sin_phi_vec(abs(q)) * SQRT_2N_2(p)
        else if (q == 0) then
          zn(i) = zn_mat(p+1, q+1) * SQRT_N_1(p)
        else
          zn(i) = zn_mat(p+1, q+1) * cos_phi_vec(abs(q)+1) * SQRT_2N_2(p)
        end if
        i = i + 1
      end do
    end do
  end subroutine calc_zn

!===============================================================================
! EVALUATE_LEGENDRE Find the value of f(x) given a set of Legendre coefficients
! and the value of x
!===============================================================================

  pure function evaluate_legendre(n, data, x) result(val) bind(C)
    integer(C_INT), intent(in) :: n
    real(C_DOUBLE), intent(in) :: data(n)
    real(C_DOUBLE), intent(in) :: x
    real(C_DOUBLE)             :: val

    integer(C_INT) :: l

    val =  HALF * data(1)
    do l = 1, n - 1
      val = val + (real(l, 8) +  HALF) * data(l + 1) * calc_pn(l,x)
    end do

  end function evaluate_legendre

!===============================================================================
! ROTATE_ANGLE rotates direction cosines through a polar angle whose cosine is
! mu and through an azimuthal angle sampled uniformly. Note that this is done
! with direct sampling rather than rejection as is done in MCNP and SERPENT.
!===============================================================================

  subroutine rotate_angle(uvw0, mu, uvw, phi) bind(C)
    real(C_DOUBLE), intent(in)  :: uvw0(3) ! directional cosine
    real(C_DOUBLE), intent(in)  :: mu      ! cosine of angle in lab or CM
    real(C_DOUBLE), intent(out) :: uvw(3)  ! rotated directional cosine
    real(C_DOUBLE), optional    :: phi     ! azimuthal angle

    real(C_DOUBLE) :: phi_   ! azimuthal angle
    real(C_DOUBLE) :: sinphi ! sine of azimuthal angle
    real(C_DOUBLE) :: cosphi ! cosine of azimuthal angle
    real(C_DOUBLE) :: a      ! sqrt(1 - mu^2)
    real(C_DOUBLE) :: b      ! sqrt(1 - w^2)
    real(C_DOUBLE) :: u0     ! original cosine in x direction
    real(C_DOUBLE) :: v0     ! original cosine in y direction
    real(C_DOUBLE) :: w0     ! original cosine in z direction

    ! Copy original directional cosines
    u0 = uvw0(1)
    v0 = uvw0(2)
    w0 = uvw0(3)

    ! Sample azimuthal angle in [0,2pi) if none provided
    if (present(phi)) then
      phi_ = phi
    else
      phi_ = TWO * PI * prn()
    end if

    ! Precompute factors to save flops
    sinphi = sin(phi_)
    cosphi = cos(phi_)
    a = sqrt(max(ZERO, ONE - mu*mu))
    b = sqrt(max(ZERO, ONE - w0*w0))

    ! Need to treat special case where sqrt(1 - w**2) is close to zero by
    ! expanding about the v component rather than the w component
    if (b > 1e-10) then
      uvw(1) = mu*u0 + a*(u0*w0*cosphi - v0*sinphi)/b
      uvw(2) = mu*v0 + a*(v0*w0*cosphi + u0*sinphi)/b
      uvw(3) = mu*w0 - a*b*cosphi
    else
      b = sqrt(ONE - v0*v0)
      uvw(1) = mu*u0 + a*(u0*v0*cosphi + w0*sinphi)/b
      uvw(2) = mu*v0 - a*b*cosphi
      uvw(3) = mu*w0 + a*(v0*w0*cosphi - u0*sinphi)/b
    end if

  end subroutine rotate_angle

!===============================================================================
! MAXWELL_SPECTRUM samples an energy from the Maxwell fission distribution based
! on a direct sampling scheme. The probability distribution function for a
! Maxwellian is given as p(x) = 2/(T*sqrt(pi))*sqrt(x/T)*exp(-x/T). This PDF can
! be sampled using rule C64 in the Monte Carlo Sampler LA-9721-MS.
!===============================================================================

  function maxwell_spectrum(T) result(E_out) bind(C)

    real(C_DOUBLE), intent(in)  :: T     ! tabulated function of incoming E
    real(C_DOUBLE)              :: E_out ! sampled energy

    real(C_DOUBLE) :: r1, r2, r3  ! random numbers
    real(C_DOUBLE) :: c           ! cosine of pi/2*r3

    r1 = prn()
    r2 = prn()
    r3 = prn()

    ! determine cosine of pi/2*r
    c = cos(PI/TWO*r3)

    ! determine outgoing energy
    E_out = -T*(log(r1) + log(r2)*c*c)

  end function maxwell_spectrum

!===============================================================================
! WATT_SPECTRUM samples the outgoing energy from a Watt energy-dependent fission
! spectrum. Although fitted parameters exist for many nuclides, generally the
! continuous tabular distributions (LAW 4) should be used in lieu of the Watt
! spectrum. This direct sampling scheme is an unpublished scheme based on the
! original Watt spectrum derivation (See F. Brown's MC lectures).
!===============================================================================

  function watt_spectrum(a, b) result(E_out) bind(C)

    real(C_DOUBLE), intent(in) :: a     ! Watt parameter a
    real(C_DOUBLE), intent(in) :: b     ! Watt parameter b
    real(C_DOUBLE)             :: E_out ! energy of emitted neutron

    real(C_DOUBLE) :: w ! sampled from Maxwellian

    w     = maxwell_spectrum(a)
    E_out = w + a*a*b/4. + (TWO*prn() - ONE)*sqrt(a*a*b*w)

  end function watt_spectrum

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

!===============================================================================
! BROADEN_WMP_POLYNOMIALS Doppler broadens the windowed multipole curvefit.  The
! curvefit is a polynomial of the form
! a/E + b/sqrt(E) + c + d sqrt(E) ...
!===============================================================================

  subroutine broaden_wmp_polynomials(E, dopp, n, factors) bind(C)
    real(C_DOUBLE), intent(in) :: E          ! Energy to evaluate at
    real(C_DOUBLE), intent(in) :: dopp       ! sqrt(atomic weight ratio / kT),
                                      !  kT given in eV.
    integer(C_INT), intent(in) :: n          ! number of components to polynomial
    real(C_DOUBLE), intent(out):: factors(n) ! output leading coefficient

    integer :: i

    real(8) :: sqrtE               ! sqrt(energy)
    real(8) :: beta                ! sqrt(atomic weight ratio * E / kT)
    real(8) :: half_inv_dopp2      ! 0.5 / dopp**2
    real(8) :: quarter_inv_dopp4   ! 0.25 / dopp**4
    real(8) :: erf_beta            ! error function of beta
    real(8) :: exp_m_beta2         ! exp(-beta**2)

    sqrtE = sqrt(E)
    beta = sqrtE * dopp
    half_inv_dopp2 = HALF / dopp**2
    quarter_inv_dopp4 = half_inv_dopp2**2

    if (beta > 6.0_8) then
      ! Save time, ERF(6) is 1 to machine precision.
      ! beta/sqrtpi*exp(-beta**2) is also approximately 1 machine epsilon.
      erf_beta = ONE
      exp_m_beta2 = ZERO
    else
      erf_beta = erf(beta)
      exp_m_beta2 = exp(-beta**2)
    end if

    ! Assume that, for sure, we'll use a second order (1/E, 1/V, const)
    ! fit, and no less.

    factors(1) = erf_beta / E
    factors(2) = ONE / sqrtE
    factors(3) = factors(1) * (half_inv_dopp2 + E) &
         + exp_m_beta2 / (beta * SQRT_PI)

    ! Perform recursive broadening of high order components
    do i = 1, n-3
      if (i /= 1) then
        factors(i+3) = -factors(i-1) * (i - ONE) * i * quarter_inv_dopp4 &
             + factors(i+1) * (E + (ONE + TWO * i) * half_inv_dopp2)
      else
        ! Although it's mathematically identical, factors(0) will contain
        ! nothing, and we don't want to have to worry about memory.
        factors(i+3) = factors(i+1)*(E + (ONE + TWO * i) * half_inv_dopp2)
      end if
    end do
  end subroutine broaden_wmp_polynomials

end module math
