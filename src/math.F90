module math

  use, intrinsic :: ISO_C_BINDING

  use constants
  use random_lcg, only: prn

  implicit none

!===============================================================================
! FADDEEVA_W evaluates the scaled complementary error function.  This
! interfaces with the MIT C library
!===============================================================================

  interface
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

  elemental function normal_percentile(p) result(z)

    real(8), intent(in) :: p ! probability level
    real(8)             :: z ! corresponding z-value

    real(8)            :: q
    real(8)            :: r
    real(8), parameter :: p_low  = 0.02425_8
    real(8), parameter :: a(6) = (/ &
         -3.969683028665376e1_8, 2.209460984245205e2_8, -2.759285104469687e2_8, &
         1.383577518672690e2_8, -3.066479806614716e1_8, 2.506628277459239e0_8 /)
    real(8), parameter :: b(5) = (/ &
         -5.447609879822406e1_8, 1.615858368580409e2_8, -1.556989798598866e2_8, &
         6.680131188771972e1_8, -1.328068155288572e1_8 /)
    real(8), parameter :: c(6) = (/ &
         -7.784894002430293e-3_8, -3.223964580411365e-1_8, -2.400758277161838_8, &
         -2.549732539343734_8, 4.374664141464968_8, 2.938163982698783_8 /)
    real(8), parameter :: d(4) = (/ &
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

  elemental function t_percentile(p, df) result(t)

    real(8), intent(in) :: p  ! probability level
    integer, intent(in) :: df ! degrees of freedom
    real(8)             :: t  ! corresponding t-value

    real(8)            :: n  ! degrees of freedom as a real(8)
    real(8)            :: k  ! n - 2
    real(8)            :: z  ! percentile of normal distribution
    real(8)            :: z2 ! z * z

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

  elemental function calc_pn(n,x) result(pnx)

    integer, intent(in) :: n   ! Legendre order requested
    real(8), intent(in) :: x   ! Independent variable the Legendre is to be
                               ! evaluated at; x must be in the domain [-1,1]
    real(8)             :: pnx ! The Legendre poly of order n evaluated at x

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
! CALC_RN calculates the n-th order spherical harmonics for a given angle
! (in terms of (u,v,w)).  All Rn,m values are provided (where -n<=m<=n)
!===============================================================================

  pure function calc_rn(n,uvw) result(rn)

    integer, intent(in) :: n      ! Order requested
    real(8), intent(in) :: uvw(3) ! Direction of travel, assumed to be on unit sphere
    real(8)             :: rn(2*n + 1)     ! The resultant R_n(uvw)

    real(8) :: phi, w ! Azimuthal and Cosine of Polar angles (from uvw)
    real(8) :: w2m1   ! (w^2 - 1), frequently used in these

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

  end function calc_rn

!===============================================================================
! CALC_ZN calculates the n-th order Zernike polynomial moment for a given angle
! (rho, theta) location in the unit disk.
!===============================================================================

  pure function calc_zn(n, rho, phi) result(zn)

    integer, intent(in) :: n           ! Order requested
    real(8), intent(in) :: rho         ! Radial location in the unit disk
    real(8), intent(in) :: phi         ! Theta (radians) location in the unit disk
    real(8)             :: zn(n + 1)   ! The resultant Z_n(uvw)

    ! n == radial degree
    ! m == azimuthal frequency

    select case(n)
    case(0)
      ! n = 0, m = 0
      zn(1) = ( ( 1.00 ) ) &
           * ( 1.000000000000 )
    case(1)
      ! n = 1, m = -1
      zn(1) = ( ( 1.00 ) * rho ) &
           * ( 2.000000000000 ) * sin(1.00 * phi)
      ! n = 1, m = 1
      zn(2) = ( ( 1.00 ) * rho ) &
           * ( 2.000000000000 ) * cos(1.00 * phi)
    case(2)
      ! n = 2, m = -2
      zn(1) = ( ( 1.00 ) *rho * rho ) &
           * ( 2.449489742783 ) * sin(2.00 * phi)
      ! n = 2, m = 0
      zn(2) = ( ( -1.00 ) + &
           ( 2.00 ) *rho * rho ) &
           * ( 1.732050807569 )
      ! n = 2, m = 2
      zn(3) = ( ( 1.00 ) *rho * rho ) &
           * ( 2.449489742783 ) * cos(2.00 * phi)
    case(3)
      ! n = 3, m = -3
      zn(1) = ( ( 1.00 ) *rho *rho * rho ) &
           * ( 2.828427124746 ) * sin(3.00 * phi)
      ! n = 3, m = -1
      zn(2) = ( ( -2.00 ) * rho + &
           ( 3.00 ) *rho *rho * rho ) &
           * ( 2.828427124746 ) * sin(3.00 * phi)
      ! n = 3, m = 1
      zn(3) = ( ( -2.00 ) * rho + &
           ( 3.00 ) *rho *rho * rho ) &
           * ( 2.828427124746 ) * cos(3.00 * phi)
      ! n = 3, m = 3
      zn(4) = ( ( 1.00 ) *rho *rho * rho ) &
           * ( 2.828427124746 ) * cos(3.00 * phi)
    case(4)
      ! n = 4, m = -4
      zn(1) = ( ( 1.00 ) *rho *rho *rho * rho ) &
           * ( 3.162277660168 ) * sin(4.00 * phi)
      ! n = 4, m = -2
      zn(2) = ( ( -3.00 ) *rho * rho + &
           ( 4.00 ) *rho *rho *rho * rho ) &
           * ( 3.162277660168 ) * sin(4.00 * phi)
      ! n = 4, m = 0
      zn(3) = ( ( 1.00 ) + &
           ( -6.00 ) *rho * rho + &
           ( 6.00 ) *rho *rho *rho * rho ) &
           * ( 2.236067977500 )
      ! n = 4, m = 2
      zn(4) = ( ( -3.00 ) *rho * rho + &
           ( 4.00 ) *rho *rho *rho * rho ) &
           * ( 3.162277660168 ) * cos(4.00 * phi)
      ! n = 4, m = 4
      zn(5) = ( ( 1.00 ) *rho *rho *rho * rho ) &
           * ( 3.162277660168 ) * cos(4.00 * phi)
    case(5)
      ! n = 5, m = -5
      zn(1) = ( ( 1.00 ) *rho *rho *rho *rho * rho ) &
           * ( 3.464101615138 ) * sin(5.00 * phi)
      ! n = 5, m = -3
      zn(2) = ( ( -4.00 ) *rho *rho * rho + &
           ( 5.00 ) *rho *rho *rho *rho * rho ) &
           * ( 3.464101615138 ) * sin(5.00 * phi)
      ! n = 5, m = -1
      zn(3) = ( ( 3.00 ) * rho + &
           ( -12.00 ) *rho *rho * rho + &
           ( 10.00 ) *rho *rho *rho *rho * rho ) &
           * ( 3.464101615138 ) * sin(5.00 * phi)
      ! n = 5, m = 1
      zn(4) = ( ( 3.00 ) * rho + &
           ( -12.00 ) *rho *rho * rho + &
           ( 10.00 ) *rho *rho *rho *rho * rho ) &
           * ( 3.464101615138 ) * cos(5.00 * phi)
      ! n = 5, m = 3
      zn(5) = ( ( -4.00 ) *rho *rho * rho + &
           ( 5.00 ) *rho *rho *rho *rho * rho ) &
           * ( 3.464101615138 ) * cos(5.00 * phi)
      ! n = 5, m = 5
      zn(6) = ( ( 1.00 ) *rho *rho *rho *rho * rho ) &
           * ( 3.464101615138 ) * cos(5.00 * phi)
    case(6)
      ! n = 6, m = -6
      zn(1) = ( ( 1.00 ) *rho *rho *rho *rho *rho * rho ) &
           * ( 3.741657386774 ) * sin(6.00 * phi)
      ! n = 6, m = -4
      zn(2) = ( ( -5.00 ) *rho *rho *rho * rho + &
           ( 6.00 ) *rho *rho *rho *rho *rho * rho ) &
           * ( 3.741657386774 ) * sin(6.00 * phi)
      ! n = 6, m = -2
      zn(3) = ( ( 6.00 ) *rho * rho + &
           ( -20.00 ) *rho *rho *rho * rho + &
           ( 15.00 ) *rho *rho *rho *rho *rho * rho ) &
           * ( 3.741657386774 ) * sin(6.00 * phi)
      ! n = 6, m = 0
      zn(4) = ( ( -1.00 ) + &
           ( 12.00 ) *rho * rho + &
           ( -30.00 ) *rho *rho *rho * rho + &
           ( 20.00 ) *rho *rho *rho *rho *rho * rho ) &
           * ( 2.645751311065 )
      ! n = 6, m = 2
      zn(5) = ( ( 6.00 ) *rho * rho + &
           ( -20.00 ) *rho *rho *rho * rho + &
           ( 15.00 ) *rho *rho *rho *rho *rho * rho ) &
           * ( 3.741657386774 ) * cos(6.00 * phi)
      ! n = 6, m = 4
      zn(6) = ( ( -5.00 ) *rho *rho *rho * rho + &
           ( 6.00 ) *rho *rho *rho *rho *rho * rho ) &
           * ( 3.741657386774 ) * cos(6.00 * phi)
      ! n = 6, m = 6
      zn(7) = ( ( 1.00 ) *rho *rho *rho *rho *rho * rho ) &
           * ( 3.741657386774 ) * cos(6.00 * phi)
    case(7)
      ! n = 7, m = -7
      zn(1) = ( ( 1.00 ) *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.000000000000 ) * sin(7.00 * phi)
      ! n = 7, m = -5
      zn(2) = ( ( -6.00 ) *rho *rho *rho *rho * rho + &
           ( 7.00 ) *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.000000000000 ) * sin(7.00 * phi)
      ! n = 7, m = -3
      zn(3) = ( ( 10.00 ) *rho *rho * rho + &
           ( -30.00 ) *rho *rho *rho *rho * rho + &
           ( 21.00 ) *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.000000000000 ) * sin(7.00 * phi)
      ! n = 7, m = -1
      zn(4) = ( ( -4.00 ) * rho + &
           ( 30.00 ) *rho *rho * rho + &
           ( -60.00 ) *rho *rho *rho *rho * rho + &
           ( 35.00 ) *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.000000000000 ) * sin(7.00 * phi)
      ! n = 7, m = 1
      zn(5) = ( ( -4.00 ) * rho + &
           ( 30.00 ) *rho *rho * rho + &
           ( -60.00 ) *rho *rho *rho *rho * rho + &
           ( 35.00 ) *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.000000000000 ) * cos(7.00 * phi)
      ! n = 7, m = 3
      zn(6) = ( ( 10.00 ) *rho *rho * rho + &
           ( -30.00 ) *rho *rho *rho *rho * rho + &
           ( 21.00 ) *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.000000000000 ) * cos(7.00 * phi)
      ! n = 7, m = 5
      zn(7) = ( ( -6.00 ) *rho *rho *rho *rho * rho + &
           ( 7.00 ) *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.000000000000 ) * cos(7.00 * phi)
      ! n = 7, m = 7
      zn(8) = ( ( 1.00 ) *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.000000000000 ) * cos(7.00 * phi)
    case(8)
      ! n = 8, m = -8
      zn(1) = ( ( 1.00 ) *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.242640687119 ) * sin(8.00 * phi)
      ! n = 8, m = -6
      zn(2) = ( ( -7.00 ) *rho *rho *rho *rho *rho * rho + &
           ( 8.00 ) *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.242640687119 ) * sin(8.00 * phi)
      ! n = 8, m = -4
      zn(3) = ( ( 15.00 ) *rho *rho *rho * rho + &
           ( -42.00 ) *rho *rho *rho *rho *rho * rho + &
           ( 28.00 ) *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.242640687119 ) * sin(8.00 * phi)
      ! n = 8, m = -2
      zn(4) = ( ( -10.00 ) *rho * rho + &
           ( 60.00 ) *rho *rho *rho * rho + &
           ( -105.00 ) *rho *rho *rho *rho *rho * rho + &
           ( 56.00 ) *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.242640687119 ) * sin(8.00 * phi)
      ! n = 8, m = 0
      zn(5) = ( ( 1.00 ) + &
           ( -20.00 ) *rho * rho + &
           ( 90.00 ) *rho *rho *rho * rho + &
           ( -140.00 ) *rho *rho *rho *rho *rho * rho + &
           ( 70.00 ) *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 3.000000000000 )
      ! n = 8, m = 2
      zn(6) = ( ( -10.00 ) *rho * rho + &
           ( 60.00 ) *rho *rho *rho * rho + &
           ( -105.00 ) *rho *rho *rho *rho *rho * rho + &
           ( 56.00 ) *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.242640687119 ) * cos(8.00 * phi)
      ! n = 8, m = 4
      zn(7) = ( ( 15.00 ) *rho *rho *rho * rho + &
           ( -42.00 ) *rho *rho *rho *rho *rho * rho + &
           ( 28.00 ) *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.242640687119 ) * cos(8.00 * phi)
      ! n = 8, m = 6
      zn(8) = ( ( -7.00 ) *rho *rho *rho *rho *rho * rho + &
           ( 8.00 ) *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.242640687119 ) * cos(8.00 * phi)
      ! n = 8, m = 8
      zn(9) = ( ( 1.00 ) *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.242640687119 ) * cos(8.00 * phi)
    case(9)
      ! n = 9, m = -9
      zn(1) = ( ( 1.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.472135955000 ) * sin(9.00 * phi)
      ! n = 9, m = -7
      zn(2) = ( ( -8.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( 9.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.472135955000 ) * sin(9.00 * phi)
      ! n = 9, m = -5
      zn(3) = ( ( 21.00 ) *rho *rho *rho *rho * rho + &
           ( -56.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( 36.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.472135955000 ) * sin(9.00 * phi)
      ! n = 9, m = -3
      zn(4) = ( ( -20.00 ) *rho *rho * rho + &
           ( 105.00 ) *rho *rho *rho *rho * rho + &
           ( -168.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( 84.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.472135955000 ) * sin(9.00 * phi)
      ! n = 9, m = -1
      zn(5) = ( ( 5.00 ) * rho + &
           ( -60.00 ) *rho *rho * rho + &
           ( 210.00 ) *rho *rho *rho *rho * rho + &
           ( -280.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( 126.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.472135955000 ) * sin(9.00 * phi)
      ! n = 9, m = 1
      zn(6) = ( ( 5.00 ) * rho + &
           ( -60.00 ) *rho *rho * rho + &
           ( 210.00 ) *rho *rho *rho *rho * rho + &
           ( -280.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( 126.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.472135955000 ) * cos(9.00 * phi)
      ! n = 9, m = 3
      zn(7) = ( ( -20.00 ) *rho *rho * rho + &
           ( 105.00 ) *rho *rho *rho *rho * rho + &
           ( -168.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( 84.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.472135955000 ) * cos(9.00 * phi)
      ! n = 9, m = 5
      zn(8) = ( ( 21.00 ) *rho *rho *rho *rho * rho + &
           ( -56.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( 36.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.472135955000 ) * cos(9.00 * phi)
      ! n = 9, m = 7
      zn(9) = ( ( -8.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( 9.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.472135955000 ) * cos(9.00 * phi)
      ! n = 9, m = 9
      zn(10) = ( ( 1.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.472135955000 ) * cos(9.00 * phi)
    case(10)
      ! n = 10, m = -10
      zn(1) = ( ( 1.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.690415759823 ) * sin(10.00 * phi)
      ! n = 10, m = -8
      zn(2) = ( ( -9.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 10.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.690415759823 ) * sin(10.00 * phi)
      ! n = 10, m = -6
      zn(3) = ( ( 28.00 ) *rho *rho *rho *rho *rho * rho + &
           ( -72.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 45.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.690415759823 ) * sin(10.00 * phi)
      ! n = 10, m = -4
      zn(4) = ( ( -35.00 ) *rho *rho *rho * rho + &
           ( 168.00 ) *rho *rho *rho *rho *rho * rho + &
           ( -252.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 120.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.690415759823 ) * sin(10.00 * phi)
      ! n = 10, m = -2
      zn(5) = ( ( 15.00 ) *rho * rho + &
           ( -140.00 ) *rho *rho *rho * rho + &
           ( 420.00 ) *rho *rho *rho *rho *rho * rho + &
           ( -504.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 210.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.690415759823 ) * sin(10.00 * phi)
      ! n = 10, m = 0
      zn(6) = ( ( -1.00 ) + &
           ( 30.00 ) *rho * rho + &
           ( -210.00 ) *rho *rho *rho * rho + &
           ( 560.00 ) *rho *rho *rho *rho *rho * rho + &
           ( -630.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 252.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 3.316624790355 )
      ! n = 10, m = 2
      zn(7) = ( ( 15.00 ) *rho * rho + &
           ( -140.00 ) *rho *rho *rho * rho + &
           ( 420.00 ) *rho *rho *rho *rho *rho * rho + &
           ( -504.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 210.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.690415759823 ) * cos(10.00 * phi)
      ! n = 10, m = 4
      zn(8) = ( ( -35.00 ) *rho *rho *rho * rho + &
           ( 168.00 ) *rho *rho *rho *rho *rho * rho + &
           ( -252.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 120.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.690415759823 ) * cos(10.00 * phi)
      ! n = 10, m = 6
      zn(9) = ( ( 28.00 ) *rho *rho *rho *rho *rho * rho + &
           ( -72.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 45.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.690415759823 ) * cos(10.00 * phi)
      ! n = 10, m = 8
      zn(10) = ( ( -9.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 10.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.690415759823 ) * cos(10.00 * phi)
      ! n = 10, m = 10
      zn(11) = ( ( 1.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.690415759823 ) * cos(10.00 * phi)
    case(11)
      ! n = 11, m = -11
      zn(1) = ( ( 1.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.898979485566 ) * sin(11.00 * phi)
      ! n = 11, m = -9
      zn(2) = ( ( -10.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 11.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.898979485566 ) * sin(11.00 * phi)
      ! n = 11, m = -7
      zn(3) = ( ( 36.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( -90.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 55.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.898979485566 ) * sin(11.00 * phi)
      ! n = 11, m = -5
      zn(4) = ( ( -56.00 ) *rho *rho *rho *rho * rho + &
           ( 252.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( -360.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 165.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.898979485566 ) * sin(11.00 * phi)
      ! n = 11, m = -3
      zn(5) = ( ( 35.00 ) *rho *rho * rho + &
           ( -280.00 ) *rho *rho *rho *rho * rho + &
           ( 756.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( -840.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 330.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.898979485566 ) * sin(11.00 * phi)
      ! n = 11, m = -1
      zn(6) = ( ( -6.00 ) * rho + &
           ( 105.00 ) *rho *rho * rho + &
           ( -560.00 ) *rho *rho *rho *rho * rho + &
           ( 1260.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( -1260.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 462.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.898979485566 ) * sin(11.00 * phi)
      ! n = 11, m = 1
      zn(7) = ( ( -6.00 ) * rho + &
           ( 105.00 ) *rho *rho * rho + &
           ( -560.00 ) *rho *rho *rho *rho * rho + &
           ( 1260.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( -1260.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 462.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.898979485566 ) * cos(11.00 * phi)
      ! n = 11, m = 3
      zn(8) = ( ( 35.00 ) *rho *rho * rho + &
           ( -280.00 ) *rho *rho *rho *rho * rho + &
           ( 756.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( -840.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 330.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.898979485566 ) * cos(11.00 * phi)
      ! n = 11, m = 5
      zn(9) = ( ( -56.00 ) *rho *rho *rho *rho * rho + &
           ( 252.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( -360.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 165.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.898979485566 ) * cos(11.00 * phi)
      ! n = 11, m = 7
      zn(10) = ( ( 36.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( -90.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 55.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.898979485566 ) * cos(11.00 * phi)
      ! n = 11, m = 9
      zn(11) = ( ( -10.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 11.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.898979485566 ) * cos(11.00 * phi)
      ! n = 11, m = 11
      zn(12) = ( ( 1.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.898979485566 ) * cos(11.00 * phi)
    case(12)
      ! n = 12, m = -12
      zn(1) = ( ( 1.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.099019513593 ) * sin(12.00 * phi)
      ! n = 12, m = -10
      zn(2) = ( ( -11.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 12.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.099019513593 ) * sin(12.00 * phi)
      ! n = 12, m = -8
      zn(3) = ( ( 45.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -110.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 66.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.099019513593 ) * sin(12.00 * phi)
      ! n = 12, m = -6
      zn(4) = ( ( -84.00 ) *rho *rho *rho *rho *rho * rho + &
           ( 360.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -495.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 220.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.099019513593 ) * sin(12.00 * phi)
      ! n = 12, m = -4
      zn(5) = ( ( 70.00 ) *rho *rho *rho * rho + &
           ( -504.00 ) *rho *rho *rho *rho *rho * rho + &
           ( 1260.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -1320.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 495.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.099019513593 ) * sin(12.00 * phi)
      ! n = 12, m = -2
      zn(6) = ( ( -21.00 ) *rho * rho + &
           ( 280.00 ) *rho *rho *rho * rho + &
           ( -1260.00 ) *rho *rho *rho *rho *rho * rho + &
           ( 2520.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -2310.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 792.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.099019513593 ) * sin(12.00 * phi)
      ! n = 12, m = 0
      zn(7) = ( ( 1.00 ) + &
           ( -42.00 ) *rho * rho + &
           ( 420.00 ) *rho *rho *rho * rho + &
           ( -1680.00 ) *rho *rho *rho *rho *rho * rho + &
           ( 3150.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -2772.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 924.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 3.605551275464 )
      ! n = 12, m = 2
      zn(8) = ( ( -21.00 ) *rho * rho + &
           ( 280.00 ) *rho *rho *rho * rho + &
           ( -1260.00 ) *rho *rho *rho *rho *rho * rho + &
           ( 2520.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -2310.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 792.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.099019513593 ) * cos(12.00 * phi)
      ! n = 12, m = 4
      zn(9) = ( ( 70.00 ) *rho *rho *rho * rho + &
           ( -504.00 ) *rho *rho *rho *rho *rho * rho + &
           ( 1260.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -1320.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 495.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.099019513593 ) * cos(12.00 * phi)
      ! n = 12, m = 6
      zn(10) = ( ( -84.00 ) *rho *rho *rho *rho *rho * rho + &
           ( 360.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -495.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 220.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.099019513593 ) * cos(12.00 * phi)
      ! n = 12, m = 8
      zn(11) = ( ( 45.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -110.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 66.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.099019513593 ) * cos(12.00 * phi)
      ! n = 12, m = 10
      zn(12) = ( ( -11.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 12.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.099019513593 ) * cos(12.00 * phi)
      ! n = 12, m = 12
      zn(13) = ( ( 1.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.099019513593 ) * cos(12.00 * phi)
    case(13)
      ! n = 13, m = -13
      zn(1) = ( ( 1.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.291502622129 ) * sin(13.00 * phi)
      ! n = 13, m = -11
      zn(2) = ( ( -12.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 13.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.291502622129 ) * sin(13.00 * phi)
      ! n = 13, m = -9
      zn(3) = ( ( 55.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -132.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 78.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.291502622129 ) * sin(13.00 * phi)
      ! n = 13, m = -7
      zn(4) = ( ( -120.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( 495.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -660.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 286.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.291502622129 ) * sin(13.00 * phi)
      ! n = 13, m = -5
      zn(5) = ( ( 126.00 ) *rho *rho *rho *rho * rho + &
           ( -840.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( 1980.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -1980.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 715.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.291502622129 ) * sin(13.00 * phi)
      ! n = 13, m = -3
      zn(6) = ( ( -56.00 ) *rho *rho * rho + &
           ( 630.00 ) *rho *rho *rho *rho * rho + &
           ( -2520.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( 4620.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -3960.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 1287.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.291502622129 ) * sin(13.00 * phi)
      ! n = 13, m = -1
      zn(7) = ( ( 7.00 ) * rho + &
           ( -168.00 ) *rho *rho * rho + &
           ( 1260.00 ) *rho *rho *rho *rho * rho + &
           ( -4200.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( 6930.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -5544.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 1716.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.291502622129 ) * sin(13.00 * phi)
      ! n = 13, m = 1
      zn(8) = ( ( 7.00 ) * rho + &
           ( -168.00 ) *rho *rho * rho + &
           ( 1260.00 ) *rho *rho *rho *rho * rho + &
           ( -4200.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( 6930.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -5544.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 1716.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.291502622129 ) * cos(13.00 * phi)
      ! n = 13, m = 3
      zn(9) = ( ( -56.00 ) *rho *rho * rho + &
           ( 630.00 ) *rho *rho *rho *rho * rho + &
           ( -2520.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( 4620.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -3960.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 1287.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.291502622129 ) * cos(13.00 * phi)
      ! n = 13, m = 5
      zn(10) = ( ( 126.00 ) *rho *rho *rho *rho * rho + &
           ( -840.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( 1980.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -1980.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 715.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.291502622129 ) * cos(13.00 * phi)
      ! n = 13, m = 7
      zn(11) = ( ( -120.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( 495.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -660.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 286.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.291502622129 ) * cos(13.00 * phi)
      ! n = 13, m = 9
      zn(12) = ( ( 55.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -132.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 78.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.291502622129 ) * cos(13.00 * phi)
      ! n = 13, m = 11
      zn(13) = ( ( -12.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 13.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.291502622129 ) * cos(13.00 * phi)
      ! n = 13, m = 13
      zn(14) = ( ( 1.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.291502622129 ) * cos(13.00 * phi)
    case(14)
      ! n = 14, m = -14
      zn(1) = ( ( 1.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.477225575052 ) * sin(14.00 * phi)
      ! n = 14, m = -12
      zn(2) = ( ( -13.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 14.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.477225575052 ) * sin(14.00 * phi)
      ! n = 14, m = -10
      zn(3) = ( ( 66.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -156.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 91.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.477225575052 ) * sin(14.00 * phi)
      ! n = 14, m = -8
      zn(4) = ( ( -165.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 660.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -858.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 364.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.477225575052 ) * sin(14.00 * phi)
      ! n = 14, m = -6
      zn(5) = ( ( 210.00 ) *rho *rho *rho *rho *rho * rho + &
           ( -1320.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 2970.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -2860.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 1001.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.477225575052 ) * sin(14.00 * phi)
      ! n = 14, m = -4
      zn(6) = ( ( -126.00 ) *rho *rho *rho * rho + &
           ( 1260.00 ) *rho *rho *rho *rho *rho * rho + &
           ( -4620.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 7920.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -6435.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 2002.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.477225575052 ) * sin(14.00 * phi)
      ! n = 14, m = -2
      zn(7) = ( ( 28.00 ) *rho * rho + &
           ( -504.00 ) *rho *rho *rho * rho + &
           ( 3150.00 ) *rho *rho *rho *rho *rho * rho + &
           ( -9240.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 13860.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -10296.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 3003.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.477225575052 ) * sin(14.00 * phi)
      ! n = 14, m = 0
      zn(8) = ( ( -1.00 ) + &
           ( 56.00 ) *rho * rho + &
           ( -756.00 ) *rho *rho *rho * rho + &
           ( 4200.00 ) *rho *rho *rho *rho *rho * rho + &
           ( -11550.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 16632.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -12012.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 3432.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 3.872983346207 )
      ! n = 14, m = 2
      zn(9) = ( ( 28.00 ) *rho * rho + &
           ( -504.00 ) *rho *rho *rho * rho + &
           ( 3150.00 ) *rho *rho *rho *rho *rho * rho + &
           ( -9240.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 13860.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -10296.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 3003.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.477225575052 ) * cos(14.00 * phi)
      ! n = 14, m = 4
      zn(10) = ( ( -126.00 ) *rho *rho *rho * rho + &
           ( 1260.00 ) *rho *rho *rho *rho *rho * rho + &
           ( -4620.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 7920.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -6435.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 2002.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.477225575052 ) * cos(14.00 * phi)
      ! n = 14, m = 6
      zn(11) = ( ( 210.00 ) *rho *rho *rho *rho *rho * rho + &
           ( -1320.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 2970.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -2860.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 1001.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.477225575052 ) * cos(14.00 * phi)
      ! n = 14, m = 8
      zn(12) = ( ( -165.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 660.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -858.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 364.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.477225575052 ) * cos(14.00 * phi)
      ! n = 14, m = 10
      zn(13) = ( ( 66.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -156.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 91.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.477225575052 ) * cos(14.00 * phi)
      ! n = 14, m = 12
      zn(14) = ( ( -13.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 14.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.477225575052 ) * cos(14.00 * phi)
      ! n = 14, m = 14
      zn(15) = ( ( 1.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.477225575052 ) * cos(14.00 * phi)
    case(15)
      ! n = 15, m = -15
      zn(1) = ( ( 1.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.656854249492 ) * sin(15.00 * phi)
      ! n = 15, m = -13
      zn(2) = ( ( -14.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 15.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.656854249492 ) * sin(15.00 * phi)
      ! n = 15, m = -11
      zn(3) = ( ( 78.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -182.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 105.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.656854249492 ) * sin(15.00 * phi)
      ! n = 15, m = -9
      zn(4) = ( ( -220.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 858.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -1092.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 455.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.656854249492 ) * sin(15.00 * phi)
      ! n = 15, m = -7
      zn(5) = ( ( 330.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( -1980.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 4290.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -4004.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 1365.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.656854249492 ) * sin(15.00 * phi)
      ! n = 15, m = -5
      zn(6) = ( ( -252.00 ) *rho *rho *rho *rho * rho + &
           ( 2310.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( -7920.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 12870.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -10010.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 3003.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.656854249492 ) * sin(15.00 * phi)
      ! n = 15, m = -3
      zn(7) = ( ( 84.00 ) *rho *rho * rho + &
           ( -1260.00 ) *rho *rho *rho *rho * rho + &
           ( 6930.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( -18480.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 25740.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -18018.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 5005.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.656854249492 ) * sin(15.00 * phi)
      ! n = 15, m = -1
      zn(8) = ( ( -8.00 ) * rho + &
           ( 252.00 ) *rho *rho * rho + &
           ( -2520.00 ) *rho *rho *rho *rho * rho + &
           ( 11550.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( -27720.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 36036.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -24024.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 6435.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.656854249492 ) * sin(15.00 * phi)
      ! n = 15, m = 1
      zn(9) = ( ( -8.00 ) * rho + &
           ( 252.00 ) *rho *rho * rho + &
           ( -2520.00 ) *rho *rho *rho *rho * rho + &
           ( 11550.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( -27720.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 36036.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -24024.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 6435.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.656854249492 ) * cos(15.00 * phi)
      ! n = 15, m = 3
      zn(10) = ( ( 84.00 ) *rho *rho * rho + &
           ( -1260.00 ) *rho *rho *rho *rho * rho + &
           ( 6930.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( -18480.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 25740.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -18018.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 5005.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.656854249492 ) * cos(15.00 * phi)
      ! n = 15, m = 5
      zn(11) = ( ( -252.00 ) *rho *rho *rho *rho * rho + &
           ( 2310.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( -7920.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 12870.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -10010.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 3003.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.656854249492 ) * cos(15.00 * phi)
      ! n = 15, m = 7
      zn(12) = ( ( 330.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( -1980.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 4290.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -4004.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 1365.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.656854249492 ) * cos(15.00 * phi)
      ! n = 15, m = 9
      zn(13) = ( ( -220.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 858.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -1092.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 455.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.656854249492 ) * cos(15.00 * phi)
      ! n = 15, m = 11
      zn(14) = ( ( 78.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -182.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 105.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.656854249492 ) * cos(15.00 * phi)
      ! n = 15, m = 13
      zn(15) = ( ( -14.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 15.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.656854249492 ) * cos(15.00 * phi)
      ! n = 15, m = 15
      zn(16) = ( ( 1.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.656854249492 ) * cos(15.00 * phi)
    case(16)
      ! n = 16, m = -16
      zn(1) = ( ( 1.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.830951894845 ) * sin(16.00 * phi)
      ! n = 16, m = -14
      zn(2) = ( ( -15.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 16.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.830951894845 ) * sin(16.00 * phi)
      ! n = 16, m = -12
      zn(3) = ( ( 91.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -210.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 120.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.830951894845 ) * sin(16.00 * phi)
      ! n = 16, m = -10
      zn(4) = ( ( -286.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 1092.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -1365.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 560.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.830951894845 ) * sin(16.00 * phi)
      ! n = 16, m = -8
      zn(5) = ( ( 495.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -2860.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 6006.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -5460.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 1820.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.830951894845 ) * sin(16.00 * phi)
      ! n = 16, m = -6
      zn(6) = ( ( -462.00 ) *rho *rho *rho *rho *rho * rho + &
           ( 3960.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -12870.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 20020.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -15015.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 4368.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.830951894845 ) * sin(16.00 * phi)
      ! n = 16, m = -4
      zn(7) = ( ( 210.00 ) *rho *rho *rho * rho + &
           ( -2772.00 ) *rho *rho *rho *rho *rho * rho + &
           ( 13860.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -34320.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 45045.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -30030.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 8008.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.830951894845 ) * sin(16.00 * phi)
      ! n = 16, m = -2
      zn(8) = ( ( -36.00 ) *rho * rho + &
           ( 840.00 ) *rho *rho *rho * rho + &
           ( -6930.00 ) *rho *rho *rho *rho *rho * rho + &
           ( 27720.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -60060.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 72072.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -45045.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 11440.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.830951894845 ) * sin(16.00 * phi)
      ! n = 16, m = 0
      zn(9) = ( ( 1.00 ) + &
           ( -72.00 ) *rho * rho + &
           ( 1260.00 ) *rho *rho *rho * rho + &
           ( -9240.00 ) *rho *rho *rho *rho *rho * rho + &
           ( 34650.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -72072.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 84084.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -51480.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 12870.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.123105625618 )
      ! n = 16, m = 2
      zn(10) = ( ( -36.00 ) *rho * rho + &
           ( 840.00 ) *rho *rho *rho * rho + &
           ( -6930.00 ) *rho *rho *rho *rho *rho * rho + &
           ( 27720.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -60060.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 72072.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -45045.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 11440.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.830951894845 ) * cos(16.00 * phi)
      ! n = 16, m = 4
      zn(11) = ( ( 210.00 ) *rho *rho *rho * rho + &
           ( -2772.00 ) *rho *rho *rho *rho *rho * rho + &
           ( 13860.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -34320.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 45045.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -30030.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 8008.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.830951894845 ) * cos(16.00 * phi)
      ! n = 16, m = 6
      zn(12) = ( ( -462.00 ) *rho *rho *rho *rho *rho * rho + &
           ( 3960.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -12870.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 20020.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -15015.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 4368.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.830951894845 ) * cos(16.00 * phi)
      ! n = 16, m = 8
      zn(13) = ( ( 495.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -2860.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 6006.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -5460.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 1820.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.830951894845 ) * cos(16.00 * phi)
      ! n = 16, m = 10
      zn(14) = ( ( -286.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 1092.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -1365.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 560.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.830951894845 ) * cos(16.00 * phi)
      ! n = 16, m = 12
      zn(15) = ( ( 91.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -210.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 120.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.830951894845 ) * cos(16.00 * phi)
      ! n = 16, m = 14
      zn(16) = ( ( -15.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 16.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.830951894845 ) * cos(16.00 * phi)
      ! n = 16, m = 16
      zn(17) = ( ( 1.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.830951894845 ) * cos(16.00 * phi)
    case(17)
      ! n = 17, m = -17
      zn(1) = ( ( 1.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 6.000000000000 ) * sin(17.00 * phi)
      ! n = 17, m = -15
      zn(2) = ( ( -16.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 17.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 6.000000000000 ) * sin(17.00 * phi)
      ! n = 17, m = -13
      zn(3) = ( ( 105.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -240.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 136.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 6.000000000000 ) * sin(17.00 * phi)
      ! n = 17, m = -11
      zn(4) = ( ( -364.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 1365.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -1680.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 680.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 6.000000000000 ) * sin(17.00 * phi)
      ! n = 17, m = -9
      zn(5) = ( ( 715.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -4004.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 8190.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -7280.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 2380.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 6.000000000000 ) * sin(17.00 * phi)
      ! n = 17, m = -7
      zn(6) = ( ( -792.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( 6435.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -20020.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 30030.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -21840.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 6188.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 6.000000000000 ) * sin(17.00 * phi)
      ! n = 17, m = -5
      zn(7) = ( ( 462.00 ) *rho *rho *rho *rho * rho + &
           ( -5544.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( 25740.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -60060.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 75075.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -48048.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 12376.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 6.000000000000 ) * sin(17.00 * phi)
      ! n = 17, m = -3
      zn(8) = ( ( -120.00 ) *rho *rho * rho + &
           ( 2310.00 ) *rho *rho *rho *rho * rho + &
           ( -16632.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( 60060.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -120120.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 135135.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -80080.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 19448.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 6.000000000000 ) * sin(17.00 * phi)
      ! n = 17, m = -1
      zn(9) = ( ( 9.00 ) * rho + &
           ( -360.00 ) *rho *rho * rho + &
           ( 4620.00 ) *rho *rho *rho *rho * rho + &
           ( -27720.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( 90090.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -168168.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 180180.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -102960.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 24310.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 6.000000000000 ) * sin(17.00 * phi)
      ! n = 17, m = 1
      zn(10) = ( ( 9.00 ) * rho + &
           ( -360.00 ) *rho *rho * rho + &
           ( 4620.00 ) *rho *rho *rho *rho * rho + &
           ( -27720.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( 90090.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -168168.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 180180.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -102960.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 24310.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 6.000000000000 ) * cos(17.00 * phi)
      ! n = 17, m = 3
      zn(11) = ( ( -120.00 ) *rho *rho * rho + &
           ( 2310.00 ) *rho *rho *rho *rho * rho + &
           ( -16632.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( 60060.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -120120.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 135135.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -80080.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 19448.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 6.000000000000 ) * cos(17.00 * phi)
      ! n = 17, m = 5
      zn(12) = ( ( 462.00 ) *rho *rho *rho *rho * rho + &
           ( -5544.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( 25740.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -60060.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 75075.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -48048.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 12376.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 6.000000000000 ) * cos(17.00 * phi)
      ! n = 17, m = 7
      zn(13) = ( ( -792.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( 6435.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -20020.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 30030.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -21840.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 6188.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 6.000000000000 ) * cos(17.00 * phi)
      ! n = 17, m = 9
      zn(14) = ( ( 715.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -4004.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 8190.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -7280.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 2380.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 6.000000000000 ) * cos(17.00 * phi)
      ! n = 17, m = 11
      zn(15) = ( ( -364.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 1365.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -1680.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 680.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 6.000000000000 ) * cos(17.00 * phi)
      ! n = 17, m = 13
      zn(16) = ( ( 105.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -240.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 136.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 6.000000000000 ) * cos(17.00 * phi)
      ! n = 17, m = 15
      zn(17) = ( ( -16.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 17.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 6.000000000000 ) * cos(17.00 * phi)
      ! n = 17, m = 17
      zn(18) = ( ( 1.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 6.000000000000 ) * cos(17.00 * phi)
    case(18)
      ! n = 18, m = -18
      zn(1) = ( ( 1.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 6.164414002969 ) * sin(18.00 * phi)
      ! n = 18, m = -16
      zn(2) = ( ( -17.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 18.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 6.164414002969 ) * sin(18.00 * phi)
      ! n = 18, m = -14
      zn(3) = ( ( 120.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -272.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 153.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 6.164414002969 ) * sin(18.00 * phi)
      ! n = 18, m = -12
      zn(4) = ( ( -455.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 1680.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -2040.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 816.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 6.164414002969 ) * sin(18.00 * phi)
      ! n = 18, m = -10
      zn(5) = ( ( 1001.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -5460.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 10920.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -9520.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 3060.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 6.164414002969 ) * sin(18.00 * phi)
      ! n = 18, m = -8
      zn(6) = ( ( -1287.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 10010.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -30030.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 43680.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -30940.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 8568.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 6.164414002969 ) * sin(18.00 * phi)
      ! n = 18, m = -6
      zn(7) = ( ( 924.00 ) *rho *rho *rho *rho *rho * rho + &
           ( -10296.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 45045.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -100100.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 120120.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -74256.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 18564.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 6.164414002969 ) * sin(18.00 * phi)
      ! n = 18, m = -4
      zn(8) = ( ( -330.00 ) *rho *rho *rho * rho + &
           ( 5544.00 ) *rho *rho *rho *rho *rho * rho + &
           ( -36036.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 120120.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -225225.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 240240.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -136136.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 31824.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 6.164414002969 ) * sin(18.00 * phi)
      ! n = 18, m = -2
      zn(9) = ( ( 45.00 ) *rho * rho + &
           ( -1320.00 ) *rho *rho *rho * rho + &
           ( 13860.00 ) *rho *rho *rho *rho *rho * rho + &
           ( -72072.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 210210.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -360360.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 360360.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -194480.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 43758.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 6.164414002969 ) * sin(18.00 * phi)
      ! n = 18, m = 0
      zn(10) = ( ( -1.00 ) + &
           ( 90.00 ) *rho * rho + &
           ( -1980.00 ) *rho *rho *rho * rho + &
           ( 18480.00 ) *rho *rho *rho *rho *rho * rho + &
           ( -90090.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 252252.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -420420.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 411840.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -218790.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 48620.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.358898943541 )
      ! n = 18, m = 2
      zn(11) = ( ( 45.00 ) *rho * rho + &
           ( -1320.00 ) *rho *rho *rho * rho + &
           ( 13860.00 ) *rho *rho *rho *rho *rho * rho + &
           ( -72072.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 210210.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -360360.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 360360.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -194480.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 43758.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 6.164414002969 ) * cos(18.00 * phi)
      ! n = 18, m = 4
      zn(12) = ( ( -330.00 ) *rho *rho *rho * rho + &
           ( 5544.00 ) *rho *rho *rho *rho *rho * rho + &
           ( -36036.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 120120.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -225225.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 240240.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -136136.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 31824.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 6.164414002969 ) * cos(18.00 * phi)
      ! n = 18, m = 6
      zn(13) = ( ( 924.00 ) *rho *rho *rho *rho *rho * rho + &
           ( -10296.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 45045.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -100100.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 120120.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -74256.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 18564.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 6.164414002969 ) * cos(18.00 * phi)
      ! n = 18, m = 8
      zn(14) = ( ( -1287.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 10010.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -30030.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 43680.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -30940.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 8568.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 6.164414002969 ) * cos(18.00 * phi)
      ! n = 18, m = 10
      zn(15) = ( ( 1001.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -5460.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 10920.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -9520.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 3060.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 6.164414002969 ) * cos(18.00 * phi)
      ! n = 18, m = 12
      zn(16) = ( ( -455.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 1680.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -2040.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 816.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 6.164414002969 ) * cos(18.00 * phi)
      ! n = 18, m = 14
      zn(17) = ( ( 120.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -272.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 153.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 6.164414002969 ) * cos(18.00 * phi)
      ! n = 18, m = 16
      zn(18) = ( ( -17.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 18.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 6.164414002969 ) * cos(18.00 * phi)
      ! n = 18, m = 18
      zn(19) = ( ( 1.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 6.164414002969 ) * cos(18.00 * phi)
    case default
      zn = ONE
    end select

  end function calc_zn

!===============================================================================
! EXPAND_HARMONIC expands a given series of real spherical harmonics
!===============================================================================

  pure function expand_harmonic(data, order, uvw) result(val)
    real(8), intent(in) :: data(:)
    integer, intent(in) :: order
    real(8), intent(in) :: uvw(3)
    real(8)             :: val

    integer :: l, lm_lo, lm_hi

    val = data(1)
    lm_lo = 2
    lm_hi = 4
    do l = 1, order - 1
      val = val + sqrt(TWO * real(l,8) + ONE) * &
           dot_product(calc_rn(l,uvw), data(lm_lo:lm_hi))
      lm_lo = lm_hi + 1
      lm_hi = lm_lo + 2 * (l + 1)
    end do

  end function expand_harmonic

!===============================================================================
! EVALUATE_LEGENDRE Find the value of f(x) given a set of Legendre coefficients
! and the value of x
!===============================================================================

  pure function evaluate_legendre(data, x) result(val)
    real(8), intent(in) :: data(:)
    real(8), intent(in) :: x
    real(8)             :: val

    integer :: l

    val =  HALF * data(1)
    do l = 1, size(data) - 1
      val = val + (real(l,8) +  HALF) * data(l + 1) * calc_pn(l,x)
    end do

  end function evaluate_legendre

!===============================================================================
! ROTATE_ANGLE rotates direction cosines through a polar angle whose cosine is
! mu and through an azimuthal angle sampled uniformly. Note that this is done
! with direct sampling rather than rejection as is done in MCNP and SERPENT.
!===============================================================================

  function rotate_angle(uvw0, mu, phi) result(uvw)
    real(8), intent(in) :: uvw0(3) ! directional cosine
    real(8), intent(in) :: mu      ! cosine of angle in lab or CM
    real(8), optional   :: phi     ! azimuthal angle
    real(8)             :: uvw(3)  ! rotated directional cosine

    real(8) :: phi_   ! azimuthal angle
    real(8) :: sinphi ! sine of azimuthal angle
    real(8) :: cosphi ! cosine of azimuthal angle
    real(8) :: a      ! sqrt(1 - mu^2)
    real(8) :: b      ! sqrt(1 - w^2)
    real(8) :: u0     ! original cosine in x direction
    real(8) :: v0     ! original cosine in y direction
    real(8) :: w0     ! original cosine in z direction

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

  end function rotate_angle

!===============================================================================
! MAXWELL_SPECTRUM samples an energy from the Maxwell fission distribution based
! on a direct sampling scheme. The probability distribution function for a
! Maxwellian is given as p(x) = 2/(T*sqrt(pi))*sqrt(x/T)*exp(-x/T). This PDF can
! be sampled using rule C64 in the Monte Carlo Sampler LA-9721-MS.
!===============================================================================

  function maxwell_spectrum(T) result(E_out)

    real(8), intent(in)  :: T     ! tabulated function of incoming E
    real(8)              :: E_out ! sampled energy

    real(8) :: r1, r2, r3  ! random numbers
    real(8) :: c           ! cosine of pi/2*r3

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

  function watt_spectrum(a, b) result(E_out)

    real(8), intent(in) :: a     ! Watt parameter a
    real(8), intent(in) :: b     ! Watt parameter b
    real(8)             :: E_out ! energy of emitted neutron

    real(8) :: w ! sampled from Maxwellian

    w     = maxwell_spectrum(a)
    E_out = w + a*a*b/4. + (TWO*prn() - ONE)*sqrt(a*a*b*w)

  end function watt_spectrum

!===============================================================================
! FADDEEVA the Faddeeva function, using Stephen Johnson's implementation
!===============================================================================

  function faddeeva(z) result(wv)
    complex(C_DOUBLE_COMPLEX), intent(in) :: z ! The point to evaluate Z at
    complex(8)     :: wv     ! The resulting w(z) value
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

  recursive function w_derivative(z, order) result(wv)
    complex(C_DOUBLE_COMPLEX), intent(in) :: z ! The point to evaluate Z at
    integer,                   intent(in) :: order
    complex(8)     :: wv     ! The resulting w(z) value

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

  subroutine broaden_wmp_polynomials(E, dopp, n, factors)
    real(8), intent(in) :: E          ! Energy to evaluate at
    real(8), intent(in) :: dopp       ! sqrt(atomic weight ratio / kT),
                                      !  kT given in eV.
    integer, intent(in) :: n          ! number of components to polynomial
    real(8), intent(out):: factors(n) ! output leading coefficient

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
