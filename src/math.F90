module math

  use constants
  use random_lcg, only: prn

  implicit none

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

end module math
