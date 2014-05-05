module math

  use constants,  only: PI, ONE, TWO
  use random_lcg, only: prn

  implicit none

contains

!===============================================================================
! NORMAL_PERCENTILE calculates the percentile of the standard normal
! distribution with a specified probability level
!===============================================================================

  function normal_percentile(p) result(z)

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
           ((((d(1)*q + d(2))*q + d(3))*q + d(4))*q + 1.)

    elseif (p <= 1. - p_low) then
      ! Rational approximation for central region

      q = p - 0.5
      r = q*q
      z = (((((a(1)*r + a(2))*r + a(3))*r + a(4))*r + a(5))*r + a(6))*q / &
           (((((b(1)*r + b(2))*r + b(3))*r + b(4))*r + b(5))*r + 1.)

    else
      ! Rational approximation for upper region

      q = sqrt(-2*log(1. - p))
      z = -(((((c(1)*q + c(2))*q + c(3))*q + c(4))*q + c(5))*q + c(6)) / &
           ((((d(1)*q + d(2))*q + d(3))*q + d(4))*q + 1.)
    endif

    ! Refinement based on Newton's method
#ifndef NO_F2008
    z = z - (0.5 * erfc(-z/sqrt(TWO)) - p) * sqrt(TWO*PI) * exp(0.5*z*z)
#endif

  end function normal_percentile

!===============================================================================
! T_PERCENTILE calculates the percentile of the Student's t distribution with a
! specified probability level and number of degrees of freedom
!===============================================================================

  function t_percentile(p, df) result(t)

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

      t = tan(PI*(p - 0.5))

    elseif (df == 2) then
      ! For two degrees of freedom, the cdf is given by 1/2 + x/(2*sqrt(x^2 +
      ! 2)). This can be directly inverted to yield the solution below

      t = TWO*sqrt(TWO)*(p - 0.5)/sqrt(ONE - 4.*(p - 0.5)**2)

    else

      ! This approximation is from E. Olusegun George and Meenakshi Sivaram, "A
      ! modification of the Fisher-Cornish approximation for the student t
      ! percentiles," Communication in Statistics - Simulation and Computation,
      ! 16 (4), pp. 1123-1132 (1987).

      n = real(df,8)
      k = 1./(n - 2.)
      z = normal_percentile(p)
      z2 = z * z
      t = sqrt(n*k) * (z + (z2 - 3.)*z*k/4. + ((5.*z2 - 56.)*z2 + &
           75.)*z*k*k/96. + (((z2 - 27.)*3.*z2 + 417.)*z2 - 315.) &
           *z*k*k*k/384.)

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
  
  pure function calc_pn(n,x) result(pnx)

    integer, intent(in) :: n   ! Legendre order requested
    real(8), intent(in) :: x   ! Independent variable the Legendre is to be 
                               ! evaluated at; x must be in the domain [-1,1]
    real(8)             :: pnx ! The Legendre poly of order n evaluated at x
    
    select case(n)
    case(1)
      pnx = x
    case(2)
      pnx = 1.5_8 * x * x - 0.5_8
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
