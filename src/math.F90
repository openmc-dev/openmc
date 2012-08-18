module math

  use constants, only: PI, ONE, TWO

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

end module math
