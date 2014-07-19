module doppler

  use constants, only: ZERO, ONE, PI, K_BOLTZMANN

  implicit none

  real(8), parameter :: sqrt_pi_inv = ONE / sqrt(PI)

contains

!===============================================================================
! BROADEN takes a microscopic cross section at a temperature T_1 and Doppler
! broadens it to a higher temperature T_2 based on a method originally developed
! by Cullen and Weisbin (see "Exact Doppler Broadening of Tabulated Cross
! Sections," Nucl. Sci. Eng. 60, 199-229 (1976)). The only difference here is
! the F functions are evaluated based on complementary error functions rather
! than error functions as is done in the BROADR module of NJOY.
!===============================================================================

  subroutine broaden(energy, xs, A_target, T, sigmaNew)

    real(8), intent(in)  :: energy(:)   ! energy grid
    real(8), intent(in)  :: xs(:)       ! unbroadened cross section
    integer, intent(in)  :: A_target    ! mass number of target
    real(8), intent(in)  :: T           ! temperature (difference)
    real(8), intent(out) :: sigmaNew(:) ! broadened cross section

    integer              :: i, k     ! loop indices
    integer              :: n        ! number of energy points
    real(8)              :: F_a(0:4) ! F(a) functions as per C&W
    real(8)              :: F_b(0:4) ! F(b) functions as per C&W
    real(8)              :: H(0:4)   ! H functions as per C&W
    real(8), allocatable :: x(:)     ! proportional to relative velocity
    real(8)              :: y        ! proportional to neutron velocity
    real(8)              :: y_sq     ! y**2
    real(8)              :: y_inv    ! 1/y
    real(8)              :: y_inv_sq ! 1/y**2
    real(8)              :: alpha    ! constant equal to A/kT
    real(8)              :: slope    ! slope of xs between adjacent points
    real(8)              :: Ak, Bk   ! coefficients at each point
    real(8)              :: a, b     ! values of x(k)-y and x(k+1)-y
    real(8)              :: sigma    ! broadened cross section at one point

    ! Determine alpha parameter -- have to convert k to MeV/K
    alpha = A_target/(K_BOLTZMANN * T)

    ! Allocate memory for x and assign values
    n = size(energy)
    allocate(x(n))
    x = sqrt(alpha * energy)

    ! Loop over incoming neutron energies
    ENERGY_NEUTRON: do i = 1, n

       sigma    = ZERO
       y        = x(i)
       y_sq     = y*y
       y_inv    = ONE / y
       y_inv_sq = y_inv / y

       ! =======================================================================
       ! EVALUATE FIRST TERM FROM x(k) - y = 0 to -4

       k = i
       a = ZERO
       call calculate_F(F_a, a)

       do while (a >= -4.0 .and. k > 1)
          ! Move to next point
          F_b = F_a
          k = k - 1
          a = x(k) - y

          ! Calculate F and H functions
          call calculate_F(F_a, a)
          H = F_a - F_b

          ! Calculate A(k), B(k), and slope terms
          Ak = y_inv_sq*H(2) + 2.0*y_inv*H(1) + H(0)
          Bk = y_inv_sq*H(4) + 4.0*y_inv*H(3) + 6.0*H(2) + 4.0*y*H(1) + y_sq*H(0)
          slope = (xs(k+1) - xs(k)) / (x(k+1)**2 - x(k)**2)

          ! Add contribution to broadened cross section
          sigma = sigma + Ak*(xs(k) - slope*x(k)**2) + slope*Bk
       end do

       ! =======================================================================
       ! EXTEND CROSS SECTION TO 0 ASSUMING 1/V SHAPE

       if (k == 1 .and. a >= -4.0) then
          ! Since x = 0, this implies that a = -y
          F_b = F_a
          a = -y

          ! Calculate F and H functions
          call calculate_F(F_a, a)
          H = F_a - F_b

          ! Add contribution to broadened cross section
          sigma = sigma + xs(k)*x(k)*(y_inv_sq*H(1) + y_inv*H(0))
       end if

       ! =======================================================================
       ! EVALUATE FIRST TERM FROM x(k) - y = 0 to 4

       k = i
       b = ZERO
       call calculate_F(F_b, b)

       do while (b <= 4.0 .and. k < n)
          ! Move to next point
          F_a = F_b
          k = k + 1
          b = x(k) - y

          ! Calculate F and H functions
          call calculate_F(F_b, b)
          H = F_a - F_b

          ! Calculate A(k), B(k), and slope terms
          Ak = y_inv_sq*H(2) + 2.0*y_inv*H(1) + H(0)
          Bk = y_inv_sq*H(4) + 4.0*y_inv*H(3) + 6.0*H(2) + 4.0*y*H(1) + y_sq*H(0)
          slope = (xs(k) - xs(k-1)) / (x(k)**2 - x(k-1)**2)

          ! Add contribution to broadened cross section
          sigma = sigma + Ak*(xs(k) - slope*x(k)**2) + slope*Bk
       end do

       ! =======================================================================
       ! EXTEND CROSS SECTION TO INFINITY ASSUMING CONSTANT SHAPE

       if (k == n .and. b <= 4.0) then
          ! Calculate F function at last energy point
          a = x(k) - y
          call calculate_F(F_a, a)

          ! Add contribution to broadened cross section
          sigma = sigma + xs(k) * (y_inv_sq*F_a(2) + 2.0*y_inv*F_a(1) + F_a(0))
       end if

       ! =======================================================================
       ! EVALUATE SECOND TERM FROM x(k) + y = 0 to +4

       if (y <= 4.0) then
          ! Swap signs on y
          y = -y
          y_inv = -y_inv
          k = 1

          ! Calculate a and b based on 0 and x(1)
          a = -y
          b = x(k) - y

          ! Calculate F and H functions
          call calculate_F(F_a, a)
          call calculate_F(F_b, b)
          H = F_a - F_b

          ! Add contribution to broadened cross section
          sigma = sigma - xs(k) * x(k) * (y_inv_sq*H(1) + y_inv*H(0))

          ! Now progress forward doing the remainder of the second term
          do while (b <= 4.0)
             ! Move to next point
             F_a = F_b
             k = k + 1
             b = x(k) - y

             ! Calculate F and H functions
             call calculate_F(F_b, b)
             H = F_a - F_b

             ! Calculate A(k), B(k), and slope terms
             Ak = y_inv_sq*H(2) + 2.0*y_inv*H(1) + H(0)
             Bk = y_inv_sq*H(4) + 4.0*y_inv*H(3) + 6.0*H(2) + 4.0*y*H(1) + y_sq*H(0)
             slope = (xs(k) - xs(k-1)) / (x(k)**2 - x(k-1)**2)

             ! Add contribution to broadened cross section
             sigma = sigma - Ak*(xs(k) - slope*x(k)**2) - slope*Bk
          end do
       end if

       ! Set broadened cross section
       sigmaNew(i) = sigma

    end do ENERGY_NEUTRON

  end subroutine broaden

!===============================================================================
! CALCULATE_F evaluates the function:
!
!    F(n,a) = 1/sqrt(pi)*int(z^n*exp(-z^2), z = a to infinity)
!
! The five values returned in a vector correspond to the integral for n = 0
! through 4. These functions are called over and over during the Doppler
! broadening routine.
!===============================================================================

  subroutine calculate_F(F, a)

    real(8), intent(inout) :: F(0:4)
    real(8), intent(in)    :: a

#ifndef NO_F2008
    F(0) = 0.5*erfc(a)
#endif
    F(1) = 0.5*sqrt_pi_inv*exp(-a*a)
    F(2) = 0.5*F(0) + a*F(1)
    F(3) = F(1)*(1.0 + a*a)
    F(4) = 0.75*F(0) + F(1)*a*(1.5 + a*a)

  end subroutine calculate_F

end module doppler
