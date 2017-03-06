!>@brief Interfaces with C functions to evaluate the Faddeeva function

!> This is a slightly modified and trimmed down version of the
!! math.F90 module taken from ClosedMC, July 2014
!! (https://github.com/cjosey/closedmc/tree/multipole).  Elements of
!! that module are taken from codes written at, and obtained from, ANL.
module URR_faddeeva

  use URR_constants, only: ZERO, HALF, ONE, ONEI
  use ISO_C_BINDING
 
  implicit none
  private
  public :: quickw,&
       faddeeva_w,&
       tabulate_w

  ! Tabulated Faddeeva evaluations for use by QUICKW
  complex(8) :: w_tabulated(-1:60,-1:60)
  

!> Evaluates the scaled complementary error function.  This
!! interfaces with the MIT C library.
  interface
    COMPLEX (C_DOUBLE_COMPLEX) FUNCTION faddeeva_w &
            (Z, RELERR) BIND(C, NAME='Faddeeva_w')
            use ISO_C_BINDING
            implicit none
            COMPLEX(C_DOUBLE_COMPLEX), value :: Z
            REAL(C_DOUBLE)           , value :: RELERR
    end function faddeeva_w
  end interface

contains


!> Calculates the Faddeeva function, also known as the complex probability
!! integral, for complex arguments. For |z| < 6, it uses a six-point
!! interpolation scheme based on pre-tabulated data that is accurate to
!! O(10^-3). For |z| > 6, it uses a three-term asymptotic approximation that is
!! accurate to O(10^-6).
  function quickw(z) result(w)

    complex(8), intent(in) :: z
    complex(8)             :: w

    real(8) :: p   ! interpolation factor on real axis
    real(8) :: q   ! interpolation factor on imaginary axis
    real(8) :: pp  ! p*p
    real(8) :: qq  ! q*q
    real(8) :: pq  ! p*q
    real(8) :: a_l ! coefficient for left point
    real(8) :: a_c ! coefficient for center point
    real(8) :: a_b ! coefficient for bottom point
    real(8) :: a_r ! coefficient for right point
    real(8) :: a_t ! coefficient for top point

    ! Asymptotic expansion parameters
    ! Expanded ahead of time to save in performance.
    real(8) :: a = 0.512424224754768462984202823134979415014943561548661637413182_8
    real(8) :: b = 0.275255128608410950901357962647054304017026259671664935783653_8
    real(8) :: c = 0.051765358792987823963876628425793170829107067780337219430904_8
    real(8) :: d = 2.724744871391589049098642037352945695982973740328335064216346_8

    integer :: l           ! interpolation index for real axis
    integer :: m           ! interpolation index for imaginary axis

    if (abs(z) < 6.) then
      ! Use interpolation for |z| < 6. The interpolation scheme uses a bivariate
      ! six-point quadrature described in Abramowitz and Stegun 25.2.67. This
      ! interpolation is accurate to O(h^3) = O(10^-3).
      !      
      !     l-1  l  l+1
      ! m+1      +   +
      !          |
      ! m    +---+---+
      !          |
      ! m-1      +

      ! Determine indices on grid for interpolation and interpolation factors --
      ! note that in previous implementations it was necessary to add/subtract
      ! two in places because of the indexing on the tabulated function. Because
      ! w_tabulated is indexed from -1 to 60, we don't need to do that here
      p = 10. * abs(real(z))
      q = 10. * aimag(z)
      l = int(p)
      m = int(q)
      p = p - l
      q = q - m
      pp = p * p
      qq = q * q
      pq = p * q

      ! Coefficients for interpolation
      a_b = HALF*(qq - q)       ! bottom
      a_l = HALF*(pp - p)       ! left
      a_c = ONE + pq - pp - qq  ! center
      a_r = HALF*(pp + p) - pq  ! right
      a_t = HALF*(qq + q) - pq  ! top

      ! Use six-point interpolation to calculate real and imaginary parts
      w =  a_b * w_tabulated(l,m-1) + a_l * w_tabulated(l-1,m) + &
           a_c * w_tabulated(l,m)   + a_r * w_tabulated(l+1,m) + &
           a_t * w_tabulated(l,m+1) + pq  * w_tabulated(l+1,m+1)

      ! Handle arguments with negative real parts
      if (real(z) < 0) w = conjg(w)

    else
      ! Use three-term asymptotic expansion for |z| > 6
      w = ONEI * z * (a / (z * z - b) + c / (z * z - d))

    end if

  end function quickw


!> Calculates the Faddeeva (W) function on a 62 x 62 grid
  subroutine tabulate_w()

    integer :: i ! real loop index
    integer :: j ! imaginary loop index
    real(C_DOUBLE), parameter :: w = 0.1_8 ! width on grid
    real(C_DOUBLE)            :: x         ! real part
    real(C_DOUBLE)            :: y         ! imaginary part
    complex(C_DOUBLE_COMPLEX) :: z         ! argument to Faddeeva function

    do j = -1, 60
      y = w * j

      do i = -1, 60
        x = w * i
        z = cmplx(x, y, 8)
        w_tabulated(i,j) = faddeeva_w(z, ZERO)

      end do
    end do

  end subroutine tabulate_w


end module URR_faddeeva
