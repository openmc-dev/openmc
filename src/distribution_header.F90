module distribution_header

  use constants, only: ZERO, HALF, ONE, HISTOGRAM, LINEAR_LINEAR
  use error, only: fatal_error
  use math, only: rotate_angle, maxwell_spectrum, watt_spectrum
  use random_lcg, only: prn

!===============================================================================
! DISTRIBUTION type defines a probability density function
!===============================================================================

  type, abstract :: Distribution
  contains
    procedure(iSample), deferred :: sample
  end type Distribution

  abstract interface
    function iSample(this) result(x)
      import Distribution
      class(Distribution), intent(in) :: this
      real(8) :: x
    end function iSample
  end interface

!===============================================================================
! Derived classes of Distribution
!===============================================================================

  ! delta function at a single point
  type, extends(Distribution) :: Delta
    real(8) :: x0
  contains
    procedure :: sample => delta_sample
  end type Delta

  ! Uniform distribution over the interval [a,b]
  type, extends(Distribution) :: Uniform
    real(8) :: a
    real(8) :: b
  contains
    procedure :: sample => uniform_sample
  end type Uniform

  ! Maxwellian distribution of form c*E*exp(-E/a)
  type, extends(Distribution) :: Maxwell
    real(8) :: theta
  contains
    procedure :: sample => maxwell_sample
  end type Maxwell

  ! Watt fission spectrum with form c*exp(-E/a)*sinh(sqrt(b*E))
  type, extends(Distribution) :: Watt
    real(8) :: a
    real(8) :: b
  contains
    procedure :: sample => watt_sample
  end type Watt

  ! Histogram or linear-linear interpolated tabular distribution
  type, extends(Distribution) :: Tabular
    integer :: interpolation
    real(8), allocatable :: x(:) ! tabulated independent variable
    real(8), allocatable :: p(:) ! tabulated probability density
    real(8), allocatable, private :: c(:) ! cumulative distribution at tabulated values
  contains
    procedure :: sample => tabular_sample
    procedure :: initialize => tabular_initialize
  end type Tabular

!===============================================================================
! AngleDistribution
!===============================================================================

  type :: AngleDistribution
    real(8) :: reference_uvw(3)
    class(Distribution), allocatable :: mu
  contains
    procedure :: sample => angle_sample
  end type AngleDistribution

contains

  function delta_sample(this) result(x)
    class(Delta), intent(in) :: this
    real(8) :: x

    x = this%x0
  end function delta_sample

  function uniform_sample(this) result(x)
    class(Uniform), intent(in) :: this
    real(8) :: x

    x = this%a + prn()*(this%b - this%a)
  end function uniform_sample

  function maxwell_sample(this) result(x)
    class(Maxwell), intent(in) :: this
    real(8) :: x

    x = maxwell_spectrum(this%theta)
  end function maxwell_sample

  function watt_sample(this) result(x)
    class(Watt), intent(in) :: this
    real(8) :: x

    x = watt_spectrum(this%a, this%b)
  end function watt_sample

  function tabular_sample(this) result(x)
    class(Tabular), intent(in) :: this
    real(8) :: x

    integer :: i
    real(8) :: c         ! sampled cumulative frequency
    real(8) :: m         ! slope of PDF
    real(8) :: x_i, x_i1 ! i-th and (i+1)th x values
    real(8) :: c_i, c_i1 ! i-th and (i+1)th cumulative distribution values
    real(8) :: p_i, p_i1 ! i-th and (i+1)th probability density values

    ! Sample value of CDF
    c = prn()

    ! Find first CDF bin which is above the sampled value
    c_i = this%c(1)
    do i = 1, size(this%c) - 1
      c_i1 = this%c(i + 1)
      if (c <= c_i1) exit
      c_i = c_i1
    end do

    ! Determine bounding PDF values
    x_i = this%x(i)
    p_i = this%p(i)

    if (this%interpolation == HISTOGRAM) then
      ! Histogram interpolation
      if (p_i > ZERO) then
        x = x_i + (c - c_i)/p_i
      else
        x = x_i
      end if
    else
      ! Linear-linear interpolation
      x_i1 = this%x(i + 1)
      p_i1 = this%p(i + 1)

      m = (p_i1 - p_i)/(x_i1 - x_i)
      if (m == ZERO) then
        x = x_i + (c - c_i)/p_i
      else
        x = x_i + (sqrt(max(ZERO, p_i*p_i + 2*m*(c - c_i))) - p_i)/m
      end if
    end if
  end function tabular_sample

  subroutine tabular_initialize(this, x, p, interp)
    class(Tabular), intent(inout) :: this
    real(8), intent(in) :: x(:)
    real(8), intent(in) :: p(:)
    integer, intent(in) :: interp

    integer :: n

    ! Check interpolation parameter

    if (interp /= HISTOGRAM .and. interp /= LINEAR_LINEAR) then
      call fatal_error('Only histogram and linear-linear interpolation for tabular &
           &distribution is supported.')
    end if

    ! Check length of x, p arrays
    if (size(x) /= size(p)) then
      call fatal_error('Tabulated probabilities not of same length as &
           &independent variable.')
    end if

    ! Copy probability density function and interpolation parameter
    n = size(x)
    allocate(this%x(n), this%p(n), this%c(n))
    this%interpolation = interp
    this%x(:) = x(:)
    this%p(:) = p(:)

    ! Calculate cumulative distribution function
    this%c(1) = ZERO
    do i = 2, n
      if (this%interpolation == HISTOGRAM) then
        this%c(i) = this%c(i-1) + this%p(i-1)*(this%x(i) - this%x(i-1))
      elseif (this%interpolation == LINEAR_LINEAR) then
        this%c(i) = this%c(i-1) + HALF*(this%p(i-1) + this%p(i)) * &
             (this%x(i) - this%x(i-1))
      end if
    end do

    ! Normalize density and distribution functions
    this%p(:) = this%p(:)/this%c(n)
    this%c(:) = this%c(:)/this%c(n)
  end subroutine tabular_initialize

  function angle_sample(this) result(uvw)
    class(AngleDistribution), intent(in) :: this
    real(8) :: uvw(3)

    real(8) :: mu

    mu = this%mu%sample()
    if (mu == ONE) then
      uvw(:) = this%reference_uvw
    else
      uvw(:) = rotate_angle(this%reference_uvw, mu)
    end if
  end function angle_sample

end module distribution_header
