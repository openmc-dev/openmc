module distribution_univariate

  use constants,  only: ZERO, ONE, HALF, HISTOGRAM, LINEAR_LINEAR, &
       MAX_LINE_LEN, MAX_WORD_LEN
  use error,      only: fatal_error
  use random_lcg, only: prn
  use math,       only: maxwell_spectrum, watt_spectrum
  use string,     only: to_lower
  use xml_interface

  implicit none

!===============================================================================
! DISTRIBUTION type defines a probability density function
!===============================================================================

  type, abstract :: Distribution
  contains
    procedure(distribution_sample_), deferred :: sample
  end type Distribution

  type DistributionContainer
    class(Distribution), allocatable :: obj
  end type DistributionContainer

  abstract interface
    function distribution_sample_(this) result(x)
      import Distribution
      class(Distribution), intent(in) :: this
      real(8) :: x
    end function distribution_sample_
  end interface

!===============================================================================
! Derived classes of Distribution
!===============================================================================

  ! Discrete distribution
  type, extends(Distribution) :: Discrete
    real(8), allocatable :: x(:)
    real(8), allocatable :: p(:)
  contains
    procedure :: sample => discrete_sample
    procedure :: initialize => discrete_initialize
  end type Discrete

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
    real(8), allocatable :: c(:) ! cumulative distribution at tabulated values
  contains
    procedure :: sample => tabular_sample
    procedure :: initialize => tabular_initialize
  end type Tabular

  type, extends(Distribution) :: Equiprobable
    real(8), allocatable :: x(:)
  contains
    procedure :: sample => equiprobable_sample
  end type Equiprobable

contains

  function discrete_sample(this) result(x)
    class(Discrete), intent(in) :: this
    real(8) :: x

    integer :: i  ! loop counter
    integer :: n  ! size of distribution
    real(8) :: c  ! cumulative frequency
    real(8) :: xi ! sampled CDF value

    n = size(this%x)
    if (n > 1) then
      xi = prn()
      c = ZERO
      do i = 1, size(this%x)
        c = c + this%p(i)
        if (xi < c) exit
      end do
      x = this%x(i)
    else
      x = this%x(1)
    end if
  end function discrete_sample

  subroutine discrete_initialize(this, x, p)
    class(Discrete), intent(inout) :: this
    real(8), intent(in) :: x(:)
    real(8), intent(in) :: p(:)

    integer :: n

    ! Check length of x, p arrays
    if (size(x) /= size(p)) then
      call fatal_error('Tabulated probabilities not of same length as &
           &independent variable.')
    end if

    ! Copy probability density function
    n = size(x)
    allocate(this%x(n), this%p(n))
    this%x(:) = x(:)
    this%p(:) = p(:)

    ! Normalize density function
    this%p(:) = this%p(:)/sum(this%p)
  end subroutine

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

    integer :: i
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

  function equiprobable_sample(this) result(x)
    class(Equiprobable), intent(in) :: this
    real(8) :: x

    integer :: i
    integer :: n
    real(8) :: r
    real(8) :: xl, xr

    n = size(this%x)

    r = prn()
    i = 1 + int((n - 1)*r)

    xl = this%x(i)
    xr = this%x(i+1)
    x = xl + ((n - 1)*r - i + ONE) * (xr - xl)
  end function equiprobable_sample

  subroutine distribution_from_xml(dist, node_dist)
    class(Distribution), allocatable, intent(inout) :: dist
    type(Node), pointer :: node_dist

    character(MAX_WORD_LEN) :: type
    character(MAX_LINE_LEN) :: temp_str
    integer :: n
    integer :: temp_int
    real(8), allocatable :: temp_real(:)

    if (check_for_node(node_dist, "type")) then
      ! Determine type of distribution
      call get_node_value(node_dist, "type", type)

      ! Determine number of parameters specified
      if (check_for_node(node_dist, "parameters")) then
        n = get_arraysize_double(node_dist, "parameters")
      else
        n = 0
      end if

      ! Allocate extension of Distribution
      select case (to_lower(type))
      case ('uniform')
        allocate(Uniform :: dist)
        if (n /= 2) then
          call fatal_error('Uniform distribution must have two &
               &parameters specified.')
        end if

      case ('maxwell')
        allocate(Maxwell :: dist)
        if (n /= 1) then
          call fatal_error('Maxwell energy distribution must have one &
               &parameter specified.')
        end if

      case ('watt')
        allocate(Watt :: dist)
        if (n /= 2) then
          call fatal_error('Watt energy distribution must have two &
               &parameters specified.')
        end if

      case ('discrete')
        allocate(Discrete :: dist)

      case ('tabular')
        allocate(Tabular :: dist)

      case default
        call fatal_error('Invalid distribution type: ' // trim(type) // '.')

      end select

      ! Read parameters and interpolation for distribution
      select type (dist)
      type is (Uniform)
        allocate(temp_real(2))
        call get_node_array(node_dist, "parameters", temp_real)
        dist%a = temp_real(1)
        dist%b = temp_real(2)
        deallocate(temp_real)

      type is (Maxwell)
        call get_node_value(node_dist, "parameters", dist%theta)

      type is (Watt)
        allocate(temp_real(2))
        call get_node_array(node_dist, "parameters", temp_real)
        dist%a = temp_real(1)
        dist%b = temp_real(2)
        deallocate(temp_real)

      type is (Discrete)
        allocate(temp_real(n))
        call get_node_array(node_dist, "parameters", temp_real)
        call dist%initialize(temp_real(1:n/2), temp_real(n/2+1:n))
        deallocate(temp_real)

      type is (Tabular)
        ! Read interpolation
        if (check_for_node(node_dist, "interpolation")) then
          call get_node_value(node_dist, "interpolation", temp_str)
          select case(to_lower(temp_str))
          case ('histogram')
            temp_int = HISTOGRAM
          case ('linear-linear')
            temp_int = LINEAR_LINEAR
          case default
            call fatal_error("Unknown interpolation type for distribution: " &
                 // trim(temp_str))
          end select
        else
          temp_int = HISTOGRAM
        end if

        ! Read and initialize tabular distribution
        allocate(temp_real(n))
        call get_node_array(node_dist, "parameters", temp_real)
        call dist%initialize(temp_real(1:n/2), temp_real(n/2+1:n), temp_int)
        deallocate(temp_real)
      end select
    end if
  end subroutine distribution_from_xml

end module distribution_univariate
