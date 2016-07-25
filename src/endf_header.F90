module endf_header

  use constants, only: ZERO, HISTOGRAM, LINEAR_LINEAR, LINEAR_LOG, &
       LOG_LINEAR, LOG_LOG
  use search, only: binary_search

implicit none

  type, abstract :: Function1D
  contains
    procedure(function1d_evaluate_), deferred :: evaluate
  end type Function1D

  abstract interface
    pure function function1d_evaluate_(this, x) result(y)
      import Function1D
      class(Function1D), intent(in) :: this
      real(8),           intent(in) :: x
      real(8)                       :: y
    end function function1d_evaluate_
  end interface

!===============================================================================
! CONSTANT1D represents a constant one-dimensional function
!===============================================================================

  type, extends(Function1D) :: Constant1D
    real(8) :: y
  contains
    procedure :: evaluate => constant1d_evaluate
  end type Constant1D

!===============================================================================
! POLYNOMIAL represents a one-dimensional function expressed as a polynomial
!===============================================================================

  type, extends(Function1D) :: Polynomial
    real(8), allocatable :: coef(:)  ! coefficients
  contains
    procedure :: evaluate => polynomial_evaluate
    procedure :: from_ace => polynomial_from_ace
  end type Polynomial

!===============================================================================
! TABULATED1D represents a one-dimensional interpolable function
!===============================================================================

  type, extends(Function1D) :: Tabulated1D
    integer :: n_regions = 0       ! # of interpolation regions
    integer, allocatable :: nbt(:) ! values separating interpolation regions
    integer, allocatable :: int(:) ! interpolation scheme
    integer :: n_pairs             ! # of pairs of (x,y) values
    real(8), allocatable :: x(:)   ! values of abscissa
    real(8), allocatable :: y(:)   ! values of ordinate
  contains
    procedure :: from_ace => tabulated1d_from_ace
    procedure :: evaluate => tabulated1d_evaluate
  end type Tabulated1D

contains

!===============================================================================
! Constant1D implementation
!===============================================================================

  pure function constant1d_evaluate(this, x) result(y)
    class(Constant1D), intent(in) :: this
    real(8),           intent(in) :: x
    real(8)                       :: y

    y = this % y
  end function constant1d_evaluate

!===============================================================================
! Polynomial implementation
!===============================================================================

  subroutine polynomial_from_ace(this, xss, idx)
    class(Polynomial), intent(inout) :: this
    real(8),           intent(in)    :: xss(:)
    integer,           intent(in)    :: idx

    integer :: nc  ! number of coefficients (order - 1)

    ! Clear space
    if (allocated(this % coef)) deallocate(this % coef)

    ! Determine number of coefficients
    nc = nint(xss(idx))

    ! Allocate space for and read coefficients
    allocate(this % coef(nc))
    this % coef(:) = xss(idx + 1 : idx + nc)
  end subroutine polynomial_from_ace

  pure function polynomial_evaluate(this, x) result(y)
    class(Polynomial), intent(in) :: this
    real(8),           intent(in) :: x
    real(8)                       :: y

    integer :: i

    ! Use Horner's rule to evaluate polynomial. Note that coefficients are
    ! ordered in increasing powers of x.
    y = ZERO
    do i = size(this % coef), 1, -1
      y = y*x + this % coef(i)
    end do
  end function polynomial_evaluate

!===============================================================================
! Tabulated1D implementation
!===============================================================================

  subroutine tabulated1d_from_ace(this, xss, idx)
    class(Tabulated1D), intent(inout) :: this
    real(8), intent(in) :: xss(:)
    integer, intent(in) :: idx

    integer :: nr, ne

    ! Clear space
    if (allocated(this % nbt)) deallocate(this % nbt)
    if (allocated(this % int)) deallocate(this % int)
    if (allocated(this % x)) deallocate(this % x)
    if (allocated(this % y)) deallocate(this % y)

    ! Determine number of regions
    nr = nint(xss(idx))
    this%n_regions = nr

    ! Read interpolation region data
    if (nr > 0) then
      allocate(this%nbt(nr))
      allocate(this%int(nr))
      this%nbt(:) = nint(xss(idx + 1 : idx + nr))
      this%int(:) = nint(xss(idx + nr + 1 : idx + 2*nr))
    end if

    ! Determine number of pairs
    ne = int(XSS(idx + 2*nr + 1))
    this%n_pairs = ne

    ! Read (x,y) pairs
    allocate(this%x(ne))
    allocate(this%y(ne))
    this%x(:) = xss(idx + 2*nr + 2 : idx + 2*nr + 1 + ne)
    this%y(:) = xss(idx + 2*nr + 2 + ne : idx + 2*nr + 1 + 2*ne)
  end subroutine tabulated1d_from_ace

  pure function tabulated1d_evaluate(this, x) result(y)
    class(Tabulated1D), intent(in) :: this
    real(8),     intent(in) :: x    ! x value to find y at
    real(8)                 :: y    ! y(x)

    integer :: i               ! bin in which to interpolate
    integer :: j               ! index for interpolation region
    integer :: n_regions       ! number of interpolation regions
    integer :: n_pairs         ! number of tabulated values
    integer :: interp          ! ENDF interpolation scheme
    real(8) :: r               ! interpolation factor
    real(8) :: x0, x1          ! bounding x values
    real(8) :: y0, y1          ! bounding y values

    ! determine number of interpolation regions and pairs
    n_regions = this % n_regions
    n_pairs   = this % n_pairs

    ! find which bin the abscissa is in -- if the abscissa is outside the
    ! tabulated range, the first or last point is chosen, i.e. no interpolation
    ! is done outside the energy range
    if (x < this % x(1)) then
      y = this % y(1)
      return
    elseif (x > this % x(n_pairs)) then
      y = this % y(n_pairs)
      return
    else
      i = binary_search(this % x, n_pairs, x)
    end if

    ! determine interpolation scheme
    if (n_regions == 0) then
      interp = LINEAR_LINEAR
    elseif (n_regions == 1) then
      interp = this % int(1)
    elseif (n_regions > 1) then
      do j = 1, n_regions
        if (i < this % nbt(j)) then
          interp = this % int(j)
          exit
        end if
      end do
    end if

    ! handle special case of histogram interpolation
    if (interp == HISTOGRAM) then
      y = this % y(i)
      return
    end if

    ! determine bounding values
    x0 = this % x(i)
    x1 = this % x(i + 1)
    y0 = this % y(i)
    y1 = this % y(i + 1)

    ! determine interpolation factor and interpolated value
    select case (interp)
    case (LINEAR_LINEAR)
      r = (x - x0)/(x1 - x0)
      y = y0 + r*(y1 - y0)
    case (LINEAR_LOG)
      r = log(x/x0)/log(x1/x0)
      y = y0 + r*(y1 - y0)
    case (LOG_LINEAR)
      r = (x - x0)/(x1 - x0)
      y = y0*exp(r*log(y1/y0))
    case (LOG_LOG)
      r = log(x/x0)/log(x1/x0)
      y = y0*exp(r*log(y1/y0))
    end select

  end function tabulated1d_evaluate

end module endf_header
