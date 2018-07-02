module endf_header

  use algorithm, only: binary_search
  use constants, only: ZERO, HISTOGRAM, LINEAR_LINEAR, LINEAR_LOG, &
       LOG_LINEAR, LOG_LOG
  use hdf5_interface

  implicit none

  type, abstract :: Function1D
  contains
    procedure(function1d_evaluate_), deferred :: evaluate
    procedure(function1d_from_hdf5_), deferred :: from_hdf5
  end type Function1D

  abstract interface
    pure function function1d_evaluate_(this, x) result(y)
      import Function1D
      class(Function1D), intent(in) :: this
      real(8),           intent(in) :: x
      real(8)                       :: y
    end function function1d_evaluate_

    subroutine function1d_from_hdf5_(this, dset_id)
      import Function1D, HID_T
      class(Function1D), intent(inout) :: this
      integer(HID_T),    intent(in)    :: dset_id
    end subroutine function1d_from_hdf5_
  end interface

!===============================================================================
! POLYNOMIAL represents a one-dimensional function expressed as a polynomial
!===============================================================================

  type, extends(Function1D) :: Polynomial
    real(8), allocatable :: coef(:)  ! coefficients
  contains
    procedure :: from_hdf5 => polynomial_from_hdf5
    procedure :: evaluate => polynomial_evaluate
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
    procedure :: from_hdf5 => tabulated1d_from_hdf5
    procedure :: evaluate => tabulated1d_evaluate
  end type Tabulated1D

contains

!===============================================================================
! Polynomial implementation
!===============================================================================

  subroutine polynomial_from_hdf5(this, dset_id)
    class(Polynomial), intent(inout) :: this
    integer(HID_T),    intent(in)    :: dset_id

    integer(HSIZE_T) :: dims(1)

    call get_shape(dset_id, dims)
    allocate(this % coef(dims(1)))
    call read_dataset(this % coef, dset_id)
  end subroutine polynomial_from_hdf5

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

  subroutine tabulated1d_from_hdf5(this, dset_id)
    class(Tabulated1D), intent(inout) :: this
    integer(HID_T), intent(in) :: dset_id

    real(8), allocatable :: xy(:,:)
    integer(HSIZE_T) :: dims(2)

    call read_attribute(this % nbt, dset_id, 'breakpoints')
    call read_attribute(this % int, dset_id, 'interpolation')
    this % n_regions = size(this % nbt)

    call get_shape(dset_id, dims)
    this % n_pairs = int(dims(1), 4)
    allocate(this % x(this % n_pairs))
    allocate(this % y(this % n_pairs))

    allocate(xy(dims(1), dims(2)))
    call read_dataset(xy, dset_id)
    this % x(:) = xy(:,1)
    this % y(:) = xy(:,2)
  end subroutine tabulated1d_from_hdf5

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
