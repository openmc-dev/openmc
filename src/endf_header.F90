module endf_header

  implicit none

!===============================================================================
! TAB1 represents a one-dimensional interpolable function
!===============================================================================

  type Tab1
    integer :: n_regions = 0       ! # of interpolation regions
    integer, allocatable :: nbt(:) ! values separating interpolation regions
    integer, allocatable :: int(:) ! interpolation scheme
    integer :: n_pairs             ! # of pairs of (x,y) values
    real(8), allocatable :: x(:)   ! values of abscissa
    real(8), allocatable :: y(:)   ! values of ordinate
  contains
    procedure :: from_ace
  end type Tab1

contains

  subroutine from_ace(this, xss, idx)
    class(Tab1), intent(inout) :: this
    real(8), intent(in) :: xss(:)
    integer, intent(in) :: idx

    integer :: nr, ne

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
  end subroutine from_ace

end module endf_header
