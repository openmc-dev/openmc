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
  end type Tab1

end module endf_header
