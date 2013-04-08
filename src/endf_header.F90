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
    
    ! Type-Bound procedures
    contains
      procedure :: clear => tab1_clear ! deallocates a Tab1 Object.
  end type Tab1
  
  contains
  
!===============================================================================
! TAB1_CLEAR deallocates the items in Tab1
!===============================================================================

    subroutine tab1_clear(this)
      
      class(Tab1), intent(inout) :: this ! The Tab1 to clear
      
      if (allocated(this % nbt)) &
           deallocate(this % nbt, this % int)
        
      if (allocated(this % x)) &
           deallocate(this % x, this % y)
        
    end subroutine tab1_clear

end module endf_header
