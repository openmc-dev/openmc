module progress_header

  use, intrinsic :: ISO_FORTRAN_ENV

  implicit none

!===============================================================================
! PROGRESSBAR
!===============================================================================

  type ProgressBar
    private
    character(len=72) :: bar="???% |                                      " // &
                             "                           |"
   contains
     procedure :: set_value => bar_set_value
  end type ProgressBar

contains
  
!===============================================================================
! BAR_SET_VALUE prints the progress bar without advancing.  The value is
! specified as percent completion, from 0 to 100.  If the value is ever set to
! 100 or above, the 
!===============================================================================

  subroutine bar_set_value(self, val)

    class(ProgressBar), intent(inout) :: self
    real(8),            intent(in)    :: val
    
    integer :: i

    ! set the percentage
    if (val >= 100.) then
      write(unit=self % bar(1:3),fmt="(i3)") 100
    else
      write(unit=self % bar(1:3),fmt="(i3)") int(val)
    end if

    ! set the bar width
    if (val >= 100.) then
      do i=1,65
        self % bar(i+6:i+6) = '='
      end do
    else
      do i=1,int(dble(65)*val/100.)
        self % bar(i+6:i+6) = '='
      end do
    end if

    open(UNIT=OUTPUT_UNIT, carriagecontrol='fortran')
    write(UNIT=OUTPUT_UNIT, FMT='(a1,a1,a72)', ADVANCE='no') '+', char(13), &
                                                             self % bar
    call flush(OUTPUT_UNIT)
    
    if (val >= 100.) then
    
      ! make new line
      write(UNIT=OUTPUT_UNIT, FMT="(A)") ""
    
      ! reset the bar in case we want to use this instance again
      self % bar = "???% |                                      " // &
                   "                           |"
      
      close(UNIT=OUTPUT_UNIT)
      
    end if
    
  end subroutine bar_set_value

end module progress_header
