module progress_header

  use, intrinsic :: ISO_FORTRAN_ENV,      only: OUTPUT_UNIT

  implicit none

  interface
    function check_isatty(fd) bind(C, name = 'isatty')
      use, intrinsic :: ISO_C_BINDING, only: c_int
      integer(c_int)        :: check_isatty
      integer(c_int), value :: fd
    end function
  end interface

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

    if (check_isatty(1) == 0) return

    ! set the percentage
    if (val >= 100.) then
      write(self % bar(1:3), "(I3)") 100
    else if (val <= 0.) then
      write(self % bar(1:3), "(I3)") 0
    else
      write(self % bar(1:3), "(I3)") int(val)
    end if

    ! set the bar width
    if (val >= 100.) then
      do i=1,65
        write(self % bar(i+6:i+6), '(A)') '='
      end do
    else
      do i=1,int(dble(65)*val/100.)
        write(self % bar(i+6:i+6), '(A)') '='
      end do
    end if

    write(OUTPUT_UNIT, '(A1,A1,A72)', ADVANCE='no') '+', char(13), self % bar
    flush(OUTPUT_UNIT)
    
    if (val >= 100.) then
    
      ! make new line
      write(OUTPUT_UNIT, "(A)") ""
      flush(OUTPUT_UNIT)
      
      ! reset the bar in case we want to use this instance again
      self % bar = "???% |                                      " // &
                   "                           |"
      
    end if
    
  end subroutine bar_set_value

end module progress_header
