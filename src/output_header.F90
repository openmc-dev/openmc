module output_header 

  use, intrinsic :: ISO_FORTRAN_ENV

  implicit none
  
  public

  ! Short names for output and error units
  integer :: ou = OUTPUT_UNIT
  integer :: eu = ERROR_UNIT

  contains
!===============================================================================
! OUTPUT_MESSAGE displays an informational message to the log file and the
! standard output stream.
!===============================================================================

    subroutine output_message(msg)

      character(len=*), intent(in) :: msg

      integer :: i_start    ! starting position
      integer :: i_end      ! ending position
      integer :: line_wrap  ! length of line
      integer :: length     ! length of message
      integer :: last_space ! index of last space (relative to start)

      ! Set length of line
      line_wrap = 80

      ! Determine length of message
      length = len_trim(msg)

      i_start = 0
      do
        if (length - i_start < line_wrap - 1) then
          ! Remainder of message will fit on line
          write(ou, fmt='(1X,A)') msg(i_start+1:length)
          exit

        else
          ! Determine last space in current line
          last_space = index(msg(i_start+1:i_start+line_wrap), &
               ' ', BACK=.true.)
          if (last_space == 0) then
            i_end = min(length + 1, i_start+line_wrap) - 1
            write(ou, fmt='(1X,A)') msg(i_start+1:i_end)
          else
            i_end = i_start + last_space
            write(ou, fmt='(1X,A)') msg(i_start+1:i_end-1)
          end if

          ! Write up to last space

          ! Advance starting position
          i_start = i_end
          if (i_start > length) exit
        end if
      end do

    end subroutine output_message
    
end module output_header 
