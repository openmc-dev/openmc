module error

  use ISO_FORTRAN_ENV
  use global

  implicit none

  integer :: ou = OUTPUT_UNIT
  integer :: eu = ERROR_UNIT

contains

!===============================================================================
! WARNING issues a warning to the user in the log file and the standard output
! stream.
!===============================================================================

  subroutine warning(msg)

    character(*), intent(in) :: msg

    integer :: n_lines
    integer :: i

    ! Only allow master to print to screen
    if (.not. master) return

    write(ou, fmt='(1X,A9)', advance='no') 'WARNING: '

    n_lines = (len_trim(msg)-1)/70 + 1
    do i = 1, n_lines
       if (i == 1) then
          write(ou, fmt='(A70)') msg(70*(i-1)+1:70*i)
       else
          write(ou, fmt='(10X,A70)') msg(70*(i-1)+1:70*i)
       end if
    end do

  end subroutine warning

!===============================================================================
! FATAL_ERROR alerts the user that an error has been encountered and displays a
! message about the particular problem. Errors are considered 'fatal' and hence
! the program is aborted.
!===============================================================================

  subroutine fatal_error(msg)

    character(*), intent(in) :: msg

    integer :: n_lines
    integer :: i

    ! Only allow master to print to screen
    if (master) then
       write(eu, fmt='(1X,A7)', advance='no') 'ERROR: '

       n_lines = (len_trim(msg)-1)/72 + 1
       do i = 1, n_lines
          if (i == 1) then
             write(eu, fmt='(A72)') msg(72*(i-1)+1:72*i)
          else
             write(eu, fmt='(7X,A72)') msg(72*(i-1)+1:72*i)
          end if
       end do
       write(eu,*)
    end if

    ! All processors abort
    call free_memory()

  end subroutine fatal_error

end module error
