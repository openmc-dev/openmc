module error

  use, intrinsic :: ISO_FORTRAN_ENV

  use global

#ifdef MPI
  use mpi
#endif

  implicit none

contains

!===============================================================================
! WARNING issues a warning to the user in the log file and the standard output
! stream.
!===============================================================================

  subroutine warning()

    integer :: i_start   ! starting position
    integer :: i_end     ! ending position
    integer :: line_wrap ! length of line
    integer :: length    ! length of message
    integer :: indent    ! length of indentation

    ! Only allow master to print to screen
    if (.not. master) return

    ! Write warning at beginning
    write(OUTPUT_UNIT, fmt='(1X,A)', advance='no') 'WARNING: '

    ! Set line wrapping and indentation
    line_wrap = 80
    indent = 10

    ! Determine length of message
    length = len_trim(message)

    i_start = 0
    do
      if (length - i_start < line_wrap - indent + 1) then
        ! Remainder of message will fit on line
        write(OUTPUT_UNIT, fmt='(A)') message(i_start+1:length)
        exit

      else
        ! Determine last space in current line
        i_end = i_start + index(message(i_start+1:i_start+line_wrap-indent+1), &
             ' ', BACK=.true.)

        if (i_end == i_start) then
          ! This is a special case where there is no space
          i_end = i_start + line_wrap - indent + 1
          write(ERROR_UNIT, fmt='(A/A)', advance='no') &
               message(i_start+1:i_end-1), repeat(' ', indent)
          i_end = i_end - 1
        else
          ! Write up to last space
          write(OUTPUT_UNIT, fmt='(A/A)', advance='no') &
               message(i_start+1:i_end-1), repeat(' ', indent)
        end if

        ! Advance starting position
        i_start = i_end
        if (i_start > length) exit
      end if
    end do

  end subroutine warning

!===============================================================================
! FATAL_ERROR alerts the user that an error has been encountered and displays a
! message about the particular problem. Errors are considered 'fatal' and hence
! the program is aborted.
!===============================================================================

  subroutine fatal_error(error_code)

    integer, optional :: error_code ! error code

    integer :: code      ! error code
    integer :: i_start   ! starting position
    integer :: i_end     ! ending position
    integer :: line_wrap ! length of line
    integer :: length    ! length of message
    integer :: indent    ! length of indentation


    ! set default error code
    if (present(error_code)) then
      code = error_code
    else
      code = -1
    end if

    ! Write error at beginning
    write(ERROR_UNIT, fmt='(1X,A)', advance='no') 'ERROR: '

    ! Set line wrapping and indentation
    line_wrap = 80
    indent = 8

    ! Determine length of message
    length = len_trim(message)

    i_start = 0
    do
      if (length - i_start < line_wrap - indent + 1) then
        ! Remainder of message will fit on line
        write(ERROR_UNIT, fmt='(A)') message(i_start+1:length)
        exit

      else
        ! Determine last space in current line
        i_end = i_start + index(message(i_start+1:i_start+line_wrap-indent+1), &
             ' ', BACK=.true.)

        if (i_end == i_start) then
          ! This is a special case where there is no space
          i_end = i_start + line_wrap - indent + 1
          write(ERROR_UNIT, fmt='(A/A)', advance='no') &
               message(i_start+1:i_end-1), repeat(' ', indent)
          i_end = i_end - 1
        else
          ! Write up to last space
          write(ERROR_UNIT, fmt='(A/A)', advance='no') &
               message(i_start+1:i_end-1), repeat(' ', indent)
        end if

        ! Advance starting position
        i_start = i_end
        if (i_start > length) exit
      end if
    end do

    ! Write information on current batch, generation, and particle
    if (current_batch > 0) then
      write(ERROR_UNIT,'(1X,A,I12) ') 'Batch:     ', current_batch
      write(ERROR_UNIT,'(1X,A,I12) ') 'Generation:', current_gen
      write(ERROR_UNIT,'(1X,A,I12)')  'Particle:  ', p % id
      write(ERROR_UNIT,'(1X,A,3ES12.4)') 'Location:  ', p % coord0 % xyz
      write(ERROR_UNIT,'(1X,A,3ES12.4)') 'Direction: ', p % coord0 % uvw
      write(ERROR_UNIT,*)
    end if

    ! Release memory from all allocatable arrays
    call free_memory()

#ifdef MPI
    ! Abort MPI
    call MPI_ABORT(MPI_COMM_WORLD, code, mpi_err)
#endif

    ! Abort program
#ifdef NO_F2008
    stop
#else
    error stop 
#endif

  end subroutine fatal_error

end module error
