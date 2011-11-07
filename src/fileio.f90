module fileio

  use constants
  use global,    only: MAX_LINE_LEN, MAX_WORD_LEN, MAX_WORDS
  use string,    only: split_string_wl, lower_case

  implicit none

contains

!===============================================================================
! READ_LINE reads a line from a file open on a unit
!===============================================================================

  subroutine read_line(unit, line, ioError)

    integer,      intent(in)  :: unit    ! unit to read from
    character(*), intent(out) :: line    ! line to return
    integer,      intent(out) :: ioError ! error status

    read(UNIT=unit, FMT='(A)', IOSTAT=ioError) line

  end subroutine read_line

!===============================================================================
! GET_NEXT_LINE reads the next line to the file connected on the specified unit
! including any continuation lines. If a line ends in an ampersand, the next
! line is read and its words are appended to the final array
!===============================================================================

  subroutine get_next_line(unit, words, n, ioError)

    integer,      intent(in)  :: unit             ! unit to read from
    character(*), intent(out) :: words(MAX_WORDS) ! words read
    integer,      intent(out) :: n                ! number of words
    integer,      intent(out) :: ioError          ! error status

    character(MAX_LINE_LEN) :: line                   ! single line
    character(MAX_WORD_LEN) :: local_words(MAX_WORDS) ! words on one line
    integer                 :: index                  ! index of words

    index = 0
    do
       ! read line from file
       read(UNIT=unit, FMT='(A100)', IOSTAT=ioError) line

       ! if we're at the end of the file, return
       if (ioError /= 0) return

       ! split a single line into words
       call split_string_wl(line, local_words, n)

       ! if there are no words, we're done
       if (n == 0) exit

       ! Check whether there is a continuation line
       if (local_words(n) == '&') then
          words(index+1:index+n-1) = local_words(1:n-1)
          index = index + n - 1
       else
          words(index+1:index+n) = local_words(1:n)
          index = index + n
          exit
       end if
    end do

    ! set total number of words
    n = index

  end subroutine get_next_line

!===============================================================================
! SKIP_LINES skips 'n_lines' lines from a file open on a unit
!===============================================================================

  subroutine skip_lines(unit, n_lines, ioError)

    integer, intent(in)  :: unit    ! unit to read from
    integer, intent(in)  :: n_lines ! number of lines to skip
    integer, intent(out) :: ioError ! error status 

    integer                 :: i   ! index for number of lines
    character(MAX_LINE_LEN) :: tmp ! single line

    do i = 1, n_lines
       read(UNIT=unit, FMT='(A)', IOSTAT=ioError) tmp
    end do

  end subroutine skip_lines

end module fileio
