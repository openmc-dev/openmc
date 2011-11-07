module fileio

  use constants
  use global,    only: MAX_LINE_LEN, MAX_WORD_LEN, MAX_WORDS
  use string,    only: split_string_wl, lower_case

  implicit none

!===============================================================================
! READ_DATA interface allows data to be read with one function regardless of
! whether it is integer or real data. E.g. NXS and JXS can be read with the
! integer version and XSS can be read with the real version
!===============================================================================

  interface read_data
     module procedure read_data_int, read_data_real
  end interface

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

!===============================================================================
! READ_DATA_INT reads integer data into an array from a file open
!===============================================================================

  subroutine read_data_int(unit, array, n, lines, words_per_line)

    integer, intent(in)  :: unit           ! unit to read from
    integer, intent(in)  :: n              ! total number of ints
    integer, intent(out) :: array(n)       ! ints read from file
    integer, intent(in)  :: lines          ! total number of lines
    integer, intent(in)  :: words_per_line ! number of words per line

    integer :: i   ! line index
    integer :: loc ! locator for array

    loc = 0
    do i = 1, lines
       if (i == lines) then
          read(UNIT=unit,FMT=*) array(loc+1:n)
       else
          read(UNIT=unit,FMT=*) array(loc+1:loc+words_per_line)
          loc = loc + words_per_line
       end if
    end do

  end subroutine read_data_int

!===============================================================================
! READ_DATA_REAL reads real(8) data into an array from a file open
!===============================================================================

  subroutine read_data_real(unit, array, n, lines, words_per_line)

    integer, intent(in)  :: unit           ! unit to read from
    integer, intent(in)  :: n              ! total number of ints
    real(8), intent(out) :: array(n)       ! real(8)s read from file
    integer, intent(in)  :: lines          ! total number of lines
    integer, intent(in)  :: words_per_line ! number of words per line

    integer :: i   ! line index
    integer :: loc ! locator for array

    loc = 0
    do i = 1, lines
       if (i == lines) then
          read(UNIT=unit,FMT=*) array(loc+1:n)
       else
          read(UNIT=unit,FMT=*) array(loc+1:loc+words_per_line)
          loc = loc + words_per_line
       end if
    end do

  end subroutine read_data_real

end module fileio
