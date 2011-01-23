module string

  use global, only: max_words
  use output, only: error

  implicit none

contains

!=====================================================================
! SPLIT_STRING takes a sentence and splits it into separate words much
! the Python string.split() method.
!
! Arguments:
!   string = input line
!   words  = array of words
!   n      = total number of words
!=====================================================================

  subroutine split_string(string, words, n)

    character(*), intent(in)  :: string
    character(*), intent(out) :: words(max_words)
    integer,      intent(out) :: n

    character(1)  :: char    ! current character
    integer       :: i       ! current index
    integer       :: i_start ! starting index of word
    integer       :: i_end   ! ending index of word

    i_start = 0
    i_end = 0
    n = 0
    do i = 1, len(string)
       char = string(i:i)

       ! Note that ACHAR(9) is a horizontal tab
       if ((i_start == 0) .and. (char /= ' ') .and. (char /= achar(9))) then
          i_start = i
       end if
       if (i_start > 0) then
          if ((char == ' ') .or. (char == achar(9))) i_end = i - 1
          if (i == len(string))   i_end = i
          if (i_end > 0) then
             n = n + 1
             words(n) = string(i_start:i_end)
             ! reset indices
             i_start = 0
             i_end = 0
          end if
       end if
    end do

  end subroutine split_string

!=====================================================================
! SPLIT_STRING_WL takes a string that includes logical expressions for
! a list of bounding surfaces in a cell and splits it into separate
! words. The characters (, ), :, and # count as separate words since
! they represent operators.
!
! Arguments:
!   string  = input line
!   words   = array of words
!   n       = number of words
!=====================================================================

  subroutine split_string_wl(string, words, n)

    character(*), intent(in)  :: string
    character(*), intent(out) :: words(max_words)
    integer,      intent(out) :: n

    character(1)  :: char    ! current character
    integer       :: i       ! current index
    integer       :: i_start ! starting index of word
    integer       :: i_end   ! ending index of word

    i_start = 0
    i_end = 0
    n = 0
    do i = 1, len_trim(string)
       char = string(i:i)

       ! Check for special characters
       if (index('():#', char) > 0) then
          if (i_start > 0 ) then
             i_end = i - 1
             n = n + 1
             words(n) = string(i_start:i_end)
          end if
          n = n + 1
          words(n) = char
          i_start = 0
          i_end = 0
          cycle
       end if

       if ((i_start == 0) .and. (char /= ' ')) then
          i_start = i
       end if
       if (i_start > 0) then
          if (char == ' ')           i_end = i - 1
          if (i == len_trim(string)) i_end = i
          if (i_end > 0) then
             n = n + 1
             words(n) = string(i_start:i_end)
             ! reset indices
             i_start = 0
             i_end = 0
          end if
       end if
    end do
  end subroutine split_string_wl

!=====================================================================
! CONCATENATE takes an array of words and concatenates them
! together in one string with a single space between words
!
! Arguments:
!   words   = array of words
!   n_words = total number of words
!   string  = concatenated string
!=====================================================================

  subroutine concatenate( words, n_words, string )

    character(*),   intent(in)  :: words(n_words)
    integer,        intent(in)  :: n_words
    character(250), intent(out) :: string

    integer :: i ! index

    string = words(1)
    if ( n_words == 1 ) return
    do i = 2, n_words
       string = trim(string) // ' ' // words(i)
    end do

  end subroutine concatenate

!=====================================================================
! LOWER_CASE converts a string to all lower case characters
!=====================================================================

  elemental subroutine lower_case(word)
    ! convert a word to lower case

    character(*), intent(inout) :: word

    integer :: i
    integer :: ic

    do i = 1,len(word)
       ic = ichar(word(i:i))
       if (ic >= 65 .and. ic < 90) word(i:i) = char(ic+32)
    end do

  end subroutine lower_case

!=====================================================================
! STR_TO_REAL converts an arbitrary string to a real(8). Generally
! this function is intended for strings for which the exact format is
! not known. If the format of the number is known a priori, the
! appropriate format descriptor should be used in lieu of this routine
! because of the extra overhead.
!
! Arguments:
!   string = character(*) containing number to convert
!=====================================================================

  function str_to_real( string )

    character(*), intent(in) :: string
    real(8) :: str_to_real

    integer :: index_decimal  ! index of decimal point
    integer :: index_exponent ! index of exponent character
    integer :: w              ! total field width
    integer :: d              ! number of digits to right of decimal point
    integer :: readError

    character(8) :: fmt   ! format for reading string
    character(250) :: msg ! error message

    ! Determine total field width
    w = len_trim(string)

    ! Determine number of digits to right of decimal point
    index_decimal = index(string, '.')
    index_exponent = max(index(string, 'd'), index(string, 'D'), &
         & index(string, 'e'), index(string, 'E'))
    if ( index_decimal > 0 ) then
       if ( index_exponent > 0 ) then
          d = index_exponent - index_decimal - 1
       else
          d = w - index_decimal
       end if
    else
       d = 0
    end if

    ! Create format specifier for reading string
    write( fmt, '("(E",I2,".",I2,")")') w, d

    ! Read string
    read( string, fmt, iostat=readError ) str_to_real
    if ( readError > 0 ) then
       msg = "Could not read value: " // string
       call error( msg )
    end if

  end function str_to_real

end module string
