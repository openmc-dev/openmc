module string

  use global
  use error, only: warning

  implicit none

contains

!===============================================================================
! SPLIT_STRING takes a sentence and splits it into separate words much like the
! Python string.split() method.
!
! Arguments:
!   string = input line
!   words  = array of words
!   n      = total number of words
!===============================================================================

  subroutine split_string(string, words, n)

    character(*), intent(in)  :: string
    character(*), intent(out) :: words(max_words)
    integer,      intent(out) :: n

    character(1)  :: char    ! current character
    integer       :: i       ! current index
    integer       :: i_start ! starting index of word
    integer       :: i_end   ! ending index of word
    character(max_line_len) :: msg

    i_start = 0
    i_end = 0
    n = 0
    do i = 1, len_trim(string)
       char = string(i:i)

       ! Note that ACHAR(9) is a horizontal tab
       if ((i_start == 0) .and. (char /= ' ') .and. (char /= achar(9))) then
          i_start = i
       end if
       if (i_start > 0) then
          if ((char == ' ') .or. (char == achar(9))) i_end = i - 1
          if (i == len_trim(string))   i_end = i
          if (i_end > 0) then
             n = n + 1
             if (i_end - i_start + 1 > len(words(n))) then
                msg = "The word '" // string(i_start:i_end) // "' is longer than " &
                     & // "the space allocated for it."
                call warning(msg)
             end if
             words(n) = string(i_start:i_end)
             ! reset indices
             i_start = 0
             i_end = 0
          end if
       end if
    end do

  end subroutine split_string

!===============================================================================
! SPLIT_STRING_WL takes a string that includes logical expressions for a list of
! bounding surfaces in a cell and splits it into separate words. The characters
! (, ), :, and # count as separate words since they represent operators.
!
! Arguments:
!   string  = input line
!   words   = array of words
!   n       = number of words
!===============================================================================

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
          if (i_start > 0) then
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

!===============================================================================
! CONCATENATE takes an array of words and concatenates them together in one
! string with a single space between words
!
! Arguments:
!   words   = array of words
!   n_words = total number of words
!   string  = concatenated string
!===============================================================================

  function concatenate(words, n_words) result(string)

    character(*),   intent(in)  :: words(n_words)
    integer,        intent(in)  :: n_words
    character(max_line_len)     :: string

    integer :: i ! index

    string = words(1)
    if (n_words == 1) return
    do i = 2, n_words
       string = trim(string) // ' ' // words(i)
    end do

  end function concatenate

!===============================================================================
! LOWER_CASE converts a string to all lower case characters
!===============================================================================

  elemental subroutine lower_case(word)

    character(*), intent(inout) :: word

    integer :: i
    integer :: ic

    do i = 1,len(word)
       ic = ichar(word(i:i))
       if (ic >= 65 .and. ic < 90) word(i:i) = char(ic+32)
    end do

  end subroutine lower_case

!===============================================================================
! UPPER_CASE converts a string to all upper case characters
!===============================================================================

  elemental subroutine upper_case(word)

    character(*), intent(inout) :: word

    integer :: i
    integer :: ic

    do i = 1,len(word)
       ic = ichar(word(i:i))
       if (ic >= 97 .and. ic < 122) word(i:i) = char(ic-32)
    end do

  end subroutine upper_case

!===============================================================================
! IS_NUMBER determines whether a string of characters is all 0-9 characters
!===============================================================================

  function is_number(word) result(number)

    character(*), intent(in) :: word
    logical                  :: number

    integer :: i
    integer :: ic

    number = .true.
    do i = 1, len_trim(word)
       ic = ichar(word(i:i))
       if (ic < 48 .or. ic >= 58) number = .false.
    end do
    
  end function is_number

end module string
