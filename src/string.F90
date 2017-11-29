module string

  use, intrinsic :: ISO_C_BINDING

  use constants, only: MAX_WORDS, MAX_LINE_LEN, ERROR_INT, ERROR_REAL, &
       OP_LEFT_PAREN, OP_RIGHT_PAREN, OP_COMPLEMENT, OP_INTERSECTION, OP_UNION
  use error,     only: fatal_error, warning
  use stl_vector, only: VectorInt

  implicit none

  interface to_str
    module procedure int4_to_str, int8_to_str, real_to_str
  end interface

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
    character(*), intent(out) :: words(MAX_WORDS)
    integer,      intent(out) :: n

    character(1)  :: chr     ! current character
    integer       :: i       ! current index
    integer       :: i_start ! starting index of word
    integer       :: i_end   ! ending index of word

    i_start = 0
    i_end = 0
    n = 0
    do i = 1, len_trim(string)
      chr = string(i:i)

      ! Note that ACHAR(9) is a horizontal tab
      if ((i_start == 0) .and. (chr /= ' ') .and. (chr /= achar(9))) then
        i_start = i
      end if
      if (i_start > 0) then
        if ((chr == ' ') .or. (chr == achar(9))) i_end = i - 1
        if (i == len_trim(string))   i_end = i
        if (i_end > 0) then
          n = n + 1
          if (i_end - i_start + 1 > len(words(n))) then
            call warning("The word '" // string(i_start:i_end) &
                 // "' is longer than the space allocated for it.")
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
! TOKENIZE takes a string that includes logical expressions for a list of
! bounding surfaces in a cell and splits it into separate tokens. The characters
! (, ), |, and ~ count as separate tokens since they represent operators.
!===============================================================================

  subroutine tokenize(string, tokens)
    character(*), intent(in) :: string
    type(VectorInt), intent(inout) :: tokens

    integer :: i       ! current index
    integer :: i_start ! starting index of word
    integer :: token
    character(len=len_trim(string)) :: string_

    ! Remove leading blanks
    string_ = adjustl(string)

    i_start = 0
    i = 1
    do while (i <= len_trim(string_))
      ! Check for special characters
      if (index('()|~ ', string_(i:i)) > 0) then
        ! If the special character appears immediately after a non-operator,
        ! create a token with the surface half-space
        if (i_start > 0) then
          call tokens%push_back(int(str_to_int(&
               string_(i_start:i - 1)), 4))
        end if

        select case (string_(i:i))
        case ('(')
          call tokens%push_back(OP_LEFT_PAREN)
        case (')')
          if (tokens%size() > 0) then
            token = tokens%data(tokens%size())
            if (token >= OP_UNION .and. token < OP_RIGHT_PAREN) then
              call fatal_error("Right parentheses cannot follow an operator in &
                   &region specification: " // trim(string))
            end if
          end if
          call tokens%push_back(OP_RIGHT_PAREN)
        case ('|')
          if (tokens%size() > 0) then
            token = tokens%data(tokens%size())
            if (.not. (token < OP_UNION .or. token == OP_RIGHT_PAREN)) then
              call fatal_error("Union cannot follow an operator in region &
                   &specification: " // trim(string))
            end if
          end if
          call tokens%push_back(OP_UNION)
        case ('~')
          call tokens%push_back(OP_COMPLEMENT)
        case (' ')
          ! Find next non-space character
          do while (string_(i+1:i+1) == ' ')
            i = i + 1
          end do

          ! If previous token is a halfspace or right parenthesis and next token
          ! is not a left parenthese or union operator, that implies that the
          ! whitespace is to be interpreted as an intersection operator
          if (i_start > 0 .or. tokens%data(tokens%size()) == OP_RIGHT_PAREN) then
            if (index(')|', string_(i+1:i+1)) == 0) then
              call tokens%push_back(OP_INTERSECTION)
            end if
          end if
        end select

        i_start = 0
      else
        ! Check for invalid characters
        if (index('-+0123456789', string_(i:i)) == 0) then
          call fatal_error("Invalid character '" // string_(i:i) // "' in &
               &region specification.")
        end if

        ! If we haven't yet reached the start of a word, start a new word
        if (i_start == 0) i_start = i
      end if

      i = i + 1
    end do

    ! If we've reached the end and we're still in a word, create a token from it
    ! and add it to the list
    if (i_start > 0) then
      call tokens%push_back(int(str_to_int(&
           string_(i_start:len_trim(string_))), 4))
    end if
  end subroutine tokenize

!===============================================================================
! CONCATENATE takes an array of words and concatenates them together in one
! string with a single space between words
!
! Arguments:
!   words   = array of words
!   n_words = total number of words
!   string  = concatenated string
!===============================================================================

  pure function concatenate(words, n_words) result(string)

    integer,        intent(in)  :: n_words
    character(*),   intent(in)  :: words(n_words)
    character(MAX_LINE_LEN)     :: string

    integer :: i ! index

    string = words(1)
    if (n_words == 1) return
    do i = 2, n_words
      string = trim(string) // ' ' // words(i)
    end do

  end function concatenate

!===============================================================================
! TO_LOWER converts a string to all lower case characters
!===============================================================================

  pure function to_lower(word) result(word_lower)
    character(*), intent(in) :: word
    character(len=len(word)) :: word_lower

    integer :: i
    integer :: ic

    do i = 1, len(word)
      ic = ichar(word(i:i))
      if (ic >= 65 .and. ic <= 90) then
        word_lower(i:i) = char(ic+32)
      else
        word_lower(i:i) = word(i:i)
      end if
    end do

  end function to_lower

!===============================================================================
! TO_UPPER converts a string to all upper case characters
!===============================================================================

  pure function to_upper(word) result(word_upper)
    character(*), intent(in) :: word
    character(len=len(word)) :: word_upper

    integer :: i
    integer :: ic

    do i = 1, len(word)
      ic = ichar(word(i:i))
      if (ic >= 97 .and. ic <= 122) then
        word_upper(i:i) = char(ic-32)
      else
        word_upper(i:i) = word(i:i)
      end if
    end do

  end function to_upper

!===============================================================================
! ZERO_PADDED returns a string of the input integer padded with zeros to the
! desired number of digits.  Do not include the sign in n_digits for negative
! integers.
!===============================================================================

  function zero_padded(num, n_digits) result(str)
    integer, intent(in) :: num
    integer, intent(in) :: n_digits
    character(11)       :: str

    character(8)        :: zp_form

    ! Make sure n_digits is reasonable. 10 digits is the maximum needed for the
    ! largest integer(4).
    if (n_digits > 10) then
      call fatal_error('zero_padded called with an unreasonably large &
           &n_digits (>10)')
    end if

    ! Write a format string of the form '(In.m)' where n is the max width and
    ! m is the min width.  If a sign is present, then n must be one greater
    ! than m.
    if (num < 0) then
      write(zp_form, '("(I", I0, ".", I0, ")")') n_digits+1, n_digits
    else
      write(zp_form, '("(I", I0, ".", I0, ")")') n_digits, n_digits
    end if

    ! Format the number.
    write(str, zp_form) num
  end function zero_padded

!===============================================================================
! IS_NUMBER determines whether a string of characters is all 0-9 characters
!===============================================================================

  pure function is_number(word) result(number)
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

!===============================================================================
! STARTS_WITH determines whether a string starts with a certain
! sequence of characters
!===============================================================================

  pure logical function starts_with(str, seq)
    character(*), intent(in) :: str ! string to check
    character(*), intent(in) :: seq ! sequence of characters

    integer :: i
    integer :: i_start
    integer :: str_len
    integer :: seq_len

    str_len = len_trim(str)
    seq_len = len_trim(seq)

    ! determine how many spaces are at beginning of string
    i_start = 0
    do i = 1, str_len
      if (str(i:i) == ' ' .or. str(i:i) == achar(9)) cycle
      i_start = i
      exit
    end do

    ! Check if string starts with sequence using INDEX intrinsic
    if (index(str(1:str_len), seq(1:seq_len)) == i_start) then
      starts_with = .true.
    else
      starts_with = .false.
    end if

  end function starts_with

!===============================================================================
! ENDS_WITH determines whether a string ends with a certain sequence
! of characters
!===============================================================================

  pure logical function ends_with(str, seq)
    character(*), intent(in) :: str ! string to check
    character(*), intent(in) :: seq ! sequence of characters

    integer :: i_start
    integer :: str_len
    integer :: seq_len

    str_len = len_trim(str)
    seq_len = len_trim(seq)

    ! determine how many spaces are at beginning of string
    i_start = str_len - seq_len + 1

    ! Check if string starts with sequence using INDEX intrinsic
    if (index(str(1:str_len), seq(1:seq_len), .true.) == i_start) then
      ends_with = .true.
    else
      ends_with = .false.
    end if

  end function ends_with

!===============================================================================
! COUNT_DIGITS returns the number of digits needed to represent the input
! integer.
!===============================================================================

  pure function count_digits(num) result(n_digits)
    integer, intent(in) :: num
    integer             :: n_digits

    n_digits = 1
    do while (num / 10**(n_digits) /= 0 .and. abs(num / 10 **(n_digits-1)) /= 1&
              &.and. n_digits /= 10)
      ! Note that 10 digits is the maximum needed to represent an integer(4) so
      ! the loop automatically exits when n_digits = 10.
      n_digits = n_digits + 1
    end do

  end function count_digits

!===============================================================================
! INT4_TO_STR converts an integer(4) to a string.
!===============================================================================

  pure function int4_to_str(num) result(str)

    integer, intent(in) :: num
    character(11) :: str

    write (str, '(I11)') num
    str = adjustl(str)

  end function int4_to_str

!===============================================================================
! INT8_TO_STR converts an integer(8) to a string.
!===============================================================================

  pure function int8_to_str(num) result(str)

    integer(8), intent(in) :: num
    character(21) :: str

    write (str, '(I21)') num
    str = adjustl(str)

  end function int8_to_str

!===============================================================================
! STR_TO_INT converts a string to an integer.
!===============================================================================

  pure function str_to_int(str) result(num)

    character(*), intent(in) :: str
    integer(8) :: num

    character(5) :: fmt
    integer      :: w
    integer      :: ioError

    ! Determine width of string
    w = len_trim(str)

    ! Create format specifier for reading string
    write(UNIT=fmt, FMT='("(I",I2,")")') w

    ! read string into integer
    read(UNIT=str, FMT=fmt, IOSTAT=ioError) num
    if (ioError > 0) num = ERROR_INT

  end function str_to_int

!===============================================================================
! STR_TO_REAL converts an arbitrary string to a real(8)
!===============================================================================

  pure function str_to_real(string) result(num)

    character(*), intent(in) :: string
    real(8)                  :: num

    integer :: ioError

    ! Read string
    read(UNIT=string, FMT=*, IOSTAT=ioError) num
    if (ioError > 0) num = ERROR_REAL

  end function str_to_real

!===============================================================================
! REAL_TO_STR converts a real(8) to a string based on how large the value is and
! how many significant digits are desired. By default, six significants digits
! are used.
!===============================================================================

  pure function real_to_str(num, sig_digits) result(string)

    real(8),           intent(in) :: num        ! number to convert
    integer, optional, intent(in) :: sig_digits ! # of significant digits
    character(15)                 :: string     ! string returned

    integer      :: decimal ! number of places after decimal
    integer      :: width   ! total field width
    real(8)      :: num2    ! absolute value of number
    character(9) :: fmt     ! format specifier for writing number

    ! set default field width
    width = 15

    ! set number of places after decimal
    if (present(sig_digits)) then
      decimal = sig_digits
    else
      decimal = 6
    end if

    ! Create format specifier for writing character
    num2 = abs(num)
    if (num2 == 0.0_8) then
      write(fmt, '("(F",I2,".",I2,")")') width, 1
    elseif (num2 < 1.0e-1_8) then
      write(fmt, '("(ES",I2,".",I2,")")') width, decimal - 1
    elseif (num2 >= 1.0e-1_8 .and. num2 < 1.0_8) then
      write(fmt, '("(F",I2,".",I2,")")') width, decimal
    elseif (num2 >= 1.0_8 .and. num2 < 10.0_8) then
      write(fmt, '("(F",I2,".",I2,")")') width, max(decimal-1, 0)
    elseif (num2 >= 10.0_8 .and. num2 < 100.0_8) then
      write(fmt, '("(F",I2,".",I2,")")') width, max(decimal-2, 0)
    elseif (num2 >= 100.0_8 .and. num2 < 1000.0_8) then
      write(fmt, '("(F",I2,".",I2,")")') width, max(decimal-3, 0)
    elseif (num2 >= 100.0_8 .and. num2 < 10000.0_8) then
      write(fmt, '("(F",I2,".",I2,")")') width, max(decimal-4, 0)
    elseif (num2 >= 10000.0_8 .and. num2 < 100000.0_8) then
      write(fmt, '("(F",I2,".",I2,")")') width, max(decimal-5, 0)
    else
      write(fmt, '("(ES",I2,".",I2,")")') width, decimal - 1
    end if

    ! Write string and left adjust
    write(string, fmt) num
    string = adjustl(string)

  end function real_to_str

!===============================================================================
! WORD_COUNT determines the number of words in a string
!===============================================================================

  function word_count(str) result(n)
    character(*), intent(in) :: str
    integer :: n

    integer :: i
    character(kind=C_CHAR) :: chr
    logical :: inword

    ! Count number of words
    inword = .false.
    n = 0
    do i = 1, len_trim(str)
      chr = str(i:i)
      select case (chr)
      case (' ', C_HORIZONTAL_TAB, C_NEW_LINE, C_CARRIAGE_RETURN)
        if (inword) then
          inword = .false.
          n = n + 1
        end if
      case default
        inword = .true.
      end select
    end do
    if (inword) n = n + 1
  end function word_count

!===============================================================================
! TO_F_STRING takes a null-terminated array of C chars and turns it into a
! deferred-length character string. Yay Fortran 2003!
!===============================================================================

  function to_f_string(c_string) result(f_string)
    character(kind=C_CHAR), intent(in) :: c_string(*)
    character(:), allocatable :: f_string

    integer :: i, n

    ! Determine length of original string
    n = 0
    do while (c_string(n + 1) /= C_NULL_CHAR)
      n = n + 1
    end do

    ! Copy C string character by character
    allocate(character(len=n) :: f_string)
    do i = 1, n
      f_string(i:i) = c_string(i)
    end do
  end function to_f_string

!===============================================================================
! TO_C_STRING takes a space-padded Fortran character and turns it into a
! null-terminated C char array. Yay Fortran 2003!
!===============================================================================

  function to_c_string(f_string) result(c_string)
    character(*), intent(in) :: f_string
    character(kind=C_CHAR) :: c_string(len_trim(f_string) + 1)

    integer :: i, n

    ! Copy Fortran string character by character
    n = len_trim(f_string)
    do i = 1, n
      c_string(i) = f_string(i:i)
    end do
    c_string(n + 1) = C_NULL_CHAR
  end function to_c_string

end module string
