module string

  use constants,     only: MAX_WORDS, MAX_LINE_LEN, ERROR_INT, ERROR_REAL, &
                           OP_LEFT_PAREN, OP_RIGHT_PAREN, OP_COMPLEMENT, &
                           OP_INTERSECTION, OP_UNION
  use error,         only: fatal_error, warning
  use global,        only: master
  use simple_string, only: str_to_int, str_to_real
  use stl_vector,    only: VectorInt

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
            if (master) call warning("The word '" // string(i_start:i_end) &
                 &// "' is longer than the space allocated for it.")
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

end module string
