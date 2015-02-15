undefine(`index')dnl
undefine(`len')dnl
undefine(`format')dnl
dnl
include(`foreach.m4')dnl
dnl
include(`common.m4')dnl
dnl
define(`TOHW_defaultdecls', `dnl
    logical :: bracketed
    integer :: i, j, ij, k, s_i, err, ios, length
    real :: r, c
')dnl
dnl
define(`TOHW_check_errors', `dnl
    num = ij
    if (verify(s(s_i:), whitespace)/=0) num = 0
    if (err/=0) num=0
')dnl
dnl
define(`TOHW_parse_strings_csv', `dnl
      if (s_i>len(s)) then
        ij = ij + 1
        exit loop
      endif
      k = verify(s(s_i:), achar(10)//achar(13))
      if (k==0) then
        ij = ij + 1
        exit loop
      elseif (s(s_i+k-1:s_i+k-1)=="""") then
        ! we have a quote-delimited string;
        s_i = s_i + k
        s2 = ""
        quote: do 
          k = index(s(s_i:), """")
          if (k==0) then
            err = 2
            exit loop
          endif
          k = s_i + k - 1
          s2(m:) = s(s_i:k)
          m = m + (k-s_i+1)
          k = k + 2
          if (k>len(s)) then
            err = 2
            exit loop
          endif
          if (s(k:k)/="""") exit
          s_i = k + 1
          if (s_i > len(s)) then
            err = 2
            exit loop
          endif
          m = m + 1
          s2(m:m) = """"
        enddo quote
        k  = scan(s(s_i:), whitespace)
        if (k==0) then
          err = 2
          exit loop
        endif
      else
        s_i = s_i + k - 1
        k = scan(s(s_i:), achar(10)//achar(13)//",")
        if (k==0) then
          eof = .true.
          k = len(s)
        else
          k = s_i + k - 2
        endif
        if (index(s(s_i:k), """")/=0) then
          err = 2
          exit loop
        endif
      endif
      ij = ij + 1
      s_i = k + 2
      if (eof) exit loop
')dnl
define(`TOHW_parse_strings', `dnl
      if (present(separator)) then
        k = index(s(s_i:), separator)
      else
        k = verify(s(s_i:), whitespace)
        if (k==0) exit loop
        s_i = s_i + k - 1
        k = scan(s(s_i:), whitespace)
      endif
      if (k==0) then
        k = len(s)
      else
        k = s_i + k - 2
      endif
      ij = ij + 1
      s_i = k + 2
      if (s_i>len(s)) exit loop
')dnl
dnl
define(`TOHW_parse_logical', `dnl
      k = verify(s(s_i:), whitespace)
      if (k==0) exit loop
      s_i = s_i + k - 1
      if (s(s_i:s_i)==",") then
        if (s_i+1>len(s)) then
          err = 2
          exit loop
        endif
        k = verify(s(s_i+1:), whitespace)
        s_i = s_i + k - 1
      endif
      k = scan(s(s_i:), whitespace//",")
      if (k==0) then
        k = len(s)
      else
        k = s_i + k - 2
      endif
      if (.not.((s(s_i:k)=="true".or.s(s_i:k)=="1").or. & 
          (s(s_i:k)=="false".or.s(s_i:k)=="0"))) then
        err = 2
        exit loop
      endif
      ij = ij + 1
      s_i = k + 2
      if (ij<length.and.s_i>len(s)) exit loop
')dnl
dnl
define(`TOHW_parse_numbers', `dnl
      k = verify(s(s_i:), whitespace)
      if (k==0) exit loop
      s_i = s_i + k - 1
      if (s(s_i:s_i)==",") then
        if (s_i+1>len(s)) then
          err = 2
          exit loop
        endif
        k = verify(s(s_i+1:), whitespace)
        s_i = s_i + k - 1
      endif
      k = scan(s(s_i:), whitespace//",")
      if (k==0) then
        k = len(s)
      else
        k = s_i + k - 2
      endif
      read(s(s_i:k), *, iostat=ios) dummy_data
      if (ios/=0) then
        err = 2
        exit loop
      endif
      ij = ij + 1
      s_i = k + 2
      if (ij<length.and.s_i>len(s)) exit loop
')dnl
dnl
define(`TOHW_parse_complex', `dnl
      bracketed = .false.
      k = verify(s(s_i:), whitespace)
      if (k==0) exit loop
      s_i = s_i + k - 1
      select case (s(s_i:s_i))
      case ("(")
        bracketed = .true.
        k = verify(s(s_i:), whitespace)
        if (k==0) then
          err = 2
          exit loop
        endif
        s_i = s_i + k
      case (",")
        k = verify(s(s_i:), whitespace)
        if (k==0) then
          err = 2
          exit loop
        endif
        s_i = s_i + k - 1
      case ("+", "-", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9")
        continue
      case default
        err = 2
        exit loop
      end select
      if (bracketed) then
        k = index(s(s_i:), ")+i(")
      else
        k = scan(s(s_i:), whitespace//",")
      endif
      if (k==0) then
        err = 2
        exit loop
      endif
      k = s_i + k - 2
      read(s(s_i:k), *, iostat=ios) r
      if (ios/=0) then
        err = 2
        exit loop
      endif
      if (bracketed) then
        s_i = k + 5
        if (s_i>len(s)) then
          err = 2
          exit loop
        endif
      else
        s_i = k + 2
      endif
      if (bracketed) then
        k = index(s(s_i:), ")")
        if (k==0) then
          err = 2
          exit loop
        endif
        k = s_i + k - 2
      else
        k = scan(s(s_i:), whitespace//",")
        if (k==0) then
          k = len(s)
        else
          k = s_i + k - 2
        endif
      endif
      read(s(s_i:k), *, iostat=ios) c
      if (ios/=0) then
        err = 2
        exit loop
      endif
      ij = ij + 1
      s_i = k + 2
      if (ij<length.and.s_i>len(s)) exit loop
')dnl
dnl
module fox_m_fsys_count_parse_input

  use fox_m_fsys_realtypes, only: sp, dp

  implicit none
  private

  character(len=1), parameter :: SPACE           = achar(32)
  character(len=1), parameter :: NEWLINE         = achar(10)
  character(len=1), parameter :: CARRIAGE_RETURN = achar(13)
  character(len=1), parameter :: TAB             = achar(9)
  character(len=*), parameter :: whitespace = &
    SPACE//NEWLINE//CARRIAGE_RETURN//TAB

  interface countrts
    module procedure countstring
    module procedure countlogical
    module procedure countinteger
    module procedure countrealsp
    module procedure countrealdp
    module procedure countcomplexsp
    module procedure countcomplexdp
  end interface

  public :: countrts

contains

define(`m4f_thisfunc', `countstring')dnl
  pure function m4f_thisfunc`'(s, datatype, separator, csv) result(num)
    character(len=*), intent(in) :: s
    character(len=*), intent(in) :: datatype
    character, intent(in), optional :: separator
    logical, intent(in), optional :: csv
    integer :: num
#ifndef DUMMYLIB
TOHW_defaultdecls
    character(len=len(s)) :: s2
    logical :: csv_, eof
    integer :: m

    csv_ = .false.
    if (present(csv)) then
      if (csv) csv_ = csv
    endif

    s_i = 1
    err = 0
    eof = .false.
    ij = 0
    loop: do
      if (csv_) then
TOHW_parse_strings_csv(`data(i)')
      else
TOHW_parse_strings(`data(i)')
      endif
    end do loop

    num = ij
    if (err/=0) num=0

#else
    num = 0
#endif
  end function m4f_thisfunc

define(`m4f_thisfunc', `countlogical')dnl
  pure function m4f_thisfunc`'(s, datatype) result(num)
    character(len=*), intent(in) :: s
    logical, intent(in) :: datatype
    logical :: dummy_data
    integer :: num
#ifndef DUMMYLIB
TOHW_defaultdecls

    s_i = 1
    err = 0
    ij  = 0
    length = 1
    loop: do 
TOHW_parse_logical(`data(i)')
    end do loop

TOHW_check_errors

#else
    num = 0
#endif
  end function m4f_thisfunc

define(`m4f_thisfunc', `countinteger')dnl
  pure function m4f_thisfunc`'(s, datatype) result(num)
    character(len=*), intent(in) :: s
    integer, intent(in) :: datatype
    integer :: dummy_data
    integer num
#ifndef DUMMYLIB
TOHW_defaultdecls
    s_i = 1
    err = 0
    ij  = 0
    length = 1
    loop: do 
TOHW_parse_numbers(`data(i)')
    end do loop

TOHW_check_errors

#else
    num = 0
#endif
end function m4f_thisfunc

define(`m4f_thisfunc', `countrealsp')dnl
  pure function m4f_thisfunc`'(s, datatype) result(num)
    character(len=*), intent(in) :: s
    real(sp), intent(in) :: datatype
    real(sp) :: dummy_data
    integer :: num
#ifndef DUMMYLIB
TOHW_defaultdecls

    s_i = 1
    err = 0
    ij  = 0
    length = 1
    loop: do
TOHW_parse_numbers(`data(i)')
    end do loop

TOHW_check_errors

#else
    num = 0
#endif
  end function m4f_thisfunc

define(`m4f_thisfunc', `countrealdp')dnl
  pure function m4f_thisfunc`'(s, datatype) result(num)
    character(len=*), intent(in) :: s
    real(dp), intent(in) :: datatype
    real(dp) :: dummy_data
    integer :: num
#ifndef DUMMYLIB
TOHW_defaultdecls

    s_i = 1
    err = 0
    ij  = 0
    length = 1
    loop: do 
TOHW_parse_numbers(`data(i)')
    end do loop

TOHW_check_errors

#else
    num - 0
#endif
  end function m4f_thisfunc

define(`m4f_thisfunc', `countcomplexsp')dnl
  pure function m4f_thisfunc`'(s, datatype) result(num)
    character(len=*), intent(in) :: s
    complex(sp), intent(in) :: datatype
    complex(sp) :: dummy_data
    integer :: num
#ifndef DUMMYLIB
TOHW_defaultdecls

    s_i = 1
    err = 0
    ij  = 0
    length = 1
    loop: do 
TOHW_parse_complex(`data(i)')
    end do loop

TOHW_check_errors

#else
    num = 0
#endif
  end function m4f_thisfunc

define(`m4f_thisfunc', `countcomplexdp')dnl
  pure function m4f_thisfunc`'(s, datatype) result(num)
    character(len=*), intent(in) :: s
    complex(dp), intent(in) :: datatype
    complex(dp) :: dummy_data
    integer :: num
#ifndef DUMMYLIB
TOHW_defaultdecls

    s_i = 1
    err = 0
    ij  = 0
    length = 1
    loop: do 
TOHW_parse_complex(`data(i)')
    end do loop

TOHW_check_errors

#else
    num = 0
#endif
  end function m4f_thisfunc

end module fox_m_fsys_count_parse_input
