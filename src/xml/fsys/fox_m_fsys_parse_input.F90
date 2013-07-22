module fox_m_fsys_parse_input

  use fox_m_fsys_realtypes, only: sp, dp

  implicit none
  private

  character(len=1), parameter :: SPACE           = achar(32)
  character(len=1), parameter :: NEWLINE         = achar(10)
  character(len=1), parameter :: CARRIAGE_RETURN = achar(13)
  character(len=1), parameter :: TAB             = achar(9)
  character(len=*), parameter :: whitespace = &
    SPACE//NEWLINE//CARRIAGE_RETURN//TAB

  interface rts
    module procedure scalartostring
    module procedure scalartological
    module procedure scalartointeger
    module procedure scalartorealsp
    module procedure scalartorealdp
    module procedure scalartocomplexsp
    module procedure scalartocomplexdp
    module procedure arraytostring
    module procedure arraytological
    module procedure arraytointeger
    module procedure arraytorealsp
    module procedure arraytorealdp
    module procedure arraytocomplexsp
    module procedure arraytocomplexdp
    module procedure matrixtostring
    module procedure matrixtological
    module procedure matrixtointeger
    module procedure matrixtorealsp
    module procedure matrixtorealdp
    module procedure matrixtocomplexsp
    module procedure matrixtocomplexdp
  end interface

  public :: rts

contains

  subroutine scalartostring(s, data, separator, csv, num, iostat)
    character(len=*), intent(in) :: s
    character(len=*), intent(out) :: data
    character, intent(in), optional :: separator
    logical, intent(in), optional :: csv
    integer, intent(out), optional :: num, iostat
#ifndef DUMMYLIB
    logical :: bracketed
    integer :: i, j, ij, k, s_i, err, ios, length
    real :: r, c

    character(len=len(s)) :: s2
    logical :: csv_, eof, sp
    integer :: m

    csv_ = .false.
    if (present(csv)) then
      csv_ = csv
    endif

    s_i = 1
    err = 0
    eof = .false.
    data = ""
    ij = 0
    length = 1
    loop: do
      if (csv_) then
              if (s_i>len(s)) then
        data = ""
        ij = ij + 1
        exit loop
      endif
      k = verify(s(s_i:), achar(10)//achar(13))
      if (k==0) then
        data = ""
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
        data = s2
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
          if (ij+1==length.and.s(s_i+k-1:s_i+k-1)==",") err = 1
          k = s_i + k - 2
        endif
        data = s(s_i:k)
        if (index(data, """")/=0) then
          err = 2
          exit loop
        endif
      endif
      ij = ij + 1
      s_i = k + 2
      if (eof) exit loop

    if (present(num)) num = ij
    if (ij<length) then
      if (err==0) err = -1
    else
      if (verify(s(s_i:), whitespace)/=0) err = 1
    endif


    if (present(iostat)) then
      iostat = err
    else
      select case (err)
      case(-1)
        write(0, *) "Error in scalartostring"
        write(0, *) "Too few elements found"
        stop
      case(1)
        write(0, *) "Error in scalartostring"
        write(0, *) "Too many elements found"
        stop
      case(2)
        write(0, *) "Error in scalartostring"
        write(0, *) "Malformed input"
        stop
      end select
    end if

      else
        sp = .true.
        do i = 1, len(s)
          if (s_i.gt.len(data)) exit ! Truncate long input... 
                                     ! should we set iostat to 1?
                                     ! probably not - it will break the tests.
          if (sp) then
            if (verify(s(i:i), whitespace)/=0) then
              data(s_i:s_i) = s(i:i)
              s_i = s_i + 1
              sp = .false.
            endif
          else
            if (verify(s(i:i), whitespace)==0) then
              data(s_i:s_i) = " "
              sp = .true.
            else
              data(s_i:s_i) = s(i:i)
            endif
            s_i = s_i + 1
          endif
        enddo
        if (present(num)) num = 1
        if (present(iostat)) iostat = 0
      endif
      exit
    enddo loop

#else
    data = ""
#endif
  end subroutine scalartostring

  subroutine scalartological(s, data, num, iostat)
    character(len=*), intent(in) :: s
    logical, intent(out) :: data
    integer, intent(out), optional :: num, iostat
#ifndef DUMMYLIB 
    logical :: bracketed
    integer :: i, j, ij, k, s_i, err, ios, length
    real :: r, c


    s_i = 1
    err = 0
    data = .false.
    ij = 0
    length = 1
    loop: do i = 1, 1
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
      if (s(s_i:k)=="true".or.s(s_i:k)=="1") then
        data = .true.
      elseif (s(s_i:k)=="false".or.s(s_i:k)=="0") then
        data = .false.
      else
        err = 2
        exit loop
      endif
      ij = ij + 1
      s_i = k + 2
      if (ij<length.and.s_i>len(s)) exit loop

    end do loop

    if (present(num)) num = ij
    if (ij<length) then
      if (err==0) err = -1
    else
      if (verify(s(s_i:), whitespace)/=0) err = 1
    endif


    if (present(iostat)) then
      iostat = err
    else
      select case (err)
      case(-1)
        write(0, *) "Error in scalartological"
        write(0, *) "Too few elements found"
        stop
      case(1)
        write(0, *) "Error in scalartological"
        write(0, *) "Too many elements found"
        stop
      case(2)
        write(0, *) "Error in scalartological"
        write(0, *) "Malformed input"
        stop
      end select
    end if

#else
    data = .false.
#endif
  end subroutine scalartological

  subroutine scalartointeger(s, data, num, iostat)
    character(len=*), intent(in) :: s
    integer, intent(out) :: data
    integer, intent(out), optional :: num, iostat
#ifndef DUMMYLIB
    logical :: bracketed
    integer :: i, j, ij, k, s_i, err, ios, length
    real :: r, c


    s_i = 1
    err = 0
    data = 0
    ij = 0
    length = 1
    loop: do i = 1, 1
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
      read(s(s_i:k), *, iostat=ios) data
      if (ios/=0) then
        err = 2
        exit loop
      endif
      ij = ij + 1
      s_i = k + 2
      if (ij<length.and.s_i>len(s)) exit loop

    end do loop

    if (present(num)) num = ij
    if (ij<length) then
      if (err==0) err = -1
    else
      if (verify(s(s_i:), whitespace)/=0) err = 1
    endif


    if (present(iostat)) then
      iostat = err
    else
      select case (err)
      case(-1)
        write(0, *) "Error in scalartointeger"
        write(0, *) "Too few elements found"
        stop
      case(1)
        write(0, *) "Error in scalartointeger"
        write(0, *) "Too many elements found"
        stop
      case(2)
        write(0, *) "Error in scalartointeger"
        write(0, *) "Malformed input"
        stop
      end select
    end if

#else
    data = 0
#endif
  end subroutine scalartointeger

  subroutine scalartorealsp(s, data, num, iostat)
    character(len=*), intent(in) :: s
    real(sp), intent(out) :: data
    integer, intent(out), optional :: num, iostat
#ifndef DUMMYLIB
    logical :: bracketed
    integer :: i, j, ij, k, s_i, err, ios, length
    real :: r, c


    s_i = 1
    err = 0
    data = 0.0_sp
    ij = 0
    length = 1
    loop: do i = 1, 1
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
      read(s(s_i:k), *, iostat=ios) data
      if (ios/=0) then
        err = 2
        exit loop
      endif
      ij = ij + 1
      s_i = k + 2
      if (ij<length.and.s_i>len(s)) exit loop

    end do loop

    if (present(num)) num = ij
    if (ij<length) then
      if (err==0) err = -1
    else
      if (verify(s(s_i:), whitespace)/=0) err = 1
    endif


    if (present(iostat)) then
      iostat = err
    else
      select case (err)
      case(-1)
        write(0, *) "Error in scalartorealsp"
        write(0, *) "Too few elements found"
        stop
      case(1)
        write(0, *) "Error in scalartorealsp"
        write(0, *) "Too many elements found"
        stop
      case(2)
        write(0, *) "Error in scalartorealsp"
        write(0, *) "Malformed input"
        stop
      end select
    end if

#else
    data = 0.0_sp
#endif
  end subroutine scalartorealsp

  subroutine scalartorealdp(s, data, num, iostat)
    character(len=*), intent(in) :: s
    real(dp), intent(out) :: data
    integer, intent(out), optional :: num, iostat
#ifndef DUMMYLIB
    logical :: bracketed
    integer :: i, j, ij, k, s_i, err, ios, length
    real :: r, c


    s_i = 1
    err = 0
    data = 0.0_dp
    ij = 0
    length = 1
    loop: do i = 1, 1
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
      read(s(s_i:k), *, iostat=ios) data
      if (ios/=0) then
        err = 2
        exit loop
      endif
      ij = ij + 1
      s_i = k + 2
      if (ij<length.and.s_i>len(s)) exit loop

    end do loop

    if (present(num)) num = ij
    if (ij<length) then
      if (err==0) err = -1
    else
      if (verify(s(s_i:), whitespace)/=0) err = 1
    endif


    if (present(iostat)) then
      iostat = err
    else
      select case (err)
      case(-1)
        write(0, *) "Error in scalartorealdp"
        write(0, *) "Too few elements found"
        stop
      case(1)
        write(0, *) "Error in scalartorealdp"
        write(0, *) "Too many elements found"
        stop
      case(2)
        write(0, *) "Error in scalartorealdp"
        write(0, *) "Malformed input"
        stop
      end select
    end if

#else
    data = 0.0_dp
#endif
  end subroutine scalartorealdp

  subroutine scalartocomplexsp(s, data, num, iostat)
    character(len=*), intent(in) :: s
    complex(sp), intent(out) :: data
    integer, intent(out), optional :: num, iostat
#ifndef DUMMYLIB
    logical :: bracketed
    integer :: i, j, ij, k, s_i, err, ios, length
    real :: r, c


    s_i = 1
    err = 0
    data = 0.0_sp
    ij = 0
    length = 1
    loop: do i = 1, 1
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
      data = cmplx(r, c)
      ij = ij + 1
      s_i = k + 2
      if (ij<length.and.s_i>len(s)) exit loop

    end do loop

    if (present(num)) num = ij
    if (ij<length) then
      if (err==0) err = -1
    else
      if (verify(s(s_i:), whitespace)/=0) err = 1
    endif


    if (present(iostat)) then
      iostat = err
    else
      select case (err)
      case(-1)
        write(0, *) "Error in scalartocomplexsp"
        write(0, *) "Too few elements found"
        stop
      case(1)
        write(0, *) "Error in scalartocomplexsp"
        write(0, *) "Too many elements found"
        stop
      case(2)
        write(0, *) "Error in scalartocomplexsp"
        write(0, *) "Malformed input"
        stop
      end select
    end if

#else
    data = (0.0_sp, 0.0_sp)
#endif
  end subroutine scalartocomplexsp

  subroutine scalartocomplexdp(s, data, num, iostat)
    character(len=*), intent(in) :: s
    complex(dp), intent(out) :: data
    integer, intent(out), optional :: num, iostat
#ifndef DUMMYLIB
    logical :: bracketed
    integer :: i, j, ij, k, s_i, err, ios, length
    real :: r, c


    s_i = 1
    err = 0
    data = 0.0_dp
    ij = 0
    length = 1
    loop: do i = 1, 1
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
      data = cmplx(r, c)
      ij = ij + 1
      s_i = k + 2
      if (ij<length.and.s_i>len(s)) exit loop

    end do loop

    if (present(num)) num = ij
    if (ij<length) then
      if (err==0) err = -1
    else
      if (verify(s(s_i:), whitespace)/=0) err = 1
    endif


    if (present(iostat)) then
      iostat = err
    else
      select case (err)
      case(-1)
        write(0, *) "Error in scalartocomplexdp"
        write(0, *) "Too few elements found"
        stop
      case(1)
        write(0, *) "Error in scalartocomplexdp"
        write(0, *) "Too many elements found"
        stop
      case(2)
        write(0, *) "Error in scalartocomplexdp"
        write(0, *) "Malformed input"
        stop
      end select
    end if

#else
    data = (0.0_dp, 0.0_dp)
#endif
  end subroutine scalartocomplexdp

  subroutine arraytostring(s, data, separator, csv, num, iostat)
    character(len=*) :: data(:)
    character, intent(in), optional :: separator
    logical, intent(in), optional :: csv
    character(len=*), intent(in) :: s
    integer, intent(out), optional :: num
    integer, intent(out), optional :: iostat

#ifndef DUMMYLIB
    logical :: bracketed
    integer :: i, j, ij, k, s_i, err, ios, length
    real :: r, c

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
    data = ""
    ij = 0
    length = size(data)
    loop: do i = 1, size(data)
      if (csv_) then
      if (s_i>len(s)) then
        data(i) = ""
        ij = ij + 1
        exit loop
      endif
      k = verify(s(s_i:), achar(10)//achar(13))
      if (k==0) then
        data(i) = ""
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
        data(i) = s2
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
          if (ij+1==length.and.s(s_i+k-1:s_i+k-1)==",") err = 1
          k = s_i + k - 2
        endif
        data(i) = s(s_i:k)
        if (index(data(i), """")/=0) then
          err = 2
          exit loop
        endif
      endif
      ij = ij + 1
      s_i = k + 2
      if (eof) exit loop

      else
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
      data(i) = s(s_i:k)
      ij = ij + 1
      s_i = k + 2
      if (ij<length.and.s_i>len(s)) exit loop

      endif
    end do loop

    if (present(num)) num = ij
    if (ij<size(data)) then
      err = -1
    else
      if (present(separator)) then
        if (len(s)-s_i>=0) &
          err = 1
      else
        if (verify(s(s_i:), whitespace)/=0) &
          err = 1
      endif
    endif

    if (present(iostat)) then
      iostat = err
    else
      select case (err)
      case(-1)
        write(0, *) "Error in arraytostring"
        write(0, *) "Too few elements found"
        stop
      case(1)
        write(0, *) "Error in arraytostring"
        write(0, *) "Too many elements found"
        stop
      case(2)
        write(0, *) "Error in arraytostring"
        write(0, *) "Malformed input"
        stop
      end select
    end if

#else
    data = ""
#endif
  end subroutine arraytostring

  subroutine matrixtostring(s, data, separator, csv, num, iostat)
    character(len=*) :: data(:,:)
    character, intent(in), optional :: separator
    logical, intent(in), optional :: csv
    character(len=*), intent(in) :: s
    integer, intent(out), optional :: num
    integer, intent(out), optional :: iostat

#ifndef DUMMYLIB
    logical :: bracketed
    integer :: i, j, ij, k, s_i, err, ios, length
    real :: r, c

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
    data = ""
    ij = 0
    length = size(data)
    loop: do j = 1, size(data, 2)
      do i = 1, size(data, 1)
        if (csv_) then
      if (s_i>len(s)) then
        data(i, j) = ""
        ij = ij + 1
        exit loop
      endif
      k = verify(s(s_i:), achar(10)//achar(13))
      if (k==0) then
        data(i, j) = ""
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
        data(i, j) = s2
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
          if (ij+1==length.and.s(s_i+k-1:s_i+k-1)==",") err = 1
          k = s_i + k - 2
        endif
        data(i, j) = s(s_i:k)
        if (index(data(i, j), """")/=0) then
          err = 2
          exit loop
        endif
      endif
      ij = ij + 1
      s_i = k + 2
      if (eof) exit loop
        else
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
      data(i, j) = s(s_i:k)
      ij = ij + 1
      s_i = k + 2
      if (ij<length.and.s_i>len(s)) exit loop
        endif
      end do
    end do loop

    if (present(num)) num = ij
    if (ij<size(data)) then
      err = -1
    else
      if (present(separator)) then
        if (len(s)-s_i>=0) &
          err = 1
      else
        if (verify(s(s_i:), whitespace)/=0) &
          err = 1
      endif
    endif

    if (present(iostat)) then
      iostat = err
    else
      select case (err)
      case(-1)
        write(0, *) "Error in matrixtostring"
        write(0, *) "Too few elements found"
        stop
      case(1)
        write(0, *) "Error in matrixtostring"
        write(0, *) "Too many elements found"
        stop
      case(2)
        write(0, *) "Error in matrixtostring"
        write(0, *) "Malformed input"
        stop
      end select
    end if

#else
    data = ""
#endif
  end subroutine matrixtostring

  subroutine arraytological(s, data, num, iostat)
    logical, intent(out) :: data(:)
    character(len=*), intent(in) :: s
    integer, intent(out), optional :: num
    integer, intent(out), optional :: iostat

#ifndef DUMMYLIB
    logical :: bracketed
    integer :: i, j, ij, k, s_i, err, ios, length
    real :: r, c


    s_i = 1
    err = 0
    data = .false.
    ij  = 0
    length = size(data)
    loop: do i = 1, size(data)
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
      if (s(s_i:k)=="true".or.s(s_i:k)=="1") then
        data(i) = .true.
      elseif (s(s_i:k)=="false".or.s(s_i:k)=="0") then
        data(i) = .false.
      else
        err = 2
        exit loop
      endif
      ij = ij + 1
      s_i = k + 2
      if (ij<length.and.s_i>len(s)) exit loop

    end do loop

    if (present(num)) num = ij
    if (ij<length) then
      if (err==0) err = -1
    else
      if (verify(s(s_i:), whitespace)/=0) err = 1
    endif


    if (present(iostat)) then
      iostat = err
    else
      select case (err)
      case(-1)
        write(0, *) "Error in arraytological"
        write(0, *) "Too few elements found"
        stop
      case(1)
        write(0, *) "Error in arraytological"
        write(0, *) "Too many elements found"
        stop
      case(2)
        write(0, *) "Error in arraytological"
        write(0, *) "Malformed input"
        stop
      end select
    end if

#else
    data = .false.
#endif
  end subroutine arraytological

  subroutine matrixtological(s, data, num, iostat)
    logical, intent(out) :: data(:,:)
    character(len=*), intent(in) :: s
    integer, intent(out), optional :: num
    integer, intent(out), optional :: iostat

#ifndef DUMMYLIB
    logical :: bracketed
    integer :: i, j, ij, k, s_i, err, ios, length
    real :: r, c


    s_i = 1
    err = 0
    data = .false.
    ij = 0
    length = size(data)
    loop: do j = 1, size(data, 2)
    do i = 1, size(data, 1)
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
      if (s(s_i:k)=="true".or.s(s_i:k)=="1") then
        data(i, j) = .true.
      elseif (s(s_i:k)=="false".or.s(s_i:k)=="0") then
        data(i, j) = .false.
      else
        err = 2
        exit loop
      endif
      ij = ij + 1
      s_i = k + 2
      if (ij<length.and.s_i>len(s)) exit loop

    end do
    end do loop

    if (present(num)) num = ij
    if (ij<length) then
      if (err==0) err = -1
    else
      if (verify(s(s_i:), whitespace)/=0) err = 1
    endif


    if (present(iostat)) then
      iostat = err
    else
      select case (err)
      case(-1)
        write(0, *) "Error in matrixtological"
        write(0, *) "Too few elements found"
        stop
      case(1)
        write(0, *) "Error in matrixtological"
        write(0, *) "Too many elements found"
        stop
      case(2)
        write(0, *) "Error in matrixtological"
        write(0, *) "Malformed input"
        stop
      end select
    end if

#else
    data = .false.
#endif
  end subroutine matrixtological

  subroutine arraytointeger(s, data, num, iostat)
    integer, intent(out) :: data(:)
    character(len=*), intent(in) :: s
    integer, intent(out), optional :: num
    integer, intent(out), optional :: iostat

#ifndef DUMMYLIB
    logical :: bracketed
    integer :: i, j, ij, k, s_i, err, ios, length
    real :: r, c

    s_i = 1
    err = 0
    data = 0
    ij  = 0
    length = size(data)
    loop: do i = 1, size(data)
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
      read(s(s_i:k), *, iostat=ios) data(i)
      if (ios/=0) then
        err = 2
        exit loop
      endif
      ij = ij + 1
      s_i = k + 2
      if (ij<length.and.s_i>len(s)) exit loop

    end do loop

    if (present(num)) num = ij
    if (ij<length) then
      if (err==0) err = -1
    else
      if (verify(s(s_i:), whitespace)/=0) err = 1
    endif


    if (present(iostat)) then
      iostat = err
    else
      select case (err)
      case(-1)
        write(0, *) "Error in arraytointeger"
        write(0, *) "Too few elements found"
        stop
      case(1)
        write(0, *) "Error in arraytointeger"
        write(0, *) "Too many elements found"
        stop
      case(2)
        write(0, *) "Error in arraytointeger"
        write(0, *) "Malformed input"
        stop
      end select
    end if

#else
    data = 0
#endif
  end subroutine arraytointeger

  subroutine matrixtointeger(s, data, num, iostat)
    integer, intent(out) :: data(:, :)
    character(len=*), intent(in) :: s
    integer, intent(out), optional :: num
    integer, intent(out), optional :: iostat

#ifndef DUMMYLIB
    logical :: bracketed
    integer :: i, j, ij, k, s_i, err, ios, length
    real :: r, c


    s_i = 1
    err = 0
    data = 0
    ij = 0
    length = size(data)
    loop: do j = 1, size(data, 2)
    do i = 1, size(data, 1)
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
      read(s(s_i:k), *, iostat=ios) data(i, j)
      if (ios/=0) then
        err = 2
        exit loop
      endif
      ij = ij + 1
      s_i = k + 2
      if (ij<length.and.s_i>len(s)) exit loop

    end do
    end do loop

    if (present(num)) num = ij
    if (ij<length) then
      if (err==0) err = -1
    else
      if (verify(s(s_i:), whitespace)/=0) err = 1
    endif


    if (present(iostat)) then
      iostat = err
    else
      select case (err)
      case(-1)
        write(0, *) "Error in matrixtointeger"
        write(0, *) "Too few elements found"
        stop
      case(1)
        write(0, *) "Error in matrixtointeger"
        write(0, *) "Too many elements found"
        stop
      case(2)
        write(0, *) "Error in matrixtointeger"
        write(0, *) "Malformed input"
        stop
      end select
    end if

#else
    data = 0
#endif
  end subroutine matrixtointeger

  subroutine arraytorealsp(s, data, num, iostat)
    real(sp), intent(out) :: data(:)
    character(len=*), intent(in) :: s
    integer, intent(out), optional :: num
    integer, intent(out), optional :: iostat

#ifndef DUMMYLIB
    logical :: bracketed
    integer :: i, j, ij, k, s_i, err, ios, length
    real :: r, c


    s_i = 1
    err = 0
    data = 0
    ij  = 0
    length = size(data)
    loop: do i = 1, size(data)
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
      read(s(s_i:k), *, iostat=ios) data(i)
      if (ios/=0) then
        err = 2
        exit loop
      endif
      ij = ij + 1
      s_i = k + 2
      if (ij<length.and.s_i>len(s)) exit loop

    end do loop

    if (present(num)) num = ij
    if (ij<length) then
      if (err==0) err = -1
    else
      if (verify(s(s_i:), whitespace)/=0) err = 1
    endif


    if (present(iostat)) then
      iostat = err
    else
      select case (err)
      case(-1)
        write(0, *) "Error in arraytorealsp"
        write(0, *) "Too few elements found"
        stop
      case(1)
        write(0, *) "Error in arraytorealsp"
        write(0, *) "Too many elements found"
        stop
      case(2)
        write(0, *) "Error in arraytorealsp"
        write(0, *) "Malformed input"
        stop
      end select
    end if

#else
    data = 0.0_sp
#endif
  end subroutine arraytorealsp

  subroutine matrixtorealsp(s, data, num, iostat)
    real(sp), intent(out) :: data(:,:)
    character(len=*), intent(in) :: s
    integer, intent(out), optional :: num
    integer, intent(out), optional :: iostat

#ifndef DUMMYLIB
    logical :: bracketed
    integer :: i, j, ij, k, s_i, err, ios, length
    real :: r, c


    s_i = 1
    err = 0
    data = 0
    ij = 0
    length = size(data)
    loop: do j = 1, size(data, 2)
    do i = 1, size(data, 1)
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
      read(s(s_i:k), *, iostat=ios) data(i, j)
      if (ios/=0) then
        err = 2
        exit loop
      endif
      ij = ij + 1
      s_i = k + 2
      if (ij<length.and.s_i>len(s)) exit loop

    end do
    end do loop

    if (present(num)) num = ij
    if (ij<length) then
      if (err==0) err = -1
    else
      if (verify(s(s_i:), whitespace)/=0) err = 1
    endif


    if (present(iostat)) then
      iostat = err
    else
      select case (err)
      case(-1)
        write(0, *) "Error in matrixtorealsp"
        write(0, *) "Too few elements found"
        stop
      case(1)
        write(0, *) "Error in matrixtorealsp"
        write(0, *) "Too many elements found"
        stop
      case(2)
        write(0, *) "Error in matrixtorealsp"
        write(0, *) "Malformed input"
        stop
      end select
    end if

#else
    data = 0.0_sp
#endif
  end subroutine matrixtorealsp

  subroutine arraytorealdp(s, data, num, iostat)
    real(dp), intent(out) :: data(:)
    character(len=*), intent(in) :: s
    integer, intent(out), optional :: num
    integer, intent(out), optional :: iostat

#ifndef DUMMYLIB
    logical :: bracketed
    integer :: i, j, ij, k, s_i, err, ios, length
    real :: r, c


    s_i = 1
    err = 0
    data = 0
    ij  = 0
    length = size(data)
    loop: do i = 1, size(data)
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
      read(s(s_i:k), *, iostat=ios) data(i)
      if (ios/=0) then
        err = 2
        exit loop
      endif
      ij = ij + 1
      s_i = k + 2
      if (ij<length.and.s_i>len(s)) exit loop

    end do loop

    if (present(num)) num = ij
    if (ij<length) then
      if (err==0) err = -1
    else
      if (verify(s(s_i:), whitespace)/=0) err = 1
    endif


    if (present(iostat)) then
      iostat = err
    else
      select case (err)
      case(-1)
        write(0, *) "Error in arraytorealdp"
        write(0, *) "Too few elements found"
        stop
      case(1)
        write(0, *) "Error in arraytorealdp"
        write(0, *) "Too many elements found"
        stop
      case(2)
        write(0, *) "Error in arraytorealdp"
        write(0, *) "Malformed input"
        stop
      end select
    end if

#else
    data = 0.0_dp
#endif
  end subroutine arraytorealdp

  subroutine matrixtorealdp(s, data, num, iostat)
    real(dp), intent(out) :: data(:,:)
    character(len=*), intent(in) :: s
    integer, intent(out), optional :: num
    integer, intent(out), optional :: iostat

#ifndef DUMMYLIB
    logical :: bracketed
    integer :: i, j, ij, k, s_i, err, ios, length
    real :: r, c


    s_i = 1
    err = 0
    data = cmplx(0,0)
    ij = 0
    length = size(data)
    loop: do j = 1, size(data, 2)
    do i = 1, size(data, 1)
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
      read(s(s_i:k), *, iostat=ios) data(i, j)
      if (ios/=0) then
        err = 2
        exit loop
      endif
      ij = ij + 1
      s_i = k + 2
      if (ij<length.and.s_i>len(s)) exit loop

    end do
    end do loop

    if (present(num)) num = ij
    if (ij<length) then
      if (err==0) err = -1
    else
      if (verify(s(s_i:), whitespace)/=0) err = 1
    endif


    if (present(iostat)) then
      iostat = err
    else
      select case (err)
      case(-1)
        write(0, *) "Error in matrixtorealdp"
        write(0, *) "Too few elements found"
        stop
      case(1)
        write(0, *) "Error in matrixtorealdp"
        write(0, *) "Too many elements found"
        stop
      case(2)
        write(0, *) "Error in matrixtorealdp"
        write(0, *) "Malformed input"
        stop
      end select
    end if

#else
    data = 0.0_dp
#endif
  end subroutine matrixtorealdp

  subroutine arraytocomplexsp(s, data, num, iostat)
    complex(sp), intent(out) :: data(:)
    character(len=*), intent(in) :: s
    integer, intent(out), optional :: num
    integer, intent(out), optional :: iostat

#ifndef DUMMYLIB
    logical :: bracketed
    integer :: i, j, ij, k, s_i, err, ios, length
    real :: r, c


    s_i = 1
    err = 0
    data = cmplx(0,0)
    ij  = 0
    length = size(data)
    loop: do i = 1, size(data)
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
      data(i) = cmplx(r, c)
      ij = ij + 1
      s_i = k + 2
      if (ij<length.and.s_i>len(s)) exit loop

    end do loop

    if (present(num)) num = ij
    if (ij<length) then
      if (err==0) err = -1
    else
      if (verify(s(s_i:), whitespace)/=0) err = 1
    endif


    if (present(iostat)) then
      iostat = err
    else
      select case (err)
      case(-1)
        write(0, *) "Error in arraytocomplexsp"
        write(0, *) "Too few elements found"
        stop
      case(1)
        write(0, *) "Error in arraytocomplexsp"
        write(0, *) "Too many elements found"
        stop
      case(2)
        write(0, *) "Error in arraytocomplexsp"
        write(0, *) "Malformed input"
        stop
      end select
    end if

#else
    data = (0.0_sp, 0.0_sp)
#endif
  end subroutine arraytocomplexsp

  subroutine matrixtocomplexsp(s, data, num, iostat)
    complex(sp), intent(out) :: data(:,:)
    character(len=*), intent(in) :: s
    integer, intent(out), optional :: num
    integer, intent(out), optional :: iostat

#ifndef DUMMYLIB
    logical :: bracketed
    integer :: i, j, ij, k, s_i, err, ios, length
    real :: r, c


    s_i = 1
    err = 0
    data = cmplx(0,0)
    ij = 0
    length = size(data)
    loop: do j = 1, size(data, 2)
    do i = 1, size(data, 1)
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
      data(i, j) = cmplx(r, c)
      ij = ij + 1
      s_i = k + 2
      if (ij<length.and.s_i>len(s)) exit loop

    end do
    end do loop

    if (present(num)) num = ij
    if (ij<length) then
      if (err==0) err = -1
    else
      if (verify(s(s_i:), whitespace)/=0) err = 1
    endif


    if (present(iostat)) then
      iostat = err
    else
      select case (err)
      case(-1)
        write(0, *) "Error in matrixtocomplexsp"
        write(0, *) "Too few elements found"
        stop
      case(1)
        write(0, *) "Error in matrixtocomplexsp"
        write(0, *) "Too many elements found"
        stop
      case(2)
        write(0, *) "Error in matrixtocomplexsp"
        write(0, *) "Malformed input"
        stop
      end select
    end if

#else
    data = (0.0_sp, 0.0_sp)
#endif
  end subroutine matrixtocomplexsp

  subroutine arraytocomplexdp(s, data, num, iostat)
    complex(dp), intent(out) :: data(:)
    character(len=*), intent(in) :: s
    integer, intent(out), optional :: num
    integer, intent(out), optional :: iostat

#ifndef DUMMYLIB
    logical :: bracketed
    integer :: i, j, ij, k, s_i, err, ios, length
    real :: r, c


    s_i = 1
    err = 0
    data = cmplx(0)
    ij  = 0
    length = size(data)
    loop: do i = 1, size(data)
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
      data(i) = cmplx(r, c)
      ij = ij + 1
      s_i = k + 2
      if (ij<length.and.s_i>len(s)) exit loop

    end do loop

    if (present(num)) num = ij
    if (ij<length) then
      if (err==0) err = -1
    else
      if (verify(s(s_i:), whitespace)/=0) err = 1
    endif


    if (present(iostat)) then
      iostat = err
    else
      select case (err)
      case(-1)
        write(0, *) "Error in arraytocomplexdp"
        write(0, *) "Too few elements found"
        stop
      case(1)
        write(0, *) "Error in arraytocomplexdp"
        write(0, *) "Too many elements found"
        stop
      case(2)
        write(0, *) "Error in arraytocomplexdp"
        write(0, *) "Malformed input"
        stop
      end select
    end if

#else
    data = (0.0_dp, 0.0_dp)
#endif
  end subroutine arraytocomplexdp

  subroutine matrixtocomplexdp(s, data, num, iostat)
    complex(dp), intent(out) :: data(:,:)
    character(len=*), intent(in) :: s
    integer, intent(out), optional :: num
    integer, intent(out), optional :: iostat

#ifndef DUMMYLIB
    logical :: bracketed
    integer :: i, j, ij, k, s_i, err, ios, length
    real :: r, c


    s_i = 1
    err = 0
    data = cmplx(0,0)
    ij = 0
    length = size(data)
    loop: do j = 1, size(data, 2)
    do i = 1, size(data, 1)
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
      data(i, j) = cmplx(r, c)
      ij = ij + 1
      s_i = k + 2
      if (ij<length.and.s_i>len(s)) exit loop

    end do
    end do loop

    if (present(num)) num = ij
    if (ij<length) then
      if (err==0) err = -1
    else
      if (verify(s(s_i:), whitespace)/=0) err = 1
    endif


    if (present(iostat)) then
      iostat = err
    else
      select case (err)
      case(-1)
        write(0, *) "Error in matrixtocomplexdp"
        write(0, *) "Too few elements found"
        stop
      case(1)
        write(0, *) "Error in matrixtocomplexdp"
        write(0, *) "Too many elements found"
        stop
      case(2)
        write(0, *) "Error in matrixtocomplexdp"
        write(0, *) "Malformed input"
        stop
      end select
    end if

#else
    data = (0.0_dp, 0.0_dp)
#endif
  end subroutine matrixtocomplexdp
  
end module fox_m_fsys_parse_input
