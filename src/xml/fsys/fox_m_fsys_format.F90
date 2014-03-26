module fox_m_fsys_format

!Note that there are several oddities to this package,
!to get round assorted compiler bugs.

!All the _matrix_ subroutines should be straight
!call-throughs to the relevant _array_ subroutine,
!but with flattened arrays. (this would allow easy
!generation of all functions up to 7 dimensions)
!but unfortunately that breaks PGI-6.1, and causes
!errors on Pathscale-2.4.

!The Logical array/matrix functions should be able
!to COUNT their length inline in the specification
!expression, but Pathscale-2.4 gives an error on that.

  use fox_m_fsys_abort_flush, only: pxfflush
  use fox_m_fsys_realtypes, only: sp, dp

  implicit none
  private

#ifndef DUMMYLIB
  integer, parameter :: sig_sp = digits(1.0_sp)/4
  integer, parameter :: sig_dp = digits(1.0_dp)/4 ! Approximate precision worth outputting of each type.

  character(len=*), parameter :: digit = "0123456789:"
  character(len=*), parameter :: hexdigit = "0123456789abcdefABCDEF"
#endif

  interface str
! This is for external use only: str should not be called within this
! file.
! All *_chk subroutines check that the fmt they are passed is valid.
    module procedure str_string, str_string_array, str_string_matrix, &
                     str_integer, str_integer_array, str_integer_matrix, &
                     str_integer_fmt, str_integer_array_fmt, str_integer_matrix_fmt, &
                     str_logical, str_logical_array, str_logical_matrix, &
                     str_real_sp, str_real_sp_fmt_chk, &
                     str_real_sp_array, str_real_sp_array_fmt_chk, &
                     str_real_sp_matrix, str_real_sp_matrix_fmt_chk, &
                     str_real_dp, str_real_dp_fmt_chk, &
                     str_real_dp_array, str_real_dp_array_fmt_chk, &
                     str_real_dp_matrix, str_real_dp_matrix_fmt_chk, &
                     str_complex_sp, str_complex_sp_fmt_chk, &
                     str_complex_sp_array, str_complex_sp_array_fmt_chk, &
                     str_complex_sp_matrix, str_complex_sp_matrix_fmt_chk, &
                     str_complex_dp, str_complex_dp_fmt_chk, &
                     str_complex_dp_array, str_complex_dp_array_fmt_chk, &
                     str_complex_dp_matrix, str_complex_dp_matrix_fmt_chk
  end interface str

#ifndef DUMMYLIB
  interface safestr
! This is for internal use only - no check is made on the validity of 
! any fmt input.
    module procedure str_string, str_string_array, str_string_matrix, &
                     str_integer, str_integer_array, str_integer_matrix, &
                     str_logical, str_logical_array, str_logical_matrix, &
                     str_real_sp, str_real_sp_fmt, &
                     str_real_sp_array, str_real_sp_array_fmt, &
                     str_real_sp_matrix, str_real_sp_matrix_fmt, &
                     str_real_dp, str_real_dp_fmt, &
                     str_real_dp_array, str_real_dp_array_fmt, &
                     str_real_dp_matrix, str_real_dp_matrix_fmt, &
                     str_complex_sp, str_complex_sp_fmt, &
                     str_complex_sp_array, str_complex_sp_array_fmt, &
                     str_complex_sp_matrix, str_complex_sp_matrix_fmt, &
                     str_complex_dp, str_complex_dp_fmt, &
                     str_complex_dp_array, str_complex_dp_array_fmt, &
                     str_complex_dp_matrix, str_complex_dp_matrix_fmt
  end interface safestr

  interface len
    module procedure str_integer_len, str_integer_array_len, str_integer_matrix_len, &
                     str_integer_fmt_len, str_integer_array_fmt_len, str_integer_matrix_fmt_len, &
                     str_logical_len, str_logical_array_len, str_logical_matrix_len, &
                     str_real_sp_len, str_real_sp_fmt_len, &
                     str_real_sp_array_len, str_real_sp_array_fmt_len, &
                     str_real_sp_matrix_len, str_real_sp_matrix_fmt_len, &
                     str_real_dp_len, str_real_dp_fmt_len, &
                     str_real_dp_array_len, str_real_dp_array_fmt_len, &
                     str_real_dp_matrix_len, str_real_dp_matrix_fmt_len, &
                     str_complex_sp_len, str_complex_sp_fmt_len, &
                     str_complex_sp_array_len, str_complex_sp_array_fmt_len, &
                     str_complex_sp_matrix_len, str_complex_sp_matrix_fmt_len, &
                     str_complex_dp_len, str_complex_dp_fmt_len, &
                     str_complex_dp_array_len, str_complex_dp_array_fmt_len, &
                     str_complex_dp_matrix_len, str_complex_dp_matrix_fmt_len
  end interface
#endif

  interface operator(//)
    module procedure concat_str_int, concat_int_str, &
      concat_str_logical, concat_logical_str, &
      concat_real_sp_str, concat_str_real_sp, &
      concat_real_dp_str, concat_str_real_dp, &
      concat_complex_sp_str, concat_str_complex_sp, &
      concat_complex_dp_str, concat_str_complex_dp
  end interface

  public :: str
  public :: operator(//)

#ifndef DUMMYLIB
  public :: str_to_int_10
  public :: str_to_int_16
#endif

contains

#ifndef DUMMYLIB
  ! NB: The len generic module procedure is used in
  !     many initialisation statments (to set the 
  !     length of the output string needed for the
  !     converted number). As of the Fortran 2008
  !     spec every specific function belonging to
  !     a generic used in this way must be defined
  !     in the module before use. This is enforced
  !     by at least version 7.4.4 of the Cray 
  !     Fortran compiler. Hence we put all the *_len
  !     functions here at the top of the file.
  pure function str_string_array_len(st) result(n)
    character(len=*), dimension(:), intent(in) :: st
    integer :: n

    integer :: k

    n = size(st) - 1
    do k = 1, size(st)
      n = n + len(st(k))
    enddo

  end function str_string_array_len

  pure function str_string_matrix_len(st) result(n)
    character(len=*), dimension(:, :), intent(in) :: st
    integer :: n

    n = len(st) * size(st) + size(st) - 1
  end function str_string_matrix_len

  pure function str_integer_len(i) result(n)
    integer, intent(in) :: i
    integer :: n
    
    n = int(log10(real(max(abs(i),1)))) + 1 + dim(-i,0)/max(abs(i),1)

  end function str_integer_len

  pure function str_integer_base_len(i, b) result(n)
    integer, intent(in) :: i, b
    integer :: n
    
    n = int(log10(real(max(abs(i),1)))/log10(real(b))) &
      + 1 + dim(-i,0)/max(abs(i),1)

  end function str_integer_base_len

  pure function str_integer_fmt_len(i, fmt) result(n)
    integer, intent(in) :: i
    character(len=*), intent(in) :: fmt
    integer :: n
    
    select case (len(fmt))
    case(0)
      n = 0
    case(1)
      if (fmt=="x") then
        n = int(log10(real(max(abs(i),1)))/log10(16.0)) + 1 + dim(-i,0)/max(abs(i),1)
      elseif (fmt=="d") then
        n = int(log10(real(max(abs(i),1)))) + 1 + dim(-i,0)/max(abs(i),1)
      else
        return
      endif
    case default
      if (fmt(1:1)/='x'.and.fmt(1:1)/='d') then
        n = 0
      elseif (verify(fmt(2:), digit)==0) then
        n = str_to_int_10(fmt(2:))
      else
        n = 0 
      endif
    end select

  end function str_integer_fmt_len

  pure function str_integer_array_len(ia) result(n)
    integer, dimension(:), intent(in) :: ia
    integer :: n
    
    integer :: j

    n = size(ia) - 1

    do j = 1, size(ia)
      n = n + len(ia(j))
    enddo

  end function str_integer_array_len

  pure function str_integer_array_fmt_len(ia, fmt) result(n)
    integer, dimension(:), intent(in) :: ia
    character(len=*), intent(in) :: fmt
    integer :: n
    
    integer :: j

    n = size(ia) - 1

    do j = 1, size(ia)
      n = n + len(ia(j), fmt)
    enddo

  end function str_integer_array_fmt_len

  pure function str_integer_matrix_len(ia) result(n)
    integer, dimension(:,:), intent(in) :: ia
    integer :: n

    integer :: j, k

    n = size(ia) - 1

    do k = 1, size(ia, 2)
      do j = 1, size(ia, 1)
        n = n + len(ia(j, k))
      enddo
    enddo

  end function str_integer_matrix_len

  pure function str_integer_matrix_fmt_len(ia, fmt) result(n)
    integer, dimension(:,:), intent(in) :: ia
    character(len=*), intent(in) :: fmt
    integer :: n

    integer :: j, k

    n = size(ia) - 1

    do k = 1, size(ia, 2)
      do j = 1, size(ia, 1)
        n = n + len(ia(j, k), fmt)
      enddo
    enddo

  end function str_integer_matrix_fmt_len

  pure function str_logical_len(l) result (n)
    logical, intent(in) :: l
    integer :: n
    
    if (l) then
      n = 4
    else
      n = 5
    endif
  end function str_logical_len

  pure function str_logical_array_len(la) result(n)
! This function should be inlined in the declarations of
! str_logical_array below but PGI and pathscale don't like it.
    logical, dimension(:), intent(in)   :: la
    integer :: n
    n = 5*size(la) - 1 + count(.not.la)
  end function str_logical_array_len

  pure function str_logical_matrix_len(la) result(n)
! This function should be inlined in the declarations of
! str_logical_matrix below but PGI and pathscale don't like it.
    logical, dimension(:,:), intent(in)   :: la
    integer :: n
    n = 5*size(la) - 1 + count(.not.la)
  end function str_logical_matrix_len

  pure function str_real_sp_fmt_len(x, fmt) result(n)
    real(sp), intent(in) :: x
    character(len=*), intent(in) :: fmt
    integer :: n

    integer :: dec, sig
    integer :: e

    if (.not.checkFmt(fmt)) then
      n = 0
      return
    endif

    if (x == 0.0_sp) then
      e = 1
    else
      e = floor(log10(abs(x)))
    endif
      
    if (x < 0.0_sp) then
      n = 1
    else
      n = 0
    endif
      
    if (len(fmt) == 0) then
      sig = sig_sp

      n = n + sig + 2 + len(e) 
      ! for the decimal point and the e

    elseif (fmt(1:1) == "s") then
      if (len(fmt) > 1) then
        sig = str_to_int_10(fmt(2:))
      else
        sig = sig_sp
      endif
      sig = max(sig, 1)
      sig = min(sig, digits(1.0_sp))

      if (sig > 1) n = n + 1 
      ! for the decimal point
      
      n = n + sig + 1 + len(e)

    elseif (fmt(1:1) == "r") then

      if (len(fmt) > 1) then
        dec = str_to_int_10(fmt(2:))
      else
        dec = sig_sp - e - 1
      endif
      dec = min(dec, digits(1.0_sp)-e)
      dec = max(dec, 0)

      if (dec > 0) n = n + 1
      if (abs(x) >= 1.0_sp) n = n + 1

      ! Need to know if there's an overflow ....
      if (e+dec+1 > 0) then
        if (index(real_sp_str(abs(x), e+dec+1), "!") == 1) &
             e = e + 1
      endif

      n = n + abs(e) + dec

    endif

  end function str_real_sp_fmt_len

  pure function str_real_sp_len(x) result(n)
    real(sp), intent(in) :: x
    integer :: n

    n = len(x, "")

  end function str_real_sp_len

  pure function str_real_sp_array_len(xa) result(n)
    real(sp), dimension(:), intent(in) :: xa
    integer :: n

    integer :: k

    n = size(xa) - 1
    do k = 1, size(xa)
      n = n + len(xa(k), "")
    enddo

  end function str_real_sp_array_len

  pure function str_real_sp_array_fmt_len(xa, fmt) result(n)
    real(sp), dimension(:), intent(in) :: xa
    character(len=*), intent(in) :: fmt
    integer :: n

    integer :: k

    n = size(xa) - 1
    do k = 1, size(xa)
      n = n + len(xa(k), fmt)
    enddo
    
  end function str_real_sp_array_fmt_len

  pure function str_real_sp_matrix_fmt_len(xa, fmt) result(n)
    real(sp), dimension(:,:), intent(in) :: xa
    character(len=*), intent(in) :: fmt
    integer :: n

    integer :: j, k

    n = size(xa) - 1
    do k = 1, size(xa, 2)
      do j = 1, size(xa, 1)
        n = n + len(xa(j,k), fmt)
      enddo
    enddo

  end function str_real_sp_matrix_fmt_len

  pure function str_real_sp_matrix_len(xa) result(n)
    real(sp), dimension(:,:), intent(in) :: xa
    integer :: n

    n = len(xa, "")
  end function str_real_sp_matrix_len

  pure function str_real_dp_fmt_len(x, fmt) result(n)
    real(dp), intent(in) :: x
    character(len=*), intent(in) :: fmt
    integer :: n

    integer :: dec, sig
    integer :: e

    if (.not.checkFmt(fmt)) then
      n = 0
      return
    endif

    if (x == 0.0_dp) then
      e = 1
    else
      e = floor(log10(abs(x)))
    endif
      
    if (x < 0.0_dp) then
      n = 1
    else
      n = 0
    endif
      
    if (len(fmt) == 0) then
      sig = sig_dp

      n = n + sig + 2 + len(e) 
      ! for the decimal point and the e

    elseif (fmt(1:1) == "s") then
      if (len(fmt) > 1) then
        sig = str_to_int_10(fmt(2:))
      else
        sig = sig_dp
      endif
      sig = max(sig, 1)
      sig = min(sig, digits(1.0_dp))

      if (sig > 1) n = n + 1 
      ! for the decimal point
      
      n = n + sig + 1 + len(e)

    elseif (fmt(1:1) == "r") then

      if (len(fmt) > 1) then
        dec = str_to_int_10(fmt(2:))
      else
        dec = sig_dp - e - 1
      endif
      dec = min(dec, digits(1.0_dp)-e)
      dec = max(dec, 0)

      if (dec > 0) n = n + 1
      if (abs(x) >= 1.0_dp) n = n + 1

      ! Need to know if there's an overflow ....
      if (e+dec+1 > 0) then
        if (index(real_dp_str(abs(x), e+dec+1), "!") == 1) &
             e = e + 1
      endif

      n = n + abs(e) + dec

    endif

  end function str_real_dp_fmt_len

  pure function str_real_dp_len(x) result(n)
    real(dp), intent(in) :: x
    integer :: n

    n = len(x, "")

  end function str_real_dp_len

  pure function str_real_dp_array_len(xa) result(n)
    real(dp), dimension(:), intent(in) :: xa
    integer :: n

    integer :: k

    n = size(xa) - 1
    do k = 1, size(xa)
      n = n + len(xa(k), "")
    enddo
    
  end function str_real_dp_array_len

  pure function str_real_dp_array_fmt_len(xa, fmt) result(n)
    real(dp), dimension(:), intent(in) :: xa
    character(len=*), intent(in) :: fmt
    integer :: n

    integer :: k

    n = size(xa) - 1
    do k = 1, size(xa)
      n = n + len(xa(k), fmt)
    enddo
    
  end function str_real_dp_array_fmt_len

  pure function str_real_dp_matrix_fmt_len(xa, fmt) result(n)
    real(dp), dimension(:,:), intent(in) :: xa
    character(len=*), intent(in) :: fmt
    integer :: n

    integer :: j, k

    n = size(xa) - 1
    do k = 1, size(xa, 2)
      do j = 1, size(xa, 1)
        n = n + len(xa(j,k), fmt)
      enddo
    enddo

  end function str_real_dp_matrix_fmt_len

  pure function str_real_dp_matrix_len(xa) result(n)
    real(dp), dimension(:,:), intent(in) :: xa
    integer :: n

    n = len(xa, "")
  end function str_real_dp_matrix_len

  pure function str_complex_sp_fmt_len(c, fmt) result(n)
    complex(sp), intent(in) :: c
    character(len=*), intent(in) :: fmt
    integer :: n

    real(sp) :: re, im
    re = real(c)
    im = aimag(c)

    n = len(re, fmt) + len(im, fmt) + 6
  end function str_complex_sp_fmt_len

  pure function str_complex_sp_len(c) result(n)
    complex(sp), intent(in) :: c
    integer :: n

    n = len(c, "")
  end function str_complex_sp_len

  pure function str_complex_sp_array_fmt_len(ca, fmt) result(n)
    complex(sp), dimension(:), intent(in) :: ca
    character(len=*), intent(in) :: fmt
    integer :: n

    integer :: i

    n = size(ca) - 1
    do i = 1, size(ca)
      n = n + len(ca(i), fmt) 
    enddo
  end function str_complex_sp_array_fmt_len

  pure function str_complex_sp_array_len(ca) result(n)
    complex(sp), dimension(:), intent(in) :: ca
    integer :: n

    n = len(ca, "")
  end function str_complex_sp_array_len

  pure function str_complex_sp_matrix_fmt_len(ca, fmt) result(n)
    complex(sp), dimension(:, :), intent(in) :: ca
    character(len=*), intent(in) :: fmt
    integer :: n

    integer :: i, j

    n = size(ca) - 1
    do i = 1, size(ca, 1)
      do j = 1, size(ca, 2)
        n = n + len(ca(i, j), fmt)
      enddo
    enddo
  end function str_complex_sp_matrix_fmt_len

  pure function str_complex_sp_matrix_len(ca) result(n)
    complex(sp), dimension(:, :), intent(in) :: ca
    integer :: n

    n = len(ca, "")
  end function str_complex_sp_matrix_len
  
  pure function str_complex_dp_fmt_len(c, fmt) result(n)
    complex(dp), intent(in) :: c
    character(len=*), intent(in) :: fmt
    integer :: n

    real(dp) :: re, im
    re = real(c)
    im = aimag(c)

    n = len(re, fmt) + len(im, fmt) + 6
  end function str_complex_dp_fmt_len

  pure function str_complex_dp_len(c) result(n)
    complex(dp), intent(in) :: c
    integer :: n

    n = len(c, "")
  end function str_complex_dp_len

  pure function str_complex_dp_array_fmt_len(ca, fmt) result(n)
    complex(dp), dimension(:), intent(in) :: ca
    character(len=*), intent(in) :: fmt
    integer :: n

    integer :: i

    n = size(ca) - 1
    do i = 1, size(ca)
      n = n + len(ca(i), fmt)
    enddo
  end function str_complex_dp_array_fmt_len

  pure function str_complex_dp_array_len(ca) result(n)
    complex(dp), dimension(:), intent(in) :: ca
    integer :: n

    n = len(ca, "")
  end function str_complex_dp_array_len

  pure function str_complex_dp_matrix_fmt_len(ca, fmt) result(n)
    complex(dp), dimension(:, :), intent(in) :: ca
    character(len=*), intent(in) :: fmt
    integer :: n

    integer :: i, j

    n = size(ca) - 1
    do i = 1, size(ca, 1)
      do j = 1, size(ca, 2)
        n = n + len(ca(i, j), fmt)
      enddo
    enddo
  end function str_complex_dp_matrix_fmt_len
     
  pure function str_complex_dp_matrix_len(ca) result(n)
    complex(dp), dimension(:, :), intent(in) :: ca
    integer :: n

    n = len(ca, "")
  end function str_complex_dp_matrix_len
#endif

#ifndef DUMMYLIB
  subroutine FoX_error(msg)
    ! Emit error message and stop.
    ! No clean up is done here, but this can
    ! be overridden to include clean-up routines
    character(len=*), intent(in) :: msg

    write(0,'(a)') 'ERROR(FoX)'
    write(0,'(a)')  msg
    call pxfflush(0)

    stop

  end subroutine FoX_error


  pure function str_to_int_10(str) result(n)
    ! Takes a string containing digits, and returns
    ! the integer representable by those digits.
    ! Does not deal with negative numbers, and
    ! presumes that the number is representable
    ! in a default integer
    ! Error is flagged by returning -1
    character(len=*), intent(in) :: str
    integer :: n

    integer :: max_power, i, j

    if (verify(str, digit) > 0) then
      n = -1
      return
    endif

    max_power = len(str) - 1

    n = 0
    do i = 0, max_power
      j = max_power - i + 1
      n = n + (index(digit, str(j:j)) - 1) * 10**i
    enddo

  end function str_to_int_10

  pure function str_to_int_16(str) result(n)
    ! Takes a string containing hexadecimal digits, and returns
    ! the integer representable by those digits.
    ! Does not deal with negative numbers, and
    ! presumes that the number is representable
    ! in a default integer
    ! Error is flagged by returning -1
    character(len=*), intent(in) :: str
    integer :: n
    
    character(len=len(str)) :: str_l
    integer :: max_power, i, j

    if (verify(str, hexdigit) == 0) then
       str_l = to_lower(str)
    else
      n = -1
      return
    endif

    max_power = len(str) - 1

    n = 0
    do i = 0, max_power
      j = max_power - i + 1
      n = n + (index(hexdigit, str_l(j:j)) - 1) * 16**i
    enddo

  contains
    pure function to_lower(s) result(s2)
      character(len=*), intent(in) :: s
      character(len=len(s)) :: s2
      character(len=*), parameter :: hex = "abcdef"
      integer :: j, k
      do j = 1, len(s)
        k = index('ABCDEF', s(j:j))
        if (k > 0) then
          s2(j:j) = hex(k:k)
        else
          s2(j:j) = s(j:j)
        endif
      enddo
    end function to_lower
         
  end function str_to_int_16
#endif

  pure function str_string(st) result(s)
    character(len=*), intent(in) :: st
#ifdef DUMMYLIB
    character(len=1) :: s
    s = " " 
#else
    character(len=len(st)) :: s
    s = st
#endif
  end function str_string

  pure function str_string_array(st, delimiter) result(s)
    character(len=*), dimension(:), intent(in) :: st
    character(len=1), intent(in), optional :: delimiter
#ifdef DUMMYLIB
    character(len=1) :: s
    s = " "
#else
    character(len=str_string_array_len(st)) :: s
    
    integer :: k, n
    character(len=1) :: d
    
    if (present(delimiter)) then
      d = delimiter
    else
      d = " "
    endif

    n = 1
    do k = 1, size(st) - 1
      s(n:n+len(st(k))) = st(k)//d
      n = n + len(st(k)) + 1
    enddo
    s(n:) = st(k)
#endif
  end function str_string_array

  pure function str_string_matrix(st, delimiter) result(s)
    character(len=*), dimension(:, :), intent(in) :: st
    character(len=1), intent(in), optional :: delimiter
#ifdef DUMMYLIB
    character(len=1) :: s
    s = " "
#else
    character(len=str_string_matrix_len(st)) :: s
    
    integer :: j, k, n
    character(len=1) :: d

    if (present(delimiter)) then
      d = delimiter
    else
      d = " "
    endif

    s(1:len(st)) = st(1,1)
    n = len(st) + 1
    do j = 2, size(st, 1)
      s(n:n+len(st)) = d//st(j,1)
        n = n + len(st) + 1
    enddo
    do k = 2, size(st, 2)
      do j = 1, size(st, 1)
        s(n:n+len(st(j,k))) = d//st(j,k)
        n = n + len(st) + 1
      enddo
    enddo
#endif
  end function str_string_matrix

  pure function str_integer(i) result(s)
    integer, intent(in) :: i
#ifdef DUMMYLIB
    character(len=1) :: s
    s = " "
#else
    character(len=str_integer_len(i)) :: s

    integer :: b, ii, j, k, n

    b = 10

    if (i < 0) then
      s(1:1) = "-"
      n = 2
    else
      n = 1
    endif
    ii = abs(i)
    do k = len(s) - n, 0, -1
      j = ii/(b**k)
      ii = ii - j*(b**k)
      s(n:n) = digit(j+1:j+1)
      n = n + 1
    enddo
#endif
  end function str_integer

  pure function str_integer_fmt(i, fmt) result(s)
    integer, intent(in) :: i
    character(len=*), intent(in):: fmt
#ifdef DUMMYLIB
    character(len=1) :: s
    s = " "
#else
    character(len=str_integer_fmt_len(i, fmt)) :: s

    character :: f
    integer :: b, ii, j, k, n, ls
 
    if (len(fmt)>0) then
      if (fmt(1:1)=="d") then
        f = 'd'
        b = 10
      elseif (fmt(1:1)=="x") then
        f = 'x'
        b = 16
      else
        ! Undefined outcome
        s = ""
        return
      endif
    else
      ! Undefined outcome
      s = ""
      return
    endif

    ls = str_integer_base_len(i, b)
    n = len(s) - ls + 1

    if (i < 0) then
      if (n>0) s(:n) = "-"//repeat("0", n-1)
      n = n + 1
    else
      if (n>1) s(:n) = repeat("0", n)
    endif

    ii = abs(i)
    do k = 1, -n + 1
      j = ii/(b**k)
      ii = ii - j*(b**k)
      n = n + 1
    enddo
    do k = len(s) - n, 0, -1
      j = ii/(b**k)
      ii = ii - j*(b**k)
      s(n:n) = hexdigit(j+1:j+1)
      n = n + 1
    enddo
#endif
  end function str_integer_fmt

  pure function str_integer_array(ia) result(s)
    integer, dimension(:), intent(in) :: ia
#ifdef DUMMYLIB
    character(len=1) :: s
#else
    character(len=len(ia, "d")) :: s

    integer :: j, k, n

    n = 1
    do k = 1, size(ia) - 1
      j = len(ia(k))
      s(n:n+j) = str(ia(k))//" "
      n = n + j + 1
    enddo
    s(n:) = str(ia(k))
#endif
  end function str_integer_array


  function str_integer_array_fmt(ia, fmt) result(s)
    integer, dimension(:), intent(in) :: ia
    character(len=*), intent(in) :: fmt
#ifdef DUMMYLIB
    character(len=1) :: s
    s = " "
#else
    character(len=len(ia, fmt)) :: s

    integer :: j, k, n

    n = 1
    do k = 1, size(ia) - 1
      j = len(ia(k), fmt)
      s(n:n+j) = str(ia(k), fmt)//" "
      n = n + j + 1
    enddo
    s(n:) = str(ia(k), fmt)
#endif
  end function str_integer_array_fmt

  pure function str_integer_matrix(ia) result(s)
    integer, dimension(:,:), intent(in) :: ia
#ifdef DUMMYLIB
    character(len=1) :: s
    s = " "
#else
    character(len=len(ia, "d")) :: s

    integer :: j, k, n

    s(:len(ia(1,1))) = str(ia(1,1))
    n = len(ia(1,1)) + 1
    do j = 2, size(ia, 1)
      s(n:n+len(ia(j,1))) = " "//str(ia(j,1))
      n = n + len(ia(j,1)) + 1
    enddo
    do k = 2, size(ia, 2) 
      do j = 1, size(ia, 1)
        s(n:n+len(ia(j,k))) = " "//str(ia(j,k))
        n = n + len(ia(j,k)) + 1
      enddo
    enddo
#endif
  end function str_integer_matrix


  pure function str_integer_matrix_fmt(ia, fmt) result(s)
    integer, dimension(:,:), intent(in) :: ia
    character(len=*), intent(in) :: fmt
#ifdef DUMMYLIB
    character(len=1) :: s
    s = " "
#else
    character(len=len(ia, fmt)) :: s

    integer :: j, k, n

    s(:len(ia(1,1), fmt)) = str(ia(1,1), fmt)
    n = len(ia(1,1), fmt) + 1
    do j = 2, size(ia, 1)
      s(n:n+len(ia(j,1), fmt)) = " "//str(ia(j,1), fmt)
      n = n + len(ia(j,1), fmt) + 1
    enddo
    do k = 2, size(ia, 2) 
      do j = 1, size(ia, 1)
        s(n:n+len(ia(j,k), fmt)) = " "//str(ia(j,k), fmt)
        n = n + len(ia(j,k), fmt) + 1
      enddo
    enddo
#endif
  end function str_integer_matrix_fmt

  pure function str_logical(l) result(s)
    logical, intent(in) :: l
#ifdef DUMMYLIB
    character(len=1) :: s
    s = " "
#else
! Pathscale 2.5 gets it wrong if we use merge here
!    character(len=merge(4,5,l)) :: s
! And g95 (sep2007) cant resolve the generic here
    character(len=str_logical_len(l)) :: s
    
    if (l) then
      s="true"
    else
      s="false"
    endif
#endif
  end function str_logical

  pure function str_logical_array(la) result(s)
    logical, dimension(:), intent(in)   :: la
#ifdef DUMMYLIB
    character(len=1) :: s
    s = " "
#else
    character(len=len(la)) :: s
    
    integer :: k, n

    n = 1
    do k = 1, size(la) - 1
      if (la(k)) then
        s(n:n+3) = "true"
        n = n + 5
      else
        s(n:n+4) = "false"
        n = n + 6
      endif
      s(n-1:n-1) = " "
    enddo
    if (la(k)) then
      s(n:) = "true"
    else
      s(n:) = "false"
    endif
#endif
  end function str_logical_array

  pure function str_logical_matrix(la) result(s)
    logical, dimension(:,:), intent(in)   :: la
#ifdef DUMMYLIB
    character(len=1) :: s
    s = " "
#else
    character(len=len(la)) :: s

    integer :: j, k, n

    if (la(1,1)) then
       s(:4) = "true"
       n = 5
    else
       s(:5) = "false"
       n = 6
    endif
    do j = 2, size(la, 1)
      s(n:n) = " "
      if (la(j,1)) then
        s(n+1:n+4) = "true"
        n = n + 5
      else
        s(n+1:n+5) = "false"
        n = n + 6
      endif
    enddo
    do k = 2, size(la, 2)
      do j = 1, size(la, 1)
        s(n:n) = " "
        if (la(j,k)) then
          s(n+1:n+4) = "true"
          n = n + 5
        else
          s(n+1:n+5) = "false"
          n = n + 6
        endif
      enddo
    enddo
#endif
  end function str_logical_matrix
  
#ifndef DUMMYLIB
  ! In order to convert real numbers to strings, we need to
  ! perform an internal write - but how long will the 
  ! resultant string be? We don't know & there is no way
  ! to discover for an arbitrary format. Therefore, 
  ! (if we have the capability; f95 or better)
  ! we assume it will be less than 100 characters, write
  ! it to a string of that length, then remove leading &
  ! trailing whitespace. (this means that if the specified
  ! format includes whitespace, this will be lost.)
  !
  ! If we are working with an F90-only compiler, then
  ! we cannot do this trick - the output string will
  ! always be 100 chars in length, though we will remove
  ! leading whitespace. 


  ! The standard Fortran format functions do not give us
  ! enough control, so we write our own real number formatting
  ! routines here. For each real type, we optionally take a
  ! format like so:
  ! "r<integer>" which will produce output without an exponent,
  ! and <integer> digits after the decimal point.
  ! or
  ! "s<integer>": which implies scientific notation, with an 
  ! exponent, with <integer> significant figures.
  ! If the integer is absent, then the precision will be
  ! half of the number of significant figures available
  ! for that real type.
  ! The absence of a format implies scientific notation, with
  ! the default precision.

  ! These routines are fairly imperfect - they are inaccurate for
  ! the lower-end bits of the number, since they work by simple
  ! multiplications by 10.
  ! Also they will probably be orders of magnitude slower than library IO.
  ! Ideally they'd be rewritten to convert from teh native format by
  ! bit-twidding. Not sure how to do that portably though.

  ! The format specification could be done more nicely - but unfortunately
  ! not in F95 due to *stupid* restrictions on specification expressions.

  ! And I wouldn't have to invent my own format specification if Fortran
  ! had a proper IO library anyway.

!FIXME Signed zero is not handled correctly; don't quite understand why.
!FIXME too much duplication between sp & dp, we should m4.

  pure function real_sp_str(x, sig) result(s)
    real(sp), intent(in) :: x
    integer, intent(in) :: sig
    character(len=sig) :: s
    ! make a string of numbers sig long of x.
    integer :: e, i, j, k, n
    real(sp) :: x_

    if (sig < 1) then
      s ="" 
      return
    endif

    if (x == 0.0_sp) then
      e = 1
    else
      e = floor(log10(abs(x)))
    endif
    x_ = abs(x)
    ! Have to do this next in a loop rather than just exponentiating in
    ! order to  avoid under/over-flow.
    do i = 1, abs(e)
      ! Have to multiply by 10^-e rather than divide by 10^e
      ! to avoid rounding errors.
      x_ = x_ * (10.0_sp**(-abs(e)/e))
    enddo
    n = 1
    do k = sig - 2, 0, -1
      ! This baroque way of taking int() ensures the optimizer 
      ! stores it in j without keeping a different value in cache.
      j = iachar(digit(int(x_)+1:int(x_)+1)) - 48
      if (j==10) then
        ! This can happen if, on the previous cycle, int(x_) in 
        ! the line above gave a result approx. 1.0 less than
        ! expected.
        ! In this case we want to quit the cycle & just get 999... to the end
        s(n:) = repeat("9", sig - n + 1)
        return
      endif
      s(n:n) = digit(j+1:j+1)
      n = n + 1
      x_ = (x_ - j) * 10.0_sp
    enddo
    j = nint(x_)
    if (j == 10) then
      ! Now round ...
      s(n:n) = "9"
      ! Are they all 9's?
      i = verify(s, "9", .true.)
      if (i == 0) then
        s(1:1) = "!"
        ! overflow
        return
      endif
      j = index(digit, s(i:i))
      s(i:i) = digit(j+1:j+1)
      s(i+1:) = repeat("0", sig - i + 1)
    else
      s(n:n) = digit(j+1:j+1)
    endif

  end function real_sp_str

#endif

  function str_real_sp_fmt_chk(x, fmt) result(s)
    real(sp), intent(in) :: x
    character(len=*), intent(in) :: fmt
#ifdef DUMMYLIB
    character(len=1) :: s
    s = " "
#else
    character(len=len(x, fmt)) :: s

    if (checkFmt(fmt)) then
      s = safestr(x, fmt)
    else
      call FoX_error("Invalid format: "//fmt)
    endif
#endif
  end function str_real_sp_fmt_chk

#ifndef DUMMYLIB
  pure function str_real_sp_fmt(x, fmt) result(s)
    real(sp), intent(in) :: x
    character(len=*), intent(in) :: fmt
    character(len=len(x, fmt)) :: s

    integer :: sig, dec
    integer :: e, n
    character(len=len(x, fmt)) :: num !this will always be enough memory.

    if (x == 0.0_sp) then
      e = 0
    else
      e = floor(log10(abs(x)))
    endif

    if (x < 0.0_sp) then
      s(1:1) = "-"
      n = 2
    else
      n = 1
    endif

    if (len(fmt) == 0) then

      sig = sig_sp

      num = real_sp_str(abs(x), sig)
      if (num(1:1) == "!") then
        e = e + 1
        num = "1"//repeat("0",len(num)-1)
      endif

      if (sig == 1) then
        s(n:n) = num
        n = n + 1
      else
        s(n:n+1) = num(1:1)//"."
        s(n+2:n+sig) = num(2:)
        n = n + sig + 1
      endif

      s(n:n) = "e"
      s(n+1:) = str(e)

    elseif (fmt(1:1) == "s") then

      if (len(fmt) > 1) then
        sig = str_to_int_10(fmt(2:))
      else
        sig = sig_sp
      endif
      sig = max(sig, 1)
      sig = min(sig, digits(1.0_sp))

      num = real_sp_str(abs(x), sig)
      if (num(1:1) == "!") then
        e = e + 1
        num = "1"//repeat("0",len(num)-1)
      endif

      if (sig == 1) then
        s(n:n) = num
        n = n + 1
      else
        s(n:n+1) = num(1:1)//"."
        s(n+2:n+sig) = num(2:)
        n = n + sig + 1
      endif

      s(n:n) = "e"
      s(n+1:) = str(e)

    elseif (fmt(1:1) == "r") then

      if (len(fmt) > 1) then
        dec = str_to_int_10(fmt(2:))
      else
        dec = sig_sp - e - 1
      endif
      dec = min(dec, digits(1.0_sp)-e-1)
      dec = max(dec, 0)

      if (e+dec+1 > 0) then
        num = real_sp_str(abs(x), e+dec+1)
      else
        num = ""
      endif
      if (num(1:1) == "!") then
        e = e + 1
        num = "1"//repeat("0",len(num)-1)
      endif

      if (abs(x) >= 1.0_sp) then
        s(n:n+e) = num(:e+1)
        n = n + e + 1
        if (dec > 0) then
          s(n:n) = "."
          n = n + 1
          s(n:) = num(e+2:)
        endif
      else
        s(n:n) = "0"
        if (dec > 0) then
          s(n+1:n+1) = "."
          n = n + 2
          if (dec < -e-1) then
            s(n:) = repeat("0", dec)
          else
            s(n:n-e-2) = repeat("0", max(-e-1,0))
            n = n - min(e,-1) - 1
            if (n <= len(s)) then
              s(n:) = num
            endif
          endif
        endif
      endif

    endif

  end function str_real_sp_fmt
#endif

  pure function str_real_sp(x) result(s)
    real(sp), intent(in) :: x
#ifdef DUMMYLIB
    character(len=1) :: s
    s = " "
#else
    character(len=len(x)) :: s

    s = safestr(x, "")
#endif
  end function str_real_sp

  pure function str_real_sp_array(xa) result(s)
    real(sp), dimension(:), intent(in) :: xa
#ifdef DUMMYLIB
    character(len=1) :: s
    s = " "
#else
    character(len=len(xa)) :: s
    
    integer :: j, k, n

    n = 1
    do k = 1, size(xa) - 1
      j = len(xa(k), "")
      s(n:n+j) = safestr(xa(k), "")//" "
      n = n + j + 1
    enddo
    s(n:) = safestr(xa(k), "")
#endif
  end function str_real_sp_array

#ifndef DUMMYLIB
  pure function str_real_sp_array_fmt(xa, fmt) result(s)
    real(sp), dimension(:), intent(in) :: xa
    character(len=*), intent(in) :: fmt
    character(len=len(xa, fmt)) :: s
    
    integer :: j, k, n

    n = 1
    do k = 1, size(xa) - 1
      j = len(xa(k), fmt)
      s(n:n+j) = safestr(xa(k), fmt)//" "
      n = n + j + 1
    enddo
    s(n:) = safestr(xa(k), fmt)

  end function str_real_sp_array_fmt
#endif

  function str_real_sp_array_fmt_chk(xa, fmt) result(s)
    real(sp), dimension(:), intent(in) :: xa
    character(len=*), intent(in) :: fmt
#ifdef DUMMYLIB
    character(len=1) :: s
    s = " "
#else
    character(len=len(xa, fmt)) :: s
    
    if (checkFmt(fmt)) then
      s = safestr(xa, fmt)
    else
      call FoX_error("Invalid format: "//fmt)
    endif
#endif
  end function str_real_sp_array_fmt_chk

#ifndef DUMMYLIB
  pure function str_real_sp_matrix_fmt(xa, fmt) result(s)
    real(sp), dimension(:,:), intent(in) :: xa
    character(len=*), intent(in) :: fmt
    character(len=len(xa,fmt)) :: s

    integer :: i, j, k, n

    i = len(xa(1,1), fmt)
    s(:i) = safestr(xa(1,1), fmt)
    n = i + 1
    do j = 2, size(xa, 1)
      i = len(xa(j,1), fmt)
      s(n:n+i) = " "//safestr(xa(j,1), fmt)
      n = n + i + 1
    enddo
    do k = 2, size(xa, 2)
      do j = 1, size(xa, 1)
        i = len(xa(j,k), fmt)
        s(n:n+i) = " "//safestr(xa(j,k), fmt)
        n = n + i + 1
      enddo
    enddo

  end function str_real_sp_matrix_fmt
#endif

  function str_real_sp_matrix_fmt_chk(xa, fmt) result(s)
    real(sp), dimension(:,:), intent(in) :: xa
    character(len=*), intent(in) :: fmt
#ifdef DUMMYLIB
    character(len=1) :: s
    s = " "
#else
    character(len=len(xa,fmt)) :: s

    if (checkFmt(fmt)) then
      s = safestr(xa, fmt)
    else
      call FoX_error("Invalid format: "//fmt)
    end if
#endif
  end function str_real_sp_matrix_fmt_chk

  pure function str_real_sp_matrix(xa) result(s)
    real(sp), dimension(:,:), intent(in) :: xa
#ifdef DUMMYLIB
    character(len=1) :: s
    s = " "
#else
    character(len=len(xa)) :: s

    s = safestr(xa, "")
#endif
  end function str_real_sp_matrix
    
#ifndef DUMMYLIB
  pure function real_dp_str(x, sig) result(s)
    real(dp), intent(in) :: x
    integer, intent(in) :: sig
    character(len=sig) :: s
    ! make a string of numbers sig long of x.
    integer :: e, i, j, k, n
    real(dp) :: x_

    if (sig < 1) then
      s ="" 
      return
    endif

    if (x == 0.0_dp) then
      e = 1
    else
      e = floor(log10(abs(x)))
    endif
    x_ = abs(x)
    ! Have to do this next in a loop rather than just exponentiating in
    ! order to  avoid under/over-flow.
    do i = 1, abs(e)
      ! Have to multiply by 10^-e rather than divide by 10^e
      ! to avoid rounding errors.
      x_ = x_ * (10.0_dp**(-abs(e)/e))
    enddo
    n = 1
    do k = sig - 2, 0, -1
      ! This baroque way of taking int() ensures the optimizer definitely
      ! stores it in j without keeping a different value in cache.
      j = iachar(digit(int(x_)+1:int(x_)+1)) - 48
      if (j==10) then
        ! This can happen if, on the previous cycle, int(x_) in 
        ! the line above gave a result almost exactly 1.0 less than
        ! expected - but FP arithmetic is not consistent.
        ! In this case we want to quit the cycle & just get 999... to the end
        s(n:) = repeat("9", sig - n + 1)
        return
      endif
      s(n:n) = digit(j+1:j+1)
      n = n + 1
      x_ = (x_ - j) * 10.0_dp
    enddo
    j = nint(x_)
    if (j == 10) then
      ! Now round ...
      s(n:n) = "9"
      i = verify(s, "9", .true.)
      if (i == 0) then
        s(1:1) = "!"
        !overflow
        return
      endif
      j = index(digit, s(i:i))
      s(i:i) = digit(j+1:j+1)
      s(i+1:) = repeat("0", sig - i + 1)
    else
      s(n:n) = digit(j+1:j+1)
    endif

  end function real_dp_str


#endif

  function str_real_dp_fmt_chk(x, fmt) result(s)
    real(dp), intent(in) :: x
    character(len=*), intent(in) :: fmt
#ifdef DUMMYLIB
    character(len=1) :: s
    s = " "
#else
    character(len=len(x, fmt)) :: s

    if (checkFmt(fmt)) then
      s = safestr(x, fmt)
    else
      call FoX_error("Invalid format: "//fmt)
    endif
#endif
  end function str_real_dp_fmt_chk

#ifndef DUMMYLIB
  pure function str_real_dp_fmt(x, fmt) result(s)
    real(dp), intent(in) :: x
    character(len=*), intent(in) :: fmt
    character(len=len(x, fmt)) :: s

    integer :: sig, dec
    integer :: e, n
    character(len=len(x, fmt)) :: num !this will always be enough memory.

    if (x == 0.0_dp) then
      e = 0
    else
      e = floor(log10(abs(x)))
    endif

    if (x < 0.0_dp) then
      s(1:1) = "-"
      n = 2
    else
      n = 1
    endif

    if (len(fmt) == 0) then

      sig = sig_dp

      num = real_dp_str(abs(x), sig)
      if (num(1:1) == "!") then
        e = e + 1
        num = "1"//repeat("0",len(num)-1)
      endif

      if (sig == 1) then
        s(n:n) = num
        n = n + 1
      else
        s(n:n+1) = num(1:1)//"."
        s(n+2:n+sig) = num(2:)
        n = n + sig + 1
      endif

      s(n:n) = "e"
      s(n+1:) = safestr(e)

    elseif (fmt(1:1) == "s") then

      if (len(fmt) > 1) then
        sig = str_to_int_10(fmt(2:))
      else
        sig = sig_dp
      endif
      sig = max(sig, 1)
      sig = min(sig, digits(1.0_dp))

      num = real_dp_str(abs(x), sig)
      if (num(1:1) == "!") then
        e = e + 1
        num = "1"//repeat("0",len(num)-1)
      endif

      if (sig == 1) then
        s(n:n) = num
        n = n + 1
      else
        s(n:n+1) = num(1:1)//"."
        s(n+2:n+sig) = num(2:)
        n = n + sig + 1
      endif

      s(n:n) = "e"
      s(n+1:) = safestr(e)

    elseif (fmt(1:1) == "r") then

      if (len(fmt) > 1) then
        dec = str_to_int_10(fmt(2:))
      else
        dec = sig_dp - e - 1
      endif
      dec = min(dec, digits(1.0_dp)-e-1)
      dec = max(dec, 0)

      if (e+dec+1 > 0) then
        num = real_dp_str(abs(x), e+dec+1)
      else
        num = ""
      endif
      if (num(1:1) == "!") then
        e = e + 1
        num = "1"//repeat("0",len(num)-1)
      endif

      if (abs(x) >= 1.0_dp) then
        s(n:n+e) = num(:e+1)
        n = n + e + 1
        if (dec > 0) then
          s(n:n) = "."
          n = n + 1
          s(n:) = num(e+2:)
        endif
      else
        s(n:n) = "0"
        if (dec > 0) then
          s(n+1:n+1) = "."
          n = n + 2
          if (dec < -e-1) then
            s(n:) = repeat("0", dec)
          else
            s(n:n-e-2) = repeat("0", max(-e-1,0))
            n = n - min(e,-1) - 1
            if (n <= len(s)) then
              s(n:) = num
            endif
          endif
        endif
      endif

    endif

  end function str_real_dp_fmt

#endif

  pure function str_real_dp(x) result(s)
    real(dp), intent(in) :: x
#ifdef DUMMYLIB
    character(len=1) :: s
    s = " "
#else
    character(len=len(x)) :: s

    s = safestr(x, "")
#endif
  end function str_real_dp

  pure function str_real_dp_array(xa) result(s)
    real(dp), dimension(:), intent(in) :: xa
#ifdef DUMMYLIB
    character(len=1) :: s
    s = " "
#else
    character(len=len(xa)) :: s
    
    integer :: j, k, n

    n = 1
    do k = 1, size(xa) - 1
      j = len(xa(k), "")
      s(n:n+j) = safestr(xa(k), "")//" "
      n = n + j + 1
    enddo
    s(n:) = safestr(xa(k))
#endif
  end function str_real_dp_array

#ifndef DUMMYLIB
  pure function str_real_dp_array_fmt(xa, fmt) result(s)
    real(dp), dimension(:), intent(in) :: xa
    character(len=*), intent(in) :: fmt
    character(len=len(xa, fmt)) :: s
    
    integer :: j, k, n

    n = 1
    do k = 1, size(xa) - 1
      j = len(xa(k), fmt)
      s(n:n+j) = safestr(xa(k), fmt)//" "
      n = n + j + 1
    enddo
    s(n:) = safestr(xa(k), fmt)

  end function str_real_dp_array_fmt
#endif

  function str_real_dp_array_fmt_chk(xa, fmt) result(s)
    real(dp), dimension(:), intent(in) :: xa
    character(len=*), intent(in) :: fmt
#ifdef DUMMYLIB
    character(len=1) :: s
    s = " "
#else
    character(len=len(xa, fmt)) :: s
    
    if (checkFmt(fmt)) then
      s = safestr(xa, fmt)
    else
      call FoX_error("Invalid format: "//fmt)
    endif
#endif
  end function str_real_dp_array_fmt_chk

#ifndef DUMMYLIB
  function str_real_dp_matrix_fmt(xa, fmt) result(s)
    real(dp), dimension(:,:), intent(in) :: xa
    character(len=*), intent(in) :: fmt
    character(len=len(xa,fmt)) :: s

    integer :: i, j, k, n

    i = len(xa(1,1), fmt)
    s(:i) = safestr(xa(1,1), fmt)
    n = i + 1
    do j = 2, size(xa, 1)
      i = len(xa(j,1), fmt)
      s(n:n+i) = " "//safestr(xa(j,1), fmt)
      n = n + i + 1
    enddo
    do k = 2, size(xa, 2)
      do j = 1, size(xa, 1)
        i = len(xa(j,k), fmt)
        s(n:n+i) = " "//safestr(xa(j,k), fmt)
        n = n + i + 1
      enddo
    enddo

  end function str_real_dp_matrix_fmt
#endif

  function str_real_dp_matrix_fmt_chk(xa, fmt) result(s)
    real(dp), dimension(:,:), intent(in) :: xa
    character(len=*), intent(in) :: fmt
#ifdef DUMMYLIB
    character(len=1) :: s
    s = " "
#else
    character(len=len(xa,fmt)) :: s

    if (checkFmt(fmt)) then
      s = safestr(xa, fmt)
    else
      call FoX_error("Invalid format: "//fmt)
    end if
#endif
  end function str_real_dp_matrix_fmt_chk

  function str_real_dp_matrix(xa) result(s)
    real(dp), dimension(:,:), intent(in) :: xa
#ifdef DUMMYLIB
    character(len=1) :: s
    s = " "
#else
    character(len=len(xa)) :: s

    s = safestr(xa, "")
#endif
  end function str_real_dp_matrix

! For complex numbers, there's not really much prior art, so
! we use the easy solution: a+bi, where a & b are real numbers
! as output above.

  function str_complex_sp_fmt_chk(c, fmt) result(s)
    complex(sp), intent(in) :: c
    character(len=*), intent(in) :: fmt
#ifdef DUMMYLIB
    character(len=1) :: s
    s = " "
#else
    character(len=len(c, fmt)) :: s
    
    if (checkFmt(fmt)) then
      s = safestr(c, fmt)
    else
      call FoX_error("Invalid format: "//fmt)
    endif
#endif
  end function str_complex_sp_fmt_chk

#ifndef DUMMYLIB
  pure function str_complex_sp_fmt(c, fmt) result(s)
    complex(sp), intent(in) :: c
    character(len=*), intent(in) :: fmt
    character(len=len(c, fmt)) :: s
    
    real(sp) :: re, im
    integer :: i
    re = real(c)
    im = aimag(c)
    i = len(re, fmt)
    s(:i+4) = "("//safestr(re, fmt)//")+i"
    s(i+5:)="("//safestr(im,fmt)//")"
  end function str_complex_sp_fmt
#endif

  pure function str_complex_sp(c) result(s)
    complex(sp), intent(in) :: c
#ifdef DUMMYLIB
    character(len=1) :: s
    s = " "
#else
    character(len=len(c, "")) :: s

    s = safestr(c, "")
#endif
  end function str_complex_sp

#ifndef DUMMYLIB
  pure function str_complex_sp_array_fmt(ca, fmt) result(s)
    complex(sp), dimension(:), intent(in) :: ca
    character(len=*), intent(in) :: fmt
    character(len=len(ca, fmt)) :: s

    integer :: i, n
 
    s(1:len(ca(1), fmt)) = safestr(ca(1), fmt)
    n = len(ca(1), fmt)+1
    do i = 2, size(ca) 
      s(n:n+len(ca(i), fmt)) = " "//safestr(ca(i), fmt)
      n = n + len(ca(i), fmt)+1
    enddo
  end function str_complex_sp_array_fmt
#endif

  function str_complex_sp_array_fmt_chk(ca, fmt) result(s)
    complex(sp), dimension(:), intent(in) :: ca
    character(len=*), intent(in) :: fmt
#ifdef DUMMYLIB
    character(len=1) :: s
    s = " "
#else
    character(len=len(ca, fmt)) :: s

    if (checkFmt(fmt)) then
      s = safestr(ca, fmt)
    else
      call FoX_error("Invalid format: "//fmt)
    endif
#endif
  end function str_complex_sp_array_fmt_chk

  pure function str_complex_sp_array(ca) result(s)
    complex(sp), dimension(:), intent(in) :: ca
#ifdef DUMMYLIB
    character(len=1) :: s
    s = " "
#else
    character(len=len(ca)) :: s

    s = safestr(ca, "")
#endif
  end function str_complex_sp_array

#ifndef DUMMYLIB
  pure function str_complex_sp_matrix_fmt(ca, fmt) result(s)
    complex(sp), dimension(:, :), intent(in) :: ca
    character(len=*), intent(in) :: fmt
    character(len=len(ca, fmt)) :: s

    integer :: i, j, k, n

    i = len(ca(1,1), fmt)
    s(:i) = safestr(ca(1,1), fmt)
    n = i + 1
    do j = 2, size(ca, 1)
      i = len(ca(j,1), fmt)
      s(n:n+i) = " "//safestr(ca(j,1), fmt)
      n = n + i + 1
    enddo
    do k = 2, size(ca, 2)
      do j = 1, size(ca, 1)
        i = len(ca(j,k), fmt)
        s(n:n+i) = " "//safestr(ca(j,k), fmt)
        n = n + i + 1
      enddo
    enddo

  end function str_complex_sp_matrix_fmt
#endif

  function str_complex_sp_matrix_fmt_chk(ca, fmt) result(s)
    complex(sp), dimension(:, :), intent(in) :: ca
    character(len=*), intent(in) :: fmt
#ifdef DUMMYLIB
    character(len=1) :: s
    s = " "
#else
    character(len=len(ca, fmt)) :: s

    if (checkFmt(fmt)) then
      s = safestr(ca, fmt)
    else
      call FoX_error("Invalid format: "//fmt)
    endif
#endif
  end function str_complex_sp_matrix_fmt_chk

  pure function str_complex_sp_matrix(ca) result(s)
    complex(sp), dimension(:, :), intent(in) :: ca
#ifdef DUMMYLIB
    character(len=1) :: s
    s = " "
#else
    character(len=len(ca)) :: s

    s = safestr(ca, "")
#endif
  end function str_complex_sp_matrix

  function str_complex_dp_fmt_chk(c, fmt) result(s)
    complex(dp), intent(in) :: c
    character(len=*), intent(in) :: fmt
#ifdef DUMMYLIB
    character(len=1) :: s
    s = " "
#else
    character(len=len(c, fmt)) :: s
    
    if (checkFmt(fmt)) then
      s = safestr(c, fmt)
    else
      call FoX_error("Invalid format: "//fmt)
    endif
#endif
  end function str_complex_dp_fmt_chk

#ifndef DUMMYLIB
  pure function str_complex_dp_fmt(c, fmt) result(s)
    complex(dp), intent(in) :: c
    character(len=*), intent(in) :: fmt
    character(len=len(c, fmt)) :: s
    
    real(dp) :: re, im
    integer :: i
    re = real(c)
    im = aimag(c)
    i = len(re, fmt)
    s(:i+4) = "("//safestr(re, fmt)//")+i"
    s(i+5:)="("//safestr(im,fmt)//")"
  end function str_complex_dp_fmt
#endif

  pure function str_complex_dp(c) result(s)
    complex(dp), intent(in) :: c
#ifdef DUMMYLIB
    character(len=1) :: s
    s = " "
#else
    character(len=len(c, "")) :: s

    s = safestr(c, "")
#endif
  end function str_complex_dp

#ifndef DUMMYLIB
  pure function str_complex_dp_array_fmt(ca, fmt) result(s)
    complex(dp), dimension(:), intent(in) :: ca
    character(len=*), intent(in) :: fmt
    character(len=len(ca, fmt)) :: s

    integer :: i, n

    s(1:len(ca(1), fmt)) = safestr(ca(1), fmt)
    n = len(ca(1), fmt)+1
    do i = 2, size(ca) 
      s(n:n+len(ca(i), fmt)) = " "//safestr(ca(i), fmt)
      n = n + len(ca(i), fmt)+1
    enddo
  end function str_complex_dp_array_fmt
#endif

  function str_complex_dp_array_fmt_chk(ca, fmt) result(s)
    complex(dp), dimension(:), intent(in) :: ca
    character(len=*), intent(in) :: fmt
#ifdef DUMMYLIB
    character(len=1) :: s
    s = " "
#else
    character(len=len(ca, fmt)) :: s

    if (checkFmt(fmt)) then
      s = safestr(ca, fmt)
    else
      call FoX_error("Invalid format: "//fmt)
    endif
#endif
  end function str_complex_dp_array_fmt_chk

  pure function str_complex_dp_array(ca) result(s)
    complex(dp), dimension(:), intent(in) :: ca
#ifdef DUMMYLIB
    character(len=1) :: s
    s = " "
#else
    character(len=len(ca)) :: s

    s = safestr(ca, "")
#endif
  end function str_complex_dp_array

#ifndef DUMMYLIB
  pure function str_complex_dp_matrix_fmt(ca, fmt) result(s)
    complex(dp), dimension(:, :), intent(in) :: ca
    character(len=*), intent(in) :: fmt
    character(len=len(ca, fmt)) :: s

    integer :: i, j, k, n

    i = len(ca(1,1), fmt)
    s(:i) = safestr(ca(1,1), fmt)
    n = i + 1
    do j = 2, size(ca, 1)
      i = len(ca(j,1), fmt)
      s(n:n+i) = " "//safestr(ca(j,1), fmt)
      n = n + i + 1
    enddo
    do k = 2, size(ca, 2)
      do j = 1, size(ca, 1)
        i = len(ca(j,k), fmt)
        s(n:n+i) = " "//safestr(ca(j,k), fmt)
        n = n + i + 1
      enddo
    enddo

  end function str_complex_dp_matrix_fmt
#endif

  function str_complex_dp_matrix_fmt_chk(ca, fmt) result(s)
    complex(dp), dimension(:, :), intent(in) :: ca
    character(len=*), intent(in) :: fmt
#ifdef DUMMYLIB
    character(len=1) :: s
    s = " "
#else
    character(len=len(ca, fmt)) :: s

    if (checkFmt(fmt)) then
      s = safestr(ca, fmt)
    else
      call FoX_error("Invalid format: "//fmt)
    endif
#endif
  end function str_complex_dp_matrix_fmt_chk

  pure function str_complex_dp_matrix(ca) result(s)
    complex(dp), dimension(:, :), intent(in) :: ca
#ifdef DUMMYLIB
    character(len=1) :: s
    s = " "
#else
    character(len=len(ca)) :: s

    s = safestr(ca, "")
#endif
  end function str_complex_dp_matrix

#ifndef DUMMYLIB
  pure function checkFmt(fmt) result(good)
    character(len=*), intent(in) :: fmt
    logical :: good

    ! should be ([rs]\d*)?

    if (len(fmt) > 0) then
      if (fmt(1:1) == "r" .or. fmt(1:1) == "s") then
        if (len(fmt) > 1) then
          good = (verify(fmt(2:), digit) == 0)
        else
          good = .true.
        endif
      else
        good = .false.
      endif
    else
      good = .true.
    endif
  end function checkFmt
#endif

  pure function concat_str_int(s1, s2) result(s3)
    character(len=*), intent(in) :: s1
    integer, intent(in) :: s2
#ifdef DUMMYLIB
    character(len=1) :: s3
    s3 = " "
#else
    character(len=len(s1)+len(s2)) :: s3
    s3 = s1//str(s2)
#endif
  end function concat_str_int
  pure function concat_int_str(s1, s2) result(s3)
    integer, intent(in) :: s1
    character(len=*), intent(in) :: s2
#ifdef DUMMYLIB
    character(len=1) :: s3
    s3 = " "
#else
    character(len=len(s1)+len(s2)) :: s3
    s3 = str(s1)//s2
#endif
  end function concat_int_str

  pure function concat_str_logical(s1, s2) result(s3)
    character(len=*), intent(in) :: s1
    logical, intent(in) :: s2
#ifdef DUMMYLIB
    character(len=1) :: s3
    s3 = " "
#else
    character(len=len(s1)+len(s2)) :: s3
    s3 = s1//str(s2)
#endif
  end function concat_str_logical
  pure function concat_logical_str(s1, s2) result(s3)
    logical, intent(in) :: s1
    character(len=*), intent(in) :: s2
#ifdef DUMMYLIB
    character(len=1) :: s3
    s3 = " "
#else
    character(len=len(s1)+len(s2)) :: s3
    s3 = str(s1)//s2
#endif
  end function concat_logical_str

  pure function concat_str_real_sp(s1, s2) result(s3)
    character(len=*), intent(in) :: s1
    real(sp), intent(in) :: s2
#ifdef DUMMYLIB
    character(len=1) :: s3
    s3 = " "
#else
    character(len=len(s1)+len(s2)) :: s3
    s3 = s1//str(s2)
#endif
  end function concat_str_real_sp
  pure function concat_real_sp_str(s1, s2) result(s3)
    real(sp), intent(in) :: s1
    character(len=*), intent(in) :: s2
#ifdef DUMMYLIB
    character(len=1) :: s3
    s3 = " "
#else
    character(len=len(s1)+len(s2)) :: s3
    s3 = str(s1)//s2
#endif
  end function concat_real_sp_str

  pure function concat_str_real_dp(s1, s2) result(s3)
    character(len=*), intent(in) :: s1
    real(dp), intent(in) :: s2
#ifdef DUMMYLIB
    character(len=1) :: s3
    s3 = " "
#else
    character(len=len(s1)+len(s2)) :: s3
    s3 = s1//str(s2)
#endif
  end function concat_str_real_dp
  pure function concat_real_dp_str(s1, s2) result(s3)
    real(dp), intent(in) :: s1
    character(len=*), intent(in) :: s2
#ifdef DUMMYLIB
    character(len=1) :: s3
    s3 = " "
#else
    character(len=len(s1)+len(s2)) :: s3
    s3 = str(s1)//s2
#endif
  end function concat_real_dp_str

  pure function concat_str_complex_sp(s1, s2) result(s3)
    character(len=*), intent(in) :: s1
    complex(sp), intent(in) :: s2
#ifdef DUMMYLIB
    character(len=1) :: s3
    s3 = " "
#else
    character(len=len(s1)+len(s2)) :: s3
    s3 = s1//str(s2)
#endif
  end function concat_str_complex_sp
  pure function concat_complex_sp_str(s1, s2) result(s3)
    complex(sp), intent(in) :: s1
    character(len=*), intent(in) :: s2
#ifdef DUMMYLIB
    character(len=1) :: s3
    s3 = " "
#else
    character(len=len(s1)+len(s2)) :: s3
    s3 = str(s1)//s2
#endif
  end function concat_complex_sp_str

  pure function concat_str_complex_dp(s1, s2) result(s3)
    character(len=*), intent(in) :: s1
    complex(dp), intent(in) :: s2
#ifdef DUMMYLIB
    character(len=1) :: s3
    s3 = " "
#else
    character(len=len(s1)+len(s2)) :: s3
    s3 = s1//str(s2)
#endif
  end function concat_str_complex_dp
  pure function concat_complex_dp_str(s1, s2) result(s3)
    complex(dp), intent(in) :: s1
    character(len=*), intent(in) :: s2
#ifdef DUMMYLIB
    character(len=1) :: s3
    s3 = " "
#else
    character(len=len(s1)+len(s2)) :: s3
    s3 = str(s1)//s2
#endif
  end function concat_complex_dp_str

end module fox_m_fsys_format
