module m_common_buffer

#ifndef DUMMYLIB
  use fox_m_fsys_format, only: str
  use m_common_charset, only: XML1_0
  use m_common_error, only: FoX_error, FoX_warning

  implicit none
  private
  
  ! At this point we use a fixed-size buffer. 
  ! Note however that buffer overflows will only be
  ! triggered by overly long *unbroken* pcdata values, or
  ! by overly long attribute values. Hopefully
  ! element or attribute names are "short enough".
  !
  ! In a forthcoming implementation it could be made dynamical...
  
  ! MAX_BUFF_SIZE cannot be bigger than the maximum available
  ! record length for a compiler. In practice, this means
  ! 1024 seems to be the biggest available size.
  
  integer, parameter :: MAX_BUFF_SIZE  = 1024
  
  type buffer_t
    private
    integer                       :: size
    character(len=MAX_BUFF_SIZE)  :: str
    integer                       :: unit
    integer                       :: xml_version
  end type buffer_t
  
  public :: buffer_t
  
  public :: add_to_buffer
  public :: print_buffer, str, char, len
  public :: buffer_to_chararray
  public :: reset_buffer
  public :: dump_buffer

  interface str
    module procedure buffer_to_str
  end interface
  
  interface char
    module procedure buffer_to_str
  end interface
  
  interface len
    module procedure buffer_length
  end interface
  
contains

  subroutine reset_buffer(buffer, unit, xml_version)
    type(buffer_t), intent(inout)  :: buffer
    integer, intent(in), optional :: unit
    integer, intent(in) :: xml_version

    buffer%size = 0
    if (present(unit)) then
      buffer%unit = unit
    else 
      buffer%unit = 6
    endif
    buffer%xml_version = xml_version
    
  end subroutine reset_buffer
  

  subroutine print_buffer(buffer)
    type(buffer_t), intent(in)  :: buffer
    
    write(unit=6,fmt="(a)") buffer%str(:buffer%size)

  end subroutine print_buffer


  function buffer_to_str(buffer) result(str)
    type(buffer_t), intent(in)          :: buffer
    character(len=buffer%size)          :: str
    
    str = buffer%str(:buffer%size)
  end function buffer_to_str


  function buffer_to_chararray(buffer) result(str)
    type(buffer_t), intent(in)               :: buffer
    character(len=1), dimension(buffer%size) :: str
    integer :: i
    
    do i = 1, buffer%size
      str(i) = buffer%str(i:i)
    enddo
  end function buffer_to_chararray


  function buffer_length(buffer) result(length)
    type(buffer_t), intent(in)          :: buffer
    integer                             :: length
    
    length = buffer%size 

  end function buffer_length

  
  subroutine dump_buffer(buffer, lf)
    type(buffer_t), intent(inout) :: buffer
    logical, intent(in), optional :: lf

    integer :: i, n
    logical :: lf_

    if (present(lf)) then
      lf_ = lf
    else
      lf_ = .true.
    endif

    i = scan(buffer%str(:buffer%size), achar(10)//achar(13))
    n = 1
    do while (i>0)
      write(buffer%unit, '(a)', advance="yes") buffer%str(n:n+i-2)
      n = n + i
      if (n>buffer%size) exit
      i = scan(buffer%str(n:), achar(10)//achar(13))
    enddo
    
    if (n<=buffer%size) then
      if (lf_) then
        write(buffer%unit, '(a)', advance="yes") buffer%str(n:buffer%size)
      else
        write(buffer%unit, '(a)', advance="no") buffer%str(n:buffer%size)
      endif
    endif

    buffer%size = 0
  end subroutine dump_buffer

  
  subroutine check_buffer(s, version)
    character(len=*), intent(in) :: s
    integer, intent(in) :: version

    integer :: i

!FIXME this is almost a duplicate of logic in wxml/m_wxml_escape.f90

    ! We have to do it this way (with achar etc) in case the native
    ! platform encoding is not ASCII

    do i = 1, len(s)
      select case (iachar(s(i:i)))
      case (0)
        call FoX_error("Tried to output a NUL character")
      case (1:8,11:12,14:31)
        if (version==XML1_0) then
          call FoX_error("Tried to output a character invalid under XML 1.0: &#"//str(iachar(s(i:i)))//";")
        endif
      case (128:)
        !TOHW we should maybe just disallow this ...
        call FoX_warning("emitting non-ASCII character. Platform-dependent result!")
      end select
    enddo 

  end subroutine check_buffer


  subroutine add_to_buffer(s, buffer, ws_significant)
    character(len=*), intent(in)   :: s
    type(buffer_t), intent(inout)  :: buffer
    logical, intent(in), optional :: ws_significant

    character(len=(buffer%size+len(s))) :: s2
    integer :: i, n, len_b
    logical :: warning, ws_

    ! Is whitespace significant in this context?
    ! We have to assume so unless told otherwise.
    if (present(ws_significant)) then
      ws_ = ws_significant
    else
      ws_ = .true.
    endif

    ! FIXME The algorithm below unilaterally forces all
    ! line feeds and carriage returns to native EOL, regardless
    ! of input document. Thus it is impossible to output a
    ! document containing a literal non-native newline character
    ! Ideally we would put this under the control of the user.

    ! We check if whitespace is significant. If not, we can
    ! adjust the buffer without worrying about it.
    ! But if we are not told, we warn about it.
    ! And if we are told it definitely is - then we error out.

    ! If we overreach our buffer size, we will be unable to
    ! output any more characters without a newline.
    ! Go through new string, insert newlines
    ! at spaces just before MAX_BUFF_SIZE chars
    ! until we have less than MAX_BUFF_SIZE left to go,
    ! then put that in the buffer.

    ! If no whitespace is found in the newly-added string, then
    ! insert a new line immediately before it (at the end of the
    ! current buffer)

    call check_buffer(s, buffer%xml_version)

    s2 = buffer%str(:buffer%size)//s


    ! output as much of this using existing newlines as possible.
    warning = .false.
    n = 1
    do while (n<=len(s2))
      ! Note this is an XML-1.0 only definition of newline
      i = scan(s2(n:), achar(10)//achar(13))
      if (i>0) then
        ! turn that newline into an output newline ...
        write(buffer%unit, '(a)') s2(n:n+i-2)
        n = n + i
      elseif (n<=len(s2)-MAX_BUFF_SIZE) then
        ! We need to insert a newline, or we'll overrun the buffer
        ! No suitable newline, so convert a space or tab into a newline.
        i = scan(s2(n:n+MAX_BUFF_SIZE-1), achar(9)//achar(32), back=.true.)
        ! If no space or tab is present, we fail.
        if (i>0.and..not.present(ws_significant)) then
          ! We can insert a newline, but we don't know whether it is    significant. Warn:
          if (.not.warning) then
            ! We only output this warning once.
            call FoX_warning( &
            "Fortran made FoX insert a newline. "// &
            "If whitespace might be significant, check your output.")
            warning = .true.
          endif
        elseif (i==0) then
          call FoX_error( &
            "Fortran made FoX insert a newline but it can't. Stopping now.")
        elseif (ws_) then
          call FoX_error( &
            "Fortran made FoX insert a newline but whitespace is  significant. Stopping now.")
        else
          continue ! without error or warning, because whitespace is not significant
        endif
        write(buffer%unit, '(a)') s2(n:n+i-1)
        n = n + i
      else
        ! We don't need to do anything, just add the remainder to the buffer.
        exit
      endif
    enddo

    len_b = len(s2) - n + 1
    buffer%str(:len_b) = s2(n:)
    buffer%size = len_b

  end subroutine add_to_buffer

#endif
end module m_common_buffer
