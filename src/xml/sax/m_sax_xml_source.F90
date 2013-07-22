module m_sax_xml_source

#ifndef DUMMYLIB
  use fox_m_fsys_array_str, only: str_vs, vs_str_alloc
  use fox_m_fsys_format, only: operator(//)
  use m_common_error,  only: error_stack, add_error, in_error
  use m_common_charset, only: XML_WHITESPACE, XML_INITIALENCODINGCHARS, &
    XML_ENCODINGCHARS, XML1_0, XML1_1, isXML1_0_NameChar, &
    isLegalChar, isUSASCII, allowed_encoding
  use m_common_io, only: io_eor, io_eof
  
  use FoX_utils, only: URI

  implicit none
  private

  type buffer_t
    character, dimension(:), pointer :: s
    integer :: pos = 1
  end type buffer_t

  type xml_source_t
    !FIXME private
    integer            :: lun = -1
    integer            :: xml_version = XML1_0
    character, pointer :: encoding(:) => null()
    logical            :: isUSASCII
    character, pointer :: filename(:) => null()
    type(URI), pointer :: baseURI => null()
    integer            :: line = 0
    integer            :: col = 0
    integer            :: startChar = 1 ! First character after XML decl
    character, pointer :: next_chars(:) => null()   ! pushback buffer
    type(buffer_t), pointer :: input_string => null()
    logical :: pe = .false. ! is this a parameter entity?
    logical :: eof = .false.! need to keep track of this at the end of pes
  end type xml_source_t

  public :: buffer_t
  public :: xml_source_t

  public :: get_char_from_file
  public :: push_file_chars
  public :: parse_declaration

contains


  function get_char_from_file(f, xv, eof, es) result(string)
    type(xml_source_t), intent(inout) :: f
    integer, intent(in) :: xv
    logical, intent(out) :: eof
    type(error_stack), intent(inout) :: es
    character(len=1) :: string

    integer :: iostat
    logical :: pending
    character :: c, c2

    pending = .false.
    eof = .false.
    c = read_single_char(f, iostat)
    if (iostat==io_eof) then
      eof = .true.
      return
    elseif (iostat/=0) then
      call add_error(es, "Error reading "//str_vs(f%filename))
      return
    endif
    if (.not.isLegalChar(c, f%isUSASCII, xv)) then
      call add_error(es, "Illegal character found at " &
        //str_vs(f%filename)//":"//f%line//":"//f%col)
      return
    endif
    if (c==achar(13)) then
      c = achar(10)
      c2 = read_single_char(f, iostat)
      if (iostat==io_eof) then
        ! the file has just ended on a single CR. Report is as a LF.
        ! Ignore the eof just now, it'll be picked up if we need to 
        ! perform another read.
        eof = .false.
      elseif (iostat/=0) then
        call add_error(es, "Error reading "//str_vs(f%filename))
        return
      elseif (c2/=achar(10)) then
        ! then we keep c2, otherwise we'd just ignore it.
        pending = .true.
      endif
    endif
    string = c

    if (pending) then
      ! we have one character left over, put in the pushback buffer
      deallocate(f%next_chars)
      allocate(f%next_chars(1))
      f%next_chars = c2
    endif

    if (c==achar(10)) then
      f%line = f%line + 1
      f%col = 0
    else
      f%col = f%col + 1
    endif

  end function get_char_from_file

  function read_single_char(f, iostat) result(c)
    type(xml_source_t), intent(inout) :: f
    integer, intent(out) :: iostat
    character :: c

    if (f%eof) then
      c = ""
      iostat = io_eof
      return
    endif
    if (f%lun==-1) then
      if (f%input_string%pos>size(f%input_string%s)) then
        c = ""
        if (f%pe) then
          iostat = 0
        else
          iostat = io_eof
        endif
        f%eof = .true.
      else
        iostat = 0
        c = f%input_string%s(f%input_string%pos)
        f%input_string%pos = f%input_string%pos + 1
      endif
    else
      read (unit=f%lun, iostat=iostat, advance="no", fmt="(a1)") c
      if (iostat==io_eor) then
        iostat = 0
#ifdef FC_EOR_LF
        c = achar(10)
#else
        c = achar(13)
#endif
      elseif (iostat==io_eof) then
        if (f%pe) iostat = 0
        c = ""
        f%eof = .true.
      endif
    endif
  end function read_single_char

  subroutine rewind_source(f)
    type(xml_source_t), intent(inout) :: f

    if (f%lun==-1) then
      f%input_string%pos = 1
    else
      rewind(f%lun)
    endif
  end subroutine rewind_source

  subroutine push_file_chars(f, s)
    type(xml_source_t), intent(inout) :: f
    character(len=*), intent(in) :: s
    character, dimension(:), pointer :: nc

    nc => vs_str_alloc(s//str_vs(f%next_chars))
    deallocate(f%next_chars)
    f%next_chars => nc

  end subroutine push_file_chars


  subroutine parse_declaration(f, eof, es, standalone)
    type(xml_source_t), intent(inout) :: f
    logical, intent(out) :: eof
    type(error_stack), intent(inout) :: es
    logical, intent(out), optional :: standalone

    integer :: parse_state, xd_par
    character :: c, q
    character, pointer :: ch(:), ch2(:)

    integer, parameter :: XD_0      = 0
    integer, parameter :: XD_START  = 1
    integer, parameter :: XD_TARGET = 2
    integer, parameter :: XD_MISC   = 3
    integer, parameter :: XD_PA     = 4
    integer, parameter :: XD_EQ     = 5
    integer, parameter :: XD_QUOTE  = 6
    integer, parameter :: XD_PV     = 7
    integer, parameter :: XD_END    = 8
    integer, parameter :: XD_SPACE  = 9

    integer, parameter :: xd_nothing = 0
    integer, parameter :: xd_version = 1
    integer, parameter :: xd_encoding = 2
    integer, parameter :: xd_standalone = 3

    f%xml_version = XML1_0
    if (present(standalone)) standalone = .false.

    f%startChar = 1

    parse_state = XD_0
    xd_par = xd_nothing
    ch => null()
    do
      c = get_char_from_file(f, XML1_0, eof, es)
      if (eof) then
        call rewind_source(f)
        exit
      elseif (in_error(es)) then
        goto 100
      endif
      f%startChar = f%startChar + 1

      select case (parse_state)

      case (XD_0)
        if (c=="<") then
          parse_state = XD_START
        else
          call rewind_source(f)
          exit
        endif

      case (XD_START)
        if (c=="?") then
          parse_state = XD_TARGET
          ch => vs_str_alloc("")
        else
          call rewind_source(f)
          exit
        endif

      case (XD_TARGET)
        if (isXML1_0_NameChar(c)) then
          ch2 => vs_str_alloc(str_vs(ch)//c)
          deallocate(ch)
          ch => ch2
        elseif (verify(c, XML_WHITESPACE)==0 &
          .and.str_vs(ch)=="xml") then
          deallocate(ch)
          parse_state = XD_MISC
        else
          call rewind_source(f)
          deallocate(ch)
          exit
        endif

      case (XD_SPACE)
        if (verify(c, XML_WHITESPACE)==0) then
          parse_state = XD_MISC
        elseif (c=="?") then
          parse_state = XD_END
        else
          call add_error(es, &
            "Missing space in XML declaration")
        endif

      case (XD_MISC)
        if (c=="?") then
          parse_state = XD_END
        elseif (isXML1_0_NameChar(c)) then
          ch => vs_str_alloc(c)
          parse_state = XD_PA
        elseif (verify(c, XML_WHITESPACE)>0) then
          call add_error(es, &
            "Unexpected character in XML declaration")
        endif

      case (XD_PA)
        if (isXML1_0_NameChar(c)) then
          ch2 => vs_str_alloc(str_vs(ch)//c)
          deallocate(ch)
          ch => ch2
        elseif (verify(c, XML_WHITESPACE//"=")==0) then
          select case (str_vs(ch))

          case ("version")
            select case (xd_par)
            case (xd_nothing)
              xd_par = xd_version
            case default
              call add_error(es, &
                "Cannot specify version twice in XML declaration")
            end select

          case ("encoding")
            select case (xd_par)
            case (xd_nothing)
              if (present(standalone)) then
                call add_error(es, &
                  "Must specify version before encoding in XML declaration")
              else
                xd_par = xd_encoding
              endif
            case (xd_version)
              xd_par = xd_encoding
            case (xd_encoding)
              call add_error(es, &
                "Cannot specify encoding twice in XML declaration")
            case (xd_standalone)
              call add_error(es, &
                "Cannot specify encoding after standalone in XML declaration")
            end select

          case ("standalone")
            if (.not.present(standalone)) &
              call add_error(es, &
              "Cannot specify standalone in text declaration")
            select case (xd_par)
            case (xd_nothing)
              call add_error(es, &
                "Must specify version before standalone in XML declaration")
            case (xd_version, xd_encoding)
              xd_par = xd_standalone
            case (xd_standalone)
              call add_error(es, &
                "Cannot specify standalone twice in XML declaration")
            end select

          case default
            call add_error(es, &
              "Unknown parameter "//str_vs(ch)//" in XML declaration, "//&
              "expecting version, encoding or standalone")

          end select

          deallocate(ch)
          if (c=="=") then
            parse_state = XD_QUOTE
          else
            parse_state = XD_EQ
          endif
        else
          call add_error(es, &
            "Unexpected character found in XML declaration")
        endif

      case (XD_EQ)
        if (c=="=") then
          parse_state = XD_QUOTE
        elseif (verify(c, XML_WHITESPACE)>0) then
          call add_error(es, &
            "Unexpected character found in XML declaration; expecting ""=""")
        endif

      case (XD_QUOTE)
        if (verify(c, "'""")==0) then
          q = c
          parse_state = XD_PV
          ch => vs_str_alloc("")
        elseif (verify(c, XML_WHITESPACE)>0) then
          call add_error(es, &
            "Unexpected character found in XML declaration; expecting "" or '")
        endif

      case (XD_PV)
        if (c==q) then
          select case (xd_par)
          case (xd_version)
            if (str_vs(ch)//"x"=="1.0x") then
              f%xml_version = XML1_0
              deallocate(ch)
            elseif (str_vs(ch)//"x"=="1.1x") then
              f%xml_version = XML1_1
              deallocate(ch)
            else
              call add_error(es, &
                "Unknown version number "//str_vs(ch)//" found in XML declaration; expecting 1.0 or 1.1")
            endif
          case (xd_encoding)
            if (size(ch)==0) then
              call add_error(es, &
                "Empty value for encoding not allowed in XML declaration")
            elseif (size(ch)==1.and.verify(ch(1), XML_INITIALENCODINGCHARS)>0) then
              call add_error(es, &
                "Invalid encoding found in XML declaration; illegal characters in encoding name")
            elseif (size(ch)>1.and. &
              (verify(ch(1), XML_INITIALENCODINGCHARS)>0 &
              .or.verify(str_vs(ch(2:)), XML_ENCODINGCHARS)>0)) then
              call add_error(es, &
                "Invalid encoding found in XML declaration; illegal characters in encoding name")
            elseif (.not.allowed_encoding(str_vs(ch))) then
              call add_error(es, "Unknown character encoding in XML declaration")
            else
              f%encoding => ch
              f%isUSASCII = isUSASCII(str_vs(ch))
              ch => null()
            endif
          case (xd_standalone)
            if (str_vs(ch)//"x"=="yesx") then
              standalone = .true.
              deallocate(ch)
            elseif (str_vs(ch)//"x"=="nox") then
              standalone = .false.
              deallocate(ch)
            else
              call add_error(es, &
                "Invalid value for standalone found in XML declaration; expecting yes or no")

            endif
          end select
          parse_state = XD_SPACE
        else
          ch2 => vs_str_alloc(str_vs(ch)//c)
          deallocate(ch)
          ch => ch2
        endif

      case (XD_END)
        if (c==">") then
          exit
        else
          call add_error(es, &
            "Unexpected character found in XML declaration; expecting >")
        endif

      end select

    end do

    if (.not.associated(f%encoding)) then
      if (present(standalone).or.parse_state/=XD_END) then
        f%encoding => vs_str_alloc("utf-8")
      else
        call add_error(es, "Missing encoding in text declaration")
      endif
    endif
    
100 if (associated(ch)) deallocate(ch)
    ! if there is no XML declaraion, or if parsing caused an error, then
    if (parse_state/=XD_END.or.in_error(es)) f%startChar = 1

  end subroutine parse_declaration
#endif

end module m_sax_xml_source
