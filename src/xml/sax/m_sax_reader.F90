module m_sax_reader
#ifndef DUMMYLIB

  use fox_m_fsys_array_str, only: str_vs, vs_str_alloc, vs_vs_alloc
  use fox_m_fsys_format, only: operator(//)
  use m_common_charset, only: XML1_0
  use m_common_error,  only: error_stack, FoX_error, in_error, add_error
  use m_common_io, only: setup_io, get_unit, io_err

  use FoX_utils, only: URI, parseURI, copyURI, destroyURI, &
    hasScheme, getScheme, getPath

  use m_sax_xml_source, only: xml_source_t, &
    get_char_from_file, push_file_chars, parse_declaration

  implicit none
  private

  type file_buffer_t
    !FIXME private
    type(xml_source_t), pointer :: f(:) => null()
    logical                     :: standalone = .false.
    integer                     :: xml_version = XML1_0
  end type file_buffer_t

  public :: file_buffer_t
  public :: line
  public :: column

  public :: open_file
  public :: close_file

  public :: open_new_file

  public :: push_chars

  public :: get_character
  public :: get_all_characters

  public :: open_new_string
  public :: pop_buffer_stack

  public :: parse_xml_declaration
  public :: parse_text_declaration

  public :: reading_main_file
  public :: reading_first_entity

contains

  subroutine open_file(fb, iostat, file, lun, string, es)
    type(file_buffer_t), intent(out)  :: fb
    character(len=*), intent(in), optional :: file
    integer, intent(out)              :: iostat
    integer, intent(in), optional     :: lun
    character(len=*), intent(in), optional :: string
    type(error_stack), intent(inout) :: es

    type(URI), pointer :: fileURI

    iostat = 0

    call setup_io()
    if (present(string)) then
      if (present(file)) then
        call FoX_error("Cannot specify both file and string input to open_xml")
      elseif (present(lun)) then
        call FoX_error("Cannot specify lun for string input to open_xml")
      endif
      fileURI => parseURI("")
      call open_new_string(fb, string, "", baseURI=fileURI)
    else
      fileURI => parseURI(file)
      if (.not.associated(fileURI)) then
        call add_error(es, "Could not open file "//file//" - not a valid URI")
        return
      endif
      call open_new_file(fb, fileURI, iostat, lun)
    endif
    call destroyURI(fileURI)

  end subroutine open_file


  subroutine open_new_file(fb, baseURI, iostat, lun, pe)
    type(file_buffer_t), intent(inout)  :: fb
    integer, intent(out) :: iostat
    type(URI), pointer :: baseURI
    integer, intent(in), optional :: lun
    logical, intent(in), optional :: pe

    integer :: i
    type(xml_source_t) :: f
    type(xml_source_t), pointer :: temp(:)
    logical :: pe_

    if (present(pe)) then
      pe_ = pe
    else
      pe_ = .false.
    endif

    if (hasScheme(baseURI)) then
      if (getScheme(baseURI)/="file") then
        iostat = io_err
        return
      endif
    endif

    call open_actual_file(f, getPath(baseURI), iostat, lun)
    if (iostat==0) then
      if (.not.associated(fb%f)) allocate(fb%f(0))
      ! First file

      temp => fb%f
      allocate(fb%f(size(temp)+1))
      do i = 1, size(temp)
        fb%f(i+1)%lun = temp(i)%lun
        fb%f(i+1)%xml_version = temp(i)%xml_version
        fb%f(i+1)%encoding => temp(i)%encoding
        fb%f(i+1)%filename => temp(i)%filename
        fb%f(i+1)%line = temp(i)%line
        fb%f(i+1)%col = temp(i)%col
        fb%f(i+1)%startChar = temp(i)%startChar
        fb%f(i+1)%next_chars => temp(i)%next_chars
        fb%f(i+1)%input_string => temp(i)%input_string
        fb%f(i+1)%baseURI => temp(i)%baseURI
        fb%f(i+1)%pe = temp(i)%pe
      enddo
      deallocate(temp)
      fb%f(1)%lun = f%lun
      fb%f(1)%filename => f%filename
      if (pe_) then
        fb%f(1)%next_chars => vs_str_alloc(" ")
      else
        fb%f(1)%next_chars => vs_str_alloc("")
      endif
      fb%f(1)%pe = pe_
      fb%f(1)%baseURI => copyURI(baseURI)
    endif

  end subroutine open_new_file

  subroutine open_actual_file(f, file, iostat, lun)
    type(xml_source_t), intent(out)    :: f
    character(len=*), intent(in)     :: file
    integer, intent(out)             :: iostat
    integer, intent(in), optional    :: lun

    if (present(lun)) then
      f%lun = lun
    else
      call get_unit(f%lun, iostat)
      if (iostat/=0) return
    endif
    open(unit=f%lun, file=file, form="formatted", status="old", &
      action="read", position="rewind", iostat=iostat)
    if (iostat/=0) return
    f%filename => vs_str_alloc(file)

  end subroutine open_actual_file

  subroutine close_file(fb)
    type(file_buffer_t), intent(inout)  :: fb

    integer :: i

    do i = 1, size(fb%f)
      call close_actual_file(fb%f(i))
    enddo

    deallocate(fb%f)

  end subroutine close_file


  subroutine close_actual_file(f)
    type(xml_source_t), intent(inout)    :: f

    deallocate(f%filename)
    
    if (f%lun>0) then
      close(f%lun)
    else
      deallocate(f%input_string%s)
      deallocate(f%input_string)
    endif

    if (associated(f%encoding)) deallocate(f%encoding)
    f%line = 0
    f%col = 0
    deallocate(f%next_chars)
    call destroyURI(f%baseURI)
  end subroutine close_actual_file


  subroutine open_new_string(fb, string, name, baseURI, pe)
    type(file_buffer_t), intent(inout) :: fb
    character(len=*), intent(in) :: string
    character(len=*), intent(in) :: name
    type(URI), pointer :: baseURI
    logical, intent(in), optional :: pe

    integer :: i
    type(xml_source_t), pointer :: temp(:)
    logical :: pe_

    if (present(pe)) then
      pe_ = pe
    else
      pe_ = .false.
    endif

    if (.not.associated(fb%f)) allocate(fb%f(0))

    temp => fb%f
    allocate(fb%f(size(temp)+1))
    do i = 1, size(temp)
      fb%f(i+1)%lun = temp(i)%lun
      fb%f(i+1)%xml_version = temp(i)%xml_version
      fb%f(i+1)%encoding => temp(i)%encoding
      fb%f(i+1)%filename => temp(i)%filename
      fb%f(i+1)%line = temp(i)%line
      fb%f(i+1)%col = temp(i)%col
      fb%f(i+1)%startChar = temp(i)%startChar
      fb%f(i+1)%next_chars => temp(i)%next_chars
      fb%f(i+1)%input_string => temp(i)%input_string
      fb%f(i+1)%baseURI => temp(i)%baseURI
      fb%f(i+1)%pe = temp(i)%pe
    enddo
    deallocate(temp)

    allocate(fb%f(1)%input_string)
    fb%f(1)%filename => vs_str_alloc(name)
    fb%f(1)%input_string%s => vs_str_alloc(string)
    if (pe_) then
      fb%f(1)%next_chars => vs_str_alloc(" ")
    else
      fb%f(1)%next_chars => vs_str_alloc("")
    endif
    fb%f(1)%pe = pe_
    if (associated(baseURI)) then
      fb%f(1)%baseURI => copyURI(baseURI)
    else
      fb%f(1)%baseURI => copyURI(fb%f(2)%baseURI)
    endif

  end subroutine open_new_string 

  subroutine pop_buffer_stack(fb)
    type(file_buffer_t), intent(inout) :: fb

    integer :: i
    type(xml_source_t), pointer :: temp(:)
    
    call close_actual_file(fb%f(1))

    temp => fb%f
    allocate(fb%f(size(temp)-1))
    do i = 1, size(temp)-1
      fb%f(i)%lun = temp(i+1)%lun
      fb%f(i)%xml_version = temp(i+1)%xml_version
      fb%f(i)%encoding => temp(i+1)%encoding
      fb%f(i)%filename => temp(i+1)%filename
      fb%f(i)%line = temp(i+1)%line
      fb%f(i)%col = temp(i+1)%col
      fb%f(i)%startChar = temp(i+1)%startChar
      fb%f(i)%next_chars => temp(i+1)%next_chars
      fb%f(i)%input_string => temp(i+1)%input_string
      fb%f(i)%baseURI => temp(i+1)%baseURI
      fb%f(i)%pe = temp(i+1)%pe
    enddo
    deallocate(temp)

  end subroutine pop_buffer_stack


  subroutine push_chars(fb, s)
    type(file_buffer_t), intent(inout) :: fb
    character(len=*), intent(in) :: s

    call push_file_chars(fb%f(1), s)

  end subroutine push_chars

  function get_character(fb, eof, es) result(string)
    type(file_buffer_t), intent(inout) :: fb
    logical, intent(out) :: eof
    type(error_stack), intent(inout) :: es
    character(len=1) :: string

    type(xml_source_t), pointer :: f
    character, pointer :: temp(:)

    f => fb%f(1)

    if (size(f%next_chars)>0) then
      eof = .false.
      string = f%next_chars(1)
      if (size(f%next_chars)>1) then
        temp => vs_str_alloc(str_vs(f%next_chars(2:)))
      else
        temp => vs_str_alloc("")
      endif
      deallocate(f%next_chars)
      f%next_chars => temp
    else
      string = get_char_from_file(f, fb%xml_version, eof, es)
    endif

  end function get_character

  function get_all_characters(fb, es) result(s)
    type(file_buffer_t), intent(inout) :: fb
    type(error_stack), intent(inout) :: es
    character, pointer :: s(:)

    logical :: eof
    character :: c
    character, pointer :: temp(:)

    eof = .false.
    s => vs_str_alloc("")
    do while (.not.eof)
      c = get_character(fb, eof, es)
      if (eof.or.in_error(es)) return
      temp => vs_str_alloc(str_vs(s)//c)
      deallocate(s)
      s => temp
    enddo
  end function get_all_characters

  function line(fb) result(n)
    type(file_buffer_t), intent(in) :: fb
    integer :: n

    n = fb%f(1)%line
  end function line

  function column(fb) result(n)
    type(file_buffer_t), intent(in) :: fb
    integer :: n

    n = fb%f(1)%col
  end function column

  subroutine parse_xml_declaration(fb, xv, enc, sa, es)
    type(file_buffer_t), intent(inout) :: fb
    integer, intent(out) :: xv
    character, pointer :: enc(:)
    logical, intent(out) :: sa
    type(error_stack), intent(inout) :: es

    logical :: eof

    call parse_declaration(fb%f(1), eof, es, sa)
    if (eof.or.in_error(es)) then
      call add_error(es, "Error parsing XML declaration")
    else
      fb%xml_version = fb%f(1)%xml_version
      xv = fb%xml_version
      enc => vs_vs_alloc(fb%f(1)%encoding)
    endif
  end subroutine parse_xml_declaration

  subroutine parse_text_declaration(fb, es)
    type(file_buffer_t), intent(inout) :: fb
    type(error_stack), intent(inout) :: es

    logical :: eof
    integer :: xv

    xv = fb%f(size(fb%f))%xml_version

    call parse_declaration(fb%f(1), eof, es)
    if (in_error(es)) then
      call add_error(es, "Error parsing text declaration")
      return
    elseif (xv==XML1_0.and.fb%f(1)%xml_version/=XML1_0) then
      call add_error(es, "XML 1.0 document cannot reference entities with higher version numbers")
      return
    endif

  end subroutine parse_text_declaration


  function reading_main_file(fb) result(p)
    type(file_buffer_t), intent(in) :: fb
    logical :: p

    p = (size(fb%f)==1)
  end function reading_main_file

 function reading_first_entity(fb) result(p)
    type(file_buffer_t), intent(in) :: fb
    logical :: p

    p = (size(fb%f)==2)
  end function reading_first_entity
#endif

end module m_sax_reader
