module m_handlers

  use FoX_common
  use FoX_sax

  implicit none
  private

  ! A prototype of a specific language processor.

  ! It defines the routines that are called from xml_parser in response
  ! to particular events.

  ! In this particular example we just print the names of the elements
  ! and the content of the pcdata chunks, as well as any comments, XML
  ! and SGML declarations, etc.

  ! A module such as this could use "utility routines" to convert pcdata
  ! to numerical arrays, and to populate specific data structures.
  
  public :: start_document_handler, end_document_handler
  public :: begin_element_handler, end_element_handler
  public :: characters_handler
  public :: comment_handler, processing_instruction_handler
  public :: start_prefix_handler, end_prefix_handler

contains  !=============================================================

  subroutine start_document_handler()
    print*,'Document begun!!'
  end subroutine start_document_handler

  subroutine end_document_handler()
    print*,'Document finished!!'
  end subroutine end_document_handler

  subroutine begin_element_handler(URI, localname, name,attributes)
    character(len=*), intent(in)   :: URI
    character(len=*), intent(in)   :: localname
    character(len=*), intent(in)   :: name
    type(dictionary_t), intent(in) :: attributes

    write(unit=*,fmt="(4a)") ">>Begin Element: {", URI, "}", localname
    write(unit=*,fmt="(a,i2,a)") "--- ", len(attributes), " attributes:"
    call print_dict(attributes)
  end subroutine begin_element_handler

  !--------------------------------------------------
  subroutine end_element_handler(URI, localname, name)
    character(len=*), intent(in)     :: URI
    character(len=*), intent(in)     :: localname
    character(len=*), intent(in)     :: name

    write(unit=*,fmt="(4a)") ">>End Element: {", URI, "}", localname

  end subroutine end_element_handler

  !--------------------------------------------------
  subroutine characters_handler(chunk)
    character(len=*), intent(in) :: chunk

    write(*,'(a)') "PCDATA:"
    write(unit=*,fmt="(a)",advance="no") chunk

  end subroutine characters_handler

  !--------------------------------------------------
  subroutine comment_handler(comment)
    character(len=*), intent(in) :: comment

    write(unit=*,fmt="(a)") ">>Comment: "
    write(unit=*,fmt="(a)") comment

  end subroutine comment_handler


  subroutine processing_instruction_handler(name, content, attributes)
    character(len=*), intent(in)   :: name
    character(len=*), intent(in)   :: content
    type(dictionary_t), intent(in) :: attributes

    write(unit=*,fmt="(2a)") ">>Processing Instruction: ", name
    write(unit=*, fmt="(a)") content
    call print_dict(attributes)

  end subroutine processing_instruction_handler


  subroutine start_prefix_handler(URI, prefix)
    character(len=*), intent(in) :: URI
    character(len=*), intent(in) :: prefix

    write(unit=*,fmt='(a)') "START NAMESPACE MAPPING"
    write(unit=*,fmt='(2a)') "PREFIX:", prefix
    write(unit=*,fmt='(2a)') "URI:", uri

  end subroutine start_prefix_handler


  subroutine end_prefix_handler(prefix)
    character(len=*), intent(in) :: prefix

    write(unit=*,fmt='(a)') "END NAMESPACE MAPPING"
    write(unit=*,fmt='(2a)') "PREFIX:", prefix

  end subroutine end_prefix_handler

end module m_handlers












