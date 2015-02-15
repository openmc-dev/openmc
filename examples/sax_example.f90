module m_handlers

  use FoX_common
  use FoX_sax

  implicit none
  private

  ! This module shows massive overkill with the quantity
  ! of callbacks specified. You will never need all of these
  ! unless you are implementing an alternative interface on
  ! top of SAX. (For example, the DOM interface uses all of 
  ! these.

  ! This is mainly for reference. Please see sax_example_2.f90
  ! for a more realistic example

  public :: characters_handler
  public :: endDocument_handler
  public :: endElement_handler
  public :: endPrefixMapping_handler
  public :: ignorableWhitespace_handler
  public :: processingInstruction_handler
  public :: skippedEntity_handler
  public :: startElement_handler
  public :: startDocument_handler
  public :: startPrefixMapping_handler
  public :: notationDecl_handler
  public :: unparsedEntityDecl_handler
  public :: error_handler
  public :: fatalError_handler
  public :: warning_handler
  public :: attributeDecl_handler
  public :: elementDecl_handler
  public :: externalEntityDecl_handler
  public :: internalEntityDecl_handler
  public :: comment_handler
  public :: endCdata_handler
  public :: endEntity_handler
  public :: endDTD_handler
  public :: startCdata_handler
  public :: startDTD_handler
  public :: startEntity_handler

contains

  subroutine characters_handler(chunk)
    character(len=*), intent(in) :: chunk

    write(*,'(a)') "PCDATA:"
    write(unit=*,fmt="(a)") chunk

  end subroutine characters_handler

  subroutine endDocument_handler
    write(*,'(a)') 'END DOCUMENT'
  end subroutine endDocument_handler

  subroutine endElement_handler(URI, localname, name)
    character(len=*), intent(in)     :: URI
    character(len=*), intent(in)     :: localname
    character(len=*), intent(in)     :: name

    write(unit=*,fmt="(4a)") ">>End Element: {", URI, "}", localname

  end subroutine endElement_handler
  
  subroutine endPrefixMapping_handler(prefix)
    character(len=*), intent(in) :: prefix

    write(unit=*,fmt='(a)') "END NAMESPACE MAPPING"
    write(unit=*,fmt='(2a)') "PREFIX:", prefix

  end subroutine endPrefixMapping_handler

  subroutine ignorableWhitespace_handler(chars)
    character(len=*), intent(in) :: chars
    write(*,'(2a)') str(len(chars))//" characters of ignorable whitespace found"
  end subroutine ignorableWhitespace_handler

  subroutine processingInstruction_handler(name, content)
    character(len=*), intent(in)   :: name
    character(len=*), intent(in)   :: content

    write(unit=*,fmt="(2a)") ">>Processing Instruction: ", name
    write(unit=*, fmt="(a)") content

  end subroutine processingInstruction_handler

  subroutine skippedEntity_handler(name)
    character(len=*), intent(in) :: name
    write(*,'(2a)') "Skipped entity, named "//name
  end subroutine skippedEntity_handler

  subroutine startDocument_handler
    write(*,'(a)') "Document started"
  end subroutine startDocument_handler
  
  subroutine startElement_handler(URI, localname, name,attributes)
    character(len=*), intent(in)   :: URI
    character(len=*), intent(in)   :: localname
    character(len=*), intent(in)   :: name
    type(dictionary_t), intent(in) :: attributes

    write(unit=*,fmt="(4a)") ">>Begin Element: {", URI, "}", localname
    write(unit=*,fmt="(a,i2,a)") "--- ", getLength(attributes), " attributes:"
    !call print_dict(attributes)
  end subroutine startElement_handler

  subroutine startPrefixMapping_handler(URI, prefix)
    character(len=*), intent(in) :: URI
    character(len=*), intent(in) :: prefix

    write(unit=*,fmt='(a)') "START NAMESPACE MAPPING"
    write(unit=*,fmt='(2a)') "PREFIX:", prefix
    write(unit=*,fmt='(2a)') "URI:", uri

  end subroutine startPrefixMapping_handler

  subroutine notationDecl_handler(name, publicId, systemId)
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: publicId
    character(len=*), intent(in) :: systemId

    write(*,'(2a)') "Notation declaration, "//name
    write(*,'(2a)') "Public ID: ", publicId
    write(*,'(2a)') "System ID: ", systemId
  end subroutine notationDecl_handler

  subroutine unparsedEntityDecl_handler(name, publicId, systemId, notation)
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: publicId
    character(len=*), intent(in) :: systemId
    character(len=*), intent(in) :: notation

    write(*,'(2a)') "Unparsed entity declaration, "//name
    write(*,'(2a)') "Public ID: ", publicId
    write(*,'(2a)') "System ID: ", systemId
    write(*,'(2a)') "Notation: ", notation
  end subroutine unparsedEntityDecl_handler

  subroutine error_handler(msg)
    character(len=*), intent(in) :: msg
    write(*,'(a)') "Error encountered and caught"
    write(*,'(a)') msg
  end subroutine error_handler

  subroutine fatalError_handler(msg)
    character(len=*), intent(in) :: msg
    write(*,'(a)') "Fatal error encountered and caught"
    write(*,'(a)') msg
  end subroutine fatalError_handler

  subroutine warning_handler(msg)
    character(len=*), intent(in) :: msg
    write(*,'(a)') "Warning encountered and caught"
    write(*,'(a)') msg
  end subroutine warning_handler

  subroutine attributeDecl_handler(eName, aName, type, mode, value)
    character(len=*), intent(in) :: eName
    character(len=*), intent(in) :: aName
    character(len=*), intent(in) :: type
    character(len=*), intent(in), optional :: mode
    character(len=*), intent(in), optional :: value
    write(*,'(a)') "Attribute declared"
    write(*,'(a)') eName
    write(*,'(a)') aName
    write(*,'(a)') type
    if (present(mode)) write(*,'(a)') mode
    if (present(value)) write(*,'(a)') value
  end subroutine attributeDecl_handler

  subroutine elementDecl_handler(name, model)
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: model
    write(*,'(a)') "Element declared"
    write(*,'(a)') name
    write(*,'(a)') model
  end subroutine elementDecl_handler

  subroutine externalEntityDecl_handler(name, publicId, systemId)
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: publicId
    character(len=*), intent(in) :: systemId
    write(*,'(a)') "External entity declared: "//name
    write(*,'(a)') "Public Id "//publicId
    write(*,'(a)') "System Id "//systemId
  end subroutine externalEntityDecl_handler
  
  subroutine internalEntityDecl_handler(name, value)
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: value
    write(*,'(a)') 'Internal entity declared: '//name
    write(*,'(a)') value
  end subroutine internalEntityDecl_handler

  subroutine comment_handler(comment)
    character(len=*), intent(in) :: comment
    write(unit=*,fmt="(a)") ">>Comment: "
    write(unit=*,fmt="(a)") comment
  end subroutine comment_handler

  subroutine endCdata_handler()
    write(*,'(a)') 'End CDATA section'
  end subroutine endCdata_handler

  subroutine endDTD_handler()
    write(*,'(a)') 'End DTD'
  end subroutine endDTD_handler
  
  subroutine endEntity_handler(name)
    character(len=*), intent(in) :: name
    write(*,'(a)') 'End of entity expansion: '//name
  end subroutine endEntity_handler

  subroutine startCdata_handler()
    write(*,'(a)') 'Start CDATA section'
  end subroutine startCdata_handler

  subroutine startDTD_handler(name, publicId, systemId)
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: publicId
    character(len=*), intent(in) :: systemId
    write(*,'(a)') 'DTD started: '//name
    write(*,'(a)') 'Public Id: '//publicId
    write(*,'(a)') 'System Id: '//systemId
  end subroutine startDTD_handler

  subroutine startEntity_handler(name)
    character(len=*), intent(in) :: name
    write(*,'(a)') 'Start of entity expansion: '//name
  end subroutine startEntity_handler

end module m_handlers












program sax_example
  !
  ! Example driver for a stand-alone parsing of an xml document
  !
  use FoX_sax

  use m_handlers      ! Defines begin_element, end_element, pcdata_chunk, etc

  implicit none

  integer :: iostat
  type(xml_t)  :: fxml

  call open_xml_file(fxml, "staffNS.xml", iostat=iostat)
  if (iostat /= 0) then
    write(*,*) "Cannot open file."
    stop
  endif

  call parse(fxml,&
    characters_handler,            &
    endDocument_handler,           &
    endElement_handler,            &
    endPrefixMapping_handler,      &
    ignorableWhitespace_handler,   &
    processingInstruction_handler, &
    ! setDocumentLocator
    skippedEntity_handler,         &
    startDocument_handler,         & 
    startElement_handler,          &
    startPrefixMapping_handler,    &
    notationDecl_handler,          &
    unparsedEntityDecl_handler,    &
    error_handler,                 &
    fatalError_handler,            &
    warning_handler,               &
    attributeDecl_handler,         &
    elementDecl_handler,           &
    externalEntityDecl_handler,    &
    internalEntityDecl_handler,    &
    comment_handler,               &
    endCdata_handler,              &
    endDTD_handler,                &
    endEntity_handler,             &
    startCdata_handler,            &
    startDTD_handler,              &
    startEntity_handler, validate=.true.)

  call close_xml_t(fxml)

end program sax_example
