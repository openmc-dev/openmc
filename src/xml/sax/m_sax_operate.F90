module m_sax_operate

#ifndef DUMMYLIB
  use m_common_error, only: FoX_error, in_error
  use FoX_common, only : str_vs

  use m_sax_reader, only: open_file, close_file
  use m_sax_parser, only: sax_parser_init, sax_parser_destroy, sax_parse
  use m_sax_types, only: ST_STOP
#endif
  use m_sax_types, only: xml_t

  implicit none
  private

  integer, parameter :: SAX_OPEN_ERROR = 1001

  public :: xml_t
  public :: open_xml_file
  public :: open_xml_string
  public :: close_xml_t
  public :: parse
  public :: stop_parser
  public :: SAX_OPEN_ERROR

contains

  subroutine open_xml_file(xt, file, iostat, lun)
    type(xml_t), intent(out) :: xt
    character(len=*), intent(in) :: file
    integer, intent(out), optional :: iostat
    integer, intent(in), optional :: lun
#ifdef DUMMYLIB
    if (present(iostat)) iostat = 0
#else
    integer :: i

    call open_file(xt%fb, file=trim(file), iostat=i, lun=lun, es=xt%fx%error_stack)
    if (present(iostat)) then
      if (in_error(xt%fx%error_stack)) i = SAX_OPEN_ERROR
      iostat = i
      if (i/=0) return
    else
      if (i/=0) &
        call FoX_error("Error opening file in open_xml_file")
      if (in_error(xt%fx%error_stack)) & 
        call FoX_error(str_vs(xt%fx%error_stack%stack(1)%msg))
    endif

    if (i==0) call sax_parser_init(xt%fx, xt%fb)
#endif
  end subroutine open_xml_file

  subroutine open_xml_string(xt, string)
    type(xml_t), intent(out) :: xt
    character(len=*), intent(in) :: string
#ifndef DUMMYLIB
    integer :: iostat

    call open_file(xt%fb, string=string, iostat=iostat, es=xt%fx%error_stack)
    call sax_parser_init(xt%fx, xt%fb)
#endif
  end subroutine open_xml_string

  subroutine close_xml_t(xt)
    type(xml_t), intent(inout) :: xt
#ifndef DUMMYLIB
    call close_file(xt%fb)
    call sax_parser_destroy(xt%fx)
#endif
  end subroutine close_xml_t


  subroutine parse(xt,             &
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
    startEntity_handler,           &
! Features / properties
    namespaces,                    &
    namespace_prefixes,            &
    validate,                      &
    xmlns_uris)

    type(xml_t), intent(inout) :: xt
    optional :: characters_handler
    optional :: endDocument_handler
    optional :: endElement_handler
    optional :: endPrefixMapping_handler
    optional :: ignorableWhitespace_handler
    optional :: processingInstruction_handler
    optional :: skippedEntity_handler
    optional :: startElement_handler
    optional :: startDocument_handler
    optional :: startPrefixMapping_handler
    optional :: notationDecl_handler
    optional :: unparsedEntityDecl_handler
    optional :: error_handler
    optional :: fatalError_handler
    optional :: warning_handler
    optional :: attributeDecl_handler
    optional :: elementDecl_handler
    optional :: externalEntityDecl_handler
    optional :: internalEntityDecl_handler
    optional :: comment_handler
    optional :: endCdata_handler
    optional :: endEntity_handler
    optional :: endDTD_handler
    optional :: startCdata_handler
    optional :: startDTD_handler
    optional :: startEntity_handler

    logical, intent(in), optional :: namespaces
    logical, intent(in), optional :: namespace_prefixes
    logical, intent(in), optional :: validate
    logical, intent(in), optional :: xmlns_uris

    interface

      subroutine characters_handler(chunk)
        character(len=*), intent(in) :: chunk
      end subroutine characters_handler

      subroutine endDocument_handler()     
      end subroutine endDocument_handler

      subroutine endElement_handler(namespaceURI, localName, name)
        character(len=*), intent(in)     :: namespaceURI
        character(len=*), intent(in)     :: localName
        character(len=*), intent(in)     :: name
      end subroutine endElement_handler
      
      subroutine endPrefixMapping_handler(prefix)
        character(len=*), intent(in) :: prefix
      end subroutine endPrefixMapping_handler

      subroutine ignorableWhitespace_handler(chars)
        character(len=*), intent(in) :: chars
      end subroutine ignorableWhitespace_handler

      subroutine processingInstruction_handler(name, content)
        character(len=*), intent(in)     :: name
        character(len=*), intent(in)     :: content
      end subroutine processingInstruction_handler

      subroutine skippedEntity_handler(name)
        character(len=*), intent(in) :: name
      end subroutine skippedEntity_handler

      subroutine startDocument_handler()   
      end subroutine startDocument_handler

      subroutine startElement_handler(namespaceURI, localName, name, attributes)
        use FoX_common
        character(len=*), intent(in)     :: namespaceUri
        character(len=*), intent(in)     :: localName
        character(len=*), intent(in)     :: name
        type(dictionary_t), intent(in)   :: attributes
      end subroutine startElement_handler

      subroutine startPrefixMapping_handler(namespaceURI, prefix)
        character(len=*), intent(in) :: namespaceURI
        character(len=*), intent(in) :: prefix
      end subroutine startPrefixMapping_handler

      subroutine notationDecl_handler(name, publicId, systemId)
        character(len=*), intent(in) :: name
        character(len=*), intent(in) :: publicId
        character(len=*), intent(in) :: systemId
      end subroutine notationDecl_handler

      subroutine unparsedEntityDecl_handler(name, publicId, systemId, notation)
        character(len=*), intent(in) :: name
        character(len=*), intent(in) :: publicId
        character(len=*), intent(in) :: systemId
        character(len=*), intent(in) :: notation
      end subroutine unparsedEntityDecl_handler

      subroutine error_handler(msg)
        character(len=*), intent(in)     :: msg
      end subroutine error_handler

      subroutine fatalError_handler(msg)
        character(len=*), intent(in)     :: msg
      end subroutine fatalError_handler
      
      subroutine warning_handler(msg)
        character(len=*), intent(in)     :: msg
      end subroutine warning_handler

      subroutine attributeDecl_handler(eName, aName, type, mode, value)
        character(len=*), intent(in) :: eName
        character(len=*), intent(in) :: aName
        character(len=*), intent(in) :: type
        character(len=*), intent(in), optional :: mode
        character(len=*), intent(in), optional :: value
      end subroutine attributeDecl_handler

      subroutine elementDecl_handler(name, model)
        character(len=*), intent(in) :: name
        character(len=*), intent(in) :: model
      end subroutine elementDecl_handler

      subroutine externalEntityDecl_handler(name, publicId, systemId)
        character(len=*), intent(in) :: name
        character(len=*), intent(in) :: publicId
        character(len=*), intent(in) :: systemId
      end subroutine externalEntityDecl_handler

      subroutine internalEntityDecl_handler(name, value)
        character(len=*), intent(in) :: name
        character(len=*), intent(in) :: value
      end subroutine internalEntityDecl_handler

      subroutine comment_handler(comment)
        character(len=*), intent(in) :: comment
      end subroutine comment_handler

      subroutine endCdata_handler()
      end subroutine endCdata_handler

      subroutine endDTD_handler()
      end subroutine endDTD_handler

      subroutine endEntity_handler(name)
        character(len=*), intent(in) :: name
      end subroutine endEntity_handler

      subroutine startCdata_handler()
      end subroutine startCdata_handler

      subroutine startDTD_handler(name, publicId, systemId)
        character(len=*), intent(in) :: name
        character(len=*), intent(in) :: publicId
        character(len=*), intent(in) :: systemId
      end subroutine startDTD_handler

      subroutine startEntity_handler(name)
        character(len=*), intent(in) :: name
      end subroutine startEntity_handler

    end interface
#ifndef DUMMYLIB
    ! FIXME check xt is initialized

    call sax_parse(xt%fx, xt%fb,     &
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
      startEntity_handler,           &
      namespaces=namespaces,                 &
      namespace_prefixes=namespace_prefixes, &
      validate=validate,                     &
      xmlns_uris=xmlns_uris)
#endif
  end subroutine parse

  subroutine stop_parser(xt)
    ! To be called from within a callback function;
    ! this will stop the parser.
    type(xml_t), intent(inout) :: xt
#ifndef DUMMYLIB
    xt%fx%state = ST_STOP
#endif
  end subroutine stop_parser

end module m_sax_operate
