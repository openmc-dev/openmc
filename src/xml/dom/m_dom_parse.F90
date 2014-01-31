module m_dom_parse

  use fox_m_fsys_array_str, only: str_vs, vs_str_alloc
  use fox_m_utils_uri, only: URI, parseURI, rebaseURI, expressURI, destroyURI
  use m_common_attrs, only: hasKey, getValue, getIndex, getIsId, getBase,      &
    add_item_to_dict
  use m_common_entities, only: entity_t, size, getEntityByIndex
  use m_common_error, only: FoX_error, in_error
  use m_common_struct, only: xml_doc_state
  use FoX_common, only: dictionary_t, getLength
  use FoX_common, only: getQName, getValue, getURI, isSpecified
  use m_sax_parser, only: sax_parse
  use FoX_sax, only: xml_t
  use FoX_sax, only: open_xml_file, open_xml_string, close_xml_t

  ! Public interfaces
  use m_dom_dom, only: DOMConfiguration, Node, NamedNodeMap,                   &
    TEXT_NODE,                                                                 &
    getAttributes, getData, getDocType, getEntities, getImplementation,        &
    getLastChild, getNodeType,         &
    getNotations, getParameter, getParentNode, getXmlVersion,                  &
    setAttribute, setAttributeNS, setData, setValue,                           &
    appendChild, createAttribute, createAttributeNS, createCdataSection,       &
    createComment, createDocumentType, createElement, createElementNS,         &
    createEntityReference, createProcessingInstruction, createTextNode,        &
    getNamedItem, setAttributeNode, setAttributeNodeNS, setNamedItem,          &
    getFoX_checks

  ! Private interfaces
  use m_dom_dom, only: copyDOMConfig, createEmptyDocument, setDocumentElement, &
    createEmptyEntityReference, createEntity, createNotation,    &
    getReadOnly, getStringValue, getXds, destroy, destroyAllNodesRecursively,  &
    namespaceFixup, setDocType, setDomConfig, setGCstate, setIllFormed,        &
    setIsElementContentWhitespace, setIsId, setReadOnlyMap, setReadonlyNode,   &
    setSpecified, setXds, setStringValue
    
  use m_dom_error, only: DOMException, inException, throw_exception,           &
    getExceptionCode, PARSE_ERR

  implicit none
  private

  public :: parsefile
  public :: parsestring

  type(xml_t), target, save :: fxml

  type(Node), pointer, save  :: mainDoc => null()
  type(Node), pointer, save  :: current => null()

  type(DOMConfiguration), pointer :: domConfig
  
  logical :: cdata
  character, pointer :: error(:) => null()
  character, pointer :: inEntity(:) => null()

contains

  subroutine startElement_handler(nsURI, localname, name, attrs)
    character(len=*),   intent(in) :: nsURI
    character(len=*),   intent(in) :: localname
    character(len=*),   intent(in) :: name
    type(dictionary_t), intent(in) :: attrs
   
    type(URI), pointer :: URIref, URIbase, newURI
    type(Node), pointer :: el, attr, dummy
    character, pointer :: baseURI(:)
    integer :: i

    if (getParameter(domConfig, "namespaces")) then
      el => createElementNS(mainDoc, nsURI, name)
    else
      el => createElement(mainDoc, name)
    endif

    if (getBase(attrs)/="") then
      i = getIndex(attrs, "xml:base")
      if (i>0) then
        URIbase => parseURI(getBase(attrs))
        URIref => parseURI(getValue(attrs, i))
        newURI => rebaseURI(URIbase, URIref)
        call destroyURI(URIbase)
        call destroyURI(URIref)
        baseURI => vs_str_alloc(expressURI(newURI))
        call destroyURI(newURI)
      else
        baseURI => vs_str_alloc(getBase(attrs))
      endif
      if (getParameter(domConfig, "namespaces")) then
        attr => createAttributeNS(mainDoc, &
          "http://www.w3.org/XML/1998/namespace", "xml:base")
      else
        attr => createAttribute(mainDoc, "xml:base")
      endif
      call setValue(attr, str_vs(baseURI))
      deallocate(baseURI)
      if (i>0) then
        call setSpecified(attr, isSpecified(attrs, i))
        call setIsId(attr, getIsId(attrs, i))
      endif
      if (getParameter(domConfig, "namespaces")) then
        dummy => setAttributeNodeNS(el, attr)
      else
        dummy => setAttributeNode(el, attr)
      endif
    endif

    do i = 1, getLength(attrs)
      if (getQName(attrs, i)=="xml:base") cycle
      if (getParameter(domConfig, "namespaces")) then
        attr => createAttributeNS(mainDoc, getURI(attrs, i), getQName(attrs, i))
      else
        attr => createAttribute(mainDoc, getQName(attrs, i))
      endif
      call setValue(attr, getValue(attrs, i))
      call setSpecified(attr, isSpecified(attrs, i))
      call setIsId(attr, getIsId(attrs, i))
      if (getParameter(domConfig, "namespaces")) then
        dummy => setAttributeNodeNS(el, attr)
      else
        dummy => setAttributeNode(el, attr)
      endif
      if (associated(inEntity)) call setReadOnlyNode(attr, .true., .true.)
    enddo

    if (associated(current, mainDoc)) then
      current => appendChild(current,el)
      call setDocumentElement(mainDoc, current)
    else
      current => appendChild(current,el)
    endif
    if (getParameter(domConfig, "namespaces")) &
       call namespaceFixup(current, .false.)

    if (associated(inEntity)) &
      call setReadOnlyMap(getAttributes(current), .true.)

    cdata = .false.

  end subroutine startElement_handler

  subroutine endElement_handler(URI, localName, name)
    character(len=*), intent(in)     :: URI
    character(len=*), intent(in)     :: localname
    character(len=*), intent(in)     :: name

    if (associated(inEntity)) call setReadOnlyNode(current, .true., .false.)

    current => getParentNode(current)
  end subroutine endElement_handler

  ! FIXME to pick up entity references within attribute values, we need
  ! separate just_the_element, start_attribute, attribute_text etc. calls.

  subroutine characters_handler(chunk)
    character(len=*), intent(in) :: chunk

    type(Node), pointer :: temp
    logical :: readonly

    temp => getLastChild(current)
    if (associated(temp)) then
      if (.not.cdata.and.getNodeType(temp)==TEXT_NODE) then
        readonly = getReadOnly(temp) ! Reset readonly status quickly
        call setReadOnlyNode(temp, .false., .false.)
        call setData(temp, getData(temp)//chunk)
        call setReadOnlyNode(temp, readonly, .false.)
        return
      endif
    endif
    if (cdata) then
      temp => createCdataSection(mainDoc, chunk)
      temp => appendChild(current, temp)
    else
      temp => createTextNode(mainDoc, chunk)
      temp => appendChild(current, temp)
    endif

    if (associated(inEntity)) call setReadOnlyNode(temp, .true., .false.)

  end subroutine characters_handler

  subroutine ignorableWhitespace_handler(chunk)
    character(len=*), intent(in) :: chunk

    type(Node), pointer :: temp
    logical :: readonly

    if (getParameter(domConfig, "element-content-whitespace")) then
      temp => getLastChild(current)
      if (associated(temp)) then
        if (getNodeType(temp)==TEXT_NODE) then
          readonly = getReadOnly(temp) ! Reset readonly status quickly
          call setReadOnlyNode(temp, .false., .false.)
          call setData(temp, getData(temp)//chunk)
          call setReadOnlyNode(temp, readonly, .false.)
          call setIsElementContentWhitespace(temp, .true.)
          return
        endif
      endif
      temp => createTextNode(mainDoc, chunk)
      temp => appendChild(current, temp)
      call setIsElementContentWhitespace(temp, .true.)
      if (associated(inEntity)) call setReadOnlyNode(temp, .true., .false.)
    endif

  end subroutine ignorableWhitespace_handler

  subroutine comment_handler(comment)
    character(len=*), intent(in) :: comment

    type(Node), pointer :: temp

    if (getParameter(domConfig, "comments")) then
      temp => appendChild(current, createComment(mainDoc, comment))
      if (associated(inEntity)) call setReadOnlyNode(temp, .true., .false.)
    endif

  end subroutine comment_handler

  subroutine processingInstruction_handler(target, data)
    character(len=*), intent(in) :: target
    character(len=*), intent(in) :: data

    type(Node), pointer :: temp

    temp => appendChild(current, &
      createProcessingInstruction(mainDoc, target, data))

    if (associated(inEntity)) call setReadOnlyNode(temp, .true., .false.)
  end subroutine processingInstruction_handler

  subroutine startDocument_handler
    mainDoc => createEmptyDocument()
    current => mainDoc
    call setGCstate(mainDoc, .false.)
    call setDomConfig(mainDoc, domConfig)
  end subroutine startDocument_handler

  subroutine endDocument_Handler
    call setGCstate(mainDoc, .true.)
  end subroutine endDocument_Handler

  subroutine startDTD_handler(name, publicId, systemId)
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: publicId
    character(len=*), intent(in) :: systemId

    type(Node), pointer :: np

    np => createDocumentType(getImplementation(mainDoc), name, publicId=publicId, systemId=systemId)
    np => appendChild(mainDoc, np)
    call setDocType(mainDoc, np)

  end subroutine startDTD_handler

  subroutine endDTD_handler

    type(Node), pointer :: np, oldcurrent
    type(NamedNodeMap), pointer :: entities
    type(xml_t) :: subsax
    type(xml_doc_state), pointer :: xds
    type(entity_t), pointer :: ent
    integer :: i, ios
    logical :: ok

    entities => getEntities(getDocType(mainDoc))
    xds => getXds(mainDoc)

    do i = 1, size(xds%entityList)
      ent => getEntityByIndex(xds%entityList, i)
      np => getNamedItem(entities, str_vs(ent%name))

      ok = .false.
      if (ent%external) then
        if (size(ent%notation)==0) then
          call open_xml_file(subsax, expressURI(ent%baseURI), iostat=ios)
          if (ios/=0) then
            call setIllFormed(np, .true.)
          else
            ok = .true.
          endif
        endif
      else
        call open_xml_string(subsax, getStringValue(np))
        ok = .true.
      endif
      if (ok) then
        oldcurrent => current
        current => np
        ! Run the parser over value
        ! We do this with all entities already declared.
        call sax_parse(subsax%fx, subsax%fb,                           &
          startElement_handler=startElement_handler,                   &
          endElement_handler=endElement_handler,                       &
          characters_handler=characters_handler,                       &
          startCdata_handler=startCdata_handler,                       &
          endCdata_handler=endCdata_handler,                           &
          comment_handler=comment_handler,                             &
          processingInstruction_handler=processingInstruction_handler, &
          fatalError_handler=entityErrorHandler,                       &
          startInCharData = .true.,                                    &
          externalEntity = ent%external,                               &
          xmlVersion = getXmlVersion(mainDoc),                         &
          namespaces=getParameter(domConfig, "namespaces"),            &
          initial_entities = xds%entityList)
        call close_xml_t(subsax)

        current => oldcurrent
      endif
    enddo

    if (associated(getDocType(mainDoc))) then
      call setReadonlyMap(getEntities(getDocType(mainDoc)), .true.)
      call setReadonlyMap(getNotations(getDocType(mainDoc)), .true.)
    endif

  end subroutine endDTD_handler

  subroutine FoX_endDTD_handler(state)
    type(xml_doc_state), pointer :: state

    call setXds(mainDoc, state)

  end subroutine FoX_endDTD_handler

  subroutine notationDecl_handler(name, publicId, systemId)
    character(len=*), intent(in) :: name
    character(len=*), intent(in) ::  publicId
    character(len=*), intent(in) :: systemId
    
    type(Node), pointer :: np

    np => createNotation(mainDoc, name, publicId=publicId, systemId=systemId)
    np => setNamedItem(getNotations(getDocType(mainDoc)), np)
    ! The SAX parser will never give us duplicate entities,
    ! so there is no need to check

  end subroutine notationDecl_handler

  subroutine startCdata_handler()
    if (getParameter(domConfig, "cdata-sections")) cdata = .true.
  end subroutine startCdata_handler
  subroutine endCdata_handler()
    cdata = .false.
  end subroutine endCdata_handler

  subroutine internalEntityDecl_handler(name, value)
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: value

    type(Node), pointer :: np
    
    if (name(1:1)=="%") return
    ! Do nothing with parameter entities

    ! We only note that these exist here.
    ! A second parsing stage is triggered at the end
    ! of the DTD, in order to resolve entity references (which
    ! need not be declared in order)

    np => createEntity(mainDoc, name, "", "", "")
    call setStringValue(np, value)
    np => setNamedItem(getEntities(getDocType(mainDoc)), np)

  end subroutine internalEntityDecl_handler

  subroutine normalErrorHandler(msg)
    character(len=*), intent(in) :: msg
    ! This is called if the main parsing routine fails
    error => vs_str_alloc(msg)
  end subroutine normalErrorHandler

  subroutine entityErrorHandler(msg)
    character(len=*), intent(in) :: msg

    !This gets called if parsing of an entity failed. If so,
    !then we need to destroy all nodes which were being generated as
    !children of this entity, then mark the entity as ill-formed - but
    !otherwise carry on parsing the document, and only throw an error
    !if a reference is made to it.

    call destroyAllNodesRecursively(current, except=.true.)
    call setIllFormed(current, .true.)
  end subroutine entityErrorHandler

  subroutine externalEntityDecl_handler(name, publicId, systemId)
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: publicId
    character(len=*), intent(in) :: systemId
    type(Node), pointer :: np

    if (name(1:1)=="%") return
    ! Do nothing with parameter entities

    np => createEntity(mainDoc, name, &
      publicId=publicId, systemId=systemId, notationName="")
    np => setNamedItem(getEntities(getDocType(mainDoc)), np)

  end subroutine externalEntityDecl_handler

  subroutine unparsedEntityDecl_handler(name, publicId, systemId, notationName)
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: publicId
    character(len=*), intent(in) :: systemId
    character(len=*), intent(in) :: notationName
    type(Node), pointer :: np

    np => getNamedItem(getEntities(getDocType(mainDoc)), name)
    if (.not.associated(np)) then
      np => createEntity(mainDoc, name, publicId=publicId, systemId=systemId, notationName=notationName)
      np => setNamedItem(getEntities(getDocType(mainDoc)), np)
    endif

  end subroutine unparsedEntityDecl_handler

  subroutine startEntity_handler(name)
    character(len=*), intent(in) :: name

    if (name(1:1)=="%") return
    ! Do nothing with parameter entities

    if (getParameter(domConfig, "entities")) then
      if (.not.associated(inEntity)) then
        inEntity => vs_str_alloc(name)
      endif
      current => appendChild(current, createEmptyEntityReference(mainDoc, name))
    endif
  end subroutine startEntity_handler

  subroutine endEntity_handler(name)
    character(len=*), intent(in) :: name

    if (name(1:1)=="%") return
    ! Do nothing with parameter entities
    
    if (getParameter(domConfig, "entities")) then
      call setReadOnlyNode(current, .true., .false.)
      if (str_vs(inEntity)==name) deallocate(inEntity)
      current => getParentNode(current)
    endif

  end subroutine endEntity_handler

  subroutine skippedEntity_handler(name)
    character(len=*), intent(in) :: name
    
    type(Node), pointer :: temp

    if (name(1:1)=="%") return
    ! Do nothing with parameter entities

    temp => appendChild(current, createEntityReference(mainDoc, name))
    if (associated(inEntity)) call setReadonlyNode(temp, .true., .false.)
  end subroutine skippedEntity_handler


  subroutine runParser(fxml, configuration, ex)
    type(DOMException), intent(out), optional :: ex
    type(xml_t), intent(inout) :: fxml
    type(DOMConfiguration), pointer, optional :: configuration

    allocate(DOMConfig)
    if (present(configuration)) call copyDOMConfig(DOMConfig, configuration)

! We use internal sax_parse rather than public interface in order
! to use internal callbacks to get extra info.
    call sax_parse(fx=fxml%fx, fb=fxml%fb,&
      characters_handler=characters_handler,            &
      endDocument_handler=endDocument_handler,           &
      endElement_handler=endElement_handler,            &
      !endPrefixMapping_handler,      &
      ignorableWhitespace_handler=ignorableWhitespace_handler,   &
      processingInstruction_handler=processingInstruction_handler, &
      ! setDocumentLocator
      skippedEntity_handler=skippedEntity_handler,         &
      startDocument_handler=startDocument_handler,         & 
      startElement_handler=startElement_handler,          &
      !startPrefixMapping_handler,    &
      notationDecl_handler=notationDecl_handler,          &
      unparsedEntityDecl_handler=unparsedEntityDecl_handler, &
      !error_handler,            &
      fatalError_handler=normalErrorHandler,                 &
      !warning_handler,               &
      !attributeDecl_handler,         &
      !elementDecl_handler,           &
      externalEntityDecl_handler=externalEntityDecl_handler, &
      internalEntityDecl_handler=internalEntityDecl_handler,    &
      comment_handler=comment_handler,              &
      endCdata_handler=endCdata_handler,             &
      endDTD_handler=endDTD_handler,                &
      endEntity_handler=endEntity_handler,             &
      startCdata_handler=startCdata_handler,    &
      startDTD_handler=startDTD_handler,          &
      startEntity_handler=startEntity_handler, &
      FoX_endDTD_handler=FoX_endDTD_handler, &
      namespaces = getParameter(domConfig, "namespaces"),     &
      namespace_prefixes = .true., &
      validate = getParameter(domConfig, "validate"), & ! FIXME what about validate-if-present ...
      xmlns_uris = .true.)

    call close_xml_t(fxml)

    if (associated(error)) then
      if (associated(inEntity)) deallocate(inEntity)
      ! FIXME pass the value of the error through
      ! when we let exceptions do that
      deallocate(error)
      call destroy(mainDoc)
      if (getFoX_checks().or.PARSE_ERR<200) then
  call throw_exception(PARSE_ERR, "runParser", ex)
  if (present(ex)) then
    if (inException(ex)) then
       return
    endif
  endif
endif

    endif

  end subroutine runParser


  function parsefile(filename, configuration, iostat, ex) 
    type(DOMException), intent(out), optional :: ex
    character(len=*), intent(in) :: filename
    type(DOMConfiguration), pointer, optional :: configuration
    integer, intent(out), optional :: iostat
    type(Node), pointer :: parsefile

    type(DOMException) :: ex_
    integer :: iostat_

    call open_xml_file(fxml, filename, iostat_)
    if (present(iostat)) then
      iostat = iostat_
      if (iostat/=0) return
    elseif (in_error(fxml%fx%error_stack)) then
      call FoX_error(str_vs(fxml%fx%error_stack%stack(1)%msg))
    elseif (iostat_/=0) then
      call FoX_error("Cannot open file")
    endif

    if (present(ex)) then
      call runParser(fxml, configuration, ex)
    elseif (present(iostat)) then
      call runParser(fxml, configuration, ex_)
    else
      call runParser(fxml, configuration)
    endif

    if (present(iostat).and.inException(ex_)) then
      iostat = getExceptionCode(ex_)
    endif

    parsefile => mainDoc
    mainDoc => null()

  end function parsefile


  function parsestring(string, configuration, ex) 
    type(DOMException), intent(out), optional :: ex
    character(len=*), intent(in) :: string
    type(DOMConfiguration), pointer, optional :: configuration
    type(Node), pointer :: parsestring

    call open_xml_string(fxml, string)

    call runParser(fxml, configuration, ex)

    parsestring => mainDoc
    mainDoc => null()
    
  end function parsestring

end module m_dom_parse
