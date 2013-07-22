module m_sax_parser

#ifndef DUMMYLIB
  use fox_m_fsys_array_str, only: str_vs, vs_str_alloc, vs_vs_alloc
  use fox_m_fsys_string_list, only: string_list, destroy_string_list, &
    tokenize_to_string_list, registered_string, init_string_list, &
    add_string, tokenize_and_add_strings, destroy
  use m_common_attrs, only: init_dict, destroy_dict, reset_dict, &
    add_item_to_dict, has_key, get_value, get_att_index_pointer, &
    getLength, setIsId, setBase
  use m_common_charset, only: XML1_0, XML1_1, XML_WHITESPACE
  use m_common_element, only: element_t, existing_element, add_element, &
    get_element, parse_dtd_element, parse_dtd_attlist, report_declarations, &
    declared_element, attribute_t, att_value_normalize, &
    get_attribute_declaration, &
    ATT_CDATA, ATT_ID, ATT_IDREF, ATT_IDREFS, ATT_ENTITY, ATT_ENTITIES, &
    ATT_NMTOKEN, ATT_NMTOKENS, ATT_NOTATION, ATT_ENUM, &
    ATT_REQUIRED, ATT_IMPLIED, ATT_DEFAULT, ATT_FIXED
  use m_common_elstack, only: push_elstack, pop_elstack, init_elstack, &
    destroy_elstack, is_empty, len, get_top_elstack, checkContentModel, &
    elementContent, emptyContent, checkContentModelToEnd
  use m_common_entities, only: existing_entity, init_entity_list, &
    destroy_entity_list, add_internal_entity, is_unparsed_entity, &
    expand_entity, expand_char_entity, pop_entity_list, size, &
    entity_t, getEntityByIndex, getEntityByName
  use m_common_entity_expand, only: expand_entity_value_alloc
  use m_common_error, only: FoX_error, add_error, &
    init_error_stack, destroy_error_stack, in_error
  use m_common_namecheck, only: checkName, checkPublicId, &
    checkCharacterEntityReference, likeCharacterEntityReference, &
    checkQName, checkNCName, checkPITarget, checkNmtoken, checkNmtokens, &
    checkRepCharEntityReference, checkNames, checkNCNames
  use m_common_namespaces, only: getnamespaceURI, invalidNS, &
    checkNamespaces, checkEndNamespaces, namespaceDictionary, &
    initNamespaceDictionary, destroyNamespaceDictionary
  use m_common_notations, only: init_notation_list, destroy_notation_list, &
    add_notation, notation_exists
  use m_common_struct, only: init_xml_doc_state, &
    destroy_xml_doc_state, register_internal_PE, register_external_PE, &
    register_internal_GE, register_external_GE

  use FoX_utils, only: URI, parseURI, rebaseURI, copyURI, destroyURI, &
    hasFragment, expressURI

  use m_sax_reader, only: file_buffer_t, pop_buffer_stack, open_new_string, &
    open_new_file, parse_xml_declaration, parse_text_declaration, &
    reading_main_file, reading_first_entity
  use m_sax_tokenizer, only: sax_tokenize, normalize_attribute_text, &
    expand_pe_text
  use m_sax_types ! everything, really

  implicit none
  private

  public :: getNSDict

  public :: sax_parser_init
  public :: sax_parser_destroy
  public :: sax_parse

contains

  function getNSDict(fx) result(ns)
    type(sax_parser_t), target :: fx
    type(namespaceDictionary), pointer :: ns

    ns => fx%nsDict
  end function getNSDict

  subroutine sax_parser_init(fx, fb)
    type(sax_parser_t), intent(out) :: fx
    type(file_buffer_t), intent(in) :: fb
#ifdef PGF90
    type(URI), pointer :: nullURI

    nullURI => null()
#endif

    allocate(fx%token(0))

    call init_error_stack(fx%error_stack)
    call init_elstack(fx%elstack)
    call init_dict(fx%attributes)

    call initNamespaceDictionary(fx%nsdict)
    call init_notation_list(fx%nlist)
    ! FIXME do we copy correctly from fx%nlist to fx%xds%nlist?
    allocate(fx%xds)
    call init_xml_doc_state(fx%xds)
    deallocate(fx%xds%inputEncoding)
    fx%xds%inputEncoding => vs_str_alloc("us-ascii")
    ! because it always is ...
    if (fb%f(1)%lun>0) then
      fx%xds%documentURI => vs_vs_alloc(fb%f(1)%filename)
    else
      fx%xds%documentURI => vs_str_alloc("")
    endif

    fx%xds%standalone = fb%standalone

    call init_entity_list(fx%forbidden_ge_list)
    call init_entity_list(fx%forbidden_pe_list)
    call init_entity_list(fx%predefined_e_list)

#ifdef PGF90
    call add_internal_entity(fx%predefined_e_list, 'amp', '&', nullURI, .false.)
    call add_internal_entity(fx%predefined_e_list, 'lt', '<', nullURI, .false.)
    call add_internal_entity(fx%predefined_e_list, 'gt', '>', nullURI, .false.)
    call add_internal_entity(fx%predefined_e_list, 'apos', "'", nullURI, .false.)
    call add_internal_entity(fx%predefined_e_list, 'quot', '"', nullURI, .false.)
#else
    call add_internal_entity(fx%predefined_e_list, 'amp', '&', null(), .false.)
    call add_internal_entity(fx%predefined_e_list, 'lt', '<', null(), .false.)
    call add_internal_entity(fx%predefined_e_list, 'gt', '>', null(), .false.)
    call add_internal_entity(fx%predefined_e_list, 'apos', "'", null(), .false.)
    call add_internal_entity(fx%predefined_e_list, 'quot', '"', null(), .false.)
#endif
  end subroutine sax_parser_init

  subroutine sax_parser_destroy(fx)
    type(sax_parser_t), intent(inout) :: fx

    fx%context = CTXT_NULL
    fx%state = ST_NULL

    if (associated(fx%token)) deallocate(fx%token)
    if (associated(fx%root_element)) deallocate(fx%root_element)

    call destroy_error_stack(fx%error_stack)
    call destroy_elstack(fx%elstack)
    call destroy_dict(fx%attributes)
    call destroyNamespaceDictionary(fx%nsdict)
    call destroy_notation_list(fx%nlist)
    if (.not.fx%xds_used) then
      call destroy_xml_doc_state(fx%xds)
      deallocate(fx%xds)
    endif

    call destroy_entity_list(fx%forbidden_ge_list)
    call destroy_entity_list(fx%forbidden_pe_list)
    call destroy_entity_list(fx%predefined_e_list)

    if (associated(fx%token)) deallocate(fx%token)
    if (associated(fx%content)) deallocate(fx%content)
    if (associated(fx%name)) deallocate(fx%name)
    if (associated(fx%attname)) deallocate(fx%attname)
    if (associated(fx%publicId)) deallocate(fx%publicId)
    if (associated(fx%systemId)) deallocate(fx%systemId)
    if (associated(fx%Ndata)) deallocate(fx%Ndata)

  end subroutine sax_parser_destroy

  recursive subroutine sax_parse(fx, fb, &
                                ! org.xml.sax
                                ! SAX ContentHandler
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
                                ! SAX DTDHandler
    notationDecl_handler,          &
    unparsedEntityDecl_handler,    &
                                ! SAX ErrorHandler
    error_handler,                 &
    fatalError_handler,            &
    warning_handler,               &
                                ! org.xml.sax.ext
                                ! SAX DeclHandler
    attributeDecl_handler,         &
    elementDecl_handler,           &
    externalEntityDecl_handler,    &
    internalEntityDecl_handler,    &
                                ! SAX LexicalHandler
    comment_handler,               &
    endCdata_handler,              &
    endDTD_handler,                &
    endEntity_handler,             &
    startCdata_handler,            &
    startDTD_handler,              &
    startEntity_handler,           &
    namespaces,                    &
    namespace_prefixes,            &
    xmlns_uris,                    &
    validate,                      &
    FoX_endDTD_handler,            &
    startInCharData,               &
    externalEntity,                &
    xmlVersion,                    &
    initial_entities)

    type(sax_parser_t), intent(inout) :: fx
    type(file_buffer_t), intent(inout) :: fb
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
    optional :: FoX_endDTD_handler
    optional :: startCdata_handler
    optional :: startDTD_handler
    optional :: startEntity_handler

    logical, intent(in), optional :: namespaces
    logical, intent(in), optional :: namespace_prefixes
    logical, intent(in), optional :: xmlns_uris

    logical, intent(in), optional :: validate
    logical, intent(in), optional :: startInCharData
    logical, intent(in), optional :: externalEntity
    character(len=*), intent(in), optional :: xmlVersion

    type(entity_list), optional :: initial_entities
#ifdef PGF90
    type(URI), pointer :: nullURI
#endif

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

      subroutine FoX_endDTD_handler(state)
        use m_common_struct, only: xml_doc_state
        type(xml_doc_state), pointer :: state
      end subroutine FoX_endDTD_handler

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

    logical :: validCheck, startInCharData_, processDTD, pe, nameOK, eof
    logical :: namespaces_, namespace_prefixes_, xmlns_uris_, externalEntity_
    integer :: i, iostat, temp_i, nextState, ignoreDepth, declSepValue
    character, pointer :: tempString(:)
    character :: dummy
    type(element_t), pointer :: elem
    type(attribute_t), pointer :: attDecl
    type(entity_t), pointer :: ent
    type(URI), pointer :: extSubsetURI, URIref, newURI
    integer, pointer :: wf_stack(:), temp_wf_stack(:), extEntStack(:)
    logical :: inExtSubset
    type(string_list) :: id_list, idref_list

#ifdef PGF90
    nullURI => null()
#endif
    tempString => null()
    elem => null()
    attDecl => null()

    if (present(namespaces)) then
      namespaces_ = namespaces
    else
      namespaces_ = .true.
    endif
    if (present(namespace_prefixes)) then
      namespace_prefixes_ = namespace_prefixes
    else
      namespace_prefixes_ = .false.
    endif
    if (present(xmlns_uris)) then
      xmlns_uris_ = xmlns_uris
    else
      xmlns_uris_ = .false.
    endif
    if (present(validate)) then
      validCheck = validate
    else
      validCheck = .false.
    endif
    if (present(startInCharData)) then
      startInCharData_ = startInCharData
    else
      startInCharData_ = .false.
    endif
    if (present(externalEntity)) then
      externalEntity_ = externalEntity
    else
      externalEntity_ = .false.
    endif
    if (present(initial_entities)) then
      do i = 1, size(initial_entities)
        ent => getEntityByIndex(initial_entities, i)
        if (ent%external) then
          call register_external_GE(fx%xds, &
            name=str_vs(ent%name), systemId=str_vs(ent%systemId), &
            publicId=str_vs(ent%publicId), &
            wfc=ent%wfc, baseURI=copyURI(ent%baseURI))
        else
          call register_internal_GE(fx%xds, &
            name=str_vs(ent%name), text=str_vs(ent%text), &
            wfc=ent%wfc, baseURI=copyURI(ent%baseURI))
        endif
      enddo
    endif

    allocate(wf_stack(1))
    wf_stack(1) = 0
    allocate(extEntStack(0))
    fx%inIntSubset = .false.
    extSubsetURI => null()
    inExtSubset = .false.
    declSepValue = 0
    processDTD = .true.
    iostat = 0

    if (startInCharData_) then
      fx%context = CTXT_IN_CONTENT
      fx%state = ST_CHAR_IN_CONTENT
      fx%well_formed = .true.
      if (externalEntity_) call parse_text_declaration(fb, fx%error_stack)
      if (in_error(fx%error_stack)) goto 100
      if (present(xmlVersion)) then
        if (xmlVersion=="1.0") then
          fx%xds%xml_version = XML1_0
        elseif (xmlVersion=="1.1") then
          fx%xds%xml_version = XML1_1
        endif
      endif
    elseif (reading_main_file(fb)) then
      fx%context = CTXT_BEFORE_DTD
      fx%state = ST_MISC
      if (present(startDocument_handler)) then
        call startDocument_handler()
        if (fx%state==ST_STOP) goto 100
      endif
      call parse_xml_declaration(fb, fx%xds%xml_version, fx%xds%encoding, fx%xds%standalone, fx%error_stack)
      if (in_error(fx%error_stack)) goto 100
      call init_string_list(id_list)
      call init_string_list(idref_list)
    endif

    do
      call sax_tokenize(fx, fb, eof)
      if (in_error(fx%error_stack)) then
        ! Any error, we want to quit sax_tokenizer
        call add_error(fx%error_stack, 'Error getting token')
        goto 100
      elseif (eof.and..not.reading_main_file(fb)) then
        if (inExtSubset.and.reading_first_entity(fb)) then
          if (wf_stack(1)>0) then
            call add_error(fx%error_stack, &
              "Unclosed conditional section or markup in external subset")
            goto 100
          elseif (fx%state_dtd/=ST_DTD_SUBSET) then
            call add_error(fx%error_stack, &
              "Markup not terminated in external subset")
            goto 100
          endif
          call endDTDchecks
          if (in_error(fx%error_stack)) goto 100
          if (fx%state==ST_STOP) goto 100
          inExtSubset = .false.
          fx%state = ST_MISC
          fx%context = CTXT_BEFORE_CONTENT
        elseif (fx%context==CTXT_IN_DTD) then
          if (validCheck) then
            if (wf_stack(1)/=0) then
              call add_error(fx%error_stack, &
                "Markup not terminated in parameter entity")
              goto 100
            endif
          endif
          if (declSepValue==size(wf_stack)) then
            if (wf_stack(1)/=0) then
              call add_error(fx%error_stack, &
                "Markup not terminated in parameter entity")
              goto 100
            else
              declSepValue = 0
            endif
          endif
          if (present(endEntity_handler)) then
            call endEntity_handler('%'//pop_entity_list(fx%forbidden_pe_list))
            if (fx%state==ST_STOP) goto 100
          else
            dummy = pop_entity_list(fx%forbidden_pe_list)
          endif
          if (fx%state_dtd==ST_DTD_ATTLIST_CONTENTS &
            .or.fx%state_dtd==ST_DTD_ELEMENT_CONTENTS) then
            ! stick the token back in contents ...
            fx%content => fx%token
            fx%token => vs_str_alloc("")
          endif
          if (reading_main_file(fb)) &
            fx%inIntSubset = .true.
        elseif (fx%context==CTXT_IN_CONTENT) then
          if (fx%state==ST_TAG_IN_CONTENT) fx%state = ST_CHAR_IN_CONTENT
          ! because CHAR_IN_CONTENT *always* leads to TAG_IN_CONTENT
          ! *except* when it is the end of an entity expansion
          if (present(endEntity_handler)) then
            call endEntity_handler(pop_entity_list(fx%forbidden_ge_list))
            if (fx%state==ST_STOP) goto 100
          else
            dummy = pop_entity_list(fx%forbidden_ge_list)
          endif
          if (fx%state/=ST_CHAR_IN_CONTENT.or.wf_stack(1)/=0) then
            call add_error(fx%error_stack, 'Ill-formed entity')
            goto 100
          endif
        endif
        temp_wf_stack => wf_stack
        allocate(wf_stack(size(temp_wf_stack)-1))
        wf_stack = temp_wf_stack(2:)
        ! If we are not doing validity checking, we might have 
        ! finished PE expansion with wf_stack(1) non-zero
        wf_stack(1) = wf_stack(1) + temp_wf_stack(1)
        deallocate(temp_wf_stack)
        temp_wf_stack => extEntStack
        allocate(extEntStack(size(temp_wf_stack)-1))
        extEntStack = temp_wf_stack(2:)
        deallocate(temp_wf_stack)
        call pop_buffer_stack(fb)
        cycle
      endif
      if (fx%tokenType==TOK_NULL) then
        call add_error(fx%error_stack, 'Internal error! No token found!')
        goto 100
      endif

      nextState = ST_NULL

      select case (fx%state)

      case (ST_MISC)
        !write(*,*) 'ST_MISC', str_vs(fx%token)
        select case (fx%tokenType)
        case (TOK_PI_TAG)
          wf_stack(1) = wf_stack(1) + 1
          nextState = ST_START_PI
        case (TOK_BANG_TAG)
          wf_stack(1) = wf_stack(1) + 1
          nextState = ST_BANG_TAG
        case (TOK_OPEN_TAG)
          wf_stack(1) = wf_stack(1) + 1
          nextState = ST_START_TAG
        end select


      case (ST_BANG_TAG)
        !write(*,*) 'ST_BANG_TAG'
        select case (fx%tokenType)
        case (TOK_OPEN_SB)
          nextState = ST_START_CDATA_DECLARATION
        case (TOK_OPEN_COMMENT)
          nextState = ST_START_COMMENT
        case (TOK_NAME)
          if (str_vs(fx%token)=='DOCTYPE') then
            fx%context = CTXT_IN_DTD
            nextState = ST_IN_DOCTYPE
          endif
        end select


      case (ST_START_PI)
        !write(*,*)'ST_START_PI'
        select case (fx%tokenType)
        case (TOK_NAME)
          if (namespaces_) then
            nameOk = checkNCName(str_vs(fx%token), fx%xds%xml_version)
          else
            nameOk = checkName(str_vs(fx%token), fx%xds%xml_version)
          endif
          if (nameOk) then
            if (str_vs(fx%token)=='xml') then
              call add_error(fx%error_stack, "XML declaration must be at start of document")
              goto 100
            elseif (checkPITarget(str_vs(fx%token), fx%xds%xml_version)) then
              nextState = ST_PI_CONTENTS
              fx%name => fx%token
              fx%token => null()
            else
              call add_error(fx%error_stack, "Invalid PI target name")
              goto 100
            endif
          endif
        end select

      case (ST_PI_CONTENTS)
        !write(*,*)'ST_PI_CONTENTS'
        if (validCheck) then
          if (emptyContent(fx%elstack)) then
            call add_error(fx%error_stack, "Content inside empty element")
            goto 100
          endif
        endif
        wf_stack(1) = wf_stack(1) - 1
          
        select case(fx%tokenType)
        case (TOK_CHAR)
          if (present(processingInstruction_handler)) then
            call processingInstruction_handler(str_vs(fx%name), str_vs(fx%token))
            if (fx%state==ST_STOP) goto 100
          endif
          deallocate(fx%name)
          nextState = ST_PI_END
        case (TOK_PI_END)
          if (present(processingInstruction_handler)) then
            call processingInstruction_handler(str_vs(fx%name), "")
            if (fx%state==ST_STOP) goto 100
          endif
          deallocate(fx%name)
          if (fx%context==CTXT_IN_CONTENT) then
            nextState = ST_CHAR_IN_CONTENT
          else
            nextState = ST_MISC
          endif
        end select

      case (ST_PI_END)
        !write(*,*)'ST_PI_END'
        select case(fx%tokenType)
        case (TOK_PI_END)
          if (fx%context==CTXT_IN_CONTENT) then
            nextState = ST_CHAR_IN_CONTENT
          else
            nextState = ST_MISC
          endif
        end select

      case (ST_START_COMMENT)
        !write(*,*)'ST_START_COMMENT'
        select case (fx%tokenType)
        case (TOK_CHAR)
          fx%name => fx%token
          nullify(fx%token)
          nextState = ST_COMMENT_END
        end select

      case (ST_COMMENT_END)
        !write(*,*)'ST_COMMENT_END'
        if (validCheck) then
          if (emptyContent(fx%elstack)) then
            call add_error(fx%error_stack, "Content inside empty element")
            goto 100
          endif
        endif
        wf_stack(1) = wf_stack(1) - 1

        select case (fx%tokenType)
        case (TOK_COMMENT_END)
          if (present(comment_handler)) then
            call comment_handler(str_vs(fx%name))
            if (fx%state==ST_STOP) goto 100
          endif
          deallocate(fx%name)
          if (fx%context==CTXT_IN_CONTENT) then
            nextState = ST_CHAR_IN_CONTENT
          else
            nextState = ST_MISC
          endif
        end select

      case (ST_START_TAG)
        !write(*,*)'ST_START_TAG', fx%context
        select case (fx%tokenType)
        case (TOK_NAME)
          if (fx%context==CTXT_BEFORE_DTD &
            .or. fx%context==CTXT_BEFORE_CONTENT &
            .or. fx%context==CTXT_IN_CONTENT) then
            if (namespaces_) then
              nameOk = checkQName(str_vs(fx%token), fx%xds%xml_version)
            else
              nameOk = checkName(str_vs(fx%token), fx%xds%xml_version)
            endif
            if (.not.nameOk) then
              call add_error(fx%error_stack, "Illegal element name")
              goto 100
            endif
            fx%name => fx%token
            nullify(fx%token)
            nextState = ST_IN_TAG
          elseif (fx%context == CTXT_AFTER_CONTENT) then
            call add_error(fx%error_stack, "Cannot open second root element")
            goto 100
          elseif (fx%context == CTXT_IN_DTD) then
            call add_error(fx%error_stack, "Cannot open root element before DTD is finished")
            goto 100
          endif
        end select

      case (ST_START_CDATA_DECLARATION)
        !write(*,*) "ST_START_CDATA_DECLARATION"
        select case (fx%tokenType)
        case (TOK_NAME)
          if (str_vs(fx%token)=="CDATA") then
            if (fx%context/=CTXT_IN_CONTENT) then
              call add_error(fx%error_stack, "CDATA section only allowed in text content.")
              goto 100
            else
              nextState = ST_FINISH_CDATA_DECLARATION
            endif
          else
            call add_error(fx%error_stack, "Unknown keyword found in marked section declaration.")
          endif
        end select

      case (ST_FINISH_CDATA_DECLARATION)
        !write(*,*) "ST_FINISH_CDATA_DECLARATION"
        select case (fx%tokenType)
        case (TOK_OPEN_SB)
          nextState = ST_CDATA_CONTENTS
        end select


      case (ST_CDATA_CONTENTS)
        !write(*,*)'ST_CDATA_CONTENTS'
        select case (fx%tokenType)
        case (TOK_CHAR)
          fx%name => fx%token
          nullify(fx%token)
          nextState = ST_CDATA_END
        end select

      case (ST_CDATA_END)
        !write(*,*)'ST_CDATA_END'
        if (validCheck) then
          if (emptyContent(fx%elstack).or.elementContent(fx%elstack)) then
            call add_error(fx%error_stack, "Content inside empty element")
            goto 100
          endif
        endif
        wf_stack(1) = wf_stack(1) - 1

        select case(fx%tokenType)
        case (TOK_SECTION_END)
          if (present(startCdata_handler)) then
            call startCdata_handler
            if (fx%state==ST_STOP) goto 100
          endif
          if (size(fx%name)>0) then
            if (present(characters_handler)) then
              call characters_handler(str_vs(fx%name))
              if (fx%state==ST_STOP) goto 100
            endif
          endif
          if (present(endCdata_handler)) then
            call endCdata_handler
            if (fx%state==ST_STOP) goto 100
          endif
          deallocate(fx%name)
          nextState = ST_CHAR_IN_CONTENT
        end select

      case (ST_IN_TAG)
        !write(*,*)'ST_IN_TAG'
        select case (fx%tokenType)
        case (TOK_END_TAG)
          if (fx%context /= CTXT_IN_CONTENT) then
            if (associated(fx%root_element)) then
              if (validCheck) then
                if (str_vs(fx%name)/=str_vs(fx%root_element)) then
                  call add_error(fx%error_stack, "Root element name does not match document name")
                  goto 100
                endif
              endif
              deallocate(fx%root_element)
            elseif (validCheck) then
              call add_error(fx%error_stack, "No DTD defined")
              goto 100
            else
              ! We havent had a DTD, so we havent handed xds 
              ! over to the DOM.
              if (present(FoX_endDTD_handler)) then
                fx%xds_used = .true.
                call FoX_endDTD_handler(fx%xds)
              endif
            endif
            fx%context = CTXT_IN_CONTENT
          endif
          call open_tag
          if (in_error(fx%error_stack)) goto 100
          if (fx%state==ST_STOP) goto 100
          deallocate(fx%name)
          nextState = ST_CHAR_IN_CONTENT

        case (TOK_END_TAG_CLOSE)
          if (fx%context==CTXT_IN_CONTENT) then
            nextState = ST_CHAR_IN_CONTENT
          else
            ! only a single element in this doc
            if (associated(fx%root_element)) then
              if (validCheck) then
                if (str_vs(fx%name)/=str_vs(fx%root_element)) then
                  call add_error(fx%error_stack, "Root element name does not match document name")
                  goto 100
                endif
              endif
              deallocate(fx%root_element)
            elseif (validCheck) then
              call add_error(fx%error_stack, "No DTD defined")
              goto 100
            else
              ! No DTD, so we havent handed over xds
              if (present(FoX_endDTD_handler)) then
                fx%xds_used = .true.
                call FoX_endDTD_handler(fx%xds)
              endif
            endif
          endif
          call open_tag
          if (in_error(fx%error_stack)) goto 100
          if (fx%state==ST_STOP) goto 100
          call close_tag
          if (in_error(fx%error_stack)) goto 100
          if (fx%state==ST_STOP) goto 100
          deallocate(fx%name)
          if (fx%context/=CTXT_IN_CONTENT) then
            fx%well_formed = .true.
            fx%context = CTXT_AFTER_CONTENT
            nextState = ST_MISC
          endif

        case (TOK_NAME)
          if (namespaces_) then
            nameOk = checkQName(str_vs(fx%token), fx%xds%xml_version)
          else
            nameOk = checkName(str_vs(fx%token), fx%xds%xml_version)
          endif
          if (.not.nameOk) then
            call add_error(fx%error_stack, "Illegal attribute name")
            goto 100
          endif
          !Have we already had this dictionary item?
          if (has_key(fx%attributes, str_vs(fx%token))) then
            call add_error(fx%error_stack, "Duplicate attribute name")
            goto 100
          endif
          fx%attname => fx%token
          nullify(fx%token)
          if (associated(elem)) then
            attDecl => get_attribute_declaration(elem, str_vs(fx%attname))
          else
            attDecl => null()
          endif
          nextState = ST_ATT_NAME
        end select

      case (ST_ATT_NAME)
        !write(*,*)'ST_ATT_NAME'
        select case (fx%tokenType)
        case (TOK_EQUALS)
          nextState = ST_ATT_EQUALS
        end select

      case (ST_ATT_EQUALS)
        !write(*,*)'ST_ATT_EQUALS'
        ! token is pre-processed attribute value.
        ! fx%name still contains attribute name
        select case (fx%tokenType)
        case (TOK_CHAR)
          !First, expand all entities:
          tempString => normalize_attribute_text(fx, fx%token)
          deallocate(fx%token)
          fx%token => tempString
          tempString => null()
          !If this attribute is not CDATA, we must process further;
          if (associated(attDecl)) then
            temp_i = attDecl%attType
          else
            temp_i = ATT_CDATA
          endif
          if (temp_i==ATT_CDATA) then
            call add_item_to_dict(fx%attributes, str_vs(fx%attname), &
              str_vs(fx%token), itype=ATT_CDATA, declared=associated(attDecl))
          else
            if (validCheck) then
              if (fx%xds%standalone.and..not.attDecl%internal &
                .and.(str_vs(fx%token)//"x"/=att_value_normalize(str_vs(fx%token))//"x")) then
                call add_error(fx%error_stack,  &
                  "Externally-declared attribute value normalization results in changed value "// &
                  "in standalone document")
                goto 100
              endif
            endif
            call add_item_to_dict(fx%attributes, str_vs(fx%attname), &
              att_value_normalize(str_vs(fx%token)), itype=temp_i, &
              declared=.true.)
          endif
          deallocate(fx%attname)
          nextState = ST_IN_TAG
        end select

      case (ST_CHAR_IN_CONTENT)
        !write(*,*)'ST_CHAR_IN_CONTENT'
        select case (fx%tokenType)
        case (TOK_CHAR)
          if (size(fx%token)>0) then
            if (validCheck) then
              if (elementContent(fx%elstack)) then
                if (verify(str_vs(fx%token), XML_WHITESPACE)==0) then
                  if (fx%xds%standalone.and..not.elem%internal) then
                    call add_error(fx%error_stack, &
                      "Externally-specified ignorable whitespace used in standalone document")
                    goto 100
                  endif
                  if (present(ignorableWhitespace_handler)) then
                    call ignorableWhitespace_handler(str_vs(fx%token))
                    if (fx%state==ST_STOP) goto 100
                  endif
                else
                  call add_error(fx%error_stack, "Forbidden content inside elementc: "//get_top_elstack(fx%elstack))
                  goto 100
                endif
              elseif (emptyContent(fx%elstack)) then
                call add_error(fx%error_stack, "Forbidden content inside element: "//get_top_elstack(fx%elstack))
                goto 100
              else
                if (present(characters_handler)) then
                  call characters_handler(str_vs(fx%token))
                  if (fx%state==ST_STOP) goto 100
                endif
              endif
            else
              if (present(characters_handler)) then
                call characters_handler(str_vs(fx%token))
                if (fx%state==ST_STOP) goto 100
              endif
            endif
          endif
          nextState = ST_TAG_IN_CONTENT
        end select

      case (ST_TAG_IN_CONTENT)
        !write(*,*) 'ST_TAG_IN_CONTENT'
        select case (fx%tokenType)
        case (TOK_ENTITY)
          nextState = ST_START_ENTITY
        case (TOK_PI_TAG)
          wf_stack(1) = wf_stack(1) + 1
          nextState = ST_START_PI
        case (TOK_BANG_TAG)
          wf_stack(1) = wf_stack(1) + 1
          nextState = ST_BANG_TAG
        case (TOK_CLOSE_TAG)
          nextState = ST_CLOSING_TAG
        case (TOK_OPEN_TAG)
          wf_stack(1) = wf_stack(1) + 1
          nextState = ST_START_TAG
        end select

      case (ST_START_ENTITY)
        !write(*,*) 'ST_START_ENTITY'
        select case (fx%tokenType)
        case (TOK_NAME)
          if (validCheck) then
            elem => get_element(fx%xds%element_list, get_top_elstack(fx%elstack))
            if (associated(elem)) then
              if (elem%empty) then
                call add_error(fx%error_stack, &
                  "Forbidden content inside element")
                goto 100
              endif
            else
              call add_error(fx%error_stack, &
                'Encountered reference to undeclared entity')
            endif
          endif
          ent => getEntityByName(fx%forbidden_ge_list, str_vs(fx%token))
          if (associated(ent)) then
            call add_error(fx%error_stack, 'Recursive entity reference')
            goto 100
          endif
          ent => getEntityByName(fx%predefined_e_list, str_vs(fx%token))
          if (associated(ent)) then
            if (present(startEntity_handler)) then
              call startEntity_handler(str_vs(fx%token))
              if (fx%state==ST_STOP) goto 100
            endif
            if (validCheck) then
              if (associated(elem)) then
                if (.not.elem%mixed.and..not.elem%any) then
                  call add_error(fx%error_stack, &
                    "Forbidden content inside element")
                  goto 100
                endif
              endif
            endif
            if (present(characters_handler)) then
              call characters_handler(expand_entity(fx%predefined_e_list, str_vs(fx%token)))
              if (fx%state==ST_STOP) goto 100
            endif
            if (present(endEntity_handler)) then
              call endEntity_handler(str_vs(fx%token))
              if (fx%state==ST_STOP) goto 100
            endif
          elseif (likeCharacterEntityReference(str_vs(fx%token))) then
            if (checkRepCharEntityReference(str_vs(fx%token), fx%xds%xml_version)) then
              if (validCheck) then
                if (associated(elem)) then
                  if (.not.elem%mixed.and..not.elem%any) then
                    call add_error(fx%error_stack, &
                      "Forbidden content inside element")
                    goto 100
                  endif
                endif
              endif
              if (present(characters_handler)) then
                call characters_handler(expand_char_entity(str_vs(fx%token)))
                if (fx%state==ST_STOP) goto 100
              endif
            elseif (checkCharacterEntityReference(str_vs(fx%token), fx%xds%xml_version)) then
              call add_error(fx%error_stack, "Unable to digest character entity reference in content, sorry.")
              goto 100
            else
              call add_error(fx%error_stack, "Illegal character reference")
              goto 100
            endif
          elseif (existing_entity(fx%xds%entityList, str_vs(fx%token))) then
            ent => getEntityByName(fx%xds%entityList, str_vs(fx%token))
            if (ent%wfc.and.fx%xds%standalone) then
              call add_error(fx%error_stack, &
                'Externally declared entity referenced in standalone document')
              goto 100
            elseif (str_vs(ent%notation)/="") then
              call add_error(fx%error_stack, &
                'Cannot reference unparsed entity in content')
              goto 100
            elseif (ent%external) then
              call open_new_file(fb, ent%baseURI, iostat)
              if (iostat/=0) then
                if (validCheck) then
                  call add_error(fx%error_stack, &
                    "Unable to retrieve external entity "//str_vs(fx%token))
                  goto 100
                endif
                if (present(skippedEntity_handler)) then
                  call skippedEntity_handler(str_vs(fx%token))
                  if (fx%state==ST_STOP) goto 100
                endif
              else
                if (present(startEntity_handler)) then
                  call startEntity_handler(str_vs(fx%token))
                  if (fx%state==ST_STOP) goto 100
                endif
#ifdef PGF90
                call add_internal_entity(fx%forbidden_ge_list, str_vs(fx%token), "", nullURI, .false.)
#else
                call add_internal_entity(fx%forbidden_ge_list, str_vs(fx%token), "", null(), .false.)
#endif
                temp_wf_stack => wf_stack
                allocate(wf_stack(size(temp_wf_stack)+1))
                wf_stack = (/0, temp_wf_stack/)
                deallocate(temp_wf_stack)
                temp_wf_stack => extEntStack
                allocate(extEntStack(size(temp_wf_stack)+1))
                extEntStack = (/len(fx%elstack), temp_wf_stack/)
                deallocate(temp_wf_stack)
                call parse_text_declaration(fb, fx%error_stack)
                if (in_error(fx%error_stack)) goto 100
              endif
            else
              if (validCheck.and.associated(elem)) then
                if (elem%empty) then
                  call add_error(fx%error_stack, &
                    "Forbidden content inside element")
                  goto 100
                endif
              endif
              if (present(startEntity_handler)) then
                call startEntity_handler(str_vs(fx%token))
                if (fx%state==ST_STOP) goto 100
              endif
#ifdef PGF90
              call add_internal_entity(fx%forbidden_ge_list, str_vs(fx%token), "", nullURI, .false.)
#else
              call add_internal_entity(fx%forbidden_ge_list, str_vs(fx%token), "", null(), .false.)
#endif
              call open_new_string(fb, expand_entity(fx%xds%entityList, str_vs(fx%token)), str_vs(fx%token), baseURI=ent%baseURI)
              temp_wf_stack => wf_stack
              allocate(wf_stack(size(temp_wf_stack)+1))
              wf_stack = (/0, temp_wf_stack/)
              deallocate(temp_wf_stack)
            endif
          else
            ! Unknown entity check standalone etc
            if (fx%skippedExternal.and..not.fx%xds%standalone) then
              if (present(skippedEntity_handler)) then
                call skippedEntity_handler(str_vs(fx%token))
                if (fx%state==ST_STOP) goto 100
              endif
            else
              call add_error(fx%error_stack, &
                'Encountered reference to undeclared entity')
            endif
          endif
          nextState = ST_CHAR_IN_CONTENT
        end select

      case (ST_CLOSING_TAG)
        !write(*,*)'ST_CLOSING_TAG'
        select case (fx%tokenType)
        case (TOK_NAME)
          if (checkName(str_vs(fx%token), fx%xds%xml_version)) then
            fx%name => fx%token
            nullify(fx%token)
            nextState = ST_IN_CLOSING_TAG
          else
            call add_error(fx%error_stack, "Closing tag: expecting a Name")
            goto 100
          end if
        end select

      case (ST_IN_CLOSING_TAG)
        !write(*,*)'ST_IN_CLOSING_TAG'
        select case (fx%tokenType)
        case (TOK_END_TAG)
          call close_tag
          if (in_error(fx%error_stack)) goto 100
          if (fx%state==ST_STOP) goto 100
          deallocate(fx%name)
          if (is_empty(fx%elstack)) then
            if (startInCharData_) then
              fx%well_formed = .true.
              nextState = ST_CHAR_IN_CONTENT
            else
              !we're done
              if (validCheck) then
                call checkIdRefs
                if (in_error(fx%error_stack)) goto 100
              endif
              fx%well_formed = .true.
              nextState = ST_MISC
              fx%context = CTXT_AFTER_CONTENT
            endif
          else
            nextState = ST_CHAR_IN_CONTENT
          endif
        end select

      case (ST_IN_DOCTYPE)
        !write(*,*)'ST_IN_DOCTYPE'
        select case (fx%tokenType)
        case (TOK_NAME)
          if (namespaces_) then
            nameOk = checkQName(str_vs(fx%token), fx%xds%xml_version)
          else
            nameOk = checkName(str_vs(fx%token), fx%xds%xml_version)
          endif
          if (.not.nameOk) then
            call add_error(fx%error_stack, "Invalid document name")
            goto 100
          endif
          fx%root_element => fx%token
          nullify(fx%token)
          nextState = ST_DOC_NAME
        end select

      case (ST_DOC_NAME)
        !write(*,*) 'ST_DOC_NAME ', str_vs(fx%token)
        select case (fx%tokenType)
        case (TOK_NAME)
          if (str_vs(fx%token)=='SYSTEM') then
            nextState = ST_DOC_SYSTEM
          elseif (str_vs(fx%token)=='PUBLIC') then
            nextState = ST_DOC_PUBLIC
          endif
        case (TOK_OPEN_SB)
          if (present(startDTD_handler)) then
            call startDTD_handler(str_vs(fx%root_element), "", "")
            if (fx%state==ST_STOP) goto 100
          endif
          wf_stack(1) = wf_stack(1) + 1
          nextState = ST_IN_SUBSET
          fx%inIntSubset = .true.
        case (TOK_END_TAG)
          if (present(startDTD_handler)) then
            call startDTD_handler(str_vs(fx%root_element), "", "")
            if (fx%state==ST_STOP) goto 100
          endif
          wf_stack(1) = wf_stack(1) - 1
          call endDTDchecks
          if (in_error(fx%error_stack)) goto 100
          if (fx%state==ST_STOP) goto 100
          fx%context = CTXT_BEFORE_CONTENT
          nextState = ST_MISC
        case default
          call add_error(fx%error_stack, "Unexpected token")
          goto 100
        end select

      case (ST_DOC_PUBLIC)
        !write(*,*) 'ST_DOC_PUBLIC'
        select case (fx%tokenType)
        case (TOK_CHAR)
          if (checkPublicId(str_vs(fx%token))) then
            fx%publicId => fx%token
            fx%token => null()
            nextState = ST_DOC_SYSTEM
          else
            call add_error(fx%error_stack, "Invalid document public id")
            goto 100
          endif
        end select

      case (ST_DOC_SYSTEM)
        !write(*,*) 'ST_DOC_SYSTEM'
        select case (fx%tokenType)
        case (TOK_CHAR)
          fx%systemId => fx%token
          fx%token => null()
          nextState = ST_DOC_DECL
        end select

      case (ST_DOC_DECL)
        !write(*,*) 'ST_DOC_DECL'
        select case (fx%tokenType)
        case (TOK_OPEN_SB)
          if (present(startDTD_handler)) then
            if (associated(fx%publicId)) then
              call startDTD_handler(str_vs(fx%root_element), &
                publicId=str_vs(fx%publicId), systemId=str_vs(fx%systemId))
            elseif (associated(fx%systemId)) then
              call startDTD_handler(str_vs(fx%root_element), &
                publicId="", systemId=str_vs(fx%systemId))
            else
              call startDTD_handler(str_vs(fx%root_element), "", "")
            endif
            if (fx%state==ST_STOP) goto 100
          endif
          if (associated(fx%systemId)) then
            extSubsetURI => parseURI(str_vs(fx%systemId))
            deallocate(fx%systemId)
          endif
          if (associated(fx%publicId)) deallocate(fx%publicId)
          fx%inIntSubset = .true.
          wf_stack(1) = wf_stack(1) + 1
          nextState = ST_IN_SUBSET
        case (TOK_END_TAG)
          if (present(startDTD_handler)) then
            if (associated(fx%publicId)) then
              call startDTD_handler(str_vs(fx%root_element), publicId=str_vs(fx%publicId), systemId=str_vs(fx%systemId))
              deallocate(fx%publicId)
            elseif (associated(fx%systemId)) then
              call startDTD_handler(str_vs(fx%root_element), publicId="", systemId=str_vs(fx%systemId))
            else
              call startDTD_handler(str_vs(fx%root_element), "", "")
            endif
            if (fx%state==ST_STOP) goto 100
          endif
          if (associated(fx%systemId)) then
            extSubsetURI => parseURI(str_vs(fx%systemId))
            if (.not.associated(extSubsetURI)) then
              call add_error(fx%error_stack, "Invalid URI specified for DTD SYSTEM")
              goto 100
            endif
            call open_new_file(fb, extSubsetURI, iostat)
            if (iostat==0) then
              fx%inIntSubset=.false.
              call parse_text_declaration(fb, fx%error_stack)
              if (in_error(fx%error_stack)) goto 100
              temp_wf_stack => wf_stack
              allocate(wf_stack(size(temp_wf_stack)+1))
              wf_stack = (/0, temp_wf_stack/)
              deallocate(temp_wf_stack)
              inExtSubset = .true.
              nextState = ST_IN_SUBSET
            else
              if (validCheck) then
                call add_error(fx%error_stack, &
                  "Unable to retrieve external subset "//str_vs(fx%systemId))
                goto 100
              endif
              call endDTDchecks
              if (in_error(fx%error_stack)) goto 100
              if (fx%state==ST_STOP) goto 100
              fx%context = CTXT_BEFORE_CONTENT
              nextState = ST_MISC
            endif
            call destroyURI(extSubsetURI)
          else
            call endDTDchecks
            if (in_error(fx%error_stack)) goto 100
            if (fx%state==ST_STOP) goto 100
            fx%context = CTXT_BEFORE_CONTENT
            nextState = ST_MISC
          endif
          if (associated(fx%systemId)) deallocate(fx%systemId)
          if (associated(fx%publicId)) deallocate(fx%publicId)
        case default
          call add_error(fx%error_stack, "Unexpected token in DTD")
          goto 100
        end select


      case (ST_CLOSE_DOCTYPE)
        !write(*,*) "ST_CLOSE_DOCTYPE"
        select case (fx%tokenType)
        case (TOK_END_TAG)
          if (wf_stack(1)>1) then
            call add_error(fx%error_stack, "Cannot end DTD while conditional section is still open")
            goto 100
          endif
          if (associated(extSubsetURI)) then
            call open_new_file(fb, extSubsetURI, iostat)
            call destroyURI(extSubsetURI)
            if (iostat==0) then
              fx%inIntSubset = .false.
              call parse_text_declaration(fb, fx%error_stack)
              if (in_error(fx%error_stack)) goto 100
              temp_wf_stack => wf_stack
              allocate(wf_stack(size(temp_wf_stack)+1))
              wf_stack = (/0, temp_wf_stack/)
              deallocate(temp_wf_stack)
              inExtSubset = .true.
              nextState = ST_IN_SUBSET
            else
              if (validCheck) then
                call add_error(fx%error_stack, &
                  "Unable to retrieve external subset")
                goto 100
              endif
              call endDTDchecks
              if (in_error(fx%error_stack)) goto 100
              if (fx%state==ST_STOP) goto 100
              wf_stack(1) = wf_stack(1) - 1
              nextState = ST_MISC
              fx%context = CTXT_BEFORE_CONTENT
            endif
          else
            call endDTDchecks
            if (in_error(fx%error_stack)) goto 100
            if (fx%state==ST_STOP) goto 100
            wf_stack(1) = wf_stack(1) - 1
            nextState = ST_MISC
            fx%context = CTXT_BEFORE_CONTENT
          endif
        end select

      case (ST_IN_SUBSET)
        select case(fx%tokenType)
        case (TOK_ENTITY)
          nextState = ST_START_PE
        case default
          call parseDTD
          if (in_error(fx%error_stack)) goto 100
          if (fx%state==ST_STOP) goto 100
          if (fx%state_dtd==ST_DTD_DONE) &
            fx%state_dtd = ST_DTD_SUBSET
        end select

      case (ST_START_PE)
        !write(*,*) 'ST_START_PE'
        select case (fx%tokenType)
        case (TOK_NAME)
          if (existing_entity(fx%forbidden_pe_list, str_vs(fx%token))) then
            call add_error(fx%error_stack, &
              'Recursive entity reference')
            goto 100
          endif
          ent => getEntityByName(fx%xds%PEList, str_vs(fx%token))
          if (associated(ent)) then
            if (ent%wfc.and.fx%xds%standalone) then
              call add_error(fx%error_stack, &
                "Externally declared entity used in standalone document")
              goto 100
            elseif (ent%external) then
              if (present(startEntity_handler)) then
                call startEntity_handler('%'//str_vs(fx%token))
                if (fx%state==ST_STOP) goto 100
              endif
#ifdef PGF90
              call add_internal_entity(fx%forbidden_pe_list, &
                str_vs(fx%token), "", nullURI, .false.)
#else
              call add_internal_entity(fx%forbidden_pe_list, &
                str_vs(fx%token), "", null(), .false.)
#endif
              call open_new_file(fb, ent%baseURI, iostat, pe=.true.)
              if (iostat/=0) then
                if (validCheck) then
                  call add_error(fx%error_stack, &
                    "Unable to retrieve external parameter entity "//str_vs(fx%token))
                  goto 100
                endif
                if (present(skippedEntity_handler)) then
                  call skippedEntity_handler('%'//str_vs(fx%token))
                  if (fx%state==ST_STOP) goto 100
                endif
                ! having skipped a PE, we must now not process
                ! declarations any further (unless we are declared standalone)
                ! (XML section 5.1)
                fx%skippedExternal = .true.
                processDTD = fx%xds%standalone
              else
                fx%inIntSubset = .false.
                if (present(startEntity_handler)) then
                  call startEntity_handler('%'//str_vs(fx%token))
                  if (fx%state==ST_STOP) goto 100
                endif
#ifdef PGF90
                call add_internal_entity(fx%forbidden_pe_list, &
                  str_vs(fx%token), "", nullURI, .false.)
#else
                call add_internal_entity(fx%forbidden_pe_list, &
                  str_vs(fx%token), "", null(), .false.)
#endif
                call parse_text_declaration(fb, fx%error_stack)
                if (in_error(fx%error_stack)) goto 100
                allocate(temp_wf_stack(size(wf_stack)+1))
                temp_wf_stack = (/0, wf_stack/)
                deallocate(wf_stack)
                wf_stack => temp_wf_stack
                if (fx%state_dtd==ST_DTD_SUBSET) &
                  declSepValue = size(wf_stack)
              endif
            else
              ! Expand the entity,
              if (present(startEntity_handler)) then
                call startEntity_handler('%'//str_vs(fx%token))
                if (fx%state==ST_STOP) goto 100
              endif
#ifdef PGF90
              call add_internal_entity(fx%forbidden_pe_list, &
                str_vs(fx%token), "", nullURI, .false.)
              call open_new_string(fb, &
                expand_entity(fx%xds%PEList, str_vs(fx%token)), str_vs(fx%token), baseURI=nullURI, pe=.true.)
#else
              call add_internal_entity(fx%forbidden_pe_list, &
                str_vs(fx%token), "", null(), .false.)
              call open_new_string(fb, &
                expand_entity(fx%xds%PEList, str_vs(fx%token)), str_vs(fx%token), baseURI=null(), pe=.true.)
#endif
              ! NB because we are just expanding a string here, anything
              ! evaluated as a result of this string is evaluated in the
              ! context of the currently open file, so has a baseURI of
              ! this file
              allocate(temp_wf_stack(size(wf_stack)+1))
              temp_wf_stack = (/0, wf_stack/)
              deallocate(wf_stack)
              wf_stack => temp_wf_stack
              if (fx%state_dtd==ST_DTD_SUBSET) &
                declSepValue = size(wf_stack)
            endif
            ! and do nothing else, carry on ...
          else
            ! Have we previously skipped an external entity?
            if (fx%skippedExternal.and..not.fx%xds%standalone) then
              if (processDTD) then
                if (present(skippedEntity_handler)) then
                  call skippedEntity_handler('%'//str_vs(fx%token))
                  if (fx%state==ST_STOP) goto 100
                endif
              endif
            else
              ! If not, 
              call add_error(fx%error_stack, &
                "Reference to undeclared parameter entity.")
              goto 100
            endif
          endif
          nextState = ST_IN_SUBSET

        end select

      end select

      if (nextState/=ST_NULL) then
        fx%state = nextState
      else
        call add_error(fx%error_stack, "Internal error in parser - no suitable token found")
        goto 100
      endif

    end do

100 if (associated(tempString)) deallocate(tempString)
    if (associated(extSubsetURI)) call destroyURI(extSubsetURI)
    call destroy_string_list(id_list)
    call destroy_string_list(idref_list)
    deallocate(wf_stack)
    if (associated(extEntStack)) deallocate(extEntStack)

    if (fx%state==ST_STOP) return
    if (.not.eof) then
      ! We have encountered an error before the end of a file
      if (.not.reading_main_file(fb)) then !we are parsing an entity
        if (inExtSubset) then
          call add_error(fx%error_stack, "Error encountered processing external subset.")
        else
          call add_error(fx%error_stack, "Error encountered processing entity.")
        endif
        call sax_error(fx, fatalError_handler)
      else
        call sax_error(fx, fatalError_handler)
      endif
    else
      ! EOF of main file
      if (startInChardata_) then
        if (fx%well_formed) then
          if (fx%state==ST_CHAR_IN_CONTENT.and.associated(fx%token)) then
            if (size(fx%token)>0.and.present(characters_handler)) &
              call characters_handler(str_vs(fx%token))
          endif
        else
          if (present(fatalError_handler)) &
            call fatalError_handler("Ill-formed XML fragment")
        endif
      elseif (fx%well_formed.and.fx%state==ST_MISC) then
        if (present(endDocument_handler)) &
          call endDocument_handler()
      else
        call add_error(fx%error_stack, "File is not well-formed")
        call sax_error(fx, fatalError_handler)
      endif
    endif

  contains

    subroutine parseDTD

      integer :: nextDTDState
#ifdef PGF90
      type(element_t), pointer :: nullElement

      nullElement => null()
#endif

      nextDTDState = ST_DTD_NULL

      select case (fx%state_dtd)

      case (ST_DTD_SUBSET)
        !write(*,*) "ST_DTD_SUBSET"
        select case (fx%tokenType)
        case (TOK_CLOSE_SB)
          if (.not.reading_main_file(fb)) then
            call add_error(fx%error_stack, "Cannot close DOCTYPE in external entity")
            return
          endif
          wf_stack(1) = wf_stack(1) - 1
          fx%inIntSubset = .false.
          nextState = ST_CLOSE_DOCTYPE
          nextDTDState = ST_DTD_DONE
        case (TOK_SECTION_END)
          if (wf_stack(1)==0) then
            call add_error(fx%error_stack, "Unbalanced conditional section in parameter entity")
            return
          endif
          wf_stack(1) = wf_stack(1) - 2
          nextDTDState = ST_DTD_SUBSET
        case (TOK_PI_TAG)
          wf_stack(1) = wf_stack(1) + 1
          nextDTDState = ST_DTD_START_PI
        case (TOK_BANG_TAG)
          wf_stack(1) = wf_stack(1) + 1
          nextDTDState = ST_DTD_BANG_TAG
        case default
          call add_error(fx%error_stack, "Unexpected token in document subset")
          return
        end select


      case (ST_DTD_BANG_TAG)
        !write(*,*) 'ST_DTD_BANG_TAG'
        select case (fx%tokenType)
        case (TOK_OPEN_SB)
          nextDTDState = ST_DTD_START_SECTION_DECL
        case (TOK_OPEN_COMMENT)
          nextDTDState = ST_DTD_START_COMMENT
        case (TOK_NAME)
          if (str_vs(fx%token)=='ATTLIST') then
            nextDTDState = ST_DTD_ATTLIST
          elseif (str_vs(fx%token)=='ELEMENT') then
            nextDTDState = ST_DTD_ELEMENT
          elseif (str_vs(fx%token)=='ENTITY') then
            nextDTDState = ST_DTD_ENTITY
          elseif (str_vs(fx%token)=='NOTATION') then
            nextDTDState = ST_DTD_NOTATION
          endif
        end select

      case (ST_DTD_START_SECTION_DECL)
        !write(*,*) "ST_DTD_START_SECTION_DECL"
        select case (fx%tokenType)
        case (TOK_NAME)
          if (str_vs(fx%token)=="IGNORE") then
            if (fx%context/=CTXT_IN_DTD.or.reading_main_file(fb)) then
              call add_error(fx%error_stack, "IGNORE section only allowed in external subset.")
              return
            else
              ignoreDepth = 0
              fx%context = CTXT_IGNORE
              nextDTDState = ST_DTD_FINISH_SECTION_DECL
            endif
          elseif (str_vs(fx%token)=="INCLUDE") then
            if (fx%context/=CTXT_IN_DTD.or.reading_main_file(fb)) then
              call add_error(fx%error_stack, "INCLUDE section only allowed in external subset.")
              return
            else
              nextDTDState = ST_DTD_FINISH_SECTION_DECL
            endif
          else
            call add_error(fx%error_stack, "Unknown keyword found in marked section declaration.")
          endif
        end select


      case (ST_DTD_FINISH_SECTION_DECL)
        !write(*,*) "ST_FINISH_SECTION_DECL"
        select case (fx%tokenType)
        case (TOK_OPEN_SB)
          wf_stack(1) = wf_stack(1) + 1
          if (fx%context==CTXT_IGNORE) then
            nextDTDState = ST_DTD_IN_IGNORE_SECTION
            ignoreDepth = ignoreDepth + 1
          else
            nextDTDState = ST_DTD_SUBSET
          endif
        end select


      case (ST_DTD_IN_IGNORE_SECTION)
        !write(*,*) "ST_IN_IGNORE_SECTION"
        select case (fx%tokenType)
        case (TOK_SECTION_START)
          wf_stack(1) = wf_stack(1) + 2
          ignoreDepth = ignoreDepth + 1
          nextDTDState = ST_DTD_IN_IGNORE_SECTION
        case (TOK_SECTION_END)
          wf_stack(1) = wf_stack(1) - 2
          ignoreDepth = ignoreDepth - 1
          if (ignoreDepth==0) then
            fx%context = CTXT_IN_DTD
            nextDTDState = ST_DTD_SUBSET
          else
            nextDTDState = ST_DTD_IN_IGNORE_SECTION
          endif
        end select

      case (ST_DTD_START_PI)
        !write(*,*)'ST_DTD_START_PI'
        select case (fx%tokenType)
        case (TOK_NAME)
          if (namespaces_) then
            nameOk = checkNCName(str_vs(fx%token), fx%xds%xml_version)
          else
            nameOk = checkName(str_vs(fx%token), fx%xds%xml_version)
          endif
          if (nameOk) then
            if (str_vs(fx%token)=='xml') then
              call add_error(fx%error_stack, "XML declaration must be at start of document")
              return
            elseif (checkPITarget(str_vs(fx%token), fx%xds%xml_version)) then
              nextDTDState = ST_DTD_PI_CONTENTS
              fx%name => fx%token
              fx%token => null()
            else
              call add_error(fx%error_stack, "Invalid PI target name")
              return
            endif
          endif
        end select

      case (ST_DTD_PI_CONTENTS)
        !write(*,*)'ST_DTD_PI_CONTENTS'
        if (validCheck) then
          if (fx%context==CTXT_IN_DTD.and.wf_stack(1)==0) then
            call add_error(fx%error_stack, &
              "PI not balanced in parameter entity")
            return
          endif
          if (len(fx%elstack)>0) then
            elem => &
              get_element(fx%xds%element_list, get_top_elstack(fx%elstack))
            if (associated(elem)) then
              if (elem%empty) then
                call add_error(fx%error_stack, "Content inside empty element")
              endif
            endif
          endif
        endif
        wf_stack(1) = wf_stack(1) - 1
          
        select case(fx%tokenType)
        case (TOK_CHAR)
          if (present(processingInstruction_handler)) then
            call processingInstruction_handler(str_vs(fx%name), str_vs(fx%token))
            if (fx%state==ST_STOP) return
          endif
          deallocate(fx%name)
          nextDTDState = ST_DTD_PI_END
        case (TOK_PI_END)
          if (present(processingInstruction_handler)) then
            call processingInstruction_handler(str_vs(fx%name), '')
            if (fx%state==ST_STOP) return
          endif
          deallocate(fx%name)
          nextDTDState = ST_DTD_SUBSET
        end select

      case (ST_DTD_PI_END)
        !write(*,*)'ST_DTD_PI_END'
        select case(fx%tokenType)
        case (TOK_PI_END)
          nextDTDState = ST_DTD_SUBSET
        end select

      case (ST_DTD_START_COMMENT)
        !write(*,*)'ST_DTD_START_COMMENT'
        select case (fx%tokenType)
        case (TOK_CHAR)
          fx%name => fx%token
          nullify(fx%token)
          nextDTDState = ST_DTD_COMMENT_END
        end select

      case (ST_DTD_COMMENT_END)
        !write(*,*)'ST_DTD_COMMENT_END'
        if (validCheck) then
          if (wf_stack(1)==0) then
            call add_error(fx%error_stack, &
              "Comment not balanced in entity")
            return
          endif
          if (len(fx%elstack)>0) then
            elem => &
              get_element(fx%xds%element_list, get_top_elstack(fx%elstack))
            if (associated(elem)) then
              if (elem%empty) then
                call add_error(fx%error_stack, "Content inside empty element")
              endif
            endif
          endif
        endif
        wf_stack(1) = wf_stack(1) - 1

        select case (fx%tokenType)
        case (TOK_COMMENT_END)
          if (present(comment_handler)) then
            call comment_handler(str_vs(fx%name))
            if (fx%state==ST_STOP) return
          endif
          deallocate(fx%name)
          nextDTDState = ST_DTD_SUBSET
        end select

      case (ST_DTD_ATTLIST)
        !write(*,*) 'ST_DTD_ATTLIST'
        select case (fx%tokenType)
        case (TOK_NAME)
          if (namespaces_) then
            nameOk = checkQName(str_vs(fx%token), fx%xds%xml_version)
          else
            nameOk = checkName(str_vs(fx%token), fx%xds%xml_version)
          endif
          if (.not.nameOk) then
            call add_error(fx%error_stack, "Invalid element name for ATTLIST")
            return
          endif
          if (existing_element(fx%xds%element_list, str_vs(fx%token))) then
            elem => get_element(fx%xds%element_list, str_vs(fx%token))
          else
            elem => add_element(fx%xds%element_list, str_vs(fx%token))
          endif
          nextDTDState = ST_DTD_ATTLIST_CONTENTS
        end select

      case (ST_DTD_ATTLIST_CONTENTS)
        !write(*,*) 'ST_DTD_ATTLIST_CONTENTS'
        select case (fx%tokenType)
        case (TOK_ENTITY)
          !Weve found a PEref in the middle of the element contents
          ! Leave DTD state as it is & expand the entity ...
          nextState = ST_START_PE
        case (TOK_DTD_CONTENTS)
          if (processDTD) then
            call parse_dtd_attlist(str_vs(fx%token), fx%xds%xml_version, &
              namespaces_, validCheck, fx%error_stack, elem, &
              internal=reading_main_file(fb))
          else
#ifdef PGF90
            call parse_dtd_attlist(str_vs(fx%token), fx%xds%xml_version, &
              namespaces_, validCheck, fx%error_stack, nullElement, &
              internal=reading_main_file(fb))
#else
            call parse_dtd_attlist(str_vs(fx%token), fx%xds%xml_version, &
              namespaces_, validCheck, fx%error_stack, null(), &
              internal=reading_main_file(fb))
#endif
          endif
          if (in_error(fx%error_stack)) return
          ! Normalize attribute values in attlist
          if (processDTD) then
            do i = 1, size(elem%attlist%list)
              if (associated(elem%attlist%list(i)%default)) then
                tempString => elem%attlist%list(i)%default
                elem%attlist%list(i)%default => &
                  normalize_attribute_text(fx, tempString)
                deallocate(tempString)
                if (in_error(fx%error_stack)) return
              endif
            enddo
          endif
          nextDTDState = ST_DTD_ATTLIST_END
        case (TOK_END_TAG)
          if (validCheck) then
            if (wf_stack(1)==0) then
              call add_error(fx%error_stack, &
                "ATTLIST not balanced in parameter entity")
              return
            endif
          endif
          wf_stack(1) = wf_stack(1) - 1
          if (processDTD) then
            call parse_dtd_attlist("", fx%xds%xml_version, &
              namespaces_, validCheck, fx%error_stack, elem, &
              internal=reading_main_file(fb))
          else
#ifdef PGF90
            call parse_dtd_attlist("", fx%xds%xml_version, &
              namespaces_, validCheck, fx%error_stack, nullElement, &
              internal=reading_main_file(fb))
#else
            call parse_dtd_attlist("", fx%xds%xml_version, &
              namespaces_, validCheck, fx%error_stack, null(), &
              internal=reading_main_file(fb))
#endif
          endif
          if (in_error(fx%error_stack)) return
          if (processDTD) then
            if (present(attributeDecl_handler)) then
              call report_declarations(elem, attributeDecl_handler)
              if (fx%state==ST_STOP) return
            endif
          endif
          nextDTDState = ST_DTD_SUBSET
        end select

      case (ST_DTD_ATTLIST_END)
        !write(*,*) 'ST_DTD_ATTLIST_END'
        select case (fx%tokenType)
        case (TOK_END_TAG)
          if (validCheck) then
            if (wf_stack(1)==0) then
              call add_error(fx%error_stack, &
                "ATTLIST not balanced in parameter entity")
              return
            endif
          endif
          wf_stack(1) = wf_stack(1) - 1
          if (processDTD) then
            if (present(attributeDecl_handler)) then
              call report_declarations(elem, attributeDecl_handler)
              if (fx%state==ST_STOP) return
            endif
          endif
          nextDTDState = ST_DTD_SUBSET
        end select

      case (ST_DTD_ELEMENT)
        select case (fx%tokenType)
        case (TOK_NAME)
          if (namespaces_) then
            nameOk = checkQName(str_vs(fx%token), fx%xds%xml_version)
          else
            nameOk = checkName(str_vs(fx%token), fx%xds%xml_version)
          endif
          if (.not.nameOk) then
            call add_error(fx%error_stack, "Invalid name for ELEMENT")
            return
          endif
          fx%name => fx%token
          fx%token => null()
          nextDTDState = ST_DTD_ELEMENT_CONTENTS
        end select

      case (ST_DTD_ELEMENT_CONTENTS)
        !write(*,*)'ST_DTD_ELEMENT_CONTENTS'
        select case (fx%tokenType)
        case (TOK_OPEN_PAR)
          ! increment well_formedness
          wf_stack(1) = wf_stack(1) + 1
          nextDTDState = ST_DTD_ELEMENT_CONTENTS
        case (TOK_CLOSE_PAR)
          ! increment well_formedness
          wf_stack(1) = wf_stack(1) - 1
          nextDTDState = ST_DTD_ELEMENT_CONTENTS
        case (TOK_ENTITY)
          !Weve found a PEref in the middle of the element contents
          ! Leave DTD state as it is & expand the entity ...
          nextState = ST_START_PE
        case (TOK_DTD_CONTENTS)
          if (declared_element(fx%xds%element_list, str_vs(fx%name))) then
            if (validCheck) then
              call add_error(fx%error_stack, "Duplicate Element declaration")
              return
            else
              ! Ignore contents ...
              elem => null()
            endif
          elseif (processDTD) then
            if (existing_element(fx%xds%element_list, str_vs(fx%name))) then
              elem => get_element(fx%xds%element_list, str_vs(fx%name))
            else
              elem => add_element(fx%xds%element_list, str_vs(fx%name))
            endif
          else
            elem => null()
          endif
          call parse_dtd_element(str_vs(fx%token), fx%xds%xml_version, fx%error_stack, &
            elem, reading_main_file(fb))
          if (in_error(fx%error_stack)) return
          nextDTDState = ST_DTD_ELEMENT_END
        end select

      case (ST_DTD_ELEMENT_END)
        !write(*,*)'ST_DTD_ELEMENT_END'
        select case (fx%tokenType)
        case (TOK_END_TAG)
          if (validCheck) then
            if (wf_stack(1)==0) then
              call add_error(fx%error_stack, &
                "ELEMENT not balanced in parameter entity")
              return
            endif
          endif
          wf_stack(1) = wf_stack(1) - 1
          if (processDTD.and.associated(elem)) then
            if (present(elementDecl_handler)) then
              call elementDecl_handler(str_vs(fx%name), str_vs(elem%model))
              if (fx%state==ST_STOP) return
            endif
          endif
          deallocate(fx%name)
          nextDTDState = ST_DTD_SUBSET
        end select

      case (ST_DTD_ENTITY)
        !write(*,*) 'ST_DTD_ENTITY'
        select case (fx%tokenType)
        case (TOK_NAME)
          if (str_vs(fx%token)=="%") then
            pe = .true.
            ! this will be a PE
            nextDTDState = ST_DTD_ENTITY_PE
          else
            pe = .false.
            if (namespaces_) then
              nameOk = checkNCName(str_vs(fx%token), fx%xds%xml_version)
            else
              nameOk = checkName(str_vs(fx%token), fx%xds%xml_version)
            endif
            if (.not.nameOk) then
              call add_error(fx%error_stack, &
                "Illegal name for general entity")
              return
            endif
            fx%name => fx%token
            fx%token => null()
            nextDTDState = ST_DTD_ENTITY_ID
          endif
        end select

      case (ST_DTD_ENTITY_PE)
        !write(*,*) 'ST_DTD_ENTITY_PE'
        select case (fx%tokenType)
        case (TOK_NAME)
          if (namespaces_) then
            nameOk = checkNCName(str_vs(fx%token), fx%xds%xml_version)
          else
            nameOk = checkName(str_vs(fx%token), fx%xds%xml_version)
          endif
          if (.not.nameOk) then
            call add_error(fx%error_stack, &
              "Illegal name for parameter entity")
            return
          endif
          fx%name => fx%token
          fx%token => null()
          nextDTDState = ST_DTD_ENTITY_ID
        end select

      case (ST_DTD_ENTITY_ID)
        !write(*,*) 'ST_DTD_ENTITY_ID'
        select case (fx%tokenType)
        case (TOK_NAME)
          if (str_vs(fx%token) == "PUBLIC") then
            nextDTDState = ST_DTD_ENTITY_PUBLIC
          elseif (str_vs(fx%token) == "SYSTEM") then
            nextDTDState = ST_DTD_ENTITY_SYSTEM
          else
            call add_error(fx%error_stack, "Unexpected token in ENTITY")
            return
          endif
        case (TOK_CHAR)
          if (reading_main_file(fb)) then
            tempString => fx%token
          else
            tempString => expand_pe_text(fx, fx%token, fb)
          endif
          fx%attname => expand_entity_value_alloc(tempString, fx%xds, fx%error_stack)
          if (reading_main_file(fb)) then
            tempString => null()
          else
            deallocate(tempString)
          endif
          if (in_error(fx%error_stack)) return
          nextDTDState = ST_DTD_ENTITY_END
        case default
          call add_error(fx%error_stack, "Unexpected token in ENTITY")
          return
        end select

      case (ST_DTD_ENTITY_PUBLIC)
        !write(*,*) 'ST_DTD_ENTITY_PUBLIC'
        select case (fx%tokenType)
        case (TOK_CHAR)
          if (checkPublicId(str_vs(fx%token))) then
            fx%publicId => fx%token
            fx%token => null()
            nextDTDState = ST_DTD_ENTITY_SYSTEM
          else
            call add_error(fx%error_stack, "Invalid PUBLIC id in ENTITY")
            return
          endif
        case default
          call add_error(fx%error_stack, "Unexpected token in ENTITY")
          return
        end select

      case (ST_DTD_ENTITY_SYSTEM)
        !write(*,*) 'ST_DTD_ENTITY_SYSTEM'
        select case (fx%tokenType)
        case (TOK_CHAR)
          fx%systemId => fx%token
          fx%token => null()
          nextDTDState = ST_DTD_ENTITY_NDATA
        case default
          call add_error(fx%error_stack, "Unexpected token in ENTITY")
          return
        end select

      case (ST_DTD_ENTITY_NDATA)
        !write(*,*) 'ST_DTD_ENTITY_NDATA'
        select case (fx%tokenType)
        case (TOK_END_TAG)
          if (validCheck) then
            if (wf_stack(1)==0) then
              call add_error(fx%error_stack, &
                "ENTITY not balanced in parameter entity")
              return
            endif
          endif
          wf_stack(1) = wf_stack(1) - 1
          if (processDTD) then
            call add_entity
            if (in_error(fx%error_stack)) return
            if (fx%state==ST_STOP) return
          endif
          deallocate(fx%name)
          if (associated(fx%attname)) deallocate(fx%attname)
          if (associated(fx%systemId)) deallocate(fx%systemId)
          if (associated(fx%publicId)) deallocate(fx%publicId)
          if (associated(fx%Ndata)) deallocate(fx%Ndata)
          nextDTDState = ST_DTD_SUBSET
        case (TOK_NAME)
          if (str_vs(fx%token)=='NDATA') then
            if (pe) then
              call add_error(fx%error_stack, "Parameter entity cannot have NDATA declaration")
              return
            endif
            nextDTDState = ST_DTD_ENTITY_NDATA_VALUE
          else
            call add_error(fx%error_stack, "Unexpected token in ENTITY")
            return
          endif
        case default
          call add_error(fx%error_stack, "Unexpected token in ENTITY")
          return
        end select

      case (ST_DTD_ENTITY_NDATA_VALUE)
        !write(*,*) 'ST_DTD_ENTITY_NDATA_VALUE'
        !check is a name and exists in notationlist
        select case (fx%tokenType)
        case (TOK_NAME)
          if (namespaces_) then
            nameOk = checkNCName(str_vs(fx%token), fx%xds%xml_version)
          else
            nameOk = checkName(str_vs(fx%token), fx%xds%xml_version)
          endif
          if (.not.nameOk) then
            call add_error(fx%error_stack, "Invalid name for Notation")
            return
          endif
          fx%Ndata => fx%token
          fx%token => null()
          nextDTDState = ST_DTD_ENTITY_END
        case default
          call add_error(fx%error_stack, "Unexpected token in ENTITY")
          return
        end select

      case (ST_DTD_ENTITY_END)
        !write(*,*) 'ST_DTD_ENTITY_END'
        select case (fx%tokenType)
        case (TOK_END_TAG)
          if (validCheck) then
            if (wf_stack(1)==0) then
              call add_error(fx%error_stack, &
                "ENTITY not balanced in parameter entity")
              return
            endif
          endif
          wf_stack(1) = wf_stack(1) - 1
          if (processDTD) then
            call add_entity
            if (in_error(fx%error_stack)) return
            if (fx%state==ST_STOP) return
          endif
          deallocate(fx%name)
          if (associated(fx%attname)) deallocate(fx%attname)
          if (associated(fx%systemId)) deallocate(fx%systemId)
          if (associated(fx%publicId)) deallocate(fx%publicId)
          if (associated(fx%Ndata)) deallocate(fx%Ndata)
          nextDTDState = ST_DTD_SUBSET
        case default
          call add_error(fx%error_stack, "Unexpected token at end of ENTITY")
          return
        end select

      case (ST_DTD_NOTATION)
        !write(*,*) 'ST_DTD_NOTATION'
        select case (fx%tokenType)
        case (TOK_NAME)
          if (namespaces_) then
            nameOk = checkNCName(str_vs(fx%token), fx%xds%xml_version)
          else
            nameOk = checkName(str_vs(fx%token), fx%xds%xml_version)
          endif
          if (.not.nameOk) then
            call add_error(fx%error_stack, "Invalid name for Notation")
            return
          endif
          fx%name => fx%token
          fx%token => null()
          nextDTDState = ST_DTD_NOTATION_ID
        case default
          call add_error(fx%error_stack, "Unexpected token in NOTATION")
          return
        end select

      case (ST_DTD_NOTATION_ID)
        !write(*,*)'ST_DTD_NOTATION_ID'
        select case (fx%tokenType)
        case (TOK_NAME)
          if (str_vs(fx%token)=='SYSTEM') then
            nextDTDState = ST_DTD_NOTATION_SYSTEM
          elseif (str_vs(fx%token)=='PUBLIC') then
            nextDTDState = ST_DTD_NOTATION_PUBLIC
          else
            call add_error(fx%error_stack, "Unexpected token after NOTATION")
            return
          endif
        case default
          call add_error(fx%error_stack, "Unexpected token after NOTATION")
          return
        end select

      case (ST_DTD_NOTATION_SYSTEM)
        !write(*,*)'ST_DTD_NOTATION_SYSTEM'
        select case (fx%tokenType)
        case (TOK_CHAR)
          fx%systemId => fx%token
          fx%token => null()
          nextDTDState = ST_DTD_NOTATION_END
        case default
          call add_error(fx%error_stack, "Unexpected token in NOTATION")
          return
        end select

      case (ST_DTD_NOTATION_PUBLIC)
        !write(*,*)'ST_DTD_NOTATION_PUBLIC'
        select case (fx%tokenType)
        case (TOK_CHAR)
          if (checkPublicId(str_vs(fx%token))) then
            fx%publicId => fx%token
            fx%token => null()
            nextDTDState = ST_DTD_NOTATION_PUBLIC_2
          else
            call add_error(fx%error_stack, "Invalid PUBLIC id in NOTATION")
            return
          endif
        case default
          call add_error(fx%error_stack, "Unexpected token in NOTATION")
          return
        end select

      case (ST_DTD_NOTATION_PUBLIC_2)
        !write(*,*)'ST_DTD_NOTATION_PUBLIC_2'
        select case (fx%tokenType)
        case (TOK_END_TAG)
          if (validCheck) then
            if (wf_stack(1)==0) then
              call add_error(fx%error_stack, &
                "NOTATION not balanced in parameter entity")
              return
            endif
            if (notation_exists(fx%nlist, str_vs(fx%name))) then
              call add_error(fx%error_stack, "Duplicate notation declaration")
              return
            endif
          endif
          wf_stack(1) = wf_stack(1) - 1
          if (processDTD) then
            call add_notation(fx%nlist, str_vs(fx%name), publicId=str_vs(fx%publicId))
            if (present(notationDecl_handler)) then
              call notationDecl_handler(str_vs(fx%name), publicId=str_vs(fx%publicId), systemId="")
              if (fx%state==ST_STOP) return
            endif
          endif
          deallocate(fx%name)
          deallocate(fx%publicId)
          nextDTDState = ST_DTD_SUBSET
        case (TOK_CHAR)
          fx%systemId => fx%token
          fx%token => null()
          nextDTDState = ST_DTD_NOTATION_END
        end select

      case (ST_DTD_NOTATION_END)
        !write(*,*)'ST_DTD_NOTATION_END'
        select case (fx%tokenType)
        case (TOK_END_TAG)
          if (validCheck) then
            if (wf_stack(1)==0) then
              call add_error(fx%error_stack, &
                "NOTATION not balanced in parameter entity")
              return
            endif
            if (notation_exists(fx%nlist, str_vs(fx%name))) then
              call add_error(fx%error_stack, "Duplicate notation declaration")
              return
            endif
          endif
          wf_stack(1) = wf_stack(1) - 1
          if (processDTD) then
            URIref => parseURI(str_vs(fx%systemId))
            if (.not.associated(URIref)) then
              call add_error(fx%error_stack, "Invalid SYSTEM literal")
              return
            endif
            if (hasFragment(URIref)) then
              call destroyURI(URIref)
              call add_error(fx%error_stack, "SYSTEM literal may not contain fragment")
              return
            endif
            ! We aren't ever going to do anything with this URI,
            ! since we don't do NOTATIONs.
            ! Throw it away again
            call destroyURI(URIref)
            if (associated(fx%publicId)) then
              call add_notation(fx%nlist, str_vs(fx%name), &
                publicId=str_vs(fx%publicId), systemId=str_vs(fx%systemId))
              if (present(notationDecl_handler)) then
                call notationDecl_handler(str_vs(fx%name), &
                publicId=str_vs(fx%publicId), systemId=str_vs(fx%systemId))
                if (fx%state==ST_STOP) return
              endif
            else
              call add_notation(fx%nlist, str_vs(fx%name), &
                systemId=str_vs(fx%systemId))
              if (present(notationDecl_handler)) then
                call notationDecl_handler(str_vs(fx%name), &
                publicId="", systemId=str_vs(fx%systemId))
                if (fx%state==ST_STOP) return
              endif
            endif
          endif
          if (associated(fx%publicId)) deallocate(fx%publicId)
          deallocate(fx%systemId)
          deallocate(fx%name)
          nextDTDState = ST_DTD_SUBSET
        case default
          call add_error(fx%error_stack, "Unexpected token in NOTATION")
          return
        end select

      end select

      if (nextDTDState==ST_DTD_NULL) then
        call add_error(fx%error_stack, &
          "Bad token found in DTD parsing")
      else
        fx%state_dtd = nextDTDState
      endif
      if (nextState==ST_NULL) &
        nextState = ST_IN_SUBSET

    end subroutine parseDTD

    subroutine open_tag
      elem => get_element(fx%xds%element_list, str_vs(fx%name))
      if (associated(elem)) then
        if (validCheck) then
          call checkAttributes(elem, fx%attributes)
          if (.not.checkContentModel(fx%elstack, str_vs(fx%name))) then
            call add_error(fx%error_stack, &
              "Element "//str_vs(fx%name)//" not permitted in this context")
            return
          endif
        else
          call getDefaultAttributes(elem, fx%attributes)
        endif
      else
        if (validCheck) then
          call add_error(fx%error_stack, &
          "Trying to use an undeclared element")
          return
        endif
      endif
      ! Check for namespace changes
      if (namespaces_) then
        call checkNamespaces(fx%attributes, fx%nsDict, &
        len(fx%elstack), fx%xds, namespace_prefixes_, xmlns_uris_, &
        fx%error_stack, startInCharData_, &
        startPrefixMapping_handler, endPrefixMapping_handler)
        if (fx%state==ST_STOP) return
      endif
      if (in_error(fx%error_stack)) return
      call checkXmlAttributes
      if (in_error(fx%error_stack)) return
      if (size(extEntStack)>0) then
        if (len(fx%elstack)==extEntStack(1)) &
          ! This is a top-level element in the current entity
          call setBase(fx%attributes, expressURI(fb%f(1)%baseURI))
      endif
      if (namespaces_.and.getURIofQName(fx,str_vs(fx%name))==invalidNS) then
        ! no namespace was found for the current element
        if (.not.startInCharData_) then
          ! but we ignore this if we are parsing an entity through DOM
          call add_error(fx%error_stack, "No namespace found for current element")
          return
        elseif (present(startElement_handler)) then
          ! Record it as having an empty URI
          call startElement_handler("", &
            getlocalNameofQName(str_vs(fx%name)), &
            str_vs(fx%name), fx%attributes)
          if (fx%state==ST_STOP) return
        endif
      elseif (namespaces_) then
        ! Normal state of affairs
        if (present(startElement_handler)) then
          call startElement_handler(getURIofQName(fx, str_vs(fx%name)), &
          getlocalNameofQName(str_vs(fx%name)), &
          str_vs(fx%name), fx%attributes)
          if (fx%state==ST_STOP) return
        endif
      else
        ! Non-namespace aware processing
        if (present(startElement_handler)) then
          call startElement_handler("", "", &
          str_vs(fx%name), fx%attributes)
          if (fx%state==ST_STOP) return
        endif
      endif
      if (validCheck) then
        call push_elstack(fx%elstack, str_vs(fx%name), elem%cp)
      else
        call push_elstack(fx%elstack, str_vs(fx%name))
      endif
      call reset_dict(fx%attributes)
    end subroutine open_tag

    subroutine close_tag
      character :: dummy
      wf_stack(1) = wf_stack(1) - 1
      if (wf_stack(1)<0) then
        call add_error(fx%error_stack, &
          'Ill-formed entity')
        return
      endif
      if (str_vs(fx%name)/=get_top_elstack(fx%elstack)) then
        call add_error(fx%error_stack, &
          "Mismatching close tag: trying to close "//get_top_elstack(fx%elstack) &
          //" with "//str_vs(fx%name))
        return
      endif
      if (validCheck) then
        if (.not.checkContentModelToEnd(fx%elstack)) then
          call add_error(fx%error_stack, &
            "Failed to fulfil content model for "//str_vs(fx%name))
          return
        endif
      endif
      dummy = pop_elstack(fx%elstack)
      if (present(endElement_handler)) then
        if (namespaces_.and.getURIofQName(fx,str_vs(fx%name))==invalidNS) then
          ! no namespace was found for the current element, we must be
          ! closing inside a DOM entity.
          ! Record it as having an empty URI
          call endElement_handler("", &
            getlocalNameofQName(str_vs(fx%name)), &
            str_vs(fx%name))
        elseif (namespaces_) then
          ! Normal state of affairs
          call endElement_handler(getURIofQName(fx, str_vs(fx%name)), &
            getlocalnameofQName(str_vs(fx%name)), &
            str_vs(fx%name))
        else
          ! Non-namespace-aware processing:
          call endElement_handler("", "", &
            str_vs(fx%name))
        endif
        if (fx%state==ST_STOP) return
      endif
      if (namespaces_) then
        call checkEndNamespaces(fx%nsDict, len(fx%elstack), &
        endPrefixMapping_handler)
        if (fx%state==ST_STOP) return
      endif
    end subroutine close_tag

    subroutine add_entity
      !Parameter or General Entity?
      logical :: wfc
      wfc = fb%f(1)%pe.or.inExtSubset
      if (pe) then
        !Does entity with this name exist?
        if (.not.existing_entity(fx%xds%PEList, str_vs(fx%name))) then
          ! Internal or external?
          if (associated(fx%attname)) then ! it's internal
            call register_internal_PE(fx%xds, &
              name=str_vs(fx%name), text=str_vs(fx%attname), &
              wfc=wfc, baseURI=copyURI(fb%f(1)%baseURI))
            ! FIXME need to expand value here before reporting ...
            if (present(internalEntityDecl_handler)) then
              call internalEntityDecl_handler('%'//str_vs(fx%name), str_vs(fx%attname))
              if (fx%state==ST_STOP) return
            endif
          else ! PE can't have Ndata declaration
            URIref => parseURI(str_vs(fx%systemId))
            if (.not.associated(URIref)) then
              call add_error(fx%error_stack, "Invalid URI specified for SYSTEM")
            elseif (hasFragment(URIref)) then
              call add_error(fx%error_stack, "Fragment not permitted on SYSTEM URI")
              call destroyURI(URIref)
            else
              newURI => rebaseURI(fb%f(1)%baseURI, URIref)
              call destroyURI(URIref)
              if (associated(fx%publicId)) then
                call register_external_PE(fx%xds, name=str_vs(fx%name), &
                  systemId=str_vs(fx%systemId), publicId=str_vs(fx%publicId), &
                  wfc=wfc, baseURI=newURI)
                if (present(externalEntityDecl_handler)) &
                  call externalEntityDecl_handler('%'//str_vs(fx%name), &
                  systemId=str_vs(fx%systemId), publicId=str_vs(fx%publicId))
              else
                call register_external_PE(fx%xds, name=str_vs(fx%name), &
                  systemId=str_vs(fx%systemId), &
                  wfc=wfc, baseURI=newURI)
                if (present(externalEntityDecl_handler)) &
                  call externalEntityDecl_handler('%'//str_vs(fx%name), &
                  systemId=str_vs(fx%systemId), publicId="")
              endif
            endif
          endif
          ! else we ignore it
        endif
      else !It's a general entity
        if (.not.existing_entity(fx%xds%entityList, str_vs(fx%name))) then
          ! Internal or external?
          if (associated(fx%attname)) then ! it's internal
            call register_internal_GE(fx%xds, name=str_vs(fx%name), &
              text=str_vs(fx%attname), &
              wfc=wfc, baseURI=copyURI(fb%f(1)%baseURI))
            if (present(internalEntityDecl_handler)) then
              call internalEntityDecl_handler(str_vs(fx%name),&
              str_vs(fx%attname))
              if (fx%state==ST_STOP) return
            endif
          else
            URIref => parseURI(str_vs(fx%systemId))
            if (.not.associated(URIref)) then
              call add_error(fx%error_stack, "Invalid URI specified for SYSTEM")
            elseif (hasFragment(URIref)) then
              call add_error(fx%error_stack, "Fragment not permitted on SYSTEM URI")
              call destroyURI(URIref)
            else
              newURI => rebaseURI(fb%f(1)%baseURI, URIref)
              call destroyURI(URIref)
              if (associated(fx%publicId).and.associated(fx%Ndata)) then
                call register_external_GE(fx%xds, name=str_vs(fx%name), &
                  systemId=str_vs(fx%systemId), publicId=str_vs(fx%publicId), &
                  notation=str_vs(fx%Ndata), &
                  wfc=wfc, baseURI=newURI)
                if (present(unparsedEntityDecl_handler)) &
                  call unparsedEntityDecl_handler(str_vs(fx%name), &
                  systemId=str_vs(fx%systemId), publicId=str_vs(fx%publicId), &
                  notation=str_vs(fx%Ndata))
              elseif (associated(fx%Ndata)) then
                call register_external_GE(fx%xds, name=str_vs(fx%name), &
                  systemId=str_vs(fx%systemId), notation=str_vs(fx%Ndata), &
                  wfc=wfc, baseURI=newURI)
                if (present(unparsedEntityDecl_handler)) &
                  call unparsedEntityDecl_handler(str_vs(fx%name), publicId="", &
                  systemId=str_vs(fx%systemId), notation=str_vs(fx%Ndata))
              elseif (associated(fx%publicId)) then
                call register_external_GE(fx%xds, name=str_vs(fx%name), &
                  systemId=str_vs(fx%systemId), publicId=str_vs(fx%publicId), &
                  wfc=wfc, baseURI=newURI)
                if (present(externalEntityDecl_handler)) &
                  call externalEntityDecl_handler(str_vs(fx%name), &
                  systemId=str_vs(fx%systemId), publicId=str_vs(fx%publicId))
              else
                call register_external_GE(fx%xds, name=str_vs(fx%name), &
                  systemId=str_vs(fx%systemId), wfc=wfc, baseURI=newURI)
                if (present(externalEntityDecl_handler)) &
                  call externalEntityDecl_handler(str_vs(fx%name), &
                  systemId=str_vs(fx%systemId), publicId="")
              endif
            endif
          endif
        endif
      endif
    end subroutine add_entity

    subroutine getDefaultAttributes(el, dict)
      type(element_t), pointer :: el
      type(dictionary_t), intent(inout) :: dict

      type(attribute_t), pointer :: att
      integer :: i, ind
      character, pointer :: attValue(:)

      do i = 1, size(el%attlist%list)
        att => el%attlist%list(i)
        call get_att_index_pointer(dict, str_vs(att%name), ind, attValue)
        if (.not.associated(attValue)) then
          if (att%attDefault==ATT_DEFAULT &
            .or.att%attDefault==ATT_FIXED) then
            call add_item_to_dict(dict, &
              str_vs(att%name), str_vs(att%default), &
              specified=.false., declared=.true.)
          endif
        endif
      end do
    end subroutine getDefaultAttributes

    subroutine checkAttributes(el, dict)
      type(element_t), pointer :: el
      type(dictionary_t), intent(inout) :: dict

      integer :: i, j

      type(attribute_t), pointer :: att
      type(string_list) :: s_list
      character, pointer :: attValue(:), s(:)

      integer :: ind
      logical, allocatable :: attributesLeft(:)
      allocate(attributesLeft(getLength(dict)))
      attributesLeft = .true.

      do i = 1, size(el%attlist%list)
        att => el%attlist%list(i)
        call get_att_index_pointer(dict, str_vs(att%name), ind, attValue)
        if (associated(attValue)) attributesLeft(ind) = .false.
        select case(att%attDefault)
        case (ATT_REQUIRED, ATT_IMPLIED, ATT_DEFAULT)
          if (.not.associated(attValue)) then
            if (att%attDefault==ATT_REQUIRED) then
              ! Validity Constraint: Required Attribute
              call add_error(fx%error_stack, &
                "REQUIRED attribute "//str_vs(att%name)//" not present")
              return
            elseif (att%attDefault==ATT_DEFAULT) then
              if (fx%xds%standalone.and..not.att%internal) then
                ! VC: Standalone document declaration"
                call add_error(fx%error_stack, &
                  "Externally-specifid default attribute used in non-standalone document")
                return
              else
                call add_item_to_dict(dict, &
                  str_vs(att%name), str_vs(att%default), &
                  specified=.false., declared=.true.)
              endif
            endif
          else
            select case(att%attType)
            case (ATT_ID)
              ! VC: ID
              if (namespaces_) then
                nameOk = checkNCName(str_vs(attValue), fx%xds%xml_version)
              else
                nameOk = checkName(str_vs(attValue), fx%xds%xml_version)
              endif
              if (.not.nameOk) then
                call add_error(fx%error_stack, &
                  "Attributes of type ID must have a value which is an XML Name")
                return
              endif
              if (registered_string(id_list, str_vs(attValue))) then
                call add_error(fx%error_stack, &
                  "Cannot declare the same ID twice")
                return
              endif
              call add_string(id_list, str_vs(attValue))
              call setIsId(dict, ind, .true.)
              ! We don't need to check for duplicate ID & xml:ids on the same
              ! element - if we are validating we'd already have noticed,
              ! if we're not, it's not possible
            case (ATT_IDREF)
              ! VC: IDREF
              if (namespaces_) then
                nameOk = checkNCName(str_vs(attValue), fx%xds%xml_version)
              else
                nameOk = checkName(str_vs(attValue), fx%xds%xml_version)
              endif
              if (.not.nameOk) then
                call add_error(fx%error_stack, &
                  "Attributes of type IDREF must have a value which is an XML Name")
                return
              endif
              ! FIXME in a namespaced document they must match QName
              if (.not.registered_string(idref_list, str_vs(attValue))) &
                call add_string(idref_list, str_vs(attValue))
            case (ATT_IDREFS)
              ! VC: IDREF
              if (namespaces_) then
                nameOk = checkNCNames(str_vs(attValue), fx%xds%xml_version)
              else
                nameOk = checkNames(str_vs(attValue), fx%xds%xml_version)
              endif
              if (.not.nameOk) then
                call add_error(fx%error_stack, &
                  "Attributes of type IDREFS must have a value which contains only XML Names")
                return
              endif
              call tokenize_and_add_strings(idref_list, str_vs(attValue), .true.)
            case (ATT_ENTITY)
              ! VC: Entity Name
              if (namespaces_) then
                nameOk = checkNCName(str_vs(attValue), fx%xds%xml_version)
              else
                nameOk = checkName(str_vs(attValue), fx%xds%xml_version)
              endif
              if (.not.nameOk) then
                call add_error(fx%error_stack, &
                  "Attributes of type ENTITY must have a value which is an XML Name")
                return
              endif
              ent => getEntityByName(fx%xds%entityList, str_vs(attValue))
              if (associated(ent)) then
                if (.not.is_unparsed_entity(ent)) then
                  ! Validity Constraint: Entity Name
                  call add_error(fx%error_stack, &
                    "Attribute "//str_vs(att%name) &
                    //" of element "//str_vs(el%name) &
                    //" declared as ENTITY refers to parsed entity")
                  return
                endif
              else
                ! Validity Constraint: Entity Name
                call add_error(fx%error_stack, &
                  "Attribute "//str_vs(att%name) &
                  //" of element "//str_vs(el%name) &
                  //" declared as ENTITY refers to non-existent entity")
                return
              endif
            case (ATT_ENTITIES)
              ! VC: Entity Name
              if (namespaces_) then
                nameOk = checkNCNames(str_vs(attValue), fx%xds%xml_version)
              else
                nameOk = checkNames(str_vs(attValue), fx%xds%xml_version)
              endif
              if (.not.nameOk) then
                call add_error(fx%error_stack, &
                  "Attributes of type ENTITIES must have a value which contains only XML Names")
                return
              endif
              s_list = tokenize_to_string_list(str_vs(attValue))
              do j = 1, size(s_list%list)
                s => s_list%list(j)%s
                ent => getEntityByName(fx%xds%entityList, str_vs(s))
                if (associated(ent)) then
                  if (.not.is_unparsed_entity(ent)) then
                    ! Validity Constraint: Entity Name
                    call add_error(fx%error_stack, &
                      "Attribute "//str_vs(att%name) &
                      //" of element "//str_vs(el%name) &
                      //" declared as ENTITIES refers to parsed entity "&
                      //str_vs(s))
                    call destroy_string_list(s_list)
                    return
                  endif
                else
                  ! Validity Constraint: Entity Name
                  call add_error(fx%error_stack, &
                    "Attribute "//str_vs(att%name) &
                    //" of element "//str_vs(el%name) &
                    //" declared as ENTITIES refers to non-existent entity "&
                    //str_vs(s))
                  call destroy_string_list(s_list)
                  return
                endif
              enddo
              call destroy_string_list(s_list)
            case (ATT_NMTOKEN)
              ! VC Name Token
              if (.not.checkNmtoken(str_vs(attValue), fx%xds%xml_version)) then
                call add_error(fx%error_stack, &
                  "Attributes of type NMTOKEN must have a value which is a NMTOKEN")
                return
              endif
            case (ATT_NMTOKENS)
              ! VC: Name Token
              if (.not.checkNmtokens(str_vs(attValue), fx%xds%xml_version)) then
                call add_error(fx%error_stack, &
                  "Attributes of type NMTOKENS must have a value which contain only NMTOKENs")
                return
              endif
            case (ATT_NOTATION)
              ! VC: Notation Attributes
              if (namespaces_) then
                nameOk = checkNCName(str_vs(attValue), fx%xds%xml_version)
              else
                nameOk = checkName(str_vs(attValue), fx%xds%xml_version)
              endif
              if (.not.nameOk) then
                call add_error(fx%error_stack, &
                  "Attributes of type NOTATION must have a value which is an XML Name")
                return
              endif
              if (.not.notation_exists(fx%nlist, str_vs(attValue))) then
                ! Validity Constraint: Notation Attributes
                call add_error(fx%error_stack, &
                  "Attribute "//str_vs(att%name) &
                  //" declared as NOTATION refers to non-existent notation "&
                  //str_vs(attValue))
                return
              endif
              if (att%attDefault==ATT_REQUIRED) then
                if (.not.registered_string(att%enumerations, str_vs(attValue))) &
                  call add_error(fx%error_stack, &
                  "NOTATION attribute is not among declared values.")
              endif
            case (ATT_ENUM)
              ! VC: Notation Attributes
              if (.not.checkNmtoken(str_vs(attValue), fx%xds%xml_version)) then
                call add_error(fx%error_stack, &
                  "Attributes of type ENUM must have a value which is an NMTOKEN")
                return
              endif
              if (.not.registered_string(att%enumerations, str_vs(attValue))) then
                ! Validity Constraint: Enumeration
                call add_error(fx%error_stack, &
                  "Attribute "//str_vs(att%name) &
                  //" of element "//str_vs(el%name) &
                  //" declared as ENUM refers to undeclared enumeration "&
                  //str_vs(attValue))
                return
              endif
            end select
          endif
        case (ATT_FIXED)
          if (associated(attValue)) then
            if (str_vs(att%default)//"x"/=str_vs(attValue)//"x") then
              ! Validity Constraint: Fixed Attribute Default
              call add_error(fx%error_stack, &
                "FIXED attribute has wrong value")
              return
            endif
          else
            if (fx%xds%standalone.and..not.att%internal) then
              ! VC: Standalone document declaration"
              call add_error(fx%error_stack, &
                "Externally-specified default attribute used in non-standalone document")
              return
            else
              call add_item_to_dict(dict, &
                str_vs(att%name), str_vs(att%default), &
                specified=.false., declared=.true.)
            endif
          endif
        end select
      enddo

      if (any(attributesLeft)) &
        call add_error(fx%error_stack, &
        "Undeclared attributes forbidden")

    end subroutine checkAttributes

    subroutine checkXMLAttributes
      integer :: ind
      character, pointer :: attValue(:)
      ! These must all be done with the name of the attribute,
      ! not the nsURI/localname pair, in case we are
      ! processing for a non-namespace aware application
      if (has_key(fx%attributes, 'xml:space')) then
        if (get_value(fx%attributes, 'xml:space')/='default' &
          .and. get_value(fx%attributes, 'xml:space')/='preserve') then
          call add_error(fx%error_stack, 'Illegal value of xml:space attribute')
          return
        endif
      endif
      call get_att_index_pointer(fx%attributes, "xml:id", ind, attValue)
      if (associated(attValue)) then
        ! Per xml:id spec, NCName even in non-namespace aware document
        if (.not.checkNCName(str_vs(attValue), fx%xds%xml_version)) then
          call add_error(fx%error_stack, &
            "xml:id attributes must have values which are NCNames")
          return
        elseif (registered_string(id_list, str_vs(attValue))) then
          call add_error(fx%error_stack, &
            "xml:id attributes must be unique within a document")
          return
        endif
        call add_string(id_list, str_vs(attValue))
        call setIsId(fx%attributes, ind, .true.)
      endif
      if (has_key(fx%attributes, "xml:base")) then
        URIref => parseURI(get_value(fx%attributes,"xml:base"))
        if (.not.associated(URIref)) then
          call add_error(fx%error_stack, &
            "Invalid URI reference specified for xml:base attribute")
        else
          call destroyURI(URIref)
        endif
      endif
      !if (has_key(fx%attributes, 'xml:lang')) then
      !   We never care about this at the SAX level.
      !endif
    end subroutine checkXMLAttributes

    subroutine endDTDchecks
      type(element_t), pointer :: el
      type(attribute_t), pointer :: att
      type(entity_t), pointer :: ent
      type(string_list) :: s_list
      character, pointer :: s(:)
      integer :: i, j, k

      if (present(FoX_endDTD_handler)) then
        fx%xds_used = .true.
        call FoX_endDTD_handler(fx%xds)
      endif
      if (present(endDTD_handler)) then
        call endDTD_handler()
        if (fx%state==ST_STOP) return
      endif
      ! Check that all notations used have been declared:
      if (validCheck) then
        do i = 1, size(fx%xds%entityList)
          ent => getEntityByIndex(fx%xds%entityList, i)
          if (str_vs(ent%notation)/="" &
            .and..not.notation_exists(fx%nlist, str_vs(ent%notation))) then
            ! Validity Constraint: Notation Declared
            call add_error(fx%error_stack, "Attempt to use undeclared notation")
            exit
          endif
        enddo
        validLoop: do i = 1, size(fx%xds%element_list%list)
          el => fx%xds%element_list%list(i)
          do j = 1, size(el%attlist%list)
            att => el%attlist%list(j)
            ! For NOTATION, need to check enumerated as well as default ...
            if (att%attType==ATT_NOTATION) then
              do k = 1, size(att%enumerations%list)
                s => att%enumerations%list(k)%s
                if (.not.notation_exists(fx%nlist, str_vs(s))) then
                  ! Validity Constraint: Notation Attributes
                  call add_error(fx%error_stack, &
                    "Enumerated NOTATION in "//str_vs(att%name) &
                    //" of element "//str_vs(el%name) &
                    //" refers to non-existent notation")
                  call destroy(s_list)
                  exit validLoop
                endif
              enddo
              if (associated(att%default)) then
                s_list = tokenize_to_string_list(str_vs(att%default))
                do k = 1, size(s_list%list)
                  s => s_list%list(k)%s
                  if (.not.notation_exists(fx%nlist, str_vs(s))) then
                    ! Validity Constraint: Notation Attributes
                    call add_error(fx%error_stack, &
                      "Attribute "//str_vs(att%name) &
                      //" of element "//str_vs(el%name) &
                      //" declared as NOTATION refers to non-existent notation")
                    call destroy(s_list)
                    exit validLoop
                  endif
                enddo
                call destroy(s_list)
              endif
            endif
          enddo
        enddo validLoop
      endif
    end subroutine endDTDchecks

    subroutine checkIdRefs
      integer :: i
      character, pointer :: s(:)
      do i = 1, size(idRef_list%list)
        s => idRef_list%list(i)%s
        if (.not.registered_string(id_list, str_vs(s))) then
          call add_error(fx%error_stack, &
            "Reference to undeclared ID "//str_vs(s))
          return
        endif
      enddo
    end subroutine checkIdRefs
  end subroutine sax_parse


  subroutine sax_error(fx, error_handler)
    type(sax_parser_t), intent(inout) :: fx
    optional :: error_handler
    interface
      subroutine error_handler(msg)
        character(len=*), intent(in)     :: msg
      end subroutine error_handler
    end interface

    character, dimension(:), pointer :: errmsg

    integer :: i, m, n, n_err
    n = size(fx%error_stack%stack)
    n_err = n

    do i = 1, n
      n_err = n_err + size(fx%error_stack%stack(i)%msg) ! + spaces + size of entityref
    enddo
    allocate(errmsg(n_err))
    errmsg = ''
    n = 1
    do i = 1, size(fx%error_stack%stack)
      m = size(fx%error_stack%stack(i)%msg)
      errmsg(n:n+m-1) = fx%error_stack%stack(i)%msg
      errmsg(n+m:n+m) = " "
      n = n + m + 1
    enddo
    ! FIXME put location information in here
    if (present(error_handler)) then
      call error_handler(str_vs(errmsg))
      deallocate(errmsg)
      if (fx%state==ST_STOP) return
    else
      call FoX_error(str_vs(errmsg))
    endif

  end subroutine sax_error

  pure function URIlength(fx, qname) result(l_u)
    type(sax_parser_t), intent(in) :: fx
    character(len=*), intent(in) :: qName
    integer :: l_u
    integer :: n
    n = index(QName, ':')
    if (n > 0) then
      l_u = len(getnamespaceURI(fx%nsDict, QName(1:n-1)))
    else
      l_u = len(getnamespaceURI(fx%nsDict))
    endif
  end function URIlength

  pure function getURIofQName(fx, qname) result(URI)
    type(sax_parser_t), intent(in) :: fx
    character(len=*), intent(in) :: qName
    character(len=URIlength(fx, qname)) :: URI

    integer :: n
    n = index(QName, ':')
    if (n > 0) then
      URI = getnamespaceURI(fx%nsDict, QName(1:n-1))
    else
      URI = getnamespaceURI(fx%nsDict)
    endif

  end function getURIofQName

  pure function getLocalNameofQName(qname) result(localName)
    character(len=*), intent(in) :: qName
    character(len=len(QName)-index(QName,':')) :: localName

    localName = QName(index(QName,':')+1:)
  end function getLocalNameofQName
#endif

end module m_sax_parser
