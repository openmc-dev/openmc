TOHW_m_dom_publics(`
  
  public :: getTagName
  public :: getAttribute
  public :: setAttribute
  public :: removeAttribute
  public :: getAttributeNode
  public :: setAttributeNode
  public :: removeAttributeNode
  public :: getAttributeNS
  public :: setAttributeNS
  public :: removeAttributeNS
  public :: getAttributeNodeNS
  public :: setAttributeNodeNS
  public :: removeAttributeNodeNS
  public :: hasAttribute
  public :: hasAttributeNS
  public :: setIdAttribute
  public :: setIdAttributeNS
  public :: setIdAttributeNode

')`'dnl
dnl
TOHW_m_dom_contents(`

TOHW_m_dom_get(DOMString, tagName, np%nodeName, (ELEMENT_NODE))

  pure function getAttribute_len(arg, p, name) result(n)
    type(Node), intent(in) :: arg
    logical, intent(in) :: p
    character(len=*), intent(in) :: name
    integer :: n
    
    integer :: i
    
    n = 0
    if (.not.p) return
    if (arg%nodeType/=ELEMENT_NODE) return

    do i = 1, arg%elExtras%attributes%length
      if (str_vs(arg%elExtras%attributes%nodes(i)%this%nodeName)==name) then
        n = getTextContent_len(arg%elExtras%attributes%nodes(i)%this, .true.)
        exit
      endif
    enddo

  end function getAttribute_len

  TOHW_function(getAttribute, (arg, name), c)
    type(Node), pointer :: arg
    character(len=*), intent(in) :: name
#ifdef RESTRICTED_ASSOCIATED_BUG
    character(len=getAttribute_len(arg, .true., name)) :: c
#else
    character(len=getAttribute_len(arg, associated(arg), name)) :: c
#endif

    integer :: i

    if (.not.associated(arg)) then
      TOHW_m_dom_throw_error(FoX_NODE_IS_NULL)
    endif

    if (getNodeType(arg) /= ELEMENT_NODE) then
      TOHW_m_dom_throw_error(FoX_INVALID_NODE)
    endif

    if (len(c)>0) then
      do i = 1, arg%elExtras%attributes%length
        if (str_vs(arg%elExtras%attributes%nodes(i)%this%nodeName)==name) then
          c = getTextContent(arg%elExtras%attributes%nodes(i)%this)
          exit
        endif
      enddo
    else
      c = ""
    endif
        
  end function getAttribute


  TOHW_subroutine(setAttribute, (arg, name, value))
    type(Node), pointer :: arg
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: value

    type(Node), pointer :: nn, dummy
    logical :: quickFix

    if (.not.associated(arg)) then
      TOHW_m_dom_throw_error(FoX_NODE_IS_NULL)
    endif

    if (getNodetype(arg)/=ELEMENT_NODE) then
      TOHW_m_dom_throw_error(FoX_INVALID_NODE)
    elseif (arg%readonly) then
      TOHW_m_dom_throw_error(NO_MODIFICATION_ALLOWED_ERR)
    elseif (.not.checkName(name, getXmlVersionEnum(getOwnerDocument(arg)))) then
      TOHW_m_dom_throw_error(INVALID_CHARACTER_ERR)
    elseif (.not.checkChars(value, getXmlVersionEnum(getOwnerDocument(arg)))) then
      TOHW_m_dom_throw_error(FoX_INVALID_CHARACTER)
    endif

    quickFix = getGCstate(getOwnerDocument(arg)) &
      .and. arg%inDocument

    if (quickFix) call setGCstate(getOwnerDocument(arg), .false.)
    ! then the created attribute is going straight into the document,
    ! so dont faff with hanging-node lists.

    nn => createAttribute(arg%ownerDocument, name)
    call setValue(nn, value)
    dummy => setNamedItem(getAttributes(arg), nn)
    if (associated(dummy)) then
      if (getGCstate(getOwnerDocument(arg)).and..not.dummy%inDocument) &
        call putNodesInDocument(getOwnerDocument(arg), dummy) 
      ! ... so that dummy & children are removed from hangingNodes list.
      call destroyAllNodesRecursively(dummy)
    endif

    if (quickFix) call setGCstate(getOwnerDocument(arg), .true.)

  end subroutine setAttribute


  TOHW_subroutine(removeAttribute, (arg, name))
    type(Node), pointer :: arg
    character(len=*), intent(in) :: name

    type(DOMException) :: ex2
    type(Node), pointer :: dummy
    integer :: e

    if (.not.associated(arg)) then
      TOHW_m_dom_throw_error(FoX_NODE_IS_NULL)
    endif

    if (getNodetype(arg)/=ELEMENT_NODE) then
      TOHW_m_dom_throw_error(FoX_INVALID_NODE)
    elseif (arg%readonly) then
      TOHW_m_dom_throw_error(NO_MODIFICATION_ALLOWED_ERR)
    endif
    
    if (arg%inDocument) &
      call setGCstate(getOwnerDocument(arg), .false.)

    dummy => removeNamedItem(getAttributes(arg), name, ex2)
    ! removeNamedItem took care of any default attributes
    if (inException(ex2)) then
      e = getExceptionCode(ex2)
      if (e/=NOT_FOUND_ERR) then
        TOHW_m_dom_throw_error(e)
      endif
    else
      if (.not.arg%inDocument) then
        ! dummy was not in the doc, so was on hangingNode list.
        ! To remove it from the list:
        call putNodesInDocument(arg%ownerDocument, dummy)
      endif
      call destroyAllNodesRecursively(dummy)
    endif
      
    if (arg%inDocument) &
      call setGCstate(arg%ownerDocument, .true.)

  end subroutine removeAttribute


  TOHW_function(getAttributeNode, (arg, name), attr)
    type(Node), pointer :: arg
    character(len=*), intent(in) :: name
    type(Node), pointer :: attr

    if (.not.associated(arg)) then
      TOHW_m_dom_throw_error(FoX_NODE_IS_NULL)
    endif

    if (arg%nodeType /= ELEMENT_NODE) then
      TOHW_m_dom_throw_error(FoX_INVALID_NODE)
    endif

    attr => getNamedItem(getAttributes(arg), name)

  end function getAttributeNode
  

  TOHW_function(setAttributeNode, (arg, newattr), attr)
    type(Node), pointer :: arg
    type(Node), pointer :: newattr
    type(Node), pointer :: attr
    type(Node), pointer :: dummy

    if (.not.associated(arg)) then
      TOHW_m_dom_throw_error(FoX_NODE_IS_NULL)
    endif

    if (arg%nodeType /= ELEMENT_NODE) then
      TOHW_m_dom_throw_error(FoX_INVALID_NODE)
    elseif (.not.associated(arg%ownerDocument, newattr%ownerDocument)) then
      TOHW_m_dom_throw_error(WRONG_DOCUMENT_ERR)
    elseif (arg%readonly) then
      TOHW_m_dom_throw_error(NO_MODIFICATION_ALLOWED_ERR)
    endif

    if (associated(getOwnerElement(newattr), arg)) then
      attr => newattr
      return
      ! Nothing to do, this attribute is already in this element
    elseif (associated(getOwnerElement(newattr))) then
      TOHW_m_dom_throw_error(INUSE_ATTRIBUTE_ERR)
    endif

    ! this checks if attribute exists already
    ! It also does any adding/removing of hangingnodes
    ! and sets ownerElement appropriately
    dummy => setNamedItem(getAttributes(arg), newattr, ex)
    attr => dummy

  end function setAttributeNode


  TOHW_function(removeAttributeNode, (arg, oldattr), attr)
    type(Node), pointer :: arg
    type(Node), pointer :: oldattr
    type(Node), pointer :: attr

    if (.not.associated(arg)) then
      TOHW_m_dom_throw_error(FoX_NODE_IS_NULL)
    endif

    if (arg%nodeType /= ELEMENT_NODE) then
      TOHW_m_dom_throw_error(FoX_INVALID_NODE)
    endif

    if (.not.associated(arg, getOwnerElement(oldattr))) then
      TOHW_m_dom_throw_error(NOT_FOUND_ERR)
    endif

    attr => removeNamedItem(getAttributes(arg), &
      getNodeName(oldattr), ex)

  end function removeAttributeNode


!  function getElementsByTagName - see m_dom_document


  pure function getAttributesNS_len(arg, p, localname, namespaceURI) result(n)
    type(Node), intent(in) :: arg
    logical, intent(in) :: p
    character(len=*), intent(in) :: localname
    character(len=*), intent(in) :: namespaceURI
    integer :: n
    
    integer :: i
    
    n = 0
    if (.not.p) return
    if (arg%nodeType/=ELEMENT_NODE) return

    do i = 1, arg%elExtras%attributes%length
      if ((str_vs(arg%elExtras%attributes%nodes(i)%this%elExtras%localName)==localname &
        .and. str_vs(arg%elExtras%attributes%nodes(i)%this%elExtras%namespaceURI)==namespaceURI) &
        .or. (namespaceURI=="".and.str_vs(arg%elExtras%attributes%nodes(i)%this%nodeName)==localname)) then
        n = getTextContent_len(arg%elExtras%attributes%nodes(i)%this, .true.)
        exit
      endif
    enddo

  end function getAttributesNS_len

  TOHW_function(getAttributeNS, (arg, namespaceURI, localName), c)
    type(Node), pointer :: arg
    character(len=*), intent(in) :: namespaceURI
    character(len=*), intent(in) :: localName
#ifdef RESTRICTED_ASSOCIATED_BUG
    character(len=getAttributesNS_len(arg, .true., localname, namespaceURI)) :: c
#else
    character(len=getAttributesNS_len(arg, associated(arg), localname, namespaceURI)) :: c
#endif

    integer :: i

    if (.not.associated(arg)) then
      TOHW_m_dom_throw_error(FoX_NODE_IS_NULL)
    endif

    if (arg%nodeType /= ELEMENT_NODE) then
      TOHW_m_dom_throw_error(FoX_INVALID_NODE)
    endif

    if (len(c)>0) then
      do i = 1, arg%elExtras%attributes%length
        if ((str_vs(arg%elExtras%attributes%nodes(i)%this%elExtras%localName)==localname &
          .and. str_vs(arg%elExtras%attributes%nodes(i)%this%elExtras%namespaceURI)==namespaceURI) &
          .or. (namespaceURI=="".and.str_vs(arg%elExtras%attributes%nodes(i)%this%nodeName)==localname)) then
          c = getTextContent(arg%elExtras%attributes%nodes(i)%this)
          exit
        endif
      enddo
    else
      c = ""
    endif

  end function getAttributeNS


  TOHW_subroutine(setAttributeNS, (arg, namespaceURI, qualifiedname, value))
    type(Node), pointer :: arg
    character(len=*), intent(in) :: namespaceURI
    character(len=*), intent(in) :: qualifiedName
    character(len=*), intent(in) :: value

    type(Node), pointer :: nn, dummy
    logical :: quickfix

    if (.not.associated(arg)) then
      TOHW_m_dom_throw_error(FoX_NODE_IS_NULL)
    endif

    if (arg%nodeType /= ELEMENT_NODE) then
      TOHW_m_dom_throw_error(FoX_INVALID_NODE)
    elseif (arg%readonly) then
      TOHW_m_dom_throw_error(NO_MODIFICATION_ALLOWED_ERR)
    elseif (.not.checkName(qualifiedname, getXmlVersionEnum(getOwnerDocument(arg)))) then
      TOHW_m_dom_throw_error(INVALID_CHARACTER_ERR)
    endif
    if (.not.arg%ownerDocument%docExtras%brokenNS) then
      if (.not.checkQName(qualifiedname, getXmlVersionEnum(getOwnerDocument(arg)))) then
        TOHW_m_dom_throw_error(NAMESPACE_ERR)
      elseif (prefixOfQName(qualifiedName)/="" &
        .and. namespaceURI=="") then
        TOHW_m_dom_throw_error(NAMESPACE_ERR)
      elseif (prefixOfQName(qualifiedName)=="xml" .neqv. & 
        namespaceURI=="http://www.w3.org/XML/1998/namespace") then
        TOHW_m_dom_throw_error(NAMESPACE_ERR)
      elseif (namespaceURI=="http://www.w3.org/2000/xmlns/" .neqv. &
        (qualifiedName=="xmlns" .or. prefixOfQName(qualifiedName)=="xmlns")) then
        TOHW_m_dom_throw_error(NAMESPACE_ERR)
      endif
    endif

! FIXME what if namespace is undeclared? Throw an error *only* if FoX_errors is on, otherwise its taken care of by namespace fixup on serialization

    quickFix = getGCstate(getOwnerDocument(arg)) &
      .and. arg%inDocument

    if (quickFix) call setGCstate(getOwnerDocument(arg), .false.)
    ! then the created attribute is going straight into the document,
    ! so dont faff with hanging-node lists.

    nn => createAttributeNS(arg%ownerDocument, namespaceURI, qualifiedname)
    call setValue(nn, value)
    dummy => setNamedItemNS(getAttributes(arg), nn)

    if (associated(dummy)) then
      if (getGCstate(getOwnerDocument(arg)).and..not.dummy%inDocument) &
        call putNodesInDocument(getOwnerDocument(arg), dummy) 
      ! ... so that dummy & children are removed from hangingNodes list.
      call destroyAllNodesRecursively(dummy)
    endif

    if (quickFix) call setGCstate(getOwnerDocument(arg), .true.)

  end subroutine setAttributeNS


  TOHW_subroutine(removeAttributeNS, (arg, namespaceURI, localName))
    type(Node), pointer :: arg
    character(len=*), intent(in) :: namespaceURI
    character(len=*), intent(in) :: localName

    type(DOMException) :: ex2
    type(Node), pointer :: dummy
    integer :: e

    if (.not.associated(arg)) then
      TOHW_m_dom_throw_error(FoX_NODE_IS_NULL)
    endif

    if (arg%nodeType /= ELEMENT_NODE) then
      TOHW_m_dom_throw_error(FoX_INVALID_NODE)
    elseif (arg%readonly) then
      TOHW_m_dom_throw_error(NO_MODIFICATION_ALLOWED_ERR)
    endif

    if (arg%inDocument) &
      call setGCstate(getOwnerDocument(arg), .false.)
    ! So we dont add the removed nodes to the hanging node list

    dummy => removeNamedItemNS(getAttributes(arg), namespaceURI, localName, ex2)
    ! removeNamedItemNS took care of any default attributes
    if (inException(ex2)) then
      e = getExceptionCode(ex2)
      if (e/=NOT_FOUND_ERR) then
        TOHW_m_dom_throw_error(e)
      endif
    else
      if (.not.arg%inDocument) then
        ! dummy was not in the doc, so was already on hangingNode list.
        ! To remove it from the list:
        call putNodesInDocument(arg%ownerDocument, dummy)
      endif
      call destroyAllNodesRecursively(dummy)
    endif
      
    if (arg%inDocument) &
      call setGCstate(arg%ownerDocument, .true.)

  end subroutine removeAttributeNS


  TOHW_function(getAttributeNodeNS, (arg, namespaceURI, localName), attr)
    type(Node), pointer :: arg
    character(len=*), intent(in) :: namespaceURI
    character(len=*), intent(in) :: localName
    type(Node), pointer :: attr

    if (.not.associated(arg)) then
      TOHW_m_dom_throw_error(FoX_NODE_IS_NULL)
    endif

    if (arg%nodeType /= ELEMENT_NODE) then
      TOHW_m_dom_throw_error(FoX_INVALID_NODE)
    endif

    attr => null()     ! as per specs, if not found
    attr => getNamedItemNS(getAttributes(arg), namespaceURI, localname)
  end function getAttributeNodeNS
  

  TOHW_function(setAttributeNodeNS, (arg, newattr), attr)
    type(Node), pointer :: arg
    type(Node), pointer :: newattr
    type(Node), pointer :: attr
    type(Node), pointer :: dummy

    if (.not.associated(arg)) then
      TOHW_m_dom_throw_error(FoX_NODE_IS_NULL)
    endif

    if (arg%nodeType /= ELEMENT_NODE) then
      TOHW_m_dom_throw_error(FoX_INVALID_NODE)
    elseif (.not.associated(arg%ownerDocument, newattr%ownerDocument)) then
      TOHW_m_dom_throw_error(WRONG_DOCUMENT_ERR)
    elseif (arg%readonly) then
      TOHW_m_dom_throw_error(NO_MODIFICATION_ALLOWED_ERR)
    endif

    if (associated(getOwnerElement(newattr), arg)) then
      attr => newattr
      return
      ! Nothing to do, this attribute is already in this element
    elseif (associated(getOwnerElement(newattr))) then
      TOHW_m_dom_throw_error(INUSE_ATTRIBUTE_ERR)
    endif

    ! this checks if attribute exists already
    ! It also does any adding/removing of hangingnodes
    ! and sets ownerElement appropriately
    dummy => setNamedItemNS(getAttributes(arg), newattr, ex)
    attr => dummy

  end function setAttributeNodeNS


  TOHW_function(removeAttributeNodeNS, (arg, oldattr), attr)
    type(Node), pointer :: arg
    type(Node), pointer :: oldattr
    type(Node), pointer :: attr

    if (.not.associated(arg)) then
      TOHW_m_dom_throw_error(FoX_NODE_IS_NULL)
    endif

    if (arg%nodeType /= ELEMENT_NODE) then
      TOHW_m_dom_throw_error(FoX_INVALID_NODE)
    endif

    if (.not.associated(arg, getOwnerElement(oldattr))) then
      TOHW_m_dom_throw_error(NOT_FOUND_ERR)
    endif

    attr => removeNamedItemNS(getAttributes(arg), &
      getNamespaceURI(oldattr), getLocalName(oldattr), ex)

  end function removeAttributeNodeNS


!  function getElementsByTagNameNS - see m_dom_document


  TOHW_function(hasAttribute, (arg, name), p)
    type(Node), pointer :: arg
    character(len=*), intent(in) :: name
    logical :: p

    integer :: i
    type(Node), pointer :: attr

    if (.not.associated(arg)) then
      TOHW_m_dom_throw_error(FoX_NODE_IS_NULL)
    endif
 
   if (arg%nodeType /= ELEMENT_NODE) then
      TOHW_m_dom_throw_error(FoX_INVALID_NODE)
    endif

    p = .false.
    do i = 0, getLength(getAttributes(arg)) - 1
      attr => item(getAttributes(arg), i)
      if (getNodeName(attr)==name) then
        p = .true.
        exit
      endif
    enddo

  end function hasAttribute


  TOHW_function(hasAttributeNS, (arg, namespaceURI, localName), p)
    type(Node), pointer :: arg
    character(len=*), intent(in) :: namespaceURI
    character(len=*), intent(in) :: localName
    logical :: p

    integer :: i
    type(Node), pointer :: attr

    if (.not.associated(arg)) then
      TOHW_m_dom_throw_error(FoX_NODE_IS_NULL)
    endif
 
   if (arg%nodeType /= ELEMENT_NODE) then
      TOHW_m_dom_throw_error(FoX_INVALID_NODE)
    endif

    p = .false.
    do i = 0, getLength(getAttributes(arg))-1
      attr => item(getAttributes(arg), i)
      if (getNamespaceURI(attr)==namespaceURI &
        .and. getLocalName(attr)==localName) then
        p = .true.
        exit
      endif
    enddo

  end function hasAttributeNS

  TOHW_subroutine(setIdAttribute, (arg, name, isId))
    type(Node), pointer :: arg
    character(len=*), intent(in) :: name
    logical, intent(in) :: isId

    type(Node), pointer :: np

    if (arg%readonly) then
      TOHW_m_dom_throw_error(NO_MODIFICATION_ALLOWED_ERR)
    endif

    np => getAttributeNode(arg, name)
    if (associated(np)) then
      call setIsId(np, isId)
    else
      TOHW_m_dom_throw_error(NOT_FOUND_ERR)
    endif

  end subroutine setIdAttribute

  TOHW_subroutine(setIdAttributeNS, (arg, namespaceURI, localname, isId))
    type(Node), pointer :: arg
    character(len=*), intent(in) :: namespaceURI
    character(len=*), intent(in) :: localName
    logical, intent(in) :: isId

    type(Node), pointer :: np

    if (arg%readonly) then
      TOHW_m_dom_throw_error(NO_MODIFICATION_ALLOWED_ERR)
    endif

    np => getAttributeNodeNS(arg, namespaceURI, localname)
    if (associated(np)) then
      call setIsId(np, isId)
    else
      TOHW_m_dom_throw_error(NOT_FOUND_ERR)
    endif

  end subroutine setIdAttributeNS

  TOHW_subroutine(setIdAttributeNode, (arg, idAttr, isId))
    type(Node), pointer :: arg
    type(Node), pointer :: idAttr
    logical, intent(in) :: isId

    if (arg%readonly) then
      TOHW_m_dom_throw_error(NO_MODIFICATION_ALLOWED_ERR)
    elseif (.not.associated(arg, getOwnerElement(idAttr))) then
      TOHW_m_dom_throw_error(NOT_FOUND_ERR)
    endif

    call setIsId(idAttr, isId)

  end subroutine setIdAttributeNode

')`'dnl
