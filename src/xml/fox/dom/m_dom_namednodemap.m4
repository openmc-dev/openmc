TOHW_m_dom_publics(`

  public :: getNamedItem
  public :: setNamedItem
  public :: removeNamedItem
!  public :: item
!  public :: getLength
  public :: getNamedItemNS
  public :: setNamedItemNS
  public :: removeNamedItemNS

!  public :: append
  public :: setReadOnlyMap
  public :: destroyNamedNodeMap


  interface item
    module procedure item_nnm
  end interface

  interface getLength
    module procedure getLength_nnm
  end interface

')`'dnl
dnl
TOHW_m_dom_contents(`

  TOHW_function(getNamedItem, (map, name), np)
    type(NamedNodeMap), pointer :: map
    character(len=*), intent(in) :: name
    type(Node), pointer :: np

    integer :: i

    if (.not.associated(map)) then
      TOHW_m_dom_throw_error(FoX_MAP_IS_NULL)
    endif

    do i = 1, map%length
      if (str_vs(map%nodes(i)%this%nodeName)==name) then
        np => map%nodes(i)%this
        return
      endif
    enddo

    np => null()

  end function getNamedItem


  TOHW_function(setNamedItem, (map, arg), np)
    type(NamedNodeMap), pointer :: map
    type(Node), pointer :: arg
    type(Node), pointer :: np

    integer :: i

    if (.not.associated(map)) then
      TOHW_m_dom_throw_error(FoX_MAP_IS_NULL)
    endif

    if (.not.associated(arg)) then
      TOHW_m_dom_throw_error(FoX_NODE_IS_NULL)
    endif

    if (map%readonly) then
      TOHW_m_dom_throw_error(NO_MODIFICATION_ALLOWED_ERR)
    elseif (map%ownerElement%nodeType==ELEMENT_NODE) then
      if (.not.associated(map%ownerElement%ownerDocument, arg%ownerDocument)) then
        TOHW_m_dom_throw_error(WRONG_DOCUMENT_ERR)
      elseif (getNodeType(arg)/=ATTRIBUTE_NODE) then
        !Additional check from DOM 3
        TOHW_m_dom_throw_error(HIERARCHY_REQUEST_ERR)
      endif
    endif

    if (getNodeType(arg)==ATTRIBUTE_NODE) then
      if (associated(map%ownerElement, getOwnerElement(arg))) then
        ! we are looking at literally the same node
        np => arg
        return
      elseif (associated(getOwnerElement(arg))) then
        TOHW_m_dom_throw_error(INUSE_ATTRIBUTE_ERR)    
      endif
      arg%elExtras%ownerElement => map%ownerElement
    endif

    do i = 0, getLength(map)-1
      np => item(map, i)
      if (getNodeName(np)==getNodeName(arg)) then
        map%nodes(i+1)%this => arg
        exit
      endif
    enddo

    if (i<getLength(map)) then
      if (getGCstate(getOwnerDocument(map%ownerElement)).and.np%inDocument) then
        call removeNodesFromDocument(getOwnerDocument(map%ownerElement), np)
        np%inDocument = .false.
      endif
    else
      !   If not found, insert it at the end of the linked list
      np => null()
      call append_nnm(map, arg)
    endif

    if (map%ownerElement%nodeType==ELEMENT_NODE) then
      if (getGCstate(getOwnerDocument(map%ownerElement))) then
        ! We need to worry about importing this node
        if (map%ownerElement%inDocument) then
          if (.not.arg%inDocument) &
            call putNodesInDocument(getOwnerDocument(map%ownerElement), arg)
        else
          if (arg%inDocument) &
            call removeNodesFromDocument(getOwnerDocument(map%ownerElement), arg)
          endif
      endif
    endif
    ! Otherwise we only ever setNNM when building the doc, so we know this
    ! does not matter

  end function setNamedItem


  TOHW_function(removeNamedItem, (map, name), np)
    type(NamedNodeMap), pointer :: map
    character(len=*), intent(in) :: name
    type(Node), pointer :: np

    type(xml_doc_state), pointer :: xds
    type(element_t), pointer :: elem
    type(attribute_t), pointer :: att
    type(ListNode), pointer :: temp_nl(:)
    integer :: i, i2

    if (.not.associated(map)) then
      TOHW_m_dom_throw_error(FoX_MAP_IS_NULL)
    endif

    if (map%readonly) then
      TOHW_m_dom_throw_error(NO_MODIFICATION_ALLOWED_ERR)
    endif

    do i = 0, map%length-1
      np => item(map, i)
      if (getNodeName(np)==name) then
        xds => getXds(getOwnerDocument(map%ownerElement))
        elem => get_element(xds%element_list, getNodeName(map%ownerElement))
        att => get_attribute_declaration(elem, name)
        if (associated(att)) then
          if (attribute_has_default(att)) then ! there is a default value
            ! Well swap the old one out & put a new one in.
            ! Do *nothing* about namespace handling at this stage,
            ! wait until we are asked for namespace normalization
            if (getParameter( &
              getDomConfig(getOwnerDocument(map%ownerElement)), &
                           "namespaces")) then
              np => createAttributeNS(getOwnerDocument(map%ownerElement), "", name)
            else
              np => createAttribute(getOwnerDocument(map%ownerElement), name)
            endif
            call setValue(np, str_vs(att%default))
            call setSpecified(np, .false.)
            np => setNamedItem(map, np)
            call setSpecified(np, .true.)
            return
          endif
        endif
        ! Otherwise there was no default value, so we just remove the node.
        ! Grab this node
        if (getNodeType(np)==ATTRIBUTE_NODE) np%elExtras%ownerElement => null()
        ! and shrink the node list
        temp_nl => map%nodes
        allocate(map%nodes(size(temp_nl)-1))
        do i2 = 1, i
          map%nodes(i2)%this => temp_nl(i2)%this
        enddo
        do i2 = i + 2, map%length
          map%nodes(i2-1)%this => temp_nl(i2)%this
        enddo
        map%length = size(map%nodes)
        deallocate(temp_nl)
        if (np%inDocument.and.getGCstate(getOwnerDocument(map%ownerElement))) &
          call removeNodesFromDocument(getOwnerDocument(map%ownerElement), np)
        !otherwise we are only going to destroy these nodes anyway,
        ! and finish
        return
      endif
    enddo

    TOHW_m_dom_throw_error(NOT_FOUND_ERR)

  end function removeNamedItem


  TOHW_function(item_nnm, (map, index), np)
    type(NamedNodeMap), pointer :: map
    integer, intent(in) :: index
    type(Node), pointer :: np
    
    if (.not.associated(map)) then
      TOHW_m_dom_throw_error(FoX_MAP_IS_NULL)
    endif

    if (index<0 .or. index>map%length-1) then
      np => null()
    else
      np => map%nodes(index+1)%this
    endif

   end function item_nnm

  TOHW_function(getLength_nnm, (map), n)
    type(namedNodeMap), pointer :: map
    integer :: n

    if (.not.associated(map)) then
       TOHW_m_dom_throw_error(FoX_MAP_IS_NULL)
    endif

    n = map%length
    
  end function getLength_nnm


  TOHW_function(getNamedItemNS, (map, namespaceURI, localName), np)
    type(NamedNodeMap), pointer :: map
    character(len=*), intent(in) :: namespaceURI
    character(len=*), intent(in) :: localName
    type(Node), pointer :: np

    integer :: i

    if (.not.associated(map)) then
      TOHW_m_dom_throw_error(FoX_MAP_IS_NULL)
    elseif (map%ownerElement%nodeType/=ELEMENT_NODE) then
      np => null()
      return
    endif

    do i = 0, getLength(map) - 1
      np => item(map, i)
      if (getNamespaceURI(np)==namespaceURI &
        .and. getLocalName(np)==localName) then
        return
      endif
    enddo
    np => null()

  end function getNamedItemNS


  TOHW_function(setNamedItemNS, (map, arg), np)
    type(NamedNodeMap), pointer :: map
    type(Node), pointer :: arg
    type(Node), pointer :: np

    integer :: i

    if (.not.associated(map)) then
      TOHW_m_dom_throw_error(FoX_MAP_IS_NULL)
    endif

    if (.not.associated(arg)) then
      TOHW_m_dom_throw_error(FoX_NODE_IS_NULL)
    endif

    if (map%readonly) then
      TOHW_m_dom_throw_error(NO_MODIFICATION_ALLOWED_ERR)
    elseif (map%ownerElement%nodeType==ELEMENT_NODE) then
      if (.not.associated(map%ownerElement%ownerDocument, arg%ownerDocument)) then
        TOHW_m_dom_throw_error(WRONG_DOCUMENT_ERR)
      elseif (getNodeType(arg)/=ATTRIBUTE_NODE) then
        !Additional check from DOM 3
        TOHW_m_dom_throw_error(HIERARCHY_REQUEST_ERR)
      endif
    endif

    if (getNodeType(arg)==ATTRIBUTE_NODE) then
      if (associated(map%ownerElement, getOwnerElement(arg))) then
        ! we are looking at literally the same node, so do nothing else
        np => arg
        return
      elseif (associated(getOwnerElement(arg))) then
        TOHW_m_dom_throw_error(INUSE_ATTRIBUTE_ERR)    
      endif
      arg%elExtras%ownerElement => map%ownerElement
    endif

    do i = 0, getLength(map) - 1
      np => item(map, i)
      if ((getLocalName(arg)==getLocalName(np) &
        .and.getNamespaceURI(arg)==getNamespaceURI(np)) &
        ! Additional case to catch adding of specified attributeNS over 
        ! default (NS but unspecified URI) attribute
        .or.(getNamespaceURI(arg)=="".and.getName(arg)==getName(np))) then
        map%nodes(i+1)%this => arg
        exit
      endif
    enddo

    if (i<getLength(map)) then
      if (getGCstate(getOwnerDocument(map%ownerElement))) then
        if (np%inDocument) then
          call removeNodesFromDocument(getOwnerDocument(map%ownerElement), np)
          arg%inDocument = .false.
        endif
      endif
    else
      ! If not found, insert it at the end of the linked list
      np => null()
      call append_nnm(map, arg)
    endif

    if (map%ownerElement%nodeType==ELEMENT_NODE) then
      if (getGCstate(getOwnerDocument(map%ownerElement))) then
        ! We need to worry about importing this node
        if (map%ownerElement%inDocument) then
          if (.not.arg%inDocument) &
            call putNodesInDocument(getOwnerDocument(map%ownerElement), arg)
        else
          if (arg%inDocument) &
            call removeNodesFromDocument(getOwnerDocument(map%ownerElement), arg)
        endif
      endif
    endif

  end function setNamedItemNS


  TOHW_function(removeNamedItemNS, (map, namespaceURI, localName), np)
    type(NamedNodeMap), pointer :: map
    character(len=*), intent(in) :: namespaceURI
    character(len=*), intent(in) :: localName
    type(Node), pointer :: np

    type(xml_doc_state), pointer :: xds
    type(element_t), pointer :: elem
    type(attribute_t), pointer :: att
    type(ListNode), pointer :: temp_nl(:)
    integer :: i, i2

    if (.not.associated(map)) then
      TOHW_m_dom_throw_error(FoX_MAP_IS_NULL)
    endif

    if (map%readonly) then
      TOHW_m_dom_throw_error(NO_MODIFICATION_ALLOWED_ERR)
    endif

    do i = 0, getLength(map) - 1
      np => item(map, i)
      if (getNamespaceURI(np)==namespaceURI &
          .and. getLocalName(np)==localName) then
        ! Grab this node
        xds => getXds(getOwnerDocument(map%ownerElement))
        elem => get_element(xds%element_list, getNodeName(map%ownerElement))
        att => get_attribute_declaration(elem, getName(np))
        if (associated(att)) then
          if (attribute_has_default(att)) then ! there is a default value
            ! Well swap the old one out & put a new one in.
            ! Do *nothing* about namespace handling at this stage,
            ! wait until we are asked for namespace normalization
            np => createAttributeNS(getOwnerDocument(map%ownerElement), getNamespaceURI(np), getName(np))
            call setValue(np, str_vs(att%default))
            call setSpecified(np, .false.)
            np => setNamedItemNS(map, np)
            call setSpecified(np, .true.)
            return
          endif
        endif
        ! Otherwise there was no default value, so we just remove the node.
        ! and shrink the node list
        if (getNodeType(np)==ATTRIBUTE_NODE) np%elExtras%ownerElement => null()
        temp_nl => map%nodes
        allocate(map%nodes(size(temp_nl)-1))
        do i2 = 1, i
          map%nodes(i2)%this => temp_nl(i2)%this
        enddo
        do i2 = i + 2, map%length
          map%nodes(i2-1)%this => temp_nl(i2)%this
        enddo
        map%length = size(map%nodes)
        deallocate(temp_nl)
        if (np%inDocument.and.getGCstate(getOwnerDocument(map%ownerElement))) &
          call removeNodesFromDocument(getOwnerDocument(map%ownerElement), np)
        !otherwise we are only going to destroy these nodes anyway,
        ! and finish
        return
      endif
    enddo

    TOHW_m_dom_throw_error(NOT_FOUND_ERR)

  end function removeNamedItemNS


  subroutine append_nnm(map, arg)
    type(namedNodeMap), pointer :: map
    type(node), pointer :: arg

    type(ListNode), pointer :: temp_nl(:)
    integer :: i

    if (.not.associated(map%nodes)) then
      allocate(map%nodes(1))
      map%nodes(1)%this => arg
      map%length = 1
    else
      temp_nl => map%nodes
      allocate(map%nodes(size(temp_nl)+1))
      do i = 1, size(temp_nl)
        map%nodes(i)%this => temp_nl(i)%this
      enddo
      deallocate(temp_nl)
      map%nodes(size(map%nodes))%this => arg
      map%length = size(map%nodes)
    endif
    if (getNodeType(arg)==ATTRIBUTE_NODE) arg%elExtras%ownerElement => map%ownerElement

   end subroutine append_nnm


  subroutine setReadOnlyMap(map, r)
    type(namedNodeMap), pointer :: map
    logical, intent(in) :: r

    map%readonly = r
  end subroutine setReadOnlyMap

  subroutine destroyNamedNodeMap(map)
    type(namedNodeMap), pointer :: map

    if (associated(map%nodes)) deallocate(map%nodes)
    deallocate(map)
 end subroutine destroyNamedNodeMap

')`'dnl
