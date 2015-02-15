TOHW_m_dom_publics(`
  
  public :: getNodeName
  public :: getNodeValue
  public :: setNodeValue
  public :: getNodeType
  public :: getParentNode
  public :: getChildNodes
  public :: getFirstChild
  public :: getLastChild
  public :: getNextSibling
  public :: getPreviousSibling
  public :: getAttributes
  public :: getOwnerDocument
  public :: insertBefore
  public :: replaceChild
  public :: removeChild
  public :: appendChild
  public :: hasChildNodes
  public :: cloneNode  
  public :: normalize
  public :: isSupported
  public :: getNamespaceURI
  public :: getPrefix
  public :: setPrefix
  public :: getLocalName
  public :: hasAttributes
  public :: isEqualNode
  public :: isSameNode
  public :: isDefaultNamespace
  public :: lookupNamespaceURI
  public :: lookupPrefix
  public :: getTextContent
  public :: setTextContent

  public :: getNodePath

  public :: setStringValue
  public :: getStringValue
  public :: setReadonlyNode
  public :: getReadOnly

  public :: getBaseURI

')`'dnl
TOHW_m_dom_contents(`

TOHW_m_dom_get(DOMString, nodeName, np%nodeName)

  pure function getNodeValue_len(np, p) result(n)
    type(Node), intent(in) :: np
    logical, intent(in) :: p
    integer :: n

    n = 0 
    if (.not.p) return

    select case(np%nodeType)
    case (ATTRIBUTE_NODE)
      n = getTextContent_len(np, .true.)
    case (CDATA_SECTION_NODE, COMMENT_NODE, PROCESSING_INSTRUCTION_NODE, TEXT_NODE)
      n = size(np%nodeValue)
    end select

  end function getNodeValue_len

  TOHW_function(getNodeValue, (np), c)
    type(Node), pointer :: np
#ifdef RESTRICTED_ASSOCIATED_BUG
    character(len=getNodeValue_len(np, .true.)) :: c
#else
    character(len=getNodeValue_len(np, associated(np))) :: c
#endif

    if (.not.associated(np)) then
      TOHW_m_dom_throw_error(FoX_NODE_IS_NULL)
    endif

    select case(np%nodeType)
    case (ATTRIBUTE_NODE)
      c = getTextContent(np)
    case (CDATA_SECTION_NODE, COMMENT_NODE, PROCESSING_INSTRUCTION_NODE, TEXT_NODE)
      c = str_vs(np%nodeValue)
    case default
      c = ""
    end select
    
  end function getNodeValue

  TOHW_subroutine(setNodeValue, (arg, nodeValue))
    type(Node), pointer :: arg
    character(len=*) :: nodeValue

    if (.not.associated(arg)) then
      TOHW_m_dom_throw_error(FoX_NODE_IS_NULL)
    endif

    if (associated(getOwnerDocument(arg))) then
      if (.not.checkChars(nodeValue, getXmlVersionEnum(getOwnerDocument(arg)))) then
        TOHW_m_dom_throw_error(FoX_INVALID_CHARACTER)
      endif
    endif ! Otherwise its a document node, and nothing will happen anyway

    select case(arg%nodeType)
    case (ATTRIBUTE_NODE)
      call setValue(arg, nodeValue, ex)
    case (CDATA_SECTION_NODE, COMMENT_NODE, PROCESSING_INSTRUCTION_NODE, TEXT_NODE)
      call setData(arg, nodeValue, ex)
    end select

  end subroutine setNodeValue

TOHW_m_dom_get(integer, nodeType, np%nodeType)

TOHW_m_dom_get(Node, parentNode, np%parentNode)

TOHW_m_dom_get(NodeList, childNodes, np%childNodes)

TOHW_m_dom_get(Node, firstChild, np%firstChild)

TOHW_m_dom_get(Node, lastChild, np%lastChild)

TOHW_m_dom_get(Node, previousSibling, np%previousSibling)

TOHW_m_dom_get(Node, nextSibling, np%nextSibling)

  TOHW_function(getAttributes, (arg), nnm)
    type(Node), pointer :: arg
    type(NamedNodeMap), pointer :: nnm

    if (.not.associated(arg)) then
      TOHW_m_dom_throw_error(FoX_NODE_IS_NULL)
    endif

    if (getNodeType(arg)==ELEMENT_NODE) then
      nnm => arg%elExtras%attributes
    else
      nnm => null()
    endif
  end function getAttributes

  TOHW_function(getOwnerDocument, (arg), np)
    type(Node), pointer :: arg
    type(Node), pointer :: np

    if (.not.associated(arg)) then
      TOHW_m_dom_throw_error(FoX_NODE_IS_NULL)
    endif
    
    if (arg%nodeType==DOCUMENT_NODE) then
      np => null()
    else
      np => arg%ownerDocument
    endif
  end function getOwnerDocument

TOHW_m_dom_set(Node, ownerDocument, np%ownerDocument, (DOCUMENT_NODE))

  TOHW_function(insertBefore, (arg, newChild, refChild), np)
    type(Node), pointer :: arg
    type(Node), pointer :: newChild
    type(Node), pointer :: refChild
    type(Node), pointer :: np

    type(Node), pointer :: testChild, testParent, treeroot, this
    type(ListNode), pointer :: temp_nl(:)
    integer :: i, i2, i_t, i_tree
    logical :: doneChildren, doneAttributes

    if (.not.associated(arg).or..not.associated(newChild)) then
      TOHW_m_dom_throw_error(FoX_NODE_IS_NULL)
    endif

    if (.not.associated(refChild)) then
      np => appendChild(arg, newChild, ex)
      return
    endif

    if (arg%readonly) then
      TOHW_m_dom_throw_error(NO_MODIFICATION_ALLOWED_ERR)
    endif

    testParent => arg
    ! Check if you are allowed to put a newChild nodetype under a arg nodetype
    if (newChild%nodeType==DOCUMENT_FRAGMENT_NODE) then
      do i = 1, newChild%childNodes%length
        testChild => newChild%childNodes%nodes(i)%this
        TOHW_m_dom_hierarchy_test
      enddo
    else
      testChild => newChild
      TOHW_m_dom_hierarchy_test
      ! And then check that newChild is not arg or one of args ancestors
      ! (this would never be true if newChild is a documentFragment)
      testParent => arg
      do while (associated(testParent))
        if (associated(testParent, newChild)) then
          TOHW_m_dom_throw_error(HIERARCHY_REQUEST_ERR)
        endif
        testParent => testParent%parentNode
      enddo
    endif

    if (getNodeType(newChild)/=DOCUMENT_TYPE_NODE.and. &
      .not.(associated(arg%ownerDocument, newChild%ownerDocument) &
        .or.associated(arg, newChild%ownerDocument))) then
      TOHW_m_dom_throw_error(WRONG_DOCUMENT_ERR)
    endif

    if (newChild%nodeType==DOCUMENT_FRAGMENT_NODE &
      .and. newChild%childNodes%length==0) then
      np => newChild
      return
      ! Nothing to do
    endif
    if (associated(getParentNode(newChild))) then
      np => removeChild(getParentNode(newChild), newChild, ex) 
      newChild => np
    endif
    
    if (arg%childNodes%length==0) then
      TOHW_m_dom_throw_error(NOT_FOUND_ERR)
    elseif (newChild%nodeType==DOCUMENT_FRAGMENT_NODE) then
      allocate(temp_nl(arg%childNodes%length+newChild%childNodes%length))
    else
      allocate(temp_nl(arg%childNodes%length+1))
    endif

    i_t = 0
    np => null()
    do i = 1, arg%childNodes%length
      if (associated(arg%childNodes%nodes(i)%this, refChild)) then
        np => refChild
        if (newChild%nodeType==DOCUMENT_FRAGMENT_NODE) then
          do i2 = 1, newChild%childNodes%length
            i_t = i_t + 1
            temp_nl(i_t)%this => newChild%childNodes%nodes(i2)%this
            temp_nl(i_t)%this%parentNode => arg
!            call namespaceFixup(temp_nl(i_t)%this)
          enddo
        else
          i_t = i_t + 1
          temp_nl(i_t)%this => newChild
          temp_nl(i_t)%this%parentNode => arg
!          call namespaceFixup(temp_nl(i_t)%this)
        endif
        if (i==1) then
          arg%firstChild => temp_nl(1)%this
          !temp_nl(1)%this%previousSibling => null() ! This is a no-op
        else 
          temp_nl(i-1)%this%nextSibling => temp_nl(i)%this
          temp_nl(i)%this%previousSibling => temp_nl(i-1)%this
        endif
        arg%childNodes%nodes(i)%this%previousSibling => temp_nl(i_t)%this
        temp_nl(i_t)%this%nextSibling => arg%childNodes%nodes(i)%this
      endif
      i_t = i_t + 1
      temp_nl(i_t)%this => arg%childNodes%nodes(i)%this
    enddo

    if (.not.associated(np)) then
      TOHW_m_dom_throw_error(NOT_FOUND_ERR, (temp_nl))
    endif

    np => newChild
    if (getGCstate(arg%ownerDocument)) then
      if (arg%inDocument) then
        if (newChild%nodeType==DOCUMENT_FRAGMENT_NODE) then
          do i = 1, newChild%childNodes%length
            call putNodesInDocument(arg%ownerDocument, newChild%childNodes%nodes(i)%this)
          enddo
        else
          call putNodesInDocument(arg%ownerDocument, newChild)
        endif
        ! If newChild was originally in document, it was removed above so must be re-added
        ! Ideally we would avoid the cost of removal & readding to hanging nodelist
      endif
      ! If arg was not in the document, then newChildren were either 
      ! a) removed above in call to removeChild or
      ! b) in a document fragment and therefore not part of doc either
    endif


    if (getNodeType(newChild)==DOCUMENT_FRAGMENT_NODE) then
      deallocate(newChild%childNodes%nodes)
      allocate(newChild%childNodes%nodes(0))
      newChild%childNodes%length = 0
    endif
    deallocate(arg%childNodes%nodes)
    arg%childNodes%nodes => temp_nl
    arg%childNodes%length = size(arg%childNodes%nodes)

    call updateNodeLists(arg%ownerDocument)

    call updateTextContentLength(arg, newChild%textContentLength)

  end function insertBefore


  TOHW_function(replaceChild, (arg, newChild, oldChild), np)
    type(Node), pointer :: arg
    type(Node), pointer :: newChild
    type(Node), pointer :: oldChild
    type(Node), pointer :: np

    type(Node), pointer :: testChild, testParent, treeroot, this
    type(ListNode), pointer :: temp_nl(:)
    integer :: i, i2, i_t, i_tree
    logical :: doneChildren, doneAttributes

    if (.not.associated(arg).or..not.associated(newChild).or..not.associated(oldChild)) then
      TOHW_m_dom_throw_error(FoX_NODE_IS_NULL)
    endif

    if (arg%readonly) then
      TOHW_m_dom_throw_error(NO_MODIFICATION_ALLOWED_ERR)
    endif

    testParent => arg
    ! Check if you are allowed to put a newChild nodetype under a arg nodetype
    if (newChild%nodeType==DOCUMENT_FRAGMENT_NODE) then
      do i = 1, newChild%childNodes%length
        testChild => newChild%childNodes%nodes(i)%this
        TOHW_m_dom_hierarchy_test
      enddo
    else
      testChild => newChild
      TOHW_m_dom_hierarchy_test
      ! And then check that newChild is not arg or one of args ancestors
      ! (this would never be true if newChild is a documentFragment)
      testParent => arg
      do while (associated(testParent))
        if (associated(testParent, newChild)) then
          TOHW_m_dom_throw_error(HIERARCHY_REQUEST_ERR)
        endif
        testParent => testParent%parentNode
      enddo
    endif

    if (getNodeType(newChild)/=DOCUMENT_TYPE_NODE.and. &
      .not.(associated(arg%ownerDocument, newChild%ownerDocument) &
        .or.associated(arg, newChild%ownerDocument))) then
      TOHW_m_dom_throw_error(WRONG_DOCUMENT_ERR)
    endif

    if (associated(getParentNode(newChild))) &
      newChild => removeChild(getParentNode(newChild), newChild, ex) 

    if (arg%childNodes%length==0) then
      TOHW_m_dom_throw_error(NOT_FOUND_ERR)
    elseif (newChild%nodeType==DOCUMENT_FRAGMENT_NODE) then
      allocate(temp_nl(arg%childNodes%length+newChild%childNodes%length-1))
    else
      temp_nl => arg%childNodes%nodes
    endif

    i_t = 0
    np => null()
    do i = 1, arg%childNodes%length
      if (associated(arg%childNodes%nodes(i)%this, oldChild)) then
        np => oldChild
        if (newChild%nodeType==DOCUMENT_FRAGMENT_NODE) then
          do i2 = 1, newChild%childNodes%length
            i_t = i_t + 1
            temp_nl(i_t)%this => newChild%childNodes%nodes(i2)%this
            temp_nl(i_t)%this%parentNode => arg
!            call namespaceFixup(temp_nl(i_t)%this)
          enddo
        else
          i_t = i_t + 1
          temp_nl(i_t)%this => newChild
          temp_nl(i_t)%this%parentNode => arg
!          call namespaceFixup(temp_nl(i_t)%this)
        endif
        if (i==1) then
          arg%firstChild => temp_nl(1)%this
          !temp_nl(1)%this%previousSibling => null() ! This is a no-op
        else 
          temp_nl(i-1)%this%nextSibling => temp_nl(i)%this
          temp_nl(i)%this%previousSibling => temp_nl(i-1)%this
        endif
        if (i==arg%childNodes%length) then
          arg%lastChild => temp_nl(i_t)%this
          !temp_nl(i_t)%this%nextSibling => null() ! This is a no-op
        else
          arg%childNodes%nodes(i+1)%this%previousSibling => temp_nl(i_t)%this
          temp_nl(i_t)%this%nextSibling => arg%childNodes%nodes(i+1)%this
        endif
      else
        i_t = i_t + 1
        temp_nl(i_t)%this => arg%childNodes%nodes(i)%this
      endif
    enddo

    if (.not.associated(np)) then
      TOHW_m_dom_throw_error(NOT_FOUND_ERR)
    endif
    np%parentNode => null()
    np%previousSibling => null()
    np%nextSibling => null()

!    call namespaceFixup(np)

    if (getGCstate(arg%ownerDocument)) then
      if (arg%inDocument) then
        call removeNodesFromDocument(arg%ownerDocument, oldChild)
        if (newChild%nodeType==DOCUMENT_FRAGMENT_NODE) then
          do i = 1, newChild%childNodes%length
            call putNodesInDocument(arg%ownerDocument, newChild%childNodes%nodes(i)%this)
          enddo
        else
          call putNodesInDocument(arg%ownerDocument, newChild)
        endif
        ! If newChild was originally in document, it was removed above so must be re-added
        ! Ideally we would avoid the cost of removing & re-adding to hangingnodelist
      endif
      ! If arg was not in the document, then newChildren were either 
      ! a) removed above in call to removeChild or
      ! b) in a document fragment and therefore not part of doc either
    endif

    if (newChild%nodeType==DOCUMENT_FRAGMENT_NODE) then
      deallocate(newChild%childNodes%nodes)
      allocate(newChild%childNodes%nodes(0))
      newChild%childNodes%length = 0
      deallocate(arg%childNodes%nodes)
      arg%childNodes%nodes => temp_nl
      arg%childNodes%length = size(arg%childNodes%nodes)
    endif

    call updateNodeLists(arg%ownerDocument)

    call updateTextContentLength(arg, newChild%textContentLength-oldChild%textContentLength)

  end function replaceChild


  TOHW_function(removeChild, (arg, oldChild), np)
    type(Node), pointer :: arg
    type(Node), pointer :: oldChild
    type(Node), pointer :: np

    type(ListNode), pointer :: temp_nl(:)
    integer :: i, i_t

    if (.not.associated(arg).or..not.associated(oldChild)) then
      TOHW_m_dom_throw_error(FoX_NODE_IS_NULL)
    endif

    if (arg%readonly) then
      TOHW_m_dom_throw_error(NO_MODIFICATION_ALLOWED_ERR)
    endif

    allocate(temp_nl(size(arg%childNodes%nodes)-1))
    i_t = 1
    do i = 1, size(arg%childNodes%nodes)
      if (associated(arg%childNodes%nodes(i)%this, oldChild)) then 
        if (associated(arg%firstChild, arg%lastChild)) then
          ! There is only one child, we are removing it.
          arg%firstChild => null()
          arg%lastChild => null()
        elseif (i==1) then
          ! We are removing the first child, but there is a second
          arg%firstChild => arg%childNodes%nodes(2)%this
          arg%childNodes%nodes(2)%this%previousSibling => null()
        elseif (i==size(arg%childNodes%nodes)) then
          ! We are removing the last child, but there is a second-to-last
          arg%lastChild => arg%childNodes%nodes(i-1)%this
          arg%childNodes%nodes(i-1)%this%nextSibling => null()
        else
          ! We are removing a child in the middle
          arg%childNodes%nodes(i-1)%this%nextSibling => arg%childNodes%nodes(i+1)%this
          arg%childNodes%nodes(i+1)%this%previousSibling => arg%childNodes%nodes(i-1)%this
        endif
      else
        if (i_t==size(arg%childNodes%nodes)) exit ! We have failed to find the child
        temp_nl(i_t)%this => arg%childNodes%nodes(i)%this
        i_t = i_t + 1     
      endif
    enddo

    deallocate(arg%childNodes%nodes)
    arg%childNodes%nodes => temp_nl
    arg%childNodes%length = size(temp_nl)
    if (i==i_t) then
      TOHW_m_dom_throw_error(NOT_FOUND_ERR)
    endif
    oldChild%parentNode => null()
    oldChild%previousSibling => null()
    oldChild%nextSibling => null()

!    call namespaceFixup(oldChild)

    if (getGCstate(arg%ownerDocument)) then
      if (arg%inDocument) then
        call removeNodesFromDocument(arg%ownerDocument, oldChild)
      endif
    endif

    np => oldChild

    call updateNodeLists(arg%ownerDocument)

    call updateTextContentLength(arg, -oldChild%textContentLength)

  end function removeChild


  TOHW_function(appendChild, (arg, newChild), np)
    type(Node), pointer :: arg
    type(Node), pointer :: newChild
    type(Node), pointer :: np
    
    type(Node), pointer :: testChild, testParent, treeroot, this
    type(ListNode), pointer :: temp_nl(:)
    integer :: i, i_t, i_tree
    logical :: doneChildren, doneAttributes

    if (.not.associated(arg).or..not.associated(newChild)) then
      TOHW_m_dom_throw_error(FoX_NODE_IS_NULL)
    endif

    if (arg%readonly) then
      TOHW_m_dom_throw_error(NO_MODIFICATION_ALLOWED_ERR)
    endif

    testParent => arg
    ! Check if you are allowed to put a newChild nodetype under a arg nodetype
    if (newChild%nodeType==DOCUMENT_FRAGMENT_NODE) then
      do i = 1, newChild%childNodes%length
        testChild => newChild%childNodes%nodes(i)%this
        TOHW_m_dom_hierarchy_test
      enddo
    else
      testChild => newChild
      TOHW_m_dom_hierarchy_test
      ! And then check that newChild is not arg or one of args ancestors
      ! (this would never be true if newChild is a documentFragment)
      testParent => arg
      do while (associated(testParent))
        if (associated(testParent, newChild)) then
          TOHW_m_dom_throw_error(HIERARCHY_REQUEST_ERR)
        endif
        testParent => testParent%parentNode
      enddo
    endif

    if (getNodeType(newChild)/=DOCUMENT_TYPE_NODE.and. &
      .not.(associated(arg%ownerDocument, newChild%ownerDocument) &
            .or.associated(arg, newChild%ownerDocument))) then
      TOHW_m_dom_throw_error(WRONG_DOCUMENT_ERR)
    endif

    if (newChild%nodeType==DOCUMENT_FRAGMENT_NODE &
      .and. newChild%childNodes%length==0) then
      np => newChild
      return
      ! Nothing to do
    endif

    if (associated(getParentNode(newChild))) &
      newChild => removeChild(getParentNode(newChild), newChild, ex) 

    if (newChild%nodeType==DOCUMENT_FRAGMENT_NODE) then
      allocate(temp_nl(arg%childNodes%length+newChild%childNodes%length))
    else
      allocate(temp_nl(arg%childNodes%length+1))
    endif

    do i = 1, arg%childNodes%length
      temp_nl(i)%this => arg%childNodes%nodes(i)%this
    enddo
    
    if (newChild%nodeType==DOCUMENT_FRAGMENT_NODE) then
      i_t = arg%childNodes%length
      do i = 1, newChild%childNodes%length
        i_t = i_t + 1
        temp_nl(i_t)%this => newChild%childNodes%nodes(i)%this
        if (arg%inDocument) &
          call putNodesInDocument(arg%ownerDocument, temp_nl(i_t)%this)
        temp_nl(i_t)%this%parentNode => arg
!        call namespaceFixup(temp_nl(i_t)%this)
      enddo
      if (arg%childNodes%length==0) then
        arg%firstChild => newChild%firstChild
      else
        newChild%firstChild%previousSibling => arg%lastChild
        arg%lastChild%nextSibling => newChild%firstChild
      endif
      arg%lastChild => newChild%lastChild
      newChild%firstChild => null()
      newChild%lastChild => null()
      deallocate(newChild%childNodes%nodes)
      allocate(newChild%childNodes%nodes(0))
      newChild%childNodes%length = 0
    else
      temp_nl(i)%this => newChild
      if (i==1) then
        arg%firstChild => newChild
        newChild%previousSibling => null()
      else
        temp_nl(i-1)%this%nextSibling => newChild
        newChild%previousSibling => temp_nl(i-1)%this     
      endif
      if (getGCstate(arg%ownerDocument)) then
        if (arg%inDocument.and..not.newChild%inDocument) then
          call putNodesInDocument(arg%ownerDocument, newChild)
        endif
      endif
      newChild%nextSibling => null()
      arg%lastChild => newChild
      newChild%parentNode => arg
!      call namespaceFixup(newChild)
    endif

    deallocate(arg%childNodes%nodes)
    arg%childNodes%nodes => temp_nl
    arg%childNodes%length = size(temp_nl)

    np => newChild

    call updateNodeLists(arg%ownerDocument)

    call updateTextContentLength(arg, newChild%textContentLength)

  end function appendChild


  TOHW_function(hasChildNodes, (arg))
    type(Node), pointer :: arg
    logical :: hasChildNodes
    
    if (.not.associated(arg)) then
      TOHW_m_dom_throw_error(FoX_NODE_IS_NULL)
    endif

    hasChildNodes = associated(arg%firstChild)
    
  end function hasChildNodes

  recursive TOHW_function(cloneNode, (arg, deep), np)
    ! Needs to be recursive in case of entity-references within each other.
    type(Node), pointer :: arg
    logical, intent(in) :: deep
    type(Node), pointer :: np

    type(Node), pointer :: doc, treeroot, thatParent, this, new, ERchild

    logical :: doneAttributes, doneChildren, readonly, brokenNS
    integer :: i_tree

    if (.not.associated(arg)) then
      TOHW_m_dom_throw_error(FoX_NODE_IS_NULL)
    endif

    thatParent => null()
    ERchild => null()
    doc => getOwnerDocument(arg)
    if (.not.associated(doc)) return
    np => null()
    brokenNS = doc%docExtras%brokenNS
    doc%docExtras%brokenNS = .true. ! May need to do stupid NS things
    readonly = .false.

    treeroot => arg
TOHW_m_dom_treewalk(`

      new => null()
      select case(getNodeType(this))
      case (ELEMENT_NODE)
        if (getParameter(getDomConfig(doc), "namespaces")) then
          new => createEmptyElementNS(doc, getNamespaceURI(this), getTagName(this))
        else
          new => createEmptyElement(doc, getTagName(this))
        endif
      case (ATTRIBUTE_NODE)
        if (getParameter(getDomConfig(doc), "namespaces")) then
          new => createAttributeNS(doc, getNamespaceURI(this), getName(this))
        else
          new => createAttribute(doc, getName(this))
        endif
        if (associated(this, arg)) then
          call setSpecified(new, .true.)
        else
          call setSpecified(new, getSpecified(this))
        endif
      case (TEXT_NODE)
        new => createTextNode(doc, getData(this))
      case (CDATA_SECTION_NODE)
        new => createCDataSection(doc, getData(this))
      case (ENTITY_REFERENCE_NODE)
        ERchild => this
        readonly = .true.
        new => createEntityReference(doc, getNodeName(this))
        doneChildren = .true.
      case (ENTITY_NODE)
        return
      case (PROCESSING_INSTRUCTION_NODE)
        new => createProcessingInstruction(doc, getTarget(this), getData(this))
      case (COMMENT_NODE)
        new => createComment(doc, getData(this))
      case (DOCUMENT_NODE)
        return
      case (DOCUMENT_TYPE_NODE)
        return
      case (DOCUMENT_FRAGMENT_NODE)
        new => createDocumentFragment(doc)
      case (NOTATION_NODE)
        return
      end select

      if (.not.associated(thatParent)) then
        thatParent => new
      elseif (associated(new)) then
        if (this%nodeType==ATTRIBUTE_NODE) then
          new => setAttributeNode(thatParent, new)
        else
          new => appendChild(thatParent, new)
        endif
      endif

      if (.not.deep) then
        if (getNodeType(arg)==ATTRIBUTE_NODE.or.getNodeType(arg)==ELEMENT_NODE) then
          continue
        else
          exit
        endif
      endif
', `'`

      if (getNodeType(this)==ENTITY_REFERENCE_NODE &
        .and.associated(ERchild, this)) then
          ERchild => null()
          readonly = .false.
      endif
      this%readonly = readonly
      
', `parentNode')

    np => thatParent
    doc%docExtras%brokenNS = brokenNS

  end function cloneNode

  
  TOHW_function(hasAttributes, (arg))
    type(Node), pointer :: arg
    logical :: hasAttributes

    if (.not.associated(arg)) then
      TOHW_m_dom_throw_error(FoX_NODE_IS_NULL)
    endif
    
    if (arg%nodeType == ELEMENT_NODE) then
      hasAttributes = (getLength(getAttributes(arg)) > 0)
    else
      hasAttributes = .false.
    endif
    
  end function hasAttributes

!  function getBaseURI FIXME

!  function compareDocumentPosition FIXME

  TOHW_subroutine(normalize, (arg))
    type(Node), pointer :: arg
    type(Node), pointer :: this, tempNode, oldNode, treeroot
    integer :: i_tree, i_t
    logical :: doneChildren, doneAttributes
    character, pointer :: temp(:)

    if (.not.associated(arg)) then
      TOHW_m_dom_throw_error(FoX_NODE_IS_NULL)
    endif

! DOM standard requires we ignore readonly status
    treeroot => arg
TOHW_m_dom_treewalk(`

      if (getNodeType(this)==TEXT_NODE) then
        if (associated(this, arg)) exit ! If we are called on a text node itself, then do nothing.
        i_t = getLength(this)
        tempNode => getNextSibling(this)
        do while (associated(tempNode))
          if (getNodeType(tempNode)/=TEXT_NODE) exit
          i_t = i_t + getLength(tempNode)
          tempNode => getNextSibling(tempNode)
        enddo
        if (.not.associated(tempNode, getNextSibling(this))) then
          allocate(temp(i_t))
          temp(:getLength(this)) = vs_str(getData(this))
          i_t = getLength(this)
          tempNode => getNextSibling(this)
          do while (associated(tempNode))
            if (getNodeType(tempNode)/=TEXT_NODE) exit
            temp(i_t+1:i_t+getLength(tempNode)) = vs_str(getData(tempNode))
            i_t = i_t + getLength(tempNode)
            oldNode => tempNode
            tempNode => getNextSibling(tempNode)
            oldNode => removeChild(getParentNode(oldNode), oldNode)
            call remove_node_nl(arg%ownerDocument%docExtras%hangingNodes, oldNode)
            call destroy(oldNode)
          enddo
          deallocate(this%nodeValue)
          this%nodeValue => temp
        endif
      end if
',`')


  end subroutine normalize

  TOHW_function(isSupported, (arg, feature, version), p)
    type(Node), pointer :: arg
    character(len=*), intent(in) :: feature
    character(len=*), intent(in) :: version
    logical :: p

    if (.not.associated(arg)) then
      TOHW_m_dom_throw_error(FoX_NODE_IS_NULL)
    endif

    p = hasFeature(getImplementation(arg%ownerDocument), feature, version)
  end function isSupported

  pure function getNamespaceURI_len(arg, p) result(n)
    type(Node), intent(in) :: arg
    logical, intent(in) :: p
    integer :: n

    n = 0
    if (p) then
      if (arg%nodeType==ELEMENT_NODE &
        .or. arg%nodeType==ATTRIBUTE_NODE &
        .or. arg%nodeType==XPATH_NAMESPACE_NODE) then
        n = size(arg%elExtras%namespaceURI)
      endif
    endif

  end function getNamespaceURI_len

  TOHW_function(getNamespaceURI, (arg), c)
    type(Node), pointer :: arg
#ifdef RESTRICTED_ASSOCIATED_BUG
    character(len=getNamespaceURI_len(arg, .true.)) :: c
#else
    character(len=getNamespaceURI_len(arg, associated(arg))) :: c
#endif

    if (.not.associated(arg)) then
      TOHW_m_dom_throw_error(FoX_NODE_IS_NULL)
    endif

    c = ""
    if (arg%nodeType==ELEMENT_NODE &
      .or. arg%nodeType==ATTRIBUTE_NODE &
      .or. arg%nodeType==XPATH_NAMESPACE_NODE) then
      c = str_vs(arg%elExtras%namespaceURI)
    endif
  end function getNamespaceURI

TOHW_m_dom_set(DOMString, namespaceURI, np%elExtras%namespaceURI, (XPATH_NAMESPACE_NODE))

  pure function getPrefix_len(arg, p) result(n)
    type(Node), intent(in) :: arg
    logical, intent(in) :: p
    integer :: n

    n = 0
    if (p) then
      if (arg%nodeType==ELEMENT_NODE &
        .or. arg%nodeType==ATTRIBUTE_NODE &
        .or. arg%nodeType==XPATH_NAMESPACE_NODE) then
        n = size(arg%elExtras%prefix)
      endif
    endif

  end function getPrefix_len

  TOHW_function(getPrefix, (arg), c)
    type(Node), pointer :: arg
#ifdef RESTRICTED_ASSOCIATED_BUG
    character(len=getPrefix_len(arg, .true.)) :: c
#else
    character(len=getPrefix_len(arg, associated(arg))) :: c
#endif

    if (.not.associated(arg)) then
      TOHW_m_dom_throw_error(FoX_NODE_IS_NULL)
    endif

    c = ""
    if (arg%nodeType==ELEMENT_NODE &
      .or. arg%nodeType==ATTRIBUTE_NODE &
      .or. arg%nodeType==XPATH_NAMESPACE_NODE) then
      c = str_vs(arg%elExtras%prefix)
    endif

  end function getPrefix
  
  TOHW_subroutine(setPrefix, (arg, prefix))
    type(Node), pointer :: arg
    character(len=*) :: prefix

    character, pointer :: tmp(:)
    integer :: i

    if (.not.associated(arg)) then
      TOHW_m_dom_throw_error(FoX_NODE_IS_NULL)
    endif

    if (arg%nodeType==ELEMENT_NODE &
      .or. arg%nodeType==ATTRIBUTE_NODE &
      .or. arg%nodeType==XPATH_NAMESPACE_NODE) then
      if (arg%readonly) then
        TOHW_m_dom_throw_error(NO_MODIFICATION_ALLOWED_ERR)
      elseif (.not.checkName(prefix, getXmlVersionEnum(getOwnerDocument(arg)))) then
        TOHW_m_dom_throw_error(INVALID_CHARACTER_ERR)
      elseif (.not.checkNCName(prefix, getXmlVersionEnum(getOwnerDocument(arg)))) then
        TOHW_m_dom_throw_error(NAMESPACE_ERR)
      elseif (size(arg%elExtras%namespaceURI)==0) then
        TOHW_m_dom_throw_error(NAMESPACE_ERR)
      elseif (prefix=="xml" .and. &
        str_vs(arg%elExtras%namespaceURI)/="http://www.w3.org/XML/1998/namespace") then
        TOHW_m_dom_throw_error(NAMESPACE_ERR)
      elseif (prefix=="xmlns" .and. (getNodeType(arg)/=ATTRIBUTE_NODE &
        .or. str_vs(arg%elExtras%namespaceURI)/="http://www.w3.org/2000/xmlns/")) then
        TOHW_m_dom_throw_error(NAMESPACE_ERR)
      elseif (getNodeType(arg)==ATTRIBUTE_NODE.and.getName(arg)=="xmlns") then
        TOHW_m_dom_throw_error(NAMESPACE_ERR)
      endif
! FIXME check if prefix is declared and already points to same namespace
! but only if we ever get full error-checking up and running.
      deallocate(arg%elExtras%prefix)
      arg%elExtras%prefix => vs_str_alloc(prefix)
      tmp => arg%nodeName
      i = index(str_vs(arg%nodeName), ":")
      if (i==0) then
        arg%nodeName => vs_str_alloc(prefix//":"//str_vs(tmp))
      else
        arg%nodeName => vs_str_alloc(prefix//str_vs(tmp(i:)))
      endif
      deallocate(tmp)
    endif

    call updateNodeLists(arg%ownerDocument)

  end subroutine setPrefix

  pure function getLocalName_len(arg, p) result(n)
    type(Node), intent(in) :: arg
    logical, intent(in) :: p
    integer :: n

    n = 0
    if (p) then
      if (arg%nodeType==ELEMENT_NODE &
        .or. arg%nodeType==ATTRIBUTE_NODE &
        .or. arg%nodeType==XPATH_NAMESPACE_NODE) then
        n = size(arg%elExtras%localName)
      endif
    endif

  end function getLocalName_len

  TOHW_function(getLocalName, (arg), c)
    type(Node), pointer :: arg
#ifdef RESTRICTED_ASSOCIATED_BUG
    character(len=getLocalName_len(arg, .true.)) :: c
#else
    character(len=getLocalName_len(arg, associated(arg))) :: c
#endif

    if (.not.associated(arg)) then
      TOHW_m_dom_throw_error(FoX_NODE_IS_NULL)
    endif

    c = ""
    if (arg%nodeType==ELEMENT_NODE &
      .or. arg%nodeType==ATTRIBUTE_NODE &
      .or. arg%nodeType==XPATH_NAMESPACE_NODE) then
      c = str_vs(arg%elExtras%localName)
    endif

  end function getLocalName

  recursive TOHW_function(isEqualNode, (arg, other), p)
    ! We only have one level of recursion, in case of element attributes
    type(Node), pointer :: arg
    type(Node), pointer :: other
    logical :: p

    type(Node), pointer :: this, that, treeroot, treeroot2, att1, att2
    type(NodeList), pointer :: children1, children2
    type(NamedNodeMap), pointer :: atts1, atts2

    integer :: i_tree, i
    logical :: doneChildren, doneAttributes, equal

    if (.not.associated(arg)) then
      TOHW_m_dom_throw_error(FoX_NODE_IS_NULL)
    endif

    if (isSameNode(arg, other)) then
      ! Shortcut the treewalking
      p = .true.
      return
    else
      p = .false.
    endif

    treeroot => arg
    treeroot2 => other
TOHW_m_dom_treewalk(`

      if (getNodeType(this)/=getNodeType(that)) return
      ! Check necessary equal attributes ...
      if (getNodeName(this)/=getNodeName(that) &
        .or. getLocalName(this)/=getLocalName(that) &
        .or. getNamespaceURI(this)/=getNamespaceURI(that) &
        .or. getPrefix(this)/=getPrefix(that) &
        .or. getNodeValue(this)/=getNodeValue(that)) &
        return
      children1 => getChildNodes(this)
      children2 => getChildNodes(that)
      if (getLength(children1)/=getLength(children2)) return
      ! Well get to the contents of the children later on anyway.
      if (getNodeType(this)==ELEMENT_NODE) then
        ! We must treat attributes specially here (rather than relying on
        ! treewalk) since the order can legitimately change.
        atts1 => getAttributes(this)
        atts2 => getAttributes(that)
        if (getLength(atts1)/=getLength(atts2)) return
        do i = 0, getLength(atts1)-1
          att1 => item(atts1, i)
          if (getNamespaceURI(att1)=="") then
            att2 => getNamedItem(atts2, getNodeName(att1))
          else
            att2 => getNamedItemNS(atts2, getLocalName(att1), getNamespaceURI(att1))
          endif
          if (.not.associated(att2)) return
          if (.not.isEqualNode(att1, att2)) return
        enddo
        doneAttributes = .true.
      elseif (getNodeType(this)==DOCUMENT_TYPE_NODE) then
        if (getPublicId(this)/=getPublicId(that) &
          .or. getSystemId(this)/=getSystemId(that) &
          .or. getInternalSubset(this)/=getInternalSubset(that)) return
        atts1 => getEntities(this)
        atts2 => getEntities(that)
        if (getLength(atts1)/=getLength(atts2)) return
        do i = 0, getLength(atts1)-1
          att1 => item(atts1, i)
          att2 => getNamedItem(atts2, getNodeName(att1))
          if (.not.associated(att2)) return
          if (.not.isEqualNode(att1, att2)) return
        enddo
        atts1 => getNotations(this)
        atts2 => getNotations(that)
        if (getLength(atts1)/=getLength(atts2)) return
        do i = 0, getLength(atts1)-1
          att1 => item(atts1, i)
          att2 => getNamedItem(atts2, getNodeName(att1))
          if (.not.associated(att2)) return
          if (.not.isEqualNode(att1, att2)) return
        enddo
      endif
',`',`double',`')

    p = .true.

  end function isEqualNode


  TOHW_function(isSameNode, (arg, other))
    type(Node), pointer :: arg
    type(Node), pointer :: other
    logical :: isSameNode

    if (.not.associated(arg).or..not.associated(other)) then
      TOHW_m_dom_throw_error(FoX_NODE_IS_NULL)
    endif

    isSameNode = associated(arg, other)

  end function isSameNode

  !FIXME all the lookup* functions below are out of spec,
  ! since they rely on a statically-calculated set of NSnodes
  ! which is only generated at parse time, and updated after
  ! normalize.
  ! the spec reckons it should be dynamic, but because we need
  ! to know string lengths, which must be calculated inside
  ! a pure function, we cant do the recursive walk we need to.
  ! (although isDefaultNamespace could be fixed easily enough)

  TOHW_function(isDefaultNamespace, (np, namespaceURI), p)
    type(Node), pointer :: np
    character(len=*), intent(in) :: namespaceURI
    logical :: p

    type(Node), pointer :: el
    integer :: i

    if (.not.associated(np)) then
      TOHW_m_dom_throw_error(FoX_NODE_IS_NULL)
    endif

    el => null()
    select case(getNodeType(np))
    case (ELEMENT_NODE)
      el => np
    case (ATTRIBUTE_NODE)
      el => getOwnerElement(np)
    case (DOCUMENT_NODE)
      el => getDocumentElement(np)
    end select

    p = .false.
    if (associated(el)) then
      do i = 1, el%elExtras%namespaceNodes%length
        if (size(el%elExtras%namespaceNodes%nodes(i)%this%elExtras%prefix)==0) then
          p = (str_vs(el%elExtras%namespaceNodes%nodes(i)%this%elExtras%namespaceURI)==namespaceURI)
          return
        endif
      enddo
    endif
  end function isDefaultNamespace

  pure function lookupNamespaceURI_len(np, prefix, p) result(n)
    type(Node), intent(in) :: np
    character(len=*), intent(in) :: prefix
    logical, intent(in) :: p
    integer :: n

    integer :: i

    n = 0
    if (.not.p) return
    if (np%nodeType/=ELEMENT_NODE &
      .and. np%nodeType/=ATTRIBUTE_NODE &
      .and. np%nodeType/=DOCUMENT_NODE) return

    if (prefix=="xml".or.prefix=="xmlns") then
      n = 0
      return
    endif

    select case(np%nodeType)
    case (ELEMENT_NODE)
      do i = 1, np%elExtras%namespaceNodes%length
        if (str_vs(np%elExtras%namespaceNodes%nodes(i)%this%elExtras%prefix)==prefix) then
          n = size(np%elExtras%namespaceNodes%nodes(i)%this%elExtras%namespaceURI)
          return
        endif
      enddo
    case (ATTRIBUTE_NODE)
      if (associated(np%elExtras%ownerElement)) then
        do i = 1, np%elExtras%ownerElement%elExtras%namespaceNodes%length
          if (str_vs(np%elExtras%ownerElement%elExtras%namespaceNodes%nodes(i)%this%elExtras%prefix)==prefix) then
            n = size(np%elExtras%ownerElement%elExtras%namespaceNodes%nodes(i)%this%elExtras%namespaceURI)
            return
          endif
        enddo
      endif
    case (DOCUMENT_NODE)
      if (associated(np%docExtras%documentElement)) then
        do i = 1, np%docExtras%documentElement%elExtras%namespaceNodes%length
          if (str_vs(np%docExtras%documentElement%elExtras%namespaceNodes%nodes(i)%this%elExtras%prefix)==prefix) then
            n = size(np%docExtras%documentElement%elExtras%namespaceNodes%nodes(i)%this%elExtras%namespaceURI)
            return
          endif
        enddo
      endif
    end select

  end function lookupNamespaceURI_len

  TOHW_function(lookupNamespaceURI, (np, prefix), c)
    type(Node), pointer :: np
    character(len=*), intent(in) :: prefix
#ifdef RESTRICTED_ASSOCIATED_BUG
    character(len=lookupNamespaceURI_len(np, prefix, .true.)) :: c
#else
    character(len=lookupNamespaceURI_len(np, prefix, associated(np))) :: c
#endif

    type(Node), pointer :: el
    integer :: i

    if (.not.associated(np)) then
      TOHW_m_dom_throw_error(FoX_NODE_IS_NULL)
    endif

    if (len(c)==0) then
      c = ""
      return
    endif

    el => null()
    select case(getNodeType(np))
    case (ELEMENT_NODE)
      el => np
    case (ATTRIBUTE_NODE)
      el => getOwnerElement(np)
    case (DOCUMENT_NODE)
      el => getDocumentElement(np)
    end select

    if (associated(el)) then
      do i = 1, el%elExtras%namespaceNodes%length
        if (str_vs(el%elExtras%namespaceNodes%nodes(i)%this%elExtras%prefix)==prefix) then
          c = str_vs(el%elExtras%namespaceNodes%nodes(i)%this%elExtras%namespaceURI)
          return
        endif
      enddo
    endif

  end function lookupNamespaceURI

  pure function lookupPrefix_len(np, namespaceURI, p) result(n)
    type(Node), intent(in) :: np
    character(len=*), intent(in) :: namespaceURI
    logical, intent(in) :: p
    integer :: n

    integer :: i

    n = 0
    if (.not.p) return
    if (np%nodeType/=ELEMENT_NODE &
      .and. np%nodeType/=ATTRIBUTE_NODE &
      .and. np%nodeType/=DOCUMENT_NODE) return
    
    if (namespaceURI=="" &
      .or. namespaceURI=="http://www.w3.org/XML/1998/namespace" &
      .or. namespaceURI=="http://www.w3.org/2000/xmlns/") then
      return
    endif

    select case(np%nodeType)
    case (ELEMENT_NODE)
      do i = 1, np%elExtras%namespaceNodes%length
        if (str_vs(np%elExtras%namespaceNodes%nodes(i)%this%elExtras%namespaceURI)==namespaceURI) then
          n = size(np%elExtras%namespaceNodes%nodes(i)%this%elExtras%prefix)
          return
        endif
      enddo
    case (ATTRIBUTE_NODE)
      if (associated(np%elExtras%ownerElement)) then
        do i = 1, np%elExtras%ownerElement%elExtras%namespaceNodes%length
          if (str_vs(np%elExtras%ownerElement%elExtras%namespaceNodes%nodes(i)%this%elExtras%namespaceURI)==namespaceURI) then
            n = size(np%elExtras%ownerElement%elExtras%namespaceNodes%nodes(i)%this%elExtras%prefix)
            return
          endif
        enddo
      endif
    case (DOCUMENT_NODE)
      if (associated(np%docExtras%documentElement)) then
        do i = 1, np%docExtras%documentElement%elExtras%namespaceNodes%length
          if (str_vs(np%docExtras%documentElement%elExtras%namespaceNodes%nodes(i)%this%elExtras%namespaceURI)==namespaceURI) then
            n = size(np%docExtras%documentElement%elExtras%namespaceNodes%nodes(i)%this%elExtras%prefix)
            return
          endif
        enddo
      endif
    end select

  end function lookupPrefix_len

  TOHW_function(lookupPrefix, (np, namespaceURI), c)
    type(Node), pointer :: np
    character(len=*), intent(in) :: namespaceURI
#ifdef RESTRICTED_ASSOCIATED_BUG
    character(len=lookupPrefix_len(np, namespaceURI, .true.)) :: c
#else
    character(len=lookupPrefix_len(np, namespaceURI, associated(np))) :: c
#endif

    type(Node), pointer :: el
    integer :: i

    if (.not.associated(np)) then
      TOHW_m_dom_throw_error(FoX_NODE_IS_NULL)
    endif

    if (len(c)==0) then
      c = ""
      return
    endif

    el => null()
    select case(getNodeType(np))
    case (ELEMENT_NODE)
      el => np
    case (ATTRIBUTE_NODE)
      el => getOwnerElement(np)
    case (DOCUMENT_NODE)
      el => getDocumentElement(np)
    end select

    if (associated(el)) then
      do i = 1, el%elExtras%namespaceNodes%length
        if (str_vs(el%elExtras%namespaceNodes%nodes(i)%this%elExtras%namespaceURI)==namespaceURI) then
          c = str_vs(el%elExtras%namespaceNodes%nodes(i)%this%elExtras%prefix)
          return
        endif
      enddo
    endif
  end function lookupPrefix

  ! function getUserData
  ! function setUserData
  ! will not implement ...

  subroutine updateTextContentLength(np, n)
    type(Node), pointer :: np
    integer, intent(in) :: n

    type(Node), pointer :: this

    if (n/=0) then      
      this => np
      do while (associated(this))
        this%textContentLength = this%textContentLength + n
        this => getParentNode(this)
        if (associated(this)) then
          if (getNodeType(this)==DOCUMENT_NODE) exit
        endif
      enddo
    endif
  end subroutine updateTextContentLength

  pure function getTextContent_len(arg, p) result(n)
    type(Node), intent(in) :: arg
    logical, intent(in) :: p
    integer :: n

    if (p) then
      n = arg%textContentLength
    else
      n = 0
    endif
  end function getTextContent_len

  TOHW_function(getTextContent, (arg), c)
    type(Node), pointer :: arg
#ifdef RESTRICTED_ASSOCIATED_BUG
    character(len=getTextContent_len(arg, .true.)) :: c
#else
    character(len=getTextContent_len(arg, associated(arg))) :: c
#endif

    type(Node), pointer :: this, treeroot
    integer :: i, i_tree
    logical :: doneChildren, doneAttributes

    if (.not.associated(arg)) then
      TOHW_m_dom_throw_error(FoX_NODE_IS_NULL)
    endif
    
    if (len(c) == 0) then
      c = ""
      return
    endif

    i = 1
    treeroot => arg
    TOHW_m_dom_treewalk(`
      if (associated(this, treeroot).and.isCharData(getNodeType(this))) then
        c = getData(this)
        return
      endif
      select case(getNodeType(this))
      case (ELEMENT_NODE)
        doneAttributes = .true.
        ! Ignore attributes for text content (unless this is an attribute!)
      case(TEXT_NODE, CDATA_SECTION_NODE)
        if (.not.getIsElementContentWhitespace(this)) then
          c(i:i+size(this%nodeValue)-1) = str_vs(this%nodeValue)
          i = i + size(this%nodeValue)
        endif
      end select
'`')
  end function getTextContent

  TOHW_subroutine(setTextContent, (arg, textContent))
    type(Node), pointer :: arg
    character(len=*), intent(in) :: textContent

    type(Node), pointer :: np
    integer :: i

    if (.not.associated(arg)) then
      TOHW_m_dom_throw_error(FoX_NODE_IS_NULL)
    endif

    if (.not.checkChars(textContent, getXmlVersionEnum(getOwnerDocument(arg)))) then
      TOHW_m_dom_throw_error(FoX_INVALID_CHARACTER)
    endif

    select case(getNodeType(arg))
    case (ELEMENT_NODE, ATTRIBUTE_NODE, DOCUMENT_FRAGMENT_NODE)
      if (arg%readonly) then
        TOHW_m_dom_throw_error(NO_MODIFICATION_ALLOWED_ERR)
      endif
      do i = 1, getLength(getChildNodes(arg))
        call destroyNode(arg%childNodes%nodes(i)%this)
      enddo
      deallocate(arg%childNodes%nodes)
      allocate(arg%childNodes%nodes(0))
      arg%childNodes%length = 0
      arg%firstChild => null()
      arg%lastChild => null()
      arg%textContentLength = 0
      np => createTextNode(getOwnerDocument(arg), textContent)
      np => appendChild(arg, np)
    case (TEXT_NODE, CDATA_SECTION_NODE, PROCESSING_INSTRUCTION_NODE, COMMENT_NODE)
      call setData(arg, textContent)
    case (ENTITY_NODE, ENTITY_REFERENCE_NODE)
      TOHW_m_dom_throw_error(NO_MODIFICATION_ALLOWED_ERR)
    end select
  end subroutine setTextContent

  TOHW_function(getBaseURI, (arg), baseURI)
    type(Node), pointer :: arg
    character(len=200) :: baseURI

    type(Node), pointer :: el
    type(URI), pointer :: URIref, URIbase, newURI

    select case(getNodeType(arg))
    case (ELEMENT_NODE)
      el => arg
    case (ATTRIBUTE_NODE)
      if (getName(arg)=="xml:base") then
        if (associated(getOwnerElement(arg))) then
          el => getParentNode(getOwnerElement(arg))
        else
          el => null()
        endif
      else
        el => getOwnerElement(arg)
      endif
    case (TEXT_NODE)
      ! then are we in an attribute or textContent?
      el => getParentNode(arg)
      do while (associated(el))
        if (getNodeType(el)==ELEMENT_NODE) then
          exit
        elseif (getNodeType(el)==ATTRIBUTE_NODE) then
          el => getOwnerElement(el)
          exit
        else
          el => getParentNode(el)
        endif
      enddo
    case (PROCESSING_INSTRUCTION_NODE)
      ! then are we in or out of element content?
      el => getParentNode(arg)
      do while (associated(el))
        if (getNodeType(el)==ELEMENT_NODE) then
          exit
        elseif (getNodeType(el)==DOCUMENT_NODE) then
          el => getOwnerElement(el)
          exit
        else
          el => getParentNode(el)
        endif
      enddo
    case default
      el => null()
    end select

    URIref => parseURI("")

    do while (associated(el))
      select case (getNodeType(el))
      case (ELEMENT_NODE)
        if (hasAttribute(el, "xml:base")) then
          URIbase => parseURI(getAttribute(el, "xml:base"))
          newURI => rebaseURI(URIbase, URIref)
          call destroyURI(URIbase)
          call destroyURI(URIref)
          URIref => newURI
          if (isAbsoluteURI(URIref)) exit
        endif
      case (ENTITY_REFERENCE_NODE)
        if (getSystemId(el)/="") then
          URIbase => parseURI(getSystemId(el))
          newURI => rebaseURI(URIbase, URIref)
          call destroyURI(URIbase)
          call destroyURI(URIref)
          URIref => newURI
          if (isAbsoluteURI(URIref)) exit
        endif
      case default
        exit
      end select
      el => getParentNode(el) 
    end do

    if (isAbsoluteURI(URIref)) then
      baseURI = expressURI(URIref)
    else
      baseURI = ""
    endif
    call destroyURI(URIref)

  end function getBaseURI

  recursive TOHW_function(getNodePath, (arg), c)
    ! recursive only for atts and text
    type(Node), pointer :: arg
    character(len=100) :: c
    
    type(Node), pointer :: this, this2
    character(len=len(c)) :: c2
    integer :: n

    if (.not.associated(arg)) then
      TOHW_m_dom_throw_error(FoX_NODE_IS_NULL)
    endif

    c = ""
    if (.not.arg%inDocument) return
    select case(getNodeType(arg))

    case (ELEMENT_NODE)
      this => arg
      do while (getNodeType(this)/=DOCUMENT_NODE)
        c2 = ""
        this2 => getPreviousSibling(this)
        n = 0
        do while (associated(this2))
          if (getNodeType(this2)==ELEMENT_NODE &
            .and.getNodeName(this2)==getNodeName(this)) n = n + 1
          this2 => getPreviousSibling(this2)
        enddo
        if (n==0) then
          this2 => getNextSibling(this)
          do while (associated(this2))
            if (getNodeType(this2)==ELEMENT_NODE &
              .and.getNodeName(this2)==getNodeName(this)) then
              n = 1
              exit
            endif
            this2 => getNextSibling(this2)
          enddo
        else
          n = n + 1
        endif
        if (n>0) c2 = "["//n//"]"
        ! What name to use:
        if (getNamespaceURI(this)/="".and.getPrefix(this)=="") then
          ! default namespace; need to do the * trick
          ! how many previous siblings?
          c2 = "/*"//c2
        else
          c2 = "/"//getNodeName(this)//c2
        endif
        c = trim(c2)//c
        this => getParentNode(this)
      enddo

    case (ATTRIBUTE_NODE)
      c = trim(getNodePath(getOwnerElement(arg)))//"/@"//getNodeName(arg)

    case (TEXT_NODE, CDATA_SECTION_NODE)
      ! FIXME this will give wrong answers sometimes if
      ! the tree contains entity references
      this => getParentNode(arg)
      do while (associated(this))
        if (getNodeType(this)==ELEMENT_NODE) exit
        this => getParentNode(this)
      enddo
      if (getNodeType(this)/=ELEMENT_NODE) &
        this => getOwnerElement(this)
      c = trim(getNodePath(this))//"/text()"
      this => getPreviousSibling(arg)
      n = 0
      do while (associated(this))
        if (getNodeType(this)==TEXT_NODE &
          .or.getNodeType(this)==CDATA_SECTION_NODE) n = n + 1
        this => getPreviousSibling(this)
      enddo
      if (n==0) then
        this => getNextSibling(arg)
        do while (associated(this))
          if (getNodeType(this)==COMMENT_NODE &
            .or.getNodeType(this)==CDATA_SECTION_NODE) then
            n = 1
            exit
          endif
          this => getNextSibling(this)
        enddo
      else
        n = n + 1
      endif
      if (n>0) c = trim(c)//"["//n//"]"

    case (PROCESSING_INSTRUCTION_NODE)
      this => getParentNode(arg)
      c = trim(getNodePath(this))//"/processing-instruction("//getNodeName(arg)//")"
      this => getPreviousSibling(arg)
      n = 0
      do while (associated(this))
        if (getNodeType(this)==PROCESSING_INSTRUCTION_NODE &
          .and.getNodeName(this)==getNodeName(arg)) n = n + 1
        this => getPreviousSibling(this)
      enddo
      if (n==0) then
        this => getNextSibling(arg)
        do while (associated(this))
          if (getNodeType(this)==PROCESSING_INSTRUCTION_NODE &
            .and.getNodeName(this)==getNodeName(arg)) then
            n = 1
            exit
          endif
          this => getNextSibling(this)
        enddo
      else
        n = n + 1
      endif
      if (n>0) c = trim(c)//"["//n//"]"

    case (COMMENT_NODE)
      this => getParentNode(arg)
      c = trim(getNodePath(this))//"/comment()"
      this => getPreviousSibling(arg)
      n = 0
      do while (associated(this))
        if (getNodeType(this)==COMMENT_NODE) n = n + 1
        this => getPreviousSibling(this)
      enddo
      if (n==0) then
        this => getNextSibling(arg)
        do while (associated(this))
          if (getNodeType(this)==COMMENT_NODE) then
            n = 1
            exit
          endif
          this => getNextSibling(this)
        enddo
      else
        n = n + 1
      endif
      if (n>0) c = trim(c)//"["//n//"]"

    case (DOCUMENT_NODE)
      c = "/"

    case (XPATH_NAMESPACE_NODE)
      this => getOwnerElement(arg)
      if (getPrefix(arg)=="") then
        c = trim(getNodePath(this))//"/namespace::xmlns"
      else
        c = trim(getNodePath(this))//"/namespace::"//getPrefix(arg)
      endif
      ! FIXME namespace nodes are not marked as inDocument correctly

    end select

  end function getNodePath

  subroutine putNodesInDocument(doc, arg)
    type(Node), pointer :: doc, arg
    type(Node), pointer :: this, treeroot
    logical :: doneChildren, doneAttributes
    integer :: i_tree

    treeroot => arg
TOHW_m_dom_treewalk(`
        this%inDocument = .true.
        call remove_node_nl(doc%docExtras%hangingNodes, this)
',`')


  end subroutine putNodesInDocument

  subroutine removeNodesFromDocument(doc, arg)
    type(Node), pointer :: doc, arg
    type(Node), pointer :: this, treeroot
    logical :: doneChildren, doneAttributes
    integer :: i_tree

    treeroot => arg
TOHW_m_dom_treewalk(`
        this%inDocument = .false.
        call append_nl(doc%docExtras%hangingNodes, this)
',`')

  end subroutine removeNodesFromDocument

  subroutine setReadOnlyNode(arg, p, deep)
    type(Node), pointer :: arg
    logical, intent(in) :: p
    logical, intent(in) :: deep

    type(Node), pointer :: this, treeroot
    integer :: i_tree
    logical :: doneAttributes, doneChildren

    if (deep) then
      treeroot => arg
TOHW_m_dom_treewalk(`
      this%readonly = p
      if (this%nodeType==ELEMENT_NODE) &
        this%elExtras%attributes%readonly = p
',`')
    else
      arg%readonly = p
      if (arg%nodeType==ELEMENT_NODE) &
        arg%elExtras%attributes%readonly = p
    endif

  end subroutine setReadOnlyNode

TOHW_m_dom_get(logical, readonly, np%readonly)

')`'dnl
