TOHW_m_dom_publics(`
  integer, parameter ::     ELEMENT_NODE                   = 1
  integer, parameter ::     ATTRIBUTE_NODE                 = 2
  integer, parameter ::     TEXT_NODE                      = 3
  integer, parameter ::     CDATA_SECTION_NODE             = 4
  integer, parameter ::     ENTITY_REFERENCE_NODE          = 5
  integer, parameter ::     ENTITY_NODE                    = 6
  integer, parameter ::     PROCESSING_INSTRUCTION_NODE    = 7
  integer, parameter ::     COMMENT_NODE                   = 8
  integer, parameter ::     DOCUMENT_NODE                  = 9
  integer, parameter ::     DOCUMENT_TYPE_NODE             = 10
  integer, parameter ::     DOCUMENT_FRAGMENT_NODE         = 11
  integer, parameter ::     NOTATION_NODE                  = 12
  integer, parameter ::     XPATH_NAMESPACE_NODE           = 13

  type DOMImplementation
    private
    character(len=7) :: id = "FoX_DOM"
    logical :: FoX_checks = .true. ! Do extra checks not mandated by DOM
  end type DOMImplementation

  type ListNode
    private
    type(Node), pointer :: this => null()
  end type ListNode 

  type NodeList
    private
    character, pointer :: nodeName(:) => null() ! What was getByTagName run on?
    character, pointer :: localName(:) => null() ! What was getByTagNameNS run on?
    character, pointer :: namespaceURI(:) => null() ! What was getByTagNameNS run on?
    type(Node), pointer :: element => null() ! which element or document was the getByTagName run from?
    type(ListNode), pointer :: nodes(:) => null()
    integer :: length = 0
  end type NodeList

  type NodeListptr
    private
    type(NodeList), pointer :: this
  end type NodeListptr

  type NamedNodeMap
    private
    logical :: readonly = .false.
    type(Node), pointer :: ownerElement => null()
    type(ListNode), pointer :: nodes(:) => null()
    integer :: length = 0
  end type NamedNodeMap

  type documentExtras
    type(DOMImplementation), pointer :: implementation => null() ! only for doctype
    type(Node), pointer :: docType => null()
    type(Node), pointer :: documentElement => null()
    character, pointer :: inputEncoding(:) => null()
    character, pointer :: xmlEncoding(:) => null()
    type(NodeListPtr), pointer :: nodelists(:) => null() ! document
    ! In order to keep track of all nodes not connected to the document
    logical :: liveNodeLists ! For the document, are nodelists live?
    type(NodeList) :: hangingNodes ! For the document, list of nodes not associated with doc
    type(xml_doc_state), pointer :: xds => null()
    logical :: strictErrorChecking = .true.
    logical :: brokenNS = .false. ! FIXME consolidate these logical variables into bitmask
    type(DOMConfiguration), pointer :: domConfig => null()
  end type documentExtras

  type elementOrAttributeExtras
    ! Needed for all:
    character, pointer, dimension(:) :: namespaceURI => null()
    character, pointer, dimension(:) :: prefix => null()
    character, pointer, dimension(:) :: localName => null()
    ! Needed for elements:
    type(NamedNodeMap) :: attributes
    type(NodeList) :: namespaceNodes
    ! Needed for attributes:
    type(Node), pointer :: ownerElement => null()
    logical :: specified = .true.
    logical :: isId = .false.
    logical :: dom1 = .false.
  end type elementOrAttributeExtras

  type docTypeExtras
    character, pointer :: publicId(:) => null() ! doctype, entity, notation 
    character, pointer :: systemId(:) => null() ! doctype, entity, notation
    character, pointer :: notationName(:) => null() ! entity
    logical :: illFormed = .false. ! entity
    type(namedNodeMap) :: entities ! doctype
    type(namedNodeMap) :: notations ! doctype
  end type docTypeExtras

  type Node
    private
    logical :: readonly = .false.
    character, pointer, dimension(:)         :: nodeName => null()
    character, pointer, dimension(:)         :: nodeValue => null()
    integer             :: nodeType        = 0
    type(Node), pointer :: parentNode      => null()
    type(Node), pointer :: firstChild      => null()
    type(Node), pointer :: lastChild       => null()
    type(Node), pointer :: previousSibling => null()
    type(Node), pointer :: nextSibling     => null()
    type(Node), pointer :: ownerDocument   => null()
    type(NodeList) :: childNodes ! not for text, cdata, PI, comment, notation, docType, XPath
    logical :: inDocument = .false.! For a node, is this node associated to the doc?
    logical :: ignorableWhitespace = .false. ! Text nodes only
    type(documentExtras), pointer :: docExtras => null()
    type(elementOrAttributeExtras), pointer :: elExtras => null()
    type(docTypeExtras), pointer :: dtdExtras => null()
    integer :: textContentLength = 0
  end type Node

  type(DOMImplementation), save, target :: FoX_DOM

  interface destroy
    module procedure destroyNode
    module procedure destroyNodeList
    module procedure destroyNamedNodeMap
    module procedure destroyDOMConfig
  end interface destroy

  public :: ELEMENT_NODE
  public :: ATTRIBUTE_NODE
  public :: TEXT_NODE
  public :: CDATA_SECTION_NODE
  public :: ENTITY_REFERENCE_NODE
  public :: ENTITY_NODE
  public :: PROCESSING_INSTRUCTION_NODE
  public :: COMMENT_NODE
  public :: DOCUMENT_NODE
  public :: DOCUMENT_TYPE_NODE
  public :: DOCUMENT_FRAGMENT_NODE
  public :: NOTATION_NODE

  public :: DOMImplementation
  public :: DOMConfiguration
  public :: Node

  public :: ListNode
  public :: NodeList
  public :: NamedNodeMap

  public :: destroy
  public :: destroyAllNodesRecursively

')`'dnl
dnl
TOHW_m_dom_contents(`

  TOHW_function(createNode, (arg, nodeType, nodeName, nodeValue), np)
    type(Node), pointer :: arg
    integer, intent(in) :: nodeType
    character(len=*), intent(in) :: nodeName
    character(len=*), intent(in) :: nodeValue
    type(Node), pointer :: np

    allocate(np)
    np%ownerDocument => arg
    np%nodeType = nodeType
    np%nodeName => vs_str_alloc(nodeName)
    np%nodeValue => vs_str_alloc(nodeValue)

    allocate(np%childNodes%nodes(0))

  end function createNode

  recursive TOHW_subroutine(destroyNode, (np))
    type(Node), pointer :: np

    if (.not.associated(np)) return

    select case(np%nodeType)
    case (ELEMENT_NODE, ATTRIBUTE_NODE, XPATH_NAMESPACE_NODE)
      call destroyElementOrAttribute(np,ex)
    case (DOCUMENT_TYPE_NODE)
      call destroyDocumentType(np,ex)
    case (ENTITY_NODE, NOTATION_NODE)
      call destroyEntityOrNotation(np,ex)
    case (DOCUMENT_NODE)
      call destroyDocument(np,ex)
    end select
    call destroyNodeContents(np)
    deallocate(np)

  end subroutine destroyNode

  recursive TOHW_subroutine(destroyElementOrAttribute, (np))
    type(Node), pointer :: np

    integer :: i

    if (np%nodeType /= ELEMENT_NODE &
      .and. np%nodeType /= ATTRIBUTE_NODE &
      .and. np%nodeType /= XPATH_NAMESPACE_NODE) then
       TOHW_m_dom_throw_error(FoX_INTERNAL_ERROR)
    endif

    if (associated(np%elExtras%attributes%nodes)) deallocate(np%elExtras%attributes%nodes)
    do i = 1, np%elExtras%namespaceNodes%length
      call destroyNode(np%elExtras%namespaceNodes%nodes(i)%this)
    enddo
    if (associated(np%elExtras%namespaceNodes%nodes)) deallocate(np%elExtras%namespaceNodes%nodes)
    if (associated(np%elExtras%namespaceURI)) deallocate(np%elExtras%namespaceURI)
    if (associated(np%elExtras%prefix)) deallocate(np%elExtras%prefix)
    if (associated(np%elExtras%localName)) deallocate(np%elExtras%localName)
    deallocate(np%elExtras)

  end subroutine destroyElementOrAttribute

  TOHW_subroutine(destroyEntityOrNotation, (np))
    type(Node), pointer :: np

    if (np%nodeType /= ENTITY_NODE &
      .and. np%nodeType /= NOTATION_NODE) then
       TOHW_m_dom_throw_error(FoX_INTERNAL_ERROR)
    endif

    if (associated(np%dtdExtras%publicId)) deallocate(np%dtdExtras%publicId)
    if (associated(np%dtdExtras%systemId)) deallocate(np%dtdExtras%systemId)
    if (associated(np%dtdExtras%notationName)) deallocate(np%dtdExtras%notationName)

    deallocate(np%dtdExtras)

  end subroutine destroyEntityOrNotation

  TOHW_subroutine(destroyDocumentType, (np))
    type(Node), pointer :: np

    integer :: i

    if (np%nodeType /= DOCUMENT_TYPE_NODE) then
       TOHW_m_dom_throw_error(FoX_INTERNAL_ERROR)
    endif

    if (associated(np%dtdExtras%publicId)) deallocate(np%dtdExtras%publicId)
    if (associated(np%dtdExtras%systemId)) deallocate(np%dtdExtras%systemId)

    ! Destroy all entities & notations (docType only)
    if (associated(np%dtdExtras%entities%nodes)) then
      do i = 1, size(np%dtdExtras%entities%nodes)
        call destroyAllNodesRecursively(np%dtdExtras%entities%nodes(i)%this)
      enddo
      deallocate(np%dtdExtras%entities%nodes)
    endif
    if (associated(np%dtdExtras%notations%nodes)) then
      do i = 1, size(np%dtdExtras%notations%nodes)
        call destroy(np%dtdExtras%notations%nodes(i)%this)
      enddo
      deallocate(np%dtdExtras%notations%nodes)
    endif

    deallocate(np%dtdExtras)

  end subroutine destroyDocumentType

  recursive subroutine destroyAllNodesRecursively(arg, except)
    ! Only recurses once into destroyDocumentType
    type(Node), pointer :: arg
    logical, intent(in), optional :: except
    
    type(Node), pointer :: this, deadNode, treeroot
    logical :: doneChildren, doneAttributes
    integer :: i_tree

    if (.not.associated(arg)) return

    treeroot => arg
TOHW_m_dom_treewalk(`',`',`deadNode', `')

    deallocate(arg%childNodes%nodes)
    allocate(arg%childNodes%nodes(0))
    arg%firstChild => null()
    arg%lastChild => null()

    if (.not.present(except)) call destroyNode(arg)

  end subroutine destroyAllNodesRecursively

  subroutine destroyNodeContents(np)
    type(Node), intent(inout) :: np
    
    if (associated(np%nodeName)) deallocate(np%nodeName)
    if (associated(np%nodeValue)) deallocate(np%nodeValue)

    deallocate(np%childNodes%nodes)

  end subroutine destroyNodeContents

')`'dnl
