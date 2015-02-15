module FoX_dom

  use fox_m_fsys_array_str
  use fox_m_fsys_format

  use m_dom_dom
  use m_dom_error
  use m_dom_extras
  use m_dom_parse
  use m_dom_utils

  implicit none
  private

  public :: str_vs, vs_vs_alloc, vs_str_alloc
  public :: str, operator(//)

  public :: DOMImplementation
  public :: Node
  public :: NodeList
  public :: NamedNodeMap

  ! DOM DOMString
  ! no

  ! DOM DOMTimestamp
  !  public :: DOMTimestamp

  ! DOM Exceptions
  public :: DOMException
  public :: inException
  public :: getExceptionCode

  public :: INDEX_SIZE_ERR
  public :: DOMSTRING_SIZE_ERR
  public :: HIERARCHY_REQUEST_ERR
  public :: WRONG_DOCUMENT_ERR
  public :: INVALID_CHARACTER_ERR
  public :: NO_DATA_ALLOWED_ERR
  public :: NO_MODIFICATION_ALLOWED_ERR
  public :: NOT_FOUND_ERR
  public :: NOT_SUPPORTED_ERR
  public :: INUSE_ATTRIBUTE_ERR
  public :: INVALID_STATE_ERR
  public :: SYNTAX_ERR
  public :: INVALID_MODIFICATION_ERR
  public :: NAMESPACE_ERR
  public :: INVALID_ACCESS_ERR
  public :: VALIDATION_ERR
  public :: TYPE_MISMATCH_ERR
  ! XPath
  public :: INVALID_EXPRESSION_ERR
  public :: TYPE_ERR
  ! LS
  public :: PARSE_ERR
  public :: SERIALIZE_ERR

  ! DOM Implementation
  public :: hasFeature
  public :: createDocumentType
  public :: createDocument

  ! DOM Document
  public :: getDocumentElement
  public :: getDocType
  public :: getImplementation
  public :: createDocumentFragment
  public :: createElement
  public :: createTextNode
  public :: createComment
  public :: createCDATASection
  public :: createProcessingInstruction
  public :: createAttribute
  public :: createEntityReference
  public :: getElementsByTagName
  public :: getElementById

  public :: importNode
  public :: createElementNS
  public :: createAttributeNS
  public :: getElementsByTagNameNS

  public :: getXmlStandalone
  public :: setXmlStandalone
  public :: getXmlVersion
  public :: setXmlVersion
  public :: getXmlEncoding
  public :: getInputEncoding
  public :: getDocumentURI
  public :: setDocumentURI
  public :: getStrictErrorChecking
  public :: setStrictErrorChecking
  public :: getDomConfig
  public :: normalizeDocument
  public :: renameNode
  public :: adoptNode

  ! DOM Node
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

  public :: getNodeName
  public :: getNodeValue
  public :: setNodeValue
  public :: getNodeType
  public :: getFirstChild
  public :: getLastChild
  public :: getAttributes
  public :: getNextSibling
  public :: getPreviousSibling
  public :: getParentNode
  public :: getChildNodes
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

  public :: getTextContent
  public :: setTextContent
  public :: isEqualNode
  public :: isSameNode
  public :: isDefaultNamespace
  public :: lookupNamespaceURI
  public :: lookupPrefix

  ! DOM NodeList
  public :: item
  public :: append

  ! DOM NamedNodeMap
  public :: getLength
  public :: getNamedItem
  public :: setNamedItem
  public :: removeNamedItem
!  public :: item
  public :: getNamedItemNS
  public :: setNamedItemNS
  public :: removeNamedItemNS

  ! DOM CharacterData
  ! NB We use the native Fortran string type here
  ! rather than inventing a DOM String, thus no
  ! string type to make public
! public :: getData
! public :: setData
  public :: substringData
  public :: appendData
  public :: insertData
  public :: deleteData
  public :: replaceData

  ! DOM Attr
!  public :: getName
  public :: getSpecified
  public :: getValue
  public :: setValue
  public :: getOwnerElement
  public :: getIsId

  ! DOM Element
  public :: getTagName
  public :: getAttribute
  public :: setAttribute
  public :: removeAttribute
  public :: getAttributeNode
  public :: setAttributeNode
  public :: removeAttributeNode
!  public :: getElementsByTagName
  public :: getAttributeNS
  public :: setAttributeNS
  public :: removeAttributeNS
  public :: getAttributeNodeNS
  public :: setAttributeNodeNS
!  public :: getElementsByTagNameNS
  public :: hasAttribute
  public :: hasAttributeNS
  public :: setIdAttribute
  public :: setIdAttributeNS
  public :: setIdAttributeNode

  !DOM Text
  public :: splitText
  public :: getIsElementContentWhitespace

  !DOM CData
! public :: getData
! public :: setData

  !DOM DocumentType
  public :: getEntities
  public :: getNotations
  public :: getInternalSubset
  
  !DOM Notation

  !DOM Entity
  public :: getNotationName

  !DOM EntityReference

  !DOM ProcessingInstruction
!  public :: getData
!  public :: setData
  public :: getTarget

  !DOM common
  public :: getData
  public :: setData
  public :: getName
  public :: getPublicId
  public :: getSystemId

  !DOM Configuration
  public :: DOMConfiguration
  public :: getParameter
  public :: setParameter
  public :: canSetParameter
  public :: getParameterNames
  
  ! FoX-only interfaces
  public :: newDOMConfig

  public :: getNodePath

  public :: extractDataContent
  public :: extractDataAttribute
  public :: extractDataAttributeNS

  public :: parseFile
  public :: parseString
  public :: serialize

  public :: destroy
  public :: getFoX_checks
  public :: setFoX_checks
  public :: getLiveNodeLists
  public :: setLiveNodeLists
  public :: getNamespaceNodes

end module FoX_dom
