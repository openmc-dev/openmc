TOHW_m_dom_publics(`

! Assorted functions with identical signatures despite belonging to different types.

  public :: getData
  public :: setData
  public :: getName
  public :: getPublicId
  public :: getSystemId

')`'dnl
dnl
TOHW_m_dom_contents(`

  TOHW_m_dom_get(DOMString, data, np%nodeValue, (TEXT_NODE, COMMENT_NODE, CDATA_SECTION_NODE, PROCESSING_INSTRUCTION_NODE))

  TOHW_subroutine(setData, (arg, data))
    type(Node), pointer :: arg
    character(len=*) :: data
    
    integer :: n

    if (.not.associated(arg)) then
      TOHW_m_dom_throw_error(FoX_NODE_IS_NULL)
    endif
    
!NB special case in order to check readonly correctly
    if (arg%nodeType==TEXT_NODE .or. &
      arg%nodeType==COMMENT_NODE .or. &
      arg%nodeType==CDATA_SECTION_NODE .or. &
      arg%nodeType==PROCESSING_INSTRUCTION_NODE) then
      if (arg%readonly) then
        TOHW_m_dom_throw_error(NO_MODIFICATION_ALLOWED_ERR)
      endif
    else
      TOHW_m_dom_throw_error(FoX_INVALID_NODE)
    endif

    select case (arg%nodeType)
    case (CDATA_SECTION_NODE)
      if (index(data,"]]>")>0) then
        TOHW_m_dom_throw_error(FoX_INVALID_CDATA_SECTION)
      endif
    case (COMMENT_NODE)
      if (index(data,"--")>0) then
        TOHW_m_dom_throw_error(FoX_INVALID_COMMENT)
      endif
    case (PROCESSING_INSTRUCTION_NODE)
      if (index(data,"?>")>0) then
        TOHW_m_dom_throw_error(FoX_INVALID_PI_DATA)
      endif
    end select

    deallocate(arg%nodeValue)
    arg%nodeValue => vs_str_alloc(data)

    if (arg%nodeType==TEXT_NODE .or. &
      arg%nodeType==CDATA_SECTION_NODE) then
      n = len(data) - arg%textContentLength
      call updateTextContentLength(arg, n)
    endif

  end subroutine setData
  
  TOHW_m_dom_get(DOMString, name, np%nodeName, (DOCUMENT_TYPE_NODE, ATTRIBUTE_NODE))

  TOHW_m_dom_get(DOMString, publicId, np%dtdExtras%publicId, (DOCUMENT_TYPE_NODE, NOTATION_NODE, ENTITY_NODE))

  TOHW_m_dom_get(DOMString, systemId, np%dtdExtras%systemId, (DOCUMENT_TYPE_NODE, NOTATION_NODE, ENTITY_NODE))

')`'dnl
