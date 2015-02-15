TOHW_m_dom_publics(`  
  public :: splitText
  public :: getIsElementContentWhitespace
  public :: setIsElementContentWhitespace
')`'dnl
dnl
TOHW_m_dom_contents(`

  TOHW_function(splitText, (arg, offset), np)
    type(Node), pointer :: arg
    integer, intent(in) :: offset

    type(Node), pointer :: np

    character, pointer :: tmp(:)

    if (.not.associated(arg)) then
      TOHW_m_dom_throw_error(FoX_NODE_IS_NULL)
    endif

    if (.not.(arg%nodeType==TEXT_NODE.or.arg%nodeType==CDATA_SECTION_NODE)) then
      TOHW_m_dom_throw_error(FoX_INVALID_NODE)
    elseif (arg%readonly) then
      TOHW_m_dom_throw_error(NO_MODIFICATION_ALLOWED_ERR)
    elseif (offset<0 .or. offset>size(arg%nodeValue)) then
      TOHW_m_dom_throw_error(INDEX_SIZE_ERR)
    endif

    tmp => arg%nodeValue
    if (arg%nodeType==TEXT_NODE) then
      np => createTextNode(arg%ownerDocument, str_vs(tmp(offset+1:)))
    elseif (arg%nodeType==CDATA_SECTION_NODE) then
      np => createCdataSection(arg%ownerDocument, str_vs(tmp(offset+1:)))
    endif
    arg%nodeValue => vs_str_alloc(str_vs(tmp(:offset)))     
    deallocate(tmp)
    if (associated(arg%parentNode)) then
      if (associated(arg%nextSibling)) then
        np => insertBefore(arg%parentNode, np, arg%nextSibling)
      else
        np => appendChild(arg%parentNode, np)
      endif
    endif

  end function splitText

TOHW_m_dom_get(logical, isElementContentWhitespace, np%ignorableWhitespace, (TEXT_NODE, CDATA_SECTION_NODE))

  TOHW_subroutine(setIsElementContentWhitespace, (np, isElementContentWhitespace))
    type(Node), pointer :: np
    logical :: isElementContentWhitespace

    integer :: n

    np%ignorableWhitespace = isElementContentWhitespace

    if (isElementContentWhitespace) then
      n = -np%textContentLength
    else
      n = size(np%nodeValue)
    endif

    call updateTextContentLength(np, n)
  end subroutine setIsElementContentWhitespace

! function getWholeText
! function replaceWholeText

')`'dnl
