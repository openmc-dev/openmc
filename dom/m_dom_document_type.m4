TOHW_m_dom_publics(`

  !public :: getName
  public :: getEntities
  public :: getNotations
!  public :: getPublicId
!  public :: getSystemId
  public :: getInternalSubset

')`'dnl
dnl
TOHW_m_dom_contents(`

!  function getName(docType) result(c) See m_dom_common

  TOHW_function(getEntities, (arg), nnp)
    type(Node), pointer :: arg
    type(NamedNodeMap), pointer :: nnp

    if (.not.associated(arg)) then
      TOHW_m_dom_throw_error(FoX_NODE_IS_NULL)
    endif

    if (arg%nodeType/=DOCUMENT_TYPE_NODE) then
       TOHW_m_dom_throw_error(FoX_INVALID_NODE)
    endif

    nnp => arg%dtdExtras%entities
  end function getEntities

  TOHW_function(getNotations, (arg), nnp)
    type(Node), pointer :: arg
    type(NamedNodeMap), pointer :: nnp

    if (.not.associated(arg)) then
      TOHW_m_dom_throw_error(FoX_NODE_IS_NULL)
    endif

    if (arg%nodeType/=DOCUMENT_TYPE_NODE) then
       TOHW_m_dom_throw_error(FoX_INVALID_NODE)
    endif

    nnp => arg%dtdExtras%notations
  end function getNotations


!  function getPublicId(docType) result(c) See m_dom_common


!  function getSystemId(docType) result(c) See m_dom_common

  pure function getInternalSubset_len(arg, p) result(n)
    type(Node), pointer :: arg
    logical, intent(in) :: p
    integer :: n

    n = 0
    if (p) then
      if (associated(arg%ownerDocument)) then
        if (associated(arg%ownerDocument%docExtras%xds%intSubset)) then
          n = size(arg%ownerDocument%docExtras%xds%intSubset)
        endif
      endif
    endif
  end function getInternalSubset_len

  TOHW_function(getInternalSubset, (arg), s)
    type(Node), pointer :: arg
#ifdef RESTRICTED_ASSOCIATED_BUG
    character(len=getInternalSubset_len(arg, .true.)) :: s
#else
    character(len=getInternalSubset_len(arg, associated(arg))) :: s
#endif

    if (.not.associated(arg)) then
      TOHW_m_dom_throw_error(FoX_NODE_IS_NULL)
    endif

    if (arg%nodeType/=DOCUMENT_TYPE_NODE) then
       TOHW_m_dom_throw_error(FoX_INVALID_NODE)
    endif

    if (len(s)>0) then
      s = str_vs(arg%ownerDocument%docExtras%xds%intSubset)
    else
      s = ""
    endif
  end function getInternalSubset

')`'dnl
