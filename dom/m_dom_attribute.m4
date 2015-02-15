TOHW_m_dom_publics(`

  !public :: getName
  public :: getSpecified
  public :: setSpecified
  interface getValue
    module procedure getValue_DOM
  end interface
  public :: getValue
  public :: setValue
  public :: getOwnerElement

  public :: getIsId
  public :: setIsId
  interface getIsId
    module procedure getIsId_DOM
  end interface
  interface setIsId
    module procedure setIsId_DOM
  end interface

')`'dnl
dnl
TOHW_m_dom_contents(`
  
  ! function getName(attribute) result(c) See m_dom_common

! NB All functions manipulating attributes play with the nodelist
! directly rather than through helper functions.
! This is so that getValue_length can be pure,  and the nodeList
! can be explicitly kept up to dat.

TOHW_m_dom_get(logical, specified, np%elExtras%specified, (ATTRIBUTE_NODE))

TOHW_m_dom_set(logical, specified, np%elExtras%specified, (ATTRIBUTE_NODE))

TOHW_m_dom_get(logical, isId_DOM, np%elExtras%isId, (ATTRIBUTE_NODE))

TOHW_m_dom_set(logical, isId_DOM, np%elExtras%isId, (ATTRIBUTE_NODE))

TOHW_m_dom_get(Node, ownerElement, np%elExtras%ownerElement, (ATTRIBUTE_NODE))

  TOHW_function(getValue_DOM, (arg), c)
    type(Node), pointer :: arg
#ifdef RESTRICTED_ASSOCIATED_BUG
    character(len=getTextContent_len(arg, .true.)) :: c 
#else
    character(len=getTextContent_len(arg, associated(arg))) :: c 
#endif

    if (.not.associated(arg)) then
      TOHW_m_dom_throw_error(FoX_NODE_IS_NULL)
    endif

    if (getNodeType(arg)/=ATTRIBUTE_NODE) then
      TOHW_m_dom_throw_error(FoX_INVALID_NODE)
    endif

    c = getTextContent(arg, ex)

  end function getValue_DOM

  TOHW_subroutine(setValue, (arg, value))
    type(Node), pointer :: arg
    character(len=*), intent(in) :: value

    if (.not.associated(arg)) then
      TOHW_m_dom_throw_error(FoX_NODE_IS_NULL)
    endif   

    if (getNodeType(arg)/=ATTRIBUTE_NODE) then
      TOHW_m_dom_throw_error(FoX_INVALID_NODE)
    endif

    call setTextContent(arg, value, ex)

  end subroutine setValue

')`'dnl
