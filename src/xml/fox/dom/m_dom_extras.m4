dnl
include(`foreach.m4')dnl
dnl
include(`datatypes.m4')dnl
dnl
include(`common.m4')dnl
dnl
include(`m_dom_exception.m4')dnl
define(`TOHW_extractDataContent', `dnl
dnl
TOHW_subroutine(extractDataContent$1$2, (arg, data, dnl
ifelse(`$1',`Ch', `separator, csv, ')`'dnl
num, iostat))
    type(Node), pointer :: arg
    TOHWM4_declarationtype(`$1')`'dnl
ifelse(dnl
`$2', `Arr', `, dimension(:)',dnl
`$2', `Mat', `, dimension(:,:)')dnl
, intent(out) :: data
ifelse(`$1', `Ch', `
    logical, intent(in), optional :: csv
    character, intent(in), optional :: separator')
    integer, intent(out), optional :: num, iostat
dnl

    if (.not.associated(arg)) then
ifelse(`$1', `Ch', 
`      TOHW_m_dom_throw_error(FoX_NODE_IS_NULL, `', `data = ""')',
`      TOHW_m_dom_throw_error(FoX_NODE_IS_NULL)')
    endif
ifelse(`$1', `Ch', `dnl
    if (present(ex)) then
      call rts(getTextContent(arg, ex), data, separator, csv, num, iostat)
    else
      call rts(getTextContent(arg), data, separator, csv, num, iostat)
    endif
', `
    if (present(ex)) then
      call rts(getTextContent(arg, ex), data, num, iostat)
    else
      call rts(getTextContent(arg), data, num, iostat)
    endif
')
  end subroutine extractDataContent$1$2
')`'dnl
dnl
define(`TOHW_extractDataAttribute', `dnl
dnl
TOHW_subroutine(extractDataAttribute$1$2, (arg, name, data, dnl
ifelse(`$1',`Ch', `separator, csv, ')`'dnl
num, iostat))
    type(Node), pointer :: arg
    character(len=*), intent(in) :: name
ifelse(`$1',`Ch', `
    logical, intent(in), optional :: csv
    character, intent(in), optional :: separator')
    TOHWM4_declarationtype(`$1')`'dnl
ifelse(dnl
`$2', `Arr', `, dimension(:)',dnl
`$2', `Mat', `, dimension(:,:)')dnl
, intent(out) :: data
    integer, intent(out), optional :: num, iostat
dnl
    if (.not.associated(arg)) then
ifelse(`$1', `Ch', 
`      TOHW_m_dom_throw_error(FoX_NODE_IS_NULL, `', `data = ""')',
`      TOHW_m_dom_throw_error(FoX_NODE_IS_NULL)')
    elseif (getNodeType(arg)/=ELEMENT_NODE) then
ifelse(`$1', `Ch', 
`      TOHW_m_dom_throw_error(FoX_INVALID_NODE, `', `data = ""')',
`      TOHW_m_dom_throw_error(FoX_INVALID_NODE)')
    endif

ifelse(`$1', `Ch', `dnl
    if (present(ex)) then
      call rts(getAttribute(arg, name, ex), data, separator, csv, num, iostat)
    else
      call rts(getAttribute(arg, name), data, separator, csv, num, iostat)
    endif
', `
    if (present(ex)) then
      call rts(getAttribute(arg, name, ex), data, num, iostat)
    else
      call rts(getAttribute(arg, name), data, num, iostat)
    endif
')

  end subroutine extractDataAttribute$1$2
')`'dnl
define(`TOHW_extractDataAttributeNS', `dnl
dnl
TOHW_subroutine(extractDataAttNS$1$2, (arg, namespaceURI, localName, data, dnl
ifelse(`$1',`Ch', `separator, csv, ')`'dnl
num, iostat))
    type(Node), pointer :: arg
    character(len=*), intent(in) :: namespaceURI, localName
    TOHWM4_declarationtype(`$1')`'dnl
ifelse(dnl
`$2', `Arr', `, dimension(:)',dnl
`$2', `Mat', `, dimension(:,:)')dnl
, intent(out) :: data
ifelse(`$1',`Ch', `
    logical, intent(in), optional :: csv
    character, intent(in), optional :: separator')
    integer, intent(out), optional :: num, iostat
dnl
    if (.not.associated(arg)) then
ifelse(`$1', `Ch', 
`      TOHW_m_dom_throw_error(FoX_NODE_IS_NULL, `', `data = ""')',
`      TOHW_m_dom_throw_error(FoX_NODE_IS_NULL)')
    elseif (getNodeType(arg)/=ELEMENT_NODE) then
ifelse(`$1', `Ch', 
`      TOHW_m_dom_throw_error(FoX_INVALID_NODE, `', `data = ""')',
`      TOHW_m_dom_throw_error(FoX_INVALID_NODE)')
    endif

ifelse(`$1', `Ch', `dnl
    if (present(ex)) then
      call rts(getAttributeNS(arg, namespaceURI, localName, ex), &
        data, separator, csv, num, iostat)
    else
      call rts(getAttributeNS(arg, namespaceURI, localName), &
        data, separator, csv, num, iostat)
    endif
', `
    if (present(ex)) then
      call rts(getAttributeNS(arg, namespaceURI, localName, ex), &
        data, num, iostat)
    else
      call rts(getAttributeNS(arg, namespaceURI, localName), &
        data, num, iostat)
    endif
')

  end subroutine extractDataAttNS$1$2
')`'dnl
dnl
module m_dom_extras

  use fox_m_fsys_realtypes, only: sp, dp
  use fox_m_fsys_parse_input, only: rts

  use m_dom_error, only: DOMException, inException, throw_exception,           &
    FoX_NODE_IS_NULL, FoX_INVALID_NODE
  use m_dom_dom, only: Node, ELEMENT_NODE,                                     &
    getAttribute, getAttributeNS, getTextContent, getNodeType, getFoX_checks

  implicit none
  private

  public :: extractDataContent
  public :: extractDataAttribute
  public :: extractDataAttributeNS

  interface extractDataContent
m4_foreach(`x', TOHWM4_types, `TOHWM4_interfaceshortlist(`extractDataContent', x)')
  end interface extractDataContent

  interface extractDataAttribute
m4_foreach(`x', TOHWM4_types, `TOHWM4_interfaceshortlist(`extractDataAttribute', x)')
  end interface extractDataAttribute

  interface extractDataAttributeNS
m4_foreach(`x', TOHWM4_types, `TOHWM4_interfaceshortlist(`extractDataAttNS', x)')
  end interface extractDataAttributeNS

contains

m4_foreach(`x', TOHWM4_types, `TOHW_extractDataContent(x, `Sca')
')
m4_foreach(`x', TOHWM4_types, `TOHW_extractDataContent(x, `Arr')
')
m4_foreach(`x', TOHWM4_types, `TOHW_extractDataContent(x, `Mat')
')
m4_foreach(`x', TOHWM4_types, `TOHW_extractDataAttribute(x, `Sca')
')
m4_foreach(`x', TOHWM4_types, `TOHW_extractDataAttribute(x, `Arr')
')
m4_foreach(`x', TOHWM4_types, `TOHW_extractDataAttribute(x, `Mat')
')
m4_foreach(`x', TOHWM4_types, `TOHW_extractDataAttributeNS(x, `Sca')
')
m4_foreach(`x', TOHWM4_types, `TOHW_extractDataAttributeNS(x, `Arr')
')
m4_foreach(`x', TOHWM4_types, `TOHW_extractDataAttributeNS(x, `Mat')
')

end module m_dom_extras
