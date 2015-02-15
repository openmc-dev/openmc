include(`foreach.m4')`'dnl
dnl
dnl
define(`TOHW_m_dom_get', `dnl
dnl
dnl Provide a getter for the object:
dnl object is presumed to be a node named np
dnl $1 is attribute type
dnl $2 is attribute name
dnl $3 is location of value within node
dnl $4 is list of applicable node types
dnl
ifelse(`$1', `DOMString', `
  pure function get$2_len(np, p) result(n)
    type(Node), intent(in) :: np
    logical, intent(in) :: p
    integer :: n

ifelse(`$4', `', `dnl
    if (p) then
', `dnl
    if (p .and. ( &
dnl
m4_foreach(`x', `$4', `dnl
      np%nodeType==x .or. &
')`'dnl
      .false.)) then
')`'dnl
      n = size(`$3')
    else
      n = 0
    endif
  end function get$2_len
')dnl
dnl
dnl
TOHW_function(`get$2', (np), c)
    type(Node), pointer :: np
ifelse(`$1', `DOMString', `dnl
#ifdef RESTRICTED_ASSOCIATED_BUG
    character(len=get$2_len(np, .true.)) :: c
#else
    character(len=get$2_len(np, associated(np))) :: c
#endif
', `$1', `Node',`dnl
    type(Node), pointer :: c
', `$1', `NodeList',`dnl
    type(NodeList), pointer :: c
', `$1', `NamedNodeMap',`dnl
    type(NamedNodeMap), pointer :: c
', `$1', `DOMConfiguration',`dnl
    type(DOMConfiguration), pointer :: c
',`dnl
    $1 :: c
')

    if (.not.associated(np)) then
      TOHW_m_dom_throw_error(FoX_NODE_IS_NULL)
    endif

ifelse($4, `', `', `dnl
   if (dnl
m4_foreach(`x', `$4', `getNodeType(np)/=x .and. &
')dnl
      .true.) then
      TOHW_m_dom_throw_error(FoX_INVALID_NODE)
    endif
')`'dnl

ifelse(`$1', `DOMString', `dnl
    c = str_vs($3)
', `$1', `Node', `dnl
    c => $3
', `$1', `NodeList', `dnl
    c => $3
', `$1', `NamedNodeMap', `dnl
    c => $3
', `$1', `DOMConfiguration', `dnl
    c => $3
',`dnl
    c = $3
')
  end function get$2
')dnl
dnl
define(`TOHW_m_dom_set', `dnl
dnl
dnl Provide a setter for the object:
dnl object is presumed to be a node named np
dnl $1 is attribute type
dnl $2 is attribute name
dnl $3 is location of value within node
dnl $4 is list of applicable node types
dnl
dnl
dnl
TOHW_subroutine(`set$2', (np, c))
    type(Node), pointer :: np
ifelse(`$1', `DOMString', `dnl
    character(len=*) :: c
', `$1', `Node',`dnl
    type(Node), pointer :: c
', `$1', `NodeList',`dnl
    type(NodeList), pointer :: c
', `$1', `NamedNodeMap',`dnl
    type(NamedNodeMap), pointer :: c
', `$1', `DOMConfiguration',`dnl
    type(DOMConfiguration), pointer :: c
',`dnl
    $1 :: c
')

    if (.not.associated(np)) then
      TOHW_m_dom_throw_error(FoX_NODE_IS_NULL)
    endif

ifelse($4, `', `', `dnl
   if (dnl
m4_foreach(`x', `$4', `getNodeType(np)/=x .and. &
')dnl
      .true.) then
      TOHW_m_dom_throw_error(FoX_INVALID_NODE)
    endif
')`'dnl

ifelse(`$1', `DOMString', `dnl
    if (associated($3)) deallocate($3)
    $3 => vs_str_alloc(c)
', `$1', `Node', `dnl
    $3 => c
', `$1', `NodeList', `dnl
    $3 => c
', `$1', `NamedNodeMap', `dnl
    $3 => c
', `$1', `DOMConfiguration', `dnl
    $3 => c
',`dnl
    $3 = c
')
  end subroutine set$2
')`'dnl
