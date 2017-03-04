module pugixml

  use, intrinsic :: ISO_C_BINDING

  implicit none
  private

  interface
    function strlen(str) result(sz) bind(C)
      import C_PTR, C_SIZE_T
      type(C_PTR), value :: str
      integer(C_SIZE_T) :: sz
    end function strlen
  end interface

  !=============================================================================
  ! XMLNode C interface and derived type

  interface
    function xml_node_name(node) result(name) bind(C)
      import C_PTR
      type(C_PTR), value :: node
      type(C_PTR) :: name
    end function xml_node_name

    function xml_node_name_size(base) result(sz) bind(C)
      import C_PTR, C_SIZE_T
      type(C_PTR), value :: base
      integer(C_SIZE_T) :: sz
    end function xml_node_name_size

    function xml_node_child(node, name) result(child) bind(C)
      import C_PTR
      type(C_PTR), value :: node
      type(C_PTR), value :: name
      type(C_PTR) :: child
    end function xml_node_child

    function xml_node_next_sibling(node, name) result(next) bind(C)
      import C_PTR
      type(C_PTR), value :: node
      type(C_PTR), value :: name
      type(C_PTR) :: next
    end function xml_node_next_sibling

    function xml_node_attribute(node, name) result(attribute) bind(C)
      import C_PTR
      type(C_PTR), value :: node
      type(C_PTR), value :: name
      type(C_PTR) :: attribute
    end function xml_node_attribute

    function xml_node_child_value(node) result(val) bind(C)
      import C_PTR
      type(C_PTR), value :: node
      type(C_PTR) :: val
    end function xml_node_child_value

    function xml_node_text(node) result(text) bind(C)
      import C_PTR
      type(C_PTR), value :: node
      type(C_PTR) :: text
    end function xml_node_text
  end interface

  type :: XMLNode
    type(C_PTR) :: ptr
  contains
    procedure :: name => xmlnode_name
    procedure :: child => xmlnode_child
    procedure :: next_sibling => xmlnode_next_sibling
    procedure :: attribute => xmlnode_attribute
    procedure :: child_value => xmlnode_child_value
    procedure :: text => xmlnode_text
    procedure :: associated => xmlnode_associated
  end type XMLNode

  !=============================================================================
  ! XMLAttribute C interface and derived type

  interface
    function xml_attribute_name(attribute) result(name) bind(C)
      import C_PTR
      type(C_PTR), value :: attribute
      type(C_PTR) :: name
    end function xml_attribute_name

    function xml_attribute_name_size(base) result(sz) bind(C)
      import C_PTR, C_SIZE_T
      type(C_PTR), value :: base
      integer(C_SIZE_T) :: sz
    end function xml_attribute_name_size

    function xml_attribute_value(attribute) result(val) bind(C)
      import C_PTR
      type(C_PTR), value :: attribute
      type(C_PTR) :: val
    end function xml_attribute_value

    function xml_attribute_next_attribute(attribute) result(next) bind(C)
      import C_PTR
      type(C_PTR), value :: attribute
      type(C_PTR) :: next
    end function xml_attribute_next_attribute

    function xml_attribute_as_bool(attribute) result(x) bind(C)
      import C_PTR, C_BOOL
      type(C_PTR), value :: attribute
      logical(C_BOOL) :: x
    end function xml_attribute_as_bool

    function xml_attribute_as_int(attribute) result(x) bind(C)
      import C_PTR, C_INT
      type(C_PTR), value :: attribute
      integer(C_INT) :: x
    end function xml_attribute_as_int

    function xml_attribute_as_llong(attribute) result(x) bind(C)
      import C_PTR, C_LONG_LONG
      type(C_PTR), value :: attribute
      integer(C_LONG_LONG) :: x
    end function xml_attribute_as_llong

    function xml_attribute_as_double(attribute) result(x) bind(C)
      import C_PTR, C_DOUBLE
      type(C_PTR), value :: attribute
      real(C_DOUBLE) :: x
    end function xml_attribute_as_double
  end interface

  type :: XMLAttribute
    type(C_PTR) :: ptr
  contains
    procedure :: name => xmlattribute_name
    procedure :: value => xmlattribute_value
    procedure :: next_attribute => xmlattribute_next_attribute
    procedure :: as_bool => xmlattribute_as_bool
    procedure :: as_int => xmlattribute_as_int
    procedure :: as_llong => xmlattribute_as_llong
    procedure :: as_double => xmlattribute_as_double
    procedure :: associated => xmlattribute_associated
  end type XMLAttribute

  !=============================================================================
  ! XMLText C interface and derived type

  interface
    function xml_text_as_bool(text) result(x) bind(C)
      import C_PTR, C_BOOL
      type(C_PTR), value :: text
      logical(C_BOOL) :: x
    end function xml_text_as_bool

    function xml_text_as_int(text) result(x) bind(C)
      import C_PTR, C_INT
      type(C_PTR), value :: text
      integer(C_INT) :: x
    end function xml_text_as_int

    function xml_text_as_llong(text) result(x) bind(C)
      import C_PTR, C_LONG_LONG
      type(C_PTR), value :: text
      integer(C_LONG_LONG) :: x
    end function xml_text_as_llong

    function xml_text_as_double(text) result(x) bind(C)
      import C_PTR, C_DOUBLE
      type(C_PTR), value :: text
      real(C_DOUBLE) :: x
    end function xml_text_as_double
  end interface

  type :: XMLText
    type(C_PTR) :: ptr
  contains
    procedure :: as_bool => xmltext_as_bool
    procedure :: as_int => xmltext_as_int
    procedure :: as_llong => xmltext_as_llong
    procedure :: as_double => xmltext_as_double
  end type XMLText

  !=============================================================================
  ! XMLDocument C interface and derived type

  interface
    function xml_document_load_file(filename) result(doc) bind(C)
      import C_PTR
      type(C_PTR), value :: filename
      type(C_PTR) :: doc
    end function xml_document_load_file

    function xml_document_document_element(doc) result(node) bind(C)
      import C_PTR
      type(C_PTR), value :: doc
      type(C_PTR) :: node
    end function xml_document_document_element

    subroutine xml_document_clear(doc) bind(C)
      import C_PTR
      type(C_PTR), value :: doc
    end subroutine xml_document_clear
  end interface

  type, extends(XMLNode) :: XMLDocument
  contains
    procedure :: load_file => xmldocument_load_file
    procedure :: document_element => xmldocument_document_element
    procedure :: clear => xmldocument_clear
  end type XMLDocument

  public :: XMLNode, XMLAttribute, XMLText, XMLDocument

contains

  !=============================================================================
  ! XMLNode Implementation

  function xmlnode_name(this) result(name)
    class(XMLNode), intent(in) :: this
    character(len=:, kind=C_CHAR), allocatable :: name

    character(kind=C_CHAR), pointer :: string(:)
    integer(C_SIZE_T) :: size_string
    integer(C_SIZE_T) :: i
    type(C_PTR) :: name_ptr

    name_ptr = xml_node_name(this % ptr)
    size_string = strlen(name_ptr)
    call c_f_pointer(name_ptr, string, [size_string])
    allocate(character(len=size_string, kind=C_CHAR) :: name)
    do i = 1, size_string
      name(i:i) = string(i)
    end do
  end function xmlnode_name

  function xmlnode_child(this, name) result(child)
    class(XMLNode), intent(in) :: this
    character(len=*), optional :: name
    type(XMLNode) :: child

    integer :: i, n
    character(len=:, kind=C_CHAR), target, allocatable :: string

    if (present(name)) then
      allocate(character(len=len_trim(name) + 1, kind=C_CHAR) :: string)
      n = len_trim(name)
      do i = 1, n
        string(i:i) = name(i:i)
      end do
      string(n+1:n+1) = C_NULL_CHAR
      child % ptr = xml_node_child(this % ptr, c_loc(string))
    else
      child % ptr = xml_node_child(this % ptr, C_NULL_PTR)
    end if
  end function xmlnode_child

  function xmlnode_next_sibling(this, name) result(next)
    class(XMLNode), intent(in) :: this
    character(len=*), optional :: name
    type(XMLNode) :: next

    integer :: i, n
    character(len=:, kind=C_CHAR), target, allocatable :: string

    if (present(name)) then
      allocate(character(len=len_trim(name) + 1, kind=C_CHAR) :: string)
      n = len_trim(name)
      do i = 1, n
        string(i:i) = name(i:i)
      end do
      string(n+1:n+1) = C_NULL_CHAR
      next % ptr = xml_node_next_sibling(this % ptr, c_loc(string))
    else
      next % ptr = xml_node_next_sibling(this % ptr, C_NULL_PTR)
    end if
  end function xmlnode_next_sibling

  function xmlnode_attribute(this, name) result(attribute)
    class(XMLNode), intent(in) :: this
    character(len=*), intent(in) :: name
    type(XMLAttribute) :: attribute

    integer :: i, n
    character(kind=C_CHAR), target :: string(len_trim(name) + 1)

    n = len_trim(name)
    do i = 1, n
      string(i) = name(i:i)
    end do
    string(n+1) = C_NULL_CHAR
    attribute % ptr = xml_node_attribute(this % ptr, c_loc(string))
  end function xmlnode_attribute

  function xmlnode_child_value(this) result(val)
    class(XMLNode), intent(in) :: this
    character(len=:, kind=C_CHAR), allocatable :: val

    character(kind=C_CHAR), pointer :: string(:)
    integer(C_SIZE_T) :: size_string
    integer(C_SIZE_T) :: i
    type(C_PTR) :: text

    text = xml_node_child_value(this % ptr)
    size_string = strlen(text)
    call c_f_pointer(text, string, [size_string])
    allocate(character(len=size_string, kind=C_CHAR) :: val)
    do i = 1, size_string
      val(i:i) = string(i)
    end do
  end function xmlnode_child_value

  function xmlnode_text(this) result(text)
    class(XMLNode), intent(in) :: this
    type(XMLText) :: text

    text % ptr = xml_node_text(this % ptr)
  end function xmlnode_text

  logical function xmlnode_associated(this)
    class(XMLNode), intent(in) :: this

    xmlnode_associated = c_associated(this % ptr)
  end function xmlnode_associated

  !=============================================================================
  ! XMLAttribute Implementation

  function xmlattribute_name(this) result(name)
    class(XMLAttribute), intent(in) :: this
    character(len=:, kind=C_CHAR), allocatable :: name

    character(kind=C_CHAR), pointer :: string(:)
    integer(C_SIZE_T) :: size_string
    integer(C_SIZE_T) :: i
    type(C_PTR) :: name_ptr

    name_ptr = xml_attribute_name(this % ptr)
    size_string = strlen(name_ptr)
    call c_f_pointer(name_ptr, string, [size_string])
    allocate(character(len=size_string, kind=C_CHAR) :: name)
    do i = 1, size_string
      name(i:i) = string(i)
    end do
  end function xmlattribute_name

  function xmlattribute_value(this) result(val)
    class(XMLAttribute), intent(in) :: this
    character(len=:, kind=C_CHAR), allocatable :: val

    character(kind=C_CHAR), pointer :: string(:)
    integer(C_SIZE_T) :: size_string
    integer(C_SIZE_T) :: i
    type(C_PTR) :: val_ptr

    val_ptr = xml_attribute_value(this % ptr)
    size_string = strlen(val_ptr)
    call c_f_pointer(val_ptr, string, [size_string])
    allocate(character(len=size_string, kind=C_CHAR) :: val)
    do i = 1, size_string
      val(i:i) = string(i)
    end do
  end function xmlattribute_value

  function xmlattribute_next_attribute(this) result(next)
    class(XMLAttribute), intent(in) :: this
    type(XMLAttribute) :: next

    next % ptr = xml_attribute_next_attribute(this % ptr)
  end function xmlattribute_next_attribute

  logical(C_BOOL) function xmlattribute_as_bool(this)
    class(XMLAttribute), intent(in) :: this
    xmlattribute_as_bool = xml_attribute_as_bool(this % ptr)
  end function xmlattribute_as_bool

  integer(C_INT) function xmlattribute_as_int(this)
    class(XMLAttribute), intent(in) :: this
    xmlattribute_as_int = xml_attribute_as_int(this % ptr)
  end function xmlattribute_as_int

  integer(C_LONG_LONG) function xmlattribute_as_llong(this)
    class(XMLAttribute), intent(in) :: this
    xmlattribute_as_llong = xml_attribute_as_llong(this % ptr)
  end function xmlattribute_as_llong

  real(C_DOUBLE) function xmlattribute_as_double(this)
    class(XMLAttribute), intent(in) :: this
    xmlattribute_as_double = xml_attribute_as_double(this % ptr)
  end function xmlattribute_as_double

  logical function xmlattribute_associated(this)
    class(XMLAttribute), intent(in) :: this

    xmlattribute_associated = c_associated(this % ptr)
  end function xmlattribute_associated

  !=============================================================================
  ! XMLText Implementation

  logical(C_BOOL) function xmltext_as_bool(this)
    class(XMLText), intent(in) :: this
    xmltext_as_bool = xml_text_as_bool(this % ptr)
  end function xmltext_as_bool

  integer(C_INT) function xmltext_as_int(this)
    class(XMLText), intent(in) :: this
    xmltext_as_int = xml_text_as_int(this % ptr)
  end function xmltext_as_int

  integer(C_LONG_LONG) function xmltext_as_llong(this)
    class(XMLText), intent(in) :: this
    xmltext_as_llong = xml_text_as_llong(this % ptr)
  end function xmltext_as_llong

  real(C_DOUBLE) function xmltext_as_double(this)
    class(XMLText), intent(in) :: this
    xmltext_as_double = xml_text_as_double(this % ptr)
  end function xmltext_as_double

  !=============================================================================
  ! XMLDocument Implementation

  subroutine xmldocument_load_file(this, filename)
    class(XMLDocument), intent(inout) :: this
    character(*), intent(in) :: filename

    integer :: i, n
    character(kind=C_CHAR), target :: string(len_trim(filename) + 1)

    n = len_trim(filename)
    do i = 1, n
      string(i) = filename(i:i)
    end do
    string(n+1) = C_NULL_CHAR
    this % ptr = xml_document_load_file(c_loc(string))
  end subroutine xmldocument_load_file

  function xmldocument_document_element(this) result(node)
    class(XMLDocument), intent(in) :: this
    type(XMLNode) :: node

    node % ptr = xml_document_document_element(this % ptr)
  end function xmldocument_document_element

  subroutine xmldocument_clear(this)
    class(XMLDocument), intent(inout) :: this

    call xml_document_clear(this % ptr)
  end subroutine xmldocument_clear

end module pugixml
