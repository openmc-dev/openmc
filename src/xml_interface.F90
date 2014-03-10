module xml_interface

  use constants, only: MAX_LINE_LEN
  use error,     only: fatal_error
  use global,    only: message 
  use fox_dom

  implicit none
  private
  public :: Node
  public :: NodeList
  public :: open_xmldoc
  public :: close_xmldoc
  public :: check_for_node
  public :: get_node_ptr
  public :: get_node_list
  public :: get_list_size
  public :: get_list_item
  public :: get_node_value
  public :: get_node_array
  public :: get_arraysize_integer
  public :: get_arraysize_double
  public :: get_arraysize_string

  integer, parameter :: ATTR_NODE = 0
  integer, parameter :: ELEM_NODE = 1

  interface get_node_value
    module procedure get_node_value_integer
    module procedure get_node_value_long
    module procedure get_node_value_double
    module procedure get_node_value_string
  end interface get_node_value

  interface get_node_array
    module procedure get_node_array_integer
    module procedure get_node_array_double
    module procedure get_node_array_string
  end interface get_node_array

contains

!===============================================================================
! OPEN_XMLDOC opens and parses an XML and returns a pointer to the document
!===============================================================================

  subroutine open_xmldoc(ptr, filename)

    character(len=*) :: filename
    type(Node), pointer :: ptr

    ptr => parseFile(trim(filename))  ! Pointer to the whole XML document
    ptr => getDocumentElement(ptr)    ! Grabs root element of XML document

  end subroutine open_xmldoc

!===============================================================================
! CLOSE_XMLDOC closes and destroys all memory associated with document
!===============================================================================

  subroutine close_xmldoc(ptr)

    type(Node), pointer :: ptr

    ptr => getParentNode(ptr)  ! Go from root element ptr to document ptr
    call destroy(ptr)          ! Deallocates all nodes recursively

  end subroutine close_xmldoc

!===============================================================================
! CHECK_FOR_NODE checks for an attribute or sub-element node with the given
! node name. This should only be used for checking a single occurance of a
! sub-element node. To check for sub-element nodes that repeat, use
! get_node_list and get_list_size. This is to minimize number of calls
! to getChildrenByTagName.
!===============================================================================

  function check_for_node(ptr, node_name) result(found)

    character(len=*), intent(in) :: node_name
    logical :: found
    type(Node), pointer, intent(in) :: ptr

    type(Node), pointer :: temp_ptr
    type(NodeList), pointer :: elem_list

    ! Default that we found it and attribute
    found = .true.

    ! Get the attribute node
    temp_ptr => getAttributeNode(ptr, trim(node_name))

    ! Check if node exists
    if (associated(temp_ptr)) return

    ! Check for a sub-element 
    elem_list => getChildrenByTagName(ptr, trim(node_name))

    ! Get the length of the list
    if (getLength(elem_list) == 0) then
      found = .false.
      return
    end if

  end function check_for_node

!===============================================================================
! GET_NODE_PTR returns a Node pointer to an attribute or sub-element node.
! Note, this should only be used for sub-element nodes that occur only once.
! For repeated nodes, use get_node_list and then get_item_item.
!===============================================================================

  subroutine get_node_ptr(in_ptr, node_name, out_ptr, found)

    character(len=*), intent(in) :: node_name
    logical, intent(out), optional :: found
    type(Node), pointer, intent(in) :: in_ptr
    type(Node), pointer, intent(out) :: out_ptr

    logical :: found_
    type(NodeList), pointer :: elem_list => null()

    ! Set found to false
    found_ = .false.

    ! Check for a sub-element 
    elem_list => getChildrenByTagName(in_ptr, trim(node_name))

    ! Get the length of the list
    if (getLength(elem_list) == 0) return

    ! Point to the new element
    out_ptr => item(elem_list, 0)
    found_ = .true.

    ! Check to output found
    if (present(found)) found = found_

  end subroutine get_node_ptr

!===============================================================================
! GET_NODE_LIST is used to get a pointer to a list of sub-element nodes
!===============================================================================

  subroutine get_node_list(in_ptr, node_name, out_ptr)

    character(len=*), intent(in) :: node_name
    type(Node), pointer, intent(in) :: in_ptr
    type(NodeList), pointer, intent(out) :: out_ptr

    ! Check for a sub-element 
    out_ptr => getChildrenByTagName(in_ptr, trim(node_name))

  end subroutine get_node_list

!===============================================================================
! GET_LIST_SIZE is used to get the number of elements from a node list
!===============================================================================

  function get_list_size(in_ptr) result(n_size)

    integer :: n_size
    type(NodeList), pointer, intent(in) :: in_ptr

    ! Get the size of the list
    n_size = getLength(in_ptr)

  end function get_list_size

!===============================================================================
! GET_LIST_ITEM is used to select a specific item from a node list
!===============================================================================

  subroutine get_list_item(in_ptr, idx, out_ptr)

    integer, intent(in) :: idx
    type(NodeList), pointer, intent(in) :: in_ptr
    type(Node), pointer, intent(out) :: out_ptr

    ! Check for a sub-element 
    out_ptr => item(in_ptr, idx - 1)

  end subroutine get_list_item

!===============================================================================
! GET_NODE_VALUE_INTEGER returns a integer value from an attribute or node
!===============================================================================

  subroutine get_node_value_integer(ptr, node_name, result_int)

    character(len=*), intent(in) :: node_name
    integer :: result_int
    type(Node), pointer, intent(in) :: ptr

    integer :: node_type
    logical :: found
    type(Node), pointer :: temp_ptr

    ! Get pointer to the node
    call get_node(ptr, node_name, temp_ptr, node_type, found)

    ! Leave if it was not found
    if (.not. found) then 
      message = "Node " // node_name // " not part of Node " // &
                getNodeName(ptr) // "."
      call fatal_error()
    end if

    ! Extract value
    if (node_type == ATTR_NODE) then
      call extractDataAttribute(ptr, node_name, result_int)
    else
      call extractDataContent(temp_ptr, result_int)
    end if
    
  end subroutine get_node_value_integer

!===============================================================================
! GET_NODE_VALUE_LONG returns an 8-byte integer from attribute or node
!===============================================================================

  subroutine get_node_value_long(ptr, node_name, result_long)

    character(len=*), intent(in) :: node_name
    integer(8) :: result_long
    type(Node), pointer, intent(in) :: ptr

    integer :: node_type
    logical :: found
    type(Node), pointer :: temp_ptr

    ! Get pointer to the node
    call get_node(ptr, node_name, temp_ptr, node_type, found)

    ! Leave if it was not found
    if (.not. found) then 
      message = "Node " // node_name // " not part of Node " // &
                getNodeName(ptr) // "."
      call fatal_error()
    end if
      
    ! Extract value
    if (node_type == ATTR_NODE) then
      call extractDataAttribute(ptr, node_name, result_long)
    else
      call extractDataContent(temp_ptr, result_long)
    end if
    
  end subroutine get_node_value_long

!===============================================================================
! GET_NODE_VALUE_DOUBLE returns a double precision real from attr. or node
!===============================================================================

  subroutine get_node_value_double(ptr, node_name, result_double)

    character(len=*), intent(in) :: node_name
    real(8) :: result_double
    type(Node), pointer, intent(in) :: ptr

    integer :: node_type
    logical :: found
    type(Node), pointer :: temp_ptr

    ! Get pointer to the node
    call get_node(ptr, node_name, temp_ptr, node_type, found)

    ! Leave if it was not found
    if (.not. found) then 
      message = "Node " // node_name // " not part of Node " // &
                getNodeName(ptr) // "."
      call fatal_error()
    end if
      
    ! Extract value
    if (node_type == ATTR_NODE) then
      call extractDataAttribute(ptr, node_name, result_double)
    else
      call extractDataContent(temp_ptr, result_double)
    end if
    
  end subroutine get_node_value_double

!===============================================================================
! GET_NODE_ARRAY_INTEGER returns a 1-D array of integers
!===============================================================================

  subroutine get_node_array_integer(ptr, node_name, result_int)

    character(len=*), intent(in) :: node_name
    integer :: result_int(:)
    type(Node), pointer, intent(in) :: ptr

    integer :: node_type
    logical :: found
    type(Node), pointer :: temp_ptr

    ! Get pointer to the node
    call get_node(ptr, node_name, temp_ptr, node_type, found)

    ! Leave if it was not found
    if (.not. found) then 
      message = "Node " // node_name // " not part of Node " // &
                getNodeName(ptr) // "."
      call fatal_error()
    end if
      
    ! Extract value
    if (node_type == ATTR_NODE) then
      call extractDataAttribute(ptr, node_name, result_int)
    else
      call extractDataContent(temp_ptr, result_int)
    end if
    
  end subroutine get_node_array_integer

!===============================================================================
! GET_NODE_ARRAY_DOUBLE returns a 1-D array of double precision reals
!===============================================================================

  subroutine get_node_array_double(ptr, node_name, result_double)

    character(len=*), intent(in) :: node_name
    real(8) :: result_double(:)
    type(Node), pointer, intent(in) :: ptr

    integer :: node_type
    logical :: found
    type(Node), pointer :: temp_ptr

    ! Get pointer to the node
    call get_node(ptr, node_name, temp_ptr, node_type, found)

    ! Leave if it was not found
    if (.not. found) then 
      message = "Node " // node_name // " not part of Node " // &
                getNodeName(ptr) // "."
      call fatal_error()
    end if
      
    ! Extract value
    if (node_type == ATTR_NODE) then
      call extractDataAttribute(ptr, node_name, result_double)
    else
      call extractDataContent(temp_ptr, result_double)
    end if
    
  end subroutine get_node_array_double

!===============================================================================
! GET_NODE_ARRAY_STRING returns a 1-D array of strings
!===============================================================================

  subroutine get_node_array_string(ptr, node_name, result_string)

    character(len=*), intent(in) :: node_name
    character(len=*) :: result_string(:)
    type(Node), pointer, intent(in) :: ptr

    integer :: node_type
    logical :: found
    type(Node), pointer :: temp_ptr

    ! Get pointer to the node
    call get_node(ptr, node_name, temp_ptr, node_type, found)

    ! Leave if it was not found
    if (.not. found) then 
      message = "Node " // node_name // " not part of Node " // &
                getNodeName(ptr) // "."
      call fatal_error()
    end if
      
    ! Extract value
    if (node_type == ATTR_NODE) then
      call extractDataAttribute(ptr, node_name, result_string)
    else
      call extractDataContent(temp_ptr, result_string)
    end if
    
  end subroutine get_node_array_string

!===============================================================================
! GET_NODE_VALUE_STRING returns a single string from attr. or node
!===============================================================================

  subroutine get_node_value_string(ptr, node_name, result_str)

    character(len=*), intent(in) :: node_name
    character(len=*) :: result_str
    type(Node), pointer, intent(in) :: ptr

    integer :: node_type
    logical :: found
    type(Node), pointer :: temp_ptr

    ! Get pointer to the node
    call get_node(ptr, node_name, temp_ptr, node_type, found)

    ! Leave if it was not found
    if (.not. found) then 
      message = "Node " // node_name // " not part of Node " // &
                getNodeName(ptr) // "."
      call fatal_error()
    end if
 
    ! Extract value
    if (node_type == ATTR_NODE) then
      call extractDataAttribute(ptr, node_name, result_str)
    else
      call extractDataContent(temp_ptr, result_str)
    end if
   
  end subroutine get_node_value_string

!===============================================================================
! GET_NODE_ARRAYSIZE_INTEGER returns the size of the integer array
!===============================================================================

  function get_arraysize_integer(ptr, node_name) result(n)

    character(len=*), intent(in) :: node_name
    integer :: n
    type(Node), pointer, intent(in) :: ptr

    integer :: node_type
    logical :: found
    type(Node), pointer :: temp_ptr

    ! Get pointer to the node
    call get_node(ptr, node_name, temp_ptr, node_type, found)

    ! Leave if it was not found
    if (.not. found) then
      message = "Node " // node_name // " not part of Node " // &
                getNodeName(ptr) // "."
      call fatal_error()
    end if

    ! Get the size
    n = countrts(getTextContent(temp_ptr), 0)

  end function get_arraysize_integer

!===============================================================================
! GET_NODE_ARRAYSIZE_DOUBLE returns the size of double prec. real array
!===============================================================================

  function get_arraysize_double(ptr, node_name) result(n)

    character(len=*), intent(in) :: node_name
    integer :: n
    type(Node), pointer, intent(in) :: ptr

    integer :: node_type
    logical :: found
    type(Node), pointer :: temp_ptr

    ! Get pointer to the node
    call get_node(ptr, node_name, temp_ptr, node_type, found)

    ! Leave if it was not found
    if (.not. found) then
      message = "Node " // node_name // " not part of Node " // &
                getNodeName(ptr) // "."
      call fatal_error()
    end if

    ! Get the size
    n = countrts(getTextContent(temp_ptr), 0.0_8)

  end function get_arraysize_double

!===============================================================================
! GET_NODE_ARRAYSIZE_STRING returns the size of string array
!===============================================================================

  function get_arraysize_string(ptr, node_name) result(n)

    character(len=*), intent(in) :: node_name
    integer :: n
    type(Node), pointer, intent(in) :: ptr

    integer :: node_type
    logical :: found
    type(Node), pointer :: temp_ptr

    ! Get pointer to the node
    call get_node(ptr, node_name, temp_ptr, node_type, found)

    ! Leave if it was not found
    if (.not. found) then
      message = "Node " // node_name // " not part of Node " // &
                getNodeName(ptr) // "."
      call fatal_error()
    end if

    ! Get the size
    n = countrts(getTextContent(temp_ptr), 'a')

  end function get_arraysize_string

!===============================================================================
! GET_NODE private routine that gets a pointer to a specific node
!===============================================================================

  subroutine get_node(in_ptr, node_name, out_ptr, node_type, found)

    character(len=*), intent(in) :: node_name
    integer, intent(out) :: node_type
    logical, intent(out) :: found
    type(Node), pointer, intent(in) :: in_ptr
    type(Node), pointer, intent(out) :: out_ptr

    type(NodeList), pointer :: elem_list

    ! Default that we found it and attribute
    found = .true.
    node_type = ATTR_NODE

    ! Get the attribute node
    out_ptr => getAttributeNode(in_ptr, trim(node_name))

    ! Check if node exists
    if (associated(out_ptr)) return

    ! Check for a sub-element 
    elem_list => getChildrenByTagName(in_ptr, trim(node_name))

    ! Get the length of the list
    if (getLength(elem_list) == 0) then
      found = .false.
      return
    end if

    ! Point to the new element
    node_type = ELEM_NODE
    out_ptr => item(elem_list, 0)

  end subroutine get_node

end module xml_interface
