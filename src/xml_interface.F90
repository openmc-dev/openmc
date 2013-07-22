module xml_interface

  use constants, only: MAX_LINE_LEN
  use error,     only: fatal_error
  use global,    only: message 
  use fox_dom

  implicit none
  private
  public :: Node
  public :: open_xmldoc
  public :: close_xmldoc
  public :: check_for_node
  public :: get_node_ptr
  public :: get_node_value

  integer, parameter :: ATTR_NODE = 0
  integer, parameter :: ELEM_NODE = 1

  interface get_node_value
    module procedure get_node_value_integer
    module procedure get_node_value_long
    module procedure get_node_value_string
  end interface get_node_value

contains

!===============================================================================
! OPEN_XMLDOC
!===============================================================================

  subroutine open_xmldoc(ptr, filename)

    character(len=*) :: filename
    type(Node), pointer :: ptr

    ptr => parseFile(trim(filename))
    ptr => getDocumentElement(ptr)

  end subroutine open_xmldoc

!===============================================================================
! CLOSE_XMLDOC
!===============================================================================

  subroutine close_xmldoc(ptr)

    type(Node), pointer :: ptr

    call destroy(ptr)

  end subroutine close_xmldoc

!===============================================================================
! CHECK_FOR_NODE
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
    elem_list => getElementsByTagName(ptr, trim(node_name))

    ! Get the length of the list
    if (getLength(elem_list) == 0) then
      found = .false.
      return
    end if

  end function check_for_node

!===============================================================================
! GET_NODE_PTR
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
    elem_list => getElementsByTagName(in_ptr, trim(node_name))

    ! Get the length of the list
    if (getLength(elem_list) == 0) return

    ! Point to the new element
    out_ptr => item(elem_list, 0)
    found_ = .true.

    ! Check to output found
    if (present(found)) found = found_

  end subroutine get_node_ptr

!===============================================================================
! GET_NODE_VALUE_INTEGER
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
! GET_NODE_VALUE_INTEGER
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
! GET_NODE_VALUE_STRING
!===============================================================================

  subroutine get_node_value_string(ptr, node_name, result_str)

    character(len=*), intent(in) :: node_name
    character(len=MAX_LINE_LEN) :: result_str
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
! GET_NODE
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
    elem_list => getElementsByTagName(in_ptr, trim(node_name))

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
