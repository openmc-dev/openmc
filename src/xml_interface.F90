module xml_interface

  use, intrinsic :: ISO_C_BINDING

  use error,  only: fatal_error
  use pugixml
  use string, only: word_count

  implicit none
  private
  public :: XMLDocument, XMLNode, XMLText, XMLAttribute
  public :: check_for_node
  public :: get_node_list
  public :: get_node_value
  public :: get_node_value_bool
  public :: get_node_array
  public :: node_value_string
  public :: node_word_count

  interface get_node_value
    module procedure get_node_value_bool
    module procedure get_node_value_cbool
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
! CHECK_FOR_NODE checks for an attribute or sub-element node with the given
! node name. This should only be used for checking a single occurance of a
! sub-element node.
!===============================================================================

  function check_for_node(node, name) result(found)
    type(XMLNode), intent(in) :: node
    character(len=*), intent(in) :: name
    logical :: found

    type(XMLAttribute) :: attr
    type(XMLNode) :: child

    ! Check if an attribute or  exists
    attr = node % attribute(name)
    if (attr % associated()) then
      found = .true.
    else
      child = node % child(name)
      if (child % associated()) then
        found = .true.
      else
        found = .false.
      end if
    end if
  end function check_for_node

!===============================================================================
! GET_NODE_LIST is used to get a pointer to a list of sub-element nodes
!===============================================================================

  subroutine get_node_list(node, name, node_list)
    type(XMLNode), intent(in) :: node
    character(*), intent(in) :: name
    type(XMLNode), allocatable, intent(out) :: node_list(:)

    integer :: i, n
    type(XMLNode) :: first, current

    first = node % child(name)

    ! Determine number of nodes
    n = 0
    current = first
    do while (current % associated())
      n = n + 1
      current = current % next_sibling(name)
    end do

    ! Allocate nodes
    allocate(node_list(n))
    current = first
    do i = 1, n
      node_list(i) = current
      current = current % next_sibling(name)
    end do
  end subroutine get_node_list

!===============================================================================
! GET_NODE_VALUE_BOOL returns a logical value from an attribute or node
!===============================================================================

  subroutine get_node_value_bool(node, name, val)
    type(XMLNode), intent(in) :: node
    character(*), intent(in) :: name
    logical, intent(out) :: val

    type(XMLAttribute) :: attr
    type(XMLNode) :: child
    type(XMLText) :: text

    attr = node % attribute(name)
    if (attr % associated()) then
      val = attr % as_bool()
    else
      child = node % child(name)
      if (child % associated()) then
        text = child % text()
        val = text % as_bool()
      else
        call fatal_error("No child XML node named '" // trim(name) // "'.")
      end if
    end if
  end subroutine get_node_value_bool

  subroutine get_node_value_cbool(node, name, val)
    type(XMLNode), intent(in) :: node
    character(*), intent(in) :: name
    logical(kind=C_BOOL), intent(out) :: val

    logical :: val_
    call get_node_value_bool(node, name, val_)
    val = val_
  end subroutine get_node_value_cbool

!===============================================================================
! GET_NODE_VALUE_INTEGER returns a integer value from an attribute or node
!===============================================================================

  subroutine get_node_value_integer(node, name, val)
    type(XMLNode), intent(in) :: node
    character(*), intent(in) :: name
    integer, intent(out) :: val

    type(XMLAttribute) :: attr
    type(XMLNode) :: child
    type(XMLText) :: text

    attr = node % attribute(name)
    if (attr % associated()) then
      val = attr % as_int()
    else
      child = node % child(name)
      if (child % associated()) then
        text = child % text()
        val = text % as_int()
      else
        call fatal_error("No XML node named '" // trim(name) // "'.")
      end if
    end if
  end subroutine get_node_value_integer

!===============================================================================
! GET_NODE_VALUE_LONG returns an 8-byte integer from attribute or node
!===============================================================================

  subroutine get_node_value_long(node, name, val)
    type(XMLNode), intent(in) :: node
    character(*), intent(in) :: name
    integer(C_LONG_LONG), intent(out) :: val

    type(XMLAttribute) :: attr
    type(XMLNode) :: child
    type(XMLText) :: text

    attr = node % attribute(name)
    if (attr % associated()) then
      val = attr % as_llong()
    else
      child = node % child(name)
      if (child % associated()) then
        text = child % text()
        val = text % as_llong()
      else
        call fatal_error("No XML node named '" // trim(name) // "'.")
      end if
    end if
  end subroutine get_node_value_long

!===============================================================================
! GET_NODE_VALUE_DOUBLE returns a double precision real from attr. or node
!===============================================================================

  subroutine get_node_value_double(node, name, val)
    type(XMLNode), intent(in) :: node
    character(*), intent(in) :: name
    real(C_DOUBLE), intent(out) :: val

    type(XMLAttribute) :: attr
    type(XMLNode) :: child
    type(XMLText) :: text

    attr = node % attribute(name)
    if (attr % associated()) then
      val = attr % as_double()
    else
      child = node % child(name)
      if (child % associated()) then
        text = child % text()
        val = text % as_double()
      else
        call fatal_error("No XML node named '" // trim(name) // "'.")
      end if
    end if
  end subroutine get_node_value_double

!===============================================================================
! NODE_VALUE_STRING returns a single string from attr. or node
!===============================================================================

  function node_value_string(node, name) result(val)
    type(XMLNode), intent(in) :: node
    character(*), intent(in) :: name
    character(len=:, kind=C_CHAR), allocatable :: val

    type(XMLAttribute) :: attr
    type(XMLNode) :: child

    attr = node % attribute(name)
    if (attr % associated()) then
      val = adjustl(attr % value())
    else
      child = node % child(name)
      if (child % associated()) then
        val = adjustl(child % child_value())
      else
        call fatal_error("No child XML node named '" // trim(name) // "'.")
      end if
    end if
  end function node_value_string

!===============================================================================
! GET_NODE_VALUE_STRING returns a single string from attr. or node
!===============================================================================

  subroutine get_node_value_string(node, name, val)
    type(XMLNode), intent(in) :: node
    character(*), intent(in) :: name
    character(*), intent(out) :: val

    val = node_value_string(node, name)
  end subroutine get_node_value_string

!===============================================================================
! NODE_WORD_COUNT returns the number of words in a text node or attribute value
!===============================================================================

  integer function node_word_count(node, name)
    type(XMLNode), intent(in) :: node
    character(*), intent(in) :: name

    node_word_count = word_count(node_value_string(node, name))
  end function node_word_count

!===============================================================================
! GET_NODE_ARRAY_INTEGER returns a 1-D array of integers
!===============================================================================

  subroutine get_node_array_integer(node, name, array)
    type(XMLNode), intent(in) :: node
    character(*), intent(in) :: name
    integer, intent(out) :: array(:)

    integer :: stat
    character(len=:, kind=C_CHAR), allocatable :: str

    ! Get value of text node/attribute
    str = node_value_string(node, name)
    call whitespace_to_blanks(str)

    ! Read numbers into array
    read(UNIT=str, FMT=*, IOSTAT=stat) array
    if (stat > 0) then
      call fatal_error("Error converting XML text: " // trim(str))
    end if
  end subroutine get_node_array_integer

!===============================================================================
! GET_NODE_ARRAY_DOUBLE returns a 1-D array of double precision reals
!===============================================================================

  subroutine get_node_array_double(node, name, array)
    type(XMLNode), intent(in) :: node
    character(*), intent(in) :: name
    real(C_DOUBLE), intent(out) :: array(:)

    integer :: stat
    character(len=:, kind=C_CHAR), allocatable :: str

    ! Get value of text node/attribute
    str = node_value_string(node, name)
    call whitespace_to_blanks(str)

    ! Read numbers into array
    read(UNIT=str, FMT=*, IOSTAT=stat) array
    if (stat > 0) then
      call fatal_error("Error converting XML text: " // trim(str))
    end if
  end subroutine get_node_array_double

!===============================================================================
! GET_NODE_ARRAY_STRING returns a 1-D array of strings
!===============================================================================

  subroutine get_node_array_string(node, name, array)
    type(XMLNode), intent(in) :: node
    character(*), intent(in) :: name
    character(len=*), intent(out) :: array(:)

    integer :: i, n
    integer :: start
    logical :: inword
    character(kind=C_CHAR) :: chr
    character(len=:, kind=C_CHAR), allocatable :: str

    ! Get value of text node/attribute
    str = node_value_string(node, name)

    inword = .false.
    n = 0
    do i = 1, len(str)
      chr = str(i:i)
      select case (chr)
      case (' ', C_HORIZONTAL_TAB, C_NEW_LINE, C_CARRIAGE_RETURN)
        if (inword) then
          inword = .false.
          n = n + 1
          array(n) = str(start:i-1)
        end if
      case default
        if (.not. inword) then
          start = i
          inword = .true.
        end if
      end select
    end do
    if (inword) then
      n = n + 1
      array(n) = str(start:len(str))
    end if
  end subroutine get_node_array_string

!===============================================================================
! WHITESPACE_TO_BLANKS converts all whitespace to blanks
!===============================================================================

  subroutine whitespace_to_blanks(str)
    character(len=*, kind=C_CHAR), intent(inout) :: str

    integer :: i

    do i = 1, len(str)
      select case (str(i:i))
      case (C_NEW_LINE, C_HORIZONTAL_TAB, C_CARRIAGE_RETURN)
        str(i:i) = ' '
      end select
    end do
  end subroutine

end module xml_interface
