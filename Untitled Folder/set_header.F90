module set_header

!===============================================================================
! SET_HEADER module
!
! This module provides an implementation of sets based on the list
! implementation in list_header. The underlying datatype is a list, so adding an
! element just checks if the element is already in the list, and if not it's
! added. This results in much worse performance than an implementation based on
! hash tables or binary trees, but for our purposes, we don't expect to have
! gigantic sets where performance is critical.
!===============================================================================

  use constants, only: MAX_WORD_LEN
  use list_header

  implicit none

!===============================================================================
! SET contains a list of elements and methods to add, remove, and perform other
! basic tasks.
!===============================================================================

  type :: SetInt
    private
    type(ListInt) :: elements
  contains
    procedure :: add => set_add_int
    procedure :: clear => set_clear_int
    procedure :: contains => set_contains_int
    procedure :: get_item => set_get_item_int
    procedure :: remove => set_remove_int
    procedure :: size => set_size_int
  end type SetInt

  type :: SetChar
    private
    type(ListChar) :: elements
  contains
    procedure :: add => set_add_char
    procedure :: clear => set_clear_char
    procedure :: contains => set_contains_char
    procedure :: get_item => set_get_item_char
    procedure :: remove => set_remove_char
    procedure :: size => set_size_char
  end type SetChar

contains

!===============================================================================
! SET_ADD adds an item to a set if it is not already present in the set
!===============================================================================

  subroutine set_add_int(this, data)
    class(SetInt) :: this
    integer       :: data

    if (.not. this % elements % contains(data)) then
      call this % elements % append(data)
    end if

  end subroutine set_add_int

  subroutine set_add_char(this, data)
    class(SetChar) :: this
    character(*)   :: data

    if (.not. this % elements % contains(data)) then
      call this % elements % append(data)
    end if

  end subroutine set_add_char

!===============================================================================
! SET_CLEAR removes all items in a set
!===============================================================================

  subroutine set_clear_int(this)
    class(SetInt) :: this

    call this % elements % clear()

  end subroutine set_clear_int

  subroutine set_clear_char(this)
    class(SetChar) :: this

    call this % elements % clear()

  end subroutine set_clear_char

!===============================================================================
! SET_CONTAINS determines if a specified item is in a set
!===============================================================================

  function set_contains_int(this, data) result(in_set)
    class(SetInt) :: this
    integer       :: data
    logical       :: in_set

    in_set = this % elements % contains(data)

  end function set_contains_int

  function set_contains_char(this, data) result(in_set)
    class(SetChar) :: this
    character(*)   :: data
    logical        :: in_set

    in_set = this % elements % contains(data)

  end function set_contains_char

!===============================================================================
! SET_GET_ITEM returns the i-th item in the set
!===============================================================================

  function set_get_item_int(this, i_list) result(data)

    class(SetInt) :: this
    integer :: i_list
    integer :: data

    data = this % elements % get_item(i_list)

  end function set_get_item_int

  function set_get_item_char(this, i_list) result(data)

    class(SetChar)          :: this
    integer                 :: i_list
    character(MAX_WORD_LEN) :: data

    data = this % elements % get_item(i_list)

  end function set_get_item_char

!===============================================================================
! SET_REMOVE removes the specified item from the set. If it is not in the set,
! no action is taken.
!===============================================================================

  subroutine set_remove_int(this, data)

    class(SetInt) :: this
    integer       :: data

    call this % elements % remove(data)

  end subroutine set_remove_int

  subroutine set_remove_char(this, data)

    class(SetChar) :: this
    character(*)   :: data

    call this % elements % remove(data)

  end subroutine set_remove_char

!===============================================================================
! SET_SIZE returns the number of elements in the set
!===============================================================================

  function set_size_int(this) result(size)

    class(SetInt) :: this
    integer       :: size

    size = this % elements % size()

  end function set_size_int

  function set_size_char(this) result(size)

    class(SetChar) :: this
    integer        :: size

    size = this % elements % size()

  end function set_size_char

end module set_header
