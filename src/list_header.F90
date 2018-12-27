module list_header

!===============================================================================
! LIST_HEADER module
!
! This module contains a linked list structure with convenience methods such as
! append, contains, remove, index, get_item, size, etc. This is an updated
! implementation with type-bound procedures (F2003).
!===============================================================================

  use constants, only: ERROR_INT, ERROR_REAL, MAX_WORD_LEN

  implicit none

!===============================================================================
! LISTELEM* types hold one piece of data and a pointer to the next piece of data
!===============================================================================

  type :: ListElemChar
    character(MAX_WORD_LEN) :: data
    type(ListElemChar), pointer :: next => null()
    type(ListElemChar), pointer :: prev => null()
  end type ListElemChar

!===============================================================================
! LIST* types contain the linked list with convenience methods. We originally
! considered using unlimited polymorphism to provide a single type, but compiler
! support is still spotty, and in many cases it doesn't prevent duplication of
! code. For the time being, a separate derived type exists for each datatype.
!===============================================================================

  type, public :: ListChar
    private
    integer :: count = 0 ! Number of elements in list

    ! Used in get_item for fast sequential lookups
    integer :: last_index = huge(0)
    type(ListElemChar), pointer :: last_elem => null()

    ! Pointers to beginning and end of list
    type(ListElemChar), public, pointer :: head => null()
    type(ListElemChar), public, pointer :: tail => null()
  contains
    procedure :: append => list_append_char     ! Add item to end of list
    procedure :: clear => list_clear_char       ! Remove all items
    procedure :: contains => list_contains_char ! Does list contain?
    procedure :: get_item => list_get_item_char ! Get i-th item in list
    procedure :: index => list_index_char       ! Determine index of given item
    procedure :: insert => list_insert_char     ! Insert item in i-th position
    procedure :: remove => list_remove_char     ! Remove specified item
    procedure :: size => list_size_char         ! Size of list
  end type ListChar

contains

!===============================================================================
! LIST_APPEND appends an item to the end of the list. If the list is empty, it
! becomes the first item.
!===============================================================================

  subroutine list_append_char(this, data)
    class(ListChar) :: this
    character(*)    :: data

    type(ListElemChar), pointer :: elem

    ! Create element and set dat
    allocate(elem)
    elem % data = data

    if (.not. associated(this % head)) then
      ! If list is empty, set head and tail to new element
      this % head => elem
      this % tail => elem
    else
      ! Otherwise append element at end of list
      this % tail % next => elem
      elem % prev => this % tail
      this % tail => this % tail % next
    end if

    this % count = this % count + 1

  end subroutine list_append_char

!===============================================================================
! LIST_CLEAR removes all elements from the list
!===============================================================================

  subroutine list_clear_char(this)
    class(ListChar) :: this

    type(ListElemChar), pointer :: current => null()
    type(ListElemChar), pointer :: next => null()

    if (this % count > 0) then
      current => this % head
      do while (associated(current))
        ! Set pointer to next element
        next => current % next

        ! Deallocate memory for current element
        deallocate(current)

        ! Move to next element
        current => next
      end do

      nullify(this % head)
      nullify(this % tail)
      this % count = 0
    end if

  end subroutine list_clear_char

!===============================================================================
! LIST_CONTAINS determines whether the list contains a specified item. Since it
! relies on the index method, it is O(n).
!===============================================================================

  function list_contains_char(this, data) result(in_list)
    class(ListChar) :: this
    character(*)    :: data
    logical         :: in_list

    in_list = (this % index(data) > 0)

  end function list_contains_char

!===============================================================================
! LIST_GET_ITEM returns the item in the list at position 'i_list'. If the index
! is out of bounds, an error code is returned.
! ===============================================================================

  function list_get_item_char(this, i_list) result(data)
    class(ListChar) :: this
    integer         :: i_list
    character(MAX_WORD_LEN) :: data

    integer :: last_index

    if (i_list < 1 .or. i_list > this % count) then
      ! Check for index out of bounds
      data = ""
    elseif (i_list == 1) then
      data = this % head % data
      this % last_index = 1
      this % last_elem => this % head
    elseif (i_list == this % count) then
      data = this % tail % data
      this % last_index = this % count
      this % last_elem => this % tail
    else
      if (i_list < this % last_index) then
        this % last_index = 1
        this % last_elem => this % head
      end if

      do last_index = this % last_index + 1, i_list
        this % last_elem => this % last_elem % next
        this % last_index = last_index
      end do
      data = this % last_elem % data
    end if

  end function list_get_item_char

!===============================================================================
! LIST_INDEX determines the first index in the list that contains 'data'. If
! 'data' is not present in the list, the return value is -1.
!===============================================================================

  function list_index_char(this, data) result(i_list)

    class(ListChar) :: this
    character(*)    :: data
    integer         :: i_list

    type(ListElemChar), pointer :: elem

    i_list = 0
    elem => this % head
    do while (associated(elem))
      i_list = i_list + 1
      if (data == elem % data) exit
      elem => elem % next
    end do

    ! Check if we reached the end of the list
    if (.not. associated(elem)) i_list = -1

  end function list_index_char

!===============================================================================
! LIST_INSERT inserts 'data' at index 'i_list' within the list. If 'i_list'
! exceeds the size of the list, the data is appends at the end of the list.
!===============================================================================

  subroutine list_insert_char(this, i_list, data)

    class(ListChar) :: this
    integer         :: i_list
    character(*)    :: data

    integer :: i
    type(ListElemChar), pointer :: elem => null()
    type(ListElemChar), pointer :: new_elem => null()

    if (i_list > this % count) then
      ! Check whether specified index is greater than number of elements -- if
      ! so, just append it to the end of the list
      call this % append(data)

    else if (i_list == 1) then
      ! Check for new head element
      allocate(new_elem)
      new_elem % data = data
      new_elem % next => this % head
      this % head => new_elem
      this % count = this % count + 1

    else
      ! Default case with new element somewhere in middle of list
      if (i_list >= this % last_index) then
        i = this % last_index
        elem => this % last_elem
      else
        i = 0
        elem => this % head
      end if
      do while (associated(elem))
        if (i == i_list) then
          ! Allocate new element
          allocate(new_elem)
          new_elem % data = data

          ! Put it before the i-th element
          new_elem % prev => elem % prev
          new_elem % next => elem
          new_elem % prev % next => new_elem
          new_elem % next % prev => new_elem
          this % count = this % count + 1
          this % last_index = i_list
          this % last_elem => new_elem
          exit
        end if
        i = i + 1
        elem => elem % next
      end do
    end if

  end subroutine list_insert_char

!===============================================================================
! LIST_REMOVE removes the first item in the list that contains 'data'. If 'data'
! is not in the list, no action is taken.
!===============================================================================

  subroutine list_remove_char(this, data)

    class(ListChar) :: this
    character(*)    :: data

    type(ListElemChar), pointer :: elem => null()

    elem => this % head
    do while (associated(elem))
      ! Check for matching data
      if (elem % data == data) then

        ! Determine whether the current element is the head, tail, or a middle
        ! element
        if (associated(elem, this % head)) then
          this % head => elem % next
          if (associated(elem, this % tail)) nullify(this % tail)
          if (associated(this % head)) nullify(this % head % prev)
          deallocate(elem)
        else if (associated(elem, this % tail)) then
          this % tail => elem % prev
          deallocate(this % tail % next)
        else
          elem % prev % next => elem % next
          elem % next % prev => elem % prev
          deallocate(elem)
        end if

        ! Decrease count and exit
        this % count = this % count - 1
        exit
      end if

      ! Advance pointers
      elem => elem % next
    end do

  end subroutine list_remove_char

!===============================================================================
! LIST_SIZE returns the number of elements in the list
!===============================================================================

  function list_size_char(this) result(size)

    class(ListChar) :: this
    integer         :: size

    size = this % count

  end function list_size_char

end module list_header
