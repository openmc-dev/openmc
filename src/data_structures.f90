module data_structures

!=====================================================================
! DATA_STRUCTURES module
!
! This module implements a dictionary that has (key,value) pairs. This
! data structure is used to provide lookup features, e.g. cells and
! surfaces by name.
!
! The original version from the 'flibs' open source package
! implemented another derived type called DICT_DATA that has been
! replaced here by a simple integer in ListData. If used in the
! original form , the dictionary can store derived types (changes made
! to ListData, dict_create, dict_add_key, and dict_get_key).
!=====================================================================

  use global, only: DICT_KEY_LENGTH

  implicit none

  ! Data stored in a linked list. In this case, we store the
  ! (key,value) pair for a dictionary. Note that we need to store the
  ! key in addition to the value for collision resolution.
  type ListData
     character(len=DICT_KEY_LENGTH) :: key
     integer                        :: value
  end type ListData

  ! A simple linked list
  type LinkedList
     type(LinkedList), pointer :: next
     type(ListData)            :: data
  end type LinkedList

  type HashList
     type(LinkedList), pointer :: list
  end type HashList

  ! A dictionary of (key,value) pairs
  type Dictionary
     private
     type(HashList), pointer, dimension(:) :: table
  end type Dictionary

  ! Hide objects that do not need to be publicly accessible
  private :: ListData
  private :: HashList
  private :: LinkedList
  private :: list_create
  private :: list_destroy
  private :: list_count
  private :: list_next
  private :: list_insert
  private :: list_insert_head
  private :: list_delete_element
  private :: list_get_data
  private :: list_put_data
  private :: dict_get_elem
  private :: dict_hashkey
  
  integer, parameter, private :: hash_size  = 4993
  integer, parameter, private :: multiplier = 31
  integer, parameter :: DICT_NULL = 0

contains

!=====================================================================
! LIST_CREATE creates and initializes a list
! Arguments:
!     list       Pointer to new linked list
!     data       The data for the first element
! Note:
!     This version assumes a shallow copy is enough
!     (that is, there are no pointers within the data
!     to be stored)
!     It also assumes the argument list does not already
!     refer to a list. Use list_destroy first to
!     destroy up an old list.
!=====================================================================

  subroutine list_create( list, data )

    type(LinkedList), pointer  :: list
    type(ListData), intent(in) :: data

    allocate( list )
    list%next => null()
    list%data =  data

  end subroutine list_create

!=====================================================================
! LIST_DESTROY destroys an entire list
! Arguments:
!     list       Pointer to the list to be destroyed
! Note:
!     This version assumes that there are no
!     pointers within the data that need deallocation
!=====================================================================

  subroutine list_destroy( list )

    type(LinkedList), pointer  :: list

    type(LinkedList), pointer  :: current
    type(LinkedList), pointer  :: next

    current => list
    do while ( associated(current%next) )
       next => current%next
       deallocate( current )
       current => next
    enddo

  end subroutine list_destroy

!=====================================================================
! LIST_COUNT count the number of items in the list
! Arguments:
!     list       Pointer to the list
!=====================================================================

  integer function list_count( list )

    type(LinkedList), pointer  :: list

    type(LinkedList), pointer  :: current
    type(LinkedList), pointer  :: next

    if ( associated(list) ) then
       list_count = 1
       current => list
       do while ( associated(current%next) )
          current => current%next
          list_count = list_count + 1
       enddo
    else
       list_count = 0
    endif

  end function list_count

!=====================================================================
! LIST_NEXT returns the next element (if any)
! Arguments:
!     elem       Element in the linked list
! Result:
!=====================================================================

  function list_next( elem ) result(next)

    type(LinkedList), pointer :: elem
    type(LinkedList), pointer :: next

    next => elem%next

  end function list_next

!=====================================================================
! LIST_INSERT inserts a new element
! Arguments:
!     elem       Element in the linked list after
!                which to insert the new element
!     data       The data for the new element
!=====================================================================

  subroutine list_insert( elem, data )

    type(LinkedList), pointer  :: elem
    type(ListData), intent(in) :: data

    type(LinkedList), pointer :: next

    allocate(next)

    next%next => elem%next
    elem%next => next
    next%data =  data

  end subroutine list_insert

!=====================================================================
! LIST_INSERT_HEAD inserts a new element before the first element
! Arguments:
!     list       Start of the list
!     data       The data for the new element
!=====================================================================

  subroutine list_insert_head( list, data )

    type(LinkedList), pointer  :: list
    type(ListData), intent(in) :: data

    type(LinkedList), pointer :: elem

    allocate(elem)
    elem%data =  data

    elem%next => list
    list      => elem

  end subroutine list_insert_head

!=====================================================================
! LIST_DELETE_ELEMENT deletes an element from the list
! Arguments:
!     list       Header of the list
!     elem       Element in the linked list to be removed
!=====================================================================

  subroutine list_delete_element( list, elem )

    type(LinkedList), pointer  :: list
    type(LinkedList), pointer  :: elem

    type(LinkedList), pointer  :: current
    type(LinkedList), pointer  :: prev

    if ( associated(list,elem) ) then
       list => elem%next
       deallocate( elem )
    else
       current => list
       prev    => list
       do while ( associated(current) )
          if ( associated(current,elem) ) then
             prev%next => current%next
             deallocate( current ) ! Is also "elem"
             exit
          endif
          prev    => current
          current => current%next
       enddo
    endif

  end subroutine list_delete_element

!=====================================================================
! LIST_GET_DATA gets the data stored with a list element
! Arguments:
!     elem       Element in the linked list
!=====================================================================

  function list_get_data( elem ) result(data)

    type(LinkedList), pointer :: elem
    type(ListData)            :: data

    data = elem%data

  end function list_get_data

!=====================================================================
! LIST_PUT_DATA stores new data with a list element
! Arguments:
!     elem       Element in the linked list
!     data       The data to be stored
!=====================================================================

  subroutine list_put_data( elem, data )

    type(LinkedList), pointer  :: elem
    type(ListData), intent(in) :: data

    elem%data = data

  end subroutine list_put_data

!=====================================================================
! DICT_CREATE creates and initializes a dictionary
! Arguments:
!     dict       Pointer to new dictionary
!     key        Key for the first element
!     value      Value for the first element
! Note:
!     This version assumes a shallow copy is enough
!     (that is, there are no pointers within the data
!     to be stored)
!     It also assumes the argument list does not already
!     refer to a list. Use dict_destroy first to
!     destroy up an old list.
!=====================================================================

  subroutine dict_create( dict, key, value )

    type(Dictionary), pointer   :: dict
    character(len=*), intent(in) :: key
    integer,  intent(in) :: value

    type(ListData)              :: data
    integer                      :: i
    integer                      :: hash

    allocate( dict )
    allocate( dict%table(hash_size) )

    do i = 1,hash_size
       dict%table(i)%list => null()
    enddo

    data%key   = key
    data%value = value

    hash = dict_hashkey( trim(key ) )
    call list_create( dict%table(hash)%list, data )

  end subroutine dict_create

!=====================================================================
! DICT_DESTROY destroys an entire dictionary
! Arguments:
!     dict       Pointer to the dictionary to be destroyed
! Note:
!     This version assumes that there are no
!     pointers within the data that need deallocation
!=====================================================================

  subroutine dict_destroy( dict )

    type(Dictionary), pointer  :: dict

    integer                     :: i

    do i = 1,size(dict%table)
       if ( associated( dict%table(i)%list ) ) then
          call list_destroy( dict%table(i)%list )
       endif
    enddo
    deallocate( dict%table )
    deallocate( dict )

  end subroutine dict_destroy

!=====================================================================
! DICT_ADD_KEY add a new keys
! Arguments:
!     dict       Pointer to the dictionary
!     key        Key for the new element
!     value      Value for the new element
! Note:
!     If the key already exists, the
!     key's value is simply replaced
!=====================================================================

  subroutine dict_add_key( dict, key, value )

    type(Dictionary), pointer   :: dict
    character(len=*), intent(in) :: key
    integer,  intent(in) :: value

    type(ListData)              :: data
    type(LinkedList), pointer   :: elem
    integer                      :: hash

    elem => dict_get_elem( dict, key )

    if ( associated(elem) ) then
       elem%data%value = value
    else
       data%key   = key
       data%value = value
       hash       = dict_hashkey( trim(key) )
       if ( associated( dict%table(hash)%list ) ) then
          call list_insert( dict%table(hash)%list, data )
       else
          call list_create( dict%table(hash)%list, data )
       endif
    endif

  end subroutine dict_add_key

!=====================================================================
! DICT_DELETE_KEY deletes a key-value pair from the dictionary
! Arguments:
!     dict       Dictionary in question
!     key        Key to be removed
!=====================================================================

  subroutine dict_delete_key( dict, key )

    type(Dictionary), pointer   :: dict
    character(len=*), intent(in) :: key

    type(LinkedList), pointer   :: elem
    integer                      :: hash

    elem => dict_get_elem( dict, key )

    if ( associated(elem) ) then
       hash = dict_hashkey( trim(key) )
       call list_delete_element( dict%table(hash)%list, elem )
    endif

  end subroutine dict_delete_key

!=====================================================================
! DICT_GET_KEY gets the value belonging to a key
! Arguments:
!     dict       Pointer to the dictionary
!     key        Key for which the values are sought
!=====================================================================

  function dict_get_key( dict, key ) result(value)

    type(Dictionary), pointer   :: dict
    character(len=*), intent(in) :: key
    integer                     :: value

    type(ListData)              :: data
    type(LinkedList), pointer   :: elem

    elem => dict_get_elem( dict, key )

    if ( associated(elem) ) then
       value = elem%data%value
    else
       value = DICT_NULL
    endif

  end function dict_get_key

!=====================================================================
! DICT_HAS_KEY checks if the dictionary has a particular key
! Arguments:
!     dict       Pointer to the dictionary
!     key        Key to be sought
!=====================================================================

  function dict_has_key( dict, key ) result(has)

    type(Dictionary), pointer   :: dict
    character(len=*), intent(in) :: key
    logical                      :: has

    type(LinkedList), pointer   :: elem

    elem => dict_get_elem( dict, key )

    has = associated(elem)

  end function dict_has_key

!=====================================================================
! DICT_GET_ELEM finds the element with a particular key
! Arguments:
!     dict       Pointer to the dictionary
!     key        Key to be sought
!=====================================================================

  function dict_get_elem( dict, key ) result(elem)

    type(Dictionary), pointer   :: dict
    character(len=*), intent(in) :: key

    type(LinkedList), pointer   :: elem
    integer                      :: hash

    hash = dict_hashkey( trim(key) )
    elem => dict%table(hash)%list
    do while ( associated(elem) )
       if ( elem%data%key .eq. key ) then
          exit
       else
          elem => list_next( elem )
       endif
    enddo

  end function dict_get_elem

!=====================================================================
! DICT_HASHKEY determine the hash value from the string
! Arguments:
!     key        String to be examined
!=====================================================================

  integer function dict_hashkey( key )

    character(len=*), intent(in) :: key

    integer                      :: hash
    integer                      :: i

    dict_hashkey = 0

    do i = 1,len(key)
       dict_hashkey = multiplier * dict_hashkey + ichar(key(i:i))
    enddo

    ! Added the absolute value on dict_hashkey-1 since the sum in the
    ! do loop is susceptible to integer overflow
    dict_hashkey = 1 + mod( abs(dict_hashkey-1), hash_size )

  end function dict_hashkey

end module data_structures
