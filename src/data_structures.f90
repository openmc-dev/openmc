module data_structures

!===============================================================================
! DATA_STRUCTURES module
!
! This module implements a dictionary that has (key,value) pairs. This data
! structure is used to provide lookup features, e.g. cells and surfaces by name.
!
! The original version from the 'flibs' open source package implemented another
! derived type called DICT_DATA that has been replaced here by a simple integer
! in ListData. If used in the original form , the dictionary can store derived
! types (changes made to ListData, dict_create, dict_add_key, and dict_get_key).
!===============================================================================

  use types

  implicit none

  integer, parameter :: hash_size  = 4993
  integer, parameter :: multiplier = 31
  integer, parameter :: DICT_NULL = 0

!===============================================================================
! LIST Interfaces -- these allow one to use a single subroutine or function for
! various types of lists, namely those consisting of integers, reals, or
! (key,value) pairs.
!===============================================================================

  interface list_create
     module procedure list_real_create, list_int_create
     module procedure list_kvci_create, list_kvii_create
  end interface
  interface list_insert
     module procedure list_int_insert, list_real_insert
     module procedure list_kvci_insert, list_kvii_insert
  end interface
  interface list_delete
     module procedure list_int_delete, list_real_delete
     module procedure list_kvci_delete, list_kvii_delete
  end interface
  interface list_size
     module procedure list_int_size, list_real_size
     module procedure list_kvci_size, list_kvii_size
  end interface
  interface list_insert_head
     module procedure list_int_insert_head, list_real_insert_head
     module procedure list_kvci_insert_head, list_kvii_insert_head
  end interface
  interface list_delete_element
     module procedure list_int_del_element, list_real_del_element
     module procedure list_kvci_del_element, list_kvii_del_element
  end interface

!===============================================================================
! DICTIONARY Interfaces -- these allow one to use a single subroutine or
! function for dictionaries with varying types of (key,value) pairs
!===============================================================================

  interface dict_create
     module procedure dict_ci_create, dict_ii_create
  end interface
  interface dict_delete
     module procedure dict_ci_delete, dict_ii_delete
  end interface
  interface dict_add_key
     module procedure dict_ci_add_key, dict_ii_add_key
  end interface
  interface dict_delete_key
     module procedure dict_ci_delete_key, dict_ii_delete_key
  end interface
  interface dict_get_key
     module procedure dict_ci_get_key, dict_ii_get_key
  end interface
  interface dict_has_key
     module procedure dict_ci_has_key, dict_ii_has_key
  end interface
  interface dict_get_elem
     module procedure dict_ci_get_elem, dict_ii_get_elem
  end interface
  interface dict_hashkey
     module procedure dict_ci_hashkey, dict_ii_hashkey
  end interface
  interface dict_keys
     module procedure dict_ci_keys, dict_ii_keys
  end interface

contains

!===============================================================================
! LIST_INT_CREATE creates and initializes a list of integers
!===============================================================================

  subroutine list_int_create(list, data)

    type(ListInt), pointer    :: list
    integer,       intent(in) :: data

    allocate(list)
    list%next => null()
    list%data =  data

  end subroutine list_int_create

!===============================================================================
! LIST_REAL_CREATE creates and initializes a list of real(8)s
!===============================================================================

  subroutine list_real_create(list, data)

    type(ListReal), pointer    :: list
    real(8),        intent(in) :: data

    allocate(list)
    list%next => null()
    list%data =  data

  end subroutine list_real_create

!===============================================================================
! LIST_KVCI_CREATE creates and initializes a list of (key,value) pairs
!===============================================================================

  subroutine list_kvci_create(list, data)

    type(ListKeyValueCI), pointer    :: list
    type(KeyValueCI),     intent(in) :: data

    allocate(list)
    list%next => null()
    list%data =  data

  end subroutine list_kvci_create

!===============================================================================
! LIST_KVII_CREATE creates and initializes a list of (key,value) pairs
!===============================================================================

  subroutine list_kvii_create(list, data)

    type(ListKeyValueII), pointer    :: list
    type(KeyValueII),     intent(in) :: data

    allocate(list)
    list%next => null()
    list%data =  data

  end subroutine list_kvii_create

!===============================================================================
! LIST_INT_DELETE deallocates all data in a list of integers
!===============================================================================

  subroutine list_int_delete(list)

    type(ListInt), pointer :: list
    type(ListInt), pointer :: current
    type(ListInt), pointer :: next

    current => list
    do while (associated(current%next))
       next => current%next
       deallocate(current)
       current => next
    enddo

  end subroutine list_int_delete

!===============================================================================
! LIST_REAL_DELETE deallocates all data in a list of real(8)s
!===============================================================================

  subroutine list_real_delete(list)

    type(ListReal), pointer :: list
    type(ListReal), pointer :: current
    type(ListReal), pointer :: next

    current => list
    do while (associated(current%next))
       next => current%next
       deallocate(current)
       current => next
    enddo

  end subroutine list_real_delete

!===============================================================================
! LIST_KVCI_DELETE deallocates all data in a list of (key,value) pairs
!===============================================================================

  subroutine list_kvci_delete(list)

    type(ListKeyValueCI), pointer  :: list
    type(ListKeyValueCI), pointer  :: current
    type(ListKeyValueCI), pointer  :: next

    current => list
    do while (associated(current%next))
       next => current%next
       deallocate(current)
       current => next
    enddo

  end subroutine list_kvci_delete

!===============================================================================
! LIST_KVII_DELETE deallocates all data in a list of (key,value) pairs
!===============================================================================

  subroutine list_kvii_delete(list)

    type(ListKeyValueII), pointer  :: list
    type(ListKeyValueII), pointer  :: current
    type(ListKeyValueII), pointer  :: next

    current => list
    do while (associated(current%next))
       next => current%next
       deallocate(current)
       current => next
    enddo

  end subroutine list_kvii_delete

!===============================================================================
! LIST_INT_SIZE counts the number of items in a list of integers
!===============================================================================

  function list_int_size(list) result(n)

    type(ListInt), pointer  :: list
    type(ListInt), pointer  :: current
    integer                 :: n

    if (associated(list)) then
       n = 1
       current => list
       do while (associated(current%next))
          current => current%next
          n = n + 1
       enddo
    else
       n = 0
    endif

  end function list_int_size

!===============================================================================
! LIST_REAL_SIZE counts the number of items in a list of real(8)s
!===============================================================================

  function list_real_size(list) result(n)

    type(ListReal), pointer  :: list
    type(ListReal), pointer  :: current
    integer                  :: n

    if (associated(list)) then
       n = 1
       current => list
       do while (associated(current%next))
          current => current%next
          n = n + 1
       enddo
    else
       n = 0
    endif

  end function list_real_size

!===============================================================================
! LIST_KVCI_SIZE counts the number of items in a list of (key,value)
! pairs
!===============================================================================

  function list_kvci_size(list) result(n)

    type(ListKeyValueCI), pointer  :: list
    type(ListKeyValueCI), pointer  :: current
    integer                        :: n

    if (associated(list)) then
       n = 1
       current => list
       do while (associated(current%next))
          current => current%next
          n = n + 1
       enddo
    else
       n = 0
    endif

  end function list_kvci_size

!===============================================================================
! LIST_KVII_SIZE counts the number of items in a list of (key,value) pairs
!===============================================================================

  function list_kvii_size(list) result(n)

    type(ListKeyValueII), pointer  :: list
    type(ListKeyValueII), pointer  :: current
    integer                        :: n

    if (associated(list)) then
       n = 1
       current => list
       do while (associated(current%next))
          current => current%next
          n = n + 1
       enddo
    else
       n = 0
    endif

  end function list_kvii_size

!===============================================================================
! LIST_INT_INSERT inserts a new element
! Arguments:
!     elem       Element in the linked list after
!                which to insert the new element
!     data       The data for the new element
!===============================================================================

  subroutine list_int_insert(elem, data)

    type(ListInt), pointer    :: elem
    integer,       intent(in) :: data

    type(ListInt), pointer    :: next

    allocate(next)

    next%next => elem%next
    elem%next => next
    next%data =  data

  end subroutine list_int_insert

!===============================================================================
! LIST_REAL_INSERT inserts a new element
! Arguments:
!     elem       Element in the linked list after
!                which to insert the new element
!     data       The data for the new element
!===============================================================================

  subroutine list_real_insert(elem, data)

    type(ListReal), pointer    :: elem
    real(8),        intent(in) :: data

    type(ListReal), pointer    :: next

    allocate(next)

    next%next => elem%next
    elem%next => next
    next%data =  data

  end subroutine list_real_insert

!===============================================================================
! LIST_KVCI_INSERT inserts a new element
! Arguments:
!     elem       Element in the linked list after
!                which to insert the new element
!     data       The data for the new element
!===============================================================================

  subroutine list_kvci_insert(elem, data)

    type(ListKeyValueCI), pointer    :: elem
    type(KeyValueCI),     intent(in) :: data

    type(ListKeyValueCI), pointer    :: next

    allocate(next)

    next%next => elem%next
    elem%next => next
    next%data =  data

  end subroutine list_kvci_insert

!===============================================================================
! LIST_KVII_INSERT inserts a new element
! Arguments:
!     elem       Element in the linked list after
!                which to insert the new element
!     data       The data for the new element
!===============================================================================

  subroutine list_kvii_insert(elem, data)

    type(ListKeyValueII), pointer    :: elem
    type(KeyValueII),     intent(in) :: data

    type(ListKeyValueII), pointer    :: next

    allocate(next)

    next%next => elem%next
    elem%next => next
    next%data =  data

  end subroutine list_kvii_insert

!===============================================================================
! LIST_INT_INSERT_HEAD inserts a new element before the first element
! Arguments:
!     list       Start of the list
!     data       The data for the new element
!===============================================================================

  subroutine list_int_insert_head(list, data)

    type(ListInt), pointer    :: list
    integer,       intent(in) :: data

    type(ListInt), pointer    :: elem

    allocate(elem)
    elem%data =  data

    elem%next => list
    list      => elem

  end subroutine list_int_insert_head

!===============================================================================
! LIST_REAL_INSERT_HEAD inserts a new element before the first element
! Arguments:
!     list       Start of the list
!     data       The data for the new element
!===============================================================================

  subroutine list_real_insert_head(list, data)

    type(ListReal), pointer    :: list
    real(8),        intent(in) :: data

    type(ListReal), pointer    :: elem

    allocate(elem)
    elem%data =  data

    elem%next => list
    list      => elem

  end subroutine list_real_insert_head

!===============================================================================
! LIST_KVCI_INSERT_HEAD inserts a new element before the first element
! Arguments:
!     list       Start of the list
!     data       The data for the new element
!===============================================================================

  subroutine list_kvci_insert_head(list, data)

    type(ListKeyValueCI), pointer    :: list
    type(KeyValueCI),     intent(in) :: data

    type(ListKeyValueCI), pointer    :: elem

    allocate(elem)
    elem%data =  data

    elem%next => list
    list      => elem

  end subroutine list_kvci_insert_head

!===============================================================================
! LIST_KVII_INSERT_HEAD inserts a new element before the first element
! Arguments:
!     list       Start of the list
!     data       The data for the new element
!===============================================================================

  subroutine list_kvii_insert_head(list, data)

    type(ListKeyValueII), pointer    :: list
    type(KeyValueII),     intent(in) :: data

    type(ListKeyValueII), pointer    :: elem

    allocate(elem)
    elem%data =  data

    elem%next => list
    list      => elem

  end subroutine list_kvii_insert_head

!===============================================================================
! LIST_INT_DEL_ELEMENT deletes an element from the list
! Arguments:
!     list       Header of the list
!     elem       Element in the linked list to be removed
!===============================================================================

  subroutine list_int_del_element(list, elem)

    type(ListInt), pointer  :: list
    type(ListInt), pointer  :: elem

    type(ListInt), pointer  :: current
    type(ListInt), pointer  :: prev

    if (associated(list,elem)) then
       list => elem%next
       deallocate(elem)
    else
       current => list
       prev    => list
       do while (associated(current))
          if (associated(current,elem)) then
             prev%next => current%next
             deallocate(current) ! Is also "elem"
             exit
          endif
          prev    => current
          current => current%next
       enddo
    endif

  end subroutine list_int_del_element

!===============================================================================
! LIST_REAL_DEL_ELEMENT deletes an element from the list
! Arguments:
!     list       Header of the list
!     elem       Element in the linked list to be removed
!===============================================================================

  subroutine list_real_del_element(list, elem)

    type(ListReal), pointer  :: list
    type(ListReal), pointer  :: elem

    type(ListReal), pointer  :: current
    type(ListReal), pointer  :: prev

    if (associated(list,elem)) then
       list => elem%next
       deallocate(elem)
    else
       current => list
       prev    => list
       do while (associated(current))
          if (associated(current,elem)) then
             prev%next => current%next
             deallocate(current) ! Is also "elem"
             exit
          endif
          prev    => current
          current => current%next
       enddo
    endif

  end subroutine list_real_del_element

!===============================================================================
! LIST_KVCI_DEL_ELEMENT deletes an element from the list
! Arguments:
!     list       Header of the list
!     elem       Element in the linked list to be removed
!===============================================================================

  subroutine list_kvci_del_element(list, elem)

    type(ListKeyValueCI), pointer  :: list
    type(ListKeyValueCI), pointer  :: elem

    type(ListKeyValueCI), pointer  :: current
    type(ListKeyValueCI), pointer  :: prev

    if (associated(list,elem)) then
       list => elem%next
       deallocate(elem)
    else
       current => list
       prev    => list
       do while (associated(current))
          if (associated(current,elem)) then
             prev%next => current%next
             deallocate(current) ! Is also "elem"
             exit
          endif
          prev    => current
          current => current%next
       enddo
    endif

  end subroutine list_kvci_del_element

!===============================================================================
! LIST_KVII_DEL_ELEMENT deletes an element from the list
! Arguments:
!     list       Header of the list
!     elem       Element in the linked list to be removed
!===============================================================================

  subroutine list_kvii_del_element(list, elem)

    type(ListKeyValueII), pointer  :: list
    type(ListKeyValueII), pointer  :: elem

    type(ListKeyValueII), pointer  :: current
    type(ListKeyValueII), pointer  :: prev

    if (associated(list,elem)) then
       list => elem%next
       deallocate(elem)
    else
       current => list
       prev    => list
       do while (associated(current))
          if (associated(current,elem)) then
             prev%next => current%next
             deallocate(current) ! Is also "elem"
             exit
          endif
          prev    => current
          current => current%next
       enddo
    endif

  end subroutine list_kvii_del_element

!===============================================================================
! DICT_CI_CREATE creates and initializes a dictionary
!===============================================================================

  subroutine dict_ci_create(dict)

    type(DictionaryCI), pointer   :: dict

    integer :: i

    allocate(dict)
    allocate(dict%table(hash_size))

    do i = 1,hash_size
       dict%table(i)%list => null()
    enddo

  end subroutine dict_ci_create

!===============================================================================
! DICT_II_CREATE creates and initializes a dictionary
!===============================================================================

  subroutine dict_ii_create(dict)

    type(DictionaryII), pointer   :: dict

    integer :: i

    allocate(dict)
    allocate(dict%table(hash_size))

    do i = 1,hash_size
       dict%table(i)%list => null()
    enddo

  end subroutine dict_ii_create

!===============================================================================
! DICT_CI_DELETE destroys an entire dictionary
!===============================================================================

  subroutine dict_ci_delete(dict)

    type(DictionaryCI), pointer  :: dict

    integer :: i

    do i = 1,size(dict%table)
       if (associated(dict%table(i)%list)) then
          call list_delete(dict%table(i)%list)
       endif
    enddo
    deallocate(dict%table)
    deallocate(dict)

  end subroutine dict_ci_delete

!===============================================================================
! DICT_II_DELETE destroys an entire dictionary
!===============================================================================

  subroutine dict_ii_delete(dict)

    type(DictionaryII), pointer  :: dict

    integer :: i

    do i = 1,size(dict%table)
       if (associated(dict%table(i)%list)) then
          call list_delete(dict%table(i)%list)
       endif
    enddo
    deallocate(dict%table)
    deallocate(dict)

  end subroutine dict_ii_delete

!===============================================================================
! DICT_CI_ADD_KEY add a new keys
! Arguments:
!     dict       Pointer to the dictionary
!     key        Key for the new element
!     value      Value for the new element
! Note:
!     If the key already exists, the
!     key's value is simply replaced
!===============================================================================

  subroutine dict_ci_add_key(dict, key, value)

    type(DictionaryCI), pointer    :: dict
    character(*),       intent(in) :: key
    integer,            intent(in) :: value

    type(KeyValueCI)               :: data
    type(ListKeyValueCI), pointer  :: elem
    integer                        :: hash

    elem => dict_get_elem(dict, key)

    if (associated(elem)) then
       elem%data%value = value
    else
       data%key   = key
       data%value = value
       hash       = dict_hashkey(trim(key))
       if (associated(dict%table(hash)%list)) then
          call list_insert(dict%table(hash)%list, data)
       else
          call list_create(dict%table(hash)%list, data)
       endif
    endif

  end subroutine dict_ci_add_key

!===============================================================================
! DICT_II_ADD_KEY add a new keys
! Arguments:
!     dict       Pointer to the dictionary
!     key        Key for the new element
!     value      Value for the new element
! Note:
!     If the key already exists, the
!     key's value is simply replaced
!===============================================================================

  subroutine dict_ii_add_key(dict, key, value)

    type(DictionaryII), pointer :: dict
    integer, intent(in) :: key
    integer, intent(in) :: value

    type(KeyValueII)              :: data
    type(ListKeyValueII), pointer :: elem
    integer                       :: hash

    elem => dict_get_elem(dict, key)

    if (associated(elem)) then
       elem%data%value = value
    else
       data%key   = key
       data%value = value
       hash       = dict_hashkey(key)
       if (associated(dict%table(hash)%list)) then
          call list_insert(dict%table(hash)%list, data)
       else
          call list_create(dict%table(hash)%list, data)
       endif
    endif

  end subroutine dict_ii_add_key

!===============================================================================
! DICT_CI_DELETE_KEY deletes a key-value pair from the dictionary
! Arguments:
!     dict       Dictionary in question
!     key        Key to be removed
!===============================================================================

  subroutine dict_ci_delete_key(dict, key)

    type(DictionaryCI), pointer    :: dict
    character(*),       intent(in) :: key

    type(ListKeyValueCI), pointer  :: elem
    integer                        :: hash

    elem => dict_get_elem(dict, key)

    if (associated(elem)) then
       hash = dict_hashkey(trim(key))
       call list_delete_element(dict%table(hash)%list, elem)
    endif

  end subroutine dict_ci_delete_key

!===============================================================================
! DICT_II_DELETE_KEY deletes a key-value pair from the dictionary
! Arguments:
!     dict       Dictionary in question
!     key        Key to be removed
!===============================================================================

  subroutine dict_ii_delete_key(dict, key)

    type(DictionaryII), pointer :: dict
    integer, intent(in) :: key

    type(ListKeyValueII), pointer :: elem
    integer                       :: hash

    elem => dict_get_elem(dict, key)

    if (associated(elem)) then
       hash = dict_hashkey(key)
       call list_delete_element(dict%table(hash)%list, elem)
    endif

  end subroutine dict_ii_delete_key

!===============================================================================
! DICT_CI_GET_KEY gets the value belonging to a key
! Arguments:
!     dict       Pointer to the dictionary
!     key        Key for which the values are sought
!===============================================================================

  function dict_ci_get_key(dict, key) result(value)

    type(DictionaryCI), pointer    :: dict
    character(*),       intent(in) :: key
    integer                        :: value

    type(KeyValueCI)               :: data
    type(ListKeyValueCI), pointer  :: elem

    elem => dict_get_elem(dict, key)

    if (associated(elem)) then
       value = elem%data%value
    else
       value = DICT_NULL
    endif

  end function dict_ci_get_key

!===============================================================================
! DICT_II_GET_KEY gets the value belonging to a key
! Arguments:
!     dict       Pointer to the dictionary
!     key        Key for which the values are sought
!===============================================================================

  function dict_ii_get_key(dict, key) result(value)

    type(DictionaryII), pointer    :: dict
    integer,            intent(in) :: key
    integer                        :: value

    type(KeyValueII)               :: data
    type(ListKeyValueII), pointer  :: elem

    elem => dict_get_elem(dict, key)

    if (associated(elem)) then
       value = elem%data%value
    else
       value = DICT_NULL
    endif

  end function dict_ii_get_key

!===============================================================================
! DICT_CI_HAS_KEY checks if the dictionary has a particular key
! Arguments:
!     dict       Pointer to the dictionary
!     key        Key to be sought
!===============================================================================

  function dict_ci_has_key(dict, key) result(has)

    type(DictionaryCI), pointer    :: dict
    character(*),       intent(in) :: key
    logical                        :: has

    type(ListKeyValueCI), pointer  :: elem

    elem => dict_get_elem(dict, key)

    has = associated(elem)

  end function dict_ci_has_key

!===============================================================================
! DICT_II_HAS_KEY checks if the dictionary has a particular key
! Arguments:
!     dict       Pointer to the dictionary
!     key        Key to be sought
!===============================================================================

  function dict_ii_has_key(dict, key) result(has)

    type(DictionaryII), pointer    :: dict
    integer,            intent(in) :: key
    logical                        :: has

    type(ListKeyValueII), pointer  :: elem

    elem => dict_get_elem(dict, key)

    has = associated(elem)

  end function dict_ii_has_key

!===============================================================================
! DICT_CI_GET_ELEM finds the element with a particular key
! Arguments:
!     dict       Pointer to the dictionary
!     key        Key to be sought
!===============================================================================

  function dict_ci_get_elem(dict, key) result(elem)

    type(DictionaryCI), pointer    :: dict
    character(*),       intent(in) :: key

    type(ListKeyValueCI), pointer  :: elem
    integer                        :: hash

    hash = dict_hashkey(trim(key))
    elem => dict%table(hash)%list
    do while (associated(elem))
       if (elem%data%key .eq. key) then
          exit
       else
          elem => elem%next
       endif
    enddo

  end function dict_ci_get_elem

!===============================================================================
! DICT_II_GET_ELEM finds the element with a particular key
! Arguments:
!     dict       Pointer to the dictionary
!     key        Key to be sought
!===============================================================================

  function dict_ii_get_elem(dict, key) result(elem)

    type(DictionaryII), pointer    :: dict
    integer,            intent(in) :: key

    type(ListKeyValueII), pointer  :: elem
    integer                        :: hash

    hash = dict_hashkey(key)
    elem => dict%table(hash)%list
    do while (associated(elem))
       if (elem%data%key .eq. key) then
          exit
       else
          elem => elem%next
       endif
    enddo

  end function dict_ii_get_elem

!===============================================================================
! DICT_CI_HASHKEY determine the hash value from the string
! Arguments:
!     key        String to be examined
!===============================================================================

  function dict_ci_hashkey(key) result(val)

    character(*), intent(in) :: key

    integer :: val
    integer :: hash
    integer :: i

    val = 0

    do i = 1,len(key)
       val = multiplier * val + ichar(key(i:i))
    enddo

    ! Added the absolute val on val-1 since the sum in the do loop is
    ! susceptible to integer overflow
    val = 1 + mod(abs(val-1), hash_size)

  end function dict_ci_hashkey

!===============================================================================
! DICT_II_HASHKEY determine the hash value from the string
! Arguments:
!     key        String to be examined
!===============================================================================

  function dict_ii_hashkey(key) result(val)

    integer, intent(in) :: key

    integer :: val
    integer :: hash
    integer :: i

    val = 0

    ! Added the absolute val on val-1 since the sum in the do loop is
    ! susceptible to integer overflow
    val = 1 + mod(abs(key-1), hash_size)

  end function dict_ii_hashkey

!===============================================================================
! DICT_CI_KEYS
!===============================================================================

  function dict_ci_keys(dict) result(head)

    type(DictionaryCI),   pointer :: dict
    type(ListKeyValueCI), pointer :: head
    type(ListKeyValueCI), pointer :: current => null()
    type(ListKeyValueCI), pointer :: elem => null()

    integer :: i, j

    head => null()

    do i = 1, size(dict%table)
       elem => dict%table(i)%list
       do while (associated(elem))
          if (.not. associated(head)) then
             allocate(head)
             current => head
          else
             allocate(current%next)
             current => current%next
          end if
          current%data = elem%data
          elem => elem%next
       end do
    end do

  end function dict_ci_keys

!===============================================================================
! DICT_II_KEYS
!===============================================================================

  function dict_ii_keys(dict) result(head)

    type(DictionaryII),   pointer :: dict
    type(ListKeyValueII), pointer :: head
    type(ListKeyValueII), pointer :: current => null()
    type(ListKeyValueII), pointer :: elem => null()

    integer :: i, j

    head => null()

    do i = 1, size(dict%table)
       elem => dict%table(i)%list
       do while (associated(elem))
          if (.not. associated(head)) then
             allocate(head)
             current => head
          else
             allocate(current%next)
             current => current%next
          end if
          current%data = elem%data
          elem => elem%next
       end do
    end do

  end function dict_ii_keys

end module data_structures
