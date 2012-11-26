module datatypes_header

!===============================================================================
! DATATYPES_HEADER module contains the derived types for (key,value) pairs,
! lists, and dictionaries. It would be nice to add sets as well. This module
! could be made a LOT cleaner by using unlimited polymorphism, but we will wait
! on that since compiler support is still limited (namely gfortran which doesn't
! support it yet)
!===============================================================================

  implicit none

!===============================================================================
! KEYVALUECI stores the (key,value) pair for a dictionary where the key is a
! string and the value is an integer. Note that we need to store the key in
! addition to the value for collision resolution.
!===============================================================================

  ! Key length for dictionary
  integer, parameter :: DICT_KEY_LENGTH = 255

  type KeyValueCI
    character(len=DICT_KEY_LENGTH) :: key
    integer                        :: value
  end type KeyValueCI

!===============================================================================
! KEYVALUEII stores the (key,value) pair for a dictionary where the key is an
! integer and the value is an integer. Note that we need to store the key in
! addition to the value for collision resolution.
!===============================================================================

  type KeyValueII
    integer :: key
    integer :: value
  end type KeyValueII

!===============================================================================
! LISTKEYVALUECI stores a linked list of (key,value) pairs where the key is a
! character and the value is an integer
!===============================================================================

  type ListKeyValueCI
    type(ListKeyValueCI), pointer :: next => null()
    type(KeyValueCI)              :: data
  end type ListKeyValueCI

!===============================================================================
! LISTKEYVALUEII stores a linked list of (key,value) pairs where the key is a
! character and the value is an integer
!===============================================================================

  type ListKeyValueII
    type(ListKeyValueII), pointer :: next => null()
    type(KeyValueII)              :: data
  end type ListKeyValueII

!===============================================================================
! LISTREAL stores a linked list of real values. This is used for constructing a
! unionized energy grid.
!===============================================================================

  type ListReal
    type(ListReal), pointer :: next => null()
    real(8)                 :: data
  end type ListReal

!===============================================================================
! LISTINT stores a linked list of integer values.
!===============================================================================

  type ListInt
    type(ListInt), pointer :: next => null()
    integer                :: data
  end type ListInt

!===============================================================================
! HASHLISTCI - Since it's not possible to directly do an array of pointers, this
! derived type provides a pointer
!===============================================================================

  type HashListCI
    type(ListKeyValueCI), pointer :: list => null()
  end type HashListCI

!===============================================================================
! HASHLISTII - Since it's not possible to directly do an array of pointers, this
! derived type provides a pointer
!===============================================================================

  type HashListII
    type(ListKeyValueII), pointer :: list => null()
  end type HashListII

!===============================================================================
! DICTIONARYCI provides a dictionary data structure of (key,value) pairs where
! the keys are strings and values are integers.
!===============================================================================

  type DictionaryCI
    type(HashListCI), pointer :: table(:) => null()
  end type DictionaryCI

!===============================================================================
! DICTIONARYII provides a dictionary data structure of (key,value) pairs where
! the keys and values are both integers.
!===============================================================================

  type DictionaryII
    type(HashListII), pointer :: table(:) => null()
  end type DictionaryII

end module datatypes_header
