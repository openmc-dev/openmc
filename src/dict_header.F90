module dict_header

!===============================================================================
! DICT_HEADER module
!
! This module provides an implementation of a dictionary that has (key,value)
! pairs. This data structure is used to provide lookup features, e.g. cells and
! surfaces by name.
!
! The original version was roughly based on capabilities in the 'flibs' open
! source package. However, it was rewritten from scratch so that it could be
! used stand-alone without relying on the implementation of lists. As with
! lists, it was considered writing a single dictionary used unlimited
! polymorphism, but again compiler support is spotty and doesn't always prevent
! duplication of code.
!===============================================================================

  implicit none

  integer, parameter, private :: HASH_SIZE       = 4993
  integer, parameter, private :: HASH_MULTIPLIER = 31
  integer, parameter, private :: DICT_NULL       = -huge(0)
  integer, parameter          :: DICT_KEY_LENGTH = 255

!===============================================================================
! ELEMKEYVALUE* contains (key,value) pairs and a pointer to the next (key,value)
! pair
!===============================================================================

  type ElemKeyValueCI
    type(ElemKeyValueCI), pointer :: next => null()
    character(len=DICT_KEY_LENGTH) :: key
    integer                        :: value
  end type ElemKeyValueCI

  type ElemKeyValueII
    type(ElemKeyValueII), pointer :: next => null()
    integer :: key
    integer :: value
  end type ElemKeyValueII

!===============================================================================
! HASHLIST* types contain a single pointer to a linked list of (key,value)
! pairs. This type is necesssary so that the Dict types can be dynamically
! allocated.
!===============================================================================

  type, private :: HashListCI
    type(ElemKeyValueCI), pointer :: list => null()
  end type HashListCI

  type, private :: HashListII
    type(ElemKeyValueII), pointer :: list => null()
  end type HashListII

!===============================================================================
! DICT* is a dictionary of (key,value) pairs with convenience methods as
! type-bound procedures. DictCharInt has character(*) keys and integer values,
! and DictIntInt has integer keys and values. 
!===============================================================================

  type, public :: DictCharInt
    private
    type(HashListCI), pointer :: table(:) => null()
  contains
    procedure :: add_key => dict_add_key_ci
    procedure :: get_key => dict_get_key_ci
    procedure :: has_key => dict_has_key_ci
    procedure :: keys => dict_keys_ci
    procedure :: clear => dict_clear_ci
    procedure, private :: get_elem => dict_get_elem_ci
  end type DictCharInt

  type, public :: DictIntInt
    private
    type(HashListII), pointer :: table(:) => null()
  contains
    procedure :: add_key => dict_add_key_ii
    procedure :: get_key => dict_get_key_ii
    procedure :: has_key => dict_has_key_ii
    procedure :: keys => dict_keys_ii
    procedure :: clear => dict_clear_ii
    procedure, private :: get_elem => dict_get_elem_ii
  end type DictIntInt

contains

!===============================================================================
! DICT_ADD_KEY adds a (key,value) entry to a dictionary. If the key is already
! in the dictionary, the value is replaced by the new specified value.
!===============================================================================

  subroutine dict_add_key_ci(this, key, value)

    class(DictCharInt) :: this
    character(*), intent(in) :: key
    integer,      intent(in) :: value

    integer :: hash
    type(ElemKeyValueCI), pointer :: elem => null()
    type(ElemKeyValueCI), pointer :: new_elem => null()

    elem => this % get_elem(key)

    if (associated(elem)) then
      elem % value = value
    else
      ! Get hash 
      hash = dict_hash_key_ci(key)

      ! Create new element
      allocate(new_elem)
      new_elem % key = key
      new_elem % value = value

      ! Add element to front of list
      new_elem % next => this % table(hash) % list
      this % table(hash) % list => new_elem
    end if

  end subroutine dict_add_key_ci

  subroutine dict_add_key_ii(this, key, value)

    class(DictIntInt) :: this
    integer, intent(in) :: key
    integer, intent(in) :: value

    integer :: hash
    type(ElemKeyValueII), pointer :: elem => null()
    type(ElemKeyValueII), pointer :: new_elem => null()

    elem => this % get_elem(key)

    if (associated(elem)) then
      elem % value = value
    else
      ! Get hash 
      hash = dict_hash_key_ii(key)

      ! Create new element
      allocate(new_elem)
      new_elem % key = key
      new_elem % value = value

      ! Add element to front of list
      new_elem % next => this % table(hash) % list
      this % table(hash) % list => new_elem
    end if

  end subroutine dict_add_key_ii

!===============================================================================
! DICT_GET_ELEM returns a pointer to the (key,value) pair for a given key. This
! method is private.
!===============================================================================

  function dict_get_elem_ci(this, key) result(elem)
 
    class(DictCharInt)            :: this
    character(*), intent(in)      :: key
    type(ElemKeyValueCI), pointer :: elem

    integer :: hash
    
    ! Check for dictionary not being allocated
    if (.not. associated(this % table)) then
      allocate(this % table(HASH_SIZE))
    end if

    hash = dict_hash_key_ci(key)
    elem => this % table(hash) % list
    do while (associated(elem))
      if (elem % key == key) exit
      elem => elem % next
    end do

  end function dict_get_elem_ci

  function dict_get_elem_ii(this, key) result(elem)
 
    class(DictIntInt)             :: this
    integer, intent(in)           :: key
    type(ElemKeyValueII), pointer :: elem

    integer :: hash
    
    ! Check for dictionary not being allocated
    if (.not. associated(this % table)) then
      allocate(this % table(HASH_SIZE))
    end if

    hash = dict_hash_key_ii(key)
    elem => this % table(hash) % list
    do while (associated(elem))
      if (elem % key == key) exit
      elem => elem % next
    end do

  end function dict_get_elem_ii

!===============================================================================
! DICT_GET_KEY returns the value matching a given key. If the dictionary does
! not contain the key, the value DICT_NULL is returned.
!===============================================================================

  function dict_get_key_ci(this, key) result(value)

    class(DictCharInt)        :: this
    character(*), intent(in) :: key
    integer                  :: value

    type(ElemKeyValueCI), pointer :: elem

    elem => this % get_elem(key)

    if (associated(elem)) then
      value = elem % value
    else
      value = DICT_NULL
    end if

  end function dict_get_key_ci

  function dict_get_key_ii(this, key) result(value)

    class(DictIntInt)   :: this
    integer, intent(in) :: key
    integer             :: value

    type(ElemKeyValueII), pointer  :: elem

    elem => this % get_elem(key)

    if (associated(elem)) then
      value = elem % value
    else
      value = DICT_NULL
    end if

  end function dict_get_key_ii

!===============================================================================
! DICT_HAS_KEY determines whether a dictionary has a (key,value) pair with a
! given key.
!===============================================================================

  function dict_has_key_ci(this, key) result(has)

    class(DictCharInt)        :: this
    character(*), intent(in) :: key
    logical                  :: has

    type(ElemKeyValueCI), pointer  :: elem

    elem => this % get_elem(key)
    has = associated(elem)

  end function dict_has_key_ci

  function dict_has_key_ii(this, key) result(has)

    class(DictIntInt)   :: this
    integer, intent(in) :: key
    logical             :: has

    type(ElemKeyValueII), pointer  :: elem

    elem => this % get_elem(key)
    has = associated(elem)

  end function dict_has_key_ii

!===============================================================================
! DICT_HASH_KEY returns the hash value for a given key
!===============================================================================

  function dict_hash_key_ci(key) result(val)

    character(*), intent(in) :: key
    integer                  :: val
    
    integer :: i

    val = 0

    do i = 1, len_trim(key)
      val = HASH_MULTIPLIER * val + ichar(key(i:i))
    end do

    ! Added the absolute val on val-1 since the sum in the do loop is
    ! susceptible to integer overflow
    val = 1 + mod(abs(val-1), HASH_SIZE)

  end function dict_hash_key_ci

  function dict_hash_key_ii(key) result(val)

    integer, intent(in) :: key
    integer             :: val
    
    val = 0

    ! Added the absolute val on val-1 since the sum in the do loop is
    ! susceptible to integer overflow
    val = 1 + mod(abs(key-1), HASH_SIZE)

  end function dict_hash_key_ii

!===============================================================================
! DICT_KEYS returns a pointer to a linked list of all (key,value) pairs
!===============================================================================

  function dict_keys_ci(this) result(keys)
    class(DictCharInt)            :: this
    type(ElemKeyValueCI), pointer :: keys

    integer :: i
    type(ElemKeyValueCI), pointer :: current => null()
    type(ElemKeyValueCI), pointer :: elem => null()

    keys => null()

    do i = 1, size(this % table)
      ! Get pointer to start of bucket i
      elem => this % table(i) % list

      do while (associated(elem))
        ! Allocate (key,value) pair
        if (.not. associated(keys)) then
          allocate(keys)
          current => keys
        else
          allocate(current % next)
          current => current % next
        end if

        ! Copy (key,value) pair
        current % key   = elem % key
        current % value = elem % value

        ! Move to next element in bucket i
        elem => elem % next
      end do
    end do

  end function dict_keys_ci

  function dict_keys_ii(this) result(keys)
    class(DictIntInt)             :: this
    type(ElemKeyValueII), pointer :: keys

    integer :: i
    type(ElemKeyValueII), pointer :: current => null()
    type(ElemKeyValueII), pointer :: elem => null()

    keys => null()

    do i = 1, size(this % table)
      ! Get pointer to start of bucket i
      elem => this % table(i) % list

      do while (associated(elem))
        ! Allocate (key,value) pair
        if (.not. associated(keys)) then
          allocate(keys)
          current => keys
        else
          allocate(current % next)
          current => current % next
        end if

        ! Copy (key,value) pair
        current % key   = elem % key
        current % value = elem % value

        ! Move to next element in bucket i
        elem => elem % next
      end do
    end do

  end function dict_keys_ii

!===============================================================================
! DICT_CLEAR Deletes and deallocates the dictionary item
!===============================================================================

  subroutine dict_clear_ci(this)

    class(DictCharInt) :: this

    integer :: i
    type(ElemKeyValueCI), pointer :: current
    type(ElemKeyValueCI), pointer :: next

    if (associated(this % table)) then
      do i = 1, size(this % table)
        current => this % table(i) % list
        do while (associated(current))
          if (associated(current % next)) then
            next => current % next
          else
            nullify(next)
          end if
          deallocate(current)
          current => next
        end do
        if (associated(this % table(i) % list)) &
             nullify(this % table(i) % list)
      end do
      deallocate(this % table)
    end if

  end subroutine dict_clear_ci

  subroutine dict_clear_ii(this)

    class(DictIntInt) :: this

    integer :: i
    type(ElemKeyValueII), pointer :: current
    type(ElemKeyValueII), pointer :: next
    
    if (associated(this % table)) then
      do i = 1, size(this % table)
        current => this % table(i) % list
        do while (associated(current))
          if (associated(current % next)) then
            next => current % next
          else
            nullify(next)
          end if
          deallocate(current)
          current => next
        end do
        if (associated(this % table(i) % list)) &
             nullify(this % table(i) % list)
      end do
      deallocate(this % table)
    end if

  end subroutine dict_clear_ii
  
end module dict_header
