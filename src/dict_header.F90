module dict_header

!===============================================================================
! DICT_HEADER module
!
! This module provides an implementation of a dictionary that has (key,value)
! pairs. This data structure is used to provide lookup features, e.g. cells and
! surfaces by name. As with lists, it was considered writing a single
! dictionary used unlimited polymorphism, but again compiler support is spotty
! and doesn't always prevent duplication of code.
!===============================================================================

  implicit none

  integer, parameter          :: EMPTY           = -huge(0)
  integer, parameter, private :: DELETED         = -huge(0) + 1
  integer, parameter, private :: MIN_SIZE        = 8
  integer, parameter, private :: KEY_CHAR_LENGTH = 255
  real(8), parameter, private :: GROWTH_FACTOR   = 2
  real(8), parameter, private :: MAX_LOAD_FACTOR = 0.65

!===============================================================================
! DICTENTRY* contains (key,value) pairs
!===============================================================================

  type DictEntryCI
    character(len=KEY_CHAR_LENGTH) :: key
    integer                        :: value = EMPTY
    integer                        :: hash = EMPTY
  end type DictEntryCI

  type DictEntryII
    integer :: key = EMPTY
    integer :: value = EMPTY
  end type DictEntryII

!===============================================================================
! DICT* is a dictionary of (key,value) pairs with convenience methods as
! type-bound procedures. DictCharInt has character(*) keys and integer values,
! and DictIntInt has integer keys and values.
!===============================================================================

  type, public :: DictCharInt
    integer, private :: entries = 0
    integer, private :: capacity = 0
    type(DictEntryCI), allocatable, private :: table(:)
  contains
    procedure :: add => add_ci
    procedure :: get => get_ci
    procedure :: has => has_ci
    procedure :: remove => remove_ci
    procedure :: next_entry => next_entry_ci
    procedure :: clear => clear_ci
    procedure :: size => size_ci
    procedure, private :: get_entry => get_entry_ci
    procedure, private :: resize => resize_ci
  end type DictCharInt

  type, public :: DictIntInt
    integer, private :: entries = 0
    integer, private :: capacity = 0
    type(DictEntryII), allocatable, private :: table(:)
  contains
    procedure :: add => add_ii
    procedure :: get => get_ii
    procedure :: has => has_ii
    procedure :: remove => remove_ii
    procedure :: next_entry => next_entry_ii
    procedure :: clear => clear_ii
    procedure :: size => size_ii
    procedure, private :: get_entry => get_entry_ii
    procedure, private :: resize => resize_ii
  end type DictIntInt

contains

!===============================================================================
! GET_ENTRY returns the index of the (key,value) pair in the table for a given
! key. This method is private.
!===============================================================================

  function get_entry_ci(this, key) result(i)

    class(DictCharInt)         :: this
    character(*), intent(in)   :: key
    integer                    :: i

    integer :: hash

    if (.not. allocated(this % table)) then
      allocate(this % table(MIN_SIZE))
      this % capacity = MIN_SIZE
    end if

    hash = hash_ci(key)
    i = 1 + mod(hash, this % capacity)

    do
      if (this % table(i) % hash == hash .and. &
           this % table(i) % key == key) exit

      if (this % table(i) % hash == EMPTY) exit

      ! TODO: use more advanced probing or double hashing to update i
      i = 1 + mod(i, this % capacity)
    end do

  end function get_entry_ci

  function get_entry_ii(this, key) result(i)

    class(DictIntInt)          :: this
    integer, intent(in)        :: key
    integer                    :: i

    integer :: hash

    if (.not. allocated(this % table)) then
      allocate(this % table(MIN_SIZE))
      this % capacity = MIN_SIZE
    end if

    hash = hash_ii(key)
    i = 1 + mod(hash, this % capacity)

    do
      if (this % table(i) % key == key .or. &
           this % table(i) % key == EMPTY) exit

      ! TODO: use more advanced probing or double hashing to update i
      i = 1 + mod(i, this % capacity)
    end do

  end function get_entry_ii

!===============================================================================
! RESIZE allocates a new hash table to accomodate the number of entries and
! reinserts all of the entries into the new table. This method is private.
!===============================================================================

  subroutine resize_ci(this, new_size)

    class(DictCharInt)  :: this
    integer, intent(in) :: new_size

    type(DictEntryCI), allocatable :: table(:)
    integer :: i

    call move_alloc(this % table, table)
    allocate(this % table(new_size))
    this % capacity = new_size
    this % entries = 0

    ! Rehash each entry into the new table
    do i = 1, size(table)
      if (table(i) % hash /= EMPTY .and. table(i) % hash /= DELETED) then
        call this % add(table(i) % key, table(i) % value)
      end if
    end do

    deallocate(table)

  end subroutine resize_ci

  subroutine resize_ii(this, new_size)

    class(DictIntInt)   :: this
    integer, intent(in) :: new_size

    type(DictEntryII), allocatable :: table(:)
    integer :: i

    call move_alloc(this % table, table)
    allocate(this % table(new_size))
    this % capacity = new_size
    this % entries = 0

    ! Rehash each entry into the new table
    do i = 1, size(table)
      if (table(i) % key /= EMPTY .and. table(i) % key /= DELETED) then
        call this % add(table(i) % key, table(i) % value)
      end if
    end do

    deallocate(table)

  end subroutine resize_ii

!===============================================================================
! ADD adds a (key,value) entry to a dictionary. If the key is already in the
! dictionary, the value is replaced by the new specified value.
!===============================================================================

  subroutine add_ci(this, key, value)

    class(DictCharInt)       :: this
    character(*), intent(in) :: key
    integer, intent(in)      :: value

    integer :: hash
    integer :: i

    if (.not. allocated(this % table)) then
      allocate(this % table(MIN_SIZE))
      this % capacity = MIN_SIZE
    else if (real(this % entries, 8) / this % capacity > MAX_LOAD_FACTOR) then
      call this % resize(int(this % capacity * GROWTH_FACTOR))
    end if

    hash = hash_ci(key)
    i = 1 + mod(hash, this % capacity)

    do
      if (this % table(i) % hash == EMPTY .or. &
           this % table(i) % hash == DELETED) then
        this % table(i) % hash = hash
        this % table(i) % key = key
        this % table(i) % value = value
        this % entries = this % entries + 1
        exit
      else if (this % table(i) % hash == hash .and. &
           this % table(i) % key == key) then
        this % table(i) % value = value
        exit
      end if

      ! TODO: use more advanced probing or double hashing to update i
      i = 1 + mod(i, this % capacity)
    end do

  end subroutine add_ci

  subroutine add_ii(this, key, value)

    class(DictIntInt)   :: this
    integer, intent(in) :: key
    integer, intent(in) :: value

    integer :: hash
    integer :: i

    if (.not. allocated(this % table)) then
      allocate(this % table(MIN_SIZE))
      this % capacity = MIN_SIZE
    else if (real(this % entries, 8) / this % capacity > MAX_LOAD_FACTOR) then
      call this % resize(int(this % capacity * GROWTH_FACTOR))
    end if

    hash = hash_ii(key)
    i = 1 + mod(hash, this % capacity)

    do
      if (this % table(i) % key == EMPTY .or. &
           this % table(i) % key == DELETED) then
        this % table(i) % key = key
        this % table(i) % value = value
        this % entries = this % entries + 1
        exit
      else if (this % table(i) % key == key) then
        this % table(i) % value = value
        exit
      end if

      ! TODO: use more advanced probing or double hashing to update i
      i = 1 + mod(i, this % capacity)
    end do

  end subroutine add_ii

!===============================================================================
! GET returns the value matching a given key. If the dictionary does not contain
! the key, the value EMPTY is returned.
!===============================================================================

  function get_ci(this, key) result(value)

    class(DictCharInt)       :: this
    character(*), intent(in) :: key
    integer                  :: value

    integer :: i

    i = this % get_entry(key)
    value = this % table(i) % value

  end function get_ci

  function get_ii(this, key) result(value)

    class(DictIntInt)   :: this
    integer, intent(in) :: key
    integer             :: value

    integer :: i

    i = this % get_entry(key)
    value = this % table(i) % value

  end function get_ii

!===============================================================================
! HAS determines whether a dictionary has a (key,value) pair with a given key.
!===============================================================================

  function has_ci(this, key) result(has)

    class(DictCharInt)       :: this
    character(*), intent(in) :: key
    logical                  :: has

    integer :: i

    i = this % get_entry(key)
    has = (this % table(i) % hash /= EMPTY)

  end function has_ci

  function has_ii(this, key) result(has)

    class(DictIntInt)   :: this
    integer, intent(in) :: key
    logical             :: has

    integer :: i

    i = this % get_entry(key)
    has = (this % table(i) % key /= EMPTY)

  end function has_ii

!===============================================================================
! REMOVE deletes a (key,value) entry from a dictionary.
!===============================================================================

  subroutine remove_ci(this, key)

    class(DictCharInt)       :: this
    character(*), intent(in) :: key

    integer :: i

    i = this % get_entry(key)
    if (this % table(i) % hash /= EMPTY) then
      this % table(i) % hash = DELETED
      this % table(i) % key = ""
      this % table(i) % value = EMPTY
      this % entries = this % entries - 1
    end if

  end subroutine remove_ci

  subroutine remove_ii(this, key)

    class(DictIntInt)   :: this
    integer, intent(in) :: key

    integer :: i

    i = this % get_entry(key)
    if (this % table(i) % key /= EMPTY) then
      this % table(i) % key = DELETED
      this % table(i) % value = EMPTY
      this % entries = this % entries - 1
    end if

  end subroutine remove_ii

!===============================================================================
! NEXT_ENTRY finds the next (key,value) pair. The value of current_entry is
! updated with the (key,value) pair, and the value of i is updated with the
! index of the entry in the table. Passing in i = 0 will locate the first entry
! in the dictionary. If there are no more entries, i will be set to 0. 
!===============================================================================

  subroutine next_entry_ci(this, current_entry, i)

    class(DictCharInt)               :: this
    type(DictEntryCI), intent(inout) :: current_entry
    integer, intent(inout)           :: i

    if (.not. allocated(this % table)) return

    do
      i = i + 1
      if (i > size(this % table)) then
        i = 0
        exit
      else if (this % table(i) % hash /= EMPTY .and. &
           this % table(i) % hash /= DELETED) then
        current_entry = this % table(i)
        exit
      end if
    end do

  end subroutine next_entry_ci

  subroutine next_entry_ii(this, current_entry, i)

    class(DictIntInt)                :: this
    type(DictEntryII), intent(inout) :: current_entry
    integer, intent(inout)           :: i

    if (.not. allocated(this % table)) return

    do
      i = i + 1
      if (i > size(this % table)) then
        i = 0
        exit
      else if (this % table(i) % key /= EMPTY .and. &
           this % table(i) % key /= DELETED) then
        current_entry = this % table(i)
        exit
      end if
    end do

  end subroutine next_entry_ii

!===============================================================================
! CLEAR deletes and deallocates the dictionary item
!===============================================================================

  subroutine clear_ci(this)

    class(DictCharInt) :: this

    if (allocated(this % table)) deallocate(this % table)

  end subroutine clear_ci

  subroutine clear_ii(this)

    class(DictIntInt) :: this

    if (allocated(this % table)) deallocate(this % table)

  end subroutine clear_ii

!===============================================================================
! SIZE returns the number of entries in the dictionary
!===============================================================================

  pure function size_ci(this) result(size)

    class(DictCharInt), intent(in) :: this
    integer :: size

    size = this % entries

  end function size_ci

  pure function size_ii(this) result(size)

    class(DictIntInt), intent(in) :: this
    integer :: size

    size = this % entries

  end function size_ii

!===============================================================================
! HASH returns the hash value for a given key
!===============================================================================

  function hash_ci(key) result(hash)

    character(*), intent(in) :: key
    integer                  :: hash

    integer :: i

    hash = 0

    do i = 1, len_trim(key)
      hash = 31 * hash + ichar(key(i:i))
    end do

    hash = abs(hash - 1)

  end function hash_ci

  function hash_ii(key) result(hash)

    integer, intent(in) :: key
    integer             :: hash

    hash = abs(key - 1)

  end function hash_ii

end module dict_header
