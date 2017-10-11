module dict_header

!===============================================================================
! DICT_HEADER module
!
! This module provides an implementation of a dictionary that has (key,value)
! pairs. This data structure is used to provide lookup features, e.g. cells and
! surfaces by name.
!
! The implementation is based on Algorithm D from Knuth Vol. 3 Sec. 6.4 (open
! addressing with double hashing). Hash table sizes M are chosen such that M
! and M - 2 are twin primes, which helps reduce clustering. The sequence of
! twin primes used for the table sizes comes from
! https://github.com/anholt/hash_table/blob/master/hash_table.c. These values
! were selected so that the table would grow by approximately a factor of two
! each time the maximum load factor is exceeded. An upper limit is placed on
! the load factor to prevent exponential performance degradation as the number
! of entries approaches the number of buckets.
!===============================================================================

  implicit none

  integer, parameter          :: EMPTY           = -huge(0)
  integer, parameter, private :: DELETED         = -huge(0) + 1
  integer, parameter, private :: KEY_CHAR_LENGTH = 255
  integer, parameter, private :: TABLE_SIZES(30) = &
       [5, 7, 13, 19, 43, 73, 151, 283, 571, 1153, 2269, 4519, 9013, 18043, &
       36109, 72091, 144409, 288361, 576883, 1153459, 2307163, 4613893, &
       9227641, 18455029, 36911011, 73819861, 147639589, 295279081, &
       590559793, 1181116273]
  real(8), parameter, private :: MAX_LOAD_FACTOR = 0.65

!===============================================================================
! DICTENTRY* contains (key,value) pairs.
!===============================================================================

  type DictEntryII
    integer :: key = EMPTY
    integer :: value = EMPTY
  end type DictEntryII

  type DictEntryCI
    character(len=KEY_CHAR_LENGTH) :: key
    integer                        :: value = EMPTY
  end type DictEntryCI

!===============================================================================
! BUCKET* contains an allocatable DictEntry object for storing a (key,value)
! pair and the hash value for fast comparisons. Integer (key,value) pairs are
! stored directly in the hash table since their memory requirement is small.
!===============================================================================

  type, private :: BucketCI
    type(DictEntryCI), allocatable :: entry
    integer                        :: hash = EMPTY
  end type BucketCI

!===============================================================================
! DICT* is a dictionary of (key,value) pairs with convenience methods as
! type-bound procedures. DictIntInt has integer keys and values, and
! DictCharInt has character(*) keys and integer values.
!===============================================================================

  type, public :: DictIntInt
    integer, private :: entries = 0
    integer, private :: capacity = 0
    type(DictEntryII), allocatable, private :: table(:)
  contains
    procedure :: set => set_ii
    procedure :: get => get_ii
    procedure :: has => has_ii
    procedure :: remove => remove_ii
    procedure :: next_entry => next_entry_ii
    procedure :: clear => clear_ii
    procedure :: size => size_ii
    procedure, private :: get_entry => get_entry_ii
    procedure, private :: resize => resize_ii
  end type DictIntInt

  type, public :: DictCharInt
    integer, private :: entries = 0
    integer, private :: capacity = 0
    type(BucketCI), allocatable, private :: table(:)
  contains
    procedure :: set => set_ci
    procedure :: get => get_ci
    procedure :: has => has_ci
    procedure :: remove => remove_ci
    procedure :: next_entry => next_entry_ci
    procedure :: clear => clear_ci
    procedure :: size => size_ci
    procedure, private :: get_entry => get_entry_ci
    procedure, private :: resize => resize_ci
  end type DictCharInt

contains

!===============================================================================
! GET_ENTRY returns the index of the (key,value) pair in the table for a given
! key. This method is private.
!===============================================================================

  function get_entry_ii(this, key) result(i)

    class(DictIntInt)          :: this
    integer, intent(in)        :: key
    integer                    :: i

    integer :: hash
    integer :: c

    if (.not. allocated(this % table)) then
      allocate(this % table(TABLE_SIZES(1)))
      this % capacity = TABLE_SIZES(1)
    end if

    hash = hash_ii(key)
    i = 1 + mod(hash, this % capacity)
    c = 2 + mod(hash, this % capacity - 2)

    do
      if (this % table(i) % key == key .or. &
           this % table(i) % key == EMPTY) exit

      i = 1 + mod(i + c - 1, this % capacity)
    end do

  end function get_entry_ii

  function get_entry_ci(this, key) result(i)

    class(DictCharInt)         :: this
    character(*), intent(in)   :: key
    integer                    :: i

    integer :: hash
    integer :: c

    if (.not. allocated(this % table)) then
      allocate(this % table(TABLE_SIZES(1)))
      this % capacity = TABLE_SIZES(1)
    end if

    hash = hash_ci(key)
    i = 1 + mod(hash, this % capacity)
    c = 2 + mod(hash, this % capacity - 2)

    do
      if (this % table(i) % hash == hash) then
        if (allocated(this % table(i) % entry)) then
          if (this % table(i) % entry % key == key) exit
        end if
      end if

      if (this % table(i) % hash == EMPTY) exit

      i = 1 + mod(i + c - 1, this % capacity)
    end do

  end function get_entry_ci

!===============================================================================
! RESIZE allocates a new hash table to accomodate the number of entries and
! reinserts all of the entries into the new table. This method is private.
!===============================================================================

  subroutine resize_ii(this)

    class(DictIntInt)   :: this

    type(DictEntryII), allocatable :: table(:)
    integer :: new_size
    integer :: i

    do i = 1, size(TABLE_SIZES)
      if (TABLE_SIZES(i) > this % capacity) exit
    end do
    new_size = TABLE_SIZES(i)

    call move_alloc(this % table, table)
    allocate(this % table(new_size))
    this % capacity = new_size
    this % entries = 0

    ! Rehash each entry into the new table
    do i = 1, size(table)
      if (table(i) % key /= EMPTY .and. table(i) % key /= DELETED) then
        call this % set(table(i) % key, table(i) % value)
      end if
    end do

    deallocate(table)

  end subroutine resize_ii

  subroutine resize_ci(this)

    class(DictCharInt)  :: this

    type(BucketCI), allocatable :: table(:)
    integer :: new_size
    integer :: i

    do i = 1, size(TABLE_SIZES)
      if (TABLE_SIZES(i) > this % capacity) exit
    end do
    new_size = TABLE_SIZES(i)

    call move_alloc(this % table, table)
    allocate(this % table(new_size))
    this % capacity = new_size
    this % entries = 0

    ! Rehash each entry into the new table
    do i = 1, size(table)
      if (table(i) % hash /= EMPTY .and. table(i) % hash /= DELETED) then
        call this % set(table(i) % entry % key, table(i) % entry % value)
      end if
    end do

    deallocate(table)

  end subroutine resize_ci

!===============================================================================
! SET adds a (key,value) entry to a dictionary. If the key is already in the
! dictionary, the value is replaced by the new specified value.
!===============================================================================

  subroutine set_ii(this, key, value)

    class(DictIntInt)   :: this
    integer, intent(in) :: key
    integer, intent(in) :: value

    integer :: hash
    integer :: i
    integer :: c

    if (.not. allocated(this % table)) then
      allocate(this % table(TABLE_SIZES(1)))
      this % capacity = TABLE_SIZES(1)
    else if (real(this % entries + 1, 8) / this % capacity > MAX_LOAD_FACTOR) then
      call this % resize()
    end if

    hash = hash_ii(key)
    i = 1 + mod(hash, this % capacity)
    c = 2 + mod(hash, this % capacity - 2)

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

      i = 1 + mod(i + c - 1, this % capacity)
    end do

  end subroutine set_ii

  subroutine set_ci(this, key, value)

    class(DictCharInt)       :: this
    character(*), intent(in) :: key
    integer, intent(in)      :: value

    integer :: hash
    integer :: i
    integer :: c

    if (.not. allocated(this % table)) then
      allocate(this % table(TABLE_SIZES(1)))
      this % capacity = TABLE_SIZES(1)
    else if (real(this % entries + 1, 8) / this % capacity > MAX_LOAD_FACTOR) then
      call this % resize()
    end if

    hash = hash_ci(key)
    i = 1 + mod(hash, this % capacity)
    c = 2 + mod(hash, this % capacity - 2)

    do
      if (this % table(i) % hash == EMPTY .or. &
           this % table(i) % hash == DELETED) then
        if (.not. allocated(this % table(i) % entry)) then
          allocate(this % table(i) % entry)
        end if
        this % table(i) % hash = hash
        this % table(i) % entry % key = key
        this % table(i) % entry % value = value
        this % entries = this % entries + 1
        exit
      else if (this % table(i) % hash == hash .and. &
           this % table(i) % entry % key == key) then
        this % table(i) % entry % value = value
        exit
      end if

      i = 1 + mod(i + c - 1, this % capacity)
    end do

  end subroutine set_ci

!===============================================================================
! GET returns the value matching a given key. If the dictionary does not contain
! the key, the value EMPTY is returned.
!===============================================================================

  function get_ii(this, key) result(value)

    class(DictIntInt)   :: this
    integer, intent(in) :: key
    integer             :: value

    integer :: i

    i = this % get_entry(key)
    value = this % table(i) % value

  end function get_ii

  function get_ci(this, key) result(value)

    class(DictCharInt)       :: this
    character(*), intent(in) :: key
    integer                  :: value

    integer :: i

    i = this % get_entry(key)
    if (allocated(this % table(i) % entry)) then
      value = this % table(i) % entry % value
    else
      value = EMPTY
    end if

  end function get_ci

!===============================================================================
! HAS determines whether a dictionary has a (key,value) pair with a given key.
!===============================================================================

  function has_ii(this, key) result(has)

    class(DictIntInt)   :: this
    integer, intent(in) :: key
    logical             :: has

    integer :: i

    i = this % get_entry(key)
    has = (this % table(i) % key /= EMPTY)

  end function has_ii

  function has_ci(this, key) result(has)

    class(DictCharInt)       :: this
    character(*), intent(in) :: key
    logical                  :: has

    integer :: i

    i = this % get_entry(key)
    has = (this % table(i) % hash /= EMPTY)

  end function has_ci

!===============================================================================
! REMOVE deletes a (key,value) entry from a dictionary.
!===============================================================================

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

  subroutine remove_ci(this, key)

    class(DictCharInt)       :: this
    character(*), intent(in) :: key

    integer :: i

    i = this % get_entry(key)
    if (this % table(i) % hash /= EMPTY) then
      this % table(i) % hash = DELETED
      if (allocated(this % table(i) % entry)) then
        deallocate(this % table(i) % entry)
      end if
      this % entries = this % entries - 1
    end if

  end subroutine remove_ci

!===============================================================================
! NEXT_ENTRY finds the next (key,value) pair. The value of current_entry is
! updated with the (key,value) pair, and the value of i is updated with the
! index of the entry in the table. Passing in i = 0 will locate the first entry
! in the dictionary. If there are no more entries, i will be set to 0.
!===============================================================================

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
        current_entry = this % table(i) % entry
        exit
      end if
    end do

  end subroutine next_entry_ci

!===============================================================================
! CLEAR deletes and deallocates the dictionary item
!===============================================================================

  subroutine clear_ii(this)

    class(DictIntInt) :: this

    if (allocated(this % table)) deallocate(this % table)
    this % entries = 0
    this % capacity = 0

  end subroutine clear_ii

  subroutine clear_ci(this)

    class(DictCharInt) :: this

    if (allocated(this % table)) deallocate(this % table)
    this % entries = 0
    this % capacity = 0

  end subroutine clear_ci

!===============================================================================
! SIZE returns the number of entries in the dictionary
!===============================================================================

  pure function size_ii(this) result(size)

    class(DictIntInt), intent(in) :: this
    integer :: size

    size = this % entries

  end function size_ii

  pure function size_ci(this) result(size)

    class(DictCharInt), intent(in) :: this
    integer :: size

    size = this % entries

  end function size_ci

!===============================================================================
! HASH returns the hash value for a given key
!===============================================================================

  pure function hash_ii(key) result(hash)

    integer, intent(in) :: key
    integer             :: hash

    hash = abs(key - 1)

  end function hash_ii

  pure function hash_ci(key) result(hash)

    character(*), intent(in) :: key
    integer                  :: hash

    integer :: i

    hash = 0

    do i = 1, len_trim(key)
      hash = 31 * hash + ichar(key(i:i))
    end do

    hash = abs(hash - 1)

  end function hash_ci

end module dict_header
