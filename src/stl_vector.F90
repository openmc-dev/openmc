module stl_vector

  ! This module provides derived types that are meant to mimic the
  ! std::vector<T> type in C++. The vector type has numerous advantages over
  ! simple arrays and linked lists in that storage can grow and shrink
  ! dynamically, yet it is still contiguous in memory. Vectors can be filled
  ! element-by-element with automatic memory allocation in amortized constant
  ! time. In the implementation here, we grow the vector by a factor of 1.5 each
  ! time the capacity is exceed.
  !
  ! The member functions which have been implemented here are:
  !
  ! capacity -- Returns the size of the storage space currently allocated for
  !             the vector
  !
  ! clear -- Remove all elements from the vector, leaving it with a size of
  !          0. Note that this doesn't imply that storage is deallocated.
  !
  ! initialize -- Set the storage size of the vector and optionally fill it with
  !               a particular value.
  !
  ! pop_back -- Remove the last element of the vector, reducing the size by one.
  !
  ! push_back -- Add a new element at the end of the vector. This increases the
  !              size of the vector by one. Note that the underlying storage is
  !              reallocated only if the size exceeds the capacity.
  !
  ! reserve -- Requests that the capacity of the vector be a certain size.
  !
  ! resize -- Resize the vector so it contains n elements. If n is larger than
  !           the current size, an optional fill value can be used to set the
  !           extra elements.
  !
  ! shrink_to_fit -- Request that the capacity be reduced to fit the size.
  !
  ! size -- Returns the number of elements in the vector.

  implicit none
  private

  integer, parameter :: VECTOR_CHAR_LEN = 255
  real(8), parameter :: GROWTH_FACTOR = 1.5

  type, public :: VectorInt
    integer, private :: size_ = 0
    integer, private :: capacity_ = 0
    integer, allocatable :: data(:)
  contains
    procedure :: capacity => capacity_int
    procedure :: clear => clear_int
    generic :: initialize => &
         initialize_fill_int
    procedure, private :: initialize_fill_int
    procedure :: pop_back => pop_back_int
    procedure :: push_back => push_back_int
    procedure :: reserve => reserve_int
    procedure :: resize => resize_int
    procedure :: shrink_to_fit => shrink_to_fit_int
    procedure :: size => size_int
  end type VectorInt

  type, public :: VectorReal
    integer, private :: size_ = 0
    integer, private :: capacity_ = 0
    real(8), allocatable :: data(:)
  contains
    procedure :: capacity => capacity_real
    procedure :: clear => clear_real
    generic :: initialize => &
         initialize_fill_real
    procedure, private :: initialize_fill_real
    procedure :: pop_back => pop_back_real
    procedure :: push_back => push_back_real
    procedure :: reserve => reserve_real
    procedure :: resize => resize_real
    procedure :: shrink_to_fit => shrink_to_fit_real
    procedure :: size => size_real
  end type VectorReal

  type, public :: VectorChar
    integer, private :: size_ = 0
    integer, private :: capacity_ = 0
    character(VECTOR_CHAR_LEN), allocatable :: data(:)
  contains
    procedure :: capacity => capacity_char
    procedure :: clear => clear_char
    generic :: initialize => &
         initialize_fill_char
    procedure, private :: initialize_fill_char
    procedure :: pop_back => pop_back_char
    procedure :: push_back => push_back_char
    procedure :: reserve => reserve_char
    procedure :: resize => resize_char
    procedure :: shrink_to_fit => shrink_to_fit_char
    procedure :: size => size_char
  end type VectorChar

contains

!===============================================================================
! Implementation of VectorInt
!===============================================================================

  pure function capacity_int(this) result(capacity)
    class(VectorInt), intent(in) :: this
    integer :: capacity

    capacity = this%capacity_
  end function capacity_int

  subroutine clear_int(this)
    class(VectorInt), intent(inout) :: this

    ! Since integer is trivially destructible, we only need to set size to zero
    ! and can leave capacity as is
    this%size_ = 0
    if (allocated(this % data)) then
      this%capacity_ = size(this % data)
    else
      this%capacity_ = 0
    end if
  end subroutine clear_int

  subroutine initialize_fill_int(this, n, val)
    class(VectorInt), intent(inout) :: this
    integer, intent(in) :: n
    integer, optional, intent(in) :: val

    integer :: val_

    ! If no value given, fill the vector with zeros
    if (present(val)) then
      val_ = val
    else
      val_ = 0
    end if

    if (allocated(this%data)) deallocate(this%data)

    allocate(this%data(n), SOURCE=val_)
    this%size_ = n
    this%capacity_ = n
  end subroutine initialize_fill_int

  subroutine pop_back_int(this)
    class(VectorInt), intent(inout) :: this
    if (this%size_ > 0) this%size_ = this%size_ - 1
  end subroutine pop_back_int

  subroutine push_back_int(this, val)
    class(VectorInt), intent(inout) :: this
    integer, intent(in) :: val

    integer :: capacity
    integer, allocatable :: data(:)

    if (this%capacity_ == this%size_) then
      ! Create new data array that is GROWTH_FACTOR larger. Note that
      if (this%capacity_ == 0) then
        capacity = 8
      else
        capacity = int(GROWTH_FACTOR*this%capacity_)
      end if
      allocate(data(capacity))

      ! Copy existing elements
      if (this%size_ > 0) data(1:this%size_) = this%data

      ! Move allocation
      call move_alloc(FROM=data, TO=this%data)
      this%capacity_ = capacity
    end if

    ! Increase size of vector by one and set new element
    this%size_ = this%size_ + 1
    this%data(this%size_) = val
  end subroutine push_back_int

  subroutine reserve_int(this, n)
    class(VectorInt), intent(inout) :: this
    integer, intent(in) :: n

    integer, allocatable :: data(:)

    if (n > this%capacity_) then
      allocate(data(n))

      ! Copy existing elements
      if (this%size_ > 0) data(1:this%size_) = this%data(1:this%size_)

      ! Move allocation
      call move_alloc(FROM=data, TO=this%data)
      this%capacity_ = n
    end if
  end subroutine reserve_int

  subroutine resize_int(this, n, val)
    class(VectorInt), intent(inout) :: this
    integer, intent(in) :: n
    integer, intent(in), optional :: val

    if (n < this%size_) then
      this%size_ = n
    elseif (n > this%size_) then
      ! If requested size is greater than capacity, first reserve that many
      ! elements
      if (n > this%capacity_) call this%reserve(n)

      ! Fill added elements with specified value and increase size
      if (present(val)) this%data(this%size_ + 1 : n) = val
      this%size_ = n
    end if

  end subroutine resize_int

  subroutine shrink_to_fit_int(this)
    class(VectorInt), intent(inout) :: this

    integer, allocatable :: data(:)

    if (this%capacity_ > this%size_) then
      if (this%size_ > 0) then
        allocate(data(this%size_))
        data(:) = this%data(1:this%size_)
        call move_alloc(FROM=data, TO=this%data)
        this%capacity_ = this%size_
      else
        if (allocated(this%data)) deallocate(this%data)
      end if
    end if
  end subroutine shrink_to_fit_int

  pure function size_int(this) result(size)
    class(VectorInt), intent(in) :: this
    integer :: size

    size = this%size_
  end function size_int

!===============================================================================
! Implementation of VectorReal
!===============================================================================

  pure function capacity_real(this) result(capacity)
    class(VectorReal), intent(in) :: this
    integer :: capacity

    capacity = this%capacity_
  end function capacity_real

  subroutine clear_real(this)
    class(VectorReal), intent(inout) :: this

    ! Since real is trivially destructible, we only need to set size to zero and
    ! can leave capacity as is
    this%size_ = 0
    if (allocated(this % data)) then
      this%capacity_ = size(this % data)
    else
      this%capacity_ = 0
    end if
  end subroutine clear_real

  subroutine initialize_fill_real(this, n, val)
    class(VectorReal), intent(inout) :: this
    integer, intent(in) :: n
    real(8), optional, intent(in) :: val

    real(8) :: val_

    ! If no value given, fill the vector with zeros
    if (present(val)) then
      val_ = val
    else
      val_ = 0
    end if

    if (allocated(this%data)) deallocate(this%data)

    allocate(this%data(n), SOURCE=val_)
    this%size_ = n
    this%capacity_ = n
  end subroutine initialize_fill_real

  subroutine pop_back_real(this)
    class(VectorReal), intent(inout) :: this
    if (this%size_ > 0) this%size_ = this%size_ - 1
  end subroutine pop_back_real

  subroutine push_back_real(this, val)
    class(VectorReal), intent(inout) :: this
    real(8), intent(in) :: val

    integer :: capacity
    real(8), allocatable :: data(:)

    if (this%capacity_ == this%size_) then
      ! Create new data array that is GROWTH_FACTOR larger. Note that
      if (this%capacity_ == 0) then
        capacity = 8
      else
        capacity = int(GROWTH_FACTOR*this%capacity_)
      end if
      allocate(data(capacity))

      ! Copy existing elements
      if (this%size_ > 0) data(1:this%size_) = this%data

      ! Move allocation
      call move_alloc(FROM=data, TO=this%data)
      this%capacity_ = capacity
    end if

    ! Increase size of vector by one and set new element
    this%size_ = this%size_ + 1
    this%data(this%size_) = val
  end subroutine push_back_real

  subroutine reserve_real(this, n)
    class(VectorReal), intent(inout) :: this
    integer, intent(in) :: n

    real(8), allocatable :: data(:)

    if (n > this%capacity_) then
      allocate(data(n))

      ! Copy existing elements
      if (this%size_ > 0) data(1:this%size_) = this%data(1:this%size_)

      ! Move allocation
      call move_alloc(FROM=data, TO=this%data)
      this%capacity_ = n
    end if
  end subroutine reserve_real

  subroutine resize_real(this, n, val)
    class(VectorReal), intent(inout) :: this
    integer, intent(in) :: n
    real(8), intent(in), optional :: val

    if (n < this%size_) then
      this%size_ = n
    elseif (n > this%size_) then
      ! If requested size is greater than capacity, first reserve that many
      ! elements
      if (n > this%capacity_) call this%reserve(n)

      ! Fill added elements with specified value and increase size
      if (present(val)) this%data(this%size_ + 1 : n) = val
      this%size_ = n
    end if

  end subroutine resize_real

  subroutine shrink_to_fit_real(this)
    class(VectorReal), intent(inout) :: this

    real(8), allocatable :: data(:)

    if (this%capacity_ > this%size_) then
      if (this%size_ > 0) then
        allocate(data(this%size_))
        data(:) = this%data(1:this%size_)
        call move_alloc(FROM=data, TO=this%data)
        this%capacity_ = this%size_
      else
        if (allocated(this%data)) deallocate(this%data)
      end if
    end if
  end subroutine shrink_to_fit_real

  pure function size_real(this) result(size)
    class(VectorReal), intent(in) :: this
    integer :: size

    size = this%size_
  end function size_real

!===============================================================================
! Implementation of VectorChar
!===============================================================================

  pure function capacity_char(this) result(capacity)
    class(VectorChar), intent(in) :: this
    integer :: capacity

    capacity = this%capacity_
  end function capacity_char

  subroutine clear_char(this)
    class(VectorChar), intent(inout) :: this

    ! Since char is trivially destructible, we only need to set size to zero and
    ! can leave capacity as is
    this%size_ = 0
    if (allocated(this % data)) then
      this%capacity_ = size(this % data)
    else
      this%capacity_ = 0
    end if
  end subroutine clear_char

  subroutine initialize_fill_char(this, n, val)
    class(VectorChar), intent(inout) :: this
    integer, intent(in) :: n
    character(*), optional, intent(in) :: val

    integer :: i
    character(VECTOR_CHAR_LEN) :: val_

    ! If no value given, fill the vector with empty strings
    if (present(val)) then
      val_ = val
    else
      val_ = ''
    end if

    if (allocated(this%data)) deallocate(this%data)

    allocate(this%data(n))
    do i = 1, n
      this%data(i) = val_
    end do
    this%size_ = n
    this%capacity_ = n
  end subroutine initialize_fill_char

  subroutine pop_back_char(this)
    class(VectorChar), intent(inout) :: this
    if (this%size_ > 0) this%size_ = this%size_ - 1
  end subroutine pop_back_char

  subroutine push_back_char(this, val)
    class(VectorChar), intent(inout) :: this
    character(*), intent(in) :: val

    integer :: capacity
    character(VECTOR_CHAR_LEN), allocatable :: data(:)

    if (this%capacity_ == this%size_) then
      ! Create new data array that is GROWTH_FACTOR larger. Note that
      if (this%capacity_ == 0) then
        capacity = 8
      else
        capacity = int(GROWTH_FACTOR*this%capacity_)
      end if
      allocate(data(capacity))

      ! Copy existing elements
      if (this%size_ > 0) data(1:this%size_) = this%data

      ! Move allocation
      call move_alloc(FROM=data, TO=this%data)
      this%capacity_ = capacity
    end if

    ! Increase size of vector by one and set new element
    this%size_ = this%size_ + 1
    this%data(this%size_) = val
  end subroutine push_back_char

  subroutine reserve_char(this, n)
    class(VectorChar), intent(inout) :: this
    integer, intent(in) :: n

    character(VECTOR_CHAR_LEN), allocatable :: data(:)

    if (n > this%capacity_) then
      allocate(data(n))

      ! Copy existing elements
      if (this%size_ > 0) data(1:this%size_) = this%data(1:this%size_)

      ! Move allocation
      call move_alloc(FROM=data, TO=this%data)
      this%capacity_ = n
    end if
  end subroutine reserve_char

  subroutine resize_char(this, n, val)
    class(VectorChar), intent(inout) :: this
    integer, intent(in) :: n
    character(*), intent(in), optional :: val

    if (n < this%size_) then
      this%size_ = n
    elseif (n > this%size_) then
      ! If requested size is greater than capacity, first reserve that many
      ! elements
      if (n > this%capacity_) call this%reserve(n)

      ! Fill added elements with specified value and increase size
      if (present(val)) this%data(this%size_ + 1 : n) = val
      this%size_ = n
    end if

  end subroutine resize_char

  subroutine shrink_to_fit_char(this)
    class(VectorChar), intent(inout) :: this

    character(VECTOR_CHAR_LEN), allocatable :: data(:)

    if (this%capacity_ > this%size_) then
      if (this%size_ > 0) then
        allocate(data(this%size_))
        data(:) = this%data(1:this%size_)
        call move_alloc(FROM=data, TO=this%data)
        this%capacity_ = this%size_
      else
        if (allocated(this%data)) deallocate(this%data)
      end if
    end if
  end subroutine shrink_to_fit_char

  pure function size_char(this) result(size)
    class(VectorChar), intent(in) :: this
    integer :: size

    size = this%size_
  end function size_char

end module stl_vector
