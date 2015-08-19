module vector_header

  use constants,  only: ZERO

  implicit none
  private

  type, public :: Vector
    integer :: n        ! number of rows/cols in matrix
    real(8), allocatable :: data(:) ! where vector data is stored
    real(8), pointer :: val(:) ! pointer to vector data
  contains
    procedure :: create       => vector_create
    procedure :: destroy      => vector_destroy
    procedure :: add_value    => vector_add_value
    procedure :: copy         => vector_copy
    ! TODO: procedure :: write       => vector_write
  end type Vector

contains

!===============================================================================
! VECTOR_CREATE allocates and initializes a vector
!===============================================================================

  subroutine vector_create(self, n)

    integer, intent(in)                    :: n    ! size of vector
    class(Vector), intent(inout), target   :: self ! vector instance

    ! Preallocate vector
    if (.not.allocated(self % data)) allocate(self % data(n))
    self % val => self % data(1:n)

    ! Set n
    self % n = n

    ! Initialize to zero
    self % val = ZERO

  end subroutine vector_create

!===============================================================================
! VECTOR_DESTROY deallocates all space associated with a vector
!===============================================================================

  subroutine vector_destroy(self)

    class(Vector), intent(inout) :: self ! vector instance

    if (associated(self % val)) nullify(self % val)
    if (allocated(self % data)) deallocate(self % data)

  end subroutine vector_destroy

!===============================================================================
! VECTOR_ADD_VALUE adds a value to the vector
!===============================================================================

  subroutine vector_add_value(self, idx, val)

    integer, intent(in)          :: idx  ! index location in vector
    real(8), intent(in)          :: val  ! value to add
    class(Vector), intent(inout) :: self ! vector instance

    self % val(idx) = val

  end subroutine vector_add_value

!===============================================================================
! VECTOR_COPY allocates a separate vector and copies
!===============================================================================

  subroutine vector_copy(self, vectocopy)

    class(Vector), target, intent(inout) :: self
    type(Vector), intent(in) :: vectocopy

    ! Preallocate vector
    if (.not.allocated(self % data)) allocate(self % data(vectocopy % n))
    self % val => self % data(1:vectocopy % n)

    ! Set n
    self % n = vectocopy % n

    ! Copy values
    self % val = vectocopy % val

  end subroutine vector_copy

end module vector_header
