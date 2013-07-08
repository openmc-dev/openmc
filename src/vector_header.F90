module vector_header

  use constants,  only: ZERO

  implicit none
  private

# ifdef PETSC
#  include <finclude/petsc.h90>
# endif

  type, public :: Vector 
    integer :: n        ! number of rows/cols in matrix
    real(8), allocatable :: val(:) ! matrix value vector
#  ifdef PETSC
    Vec :: petsc_vec
#  endif
   contains
     procedure :: create       => vector_create
     procedure :: destroy      => vector_destroy
     procedure :: add_value    => vector_add_value
#   ifdef PETSC
     procedure :: setup_petsc  => vector_setup_petsc
#   endif
  end type Vector

contains

!===============================================================================
! VECTOR_CREATE allocates and initializes a vector 
!===============================================================================

  subroutine vector_create(self, n)

    integer       :: n
    class(Vector) :: self

    ! preallocate vector
    if (.not.allocated(self % val)) allocate(self % val(n))

    ! set n
    self % n = n

    ! initialize to zero
    self % val = ZERO

  end subroutine vector_create

!===============================================================================
! VECTOR_DESTROY deallocates all space associated with a vector
!===============================================================================

  subroutine vector_destroy(self)

    class(Vector) :: self

    if (allocated(self % val)) deallocate(self % val)

  end subroutine vector_destroy

!===============================================================================
! VECTOR_ADD_VALUE adds a value to the vector
!===============================================================================

  subroutine vector_add_value(self, idx, val)

    integer       :: idx
    real(8)       :: val
    class(Vector) :: self

    self % val(idx) = val

  end subroutine vector_add_value

# ifdef PETSC

!===============================================================================
! VECTOR_SETUP_PETSC
!===============================================================================

  subroutine vector_setup_petsc(self)

    class(Vector) :: self

    integer :: petsc_err

    ! link to petsc
    call VecCreateSeqWithArray(PETSC_COMM_WORLD, 1, self % n, self % val, &
         self % petsc_vec, petsc_err) 

  end subroutine vector_setup_petsc

# endif

end module vector_header
