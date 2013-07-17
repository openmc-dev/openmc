module vector_header

  use constants,  only: ZERO

  implicit none
  private

# ifdef PETSC
#  include <finclude/petsc.h90>
# endif

  type, public :: Vector 
    integer :: n        ! number of rows/cols in matrix
    real(8), allocatable :: data(:) ! where vector data is stored
    real(8), pointer :: val(:) ! pointer to vector data
#  ifdef PETSC
    Vec :: petsc_vec
#  endif
    logical :: petsc_active
   contains
     procedure :: create       => vector_create
     procedure :: destroy      => vector_destroy
     procedure :: add_value    => vector_add_value
     procedure :: setup_petsc  => vector_setup_petsc
     procedure :: write_petsc_binary => vector_write_petsc_binary
  end type Vector

  integer :: petsc_err

contains

!===============================================================================
! VECTOR_CREATE allocates and initializes a vector 
!===============================================================================

  subroutine vector_create(self, n)

    integer       :: n
    class(Vector), target :: self

    ! Preallocate vector
    if (.not.allocated(self % data)) allocate(self % data(n))
    self % val => self % data(1:n)

    ! Set n
    self % n = n

    ! Initialize to zero
    self % val = ZERO

    ! Petsc is default not active
    self % petsc_active = .false.

  end subroutine vector_create

!===============================================================================
! VECTOR_DESTROY deallocates all space associated with a vector
!===============================================================================

  subroutine vector_destroy(self)

    class(Vector) :: self

#ifdef PETSC
    if (self % petsc_active) call VecDestroy(self % petsc_vec, petsc_err)
#endif

    if (associated(self % val)) nullify(self % val)
    if (allocated(self % data)) deallocate(self % data)

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

!===============================================================================
! VECTOR_SETUP_PETSC links the data to a PETSc vector
!===============================================================================

  subroutine vector_setup_petsc(self)

    class(Vector) :: self

    ! Link to PETSc 
#ifdef PETSC
    call VecCreateSeqWithArray(PETSC_COMM_WORLD, 1, self % n, self % val, &
         self % petsc_vec, petsc_err) 
#endif

    ! Set that PETSc is now active
    self % petsc_active = .true.

  end subroutine vector_setup_petsc

!===============================================================================
! VECTOR_WRITE_PETSC_BINARY writes the PETSc vector to a binary file
!===============================================================================

  subroutine vector_write_petsc_binary(self, filename)

    character(*) :: filename
    class(Vector) :: self

#ifdef PETSC
    PetscViewer :: viewer

    call PetscViewerBinaryOpen(PETSC_COMM_WORLD, trim(filename), &
         FILE_MODE_WRITE, viewer, petsc_err)
    call VecView(self % petsc_vec, viewer, petsc_err)
    call PetscViewerDestroy(viewer, petsc_err)
#endif

  end subroutine vector_write_petsc_binary

end module vector_header
