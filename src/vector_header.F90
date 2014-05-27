module vector_header

  use constants,  only: ZERO

#ifdef PETSC
  use petscvec
#endif

  implicit none
  private

  type, public :: Vector 
    integer :: n        ! number of rows/cols in matrix
    real(8), allocatable :: data(:) ! where vector data is stored
    real(8), pointer :: val(:) ! pointer to vector data
#  ifdef PETSC
    type(vec) :: petsc_vec ! PETSc vector
#  endif
    logical :: petsc_active ! Logical if PETSc is being used
   contains
     procedure :: create       => vector_create
     procedure :: destroy      => vector_destroy
     procedure :: add_value    => vector_add_value
#ifdef PETSC
     procedure :: setup_petsc  => vector_setup_petsc
     procedure :: write_petsc_binary => vector_write_petsc_binary
#endif
  end type Vector

#ifdef PETSC
  integer :: petsc_err ! petsc error code
#endif

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

    ! Petsc is default not active
    self % petsc_active = .false.

  end subroutine vector_create

!===============================================================================
! VECTOR_DESTROY deallocates all space associated with a vector
!===============================================================================

  subroutine vector_destroy(self)

    class(Vector), intent(inout) :: self ! vector instance

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

    integer, intent(in)          :: idx  ! index location in vector
    real(8), intent(in)          :: val  ! value to add
    class(Vector), intent(inout) :: self ! vector instance

    self % val(idx) = val

  end subroutine vector_add_value

!===============================================================================
! VECTOR_SETUP_PETSC links the data to a PETSc vector
!===============================================================================

#ifdef PETSC
  subroutine vector_setup_petsc(self)

    class(Vector), intent(inout) :: self ! vector instance

    ! Link to PETSc 
    call VecCreateSeqWithArray(PETSC_COMM_WORLD, 1, self % n, self % val, &
         self % petsc_vec, petsc_err) 

    ! Set that PETSc is now active
    self % petsc_active = .true.

  end subroutine vector_setup_petsc
#endif

!===============================================================================
! VECTOR_WRITE_PETSC_BINARY writes the PETSc vector to a binary file
!===============================================================================

#ifdef PETSC
  subroutine vector_write_petsc_binary(self, filename)

    character(*), intent(in)  :: filename ! name of file to write to
    class(Vector), intent(in) :: self     ! vector instance

    type(PetscViewer) :: viewer ! PETSc viewer instance

    call PetscViewerBinaryOpen(PETSC_COMM_WORLD, trim(filename), &
         FILE_MODE_WRITE, viewer, petsc_err)
    call VecView(self % petsc_vec, viewer, petsc_err)
    call PetscViewerDestroy(viewer, petsc_err)

  end subroutine vector_write_petsc_binary
#endif

end module vector_header
