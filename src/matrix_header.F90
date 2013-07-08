module matrix_header

  implicit none
  private

# ifdef PETSC
#  include <finclude/petsc.h90>
# endif

  type, public :: matrix
    integer :: n        ! number of rows/cols in matrix
    integer :: nz       ! number of nonzeros in matrix
    integer :: n_kount  ! counter for length of matrix
    integer :: nz_kount ! counter for number of non zeros
    integer, allocatable :: row(:) ! csr row vector
    integer, allocatable :: col(:) ! column vector
    real(8), allocatable :: val(:) ! matrix value vector
#  ifdef PETSC
    Mat :: petsc_mat
#  endif
   contains
     procedure :: create       => matrix_create
     procedure :: destroy      => matrix_destroy
     procedure :: add_value    => matrix_add_value
     procedure :: new_row      => matrix_new_row
#   ifdef PETSC
     procedure :: setup_petsc  => matrix_setup_petsc
#   endif
  end type matrix

contains

!===============================================================================
! MATRIX_CREATE allocates CSR vectors
!===============================================================================

  subroutine matrix_create(self, n, nz)

    integer       :: n
    integer       :: nz
    class(matrix) :: self

    ! preallocate vectors
    if (.not.allocated(self % row)) allocate(self % row(n+1))
    if (.not.allocated(self % col)) allocate(self % col(nz))
    if (.not.allocated(self % val)) allocate(self % val(nz))

    ! set counters to 1
    self % n_kount  = 1
    self % nz_kount = 1

  end subroutine matrix_create

!===============================================================================
! MATRIX_DESTROY deallocates all space associated with a matrix
!===============================================================================

  subroutine matrix_destroy(self)

    class(matrix) :: self

    if (allocated(self % row)) deallocate(self % row)
    if (allocated(self % col)) deallocate(self % col)
    if (allocated(self % val)) deallocate(self % val)

  end subroutine matrix_destroy

!===============================================================================
! MATRIX_ADD_VALUE adds a value to the matrix
!===============================================================================

  subroutine matrix_add_value(self, col, val)

    integer       :: col
    real(8)       :: val
    class(matrix) :: self

    self % col(self % nz_kount) = col
    self % val(self % nz_kount) = val
    self % nz_kount = self % nz_kount + 1

  end subroutine matrix_add_value

!===============================================================================
! MATRIX_NEW_ROW adds a new row by saving the column position of the new row
!===============================================================================

  subroutine matrix_new_row(self)

    class(matrix) :: self

    self % row(self % n_kount) = self % nz_kount
    self % n_kount = self % n_kount + 1

  end subroutine matrix_new_row

# ifdef PETSC

!===============================================================================
! MATRIX_USE_PETSC
!===============================================================================

  subroutine matrix_use_petsc(self)

    class(matrix) :: self

    integer :: petsc_err

    ! change indices to c notation
    self % row = self % row - 1
    self % col = self % col - 1

    ! link to petsc
    call MatCreateSeqAIJWithArrays(PETSC_COMM_WORLD, self % n, self % n, &
            self % row, self % col, self % val, loss % petsc_mat, petsc_err)

  end subroutine matrix_use_petsc

# endif

end module matrix_header
