module matrix_header

  implicit none
  private

# ifdef PETSC
#  include <finclude/petsc.h90>
# endif

  type, public :: Matrix
    integer :: n        ! number of rows/cols in matrix
    integer :: nnz      ! number of nonzeros in matrix
    integer :: n_kount  ! counter for length of matrix
    integer :: nz_kount ! counter for number of non zeros
    integer, allocatable :: row(:) ! csr row vector
    integer, allocatable :: col(:) ! column vector
    real(8), allocatable :: val(:) ! matrix value vector
#  ifdef PETSC
    Mat :: petsc_mat
#  endif
    logical :: petsc_active
   contains
     procedure :: create       => matrix_create
     procedure :: destroy      => matrix_destroy
     procedure :: add_value    => matrix_add_value
     procedure :: new_row      => matrix_new_row
     procedure :: assemble     => matrix_assemble
     procedure :: setup_petsc        => matrix_setup_petsc
     procedure :: write_petsc_binary => matrix_write_petsc_binary
     procedure :: transpose          => matrix_transpose
     procedure :: vector_multiply    => matrix_vector_multiply
  end type matrix

  integer :: petsc_err

contains

!===============================================================================
! MATRIX_CREATE allocates CSR vectors
!===============================================================================

  subroutine matrix_create(self, n, nnz)

    integer       :: n
    integer       :: nnz
    class(Matrix) :: self

    ! Preallocate vectors
    if (.not.allocated(self % row)) allocate(self % row(n+1))
    if (.not.allocated(self % col)) allocate(self % col(nnz))
    if (.not.allocated(self % val)) allocate(self % val(nnz))

    ! Set counters to 1
    self % n_kount  = 1
    self % nz_kount = 1

    ! Set n and nnz
    self % n = n
    self % nnz = nnz

    ! Set petsc active by default to false
    self % petsc_active = .false.

  end subroutine matrix_create

!===============================================================================
! MATRIX_DESTROY deallocates all space associated with a matrix
!===============================================================================

  subroutine matrix_destroy(self)

    class(Matrix) :: self

#ifdef PETSC
    if (self % petsc_active) call MatDestroy(self % petsc_mat, petsc_err)
#endif

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
    class(Matrix) :: self

    self % col(self % nz_kount) = col
    self % val(self % nz_kount) = val
    self % nz_kount = self % nz_kount + 1

  end subroutine matrix_add_value

!===============================================================================
! MATRIX_NEW_ROW adds a new row by saving the column position of the new row
!===============================================================================

  subroutine matrix_new_row(self)

    class(Matrix) :: self

    self % row(self % n_kount) = self % nz_kount
    self % n_kount = self % n_kount + 1

  end subroutine matrix_new_row

!===============================================================================
! MATRIX_ASSEMBLE main rountine to sort all the columns in CSR matrix
!===============================================================================

  subroutine matrix_assemble(self)

    class(Matrix) :: self

    integer :: i
    integer :: first
    integer :: last

    ! Loop around row vector
    do i = 1, self % n

      ! Get bounds
      first = self % row(i)
      last  = self % row(i+1) - 1

      ! Sort a row
      call sort_csr(self % col, self % val, first, last)

    end do

  end subroutine matrix_assemble

!===============================================================================
! SORT_CSR main routine that performs a sort on a CSR col/val subvector
!===============================================================================

  recursive subroutine sort_csr(col, val, first, last)

    integer :: col(:)
    integer :: first
    integer :: last
    real(8) :: val(:)

    integer :: mid

    if (first < last) then
      call split(col, val, first, last, mid)
      call sort_csr(col, val, first, mid-1)
      call sort_csr(col, val, mid+1, last)
    end if

  end subroutine sort_csr

!===============================================================================
! SPLIT bisects the search space for the sorting routine
!===============================================================================

  subroutine split(col, val, low, high, mid)

    integer :: col(:)
    integer :: low
    integer :: high
    integer :: mid
    real(8) :: val(:)

    integer :: left
    integer :: right
    integer :: iswap
    integer :: pivot
    real(8) :: rswap
    real(8) :: val0

    left = low
    right = high
    pivot = col(low)
    val0 = val(low)

    ! Repeat the following while left and right havent met
    do while (left < right)

      ! Scan right to left to find element < pivot
      do while (left < right .and. col(right) >= pivot)
        right = right - 1
      end do

      ! Scan left to right to find element > pivot
      do while (left < right .and. col(left) <= pivot)
        left = left + 1
      end do

      ! If left and right havent met, exchange the items
      if (left < right) then
        iswap = col(left)
        col(left) = col(right)
        col(right) = iswap
        rswap = val(left)
        val(left) = val(right)
        val(right) = rswap
      end if

    end do

    ! Switch the element in split position with pivot
    col(low) = col(right)
    col(right) = pivot
    mid = right
    val(low) = val(right)
    val(right) = val0

  end subroutine split

!===============================================================================
! MATRIX_SETUP_PETSC configures the row/col vectors and links to a petsc object
!===============================================================================

  subroutine matrix_setup_petsc(self)

    class(Matrix) :: self

    ! change indices to c notation
    self % row = self % row - 1
    self % col = self % col - 1

    ! Link to petsc
#ifdef PETSC
    call MatCreateSeqAIJWithArrays(PETSC_COMM_WORLD, self % n, self % n, &
            self % row, self % col, self % val, self % petsc_mat, petsc_err)
#endif

    ! Petsc is now active
    self % petsc_active = .true.

  end subroutine matrix_setup_petsc

!===============================================================================
! MATRIX_WRITE_PETSC_BINARY writes a petsc matrix binary file
!===============================================================================

  subroutine matrix_write_petsc_binary(self, filename)

    character(*) :: filename
    class(Matrix) :: self

#ifdef PETSC
    PetscViewer :: viewer

    call PetscViewerBinaryOpen(PETSC_COMM_WORLD, trim(filename), &
         FILE_MODE_WRITE, viewer, petsc_err)
    call MatView(self % petsc_mat, viewer, petsc_err)
    call PetscViewerDestroy(viewer, petsc_err)
#endif

  end subroutine matrix_write_petsc_binary

!===============================================================================
! MATRIX_TRANSPOSE uses PETSc to transpose a matrix
!===============================================================================

  subroutine matrix_transpose(self)

    class(Matrix) :: self

#ifdef PETSC
    call MatTranspose(self % petsc_mat, MAT_REUSE_MATRIX, self % petsc_mat, &
         petsc_err)
#endif

  end subroutine matrix_transpose

!===============================================================================
! MATRIX_VECTOR_MULTIPLY allow a vector to multiply the matrix
!===============================================================================

  subroutine matrix_vector_multiply(self, vec_in, vec_out)

    use constants,      only: ZERO
    use vector_header,  only: Vector

    class(Matrix) :: self
    type(Vector)  :: vec_in
    type(Vector)  :: vec_out

    integer :: i
    integer :: j
    integer :: shift

    ! Set shift by default 0 and change if petsc is active
    ! This is because PETSc needs vectors to remain referenced to 0
    shift = 0
    if (self % petsc_active) shift = 1

    ! Begin loop around rows
    ROWS: do i = 1, self % n

      ! Initialize target location in vector
      vec_out % val(i) = ZERO

      ! Begin loop around columns (need to shift if petsc is active)
      COLS: do j = self % row(i) + shift, self % row(i + 1) - 1 + shift

        vec_out % val(i) = vec_out % val(i) + self % val(j) * &
                           vec_in % val(self % col(j) + shift)

      end do COLS

    end do ROWS

  end subroutine matrix_vector_multiply

end module matrix_header
