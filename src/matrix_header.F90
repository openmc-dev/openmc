module matrix_header

#ifdef PETSC
  use petscmat
#endif

  implicit none
  private

  type, public :: Matrix
    integer :: n        ! number of rows/cols in matrix
    integer :: nnz      ! number of nonzeros in matrix
    integer :: n_count  ! counter for length of matrix
    integer :: nz_count ! counter for number of non zeros
    integer, allocatable :: row(:) ! csr row vector
    integer, allocatable :: col(:) ! column vector
    real(8), allocatable :: val(:) ! matrix value vector
#  ifdef PETSC
    type(mat) :: petsc_mat
#  endif
    logical :: petsc_active
   contains
     procedure :: create       => matrix_create
     procedure :: destroy      => matrix_destroy
     procedure :: add_value    => matrix_add_value
     procedure :: new_row      => matrix_new_row
     procedure :: assemble     => matrix_assemble
     procedure :: get_row      => matrix_get_row
     procedure :: get_col      => matrix_get_col
     procedure :: vector_multiply    => matrix_vector_multiply
#ifdef PETSC
     procedure :: transpose          => matrix_transpose
     procedure :: setup_petsc        => matrix_setup_petsc
     procedure :: write_petsc_binary => matrix_write_petsc_binary
#endif
  end type matrix

#ifdef PETSC
  integer :: petsc_err
#endif

contains

!===============================================================================
! MATRIX_CREATE allocates CSR vectors
!===============================================================================

  subroutine matrix_create(self, n, nnz)

    integer, intent(in) :: n             ! dimension of matrix
    integer, intent(in) :: nnz           ! number of nonzeros
    class(Matrix), intent(inout) :: self ! matrix instance

    ! Preallocate vectors
    if (.not.allocated(self % row)) allocate(self % row(n+1))
    if (.not.allocated(self % col)) allocate(self % col(nnz))
    if (.not.allocated(self % val)) allocate(self % val(nnz))

    ! Set counters to 1
    self % n_count  = 1
    self % nz_count = 1

    ! Set n and nnz
    self % n = n
    self % nnz = nnz

    ! Set PETSc active by default to false
    self % petsc_active = .false.

  end subroutine matrix_create

!===============================================================================
! MATRIX_DESTROY deallocates all space associated with a matrix
!===============================================================================

  subroutine matrix_destroy(self)

    class(Matrix), intent(inout) :: self ! matrix instance

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

    integer, intent(in) :: col           ! col location in matrix
    real(8), intent(in) :: val           ! value to store in matrix
    class(Matrix), intent(inout) :: self ! matrix instance

    ! Record the data
    self % col(self % nz_count) = col
    self % val(self % nz_count) = val

    ! Need to adjust column indices if PETSc is active
    if (self % petsc_active) self % col(self % nz_count) = &
                             self % col(self % nz_count) - 1

    ! Increment the number of nonzeros currently stored
    self % nz_count = self % nz_count + 1

  end subroutine matrix_add_value

!===============================================================================
! MATRIX_NEW_ROW adds a new row by saving the column position of the new row
!===============================================================================

  subroutine matrix_new_row(self)

    class(Matrix), intent(inout) :: self ! matrix instance

    ! Record the current number of nonzeros
    self % row(self % n_count) = self % nz_count

    ! If PETSc is active, we have to reference indices off 0
    if (self % petsc_active) self % row(self % n_count) = &
                             self % row(self % n_count) - 1

    ! Increment the current row that we are on
    self % n_count = self % n_count + 1

  end subroutine matrix_new_row

!===============================================================================
! MATRIX_ASSEMBLE main rountine to sort all the columns in CSR matrix
!===============================================================================

  subroutine matrix_assemble(self)

    class(Matrix), intent(inout) :: self ! matrix instance

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

    integer :: col(:) ! column vector to sort
    integer :: first  ! first value in sort
    integer :: last   ! last value in sort
    real(8) :: val(:) ! value vector to be sorted like columns

    integer :: mid ! midpoint value

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

    integer :: col(:) ! column vector to sort
    integer :: low    ! low index to sort
    integer :: high   ! high index to sort
    integer :: mid    ! middle of sort
    real(8) :: val(:) ! value vector to be sorted like columns

    integer :: left  ! contains left value in sort
    integer :: right ! contains right value in sort
    integer :: iswap ! temporary interger swap
    integer :: pivot ! pivotting variable for columns
    real(8) :: rswap ! temporary real swap
    real(8) :: val0  ! pivot for value vector

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
! MATRIX_GET_ROW is a method to get row and checks for PETSc C indexing use
!===============================================================================

  function matrix_get_row(self, i) result(row)

    class(Matrix), intent(in) :: self ! matrix instance
    integer, intent(in)       :: i    ! row to get
    integer                   :: row  ! index in col where row begins

    row = self % row(i)

    if (self % petsc_active) row = row + 1

  end function matrix_get_row

!===============================================================================
! MATRIX_GET_COL is a method to get column and checks for PETSc C index use
!===============================================================================

  function matrix_get_col(self, i) result(col)

    class(Matrix), intent(in) :: self ! matrix instance
    integer, intent(in)       :: i    ! index from row vector
    integer                   :: col  ! the actual column

    col = self % col(i)

    if (self % petsc_active) col = col + 1

  end function matrix_get_col

!===============================================================================
! MATRIX_SETUP_PETSC configures the row/col vectors and links to a PETSc object
!===============================================================================

#ifdef PETSC
  subroutine matrix_setup_petsc(self)

    class(Matrix), intent(inout) :: self ! matrix instance

    ! change indices to c notation
    self % row = self % row - 1
    self % col = self % col - 1

    ! Link to petsc
    call MatCreateSeqAIJWithArrays(PETSC_COMM_WORLD, self % n, self % n, &
            self % row, self % col, self % val, self % petsc_mat, petsc_err)

    ! Petsc is now active
    self % petsc_active = .true.

  end subroutine matrix_setup_petsc
#endif

!===============================================================================
! MATRIX_WRITE_PETSC_BINARY writes a PETSc matrix binary file
!===============================================================================

#ifdef PETSC
  subroutine matrix_write_petsc_binary(self, filename)

    character(*), intent(in)  :: filename ! file name to write to
    class(Matrix), intent(in) :: self     ! matrix instance

    type(PetscViewer) :: viewer ! a petsc viewer instance

    call PetscViewerBinaryOpen(PETSC_COMM_WORLD, trim(filename), &
         FILE_MODE_WRITE, viewer, petsc_err)
    call MatView(self % petsc_mat, viewer, petsc_err)
    call PetscViewerDestroy(viewer, petsc_err)

  end subroutine matrix_write_petsc_binary
#endif

!===============================================================================
! MATRIX_TRANSPOSE uses PETSc to transpose a matrix
!===============================================================================

#ifdef PETSC
  subroutine matrix_transpose(self)

    class(Matrix), intent(inout) :: self ! matrix instance

    call MatTranspose(self % petsc_mat, MAT_REUSE_MATRIX, self % petsc_mat, &
         petsc_err)

  end subroutine matrix_transpose
#endif

!===============================================================================
! MATRIX_VECTOR_MULTIPLY allow a vector to multiply the matrix
!===============================================================================

  subroutine matrix_vector_multiply(self, vec_in, vec_out)

    use constants,      only: ZERO
    use vector_header,  only: Vector

    class(Matrix), intent(in)   :: self    ! matrix instance
    type(Vector), intent(in)    :: vec_in  ! vector to multiply matrix against
    type(Vector), intent(inout) :: vec_out ! resulting vector

    integer :: i ! row iteration counter
    integer :: j ! column iteration counter

    ! Begin loop around rows
    ROWS: do i = 1, self % n

      ! Initialize target location in vector
      vec_out % val(i) = ZERO

      ! Begin loop around columns
      COLS: do j = self % get_row(i), self % get_row(i + 1) - 1

        vec_out % val(i) = vec_out % val(i) + self % val(j) * &
                           vec_in % val(self % get_col(j))

      end do COLS

    end do ROWS

  end subroutine matrix_vector_multiply

end module matrix_header
