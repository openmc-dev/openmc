module matrix_header

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
  contains
    procedure :: create             => matrix_create
    procedure :: destroy            => matrix_destroy
    procedure :: add_value          => matrix_add_value
    procedure :: new_row            => matrix_new_row
    procedure :: assemble           => matrix_assemble
    procedure :: get_row            => matrix_get_row
    procedure :: get_col            => matrix_get_col
    procedure :: vector_multiply    => matrix_vector_multiply
    procedure :: search_indices     => matrix_search_indices
    procedure :: write              => matrix_write
    procedure :: copy               => matrix_copy
    procedure :: transpose          => matrix_transpose
  end type matrix

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

  end subroutine matrix_create

!===============================================================================
! MATRIX_DESTROY deallocates all space associated with a matrix
!===============================================================================

  subroutine matrix_destroy(self)

    class(Matrix), intent(inout) :: self ! matrix instance

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
! MATRIX_GET_ROW is a method to get row
!===============================================================================

  function matrix_get_row(self, i) result(row)

    class(Matrix), intent(in) :: self ! matrix instance
    integer, intent(in)       :: i    ! row to get
    integer                   :: row  ! index in col where row begins

    row = self % row(i)

  end function matrix_get_row

!===============================================================================
! MATRIX_GET_COL is a method to get column
!===============================================================================

  function matrix_get_col(self, i) result(col)

    class(Matrix), intent(in) :: self ! matrix instance
    integer, intent(in)       :: i    ! index from row vector
    integer                   :: col  ! the actual column

    col = self % col(i)

  end function matrix_get_col

!===============================================================================
! MATRIX_TRANSPOSE transposes a sparse matrix
!===============================================================================

  function matrix_transpose(self) result(mat)

    class(Matrix), intent(in) :: self ! matrix instance
    type(Matrix)              :: mat  ! transposed matrix

    integer :: i, j        ! loop indices for row/column
    integer :: i_transpose ! loop index for row in transposed matrix
    integer :: first, last ! indices in col() array

    call mat % create(self % n, self % nnz)
    do i_transpose = 1, mat % n
      ! Set up row in transposed matrix
      call mat % new_row()

      ROW: do i = 1, self % n
        ! Get range of columns for row i
        first = self % row(i)
        last  = self % row(i+1) - 1

        COL: do j = first, last
          if (self % col(j) == i_transpose) then
            ! If column in original matrix matches row in transposed, add value
            call mat % add_value(i, self % val(j))
          elseif (self % col(j) > i_transpose) then
            exit COL
          end if
        end do COL
      end do ROW
    end do

    call mat % new_row()

  end function matrix_transpose

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

!===============================================================================
! MATRIX_SEARCH_INDICES searches for an index in column corresponding to a row
!===============================================================================

  subroutine matrix_search_indices(self, row, col, idx, found)

    class(Matrix), intent(inout) :: self
    integer, intent(in) :: row
    integer, intent(in) :: col
    integer, intent(out) :: idx
    logical, intent(out) :: found

    integer :: j

    found = .false.

    COLS: do j = self % get_row(row), self % get_row(row + 1) - 1

      if (self % get_col(j) == col) then
        idx = j
        found = .true.
        exit
      end if

    end do COLS

  end subroutine matrix_search_indices

!===============================================================================
! MATRIX_WRITE writes a matrix to file
!===============================================================================

  subroutine matrix_write(self, filename)

    character(*), intent(in) :: filename
    class(Matrix), intent(inout) :: self

    integer :: unit_
    integer :: i
    integer :: j

    open(newunit=unit_, file=filename)

    do i = 1, self % n
      do j = self % get_row(i), self % get_row(i + 1) - 1
        write(unit_,*) i, self % get_col(j), self % val(j)
      end do
    end do

    close(unit_)

  end subroutine matrix_write

!===============================================================================
! MATRIX_COPY copies a matrix
!===============================================================================

  subroutine matrix_copy(self, mattocopy)

    class(Matrix), intent(inout) :: self
    type(Matrix), intent(in) :: mattocopy

    ! Set n and nnz
    self % n_count = mattocopy % n_count
    self % nz_count = mattocopy % nz_count
    self % n = mattocopy % n
    self % nnz = mattocopy % nnz

    ! Allocate vectors
    if (.not.allocated(self % row)) allocate(self % row(self % n + 1))
    if (.not.allocated(self % col)) allocate(self % col(self % nnz))
    if (.not.allocated(self % val)) allocate(self % val(self % nnz))

    ! Copy over data
    self % row = mattocopy % row
    self % col = mattocopy % col
    self % val = mattocopy % val

  end subroutine matrix_copy

end module matrix_header
