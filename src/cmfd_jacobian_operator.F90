module cmfd_jacobian_operator

# ifdef PETSC

  use cmfd_loss_operator, only: loss_operator,init_M_operator, &
                                build_loss_matrix,destroy_M_operator
  use cmfd_prod_operator, only: prod_operator,init_F_operator, &
                                build_prod_matrix,destroy_F_operator
  implicit none
  private
  public :: init_J_operator, build_jacobian_matrix, destroy_J_operator

# include <finclude/petsc.h90>

  integer  :: nx   ! maximum number of x cells
  integer  :: ny   ! maximum number of y cells
  integer  :: nz   ! maximum number of z cells
  integer  :: ng   ! maximum number of groups
  integer  :: ierr ! petsc error code
    
  type, public :: jacobian_operator
    Mat      :: J    ! petsc matrix for neutronic prod operator
    integer  :: n    ! dimensions of matrix
    integer  :: nnz  ! max number of nonzeros
    integer  :: localn ! local size on proc
    integer, allocatable :: d_nnz(:) ! vector of diagonal preallocation
    integer, allocatable :: o_nnz(:) ! vector of off-diagonal preallocation
  end type jacobian_operator

  type, public :: operators 
    type(loss_operator) :: loss
    type(prod_operator) :: prod
  end type operators

contains

!==============================================================================
! INIT_J_OPERATOR
!==============================================================================

  subroutine init_J_operator(this,ctx)

    type(jacobian_operator) :: this
    type(operators)         :: ctx

    ! get indices
    call get_J_indices(this)

    ! get preallocation
    call preallocate_jacobian_matrix(this,ctx)

    ! set up M operator
    call MatCreateAIJ(PETSC_COMM_WORLD, this%localn, this%localn, PETSC_DECIDE,&
         PETSC_DECIDE, PETSC_NULL_INTEGER, this%d_nnz, PETSC_NULL_INTEGER, &
         this%o_nnz, this%J,ierr)
    call MatSetOption(this%J, MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE, ierr)
    call MatSetOption(this%J, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE, ierr)

  end subroutine init_J_operator

!==============================================================================
! GET_J_INDICES
!==============================================================================

  subroutine get_J_indices(this)

    use global,  only: cmfd, cmfd_coremap

    type(jacobian_operator) :: this

    ! get maximum number of cells in each direction
    nx = cmfd%indices(1)
    ny = cmfd%indices(2)
    nz = cmfd%indices(3)
    ng = cmfd%indices(4)

    ! get number of nonzeros
    this%nnz = 7 + ng - 1

    ! calculate dimensions of matrix
    this%n = nx*ny*nz*ng

    ! calculate dimensions of matrix
    if (cmfd_coremap) then
      this%n = cmfd % mat_dim * ng
    else
      this%n = nx*ny*nz*ng
    end if

    ! add 1 for eigenvalue row
!   this%n = this%n + 1

  end subroutine get_J_indices

!===============================================================================
! PREALLOCATE_JACOBIAN_MATRIX
!===============================================================================

  subroutine preallocate_jacobian_matrix(this,ctx)

    use global,      only: cmfd, n_procs_cmfd, rank

    type(jacobian_operator) :: this
    type(operators)         :: ctx

    integer :: n             ! the extent of the matrix
    integer :: row_start     ! index of local starting row
    integer :: row_end       ! index of local final row

    ! get local problem size
    n = this%n

    ! determine local size, divide evenly between all other procs
    this%localn = n/(n_procs_cmfd)

    ! add 1 more if less proc id is less than mod
    if (rank < mod(n,n_procs_cmfd)) this%localn = this%localn + 1

    ! add another 1 on last proc
    if (rank == n_procs_cmfd - 1) this%localn = this%localn + 1

    ! determine local starting row
    row_start = 0
    if (rank < mod(n,n_procs_cmfd)) then
      row_start = rank*(n/n_procs_cmfd+1)
    else
      row_start = min(mod(n,n_procs_cmfd)*(n/n_procs_cmfd+1) + (rank - &
           mod(n,n_procs_cmfd))*(n/n_procs_cmfd),n)
    end if

    ! determine local final row
    row_end = row_start + this%localn - 1

    ! allocate counters
    if (.not. allocated(this%d_nnz)) allocate(this%d_nnz(row_start:row_end))
    if (.not. allocated(this%o_nnz)) allocate(this%o_nnz(row_start:row_end))
    this % d_nnz = 0
    this % o_nnz = 0

    ! start with pattern from loss matrix
    if (rank == n_procs_cmfd - 1) then
      this % d_nnz(row_start:row_end-1) = ctx%loss%d_nnz
      this % o_nnz(row_start:row_end-1) = ctx%loss%o_nnz
    else
      this % d_nnz = ctx%loss%d_nnz
      this % o_nnz = ctx%loss%o_nnz
    end if

    ! append -F*phi term for last processor will take care of 1 for lambda
    if (rank == n_procs_cmfd - 1) then
      this%d_nnz = this%d_nnz + 1
    else
      this%o_nnz = this%o_nnz + 1
    end if

    ! do last row which has all filled (already did lower left corner above)
    if (rank == n_procs_cmfd - 1) then
      this % d_nnz(row_end) = this % d_nnz(row_end) + (row_end - row_start)
      this % o_nnz(row_end) = this % o_nnz(row_end) + (this%n  - (row_end - &
           row_start))
    end if

  end subroutine preallocate_jacobian_matrix

!===============================================================================
! BUILD_JACOBIAN_MATRIX creates the matrix representing loss of neutrons
!===============================================================================

  subroutine build_jacobian_matrix(snes,x,jac,jac_prec,flag,ctx,ierr)

    use constants,    only: ZERO, ONE
    use global,       only: n_procs_cmfd, cmfd_write_matrices, rank

    SNES            :: snes      ! the snes context
    Vec             :: x         ! the solution vector
    Mat             :: jac       ! the jacobian matrix
    Mat             :: jac_prec  ! the jacobian preconditioner
    MatStructure    :: flag      ! not used
    type(operators) :: ctx       ! not used
    integer         :: ierr      ! petsc error flag

    Vec                  :: phi      ! flux vector
    Vec                  :: source   ! source vector
    integer              :: n        ! problem size
    integer              :: k        ! implied do loop counter
    integer              :: ncols    ! number of nonzeros in cols
    integer              :: irow     ! row counter
    integer              :: row_start! starting local row on process
    integer              :: row_end  ! ending local row on process
    integer, allocatable :: dims(:)  ! vec of starting and ending rows
    integer, allocatable :: dims1(:) ! vec of sizes on each proc
    integer, allocatable :: cols(:)  ! vector of column numbers
    real(8)              :: lambda   ! eigenvalue
    real(8), pointer     :: xptr(:)  ! pointer to solution vector
    real(8), pointer     :: sptr(:)  ! pointer to source vector
    real(8), allocatable :: vals(:)  ! vector of row values
    real(8), allocatable :: phi_tmp(:) ! temp buffer for flux

    ! create operators
    call build_loss_matrix(ctx%loss)
    call build_prod_matrix(ctx%prod)

    ! get problem size
    n = ctx%loss%n

    ! get local size on each processor
    call MatGetOwnershipRange(jac_prec, row_start, row_end, ierr)

    ! allocate cols and  initialize to zero
    if (.not. allocated(cols)) allocate(cols(&
         maxval(ctx%loss%d_nnz + ctx%loss%o_nnz)))
    if (.not. allocated(vals)) allocate(vals(&
         maxval(ctx%loss%d_nnz + ctx%loss%o_nnz)))
    cols = 0
    vals = ZERO

    ! get pointers to residual vector 
    call VecGetArrayF90(x, xptr, ierr)

    ! create petsc vector for flux
    call VecCreateMPI(PETSC_COMM_WORLD, ctx%loss%localn, PETSC_DECIDE, phi, ierr)

    ! extract flux and eigenvalue
    call VecPlaceArray(phi, xptr, ierr)
    if (rank == n_procs_cmfd - 1) lambda = xptr(size(xptr))
    call MPI_BCAST(lambda, 1, MPI_REAL8, n_procs_cmfd-1, PETSC_COMM_WORLD, ierr)

    ! compute math (M-lambda*F) M is overwritten here
    call MatAXPY(ctx%loss%M, -lambda, ctx%prod%F, &
         DIFFERENT_NONZERO_PATTERN, ierr)

    ! create tmp petsc vector for source
    call VecCreateMPI(PETSC_COMM_WORLD, ctx%loss%localn, PETSC_DECIDE, &
         source, ierr)

    ! perform math (-F*phi --> source)
    call MatMult(ctx%prod%F, phi, source, ierr)
    call VecScale(source, -ONE, ierr)

    ! get pointer to source
    call VecGetArrayF90(source, sptr, ierr)

    ! begin loop to insert things into matrix
    do irow = row_start, row_end - 1

      ! don't do last row
      if (irow == n) cycle 

      ! get row of matrix
      call MatGetRow(ctx%loss%M, irow, ncols, cols, vals, ierr)

      ! set that row to Jacobian matrix
      call MatSetValues(jac_prec, 1, (/irow/), ncols, cols(1:ncols), vals, &
           INSERT_VALUES, ierr)

      ! restore the row
      call MatRestoreRow(ctx%loss%M, irow, ncols, cols, vals, ierr)

      ! insert source value
      call MatSetValue(jac_prec, irow, n, sptr(irow-row_start+1), &
           INSERT_VALUES, ierr)

    end do

    ! allocate space for flux vector buffer
    if (rank == n_procs_cmfd - 1) then
      if (.not. allocated(phi_tmp)) allocate(phi_tmp(0:n-1))
    end if

    ! get size on each proc
    if (.not. allocated(dims)) allocate(dims(0:n_procs_cmfd))
    if (.not. allocated(dims1)) allocate(dims1(0:n_procs_cmfd-1))
    call VecGetOwnershipRanges(phi, dims, ierr)
    do k = 0, n_procs_cmfd-1
      dims1(k) = dims(k+1) - dims(k)
    end do

    ! gather data on all procs (will truncate xptr if needed for last proc)
    call MPI_GATHERV(xptr, dims1(rank), MPI_REAL8, phi_tmp, dims1, &
         dims(0:n_procs_cmfd-1), MPI_REAL8, n_procs_cmfd-1, &
         PETSC_COMM_WORLD, ierr)

    ! set values in last row of matrix
    if (rank == n_procs_cmfd - 1) then
      phi_tmp = -phi_tmp ! negate the transpose
      call MatSetValues(jac_prec, 1, (/n/), n, (/(k,k=0,n-1)/), phi_tmp, &
           INSERT_VALUES, ierr)
      call MatSetValue(jac_prec, n, n, ONE, INSERT_VALUES, ierr)
    end if

    ! assemble matrix
    call MatAssemblyBegin(jac_prec, MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyEnd(jac_prec, MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyBegin(jac, MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyEnd(jac, MAT_FINAL_ASSEMBLY, ierr)

    ! reset all vectors
    call VecResetArray(phi,ierr)

    ! restore all vectors
    call VecRestoreArrayF90(x, xptr, ierr)
    call VecRestoreArrayF90(source, sptr, ierr)

    ! destroy all temporary objects
    call VecDestroy(phi, ierr)
    call VecDestroy(source, ierr)

    ! deallocate all temporary space
    if (allocated(cols)) deallocate(cols)
    if (allocated(vals)) deallocate(vals)
    if (allocated(phi_tmp)) deallocate(phi_tmp)
    if (allocated(dims)) deallocate(dims)
    if (allocated(dims1)) deallocate(dims1)

    ! print jacobian out
    if (cmfd_write_matrices) call print_J_operator(jac_prec)

  end subroutine build_jacobian_matrix 

!===============================================================================
! PRINT_J_OPERATOR
!===============================================================================

  subroutine print_J_operator(jac)

    Mat :: jac
 
    PetscViewer :: viewer

    ! write out matrix in binary file (debugging)
    call PetscViewerBinaryOpen(PETSC_COMM_WORLD, 'jacobian.bin', &
         FILE_MODE_WRITE, viewer, ierr)
    call MatView(jac, viewer, ierr)
    call PetscViewerDestroy(viewer, ierr)

  end subroutine print_J_operator

!==============================================================================
! DESTROY_J_OPERATOR
!==============================================================================

  subroutine destroy_J_operator(this)

    type(jacobian_operator) :: this

    ! deallocate matrix
    call MatDestroy(this%J,ierr)

    ! deallocate other parameters
    if (allocated(this%d_nnz)) deallocate(this%d_nnz)
    if (allocated(this%o_nnz)) deallocate(this%o_nnz)

  end subroutine destroy_J_operator

# endif

end module cmfd_jacobian_operator
