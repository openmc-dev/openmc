module cmfd_prod_operator

  use constants,      only: CMFD_NOACCEL
  use cmfd_header,    only: cmfd, cmfd_coremap
  use matrix_header,  only: Matrix

  implicit none
  private
  public :: init_prod_matrix, build_prod_matrix

contains

!==============================================================================
! INIT_PROD_MATRIX preallocates prod matrix and initializes it
!==============================================================================

  subroutine init_prod_matrix(prod_matrix)

    type(Matrix), intent(inout) :: prod_matrix ! production matrix

    integer :: nx   ! maximum number of x cells
    integer :: ny   ! maximum number of y cells
    integer :: nz   ! maximum number of z cells
    integer :: ng   ! maximum number of groups
    integer :: n    ! total length of matrix
    integer :: nnz  ! number of nonzeros in matrix

    ! Get maximum number of cells in each direction
    nx = cmfd%indices(1)
    ny = cmfd%indices(2)
    nz = cmfd%indices(3)
    ng = cmfd%indices(4)

    ! Calculate dimensions and number of nonzeros in matrix
    if (cmfd_coremap) then
      n = cmfd % mat_dim * ng
    else
      n = nx*ny*nz*ng
    end if
    nnz = n * ng

    ! Configure prod matrix
    call prod_matrix % create(n, nnz)

  end subroutine init_prod_matrix

!===============================================================================
! BUILD_PROD_MATRIX creates the matrix representing production of neutrons
!===============================================================================

  subroutine build_prod_matrix(prod_matrix, adjoint)

    type(Matrix), intent(inout)   :: prod_matrix ! production matrix
    logical, intent(in), optional :: adjoint     ! adjoint calculation logical

    integer :: i                  ! iteration counter for x
    integer :: j                  ! iteration counter for y
    integer :: k                  ! iteration counter for z
    integer :: g                  ! iteration counter for groups
    integer :: h                  ! energy group when doing scattering
    integer :: nx                 ! maximum number of x cells
    integer :: ny                 ! maximum number of y cells
    integer :: nz                 ! maximum number of z cells
    integer :: ng                 ! maximum number of groups
    integer :: hmat_idx           ! index in matrix for energy group h
    integer :: irow               ! iteration counter over row
    logical :: adjoint_calc       ! is this a physical adjoint?
    real(8) :: nfissxs            ! nufission cross section h-->g
    real(8) :: val                ! temporary variable for nfissxs

    ! get maximum number of cells in each direction
    nx = cmfd%indices(1)
    ny = cmfd%indices(2)
    nz = cmfd%indices(3)
    ng = cmfd%indices(4)

    ! check for adjoint
    adjoint_calc = .false.
    if (present(adjoint)) adjoint_calc = adjoint

    ! begin iteration loops
    ROWS: do irow = 1, prod_matrix % n

      ! add a new row to matrix
      call prod_matrix % new_row()

      ! get indices for that row
      call matrix_to_indices(irow, g, i, j, k, ng, nx, ny, nz)

      ! check if not including reflector
      if (cmfd_coremap) then

        ! check if at a reflector
        if (cmfd % coremap(i,j,k) == CMFD_NOACCEL) then
          cycle
        end if

      end if

      ! loop around all other groups

      NFISS: do h = 1, ng

        ! get matrix column location
        call indices_to_matrix(h, i, j, k, hmat_idx, ng, nx, ny)

        ! check for adjoint and bank val
        if (adjoint_calc) then
          ! get nu-fission cross section from cell
          nfissxs = cmfd%nfissxs(g,h,i,j,k)
        else
          ! get nu-fission cross section from cell
          nfissxs = cmfd%nfissxs(h,g,i,j,k)
        end if

        ! set as value to be recorded
        val = nfissxs

        ! record value in matrix

        call prod_matrix % add_value(hmat_idx, val)

      end do NFISS

    end do ROWS

    ! CSR requires n+1 row
    call prod_matrix % new_row()

  end subroutine build_prod_matrix

!===============================================================================
! INDICES_TO_MATRIX takes (x,y,z,g) indices and computes location in matrix
!===============================================================================

  subroutine indices_to_matrix(g, i, j, k, matidx, ng, nx, ny)

    integer, intent(out) :: matidx ! the index location in matrix
    integer, intent(in)  :: i      ! current x index
    integer, intent(in)  :: j      ! current y index
    integer, intent(in)  :: k      ! current z index
    integer, intent(in)  :: g      ! current group index
    integer, intent(in)  :: nx     ! maximum number of x cells
    integer, intent(in)  :: ny     ! maximum number of y cells
    integer, intent(in)  :: ng     ! maximum number of groups

    ! check if coremap is used
    if (cmfd_coremap) then

      ! get idx from core map
      matidx = ng*(cmfd % coremap(i,j,k)) - (ng - g)

    else

      ! compute index
      matidx = g + ng*(i - 1) + ng*nx*(j - 1) + ng*nx*ny*(k - 1)

    end if

  end subroutine indices_to_matrix

!===============================================================================
! MATRIX_TO_INDICES converts matrix index to spatial and group indices
!===============================================================================

  subroutine matrix_to_indices(irow, g, i, j, k, ng, nx, ny, nz)

    integer, intent(out) :: i    ! iteration counter for x
    integer, intent(out) :: j    ! iteration counter for y
    integer, intent(out) :: k    ! iteration counter for z
    integer, intent(out) :: g    ! iteration counter for groups
    integer, intent(in)  :: irow ! iteration counter over row (0 reference)
    integer, intent(in)  :: nx   ! maximum number of x cells
    integer, intent(in)  :: ny   ! maximum number of y cells
    integer, intent(in)  :: nz   ! maximum number of z cells
    integer, intent(in)  :: ng   ! maximum number of groups

    ! check for core map
    if (cmfd_coremap) then

      ! get indices from indexmap
      g = mod(irow-1, ng) + 1
      i = cmfd % indexmap((irow-1)/ng+1,1)
      j = cmfd % indexmap((irow-1)/ng+1,2)
      k = cmfd % indexmap((irow-1)/ng+1,3)

    else

      ! compute indices
      g = mod(irow-1, ng) + 1
      i = mod(irow-1, ng*nx)/ng + 1
      j = mod(irow-1, ng*nx*ny)/(ng*nx)+ 1
      k = mod(irow-1, ng*nx*ny*nz)/(ng*nx*ny) + 1

    end if

  end subroutine matrix_to_indices

end module cmfd_prod_operator
