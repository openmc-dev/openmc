module cmfd_loss_operator

  use constants,      only: CMFD_NOACCEL, ZERO
  use cmfd_header,    only: cmfd, cmfd_coremap
  use matrix_header,  only: Matrix

  implicit none
  private
  public :: init_loss_matrix, build_loss_matrix

contains

!===============================================================================
! INIT_LOSS_MATRIX preallocates loss matrix and initializes it
!===============================================================================

  subroutine init_loss_matrix(loss_matrix)

    type(Matrix), intent(inout) :: loss_matrix ! cmfd loss matrix

    integer :: nx   ! maximum number of x cells
    integer :: ny   ! maximum number of y cells
    integer :: nz   ! maximum number of z cells
    integer :: ng   ! maximum number of groups
    integer :: n    ! total length of matrix
    integer :: nnz  ! number of nonzeros in matrix
    integer :: n_i  ! number of interior cells
    integer :: n_c  ! number of corner cells
    integer :: n_s  ! number of side cells
    integer :: n_e  ! number of edge cells
    integer :: nz_c ! number of non-zero corner cells
    integer :: nz_e ! number of non-zero edge cells
    integer :: nz_s ! number of non-zero side cells
    integer :: nz_i ! number of non-zero interior cells

    ! Get maximum number of cells in each direction
    nx = cmfd%indices(1)
    ny = cmfd%indices(2)
    nz = cmfd%indices(3)
    ng = cmfd%indices(4)

    ! Calculate dimensions of matrix
    if (cmfd_coremap) then
      n = cmfd % mat_dim * ng
    else
      n = nx*ny*nz*ng
    end if

    ! Calculate number of nonzeros, if core map -> need to determine manually
    if (cmfd_coremap) then
      nnz = preallocate_loss_matrix(nx, ny, nz, ng, n)
    else ! structured Cartesian grid
      n_c  = 8                   ! define # of corners
      n_e  = 4*(nx - 2) + 4*(ny - 2) + 4*(nz - 2) ! define # of edges
      n_s  = 2*(nx - 2)*(ny - 2) + 2*(nx - 2)*(nz - 2) &
           + 2*(ny - 2)*(nz - 2) ! define # of sides
      n_i  = nx*ny*nz - (n_c + n_e + n_s) ! define # of interiors
      nz_c = ng*n_c*(4 + ng - 1) ! define # nonzero corners
      nz_e = ng*n_e*(5 + ng - 1) ! define # nonzero edges
      nz_s = ng*n_s*(6 + ng - 1) ! define # nonzero sides
      nz_i = ng*n_i*(7 + ng - 1) ! define # nonzero interiors
      nnz  = nz_c + nz_e + nz_s + nz_i
    end if

    ! Configure loss matrix
    call loss_matrix % create(n, nnz)

  end subroutine init_loss_matrix

!===============================================================================
! PREALLOCATE_LOSS_MATRIX manually preallocates the loss matrix
!===============================================================================

  function preallocate_loss_matrix(nx, ny, nz, ng, n) result(nnz)

    integer, intent(in) :: nx  ! maximum number of x cells
    integer, intent(in) :: ny  ! maximum number of y cells
    integer, intent(in) :: nz  ! maximum number of z cells
    integer, intent(in) :: ng  ! maximum number of groups
    integer, intent(in) :: n   ! total length of matrix
    integer             :: nnz ! number of nonzeros

    integer :: i             ! iteration counter for x
    integer :: j             ! iteration counter for y
    integer :: k             ! iteration counter for z
    integer :: g             ! iteration counter for groups
    integer :: l             ! iteration counter for leakages
    integer :: h             ! energy group when doing scattering
    integer :: irow          ! row counter
    integer :: bound(6)      ! vector for comparing when looking for bound
    integer :: xyz_idx       ! index for determining if x,y or z leakage
    integer :: dir_idx       ! index for determining - or + face of cell
    integer :: neig_idx(3)   ! spatial indices of neighbour
    integer :: nxyz(3,2)     ! single vector containing bound. locations
    integer :: shift_idx     ! parameter to shift index by +1 or -1
    integer :: neig_mat_idx  ! matrix index of neighbor cell
    integer :: scatt_mat_idx ! matrix index for h-->g scattering terms

    ! Reset number of nonzeros to 0
    nnz = 0

    ! Create single vector of these indices for boundary calculation
    nxyz(1,:) = (/1,nx/)
    nxyz(2,:) = (/1,ny/)
    nxyz(3,:) = (/1,nz/)

    ! Begin loop around local rows
    ROWS: do irow = 1, n

      ! Set a nonzero for diagonal
      nnz = nnz + 1

      ! Get location indices
      call matrix_to_indices(irow, g, i, j, k, ng, nx, ny, nz)

      ! Create boundary vector
      bound = (/i,i,j,j,k,k/)

      ! Begin loop over leakages
      LEAK: do l = 1,6

        ! Define (x,y,z) and (-,+) indices
        xyz_idx = int(ceiling(real(l)/real(2)))  ! x=1, y=2, z=3
        dir_idx = 2 - mod(l,2) ! -=1, +=2

        ! Calculate spatial indices of neighbor
        neig_idx = (/i,j,k/)                ! begin with i,j,k
        shift_idx = -2*mod(l,2) +1          ! shift neig by -1 or +1
        neig_idx(xyz_idx) = shift_idx + neig_idx(xyz_idx)

        ! Check for global boundary
        if (bound(l) /= nxyz(xyz_idx,dir_idx)) then

          ! Check for coremap
          if (cmfd_coremap) then

            ! Check for neighbor that is non-acceleartred
            if (cmfd % coremap(neig_idx(1),neig_idx(2),neig_idx(3)) /= &
                 CMFD_NOACCEL) then

              ! Get neighbor matrix index
              call indices_to_matrix(g,neig_idx(1), neig_idx(2), &
                   neig_idx(3), neig_mat_idx, ng, nx, ny)

              ! Record nonzero
              nnz = nnz + 1

            end if

          else

            ! Get neighbor matrix index
            call indices_to_matrix(g, neig_idx(1), neig_idx(2), neig_idx(3), &
                 neig_mat_idx, ng, nx, ny)

            ! Record nonzero
            nnz = nnz + 1

          end if

        end if

      end do LEAK

      ! Begin loop over off diagonal in-scattering
      SCATTR: do h = 1, ng

        ! Cycle though if h=g, it was already banked in removal xs
        if (h == g) cycle

        ! Get neighbor matrix index
        call indices_to_matrix(h, i, j, k, scatt_mat_idx, ng, nx, ny)

        ! Record nonzero
        nnz = nnz + 1

      end do SCATTR

    end do ROWS

  end function preallocate_loss_matrix

!===============================================================================
! BUILD_LOSS_MATRIX creates the matrix representing loss of neutrons
!===============================================================================

  subroutine build_loss_matrix(loss_matrix, adjoint)

    type(Matrix), intent(inout)   :: loss_matrix ! cmfd loss matrix
    logical, intent(in), optional :: adjoint     ! set up the adjoint

    integer :: nxyz(3,2)     ! single vector containing bound. locations
    integer :: i             ! iteration counter for x
    integer :: j             ! iteration counter for y
    integer :: k             ! iteration counter for z
    integer :: g             ! iteration counter for groups
    integer :: l             ! iteration counter for leakages
    integer :: h             ! energy group when doing scattering
    integer :: nx            ! maximum number of x cells
    integer :: ny            ! maximum number of y cells
    integer :: nz            ! maximum number of z cells
    integer :: ng            ! maximum number of groups
    integer :: neig_mat_idx  ! matrix index of neighbor cell
    integer :: scatt_mat_idx ! matrix index for h-->g scattering terms
    integer :: bound(6)      ! vector for comparing when looking for bound
    integer :: xyz_idx       ! index for determining if x,y or z leakage
    integer :: dir_idx       ! index for determining - or + face of cell
    integer :: neig_idx(3)   ! spatial indices of neighbour
    integer :: shift_idx     ! parameter to shift index by +1 or -1
    integer :: irow          ! iteration counter over row
    logical :: adjoint_calc  ! is this a physical adjoint calculation?
    real(8) :: totxs         ! total macro cross section
    real(8) :: scattxsgg     ! scattering macro cross section g-->g
    real(8) :: scattxshg     ! scattering macro cross section h-->g
    real(8) :: dtilde(6)     ! finite difference coupling parameter
    real(8) :: dhat(6)       ! nonlinear coupling parameter
    real(8) :: hxyz(3)       ! cell lengths in each direction
    real(8) :: jn            ! direction dependent leakage coeff to neig
    real(8) :: jo(6)         ! leakage coeff in front of cell flux
    real(8) :: jnet          ! net leakage from jo
    real(8) :: val           ! temporary variable before saving to

    ! Check for adjoint
    adjoint_calc = .false.
    if (present(adjoint)) adjoint_calc = adjoint

    ! Get maximum number of cells in each direction
    nx = cmfd%indices(1)
    ny = cmfd%indices(2)
    nz = cmfd%indices(3)
    ng = cmfd%indices(4)

    ! Create single vector of these indices for boundary calculation
    nxyz(1,:) = (/1,nx/)
    nxyz(2,:) = (/1,ny/)
    nxyz(3,:) = (/1,nz/)

    ! Begin iteration loops
    ROWS: do irow = 1, loss_matrix % n

      ! Set up a new row in matrix
      call loss_matrix % new_row()

      ! Get indices for that row
      call matrix_to_indices(irow, g, i, j, k, ng, nx, ny, nz)

      ! Retrieve cell data
      totxs = cmfd%totalxs(g,i,j,k)
      scattxsgg = cmfd%scattxs(g,g,i,j,k)
      dtilde = cmfd%dtilde(:,g,i,j,k)
      hxyz = cmfd%hxyz(:,i,j,k)

      ! Check and get dhat
      if (allocated(cmfd%dhat)) then
        dhat = cmfd%dhat(:,g,i,j,k)
      else
        dhat = ZERO
      end if

      ! Create boundary vector
      bound = (/i,i,j,j,k,k/)

      ! Begin loop over leakages
      ! 1=-x, 2=+x, 3=-y, 4=+y, 5=-z, 6=+z
      LEAK: do l = 1,6

        ! Define (x,y,z) and (-,+) indices
        xyz_idx = int(ceiling(real(l)/real(2)))  ! x=1, y=2, z=3
        dir_idx = 2 - mod(l,2) ! -=1, +=2

        ! Calculate spatial indices of neighbor
        neig_idx = (/i,j,k/)                ! begin with i,j,k
        shift_idx = -2*mod(l,2) +1          ! shift neig by -1 or +1
        neig_idx(xyz_idx) = shift_idx + neig_idx(xyz_idx)

        ! Check for global boundary
        if (bound(l) /= nxyz(xyz_idx,dir_idx)) then

          ! Check for core map
          if (cmfd_coremap) then

            ! Check that neighbor is not reflector
            if (cmfd % coremap(neig_idx(1),neig_idx(2),neig_idx(3)) /= &
                 CMFD_NOACCEL) then

              ! Compute leakage coefficient for neighbor
              jn = -dtilde(l) + shift_idx*dhat(l)

              ! Get neighbor matrix index
              call indices_to_matrix(g, neig_idx(1), neig_idx(2), neig_idx(3), &
                   neig_mat_idx, ng, nx, ny)

              ! Compute value and record to bank
              val = jn/hxyz(xyz_idx)

              ! Record value in matrix
              call loss_matrix % add_value(neig_mat_idx, val)

            end if

          else

            ! Compute leakage coefficient for neighbor
            jn = -dtilde(l) + shift_idx*dhat(l)

            ! Get neighbor matrix index
            call indices_to_matrix(g, neig_idx(1), neig_idx(2), neig_idx(3), &
                 neig_mat_idx, ng, nx, ny)

            ! Compute value and record to bank
            val = jn/hxyz(xyz_idx)

            ! Record value in matrix
            call loss_matrix % add_value(neig_mat_idx, val)

          end if

        end if

        ! Compute leakage coefficient for target
        jo(l) = shift_idx*dtilde(l) + dhat(l)

      end do LEAK

      ! Calculate net leakage coefficient for target
      jnet = (jo(2) - jo(1))/hxyz(1) + (jo(4) - jo(3))/hxyz(2) + &
           (jo(6) - jo(5))/hxyz(3)

      ! Calculate loss of neutrons
      val = jnet + totxs - scattxsgg

      ! Record diagonal term
      call loss_matrix % add_value(irow, val)

      ! Begin loop over off diagonal in-scattering
      SCATTR: do h = 1, ng

        ! Cycle though if h=g, value already banked in removal xs
        if (h == g) cycle

        ! Get neighbor matrix index
        call indices_to_matrix(h, i, j, k, scatt_mat_idx, ng, nx, ny)

        ! Check for adjoint
        if (adjoint_calc) then
          ! Get scattering macro xs, transposed!
          scattxshg = cmfd%scattxs(g, h, i, j, k)
        else
          ! Get scattering macro xs
          scattxshg = cmfd%scattxs(h, g, i, j, k)
        end if

        ! Negate the scattering xs
        val = -scattxshg

        ! Record value in matrix
        call loss_matrix % add_value(scatt_mat_idx, val)

      end do SCATTR

    end do ROWS

    ! CSR requires n+1 row
    call loss_matrix % new_row()

  end subroutine build_loss_matrix

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

    ! Check if coremap is used
    if (cmfd_coremap) then

      ! Get idx from core map
      matidx = ng*(cmfd % coremap(i,j,k)) - (ng - g)

    else

      ! Compute index
      matidx = g + ng*(i - 1) + ng*nx*(j - 1) + ng*nx*ny*(k - 1)

    end if

  end subroutine indices_to_matrix

!===============================================================================
! MATRIX_TO_INDICES converts a matrix index to spatial and group indices
!===============================================================================

  subroutine matrix_to_indices(irow, g, i, j, k, ng, nx, ny, nz)

    integer, intent(out) :: i     ! iteration counter for x
    integer, intent(out) :: j     ! iteration counter for y
    integer, intent(out) :: k     ! iteration counter for z
    integer, intent(out) :: g     ! iteration counter for groups
    integer, intent(in)  :: irow  ! iteration counter over row (0 reference)
    integer, intent(in)  :: nx    ! maximum number of x cells
    integer, intent(in)  :: ny    ! maximum number of y cells
    integer, intent(in)  :: nz    ! maximum number of z cells
    integer, intent(in)  :: ng    ! maximum number of groups

    ! Check for core map
    if (cmfd_coremap) then

      ! Get indices from indexmap
      g = mod(irow-1, ng) + 1
      i = cmfd % indexmap((irow-1)/ng+1,1)
      j = cmfd % indexmap((irow-1)/ng+1,2)
      k = cmfd % indexmap((irow-1)/ng+1,3)

    else

      ! Compute indices
      g = mod(irow-1, ng) + 1
      i = mod(irow-1, ng*nx)/ng + 1
      j = mod(irow-1, ng*nx*ny)/(ng*nx)+ 1
      k = mod(irow-1, ng*nx*ny*nz)/(ng*nx*ny) + 1

    end if

  end subroutine matrix_to_indices

end module cmfd_loss_operator
