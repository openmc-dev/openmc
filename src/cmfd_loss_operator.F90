module cmfd_loss_operator

  use constants,      only: CMFD_NOACCEL, ZERO
  use global,         only: cmfd, cmfd_coremap
  use matrix_header,  only: Matrix

  implicit none
  private
  public :: init_loss_matrix, build_loss_matrix

contains

!===============================================================================
! INIT_LOSS_MATRIX
!===============================================================================

  subroutine init_loss_matrix(loss_matrix)

    type(Matrix) :: loss_matrix ! cmfd loss matrix

    integer :: nx   ! maximum number of x cells
    integer :: ny   ! maximum number of y cells
    integer :: nz   ! maximum number of z cells
    integer :: ng   ! maximum number of groups
    integer :: n    ! total length of matrix
    integer :: nnz  ! number of nonzeros in matrix
    integer :: n_i  ! number of interior cells
    integer :: n_c  ! number of corner cells
    integer :: n_s  ! number side cells
    integer :: nz_c ! number of non-zero corner cells
    integer :: nz_s ! number of non-zero side cells
    integer :: nz_i ! number of non-zero interior cells

    ! get maximum number of cells in each direction
    nx = cmfd%indices(1)
    ny = cmfd%indices(2)
    nz = cmfd%indices(3)
    ng = cmfd%indices(4)

    ! calculate dimensions of matrix
    if (cmfd_coremap) then
      n = cmfd % mat_dim * ng
    else
      n = nx*ny*nz*ng
    end if

    ! calculate number of nonzeros, if core map -> need to determine manually
    if (cmfd_coremap) then
      nnz = preallocate_loss_matrix(nx, ny, nz, ng, n)
    else ! structured Cartesian grid
      n_c  = 4                   ! define # of corners
      n_s  = 2*(nx + ny) - 8     ! define # of sides
      n_i  = nx*ny - (n_c + n_S) ! define # of interiors  
      nz_c = ng*n_c*(3 + ng - 1) ! define # nonzero corners
      nz_s = ng*n_s*(4 + ng - 1) ! define # nonzero sides
      nz_i = ng*n_i*(5 + ng - 1) ! define # nonzero interiors
      nnz  = nz_c + nz_s + nz_i
    end if

    ! configure loss matrix
    call loss_matrix % create(n, nnz)

  end subroutine init_loss_matrix

!===============================================================================
! PREALLOCATE_LOSS_MATRIX
!===============================================================================

  function preallocate_loss_matrix(nx, ny, nz, ng, n) result(nnz)

    integer, intent(in)    :: nx   ! maximum number of x cells
    integer, intent(in)    :: ny   ! maximum number of y cells
    integer, intent(in)    :: nz   ! maximum number of z cells
    integer, intent(in)    :: ng   ! maximum number of groups
    integer, intent(in)    :: n    ! total length of matrix
    integer                :: nnz  ! number of nonzeros

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

    ! reset to 0
    nnz = 0

    ! create single vector of these indices for boundary calculation
    nxyz(1,:) = (/1,nx/)
    nxyz(2,:) = (/1,ny/)
    nxyz(3,:) = (/1,nz/)

    ! begin loop around local rows
    ROWS: do irow = 1, n

      ! set a nonzero for diagonal
      nnz = nnz + 1

      ! get location indices
      call matrix_to_indices(irow, g, i, j, k, ng, nx, ny, nz)

      ! create boundary vector
      bound = (/i,i,j,j,k,k/)

      ! begin loop over leakages
      LEAK: do l = 1,6

        ! define (x,y,z) and (-,+) indices
        xyz_idx = int(ceiling(real(l)/real(2)))  ! x=1, y=2, z=3
        dir_idx = 2 - mod(l,2) ! -=1, +=2

        ! calculate spatial indices of neighbor
        neig_idx = (/i,j,k/)                ! begin with i,j,k
        shift_idx = -2*mod(l,2) +1          ! shift neig by -1 or +1
        neig_idx(xyz_idx) = shift_idx + neig_idx(xyz_idx)

        ! check for global boundary
        if (bound(l) /= nxyz(xyz_idx,dir_idx)) then

          ! check for coremap 
          if (cmfd_coremap) then

            ! check for neighbor that is non-acceleartred
            if (cmfd % coremap(neig_idx(1),neig_idx(2),neig_idx(3)) /= &
                 CMFD_NOACCEL) then

              ! get neighbor matrix index
              call indices_to_matrix(g,neig_idx(1), neig_idx(2), & 
                   neig_idx(3), neig_mat_idx, ng, nx, ny)

              ! record nonzero
              nnz = nnz + 1 

            end if

          else

            ! get neighbor matrix index
            call indices_to_matrix(g, neig_idx(1), neig_idx(2), neig_idx(3), &
                 neig_mat_idx, ng, nx, ny)

            ! record nonzero
            nnz = nnz + 1

          end if

        end if

      end do LEAK

      ! begin loop over off diagonal in-scattering
      SCATTR: do h = 1, ng

        ! cycle though if h=g
        if (h == g) cycle

        ! get neighbor matrix index
        call indices_to_matrix(h, i, j, k, scatt_mat_idx, ng, nx, ny)

        ! record nonzero
        nnz = nnz + 1

      end do SCATTR

    end do ROWS

  end function preallocate_loss_matrix

!===============================================================================
! BUILD_LOSS_MATRIX creates the matrix representing loss of neutrons
!===============================================================================

  subroutine build_loss_matrix(loss_matrix, adjoint)

    type(Matrix) :: loss_matrix     ! cmfd loss matrix
    logical, optional :: adjoint    ! set up the adjoint

    integer :: nxyz(3,2)            ! single vector containing bound. locations
    integer :: i                    ! iteration counter for x
    integer :: j                    ! iteration counter for y
    integer :: k                    ! iteration counter for z
    integer :: g                    ! iteration counter for groups
    integer :: l                    ! iteration counter for leakages
    integer :: h                    ! energy group when doing scattering
    integer :: nx                   ! maximum number of x cells
    integer :: ny                   ! maximum number of y cells
    integer :: nz                   ! maximum number of z cells
    integer :: ng                   ! maximum number of groups
    integer :: neig_mat_idx         ! matrix index of neighbor cell
    integer :: scatt_mat_idx        ! matrix index for h-->g scattering terms
    integer :: bound(6)             ! vector for comparing when looking for bound
    integer :: xyz_idx              ! index for determining if x,y or z leakage
    integer :: dir_idx              ! index for determining - or + face of cell
    integer :: neig_idx(3)          ! spatial indices of neighbour
    integer :: shift_idx            ! parameter to shift index by +1 or -1
    integer :: irow                 ! iteration counter over row
    logical :: adjoint_calc         ! is this a physical adjoint calculation?
    real(8) :: totxs                ! total macro cross section
    real(8) :: scattxsgg            ! scattering macro cross section g-->g
    real(8) :: scattxshg            ! scattering macro cross section h-->g
    real(8) :: dtilde(6)            ! finite difference coupling parameter
    real(8) :: dhat(6)              ! nonlinear coupling parameter
    real(8) :: hxyz(3)              ! cell lengths in each direction
    real(8) :: jn                   ! direction dependent leakage coeff to neig
    real(8) :: jo(6)                ! leakage coeff in front of cell flux
    real(8) :: jnet                 ! net leakage from jo
    real(8) :: val                  ! temporary variable before saving to 

    ! check for adjoint
    adjoint_calc = .false.
    if (present(adjoint)) adjoint_calc = adjoint 

    ! get maximum number of cells in each direction
    nx = cmfd%indices(1)
    ny = cmfd%indices(2)
    nz = cmfd%indices(3)
    ng = cmfd%indices(4)

    ! create single vector of these indices for boundary calculation
    nxyz(1,:) = (/1,nx/)
    nxyz(2,:) = (/1,ny/)
    nxyz(3,:) = (/1,nz/)

    ! begin iteration loops
    ROWS: do irow = 1, loss_matrix % n

      ! set up a new row in matrix
      call loss_matrix % new_row

      ! get indices for that row
      call matrix_to_indices(irow, g, i, j, k, ng, nx, ny, nz)

      ! retrieve cell data
      totxs = cmfd%totalxs(g,i,j,k)
      scattxsgg = cmfd%scattxs(g,g,i,j,k)
      dtilde = cmfd%dtilde(:,g,i,j,k)
      hxyz = cmfd%hxyz(:,i,j,k)

      ! check and get dhat
      if (allocated(cmfd%dhat)) then
        dhat = cmfd%dhat(:,g,i,j,k)
      else
        dhat = ZERO
      end if

      ! create boundary vector 
      bound = (/i,i,j,j,k,k/)

      ! begin loop over leakages
      ! 1=-x, 2=+x, 3=-y, 4=+y, 5=-z, 6=+z 
      LEAK: do l = 1,6

        ! define (x,y,z) and (-,+) indices
        xyz_idx = int(ceiling(real(l)/real(2)))  ! x=1, y=2, z=3
        dir_idx = 2 - mod(l,2) ! -=1, +=2

        ! calculate spatial indices of neighbor
        neig_idx = (/i,j,k/)                ! begin with i,j,k
        shift_idx = -2*mod(l,2) +1          ! shift neig by -1 or +1
        neig_idx(xyz_idx) = shift_idx + neig_idx(xyz_idx)

        ! check for global boundary
        if (bound(l) /= nxyz(xyz_idx,dir_idx)) then

          ! check for core map
          if (cmfd_coremap) then

            ! check that neighbor is not reflector
            if (cmfd % coremap(neig_idx(1),neig_idx(2),neig_idx(3)) /= &
                 CMFD_NOACCEL) then

              ! compute leakage coefficient for neighbor
              jn = -dtilde(l) + shift_idx*dhat(l)

              ! get neighbor matrix index
              call indices_to_matrix(g, neig_idx(1), neig_idx(2), neig_idx(3), &
                   neig_mat_idx, ng, nx, ny)

              ! compute value and record to bank
              val = jn/hxyz(xyz_idx)

              ! record value in matrix
              call loss_matrix % add_value(neig_mat_idx, val)

            end if

          else

            ! compute leakage coefficient for neighbor
            jn = -dtilde(l) + shift_idx*dhat(l)

            ! get neighbor matrix index
            call indices_to_matrix(g, neig_idx(1), neig_idx(2), neig_idx(3), &
                 neig_mat_idx, ng, nx, ny)

            ! compute value and record to bank
            val = jn/hxyz(xyz_idx)

            ! record value in matrix
            call loss_matrix % add_value(neig_mat_idx, val)

          end if

        end if

        ! compute leakage coefficient for target
        jo(l) = shift_idx*dtilde(l) + dhat(l)

      end do LEAK

      ! calate net leakage coefficient for target
      jnet = (jo(2) - jo(1))/hxyz(1) + (jo(4) - jo(3))/hxyz(2) + &
           (jo(6) - jo(5))/hxyz(3)

      ! calculate loss of neutrons
      val = jnet + totxs - scattxsgg

      ! record diagonal term
      call loss_matrix % add_value(irow, val)

      ! begin loop over off diagonal in-scattering
      SCATTR: do h = 1, ng

        ! cycle though if h=g
        if (h == g) cycle

        ! get neighbor matrix index
        call indices_to_matrix(h, i, j, k, scatt_mat_idx, ng, nx, ny)

        ! get scattering macro xs
        scattxshg = cmfd%scattxs(h, g, i, j, k)

        ! record value in matrix (negate it)
        val = -scattxshg

        ! check for adjoint and bank value
        if (adjoint_calc) then

          ! get scattering macro xs, transposed!
          scattxshg = cmfd%scattxs(g, h, i, j, k)

          ! negate the scattering xs 
          val = -scattxshg

          ! record value in matrix
          call loss_matrix % add_value(scatt_mat_idx, val)

        else

          ! get scattering macro xs
          scattxshg = cmfd%scattxs(h, g, i, j, k)

          ! negate the scattering xs 
          val = -scattxshg

          ! record value in matrix
          call loss_matrix % add_value(scatt_mat_idx, val)

        end if

      end do SCATTR

    end do ROWS 

  end subroutine build_loss_matrix

!===============================================================================
! INDICES_TO_MATRIX takes (x,y,z,g) indices and computes location in matrix
!===============================================================================

  subroutine indices_to_matrix(g, i, j, k, matidx, ng, nx, ny)

    integer :: matidx ! the index location in matrix
    integer :: i      ! current x index
    integer :: j      ! current y index
    integer :: k      ! current z index
    integer :: g      ! current group index
    integer :: nx     ! maximum number of x cells
    integer :: ny     ! maximum number of y cells
    integer :: ng     ! maximum number of groups

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
! MATRIX_TO_INDICES
!===============================================================================

  subroutine matrix_to_indices(irow, g, i, j, k, ng, nx, ny, nz)

    integer :: i     ! iteration counter for x
    integer :: j     ! iteration counter for y
    integer :: k     ! iteration counter for z
    integer :: g     ! iteration counter for groups
    integer :: irow  ! iteration counter over row (0 reference)
    integer :: nx    ! maximum number of x cells
    integer :: ny    ! maximum number of y cells
    integer :: nz    ! maximum number of z cells
    integer :: ng    ! maximum number of groups

    ! check for core map
    if (cmfd_coremap) then

      ! get indices from indexmap
      g = mod(irow, ng) + 1
      i = cmfd % indexmap(irow/ng+1,1)
      j = cmfd % indexmap(irow/ng+1,2)
      k = cmfd % indexmap(irow/ng+1,3)

    else

      ! compute indices
      g = mod(irow, ng) + 1
      i = mod(irow, ng*nx)/ng + 1
      j = mod(irow, ng*nx*ny)/(ng*nx)+ 1
      k = mod(irow, ng*nx*ny*nz)/(ng*nx*ny) + 1

    end if

  end subroutine matrix_to_indices

end module cmfd_loss_operator
