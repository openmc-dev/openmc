module cmfd_loss_operator

#ifdef PETSC

  implicit none
  private
  public :: init_M_operator,build_loss_matrix,destroy_M_operator

#include "finclude/petsc.h90"

    integer  :: nx   ! maximum number of x cells
    integer  :: ny   ! maximum number of y cells
    integer  :: nz   ! maximum number of z cells
    integer  :: ng   ! maximum number of groups
    integer  :: ierr ! petsc error code

  type, public :: loss_operator

    Mat      :: M      ! petsc matrix for neutronic loss operator
    integer  :: n      ! dimensions of matrix
    integer  :: nnz    ! max number of nonzeros
    integer  :: localn ! local size on proc

  end type loss_operator

contains

!===============================================================================
! INIT_M_OPERATOR
!===============================================================================

  subroutine init_M_operator(this)

    type(loss_operator) :: this

    ! get indices
    call get_M_indices(this)

    ! get preallocation
    call preallocate_loss_matrix(this)

    ! set up M operator
    call MatCreateMPIAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,this%n,    &
   & this%n,this%nnz,PETSC_NULL_INTEGER,this%nnz,PETSC_NULL_INTEGER,this%M,ierr)
    call MatSetOption(this%M,MAT_NEW_NONZERO_LOCATIONS,PETSC_TRUE,ierr)
    call MatSetOption(this%M,MAT_IGNORE_ZERO_ENTRIES,PETSC_TRUE,ierr)

  end subroutine init_M_operator

!===============================================================================
! GET_M_INDICES
!===============================================================================

  subroutine get_M_indices(this)

    use global, only: cmfd,cmfd_coremap

    type(loss_operator) :: this

    ! get maximum number of cells in each direction
    nx = cmfd%indices(1)
    ny = cmfd%indices(2)
    nz = cmfd%indices(3)
    ng = cmfd%indices(4)

    ! get number of nonzeros
    this%nnz = 7 + ng - 1

    ! calculate dimensions of matrix
    if (cmfd_coremap) then
      this%n = cmfd % mat_dim * ng
    else
      this%n = nx*ny*nz*ng
    end if

  end subroutine get_M_indices


!===============================================================================
! PREALLOCATE_LOSS_MATRIX
!===============================================================================

  subroutine preallocate_loss_matrix(this)

    type(loss_operator) :: this

    integer :: rank      ! rank of processor
    integer :: size      ! number of procs
    integer :: n         ! the extent of the matrix
    integer :: row_start ! index of local starting row
    integer :: row_end   ! index of local final row

    ! get rank and max rank of procs
    call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,size,ierr)

    ! get local problem size
    n = this%n

    ! determine local size, divide evenly between all other procs
    this%localn = n/(size)
        
    ! add 1 more if less proc id is less than mod
    if (rank < mod(n,size)) this%localn = this%localn + 1


    ! determine local starting row
    row_start = 0
    if (rank < mod(n,size)) then
      row_start = rank*(n/size+1)
    else
      row_start = min(mod(n,size)*(n/size+1)+(rank - mod(n,size))*(n/size),n) 
    end if

    ! determine local final row
    row_end = row_start + this%localn - 1

    print*,'Me:',rank,this%localn,row_start,row_end

  end subroutine preallocate_loss_matrix

!===============================================================================
! BUILD_LOSS_MATRIX creates the matrix representing loss of neutrons
!===============================================================================

  subroutine build_loss_matrix(this)

    use global, only: cmfd,cmfd_coremap,mpi_err

    type(loss_operator) :: this

    integer :: nxyz(3,2)            ! single vector containing bound. locations
    integer :: i                    ! iteration counter for x
    integer :: j                    ! iteration counter for y
    integer :: k                    ! iteration counter for z
    integer :: g                    ! iteration counter for groups
    integer :: l                    ! iteration counter for leakages
    integer :: h                    ! energy group when doing scattering
    integer :: neig_mat_idx         ! matrix index of neighbor cell
    integer :: scatt_mat_idx        ! matrix index for h-->g scattering terms
    integer :: bound(6)             ! vector for comparing when looking for bound
    integer :: xyz_idx              ! index for determining if x,y or z leakage
    integer :: dir_idx              ! index for determining - or + face of cell
    integer :: neig_idx(3)          ! spatial indices of neighbour
    integer :: shift_idx            ! parameter to shift index by +1 or -1
    integer :: kount                ! integer for counting values in vector
    integer :: row_start            ! the first local row on the processor
    integer :: row_finish           ! the last local row on the processor
    integer :: irow                 ! iteration counter over row
    real(8) :: totxs                ! total macro cross section
    real(8) :: scattxsgg            ! scattering macro cross section g-->g
    real(8) :: scattxshg            ! scattering macro cross section h-->g
    real(8) :: dtilde(6)            ! finite difference coupling parameter
    real(8) :: dhat(6)              ! nonlinear coupling parameter
    real(8) :: hxyz(3)              ! cell lengths in each direction
    real(8) :: jn                   ! direction dependent leakage coeff to neig
    real(8) :: jo(6)                ! leakage coeff in front of cell flux
    real(8) :: jnet                 ! net leakage from jo
    real(8) :: val                  ! temporary variable before saving to matrix 
integer :: rank
    ! create single vector of these indices for boundary calculation
    nxyz(1,:) = (/1,nx/)
    nxyz(2,:) = (/1,ny/)
    nxyz(3,:) = (/1,nz/)

    ! get row bounds for this processor
    call MatGetOwnershipRange(this%M,row_start,row_finish,ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
print *,rank,row_start,row_finish-1
    ! begin iteration loops
    ROWS: do irow = row_start,row_finish-1 

      ! get indices for that row
      call matrix_to_indices(irow,g,i,j,k)

      ! retrieve cell data
      totxs = cmfd%totalxs(g,i,j,k)
      scattxsgg = cmfd%scattxs(g,g,i,j,k)
      dtilde = cmfd%dtilde(:,g,i,j,k)
      hxyz = cmfd%hxyz(:,i,j,k)

      ! check and get dhat
      if (allocated(cmfd%dhat)) then
        dhat = cmfd%dhat(:,g,i,j,k)
      else
        dhat = 0.0_8
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
            if (cmfd % coremap(neig_idx(1),neig_idx(2),neig_idx(3)) /=         &
           &    99999) then

              ! compute leakage coefficient for neighbor
              jn = -dtilde(l) + shift_idx*dhat(l)

              ! get neighbor matrix index
              call indices_to_matrix(g,neig_idx(1),neig_idx(2),neig_idx(3),    &
             &                       neig_mat_idx)

              ! compute value and record to bank
              val = jn/hxyz(xyz_idx)

              ! record value in matrix
              call MatSetValue(this%M,irow,neig_mat_idx-1,val,                 &
             &                 INSERT_VALUES,ierr)

            end if

          else

            ! compute leakage coefficient for neighbor
            jn = -dtilde(l) + shift_idx*dhat(l)

            ! get neighbor matrix index
            call indices_to_matrix(g,neig_idx(1),neig_idx(2),neig_idx(3),      &
           &                       neig_mat_idx)

            ! compute value and record to bank
            val = jn/hxyz(xyz_idx)

            ! record value in matrix
            call MatSetValue(this%M,irow,neig_mat_idx-1,val,                   &
            &                 INSERT_VALUES,ierr)

          end if

        end if

        ! compute leakage coefficient for target
        jo(l) = shift_idx*dtilde(l) + dhat(l)

      end do LEAK

      ! calate net leakage coefficient for target
      jnet = (jo(2) - jo(1))/hxyz(1) + (jo(4) - jo(3))/hxyz(2) +               &
     &       (jo(6) - jo(5))/hxyz(3)

      ! calculate loss of neutrons
      val = jnet + totxs - scattxsgg

      ! record diagonal term
      call MatSetValue(this%M,irow,irow,val,INSERT_VALUES,ierr)

      ! begin loop over off diagonal in-scattering
      SCATTR: do h = 1,ng

        ! cycle though if h=g
        if (h == g) then
          cycle
        end if

        ! get neighbor matrix index
        call indices_to_matrix(h,i,j,k,scatt_mat_idx)

        ! get scattering macro xs
        scattxshg = cmfd%scattxs(h,g,i,j,k)

        ! record value in matrix (negate it)
        val = -scattxshg

        call MatSetValue(this%M,irow,scatt_mat_idx-1,val, INSERT_VALUES,ierr)

      end do SCATTR

    end do ROWS 

    ! assemble matrix
    call MatAssemblyBegin(this%M,MAT_FLUSH_ASSEMBLY,ierr)
    call MatAssemblyEnd(this%M,MAT_FINAL_ASSEMBLY,ierr)

    ! print out operator to file
    call print_M_operator(this)
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    stop
  end subroutine build_loss_matrix

!===============================================================================
! INDICES_TO_MATRIX takes (x,y,z,g) indices and computes location in matrix 
!===============================================================================

  subroutine indices_to_matrix(g,i,j,k,matidx)

    use global, only: cmfd,cmfd_coremap

    integer :: matidx         ! the index location in matrix
    integer :: i               ! current x index
    integer :: j               ! current y index
    integer :: k               ! current z index
    integer :: g               ! current group index

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

  subroutine matrix_to_indices(irow,g,i,j,k)

    use global, only: cmfd,cmfd_coremap

    integer :: i                    ! iteration counter for x
    integer :: j                    ! iteration counter for y
    integer :: k                    ! iteration counter for z
    integer :: g                    ! iteration counter for groups
    integer :: irow                 ! iteration counter over row (0 reference)

    ! check for core map
    if (cmfd_coremap) then

      ! get indices from indexmap
      g = mod(irow,ng) + 1
      i = cmfd % indexmap(irow+1,1)
      j = cmfd % indexmap(irow+1,2)
      k = cmfd % indexmap(irow+1,3)

    else

      ! compute indices
      g = mod(irow,ng) + 1 
      i = mod(irow,ng*nx)/ng + 1
      j = mod(irow,ng*nx*ny)/(ng*nx)+ 1
      k = mod(irow,ng*nx*ny*nz)/(ng*nx*ny) + 1 

    end if

  end subroutine matrix_to_indices

!===============================================================================
! PRINT_M_OPERATOR 
!===============================================================================

  subroutine print_M_operator(this)

    use global, only: path_input

    type(loss_operator) :: this

    PetscViewer :: viewer

    ! write out matrix in binary file (debugging)
    call PetscViewerBinaryOpen(PETSC_COMM_WORLD,trim(path_input)//'lossmat.bin' &
   &     ,FILE_MODE_WRITE,viewer,ierr)
    call MatView(this%M,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)

  end subroutine print_M_operator

!==============================================================================
! DESTROY_M_OPERATOR
!==============================================================================

  subroutine destroy_M_operator(this)

    type(loss_operator) :: this

    ! deallocate matrix
    call MatDestroy(this%M,ierr)

  end subroutine destroy_M_operator

#endif

end module cmfd_loss_operator
