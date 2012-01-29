module cmfd_prod_operator

  implicit none
  private
  public :: init_F_operator,build_prod_matrix,destroy_F_operator

#include <finclude/petsc.h90>

    integer  :: nx   ! maximum number of x cells
    integer  :: ny   ! maximum number of y cells
    integer  :: nz   ! maximum number of z cells
    integer  :: ng   ! maximum number of groups
    integer  :: ierr ! petsc error code

  type, public :: prod_operator

    Mat      :: F    ! petsc matrix for neutronic prod operator
    integer  :: n    ! dimensions of matrix
    integer  :: nnz  ! max number of nonzeros

  end type prod_operator

contains

!==============================================================================
! INIT_F_OPERATOR
!==============================================================================

  subroutine init_F_operator(this)

    type(prod_operator) :: this

    ! get indices
    call get_F_indices(this)

    ! set up M operator
    call MatCreateSeqAIJ(PETSC_COMM_SELF,this%n,this%n,this%nnz,               &
   &                     PETSC_NULL_INTEGER,this%F,ierr)
    call MatSetOption(this%F,MAT_NEW_NONZERO_LOCATIONS,PETSC_TRUE,ierr)
    call MatSetOption(this%F,MAT_IGNORE_ZERO_ENTRIES,PETSC_TRUE,ierr)
    call MatSetOption(this%F,MAT_USE_HASH_TABLE,PETSC_TRUE,ierr)

  end subroutine init_F_operator

!==============================================================================
! GET_F_INDICES
!==============================================================================

  subroutine get_F_indices(this)

    use global, only: cmfd,cmfd_coremap

    type(prod_operator) :: this

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

  end subroutine get_F_indices

!===============================================================================
! BUILD_PROD_MATRIX creates the matrix representing loss of neutrons
!===============================================================================

  subroutine build_prod_matrix(this)

    use global, only: cmfd,cmfd_coremap

    type(prod_operator) :: this

    integer :: i                  ! iteration counter for x
    integer :: j                  ! iteration counter for y
    integer :: k                  ! iteration counter for z
    integer :: g                  ! iteration counter for groups
    integer :: l                  ! iteration counter for leakages
    integer :: h                  ! energy group when doing scattering
    integer :: gmat_idx           ! index in matrix for energy group g
    integer :: hmat_idx           ! index in matrix for energy group h
    integer :: ierr               ! Petsc error code
    real(8) :: nfissxs            ! nufission cross section h-->g
    real(8) :: val                ! temporary variable for nfissxs

    ! initialize matrix for building
    call MatAssemblyBegin(this%F,MAT_FLUSH_ASSEMBLY,ierr)

    ! begin loop around energy groups and spatial indices
    ZLOOP:  do k = 1,nz

      YLOOP: do j = 1,ny

        XLOOP: do i = 1,nx

          ! check if not including reflector
          if (cmfd_coremap) then

            ! check if at a reflector
            if (cmfd % coremap(i,j,k) == 99999) then
              cycle
            end if

          end if

          GROUP: do g = 1,ng

            NFISS: do h = 1,ng

              ! get cell data
              nfissxs = cmfd%nfissxs(h,g,i,j,k)

              ! get matrix location
              gmat_idx = get_matrix_idx(g,i,j,k)
              hmat_idx = get_matrix_idx(h,i,j,k)

              ! reocrd value in matrix
              val = nfissxs
              call MatSetValue(this%F,gmat_idx-1,hmat_idx-1,val,INSERT_VALUES, &
             &                 ierr)

            end do NFISS

          end do GROUP

        end do XLOOP

      end do YLOOP

    end do ZLOOP

    ! finalize matrix assembly
    call MatAssemblyEnd(this%F,MAT_FINAL_ASSEMBLY,ierr)

    ! print out operator to file
    call print_F_operator(this)

  end subroutine build_prod_matrix

!===============================================================================
! GET_MATRIX_IDX takes (g,x,y,z) indices and computes location in matrix 
!===============================================================================

  function get_matrix_idx(g,i,j,k)

    use global, only: cmfd,cmfd_coremap

    ! arguments
    integer :: get_matrix_idx  ! the index location in matrix
    integer :: i               ! current x index
    integer :: j               ! current y index
    integer :: k               ! current z index
    integer :: g               ! current group index

    ! local variables
    integer :: nidx            ! index in matrix

    ! check if coremap is used
    if (cmfd_coremap) then

      ! get idx from core map
      nidx = ng*(cmfd % coremap(i,j,k)) - (ng - g)

    else

      ! compute index
      nidx = g + ng*(i - 1) + ng*nx*(j - 1) + ng*nx*ny*(k - 1)

    end if

    ! record value to function
    get_matrix_idx = nidx

  end function get_matrix_idx

!===============================================================================
! PRINT_F_OPERATOR 
!===============================================================================

  subroutine print_F_operator(this)

    PetscViewer :: viewer

    type(prod_operator) :: this

    ! write out matrix in binary file (debugging)
    call PetscViewerBinaryOpen(PETSC_COMM_WORLD,'prodmat.bin',FILE_MODE_WRITE, &
                               viewer,ierr)
    call MatView(this%F,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)

  end subroutine print_F_operator

!==============================================================================
! DESTROY_F_OPERATOR
!==============================================================================

  subroutine destroy_F_operator(this)

    type(prod_operator) :: this

    ! deallocate matrix
    call MatDestroy(this%F,ierr)

  end subroutine destroy_F_operator

end module cmfd_prod_operator
