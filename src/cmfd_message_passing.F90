module cmfd_message_passing

# ifdef PETSC

! This module contains MPI routines related to the CMFD calculation

  implicit none
  private
  public :: petsc_init_mpi, cmfd_bcast 

# include <finclude/petsc.h90>

contains

!===============================================================================
! PETSC_INIT_MPI
!===============================================================================

  subroutine petsc_init_mpi()

    use global,       only: n_procs_cmfd, time_cmfd, rank, mpi_err

    integer               :: new_comm   ! new communicator
    integer               :: color      ! communicator color

    ! assign color
    if (rank < n_procs_cmfd) then
      color = 1
    else
      color = 2
    end if

    ! split up procs
    call MPI_COMM_SPLIT(MPI_COMM_WORLD,color,0,new_comm,mpi_err)

    ! assign to PETSc
    PETSC_COMM_WORLD = new_comm

    ! initialize petsc on all procs
    call PetscInitialize(PETSC_NULL_CHARACTER,mpi_err)

    ! initialize timer
    call time_cmfd % reset()

  end subroutine petsc_init_mpi

!===============================================================================
! CMFD_BCAST
!===============================================================================

  subroutine cmfd_bcast()

    use global,       only: cmfd_coremap, cmfd, mpi_err
    use cmfd_header,  only: allocate_cmfd

    integer :: nx  ! number of mesh cells in x direction
    integer :: ny  ! number of mesh cells in y direction
    integer :: nz  ! number of mesh cells in z direction
    integer :: ng  ! number of energy groups

    ! extract spatial and energy indices from object
    nx = cmfd % indices(1)
    ny = cmfd % indices(2)
    nz = cmfd % indices(3)
    ng = cmfd % indices(4)

    ! initialize data
    if (.not. allocated(cmfd%flux)) call allocate_cmfd(cmfd)

    ! sync up procs
    call MPI_Barrier(PETSC_COMM_WORLD, mpi_err)

    ! broadcast all data
    call MPI_BCAST(cmfd%flux, ng*nx*ny*nz, MPI_REAL8, 0, PETSC_COMM_WORLD, mpi_err)
    call MPI_BCAST(cmfd%totalxs, ng*nx*ny*nz, MPI_REAL8, 0, PETSC_COMM_WORLD, mpi_err)
    call MPI_BCAST(cmfd%p1scattxs, ng*nx*ny*nz, MPI_REAL8, 0, PETSC_COMM_WORLD, mpi_err)
    call MPI_BCAST(cmfd%scattxs, ng*ng*nx*ny*nz, MPI_REAL8, 0, PETSC_COMM_WORLD, mpi_err)
    call MPI_BCAST(cmfd%nfissxs, ng*ng*nx*ny*nz, MPI_REAL8, 0, PETSC_COMM_WORLD, mpi_err)
    call MPI_BCAST(cmfd%diffcof, ng*nx*ny*nz, MPI_REAL8, 0, PETSC_COMM_WORLD, mpi_err)
    call MPI_BCAST(cmfd%dtilde, 6*ng*nx*ny*nz, MPI_REAL8, 0, PETSC_COMM_WORLD, mpi_err)
    call MPI_BCAST(cmfd%dhat, 6*ng*nx*ny*nz, MPI_REAL8, 0, PETSC_COMM_WORLD, mpi_err)
    call MPI_BCAST(cmfd%hxyz, 3*nx*ny*nz, MPI_REAL8, 0, PETSC_COMM_WORLD, mpi_err)
    call MPI_BCAST(cmfd%current, 12*ng*nx*ny*nz, MPI_REAL8, 0, PETSC_COMM_WORLD, mpi_err)

    ! broadcast coremap info
    if (cmfd_coremap) then
      call MPI_BCAST(cmfd%coremap, nx*ny*nz, MPI_INTEGER8, 0, PETSC_COMM_WORLD, mpi_err)
      call MPI_BCAST(cmfd%mat_dim, 1, MPI_INTEGER8, 0, PETSC_COMM_WORLD, mpi_err)
      if (.not. allocated(cmfd % indexmap)) &
           allocate(cmfd % indexmap(cmfd % mat_dim,3))
      call MPI_BCAST(cmfd%indexmap, cmfd%mat_dim*3, MPI_INTEGER8, 0, PETSC_COMM_WORLD, mpi_err)
    end if

  end subroutine cmfd_bcast

# endif

end module cmfd_message_passing
