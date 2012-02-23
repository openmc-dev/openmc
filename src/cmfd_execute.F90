module cmfd_execute

  use cmfd_data,         only: set_up_cmfd
  use cmfd_output,       only: write_cmfd_vtk
  use global,            only: cmfd,cmfd_only,time_cmfd,master,rank,mpi_err,   &
 &                             current_cycle,n_inactive,n_cycles,n_procs,n_procs_cmfd
  use timing,            only: timer_start,timer_stop,timer_reset


#ifdef PETSC
  use cmfd_power_solver, only: cmfd_power_execute
  use cmfd_slepc_solver, only: cmfd_slepc_execute
  use cmfd_snes_solver,  only: cmfd_snes_execute
#endif
 
  implicit none

#ifdef PETSC
# include <finclude/petsc.h90>
# include <finclude/slepcsys.h>
# include <finclude/slepceps.h>

contains

!==============================================================================
! EXECUTE_CMFD
!==============================================================================

  subroutine execute_cmfd()

    integer :: ierr  ! petsc error code

    ! initialize mpi communicator
    call petsc_init_mpi()

    if (rank < 6) then

    ! initialize slepc/petsc (communicates to world)
    if(current_cycle == n_inactive + 1) call SlepcInitialize                   &
   &                                       (PETSC_NULL_CHARACTER,ierr)

    ! set global variable for number of procs in cmfd calc
    call MPI_COMM_SIZE(PETSC_COMM_WORLD,n_procs_cmfd,ierr)

    ! only run if master process
    if (master) then

      ! begin timer
      call timer_start(time_cmfd)

      ! set up cmfd
      if(.not. cmfd_only) call set_up_cmfd()

    end if

    ! broadcast cmfd object to all procs
    call cmfd_bcast()

      ! execute snes solver
      call cmfd_snes_execute()

      ! stop timer
      call timer_stop(time_cmfd)

      ! write vtk file
      !if(.not. cmfd_only) call write_cmfd_vtk()

    ! finalize slepc
    if (current_cycle == n_cycles) call SlepcFinalize(ierr)

    end if

    ! sync up procs
    call MPI_Barrier(MPI_COMM_WORLD,mpi_err)

  end subroutine execute_cmfd

!==============================================================================
! EXECUTE_CMFD
!==============================================================================

  subroutine cmfd_bcast()

    use cmfd_header, only: allocate_cmfd
    use global,      only: cmfd_coremap

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
    call allocate_cmfd(cmfd)

    ! sync up procs
    call MPI_Barrier(PETSC_COMM_WORLD,mpi_err)

    ! broadcast all data
    call MPI_BCAST(cmfd%flux,ng*nx*ny*nz,MPI_REAL8,0,PETSC_COMM_WORLD,mpi_err)
    call MPI_BCAST(cmfd%totalxs,ng*nx*ny*nz,MPI_REAL8,0,PETSC_COMM_WORLD,mpi_err)
    call MPI_BCAST(cmfd%p1scattxs,ng*nx*ny*nz,MPI_REAL8,0,PETSC_COMM_WORLD,mpi_err)
    call MPI_BCAST(cmfd%scattxs,ng*ng*nx*ny*nz,MPI_REAL8,0,PETSC_COMM_WORLD,mpi_err)
    call MPI_BCAST(cmfd%nfissxs,ng*ng*nx*ny*nz,MPI_REAL8,0,PETSC_COMM_WORLD,mpi_err)
    call MPI_BCAST(cmfd%diffcof,ng*nx*ny*nz,MPI_REAL8,0,PETSC_COMM_WORLD,mpi_err)
    call MPI_BCAST(cmfd%dtilde,6*ng*nx*ny*nz,MPI_REAL8,0,PETSC_COMM_WORLD,mpi_err)
    call MPI_BCAST(cmfd%dhat,6*ng*nx*ny*nz,MPI_REAL8,0,PETSC_COMM_WORLD,mpi_err)
    call MPI_BCAST(cmfd%hxyz,3*nx*ny*nz,MPI_REAL8,0,PETSC_COMM_WORLD,mpi_err)
    call MPI_BCAST(cmfd%current,12*ng*nx*ny*nz,MPI_REAL8,0,PETSC_COMM_WORLD,mpi_err)

    ! broadcast coremap info
    if (cmfd_coremap) then
      call MPI_BCAST(cmfd%coremap,nx*ny*nz,MPI_INT,0,PETSC_COMM_WORLD,mpi_err)
      call MPI_BCAST(cmfd%mat_dim,1,MPI_INT,0,PETSC_COMM_WORLD,mpi_err)
      if (.not. allocated(cmfd % indexmap)) allocate                           &
     &           (cmfd % indexmap(cmfd % mat_dim,3))
      call MPI_BCAST(cmfd%indexmap,cmfd%mat_dim*3,MPI_INT,0,PETSC_COMM_WORLD,mpi_err)
    end if

  end subroutine cmfd_bcast

!===============================================================================
! PETSC_INIT_MPI
!===============================================================================

  subroutine petsc_init_mpi()

    integer               :: new_comm   ! new communicator
    integer               :: orig_group ! original MPI group for MPI_COMM_WORLD
    integer               :: new_group  ! new MPI group subset of orig_group
    integer,allocatable   :: ranks(:)   ! ranks to include for petsc
    integer               :: k          ! iteration counter

    ! set ranks 0-6 or min
    if (n_procs >= 6) then
      if (.not. allocated(ranks)) allocate(ranks(0:5))
      ranks = (/0,1,2,3,4,5/)
    else
      if (.not. allocated(ranks)) allocate(ranks(0:n_procs-1))
      ranks = (/(k,k=0,n_procs-1)/) 
    end if

    ! get the original mpi group
    call MPI_COMM_GROUP(MPI_COMM_WORLD,orig_group,mpi_err)

    ! new group init
    call MPI_GROUP_INCL(orig_group,size(ranks),ranks,new_group,mpi_err)

    ! create new communicator
    call MPI_COMM_CREATE(MPI_COMM_WORLD,new_group,new_comm,mpi_err)

    ! deallocate ranks
    if (allocated(ranks)) deallocate(ranks)

    ! set PETSC_COMM_WORLD to this subset
    PETSC_COMM_WORLD = new_comm

  end subroutine petsc_init_mpi

#endif

end module cmfd_execute
