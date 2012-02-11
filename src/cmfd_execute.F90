module cmfd_execute

  use cmfd_data,         only: set_up_cmfd
  use cmfd_output,       only: write_cmfd_vtk
  use global,            only: cmfd,cmfd_only,time_cmfd,master,rank,mpi_err,   &
 &                             current_cycle,n_inactive,n_cycles
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
#endif

contains

!==============================================================================
! EXECUTE_CMFD
!==============================================================================

  subroutine execute_cmfd()

    integer :: ierr  ! petsc error code

#ifdef PETSC
    ! initialize slepc/petsc (communicates to world)
    if(current_cycle == n_inactive + 1) call SlepcInitialize                   &
   &                                       (PETSC_NULL_CHARACTER,ierr)
#endif

    ! only run if master process
    if (master) then

      ! begin timer
      call timer_start(time_cmfd)

      ! set up cmfd
      if(.not. cmfd_only) call set_up_cmfd()

    end if

    ! broadcast cmfd object to all procs
    call cmfd_bcast()

#ifdef PETSC
      ! execute snes solver
      call cmfd_slepc_execute()
#endif

      ! stop timer
      call timer_stop(time_cmfd)

      ! write vtk file
      !if(.not. cmfd_only) call write_cmfd_vtk()

#ifdef PETSC
! finalize slepc
    if (current_cycle == n_cycles) call SlepcFinalize(ierr)

    ! sync up procs
    call MPI_Barrier(MPI_COMM_WORLD,mpi_err)
#endif

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
    integer rank

    ! extract spatial and energy indices from object
    nx = cmfd % indices(1)
    ny = cmfd % indices(2)
    nz = cmfd % indices(3)
    ng = cmfd % indices(4)

    ! initialize data
    call allocate_cmfd(cmfd)

    ! sync up procs
    call MPI_Barrier(MPI_COMM_WORLD,mpi_err)

    ! broadcast all data
    call MPI_BCAST(cmfd%flux,ng*nx*ny*nz,MPI_REAL8,0,MPI_COMM_WORLD,mpi_err)
    call MPI_BCAST(cmfd%totalxs,ng*nx*ny*nz,MPI_REAL8,0,MPI_COMM_WORLD,mpi_err)
    call MPI_BCAST(cmfd%p1scattxs,ng*nx*ny*nz,MPI_REAL8,0,MPI_COMM_WORLD,mpi_err)
    call MPI_BCAST(cmfd%scattxs,ng*ng*nx*ny*nz,MPI_REAL8,0,MPI_COMM_WORLD,mpi_err)
    call MPI_BCAST(cmfd%nfissxs,ng*ng*nx*ny*nz,MPI_REAL8,0,MPI_COMM_WORLD,mpi_err)
    call MPI_BCAST(cmfd%diffcof,ng*nx*ny*nz,MPI_REAL8,0,MPI_COMM_WORLD,mpi_err)
    call MPI_BCAST(cmfd%dtilde,6*ng*nx*ny*nz,MPI_REAL8,0,MPI_COMM_WORLD,mpi_err)
    call MPI_BCAST(cmfd%dhat,6*ng*nx*ny*nz,MPI_REAL8,0,MPI_COMM_WORLD,mpi_err)
    call MPI_BCAST(cmfd%hxyz,3*nx*ny*nz,MPI_REAL8,0,MPI_COMM_WORLD,mpi_err)
    call MPI_BCAST(cmfd%current,12*ng*nx*ny*nz,MPI_REAL8,0,MPI_COMM_WORLD,mpi_err)

    ! broadcast coremap info
    if (cmfd_coremap) then
      call MPI_BCAST(cmfd%coremap,nx*ny*nz,MPI_INT,0,MPI_COMM_WORLD,mpi_err)
      call MPI_BCAST(cmfd%mat_dim,1,MPI_INT,0,MPI_COMM_WORLD,mpi_err)
      if (.not. allocated(cmfd % indexmap)) allocate                           &
     &           (cmfd % indexmap(cmfd % mat_dim,3))
      call MPI_BCAST(cmfd%indexmap,cmfd%mat_dim*3,MPI_INT,0,MPI_COMM_WORLD,mpi_err)
    end if

    ! sync up procs
    call MPI_Barrier(MPI_COMM_WORLD,mpi_err)

  end subroutine cmfd_bcast

end module cmfd_execute
