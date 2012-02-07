module cmfd_execute

  use cmfd_data,         only: set_up_cmfd
  use cmfd_output,       only: write_cmfd_vtk
  use global,            only: cmfd,cmfd_only,time_cmfd,master,rank,mpi_err,   &
 &                             current_cycle,n_inactive,n_cycles
  use timing,            only: timer_start,timer_stop


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

#ifdef PETSC
      ! execute snes solver
      call cmfd_snes_execute()
#endif

      ! stop timer
      call timer_stop(time_cmfd)

      ! write vtk file
      !if(.not. cmfd_only) call write_cmfd_vtk()

    end if

#ifdef PETSC
! finalize slepc
    if (current_cycle == n_cycles) call SlepcFinalize(ierr)

    ! sync up procs
    call MPI_Barrier(MPI_COMM_WORLD,mpi_err)
#endif

  end subroutine execute_cmfd

end module cmfd_execute
