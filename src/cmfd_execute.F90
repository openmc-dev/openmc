module cmfd_execute

  use cmfd_data,         only: set_up_cmfd
  use cmfd_output,       only: write_cmfd_vtk
  use cmfd_power_solver, only: cmfd_power_execute
  use cmfd_slepc_solver, only: cmfd_slepc_execute
  use cmfd_snes_solver,  only: cmfd_snes_execute
  use global,            only: cmfd,cmfd_only,time_cmfd,master,rank,mpi_err
  use timing,            only: timer_start,timer_stop

  implicit none


#include <finclude/petsc.h90>
#include <finclude/slepcsys.h>
#include <finclude/slepceps.h>

contains

!==============================================================================
! EXECUTE_CMFD
!==============================================================================

  subroutine execute_cmfd()

    integer :: ierr  ! petsc error code

    ! initialize slepc/petsc (communicates to world)
    call SlepcInitialize(PETSC_NULL_CHARACTER,ierr)

    ! only run if master process
    if (master) then

      ! begin timer
      call timer_start(time_cmfd)

      ! set up cmfd
      if(.not. cmfd_only) call set_up_cmfd()

      ! execute snes solver
      call cmfd_snes_execute()

      ! stop timer
      call timer_stop(time_cmfd)

      ! print results
      print *,'SNES Eigenvalue is:',cmfd%keff
      print *,'CMFD Time:',time_cmfd%elapsed,' seconds.'

      ! write vtk file
      !if(.not. cmfd_only) call write_cmfd_vtk()

    end if

    ! finalize slepc
    call SlepcFinalize(ierr)     

    ! sync up procs
    call MPI_Barrier(MPI_COMM_WORLD,mpi_err)

  end subroutine execute_cmfd

end module cmfd_execute
