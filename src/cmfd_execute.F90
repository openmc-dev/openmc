module cmfd_execute

  use cmfd_data
  use cmfd_power_solver, only: cmfd_power_execute
  use cmfd_slepc_solver, only: cmfd_slepc_execute
  use cmfd_snes_solver,  only: cmfd_snes_execute
  use global

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

    ! set up cmfd
    call set_up_cmfd()    

    ! initialize slepc/petsc
    call SlepcInitialize(PETSC_NULL_CHARACTER,ierr)

    ! execute diffusion solver
    call cmfd_power_execute()

    ! print results
    print *,'Power Eigenvalue is:',cmfd%keff

    ! execute slepc solver
    call cmfd_slepc_execute()

    ! print results
    print *,'SLEPC Eigenvalue is:',cmfd%keff

    ! execute snes solver
    call cmfd_snes_execute()

    ! print results
    print *,'SNES Eigenvalue is:',cmfd%keff

    ! finalize slepc
    call SlepcFinalize(ierr)

  end subroutine execute_cmfd

end module cmfd_execute
