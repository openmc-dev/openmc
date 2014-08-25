module test_dd_comm

  use dd_comm,          only: synchronize_transfer_info
  use dd_header,        only: dd_type
  use error,            only: fatal_error, warning
  use global,           only: n_procs, rank, message
  use output,           only: write_message

  implicit none
  public

contains

!===============================================================================
! TEST_SYNCHRONIZE_INFO tests the functionality of the synchronize_info
! subroutine in dd_comm.F90
!===============================================================================

  subroutine test_synchronize_transfer_info()

    type(dd_type) :: dd

#ifdef MPI
    if (.not. n_procs == 4) then
      message = "Skipping test_synchronize_info: must be run with MPI 4 procs"
      call warning()
      return
    end if
    message = "Running test_synchronize_info..."
    call write_message(1)
#else
    message = "Skipping test_synchronize_info: requires MPI"
    call warning()
    return
#endif

    select case(rank)
      case (0)
      case (1)
      case (2)
      case (3)
    end select

  end subroutine test_synchronize_transfer_info

end module test_dd_comm
