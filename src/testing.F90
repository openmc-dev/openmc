module testing

  use global,             only: master, message
  use error,              only: fatal_error
  use output,             only: header

#ifdef TESTING
  use test_dd_comm,       only: test_distribute_source, &
                                test_synchronize_transfer_info
  use test_dd_init,       only: test_set_neighbor_meshbins, &
                                test_bins_dict
#endif

  implicit none
  private
  public :: run_tests

contains

!===============================================================================
! RUN_TESTS runs all unit tests
!===============================================================================

  subroutine run_tests()
  
#ifdef TESTING
    if (master) call header("UNIT TESTING MODE", level=1)
    
    ! Domain decomposition tests
    call test_set_neighbor_meshbins()
    call test_bins_dict()
    call test_distribute_source()
    call test_synchronize_transfer_info()
    
#else
    message = "Must be compiled with testing mode enabled to run tests."
    call fatal_error()
#endif
  
  end subroutine run_tests

end module testing
