module testing

  use global,             only: master, message, unittests
  use error,              only: fatal_error
  use output_header,      only: output_message
  use output,             only: header
  use testing_header,     only: TestSuiteClass, TestClass

#ifdef TESTING
  use test_dd_neighbor_meshbins,            only: dd_neighbor_meshbins_test
  use test_dd_bins_dict,                    only: dd_bins_dict_test
  use test_dd_distribute_source,            only: dd_distribute_source_test
  use test_dd_synchronize_transfer_info,    only: dd_sync_transfer_info_test
  use test_dd_synchronize_destination_info, only: dd_sync_destination_info_test
  use test_dd_otf_tally_allocation,         only: dd_otf_tally_allocation_test
  use test_dd_reduce_tally_results,         only: dd_reduce_tally_results_test
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
    
    unittests % master = master
    
    ! Domain decomposition tests
    call run_test(unittests, dd_neighbor_meshbins_test)
    call run_test(unittests, dd_bins_dict_test)
    call run_test(unittests, dd_distribute_source_test)
    call run_test(unittests, dd_sync_transfer_info_test)
    call run_test(unittests, dd_sync_destination_info_test)
    ! TODO: add test for DD send_recv_particles
    ! TODO: add test for synchronize_bank_dd
    call run_test(unittests, dd_otf_tally_allocation_test)
!    call run_test(unittests, dd_reduce_tally_results_test)
#else
    message = "Must be compiled with testing mode enabled to run tests."
    call fatal_error()
#endif
  
  end subroutine run_tests

!===============================================================================
! RUN_TEST
!===============================================================================

  subroutine run_test(suite, test)
    
    class(TestSuiteClass), intent(inout) :: suite
    class(TestClass),      intent(inout) :: test

    call test % init()

    if (master) call header(test % name, level=2)

    if (master) call output_message("Setting up...")
    
    call test % setup(suite)
    
    if (test % skip) return
    
    if (master) call output_message("Invoking test...")
    
    call test % execute()
    
    if (master) call output_message("Checking results...")
    
    call test % check(suite)
    
    if (master) call output_message("Tearing down...")
    
    call test % teardown()

  end subroutine run_test

end module testing
