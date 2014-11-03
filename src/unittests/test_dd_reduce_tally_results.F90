module test_dd_reduce_tally_results

  use dd_init,          only: initialize_domain_decomp
  use dd_header,        only: dd_type, deallocate_dd
  use dd_testing_setup, only: check_procs, dd_simple_four_domains
  use error,            only: warning
  use global,           only: master
  use testing_header,   only: TestSuiteClass, TestClass

  implicit none
  private

  type, extends(TestClass) :: test
    contains
      procedure         :: init     => test_init
      procedure         :: setup    => test_setup
      procedure, nopass :: execute  => test_execute
      procedure, nopass :: check    => test_check
      procedure, nopass :: teardown => test_teardown
  end type test

  type(test), public :: dd_reduce_tally_results_test
  
  type(dd_type) :: dd

contains

!===============================================================================
! INIT
!===============================================================================

  subroutine test_init(this)

    class(test), intent(inout) :: this

    this % name = "test_dd_reduce_tally_results"

  end subroutine test_init

!===============================================================================
! SETUP
!===============================================================================

  subroutine test_setup(this, suite)

    class(test),      intent(inout) :: this

    class(TestSuiteClass), intent(inout) :: suite
    
    if (check_procs(5)) then
      call suite % skip(this)
      return
    end if

    ! Get generic DD setup with 4 domains for 5 MPI ranks
    call dd_simple_four_domains(dd)

    ! Initialize dd_type
    call initialize_domain_decomp(dd)

    if (master) then
      call warning("TEST NOT IMPLEMENTED")
      call suite % skip(this)
    end if

  end subroutine test_setup

!===============================================================================
! EXECUTE
!===============================================================================

  subroutine test_execute()

    ! Call the subroutines or functions to check here
    
  end subroutine test_execute

!===============================================================================
! CHECK
!===============================================================================

  subroutine test_check(suite)

    class(TestSuiteClass), intent(inout) :: suite

    ! Add test checks here.  For example, check that the results are what you
    ! expect given whatever setup was hardcoded into test_setup
    
    ! If success, do:  call suite % pass()
    ! If failure, do:  call suite % fail()
    
  end subroutine test_check

!===============================================================================
! TEARDOWN
!===============================================================================

  subroutine test_teardown()

    call deallocate_dd(dd)
    
  end subroutine test_teardown

end module test_dd_reduce_tally_results
