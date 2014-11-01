module test_dd_otf_tally_allocation

  use dd_init,          only: initialize_domain_decomp
  use dd_header,        only: dd_type, deallocate_dd
  use dd_testing_setup, only: check_procs, dd_simple_four_domains, &
                              dd_simple_four_domain_tallies
  use global,           only: master, rank, message, free_memory, tallies
  use error,            only: warning
  use output,           only: write_message
  use particle_header,  only: Particle
  use string,           only: to_str
  use tally,            only: score_analog_tally
  use testing_header,   only: TestSuiteClass, TestClass

#ifdef MPI
  use mpi
#endif

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

  type(test), public :: dd_otf_tally_allocation_test

  type(dd_type)  :: dd
  type(Particle) :: p

contains

!===============================================================================
! INIT
!===============================================================================

  subroutine test_init(this)

    class(test), intent(inout) :: this

    this % name = "test_dd_otf_tally_allocation"

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

    ! Setup tallies and particle for scoring to them
    call dd_simple_four_domain_tallies(p)

  end subroutine test_setup

!===============================================================================
! EXECUTE
!===============================================================================

  subroutine test_execute()

    ! Simulate scoring to tally bins
    select case(rank)
      case(0, 1)
      case(2)
      case(3)
      case(4)
        p % coord % cell = 1
        call score_analog_tally(p)
        call score_analog_tally(p)
        p % coord % cell = 2
        call score_analog_tally(p)
    end select    
    
  end subroutine test_execute

!===============================================================================
! CHECK
!===============================================================================

  subroutine test_check(suite)

    class(TestSuiteClass), intent(inout) :: suite

    ! Check that the results arrays are the right size and have the right values
    select case(rank)
      case(0, 1)
        print *, rank, 'results', tallies(1) % results
      case(2)
        print *, rank, 'results', tallies(1) % results
      case(3)
        print *, rank, 'results', tallies(1) % results
      case(4)
        print *, rank, 'results', tallies(1) % results
    end select    
    

    if (.true.) then
      message = "FAILURE: Tally results arrays improperly allocated with " // &
                "OTF on rank " // trim(to_str(rank))
      call write_message()
    end if
    call suite % fail()
    
  end subroutine test_check

!===============================================================================
! TEARDOWN
!===============================================================================

  subroutine test_teardown()

    call p % clear()
    call deallocate_dd(dd)
    call free_memory()
    
  end subroutine test_teardown

end module test_dd_otf_tally_allocation
