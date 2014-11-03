module test_dd_otf_tally_allocation

  use dd_init,          only: initialize_domain_decomp
  use dd_header,        only: dd_type, deallocate_dd
  use dd_testing_setup, only: check_procs, dd_simple_four_domains, &
                              dd_simple_four_domain_tallies
  use global,           only: master, rank, free_memory, tallies, &
                              domain_decomp, dd_run
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
    call dd_simple_four_domains(domain_decomp)

    ! Initialize dd_type
    call initialize_domain_decomp(domain_decomp)
    dd_run = .true.

    ! Setup tallies and particle for scoring to them
    call dd_simple_four_domain_tallies(p)

  end subroutine test_setup

!===============================================================================
! EXECUTE
!===============================================================================

  subroutine test_execute()

    ! Simulate scoring to tally bins
    select case(rank)
      case(0)
        p % coord % cell = 1
        call score_analog_tally(p)
        p % coord % cell = 2
        call score_analog_tally(p)
        call score_analog_tally(p)
        call score_analog_tally(p)
        p % coord % cell = 4
        call score_analog_tally(p)
        call score_analog_tally(p)
        call score_analog_tally(p)
        call score_analog_tally(p)
      case(1)
        p % coord % cell = 4
        call score_analog_tally(p)
        p % coord % cell = 2
        call score_analog_tally(p)
        call score_analog_tally(p)
        call score_analog_tally(p)
        p % coord % cell = 1
        call score_analog_tally(p)
        call score_analog_tally(p)
        call score_analog_tally(p)
        call score_analog_tally(p)
      case(2)
        p % coord % cell = 3
        call score_analog_tally(p)
        p % coord % cell = 4
        call score_analog_tally(p)
        p % coord % cell = 3
        call score_analog_tally(p)
        call score_analog_tally(p)
        p % coord % cell = 4
        call score_analog_tally(p)
        call score_analog_tally(p)
        call score_analog_tally(p)
      case(3)
        p % coord % cell = 2
        call score_analog_tally(p)
        call score_analog_tally(p)
        p % coord % cell = 4
        call score_analog_tally(p)
      case(4)
        p % coord % cell = 2
        call score_analog_tally(p)
        p % coord % cell = 4
        call score_analog_tally(p)
        p % coord % cell = 3
        call score_analog_tally(p)
        call score_analog_tally(p)
        p % coord % cell = 2
        call score_analog_tally(p)
        call score_analog_tally(p)
        p % coord % cell = 4
        call score_analog_tally(p)
        call score_analog_tally(p)
    end select    
    
  end subroutine test_execute

!===============================================================================
! CHECK
!===============================================================================

  subroutine test_check(suite)

    class(TestSuiteClass), intent(inout) :: suite

    logical :: failure = .false.
#ifdef MPI
    integer :: mpi_err
    logical :: any_fail
#endif

    ! Check that the results arrays are the right size and have the right values
    select case(rank)
      case(0, 1)
        if (size(tallies(1) % results, 2) /= 3) then
          failure = .true.
        else
          if (.not. int(tallies(1) % results(1,1) % value) == 1 .and. &
              .not. int(tallies(1) % results(1,2) % value) == 3 .and. &
              .not. int(tallies(1) % results(1,3) % value) == 4) failure = .true.
        end if
      case(2)
        if (size(tallies(1) % results, 2) /= 2) then
          failure = .true.
        else
          if (.not. int(tallies(1) % results(1,1) % value) == 3 .and. &
              .not. int(tallies(1) % results(1,2) % value) == 4) failure = .true.
        end if
      case(3)
        if (size(tallies(1) % results, 2) /= 2) then
          failure = .true.
        else
          if (.not. int(tallies(1) % results(1,1) % value) == 2 .and. &
              .not. int(tallies(1) % results(1,2) % value) == 1) failure = .true.
        end if
      case(4)
        if (size(tallies(1) % results, 2) /= 3) then
          failure = .true.
        else
          if (.not. int(tallies(1) % results(1,1) % value) == 3 .and. &
              .not. int(tallies(1) % results(1,2) % value) == 3 .and. &
              .not. int(tallies(1) % results(1,3) % value) == 2) failure = .true.
        end if
    end select    
    
    if (failure) then
      call write_message("FAILURE: Tally results arrays improperly " // &
          "allocated with OTF on rank " // trim(to_str(rank)))
    end if

#ifdef MPI
    call MPI_ALLREDUCE(failure, any_fail, 1, MPI_LOGICAL, MPI_LOR, &
        MPI_COMM_WORLD, mpi_err)
    if (.not. any_fail) then
      call suite % pass()
    else
      call suite % fail()
    end if
#endif

  end subroutine test_check

!===============================================================================
! TEARDOWN
!===============================================================================

  subroutine test_teardown()

    call p % clear()
    call deallocate_dd(domain_decomp)
    call free_memory()
    
  end subroutine test_teardown

end module test_dd_otf_tally_allocation
