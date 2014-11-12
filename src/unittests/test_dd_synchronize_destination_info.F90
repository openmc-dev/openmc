module test_dd_synchronize_destination_info

  use constants
  use dd_init,          only: initialize_domain_decomp
  use dd_header,        only: dd_type, deallocate_dd
  use dd_testing_setup, only: check_procs, dd_simple_four_domains, &
                              dd_simple_four_domain_scatters
  use global,           only: master, rank
  use output,           only: write_message
  use string,           only: to_str
  use testing_header,   only: TestSuiteClass, TestClass

#ifdef MPI
  use dd_comm,          only: synchronize_transfer_info, &
                              synchronize_destination_info
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

  type(test), public :: dd_sync_destination_info_test
  
  type(dd_type) :: dd

contains

!===============================================================================
! INIT
!===============================================================================

  subroutine test_init(this)

    class(test), intent(inout) :: this

    this % name = "test_dd_synchronize_destination_info"
    
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

    ! Setup particle buffer
    dd % particle_buffer % outscatter_destination = NO_OUTSCATTER

    ! Set up local outscatter
    call dd_simple_four_domain_scatters(dd)

    ! Get dd % n_scatters_neighborhood
#ifdef MPI
    call synchronize_transfer_info(dd)
#endif

  end subroutine test_setup

!===============================================================================
! EXECUTE
!===============================================================================

  subroutine test_execute()

#ifdef MPI
    call synchronize_destination_info(dd)
#endif

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
    integer :: to_bin, pr_bin

    ! CHECK
    select case(rank)
      case(0)
        to_bin = dd % bins_dict % get_key(2)
        pr_bin = (to_bin - 1) * dd % max_domain_procs + 1
        if (.not. dd % send_rank_info(pr_bin) == 7) failure = .true.
        to_bin = dd % bins_dict % get_key(3)
        pr_bin = (to_bin - 1) * dd % max_domain_procs + 1
        if (.not. dd % send_rank_info(pr_bin) == 5) failure = .true.
      case(1)
        to_bin = dd % bins_dict % get_key(2)
        pr_bin = (to_bin - 1) * dd % max_domain_procs + 1
        if (.not. dd % send_rank_info(pr_bin) == 2) failure = .true.
        to_bin = dd % bins_dict % get_key(3)
        pr_bin = (to_bin - 1) * dd % max_domain_procs + 1
        if (.not. dd % send_rank_info(pr_bin) == 9) failure = .true.
      case(2)
        to_bin = dd % bins_dict % get_key(1)
        pr_bin = (to_bin - 1) * dd % max_domain_procs + 1
        if (.not. dd % send_rank_info(pr_bin) == 3) failure = .true.
        pr_bin = pr_bin + 1
        if (.not. dd % send_rank_info(pr_bin) == 9) failure = .true.
        to_bin = dd % bins_dict % get_key(4)
        pr_bin = (to_bin - 1) * dd % max_domain_procs + 1
        if (.not. dd % send_rank_info(pr_bin) == 3) failure = .true.
      case(3)
        to_bin = dd % bins_dict % get_key(1)
        pr_bin = (to_bin - 1) * dd % max_domain_procs + 1
        if (.not. dd % send_rank_info(pr_bin) == 7) failure = .true.
        pr_bin = pr_bin + 1
        if (.not. dd % send_rank_info(pr_bin) == 0) failure = .true.
        to_bin = dd % bins_dict % get_key(4)
        pr_bin = (to_bin - 1) * dd % max_domain_procs + 1
        if (.not. dd % send_rank_info(pr_bin) == 1) failure = .true.
      case(4)
        to_bin = dd % bins_dict % get_key(2)
        pr_bin = (to_bin - 1) * dd % max_domain_procs + 1
        if (.not. dd % send_rank_info(pr_bin) == 0) failure = .true.
        to_bin = dd % bins_dict % get_key(3)
        pr_bin = (to_bin - 1) * dd % max_domain_procs + 1
        if (.not. dd % send_rank_info(pr_bin) == 9) failure = .true.
    end select

    if (failure) then
      call write_message("FAILURE: Rank " // trim(to_str(rank)) // &
          " calculated the wrong number of particles to send to " // &
          "a process on a neighboring domain.")
    end if

#ifdef MPI
    call MPI_ALLREDUCE(failure, any_fail, 1, MPI_LOGICAL, MPI_LOR, &
        MPI_COMM_WORLD, mpi_err)
    if (.not. any_fail) then
      call suite % pass()
    else
      call suite % fail()
    end if
#else
    call suite % fail()
#endif
    
  end subroutine test_check

!===============================================================================
! TEARDOWN
!===============================================================================

  subroutine test_teardown()

    call deallocate_dd(dd)
    
  end subroutine test_teardown

end module test_dd_synchronize_destination_info
