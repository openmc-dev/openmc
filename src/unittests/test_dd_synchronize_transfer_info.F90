module test_dd_synchronize_transfer_info

  use dd_init,          only: initialize_domain_decomp
  use dd_header,        only: dd_type, deallocate_dd
  use dd_testing_setup, only: check_procs, dd_simple_four_domains, &
                              dd_simple_four_domain_scatters
  use global,           only: master, rank
  use output,           only: write_message
  use string,           only: to_str
  use testing_header,   only: TestSuiteClass, TestClass

#ifdef MPI
  use dd_comm,          only: synchronize_transfer_info
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

  type(test), public :: dd_sync_transfer_info_test
  
  type(dd_type) :: dd

contains

!===============================================================================
! INIT
!===============================================================================

  subroutine test_init(this)

    class(test), intent(inout) :: this

    this % name = "test_dd_synchronize_transfer_info"
    
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

    ! Set up local outscatter
    call dd_simple_four_domain_scatters(dd)

  end subroutine test_setup

!===============================================================================
! EXECUTE
!===============================================================================

  subroutine test_execute()

#ifdef MPI
    call synchronize_transfer_info(dd)
#endif

  end subroutine test_execute

!===============================================================================
! CHECK
!===============================================================================

  subroutine test_check(suite)

    class(TestSuiteClass), intent(inout) :: suite

    logical :: failure = .false.
    integer :: bin
#ifdef MPI
    integer :: mpi_err
    logical :: any_fail
#endif

    select case(rank)
      case(0, 1)
        bin = dd % bins_dict % get_key(2)
        if (any(dd % n_scatters_neighborhood(:, bin) /= 0)) failure = .true.
        bin = dd % bins_dict % get_key(3)
        if (dd % n_scatters_neighborhood(16, bin) /= 9) failure = .true.
      case(2)
        bin = dd % bins_dict % get_key(1)
        if (dd % n_scatters_neighborhood(20, bin) /= 7) failure = .true.
        bin = dd % bins_dict % get_key(4)
        if (dd % n_scatters_neighborhood(15, bin) /= 1) failure = .true.
      case(3)
        bin = dd % bins_dict % get_key(1)
        if (dd % n_scatters_neighborhood(10, bin) /= 12) failure = .true.
        bin = dd % bins_dict % get_key(4)
        if (dd % n_scatters_neighborhood(25, bin) /= 3) failure = .true.
      case(4)
        bin = dd % bins_dict % get_key(2)
        if (dd % n_scatters_neighborhood(9, bin) /= 9) failure = .true.
        bin = dd % bins_dict % get_key(3)
        if (dd % n_scatters_neighborhood(19, bin) /= 14) failure = .true.
    end select    

    if (failure) then
      call write_message("FAILURE: Rank " // trim(to_str(rank)) // &
          " didn't receive all the info it should have about its neighbors.")
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

end module test_dd_synchronize_transfer_info
