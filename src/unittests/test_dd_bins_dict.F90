module test_dd_bins_dict

  use dd_init,          only: initialize_domain_decomp
  use dd_header,        only: dd_type, deallocate_dd
  use dd_testing_setup, only: check_procs, dd_simple_four_domains
  use global,           only: master, rank
  use output,           only: write_message
  use string,           only: to_str
  use testing_header,   only: TestSuiteClass, TestClass

#ifdef MPI
  use mpi
#endif

  implicit none
  private

  type(dd_type), save :: dd

  type, extends(TestClass) :: test
    contains
      procedure         :: init     => test_init
      procedure         :: setup    => test_setup
      procedure, nopass :: execute  => test_execute
      procedure, nopass :: check    => test_check
      procedure, nopass :: teardown => test_teardown
  end type test

  type(test), public :: dd_bins_dict_test

contains

!===============================================================================
! INIT
!===============================================================================

  subroutine test_init(this)

    class(test), intent(inout) :: this

    this % name = "test_dd_bins_dict"

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
    
  end subroutine test_setup

!===============================================================================
! EXECUTE
!===============================================================================

  subroutine test_execute()

    call initialize_domain_decomp(dd)
    
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

    select case(rank)
      case(0, 1)
        if (.not. dd % bins_dict % get_key(2) == 4) failure = .true.
        if (.not. dd % bins_dict % get_key(3) == 2) failure = .true.
      case(2)
        if (.not. dd % bins_dict % get_key(1) == 3) failure = .true.
        if (.not. dd % bins_dict % get_key(4) == 2) failure = .true.
      case(3)
        if (.not. dd % bins_dict % get_key(1) == 1) failure = .true.
        if (.not. dd % bins_dict % get_key(4) == 4) failure = .true.
      case(4)
        if (.not. dd % bins_dict % get_key(2) == 1) failure = .true.
        if (.not. dd % bins_dict % get_key(3) == 3) failure = .true.
    end select

    if (failure) then
      call write_message("FAILURE: Domain bins_dict is incorrect on rank " // &
                trim(to_str(rank)))
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

end module test_dd_bins_dict
