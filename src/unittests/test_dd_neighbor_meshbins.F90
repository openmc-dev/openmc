module test_dd_neighbor_meshbins

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

  type(test), public :: dd_neighbor_meshbins_test

contains

!===============================================================================
! INIT
!===============================================================================

  subroutine test_init(this)

    class(test), intent(inout) :: this

    this % name = "test_dd_neighbor_meshbins"

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
        if (.not. dd % meshbin == 1) failure = .true.
        if (.not. dd % neighbor_meshbins(1) == -1) failure = .true. ! -x
        if (.not. dd % neighbor_meshbins(2) ==  3) failure = .true. ! +x
        if (.not. dd % neighbor_meshbins(3) == -1) failure = .true. ! -y
        if (.not. dd % neighbor_meshbins(4) ==  2) failure = .true. ! +y
        if (.not. dd % neighbor_meshbins(5) == -1) failure = .true. ! -z
        if (.not. dd % neighbor_meshbins(6) == -1) failure = .true. ! +z
        ! -x
        if (.not. dd % neighbor_meshbins(1 + 6) == -1) failure = .true. ! -x
        if (.not. dd % neighbor_meshbins(2 + 6) ==  1) failure = .true. ! +x
        if (.not. dd % neighbor_meshbins(3 + 6) == -1) failure = .true. ! -y
        if (.not. dd % neighbor_meshbins(4 + 6) == -1) failure = .true. ! +y
        if (.not. dd % neighbor_meshbins(5 + 6) == -1) failure = .true. ! -z
        if (.not. dd % neighbor_meshbins(6 + 6) == -1) failure = .true. ! +z
        ! +x
        if (.not. dd % neighbor_meshbins(1 + 12) ==  1) failure = .true. ! -x
        if (.not. dd % neighbor_meshbins(2 + 12) == -1) failure = .true. ! +x
        if (.not. dd % neighbor_meshbins(3 + 12) == -1) failure = .true. ! -y
        if (.not. dd % neighbor_meshbins(4 + 12) ==  4) failure = .true. ! +y
        if (.not. dd % neighbor_meshbins(5 + 12) == -1) failure = .true. ! -z
        if (.not. dd % neighbor_meshbins(6 + 12) == -1) failure = .true. ! +z
        ! -y
        if (.not. dd % neighbor_meshbins(1 + 18) == -1) failure = .true. ! -x
        if (.not. dd % neighbor_meshbins(2 + 18) == -1) failure = .true. ! +x
        if (.not. dd % neighbor_meshbins(3 + 18) == -1) failure = .true. ! -y
        if (.not. dd % neighbor_meshbins(4 + 18) ==  1) failure = .true. ! +y
        if (.not. dd % neighbor_meshbins(5 + 18) == -1) failure = .true. ! -z
        if (.not. dd % neighbor_meshbins(6 + 18) == -1) failure = .true. ! +z
        ! +y
        if (.not. dd % neighbor_meshbins(1 + 24) == -1) failure = .true. ! -x
        if (.not. dd % neighbor_meshbins(2 + 24) ==  4) failure = .true. ! +x
        if (.not. dd % neighbor_meshbins(3 + 24) ==  1) failure = .true. ! -y
        if (.not. dd % neighbor_meshbins(4 + 24) == -1) failure = .true. ! +y
        if (.not. dd % neighbor_meshbins(5 + 24) == -1) failure = .true. ! -z
        if (.not. dd % neighbor_meshbins(6 + 24) == -1) failure = .true. ! +z
        ! -z
        if (.not. dd % neighbor_meshbins(1 + 30) == -1) failure = .true. ! -x
        if (.not. dd % neighbor_meshbins(2 + 30) == -1) failure = .true. ! +x
        if (.not. dd % neighbor_meshbins(3 + 30) == -1) failure = .true. ! -y
        if (.not. dd % neighbor_meshbins(4 + 30) == -1) failure = .true. ! +y
        if (.not. dd % neighbor_meshbins(5 + 30) == -1) failure = .true. ! -z
        if (.not. dd % neighbor_meshbins(6 + 30) ==  1) failure = .true. ! +z
        ! +z
        if (.not. dd % neighbor_meshbins(1 + 36) == -1) failure = .true. ! -x
        if (.not. dd % neighbor_meshbins(2 + 36) == -1) failure = .true. ! +x
        if (.not. dd % neighbor_meshbins(3 + 36) == -1) failure = .true. ! -y
        if (.not. dd % neighbor_meshbins(4 + 36) == -1) failure = .true. ! +y
        if (.not. dd % neighbor_meshbins(5 + 36) ==  1) failure = .true. ! -z
        if (.not. dd % neighbor_meshbins(6 + 36) == -1) failure = .true. ! +z
      case(2)
        if (.not. dd % meshbin == 2) failure = .true.
        if (.not. dd % neighbor_meshbins(1) == -1) failure = .true. ! -x
        if (.not. dd % neighbor_meshbins(2) ==  4) failure = .true. ! +x
        if (.not. dd % neighbor_meshbins(3) ==  1) failure = .true. ! -y
        if (.not. dd % neighbor_meshbins(4) == -1) failure = .true. ! +y
        if (.not. dd % neighbor_meshbins(5) == -1) failure = .true. ! -z
        if (.not. dd % neighbor_meshbins(6) == -1) failure = .true. ! +z
        ! -x
        if (.not. dd % neighbor_meshbins(1 + 6) == -1) failure = .true. ! -x
        if (.not. dd % neighbor_meshbins(2 + 6) ==  2) failure = .true. ! +x
        if (.not. dd % neighbor_meshbins(3 + 6) == -1) failure = .true. ! -y
        if (.not. dd % neighbor_meshbins(4 + 6) == -1) failure = .true. ! +y
        if (.not. dd % neighbor_meshbins(5 + 6) == -1) failure = .true. ! -z
        if (.not. dd % neighbor_meshbins(6 + 6) == -1) failure = .true. ! +z
        ! +x
        if (.not. dd % neighbor_meshbins(1 + 12) ==  2) failure = .true. ! -x
        if (.not. dd % neighbor_meshbins(2 + 12) == -1) failure = .true. ! +x
        if (.not. dd % neighbor_meshbins(3 + 12) ==  3) failure = .true. ! -y
        if (.not. dd % neighbor_meshbins(4 + 12) == -1) failure = .true. ! +y
        if (.not. dd % neighbor_meshbins(5 + 12) == -1) failure = .true. ! -z
        if (.not. dd % neighbor_meshbins(6 + 12) == -1) failure = .true. ! +z
        ! -y
        if (.not. dd % neighbor_meshbins(1 + 18) == -1) failure = .true. ! -x
        if (.not. dd % neighbor_meshbins(2 + 18) ==  3) failure = .true. ! +x
        if (.not. dd % neighbor_meshbins(3 + 18) == -1) failure = .true. ! -y
        if (.not. dd % neighbor_meshbins(4 + 18) ==  2) failure = .true. ! +y
        if (.not. dd % neighbor_meshbins(5 + 18) == -1) failure = .true. ! -z
        if (.not. dd % neighbor_meshbins(6 + 18) == -1) failure = .true. ! +z
        ! +y
        if (.not. dd % neighbor_meshbins(1 + 24) == -1) failure = .true. ! -x
        if (.not. dd % neighbor_meshbins(2 + 24) == -1) failure = .true. ! +x
        if (.not. dd % neighbor_meshbins(3 + 24) ==  2) failure = .true. ! -y
        if (.not. dd % neighbor_meshbins(4 + 24) == -1) failure = .true. ! +y
        if (.not. dd % neighbor_meshbins(5 + 24) == -1) failure = .true. ! -z
        if (.not. dd % neighbor_meshbins(6 + 24) == -1) failure = .true. ! +z
        ! -z
        if (.not. dd % neighbor_meshbins(1 + 30) == -1) failure = .true. ! -x
        if (.not. dd % neighbor_meshbins(2 + 30) == -1) failure = .true. ! +x
        if (.not. dd % neighbor_meshbins(3 + 30) == -1) failure = .true. ! -y
        if (.not. dd % neighbor_meshbins(4 + 30) == -1) failure = .true. ! +y
        if (.not. dd % neighbor_meshbins(5 + 30) == -1) failure = .true. ! -z
        if (.not. dd % neighbor_meshbins(6 + 30) ==  2) failure = .true. ! +z
        ! +z
        if (.not. dd % neighbor_meshbins(1 + 36) == -1) failure = .true. ! -x
        if (.not. dd % neighbor_meshbins(2 + 36) == -1) failure = .true. ! +x
        if (.not. dd % neighbor_meshbins(3 + 36) == -1) failure = .true. ! -y
        if (.not. dd % neighbor_meshbins(4 + 36) == -1) failure = .true. ! +y
        if (.not. dd % neighbor_meshbins(5 + 36) ==  2) failure = .true. ! -z
        if (.not. dd % neighbor_meshbins(6 + 36) == -1) failure = .true. ! +z
    end select

    if (failure) then
      call write_message("FAILURE: Domain meshbin mapping is incorrect " // &
          "on rank " // trim(to_str(rank)))
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

    call deallocate_dd(dd)
    
  end subroutine test_teardown

end module test_dd_neighbor_meshbins
