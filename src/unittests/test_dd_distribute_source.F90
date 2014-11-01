module test_dd_distribute_source

  use dd_comm,          only: distribute_source
  use dd_init,          only: initialize_domain_decomp
  use dd_header,        only: dd_type, deallocate_dd
  use dd_testing_setup, only: check_procs, dd_simple_four_domains
  use global,           only: master, rank, work, message, n_particles, &
                              source_bank, size_source_bank
  use output,           only: write_message
  use random_lcg,       only: initialize_prng, set_particle_seed, prn_seed
  use string,           only: to_str
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

  type(test), public :: dd_distribute_source_test
  
  type(dd_type) :: dd

contains

!===============================================================================
! INIT
!===============================================================================

  subroutine test_init(this)

    class(test), intent(inout) :: this

    this % name = "test_dd_distribute_source"
    
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
    call initialize_domain_decomp(dd)
    
    ! Initialize some particles on all ranks, putting some on all domains
    work = 5
    size_source_bank = 5
    allocate(source_bank(5))
    source_bank(1) % wgt = real(rank*5 + 1, 8) ! set for tracking
    source_bank(1) % xyz = (/  1,  1, 0/)
    source_bank(2) % wgt = real(rank*5 + 2, 8)
    source_bank(2) % xyz = (/  1, -1, 0/)
    source_bank(3) % wgt = real(rank*5 + 3, 8)
    source_bank(3) % xyz = (/ -1,  1, 0/)
    source_bank(4) % wgt = real(rank*5 + 4, 8)
    source_bank(4) % xyz = (/ -1, -1, 0/)
    source_bank(5) % wgt = real(rank*5 + 5, 8)
    source_bank(5) % xyz = (/  1,  1, 0/)
    n_particles = 25

    ! This uses random numbers, so we need to initialize those as well
    call initialize_prng()
    call set_particle_seed(int(rank*5 + 1, 8))
    source_bank(1) % prn_seed = prn_seed
    call set_particle_seed(int(rank*5 + 2, 8))
    source_bank(2) % prn_seed = prn_seed
    call set_particle_seed(int(rank*5 + 3, 8))
    source_bank(3) % prn_seed = prn_seed
    call set_particle_seed(int(rank*5 + 4, 8))
    source_bank(4) % prn_seed = prn_seed
    call set_particle_seed(int(rank*5 + 5, 8))
    source_bank(5) % prn_seed = prn_seed

  end subroutine test_setup

!===============================================================================
! EXECUTE
!===============================================================================

  subroutine test_execute()

    ! Invoke test method
    call distribute_source(dd)
    
  end subroutine test_execute

!===============================================================================
! CHECK
!===============================================================================

  subroutine test_check(suite)

    class(TestSuiteClass), intent(inout) :: suite

    integer(8) :: i
    logical :: failure = .false.
#ifdef MPI
    integer :: mpi_err
    logical :: any_fail
#endif

    select case(rank)
      case (0)
        if (.not. work == 3) failure = .true.
        do i = 1, work
          if (.not. int(source_bank(i) % wgt) ==  4 .and. &
              .not. int(source_bank(i) % wgt) == 14 .and. &
              .not. int(source_bank(i) % wgt) == 24) failure = .true.
        end do
      case (1)
        if (.not. work == 2) failure = .true.
        do i = 1, work
          if (.not. int(source_bank(i) % wgt) ==  9 .and. &
              .not. int(source_bank(i) % wgt) == 19) failure = .true.
        end do
      case (2)
        if (.not. work == 5) failure = .true.
        do i = 1, work
          if (.not. int(source_bank(i) % wgt) == 13 .and. &
              .not. int(source_bank(i) % wgt) ==  3 .and. &
              .not. int(source_bank(i) % wgt) ==  8 .and. &
              .not. int(source_bank(i) % wgt) ==  23 .and. &
              .not. int(source_bank(i) % wgt) == 18) failure = .true.
        end do
      case (3)
        if (.not. work == 5) failure = .true.
        do i = 1, work
          if (.not. int(source_bank(i) % wgt) == 17 .and. &
              .not. int(source_bank(i) % wgt) == 12 .and. &
              .not. int(source_bank(i) % wgt) == 22 .and. &
              .not. int(source_bank(i) % wgt) ==  2 .and. &
              .not. int(source_bank(i) % wgt) ==  7) failure = .true.
        end do
      case (4)
        if (.not. work == 10) failure = .true.
        do i = 1, work
          if (.not. int(source_bank(i) % wgt) == 21 .and. &
              .not. int(source_bank(i) % wgt) == 25 .and. &
              .not. int(source_bank(i) % wgt) ==  1 .and. &
              .not. int(source_bank(i) % wgt) ==  6 .and. &
              .not. int(source_bank(i) % wgt) == 11 .and. &
              .not. int(source_bank(i) % wgt) == 15 .and. &
              .not. int(source_bank(i) % wgt) ==  5 .and. &
              .not. int(source_bank(i) % wgt) == 10 .and. &
              .not. int(source_bank(i) % wgt) == 16 .and. &
              .not. int(source_bank(i) % wgt) == 20) failure = .true.
        end do
    end select
    
    if (failure) then
      message = "FAILURE: Rank " // trim(to_str(rank)) // &
                " doesn't have all the particles it should."
      call write_message()
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

    deallocate(source_bank)
    call deallocate_dd(dd)
    
  end subroutine test_teardown

end module test_dd_distribute_source
