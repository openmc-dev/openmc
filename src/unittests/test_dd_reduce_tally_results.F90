module test_dd_reduce_tally_results

  use dd_init,          only: initialize_domain_decomp
  use dd_header,        only: dd_type, deallocate_dd
  use dd_testing_setup, only: check_procs, dd_simple_four_domains, &
                              dd_simple_four_domain_tallies, &
                              dd_score_to_four_domain_tallies
  use error,            only: warning
  use global,           only: master, domain_decomp, dd_run, free_memory, &
                              tallies, rank, total_weight
  use output,           only: write_message
  use particle_header,  only: Particle
  use string,           only: to_str
  use tally,            only: reduce_tally_results
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

  type(test), public :: dd_reduce_tally_results_test
  
  type(Particle) :: p

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
    call dd_simple_four_domains(domain_decomp)

    ! Initialize dd_type
    call initialize_domain_decomp(domain_decomp)
    dd_run = .true.

    ! Setup tallies and particle for scoring to them
    call dd_simple_four_domain_tallies(p)

    ! Score to the tallies
    call dd_score_to_four_domain_tallies(p)

  end subroutine test_setup

!===============================================================================
! EXECUTE
!===============================================================================

  subroutine test_execute()

    call reduce_tally_results()
    
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
      case(0)
        if (size(tallies(1) % results, 2) /= 3) then
          failure = .true.
        else
          if (.not. int(tallies(1) % results(1,1) % value) == 6 .and. &
              .not. int(tallies(1) % results(1,2) % value) == 4 .and. &
              .not. int(tallies(1) % results(1,3) % value) == 4) failure = .true.
        end if
        ! Check that the bins point to the proper cells
        if (.not. tallies(1) % reverse_filter_index_map % get_key(1) == 2) failure = .true.
        if (.not. tallies(1) % reverse_filter_index_map % get_key(2) == 4) failure = .true.
        if (.not. tallies(1) % reverse_filter_index_map % get_key(3) == 1) failure = .true.
      case(1)
        if (size(tallies(1) % results, 2) /= 3) then
          failure = .true.
        else
          if (.not. int(tallies(1) % results(1,1) % value) == 0 .and. &
              .not. int(tallies(1) % results(1,2) % value) == 0 .and. &
              .not. int(tallies(1) % results(1,2) % value) == 0) failure = .true.
        end if
        ! Check that the bins point to the proper cells
        if (.not. tallies(1) % reverse_filter_index_map % get_key(1) == 2) failure = .true.
        if (.not. tallies(1) % reverse_filter_index_map % get_key(2) == 4) failure = .true.
        if (.not. tallies(1) % reverse_filter_index_map % get_key(3) == 1) failure = .true.
      case(2)
        if (size(tallies(1) % results, 2) /= 2) then
          failure = .true.
        else
          if (.not. int(tallies(1) % results(1,1) % value) == 3 .and. &
              .not. int(tallies(1) % results(1,2) % value) == 4) failure = .true.
        end if
        ! Check that the bins point to the proper cells
        if (.not. tallies(1) % reverse_filter_index_map % get_key(1) == 3) failure = .true.
        if (.not. tallies(1) % reverse_filter_index_map % get_key(2) == 4) failure = .true.
      case(3)
        if (size(tallies(1) % results, 2) /= 1) then
          failure = .true.
        else
          if (.not. int(tallies(1) % results(1,1) % value) == 2) failure = .true.
        end if
        ! Check that the bins point to the proper cells
        if (.not. tallies(1) % reverse_filter_index_map % get_key(1) == 2) failure = .true.
      case(4)
        if (size(tallies(1) % results, 2) /= 2) then
          failure = .true.
        else
          if (.not. int(tallies(1) % results(1,1) % value) == 3 .and. &
              .not. int(tallies(1) % results(1,2) % value) == 4) failure = .true.
        end if
        ! Check that the bins point to the proper cells
        if (.not. tallies(1) % reverse_filter_index_map % get_key(1) == 4) failure = .true.
        if (.not. tallies(1) % reverse_filter_index_map % get_key(2) == 3) failure = .true.
    end select
    
    if (failure) then
      call write_message("FAILURE: Tally reduction failure with OTF " // &
          "tallies on rank " // trim(to_str(rank)))
    end if

    ! Check that the bins point to the proper cells
    select case(rank)
      case(0)
        if (size(tallies(1) % results, 2) /= 3) then
          failure = .true.
        else
          if (.not. tallies(1) % reverse_filter_index_map % get_key(1) == 2) failure = .true.
          if (.not. tallies(1) % reverse_filter_index_map % get_key(2) == 4) failure = .true.
          if (.not. tallies(1) % reverse_filter_index_map % get_key(3) == 1) failure = .true.
        end if
      case(1)
        if (size(tallies(1) % results, 2) /= 3) then
          failure = .true.
        else
          if (.not. tallies(1) % reverse_filter_index_map % get_key(1) == 2) failure = .true.
          if (.not. tallies(1) % reverse_filter_index_map % get_key(2) == 4) failure = .true.
          if (.not. tallies(1) % reverse_filter_index_map % get_key(3) == 1) failure = .true.
        end if
      case(2)
        if (size(tallies(1) % results, 2) /= 2) then
          failure = .true.
        else
          if (.not. tallies(1) % reverse_filter_index_map % get_key(1) == 3) failure = .true.
          if (.not. tallies(1) % reverse_filter_index_map % get_key(2) == 4) failure = .true.
        end if
      case(3)
        if (size(tallies(1) % results, 2) /= 1) then
          failure = .true.
        else
          if (.not. tallies(1) % reverse_filter_index_map % get_key(1) == 2) failure = .true.
        end if
      case(4)
        if (size(tallies(1) % results, 2) /= 2) then
          failure = .true.
        else
          if (.not. tallies(1) % reverse_filter_index_map % get_key(1) == 4) failure = .true.
          if (.not. tallies(1) % reverse_filter_index_map % get_key(2) == 3) failure = .true.
        end if
    end select

    if (failure) then
      call write_message("FAILURE: Tally reduction maps incorrect with OTF " // &
          "tallies on rank " // trim(to_str(rank)))
    end if

    ! Check total_weight
    if (.not. total_weight == 28.0) failure = .true.

    if (failure) then
      call write_message("FAILURE: Incorrect total_weight with OTF " // &
          "tallies on rank " // trim(to_str(rank)))
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

end module test_dd_reduce_tally_results
