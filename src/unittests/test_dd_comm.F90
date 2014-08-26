module test_dd_comm

  use dd_comm,          only: distribute_source, synchronize_transfer_info
  use dd_init,          only: initialize_domain_decomp
  use dd_header,        only: dd_type, deallocate_dd
  use dd_tracking,      only: cross_domain_boundary
  use error,            only: fatal_error, warning
  use global,           only: master, n_procs, rank, message, work, &
                              n_particles, source_bank, size_source_bank
  use output,           only: header, write_message
  use random_lcg,       only: initialize_prng, set_particle_seed, prn_seed
  use string,           only: to_str
  use test_dd_init,     only: check_procs, dd_simple_four_domains

#ifdef MPI
  use mpi
#endif

  implicit none
  public

  type(dd_type) :: dd

contains

!===============================================================================
! TEST_DISTRIBUTE_SOURCE
!===============================================================================

  subroutine test_distribute_source()

    integer(8) :: i
    logical :: failure = .false.
    integer :: mpi_err

    if (master) call header("test_distribute_source", level=2)

    if (check_procs()) return

    ! SETUP

    if (master) then
      message = "Setting up..."
      call write_message(1)
    end if

    ! Get generic DD setup with 4 domains for 5 MPI ranks
    call dd_simple_four_domains(dd)

    ! Initialiaze dd type
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

    ! EXECUTE

    if (master) then
      message = "Invoking test..."
      call write_message(1)
    end if

    ! Invoke test method
    call distribute_source(dd)

    ! CHECK
    
    if (master) then
      message = "Checking results..."
      call write_message(1)
    end if

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
      message = "FAILED: Rank " // trim(to_str(rank)) // &
          " doesn't have all the particles it should."
      call fatal_error()
    end if

    call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)
    
    if (master) then
      message = "PASSED"
      call write_message(1)
    end if

    ! Clean up
    deallocate(source_bank)    
    call deallocate_dd(dd)

  end subroutine test_distribute_source

!===============================================================================
! TEST_SYNCHRONIZE_INFO depends on test_distribute_source
!===============================================================================

  subroutine test_synchronize_transfer_info()

    logical :: failure = .false.
    integer :: bin
    integer :: mpi_err

    if (master) call header("test_synchronize_transfer_info", level=2)

    if (check_procs()) return

    ! EXECUTE

    if (master) then
      message = "Setting up..."
      call write_message(1)
    end if

    ! Get generic DD setup with 4 domains for 5 MPI ranks
    call dd_simple_four_domains(dd)
    ! TODO: this should test the whole neighborhood: 4 domains is not enough

    ! Initialize dd_type
    call initialize_domain_decomp(dd)

    ! Set up n_scatters_local
    dd % n_scatters_local = 0
    select case(rank)
      case (0)
        dd % n_scatters_local(dd % bins_dict % get_key(2)) = 7 ! 7 from 1 -> 2
        dd % n_scatters_local(dd % bins_dict % get_key(3)) = 5 ! 5 from 1 -> 3
      case (1)
        dd % n_scatters_local(dd % bins_dict % get_key(2)) = 2 ! 2 from 1 -> 2
        dd % n_scatters_local(dd % bins_dict % get_key(3)) = 9 ! 9 from 1 -> 3
      case (2)
        dd % n_scatters_local(dd % bins_dict % get_key(1)) = 0 ! 0 from 2 -> 1
        dd % n_scatters_local(dd % bins_dict % get_key(4)) = 3 ! 3 from 2 -> 4
      case (3)
        dd % n_scatters_local(dd % bins_dict % get_key(1)) = 2 ! 2 from 3 -> 1
        dd % n_scatters_local(dd % bins_dict % get_key(4)) = 1 ! 1 from 3 -> 4
      case (4)
        dd % n_scatters_local(dd % bins_dict % get_key(2)) = 0 ! 0 from 4 -> 2
        dd % n_scatters_local(dd % bins_dict % get_key(3)) = 9 ! 0 from 4 -> 3
    end select

    ! EXECUTE

    if (master) then
      message = "Invoking test..."
      call write_message(1)
    end if

    ! Invoke test method
    call synchronize_transfer_info(dd)

    ! CHECK

    if (master) then
      message = "Checking RESULTS"
      call write_message(1)
    end if

    select case(rank)
      case(0, 1)
        bin = dd % bins_dict % get_key(2)
        if (any(dd % n_scatters_neighborhood(:, bin) /= 0)) failure = .true.
        bin = dd % bins_dict % get_key(3)
        if (dd % n_scatters_neighborhood(16, bin) /= 9) failure = .true.
      case(2)
        bin = dd % bins_dict % get_key(1)
        if (dd % n_scatters_neighborhood(20, bin) /= 2) failure = .true.
        bin = dd % bins_dict % get_key(4)
        if (dd % n_scatters_neighborhood(15, bin) /= 1) failure = .true.
      case(3)
        bin = dd % bins_dict % get_key(1)
        if (any(dd % n_scatters_neighborhood(:, bin) /= 0)) failure = .true.
        bin = dd % bins_dict % get_key(4)
        if (dd % n_scatters_neighborhood(25, bin) /= 3) failure = .true.
      case(4)
        bin = dd % bins_dict % get_key(2)
        if (dd % n_scatters_neighborhood(9, bin) /= 9) failure = .true.
        bin = dd % bins_dict % get_key(3)
        if (dd % n_scatters_neighborhood(19, bin) /= 14) failure = .true.
    end select    

    if (failure) then
      message = "FAILED: Rank " // trim(to_str(rank)) // &
          " didn't receive all the info it should have about its neighbors."
      call fatal_error()
    end if

    call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)
    
    if (master) then
      message = "PASSED"
      call write_message(1)
    end if

    ! Clean up
    call deallocate_dd(dd)


  end subroutine test_synchronize_transfer_info

end module test_dd_comm
