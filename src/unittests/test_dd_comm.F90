module test_dd_comm

  use constants
  use dd_comm,          only: distribute_source, synchronize_transfer_info, &
                              synchronize_destination_info
  use dd_init,          only: initialize_domain_decomp
  use dd_header,        only: dd_type, deallocate_dd
  use dd_tracking,      only: cross_domain_boundary
  use error,            only: fatal_error, warning
  use global,           only: master, n_procs, rank, message, work, &
                              n_particles, source_bank, size_source_bank
  use output,           only: header
  use random_lcg,       only: initialize_prng, set_particle_seed, prn_seed
  use string,           only: to_str
  use test_dd_init,     only: check_procs, dd_simple_four_domains
  use testing_header,   only: testing_type

#ifdef MPI
  use mpi
#endif

  implicit none
  public

  type(dd_type), save :: dd

contains

!===============================================================================
! TEST_DISTRIBUTE_SOURCE
!===============================================================================

  subroutine test_distribute_source(suite)

    type(testing_type), intent(inout) :: suite
    
    integer(8) :: i
    logical :: failure = .false.
#ifdef MPI
    integer :: mpi_err
    logical :: any_fail
#endif

    if (master) call header("test_distribute_source", level=2)

    if (check_procs()) then
      if (master) call suite % skip()
      return
    end if

    ! SETUP

    if (master) call suite % setup()

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

    if (master) call suite % execute()

    ! Invoke test method
    call distribute_source(dd)

    ! CHECK
    
    if (master) call suite % check()

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
      call suite % fail("Rank " // trim(to_str(rank)) // &
        " doesn't have all the particles it should.")
    end if

#ifdef MPI
    call MPI_ALLREDUCE(failure, any_fail, 1, MPI_LOGICAL, MPI_LOR, &
        MPI_COMM_WORLD, mpi_err)
#endif
    
    if (master .and. .not. any_fail) call suite % pass()

    ! Clean up
    deallocate(source_bank)    
    call deallocate_dd(dd)

  end subroutine test_distribute_source

!===============================================================================
! DD_SIMPLE_FOUR_DOMAIN_SCATTERS hardcodes in a few particles to transfer
!===============================================================================

  subroutine dd_simple_four_domain_scatters(dd)
  
    type(dd_type), intent(inout) :: dd
    integer :: bin

    ! Initialize particle buffer
    dd % particle_buffer % outscatter_destination = NO_OUTSCATTER

    dd % n_scatters_local = 0
    select case(rank)
      case (0)
        bin = dd % bins_dict % get_key(2)
        dd % n_scatters_local(bin) =  7 !  7 from 1 -> 2
        dd % particle_buffer(2) % outscatter_destination = bin
        dd % particle_buffer(5) % outscatter_destination = bin
        dd % particle_buffer(7) % outscatter_destination = bin
        dd % particle_buffer(8) % outscatter_destination = bin
        dd % particle_buffer(9) % outscatter_destination = bin
        dd % particle_buffer(11) % outscatter_destination = bin
        dd % particle_buffer(12) % outscatter_destination = bin
        bin = dd % bins_dict % get_key(3)
        dd % n_scatters_local(bin) =  5 !  5 from 1 -> 3
        dd % particle_buffer(1) % outscatter_destination = bin
        dd % particle_buffer(3) % outscatter_destination = bin
        dd % particle_buffer(4) % outscatter_destination = bin
        dd % particle_buffer(10) % outscatter_destination = bin
        dd % particle_buffer(14) % outscatter_destination = bin
      case (1)
        bin = dd % bins_dict % get_key(2)
        dd % n_scatters_local(bin) =  2 !  2 from 1 -> 2
        dd % particle_buffer(2) % outscatter_destination = bin
        dd % particle_buffer(13) % outscatter_destination = bin
        bin = dd % bins_dict % get_key(3)
        dd % n_scatters_local(bin) =  9 !  9 from 1 -> 3
        dd % particle_buffer(1) % outscatter_destination = bin
        dd % particle_buffer(3) % outscatter_destination = bin
        dd % particle_buffer(4) % outscatter_destination = bin
        dd % particle_buffer(5) % outscatter_destination = bin
        dd % particle_buffer(7) % outscatter_destination = bin
        dd % particle_buffer(8) % outscatter_destination = bin
        dd % particle_buffer(9) % outscatter_destination = bin
        dd % particle_buffer(12) % outscatter_destination = bin
        dd % particle_buffer(14) % outscatter_destination = bin
      case (2)
        bin = dd % bins_dict % get_key(1)
        dd % n_scatters_local(bin) = 12 ! 12 from 2 -> 1
        dd % particle_buffer(3) % outscatter_destination = bin
        dd % particle_buffer(4) % outscatter_destination = bin
        dd % particle_buffer(6) % outscatter_destination = bin
        dd % particle_buffer(7) % outscatter_destination = bin
        dd % particle_buffer(13) % outscatter_destination = bin
        dd % particle_buffer(15) % outscatter_destination = bin
        dd % particle_buffer(17) % outscatter_destination = bin
        dd % particle_buffer(18) % outscatter_destination = bin
        dd % particle_buffer(19) % outscatter_destination = bin
        dd % particle_buffer(21) % outscatter_destination = bin
        dd % particle_buffer(22) % outscatter_destination = bin
        dd % particle_buffer(24) % outscatter_destination = bin
        bin = dd % bins_dict % get_key(4)
        dd % n_scatters_local(bin) =  3 !  3 from 2 -> 4
        dd % particle_buffer(8) % outscatter_destination = bin
        dd % particle_buffer(16) % outscatter_destination = bin
        dd % particle_buffer(20) % outscatter_destination = bin
      case (3)
        bin = dd % bins_dict % get_key(1)
        dd % n_scatters_local(bin) =  7 !  7 from 3 -> 1
        dd % particle_buffer(1) % outscatter_destination = bin
        dd % particle_buffer(4) % outscatter_destination = bin
        dd % particle_buffer(5) % outscatter_destination = bin
        dd % particle_buffer(7) % outscatter_destination = bin
        dd % particle_buffer(8) % outscatter_destination = bin
        dd % particle_buffer(19) % outscatter_destination = bin
        dd % particle_buffer(26) % outscatter_destination = bin
        bin = dd % bins_dict % get_key(4)
        dd % n_scatters_local(bin) =  1 !  1 from 3 -> 4
        dd % particle_buffer(14) % outscatter_destination = bin
      case (4)
        bin = dd % bins_dict % get_key(2)
        dd % n_scatters_local(bin) =  0 !  0 from 4 -> 2
        bin = dd % bins_dict % get_key(3)
        dd % n_scatters_local(bin) =  9 !  9 from 4 -> 3
        dd % particle_buffer(14) % outscatter_destination = bin
        dd % particle_buffer(15) % outscatter_destination = bin
        dd % particle_buffer(16) % outscatter_destination = bin
        dd % particle_buffer(19) % outscatter_destination = bin
        dd % particle_buffer(21) % outscatter_destination = bin
        dd % particle_buffer(23) % outscatter_destination = bin
        dd % particle_buffer(24) % outscatter_destination = bin
        dd % particle_buffer(25) % outscatter_destination = bin
        dd % particle_buffer(27) % outscatter_destination = bin
    end select

  end subroutine dd_simple_four_domain_scatters

!===============================================================================
! TEST_SYNCHRONIZE_INFO
!===============================================================================

  subroutine test_synchronize_transfer_info(suite)

    type(testing_type), intent(inout) :: suite

    logical :: failure = .false.
    integer :: bin
#ifdef MPI
    integer :: mpi_err
    logical :: any_fail
#endif

    if (master) call header("test_synchronize_transfer_info", level=2)

    if (check_procs()) then
      if (master) call suite % skip()
      return
    end if

    ! SETUP

    if (master) call suite % setup()

    ! Get generic DD setup with 4 domains for 5 MPI ranks
    call dd_simple_four_domains(dd)
    ! TODO: this should test the whole neighborhood: 4 domains is not enough

    ! Initialize dd_type
    call initialize_domain_decomp(dd)

    ! Set up local outscatter
    call dd_simple_four_domain_scatters(dd)

    ! EXECUTE

    if (master) call suite % execute()

    ! Invoke test method
    call synchronize_transfer_info(dd)

    ! CHECK

    if (master) call suite % check()

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
      call suite % fail("Rank " // trim(to_str(rank)) // &
          " didn't receive all the info it should have about its neighbors.")
    end if

#ifdef MPI
    call MPI_ALLREDUCE(failure, any_fail, 1, MPI_LOGICAL, MPI_LOR, &
        MPI_COMM_WORLD, mpi_err)
#endif
    
    if (master .and. .not. any_fail) call suite % pass()

    ! Clean up
    call deallocate_dd(dd)

  end subroutine test_synchronize_transfer_info

!===============================================================================
! TEST_SYNCHRONIZE_DESTINATION_INFO
!===============================================================================

  subroutine test_synchronize_destination_info(suite)

    type(testing_type), intent(inout) :: suite

    logical :: failure = .false.
#ifdef MPI
    integer :: mpi_err
    logical :: any_fail
#endif
    integer :: to_bin, pr_bin

    if (master) call header("test_synchronize_destination_info", level=2)

    if (check_procs()) then
      if (master) call suite % skip()
      return
    end if

    ! SETUP

    if (master) call suite % setup()

    ! Get generic DD setup with 4 domains for 5 MPI ranks
    call dd_simple_four_domains(dd)
    ! TODO: this should test the whole neighborhood: 4 domains is not enough

    ! Initialize dd_type
    call initialize_domain_decomp(dd)

    ! Setup particle buffer
    dd % particle_buffer % outscatter_destination = NO_OUTSCATTER

    ! Set up local outscatter
    call dd_simple_four_domain_scatters(dd)

    ! Get dd % n_scatters_neighborhood
    call synchronize_transfer_info(dd)

    ! EXECUTE

    if (master) call suite % execute()

    call synchronize_destination_info(dd)

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

    if (master) call suite % check()

    if (failure) then
      call suite % fail("Rank " // trim(to_str(rank)) // " calculated the " // &
                " wrong number of particles to send to a process on a " // &
                "neighboring domain.")
    end if

#ifdef MPI
    call MPI_ALLREDUCE(failure, any_fail, 1, MPI_LOGICAL, MPI_LOR, &
        MPI_COMM_WORLD, mpi_err)
#endif
    
    if (master .and. .not. any_fail) call suite % pass()

    ! Clean up
    call deallocate_dd(dd)

  end subroutine test_synchronize_destination_info

!===============================================================================
! TEST_SEND_RECV_PARTICLES
!===============================================================================

  subroutine test_send_recv_particles(suite)

    type(testing_type), intent(inout) :: suite

    logical :: failure = .false.
#ifdef MPI
    integer :: mpi_err
    logical :: any_fail
#endif

    if (master) call header("test_send_recv_particles", level=2)

    if (master) then
      message = "TEST NOT IMPLEMENTED"
      call warning(force=.true.)
      call suite % skip()
    end if
    return

    if (check_procs()) then
      if (master) call suite % skip()
      return
    end if

    ! SETUP

    if (master) call suite % setup()

    ! Get generic DD setup with 4 domains for 5 MPI ranks
    call dd_simple_four_domains(dd)
    ! TODO: this should test the whole neighborhood: 4 domains is not enough

    ! Set n_particles so the particle buffer is properly allocated
    n_particles = 97 ! we have 55 transferring

    ! Initialize dd_type
    call initialize_domain_decomp(dd)

    ! Set up local outscatter
    call dd_simple_four_domain_scatters(dd)

    ! Get dd % n_scatters_neighborhood
    call synchronize_transfer_info(dd)

    ! Set renc/recv info
    call synchronize_destination_info(dd)

    ! EXECUTE

    if (master) call suite % execute()

    ! CHECK

    if (master) call suite % check()

#ifdef MPI
    call MPI_ALLREDUCE(failure, any_fail, 1, MPI_LOGICAL, MPI_LOR, &
        MPI_COMM_WORLD, mpi_err)
#endif
    
    if (master .and. .not. any_fail) call suite % pass()

    ! Clean up
    call deallocate_dd(dd)

  end subroutine test_send_recv_particles

!===============================================================================
! TEST_SYNCHRONIZE_BANK_DD
!===============================================================================

  subroutine test_synchronize_bank_dd(suite)

    type(testing_type), intent(inout) :: suite

    logical :: failure = .false.
#ifdef MPI
    integer :: mpi_err
    logical :: any_fail
#endif

    if (master) call header("test_synchronize_bank_dd", level=2)

    if (master) then
      message = "TEST NOT IMPLEMENTED"
      call warning(force=.true.)
      call suite % skip()
    end if
    return

    if (check_procs()) then
      if (master) call suite % skip()
      return
    end if

    ! SETUP

    if (master) call suite % setup()

    ! Get generic DD setup with 4 domains for 5 MPI ranks
    call dd_simple_four_domains(dd)
    ! TODO: this should test the whole neighborhood: 4 domains is not enough

    ! Set n_particles so the particle buffer is properly allocated
    n_particles = 97 ! we have 55 transferring

    ! Initialize dd_type
    call initialize_domain_decomp(dd)

    ! Set up local outscatter
    call dd_simple_four_domain_scatters(dd)

    ! Get dd % n_scatters_neighborhood
    call synchronize_transfer_info(dd)

    ! Set renc/recv info
    call synchronize_destination_info(dd)

    ! EXECUTE

    if (master) call suite % execute()

    ! CHECK

    if (master) call suite % check()

#ifdef MPI
    call MPI_ALLREDUCE(failure, any_fail, 1, MPI_LOGICAL, MPI_LOR, &
        MPI_COMM_WORLD, mpi_err)
#endif
    
    if (master .and. .not. any_fail) call suite % pass()

    ! Clean up
    call deallocate_dd(dd)

  end subroutine test_synchronize_bank_dd

end module test_dd_comm
