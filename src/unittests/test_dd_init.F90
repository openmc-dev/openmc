module test_dd_init

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

#ifdef MPI
  use mpi
#endif

  implicit none
  public

  type(dd_type) :: dd

contains

!===============================================================================
! DD_SIMPLE_FOUR_DOMAINS hardcodes in a simple 4-domain xml input
!===============================================================================

  subroutine dd_simple_four_domains(dd)
  
    type(dd_type), intent(inout) :: dd
  
    ! Mimic everything done in input_xml for DD setup
    allocate(dd % mesh)
    allocate(dd % mesh % dimension(3))
    allocate(dd % mesh % lower_left(3))
    allocate(dd % mesh % upper_right(3))
    allocate(dd % mesh % width(3))
    dd % mesh % n_dimension = 3

    dd % mesh % dimension = (/ 2, 2, 1/)
    dd % n_domains = 4
    dd % mesh % lower_left = (/ -2, -2, -1/)
    dd % mesh % upper_right = (/ 2, 2, 1/)
    dd % mesh % width = (/ 2, 2, 2/)
  
    allocate(dd % domain_load_dist(4))
    dd % domain_load_dist = (/ 2.0, 1.0, 1.0, 1.0/)

  end subroutine dd_simple_four_domains

!===============================================================================
! CHECK_PROCS makes sure we have the right number of processors for this test
!===============================================================================

  function check_procs() result(skip)
  
    logical :: skip
    skip = .false.
  
#ifdef MPI
    if (.not. n_procs == 5) then
      message = "Skipping test_synchronize_info: must be run with MPI 5 procs"
      call warning()
      skip = .true.
    end if
#else
    message = "Skipping test_synchronize_info: requires MPI"
    call warning()
    skip = .true.
#endif

  end function check_procs

!===============================================================================
! TEST_SET_NEIGHBOR_MESHBINS
!===============================================================================

  subroutine test_set_neighbor_meshbins()

    logical :: failure = .false.
    integer :: mpi_err

    if (master) call header("test_set_neighbor_meshbins", level=2)

    if (check_procs()) return

    ! SETUP

    if (master) then
      message = "Setting up..."
      call write_message(1)
    end if

    ! Get generic DD setup with 4 domains for 5 MPI ranks
    call dd_simple_four_domains(dd)

    ! EXECUTE

    if (master) then
      message = "Invoking test..."
      call write_message(1)
    end if

    ! Invoke test method

    call initialize_domain_decomp(dd)

    ! CHECK
    
    if (master) then
      message = "Checking results..."
      call write_message(1)
    end if

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
      message = "FAILED: domain meshbin mapping is incorrect"
      call fatal_error()
    end if

#ifdef MPI
    call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)
#endif

    if (master) then
      message = "PASSED"
      call write_message(1)
    end if

    ! Clean up
    call deallocate_dd(dd)

  end subroutine test_set_neighbor_meshbins


!===============================================================================
! TEST_BINS_DICT
!===============================================================================

  subroutine test_bins_dict()

    logical :: failure = .false.
    integer :: mpi_err

    if (master) call header("test_bins_dict", level=2)

    if (check_procs()) return

    ! SETUP

    if (master) then
      message = "Setting up..."
      call write_message(1)
    end if

    ! Get generic DD setup with 4 domains for 5 MPI ranks
    call dd_simple_four_domains(dd)

    ! EXECUTE

    if (master) then
      message = "Invoking test..."
      call write_message(1)
    end if

    ! Invoke test method

    call initialize_domain_decomp(dd)

    ! CHECK
    
    if (master) then
      message = "Checking results..."
      call write_message(1)
    end if

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
      message = "FAILED: domain bins_dict is incorrect"
      call fatal_error()
    end if

#ifdef MPI
    call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)
#endif
    
    if (master) then
      message = "PASSED"
      call write_message(1)
    end if

    ! Clean up
    call deallocate_dd(dd)

  end subroutine test_bins_dict

end module test_dd_init
