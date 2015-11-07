module test_otf_material_file_write

  use constants,        only: MAX_FILE_LEN
  use dd_init,          only: initialize_domain_decomp
  use dd_header,        only: deallocate_dd
  use dd_testing_setup, only: check_procs, dd_simple_four_domains, &
                              dd_setup_four_nuc_five_comp_otf_mats, &
                              dd_create_four_nuc_five_comp_file
  use global,           only: master, rank, n_procs, n_materials, materials, &
                              domain_decomp, dd_run, free_memory
  use error,            only: warning
  use output,           only: write_message
  use state_point,      only: write_distribmat_comps
  use string,           only: to_str
  use testing_header,   only: TestSuiteClass, TestClass

  use hdf5

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

  type(test), public      :: otf_material_file_test_write

contains

!===============================================================================
! INIT
!===============================================================================

  subroutine test_init(this)

    class(test), intent(inout) :: this

    this % name = "test_otf_material_file_write"

  end subroutine test_init

!===============================================================================
! SETUP
!===============================================================================

  subroutine test_setup(this, suite)

    class(test), intent(inout) :: this
    class(TestSuiteClass), intent(inout) :: suite

    character(MAX_FILE_LEN) :: filename
    logical :: stat
#ifdef MPI
    integer :: mpi_err
#endif

    if (check_procs(5)) then
      call suite % skip(this)
      return
    end if

    filename = 'materials.h5'
    call dd_create_four_nuc_five_comp_file(filename)

    ! Check if the file exists
    INQUIRE(FILE=trim(filename), EXIST=stat)
    if (.not. stat) then
      call suite % fail()
      if (master) call write_message('FAILURE: Test file not created.')
      return
    end if

#ifdef MPI
    call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)
#endif

    ! Return to blank slate
    call free_memory()

    ! Now set up OTF DD case

    ! Get generic DD setup with 4 domains for 5 MPI ranks
    call dd_simple_four_domains(domain_decomp)

    ! Initialize dd_type
    call initialize_domain_decomp(domain_decomp)
    dd_run = .true.

    ! Initialize some OTF mats
    call dd_setup_four_nuc_five_comp_otf_mats()

    ! Manually set some materials
    select case(rank)
      case(0)
        materials(1) % next_comp_idx = 3
        call materials(1) % comp_index_map % add_key(5, 1)
        call materials(1) % reverse_comp_index_map % add_key(1, 5)
        materials(1) % otf_comp(:,1) = (/7.0_8, 7.0_8, 7.0_8, 7.0_8/)
        call materials(1) % comp_index_map % add_key(2, 2)
        call materials(1) % reverse_comp_index_map % add_key(2, 2)
        materials(1) % otf_comp(:,2) = (/2.0_8, 1.0_8, 3.0_8, 4.0_8/)
      case(1)
        materials(1) % next_comp_idx = 4
        call materials(1) % comp_index_map % add_key(3, 1)
        call materials(1) % reverse_comp_index_map % add_key(1, 3)
        materials(1) % otf_comp(:,1) = (/3.0_8, 2.0_8, 1.0_8, 4.0_8/)
        call materials(1) % comp_index_map % add_key(1, 2)
        call materials(1) % reverse_comp_index_map % add_key(2, 1)
        materials(1) % otf_comp(:,2) = (/1.0_8, 2.0_8, 3.0_8, 4.0_8/)
        call materials(1) % comp_index_map % add_key(2, 3)
        call materials(1) % reverse_comp_index_map % add_key(3, 2)
        materials(1) % otf_comp(:,3) = (/2.0_8, 1.0_8, 3.0_8, 4.0_8/)
      case(2)
        materials(1) % next_comp_idx = 1
      case(3)
        materials(1) % next_comp_idx = 3
        call materials(1) % comp_index_map % add_key(1, 1)
        call materials(1) % reverse_comp_index_map % add_key(1, 1)
        materials(1) % otf_comp(:,1) = (/1.0_8, 2.0_8, 3.0_8, 4.0_8/)
        call materials(1) % comp_index_map % add_key(4, 2)
        call materials(1) % reverse_comp_index_map % add_key(2, 4)
        materials(1) % otf_comp(:,2) = (/4.0_8, 3.0_8, 2.0_8, 1.0_8/)
      case(4)
        materials(1) % next_comp_idx = 2
        call materials(1) % comp_index_map % add_key(3, 1)
        call materials(1) % reverse_comp_index_map % add_key(1, 3)
        materials(1) % otf_comp(:,1) = (/3.0_8, 2.0_8, 1.0_8, 4.0_8/)
    end select

#ifdef MPI
    call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)
#endif

  end subroutine test_setup

!===============================================================================
! EXECUTE
!===============================================================================

  subroutine test_execute()

    character(MAX_FILE_LEN) :: filename

    ! Write materials to disk
    filename = 'materials-out.h5'
    call write_distribmat_comps(filename)

  end subroutine test_execute

!===============================================================================
! CHECK
!===============================================================================

  subroutine test_check(suite)

    class(TestSuiteClass), intent(inout) :: suite

    logical :: stat

    logical :: failure = .false.
#ifdef MPI
    integer :: mpi_err
    logical :: any_fail
#endif

    ! Check if the file exists
    INQUIRE(FILE='materials-out.h5', EXIST=stat)
    if (.not. stat) then
      call suite % fail()
      if (master) call write_message('FAILURE: Test file not created.')
      return
    end if

    select case(rank)
      case(0, 1)
        ! OTF mat writing should have sorted the composition array in real_bin order
        if (.not. materials(1) % next_comp_idx == 5) failure = .true.
        if (.not. materials(1) % otf_comp(1,1) == 1.0_8) failure = .true.
        if (.not. materials(1) % otf_comp(2,1) == 2.0_8) failure = .true.
        if (.not. materials(1) % otf_comp(3,1) == 3.0_8) failure = .true.
        if (.not. materials(1) % otf_comp(4,1) == 4.0_8) failure = .true.
        if (.not. materials(1) % otf_comp(1,2) == 2.0_8) failure = .true.
        if (.not. materials(1) % otf_comp(2,2) == 1.0_8) failure = .true.
        if (.not. materials(1) % otf_comp(3,2) == 3.0_8) failure = .true.
        if (.not. materials(1) % otf_comp(4,2) == 4.0_8) failure = .true.
        if (.not. materials(1) % otf_comp(1,3) == 3.0_8) failure = .true.
        if (.not. materials(1) % otf_comp(2,3) == 2.0_8) failure = .true.
        if (.not. materials(1) % otf_comp(3,3) == 1.0_8) failure = .true.
        if (.not. materials(1) % otf_comp(4,3) == 4.0_8) failure = .true.
        if (.not. materials(1) % otf_comp(1,4) == 7.0_8) failure = .true.
        if (.not. materials(1) % otf_comp(2,4) == 7.0_8) failure = .true.
        if (.not. materials(1) % otf_comp(3,4) == 7.0_8) failure = .true.
        if (.not. materials(1) % otf_comp(4,4) == 7.0_8) failure = .true.
      case(2)
        if (.not. materials(1) % next_comp_idx == 1) failure = .true.
      case(3)
        if (.not. materials(1) % next_comp_idx == 3) failure = .true.
        if (.not. materials(1) % otf_comp(1,1) == 1.0_8) failure = .true.
        if (.not. materials(1) % otf_comp(2,1) == 2.0_8) failure = .true.
        if (.not. materials(1) % otf_comp(3,1) == 3.0_8) failure = .true.
        if (.not. materials(1) % otf_comp(4,1) == 4.0_8) failure = .true.
        if (.not. materials(1) % otf_comp(1,2) == 4.0_8) failure = .true.
        if (.not. materials(1) % otf_comp(2,2) == 3.0_8) failure = .true.
        if (.not. materials(1) % otf_comp(3,2) == 2.0_8) failure = .true.
        if (.not. materials(1) % otf_comp(4,2) == 1.0_8) failure = .true.
      case(4)
        if (.not. materials(1) % next_comp_idx == 2) failure = .true.
        if (.not. materials(1) % otf_comp(1,1) == 3.0_8) failure = .true.
        if (.not. materials(1) % otf_comp(2,1) == 2.0_8) failure = .true.
        if (.not. materials(1) % otf_comp(3,1) == 1.0_8) failure = .true.
        if (.not. materials(1) % otf_comp(4,1) == 4.0_8) failure = .true.
    end select

    if (failure) then
      call write_message("FAILURE: OTF comp sync or ordering failure on rank " // trim(to_str(rank)))
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
    if (failure) then
      call suite % fail()
    else
      call suite % pass()
    end if
#endif

  end subroutine test_check

!===============================================================================
! TEARDOWN
!===============================================================================

  subroutine test_teardown()

    integer :: stat
#ifdef MPI
    integer :: mpi_err

    call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)
#endif

    call deallocate_dd(domain_decomp)
    call free_memory()
    if (master) then
      open(unit=1234, iostat=stat, file='materials.h5', status='old')
      if (stat.eq.0) close(1234, status='delete')
      open(unit=1234, iostat=stat, file='materials-out.h5', status='old')
      if (stat.eq.0) close(1234, status='delete')
    end if

#ifdef MPI
    call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)
#endif

  end subroutine test_teardown

end module test_otf_material_file_write
