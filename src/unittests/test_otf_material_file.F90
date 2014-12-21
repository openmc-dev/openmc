module test_otf_material_file

  use constants,        only: MAX_FILE_LEN
  use global,           only: master, rank, n_procs, n_materials, materials
  use error,            only: warning
  use input_xml,        only: init_otf_materials
  use output,           only: write_message
  use testing_header,   only: TestSuiteClass, TestClass
  use state_point,      only: write_distribmat_comps
  use string,           only: to_str

#ifdef HDF5
  use hdf5
#endif

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

  type(test), public      :: otf_material_file_test

contains

!===============================================================================
! INIT
!===============================================================================

  subroutine test_init(this)

    class(test), intent(inout) :: this

    this % name = "test_otf_material_file"

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

#ifdef MPI
    if (.not. n_procs == 5) then
      if (master) call write_message("Skipping test: must be run 5 MPI ranks")
      call suite % skip(this)
      return
    end if   
#else
    call write_message("Skipping test: must be run with MPI")
    call suite % skip(this)
    return
#endif

#ifndef HDF5
    if (master) call write_message("Skipping test: must be run with HDF5")
    call suite % skip(this)
    return
#endif

    ! Set up a material with some OTF compositions
    n_materials = 1
    allocate(materials(1))
    materials(1) % id = 1
    materials(1) % n_nuclides = 4
    materials(1) % n_comp = 5
    materials(1) % otf_compositions = .true.
    allocate(materials(1) % otf_comp(20))
    select case(rank)
      case(0)
        materials(1) % next_comp_idx = 3
        call materials(1) % reverse_comp_index_map % add_key(1, 5)
        materials(1) % otf_comp(1:4) = (/7.0_8, 7.0_8, 7.0_8, 7.0_8/)
        call materials(1) % reverse_comp_index_map % add_key(2, 2)
        materials(1) % otf_comp(5:8) = (/2.0_8, 1.0_8, 3.0_8, 4.0_8/)
      case(1)
        materials(1) % next_comp_idx = 4
        call materials(1) % reverse_comp_index_map % add_key(1, 3)
        materials(1) % otf_comp(1:4) = (/3.0_8, 2.0_8, 1.0_8, 4.0_8/)
        call materials(1) % reverse_comp_index_map % add_key(2, 1)
        materials(1) % otf_comp(5:8) = (/1.0_8, 2.0_8, 3.0_8, 4.0_8/)
        call materials(1) % reverse_comp_index_map % add_key(3, 2)
        materials(1) % otf_comp(9:12) = (/2.0_8, 1.0_8, 3.0_8, 4.0_8/)
      case(2)
        materials(1) % next_comp_idx = 1
      case(3)
        materials(1) % next_comp_idx = 3
        call materials(1) % reverse_comp_index_map % add_key(1, 1)
        materials(1) % otf_comp(1:4) = (/1.0_8, 2.0_8, 3.0_8, 4.0_8/)
        call materials(1) % reverse_comp_index_map % add_key(2, 4)
        materials(1) % otf_comp(5:8) = (/4.0_8, 3.0_8, 2.0_8, 1.0_8/)
      case(4)
        materials(1) % next_comp_idx = 2
        call materials(1) % reverse_comp_index_map % add_key(1, 3)
        materials(1) % otf_comp(1:4) = (/3.0_8, 2.0_8, 1.0_8, 4.0_8/)
    end select

#ifdef MPI
    call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)
#endif

    ! Write materials to disk
    filename = 'materials.h5'
    call write_distribmat_comps(filename)

#ifdef MPI
    call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)
#endif

    ! Check if the file exists
    INQUIRE(FILE='materials.h5', EXIST=stat)
    if (.not. stat) then
      call suite % fail()
      if (master) call write_message('FAILURE: Test file not created.')
      return
    end if

    ! Return to blank slate
    deallocate(materials)

    ! Now recreate the mat initialization that happens in input_xml
    n_materials = 1
    allocate(materials(n_materials))
    materials(1) % id = 1
    materials(1) % n_nuclides = 4
    materials(1) % n_comp = 5
    materials(1) % otf_compositions = .true.
    materials(1) % comp_file % group = 'mat-1'
    materials(1) % comp_file % n_nuclides = 4
    materials(1) % comp_file % n_instances = 5
    allocate(materials(1) % otf_comp(20))
    call init_otf_materials()

  end subroutine test_setup

!===============================================================================
! EXECUTE
!===============================================================================

  subroutine test_execute()

    real(8) :: density

    ! Read the compositions in
    select case(rank)
      case(0)
        density = materials(1) % get_density(5, 1)
        density = materials(1) % get_density(2, 1)
      case(1)
        density = materials(1) % get_density(3, 1)
        density = materials(1) % get_density(1, 1)
        density = materials(1) % get_density(2, 1)
      case(2)
        ! don't load anything
      case(3)
        density = materials(1) % get_density(1, 1)
        density = materials(1) % get_density(4, 1)
      case(4)
        density = materials(1) % get_density(3, 1)
    end select

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
      case(0)
        if (.not. materials(1) % next_comp_idx == 3) failure = .true.
        if (.not. materials(1) % otf_comp(1) == 7.0_8) failure = .true.
        if (.not. materials(1) % otf_comp(2) == 7.0_8) failure = .true.
        if (.not. materials(1) % otf_comp(3) == 7.0_8) failure = .true.
        if (.not. materials(1) % otf_comp(4) == 7.0_8) failure = .true.
        if (.not. materials(1) % otf_comp(5) == 2.0_8) failure = .true.
        if (.not. materials(1) % otf_comp(6) == 1.0_8) failure = .true.
        if (.not. materials(1) % otf_comp(7) == 3.0_8) failure = .true.
        if (.not. materials(1) % otf_comp(8) == 4.0_8) failure = .true.
      case(1)
        if (.not. materials(1) % next_comp_idx == 4) failure = .true.
        if (.not. materials(1) % otf_comp(1) == 3.0_8) failure = .true.
        if (.not. materials(1) % otf_comp(2) == 2.0_8) failure = .true.
        if (.not. materials(1) % otf_comp(3) == 1.0_8) failure = .true.
        if (.not. materials(1) % otf_comp(4) == 4.0_8) failure = .true.
        if (.not. materials(1) % otf_comp(5) == 1.0_8) failure = .true.
        if (.not. materials(1) % otf_comp(6) == 2.0_8) failure = .true.
        if (.not. materials(1) % otf_comp(7) == 3.0_8) failure = .true.
        if (.not. materials(1) % otf_comp(8) == 4.0_8) failure = .true.
        if (.not. materials(1) % otf_comp(9) == 2.0_8) failure = .true.
        if (.not. materials(1) % otf_comp(10) == 1.0_8) failure = .true.
        if (.not. materials(1) % otf_comp(11) == 3.0_8) failure = .true.
        if (.not. materials(1) % otf_comp(12) == 4.0_8) failure = .true.
      case(2)
        if (.not. materials(1) % next_comp_idx == 1) failure = .true.
      case(3)
        if (.not. materials(1) % next_comp_idx == 3) failure = .true.
        if (.not. materials(1) % otf_comp(1) == 1.0_8) failure = .true.
        if (.not. materials(1) % otf_comp(2) == 2.0_8) failure = .true.
        if (.not. materials(1) % otf_comp(3) == 3.0_8) failure = .true.
        if (.not. materials(1) % otf_comp(4) == 4.0_8) failure = .true.
        if (.not. materials(1) % otf_comp(5) == 4.0_8) failure = .true.
        if (.not. materials(1) % otf_comp(6) == 3.0_8) failure = .true.
        if (.not. materials(1) % otf_comp(7) == 2.0_8) failure = .true.
        if (.not. materials(1) % otf_comp(8) == 1.0_8) failure = .true.
      case(4)
        if (.not. materials(1) % next_comp_idx == 2) failure = .true.
        if (.not. materials(1) % otf_comp(1) == 3.0_8) failure = .true.
        if (.not. materials(1) % otf_comp(2) == 2.0_8) failure = .true.
        if (.not. materials(1) % otf_comp(3) == 1.0_8) failure = .true.
        if (.not. materials(1) % otf_comp(4) == 4.0_8) failure = .true.
    end select

    if (failure) then
      call write_message("FAILURE: Materials file reading failure on rank " // trim(to_str(rank)))
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

  end subroutine test_teardown

end module test_otf_material_file
