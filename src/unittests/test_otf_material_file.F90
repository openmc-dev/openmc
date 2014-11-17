module test_otf_material_file

  use constants,        only: MAX_LINE_LEN
  use global,           only: master
  use error,            only: warning
  use material_header,  only: Composition, CompositionFile
  use output_interface, only: BinaryOutput
  use output,           only: write_message
  use testing_header,   only: TestSuiteClass, TestClass
  use string,           only: to_str

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
  type(Composition)       :: comp(5)
  type(CompositionFile)   :: compfile
  character(MAX_LINE_LEN) :: filename

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

    integer :: i
    real(8), allocatable :: tmp(:)
    logical :: stat
    type(BinaryOutput) :: fh
#ifdef MPI
    integer :: mpi_err
#endif

    filename = 'otf_material_test'
#ifdef HDF5
    filename = trim(filename) // '.h5'
#else
    filename = trim(filename) // '.binary'
#endif

    if (master) then

      ! Create file and write header
      call fh % file_create(filename, record_len = 8)
      call fh % write_data(4, 'n_nuclides', record=1)
      call fh % write_data(5, 'n_instances', record=2)
      call fh % file_close()
      
      ! Open the file for composition writing
      call fh % file_open(filename, 'w', serial = .true., direct_access = .true., record_len = 8*4)

#ifdef HDF5
      ! For HDF5, we need to create the whole dataset first, then write to sections with hyperslabs
      allocate(tmp(20))
      call fh % write_data(tmp, "comps", length=20)
      deallocate(tmp)
      i = 1
#else
      i = 2
#endif

      ! Write the compositions
      call fh % write_data((/4.0_8, 3.0_8, 2.0_8, 1.0_8/), "comps", length=4, record=i+3)
      call fh % write_data((/7.0_8, 7.0_8, 7.0_8, 7.0_8/), "comps", length=4, record=i+4)
      call fh % write_data((/1.0_8, 2.0_8, 3.0_8, 4.0_8/), "comps", length=4, record=i)
      call fh % write_data((/3.0_8, 2.0_8, 1.0_8, 4.0_8/), "comps", length=4, record=i+2)
      call fh % write_data((/2.0_8, 1.0_8, 3.0_8, 4.0_8/), "comps", length=4, record=i+1)

      ! Close the file
      call fh % file_close()

    end if

#ifdef MPI
    call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)
#endif

    INQUIRE(FILE=filename, EXIST=stat)
    if (.not. stat) then
      call suite % skip(this)
      if (master) call write_message('Test file not created: ' // filename)
      return
    end if

  end subroutine test_setup

!===============================================================================
! EXECUTE
!===============================================================================

  subroutine test_execute()

    integer :: i
    type(BinaryOutput) :: fh

    ! Set up the composition file wrapper
    compfile % path = filename
    call fh % file_open(filename, 'r', &
                        direct_access = .true., record_len = 8)
    call fh % read_data(compfile % n_nuclides, 'n_nuclides', record = 1)
    call fh % read_data(compfile % n_instances, 'n_instances', record = 2)
    call fh % file_close()

    ! Read the compositions in
    do i = 1, compfile % n_instances
      comp(i) = compfile % load(i)
    end do

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

    if (.not. compfile % n_nuclides == 4) failure = .true.
    if (.not. compfile % n_instances == 5) failure = .true.
    
    if (.not. comp(1) % atom_density(1) == 1.0_8) failure = .true.
    if (.not. comp(1) % atom_density(2) == 2.0_8) failure = .true.
    if (.not. comp(1) % atom_density(3) == 3.0_8) failure = .true.
    if (.not. comp(1) % atom_density(4) == 4.0_8) failure = .true.

    if (.not. comp(2) % atom_density(1) == 2.0_8) failure = .true.
    if (.not. comp(2) % atom_density(2) == 1.0_8) failure = .true.
    if (.not. comp(2) % atom_density(3) == 3.0_8) failure = .true.
    if (.not. comp(2) % atom_density(4) == 4.0_8) failure = .true.
   
    if (.not. comp(3) % atom_density(1) == 3.0_8) failure = .true.
    if (.not. comp(3) % atom_density(2) == 2.0_8) failure = .true.
    if (.not. comp(3) % atom_density(3) == 1.0_8) failure = .true.
    if (.not. comp(3) % atom_density(4) == 4.0_8) failure = .true.

    if (.not. comp(4) % atom_density(1) == 4.0_8) failure = .true.
    if (.not. comp(4) % atom_density(2) == 3.0_8) failure = .true.
    if (.not. comp(4) % atom_density(3) == 2.0_8) failure = .true.
    if (.not. comp(4) % atom_density(4) == 1.0_8) failure = .true.

    if (.not. comp(5) % atom_density(1) == 7.0_8) failure = .true.
    if (.not. comp(5) % atom_density(2) == 7.0_8) failure = .true.
    if (.not. comp(5) % atom_density(3) == 7.0_8) failure = .true.
    if (.not. comp(5) % atom_density(4) == 7.0_8) failure = .true.

    if (failure) then
      call write_message("FAILURE: Materials file reading failure.")
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
