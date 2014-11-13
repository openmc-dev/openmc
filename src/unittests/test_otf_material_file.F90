module test_otf_material_file

  use global,           only: master
  use error,            only: warning
  use testing_header,   only: TestSuiteClass, TestClass
  use output_interface, only: BinaryOutput

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

  type(test), public :: otf_material_file_test

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

    real(8) :: buff(4)
    type(BinaryOutput) :: fh
    
    if (master) then
!      open(unit=666,file='/tmp/otf_material_test.binary',access='direct',recl=8,status='replace',action='write')
!      write(unit=666,rec=1) 1.0_8
!      write(unit=666,rec=2) 2.0_8
!      write(unit=666,rec=3) 3.0_8
!      write(unit=666,rec=4) 4.0_8
!      write(unit=666,rec=5) 5.0_8
!      close(unit=666)
      call fh % file_create('/tmp/otf_material_test.binary', record_len = 8)
      call fh % write_data(4, 'n_nuclides', record=1)
      call fh % write_data(5, 'n_instances', record=2)
      call fh % write_data(1.0_8, '1', record=3)
      call fh % write_data(2.0_8, '1', record=4)
      call fh % write_data(3.0_8, '1', record=5)
      call fh % write_data(4.0_8, '1', record=6)
      
!      call fh % write_data((/1.0_8, 2.0_8, 3.0_8, 4.0_8/), "1", length=4, record=2)
!      call fh % write_data((/2.0_8, 1.0_8, 3.0_8, 4.0_8/), "2", length=4, record=3)
!      call fh % write_data((/3.0_8, 2.0_8, 1.0_8, 4.0_8/), "3", length=4, record=4)
!      call fh % write_data((/4.0_8, 3.0_8, 2.0_8, 1.0_8/), "4", length=4, record=5)
!      call fh % write_data((/7.0_8, 7.0_8, 7.0_8, 7.0_8/), "5", length=4, record=6)
      call fh % file_close()
    end if

  end subroutine test_setup

!===============================================================================
! EXECUTE
!===============================================================================

  subroutine test_execute()

    real(8) :: buff(4)
    integer :: a, b

    type(BinaryOutput) :: fh

    if (master) then
!      open(unit=666,file='/tmp/otf_material_test.binary',access='direct',recl=8,action='read')
!      read(unit=666,rec=5)b
!      print *, 'buff', b
!      read(unit=666,rec=1)b
!      print *, 'buff', b
!      read(unit=666,rec=3)b
!      print *, 'buff', b
!      read(unit=666,rec=2)b
!      print *, 'buff', b
!      close(unit=666)
      call fh % file_open('/tmp/otf_material_test.binary', 'r', &
                          direct_access = .true., record_len = 8)
      call fh % read_data(a, 'n_nuclides', record = 1)
      call fh % read_data(b, 'n_instances', record = 2)
      call fh % read_data(buff(1), '1', record = 3)
      call fh % read_data(buff(2), '2', record = 4)
      call fh % read_data(buff(3), '3', record = 5)
      call fh % read_data(buff(4), '4', record = 6)
      print *,'buff',a,b,buff
!      call fh % read_data(buff, '1', length = 4, record = 1)
!      print *,'buff',buff
!      call fh % read_data(buff, '3', length = 4, record = 3)
!      print *,'buff',buff
!      call fh % file_close()
      

    end if

  end subroutine test_execute

!===============================================================================
! CHECK
!===============================================================================

  subroutine test_check(suite)

    class(TestSuiteClass), intent(inout) :: suite
    
  end subroutine test_check

!===============================================================================
! TEARDOWN
!===============================================================================

  subroutine test_teardown()

    ! Add teardown code here.  For example, deallocate all memory that was used.
    
  end subroutine test_teardown

end module test_otf_material_file
