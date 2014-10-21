module testing_header 

  use constants

  implicit none
  private
  
  type, public :: TestSuiteClass
    logical :: master
    integer :: n_passed = 0
    integer :: n_skipped = 0
    integer :: n_failed = 0
    contains
      procedure :: fail     => fail
      procedure :: skip     => skip
      procedure :: pass     => pass
  end type TestSuiteClass

  type, abstract, public :: TestClass
    character(MAX_LINE_LEN) :: name
    contains
      procedure(init_interface),             deferred :: init
      procedure(setup_interface),    nopass, deferred :: setup
      procedure(execute_interface),  nopass, deferred :: execute
      procedure(check_interface),    nopass, deferred :: check
      procedure(teardown_interface), nopass, deferred :: teardown
  end type TestClass

  abstract interface
    subroutine init_interface(self)
      import TestClass
      class(TestClass), intent(inout) :: self
    end subroutine init_interface
  end interface

  abstract interface
    subroutine setup_interface(suite)
      import TestSuiteClass
      class(TestSuiteClass), intent(inout) :: suite
    end subroutine setup_interface
  end interface
  
  abstract interface
    subroutine execute_interface()
    end subroutine execute_interface
  end interface
  
  abstract interface
    subroutine check_interface(suite)
      import TestSuiteClass
      class(TestSuiteClass), intent(inout) :: suite
    end subroutine check_interface
  end interface

  abstract interface
    subroutine teardown_interface()
    end subroutine teardown_interface
  end interface

contains

!===============================================================================
! FAIL
!===============================================================================

  subroutine fail(this)

    class(TestSuiteClass) :: this

    this % n_failed = this % n_failed + 1

  end subroutine fail

!===============================================================================
! SKIP
!===============================================================================

  subroutine skip(this)

    class(TestSuiteClass) :: this

    this % n_skipped = this % n_skipped + 1

  end subroutine skip

!===============================================================================
! PASS
!===============================================================================

  subroutine pass(this)

    class(TestSuiteClass) :: this

    this % n_passed = this % n_passed + 1

  end subroutine pass

end module testing_header
