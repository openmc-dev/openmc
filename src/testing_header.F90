module testing_header 

  use output_header, only: output_message

  implicit none
  private
  
  type, public :: testing_type
    integer :: n_passed = 0
    integer :: n_skipped = 0
    integer :: n_failed = 0
    contains
      procedure :: setup => setup
      procedure :: execute => execute
      procedure :: check => check
      procedure :: fail => fail
      procedure :: skip => skip
      procedure :: pass => pass
  end type testing_type

contains

!===============================================================================
! SETUP
!===============================================================================

  subroutine setup(this)
    
    class(testing_type) :: this
  
    call output_message("Setting up...")

  end subroutine setup

!===============================================================================
! EXECUTE
!===============================================================================

  subroutine execute(this)

    class(testing_type) :: this

    call output_message("Invoking test...")

  end subroutine execute

!===============================================================================
! CHECK
!===============================================================================

  subroutine check(this)

    class(testing_type) :: this

    call output_message("Checking results...")

  end subroutine check

!===============================================================================
! FAIL
!===============================================================================

  subroutine fail(this, msg)

    class(testing_type) :: this
    character(len=*) :: msg

    call output_message("FAILED")
    call output_message(msg)
    
    this % n_failed = this % n_failed + 1

  end subroutine fail

!===============================================================================
! SKIP
!===============================================================================

  subroutine skip(this)

    class(testing_type) :: this

    call output_message("SKIPPED")
    
    this % n_skipped = this % n_skipped + 1

  end subroutine skip

!===============================================================================
! PASS
!===============================================================================

  subroutine pass(this)

    class(testing_type) :: this

    call output_message("PASS")
    
    this % n_passed = this % n_passed + 1

  end subroutine pass

end module testing_header
