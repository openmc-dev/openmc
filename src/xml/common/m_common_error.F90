module m_common_error

#ifndef DUMMYLIB
  use fox_m_fsys_abort_flush, only: pxfabort, pxfflush
  use fox_m_fsys_array_str, only: vs_str_alloc

  implicit none
  private

  integer, parameter :: ERR_NULL = 0
  integer, parameter :: ERR_WARNING = 1
  integer, parameter :: ERR_ERROR = 2
  integer, parameter :: ERR_FATAL = 3
#endif  
  logical, save :: errors_are_fatal = .false.
  logical, save :: warnings_are_fatal = .false.

#ifndef DUMMYLIB
  type error_t
    integer :: severity = ERR_NULL
    integer :: error_code = 0
    character, dimension(:), pointer :: msg => null()
  end type error_t

  type error_stack
    type(error_t), dimension(:), pointer :: stack => null()
  end type error_stack

  interface FoX_warning
    module procedure FoX_warning_base
  end interface

  interface FoX_error
    module procedure FoX_error_base
  end interface

  interface FoX_fatal
    module procedure FoX_fatal_base
  end interface

  public :: ERR_NULL
  public :: ERR_WARNING
  public :: ERR_ERROR
  public :: ERR_FATAL

  public :: error_t
  public :: error_stack

  public :: init_error_stack
  public :: destroy_error_stack

  public :: FoX_warning
  public :: FoX_error
  public :: FoX_fatal

  public :: FoX_warning_base
  public :: FoX_error_base
  public :: FoX_fatal_base

  public :: add_error
  public :: in_error

#endif
  public :: FoX_set_fatal_errors
  public :: FoX_get_fatal_errors
  public :: FoX_set_fatal_warnings
  public :: FoX_get_fatal_warnings

contains
#ifndef DUMMYLIB

  subroutine FoX_warning_base(msg)
    ! Emit warning, but carry on.
    character(len=*), intent(in) :: msg

    if (warnings_are_fatal) then
        write(0,'(a)') 'FoX warning  made fatal'
        call FoX_fatal_base(msg)
    endif

    write(0,'(a)') 'WARNING(FoX)'
    write(0,'(a)')  msg
    call pxfflush(0)

  end subroutine FoX_warning_base


  subroutine FoX_error_base(msg)
    ! Emit error message and stop.
    ! No clean up is done here, but this can
    ! be overridden to include clean-up routines
    character(len=*), intent(in) :: msg

    if (errors_are_fatal) then
        write(0,'(a)') 'FoX error made fatal'
        call FoX_fatal_base(msg)
    endif

    write(0,'(a)') 'ERROR(FoX)'
    write(0,'(a)')  msg
    call pxfflush(0)

    stop

  end subroutine FoX_error_base

  subroutine FoX_fatal_base(msg)
    !Emit error message and abort with coredump.
    !No clean-up occurs

    character(len=*), intent(in) :: msg

    write(0,'(a)') 'ABORT(FOX)'
    write(0,'(a)')  msg
    call pxfflush(0)

    call pxfabort()

  end subroutine FoX_fatal_base


  subroutine init_error_stack(stack)
    type(error_stack), intent(inout) :: stack
    
    allocate(stack%stack(0))
  end subroutine init_error_stack


  subroutine destroy_error_stack(stack)
    type(error_stack), intent(inout) :: stack
    
    integer :: i

    do i = 1, size(stack%stack)
      deallocate(stack%stack(i)%msg)
    enddo
    deallocate(stack%stack)
  end subroutine destroy_error_stack


  subroutine add_error(stack, msg, severity, error_code)
    type(error_stack), intent(inout) :: stack
    character(len=*), intent(in) :: msg
    integer, intent(in), optional :: severity
    integer, intent(in), optional :: error_code

    integer :: i, n
    type(error_t), dimension(:), pointer :: temp_stack 

    if (.not.associated(stack%stack)) &
      call init_error_stack(stack)

    n = size(stack%stack)
    
    temp_stack => stack%stack
    allocate(stack%stack(n+1))
    do i = 1, size(temp_stack)
      stack%stack(i)%msg => temp_stack(i)%msg
      stack%stack(i)%severity = temp_stack(i)%severity
      stack%stack(i)%error_code = temp_stack(i)%error_code
    enddo
    deallocate(temp_stack)

    stack%stack(n+1)%msg => vs_str_alloc(msg)
    if (present(severity)) then
      stack%stack(n+1)%severity = severity
    else
      stack%stack(n+1)%severity = ERR_ERROR
    endif
    if (present(error_code)) then
      stack%stack(n+1)%error_code = error_code
    else
      stack%stack(n+1)%error_code = -1
    endif

  end subroutine add_error


  function in_error(stack) result(p)
    type(error_stack), intent(in) :: stack
    logical :: p

    if (associated(stack%stack)) then
      p = (size(stack%stack) > 0)
    else
      p = .false.
    endif
  end function in_error

#endif
  subroutine FoX_set_fatal_errors(newvalue)
    logical, intent(in) :: newvalue
    errors_are_fatal = newvalue
  end subroutine FoX_set_fatal_errors

  function  FoX_get_fatal_errors()
     logical :: FoX_get_fatal_errors
     FoX_get_fatal_errors = errors_are_fatal
  end function FoX_get_fatal_errors

  subroutine  FoX_set_fatal_warnings(newvalue)
    logical, intent(in) :: newvalue
    warnings_are_fatal = newvalue
  end subroutine FoX_set_fatal_warnings

  function FoX_get_fatal_warnings()
    logical :: FoX_get_fatal_warnings
    FoX_get_fatal_warnings = warnings_are_fatal
  end function FoX_get_fatal_warnings

end module m_common_error
