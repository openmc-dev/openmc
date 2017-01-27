module URR_error

  implicit none
  private
  public :: EXIT_SUCCESS,&
            EXIT_FAILURE,&
            EXIT_WARNING,&
            ERROR,&
            WARNING,&
            INFO,&
            exit_status,&
            log_message

  ! Error codes returned to application code from calls to the API
  integer, parameter :: EXIT_SUCCESS = 0
  integer, parameter :: EXIT_FAILURE = 1
  integer, parameter :: EXIT_WARNING = 2

  ! logging message codes
  integer, parameter :: DEBUG   = 0
  integer, parameter :: INFO    = 1
  integer, parameter :: WARNING = 2
  integer, parameter :: ERROR   = 3
  integer, parameter :: FATAL   = 4

  integer :: exit_code = EXIT_SUCCESS ! error code returned to application code
                                      ! from calls to the API
!$omp threadprivate(exit_code)


contains


!> Sets API exit code and writes log message
  subroutine exit_status(code, msg)

    integer :: code
    character(*), optional :: msg

    exit_code = code

    select case(exit_code)
    case(EXIT_SUCCESS)
      continue
    case(EXIT_FAILURE)
      call log_message(FATAL, msg)
    case(EXIT_WARNING)
      call log_message(WARNING, msg)
    case default
      continue
    end select

  end subroutine exit_status


!> Writes a log message
  subroutine log_message(level, msg)

    integer, intent(in) :: level
    character(*), intent(in) :: msg

    character(:), allocatable :: level_label

    select case(level)
    case(DEBUG)
      level_label = 'DEBUG: '
    case(INFO)
      level_label = 'INFO: '
    case(WARNING)
      level_label = 'WARNING: '
    case(ERROR)
      level_label = 'ERROR: '
    case(FATAL)
      level_label = 'FATAL: '
    case default
      continue
    end select

    write(*,*) level_label//msg

  end subroutine log_message


end module URR_error
