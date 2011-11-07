module logging

  use constants, only: MAX_WORD_LEN, UNIT_LOG
  use global,    only: path_input

  implicit none

contains

!===============================================================================
! CREATE_LOG creates a new log file (or overwrites the existing log)
!===============================================================================

  subroutine create_log()

    character(MAX_WORD_LEN) :: path_log     ! path of log file
    logical                 :: file_exists  ! does log file already exist?
    ! integer                 :: ioError      ! error status for file access

    ! Create filename for log file
    path_log = trim(path_input) // ".log"

    ! Check if log file already exists
    inquire(FILE=path_log, EXIST=file_exists)
    if (file_exists) then
       ! Possibly copy old log file
    end if

    ! Open log file for writing
    ! open(FILE=path_log, UNIT=UNIT_LOG, STATUS='replace', &
    !      & ACTION='write', IOSTAT=ioError)

  end subroutine create_log

!===============================================================================
! LOG_TALLIES
!===============================================================================

  subroutine log_tallies()

  end subroutine log_tallies

end module logging
