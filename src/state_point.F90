module state_point

  use constants
  use global
  use output,             only: write_message
  use string,             only: to_str
  use output_interface,   only: file_create, write_data, file_close, write_string

contains

!===============================================================================
! WRITE_STATE_POINT
!===============================================================================

  subroutine write_state_point()

    character(MAX_FILE_LEN) :: filename
    character(MAX_WORD_LEN) :: fh_state_point = "state_point"
    character(MAX_WORD_LEN) :: fh_source = "source"

    ! Set filename for state point
    filename = trim(path_output) // 'statepoint.' // &
               trim(to_str(current_batch))

    ! Write message
    message = "Creating state point " // trim(filename) // "..."
    call write_message(1)

    ! Create statepoint file 
    call file_create(filename, fh_state_point)

    if (master) then

      ! Write revision number for state point file
      call write_data(REVISION_STATEPOINT, "revision_statepoint")

      ! Write OpenMC version
      call write_data(VERSION_MAJOR, "version_major")
      call write_data(VERSION_MINOR, "version_minor")
      call write_data(VERSION_RELEASE, "version_release")

      ! Write current date and time
      call write_string(time_stamp(), "date_and_time")

    end if

    ! Close statepoint file
    call file_close(fh_state_point)

  end subroutine write_state_point

!===============================================================================
! LOAD_STATE_POINT
!===============================================================================

  subroutine load_state_point()

  end subroutine load_state_point

!===============================================================================
! WRITE_SOURCE
!===============================================================================

  subroutine write_source()

  end subroutine write_source

!===============================================================================
! READ_SOURCE
!===============================================================================

  subroutine read_source()

  end subroutine read_source

!===============================================================================
! REPLAY_BATCH_HISTORY
!===============================================================================

  subroutine replay_batch_history()

  end subroutine replay_batch_history

end module state_point
