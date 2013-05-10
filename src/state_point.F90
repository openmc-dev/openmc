module state_point

  use constants
  use global
  use output,             only: write_message, time_stamp
  use string,             only: to_str
  use output_interface

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
      call write_data(time_stamp(), "date_and_time")

      ! Write path to input
      call write_data(path_input, "path")

      ! Write out random number seed
      call write_data(seed, "seed")

      ! Write run information
      call write_data(run_mode, "run_mode")
      call write_data(n_particles, "n_particles")
      call write_data(n_batches, "n_batches")

      ! Write out current batch number
      call write_data(current_batch, "current_batch")

      ! Write out information for eigenvalue run
      if (run_mode == MODE_EIGENVALUE) then
        call write_data(n_inactive, "n_inactive")
        call write_data(gen_per_batch, "gen_per_batch")
        call write_data(k_batch, "k_batch", length=current_batch)
        call write_data(entropy, "entropy", length=current_batch)
        call write_data(k_col_abs, "k_col_abs")
        call write_data(k_col_tra, "k_col_tra")
        call write_data(k_abs_tra, "k_abs_tra")
        call write_data(k_combined, "k_combined", length=2)
      end if

      ! Write number of meshes
      call write_data(n_meshes, "n_meshes", "tallies")

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
