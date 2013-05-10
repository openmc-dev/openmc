module state_point

  use constants
  use global
  use output,             only: write_message, time_stamp
  use string,             only: to_str
  use output_interface
  use tally_header,       only: TallyObject

  implicit none

contains

!===============================================================================
! WRITE_STATE_POINT
!===============================================================================

  subroutine write_state_point()

    character(MAX_FILE_LEN) :: filename
    character(MAX_WORD_LEN) :: fh_state_point = "state_point"
    character(MAX_WORD_LEN) :: fh_source = "source"
    integer                 :: i
    integer                 :: j
    integer, allocatable    :: temp_array(:)
    type(TallyObject), pointer :: t => null()

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
      call write_data(n_meshes, "n_meshes", group="tallies")

      ! Write information for meshes
      MESH_LOOP: do i = 1, n_meshes
        call write_data(meshes(i) % id, "id", &
             group="tallies/mesh" // to_str(i))
        call write_data(meshes(i) % type, "type", &
             group="tallies/mesh" // to_str(i))
        call write_data(meshes(i) % n_dimension, "n_dimension", &
             group="tallies/mesh" // to_str(i))
        call write_data(meshes(i) % dimension, "dimension", &
             group="tallies/mesh" // to_str(i), &
             length=meshes(i) % n_dimension)
        call write_data(meshes(i) % lower_left, "lower_left", &
             group="tallies/mesh" // to_str(i), &
             length=meshes(i) % n_dimension)
        call write_data(meshes(i) % upper_right, "upper_right", &
             group="tallies/mesh" // to_str(i), &
             length=meshes(i) % n_dimension)
        call write_data(meshes(i) % width, "width", &
             group="tallies/mesh" // to_str(i), &
             length=meshes(i) % n_dimension)
      end do MESH_LOOP

      ! Write number of tallies
      call write_data(n_tallies, "n_tallies", group="tallies")

      ! Write all tally information except results
      TALLY_METADATA: do i = 1, n_tallies
        !Get pointer to tally
        t => tallies(i)

        ! Write id
        call write_data(t % id, "id", group="tallies/tally" // to_str(i))

        ! Write number of realizations
        call write_data(t % n_realizations, "n_realizations", &
             group="tallies/tally" // to_str(i))

        ! Write size of each tally
        call write_data(t % total_score_bins, "total_score_bins", &
             group="tallies/tally" // to_str(i))
        call write_data(t % total_filter_bins, "total_filter_bins", &
             group="tallies/tally" // to_str(i))

        ! Write number of filters
        call write_data(t % n_filters, "n_filters", &
             group="tallies/tally" // to_str(i))

        ! Write filter information
        FILTER_LOOP: do j = 1, t % n_filters

          ! Write type of filter
          call write_data(t % filters(j) % type, "type", &
               group="tallies/tally" // trim(to_str(i)) // "/filter" // to_str(j))

          ! Write number of bins for this filter
          call write_data(t % filters(j) % n_bins, "n_bins", &
               group="tallies/tally" // trim(to_str(i)) // "/filter" // to_str(j))

          ! Write bins
          if (t % filters(j) % type == FILTER_ENERGYIN .or. &
              t % filters(j) % type == FILTER_ENERGYOUT) then
            call write_data(t % filters(j) % real_bins, "bins", &
                 group="tallies/tally" // trim(to_str(i)) // "/filter" // to_str(j), &
                 length=size(t % filters(j) % real_bins))
          else
            call write_data(t % filters(j) % int_bins, "bins", &
                 group="tallies/tally" // trim(to_str(i)) // "/filter" // to_str(j), &
                 length=size(t % filters(j) % int_bins))
          end if

        end do FILTER_LOOP

        ! Write number of nuclide bins
        call write_data(t % n_nuclide_bins, "n_nuclide_bins", &
             group="tallies/tally" // to_str(i))

        ! Set up nuclide bin array and then write
        allocate(temp_array(t % n_nuclide_bins))
        NUCLIDE_LOOP: do j = 1, t % n_nuclide_bins
          if (t % nuclide_bins(j) > 0) then
            temp_array(j) = nuclides(t % nuclide_bins(j)) % zaid
          else
            temp_array(j) = t % nuclide_bins(j)
          end if
        end do NUCLIDE_LOOP
        call write_data(temp_array, "nuclide_bins", &
             group="tallies/tally" // to_str(i), length=t % n_nuclide_bins)
        deallocate(temp_array)

        ! Write number of score bins, score bins, and scatt order
        call write_data(t % n_score_bins, "n_score_bins", &
             group="tallies/tally" // to_str(i))
        call write_data(t % score_bins, "score_bins", &
             group="tallies/tally" // to_str(i), length=t % n_score_bins)
        call write_data(t % scatt_order, "scatt_order", &
             group="tallies/tally" // to_str(i), length=t % n_score_bins)

        ! Write number of user socre bins
        call write_data(t % n_user_score_bins, "n_user_score_bins", &
             group="tallies/tally" // to_str(i))

      end do TALLY_METADATA

    end if

    ! Check for the no-tally-reduction method
    if (.not. reduce_tallies) then

       ! TODO Call the no tally reduce routine

    elseif (master) then

      ! Write number of global realizations
      call write_data(n_realizations, "n_realizations")

      ! Write global tallies
      call write_data(N_GLOBAL_TALLIES, "n_global_tallies")
      call write_tally_result(global_tallies, "global_tallies", &
           length=N_GLOBAL_TALLIES)

      ! Write tallies
      if (tallies_on) then

        ! Indicate that tallies are on
        call write_data(1, "tallies_present")

        ! Write all tally results
        TALLY_RESULTS: do i = 1, n_tallies

          ! Set point to current tally
          t => tallies(i)

          ! Write sum and sum_sq for each bin
          call write_tally_result(t % results, "results", &
               group="tallies/tally" // to_str(i), &
               length=size(t % results, 1) * size(t % results, 2))

        end do TALLY_RESULTS

      else

        ! Indicate tallies are off
        call write_data(0, "tallies_present")

      end if

    end if

    ! Change file handles if source is separately written
#ifdef HDF5
# ifdef MPI
    source_separate = .true.
# endif
#endif

    if (source_separate) then

      ! Close statepoint file 
      call file_close(fh_state_point)

      ! Set filename for source
      filename = trim(path_output) // 'source.' // &
                 trim(to_str(current_batch))

      ! Write message
      message = "Creating source file " // trim(filename) // "..."
      call write_message(1)

      ! Create statepoint file 
      call file_create(filename, fh_source)

    end if

    ! Write out source
    call write_source_bank()

    ! Close file
    if (source_separate) then
      call file_close(fh_source)
    else
      call file_close(fh_state_point)
    end if

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
