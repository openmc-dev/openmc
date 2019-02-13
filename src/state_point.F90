module state_point

!===============================================================================
! STATE_POINT -- This module handles writing and reading state point
! files. State points are contain complete tally results, source sites, and
! various other data. They can be used to restart a run or to reconstruct
! confidence intervals for tallies (this requires post-processing via Python
! scripts).
!
! State points can be written at any batch during a simulation, or at specified
! intervals, using the <state_point ... /> tag.
!===============================================================================

  use, intrinsic :: ISO_C_BINDING

  use hdf5_interface
  use settings
  use string,             only: to_str
  use tally_header

  implicit none

contains

!===============================================================================
! OPENMC_STATEPOINT_WRITE writes an HDF5 statepoint file to disk
!===============================================================================

  subroutine statepoint_write_f(file_id) bind(C)
    integer(HID_T), value :: file_id

    integer :: i
    integer(HID_T) :: tallies_group, tally_group

    ! Open tallies group
    tallies_group = open_group(file_id, "tallies")

    if (reduce_tallies) then
      ! Write global tallies
      call write_dataset(file_id, "global_tallies", global_tallies)

      ! Write tallies
      if (active_tallies_size() > 0) then
        ! Indicate that tallies are on
        call write_attribute(file_id, "tallies_present", 1)

        ! Write all tally results
        TALLY_RESULTS: do i = 1, n_tallies
          associate (tally => tallies(i) % obj)
            ! Write sum and sum_sq for each bin
            tally_group = open_group(tallies_group, "tally " &
                 // to_str(tally % id()))
            call tally % write_results_hdf5(tally_group)
            call close_group(tally_group)
          end associate
        end do TALLY_RESULTS
      else
        ! Indicate tallies are off
        call write_attribute(file_id, "tallies_present", 0)
      end if
    end if

    call close_group(tallies_group)
  end subroutine

!===============================================================================
! LOAD_STATE_POINT
!===============================================================================

  subroutine load_state_point_f(file_id) bind(C)
    integer(HID_T), value :: file_id

    integer :: i
    integer :: temp
    integer(HID_T) :: tallies_group
    integer(HID_T) :: tally_group

    ! Read global tally data
    call read_dataset(global_tallies, file_id, "global_tallies")

    ! Check if tally results are present
    call read_attribute(temp, file_id, "tallies_present")

    ! Read in sum and sum squared
    if (temp == 1) then
      tallies_group = open_group(file_id, "tallies")

      TALLY_RESULTS: do i = 1, n_tallies
        associate (t => tallies(i) % obj)
          ! Read sum, sum_sq, and N for each bin
          tally_group = open_group(tallies_group, "tally " // &
                trim(to_str(t % id())))
          call t % read_results_hdf5(tally_group)
          call read_dataset(t % n_realizations, tally_group, &
                "n_realizations")
          call close_group(tally_group)
        end associate
      end do TALLY_RESULTS

      call close_group(tallies_group)
    end if

  end subroutine load_state_point_f

end module state_point
