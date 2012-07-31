module state_point

  use global
  use source, only: write_source_binary
  use string, only: to_str

  implicit none

contains

!===============================================================================
! CREATE_STATE_POINT creates a state point binary file that can be used for
! restarting a run
!===============================================================================

  subroutine create_state_point()

    integer :: i ! loop index

    ! Set filename for binary state point
    path_state_point = 'restart.' // trim(to_str(current_batch)) // '.binary'

    ! Open binary state point file for writing
    open(UNIT=UNIT_STATE, FILE=path_state_point, STATUS='replace', &
         ACCESS='stream')

    ! Write revision number for state point file
    write(UNIT_STATE) REVISION_STATEPOINT
    
    ! Write OpenMC version
    write(UNIT_STATE) VERSION_MAJOR, VERSION_MINOR, VERSION_RELEASE

    ! Write run information
    write(UNIT_STATE) run_mode, n_particles, n_batches, &
         n_inactive, gen_per_batch

    ! Write out current batch number
    write(UNIT_STATE) current_batch

    ! Write out global tallies sum and sum_sq
    write(UNIT_STATE) N_GLOBAL_TALLIES
    write(UNIT_STATE) global_tallies(:) % sum
    write(UNIT_STATE) global_tallies(:) % sum_sq

    ! Write out tallies sum and sum_sq
    write(UNIT_STATE) n_tallies
    do i = 1, n_tallies
       write(UNIT_STATE) size(tallies(i) % scores, 1)
       write(UNIT_STATE) size(tallies(i) % scores, 2)
       write(UNIT_STATE) tallies(i) % scores(:,:) % sum
       write(UNIT_STATE) tallies(i) % scores(:,:) % sum_sq
    end do

    ! Close binary state point file
    close(UNIT_STATE)

    ! For a criticality simulation, write the source file
    if (run_mode == MODE_CRITICALITY) call write_source_binary()

  end subroutine create_state_point

end module state_point
