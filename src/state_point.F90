module state_point

  use error,  only: warning, fatal_error
  use global
  use output, only: write_message, print_batch_keff
  use string, only: to_str
  use tally_header, only: TallyObject

  implicit none

contains

!===============================================================================
! CREATE_STATE_POINT creates a state point binary file that can be used for
! restarting a run or for getting intermediate tally results
!===============================================================================

  subroutine create_state_point()

    integer :: i, j, k ! loop indices
    type(TallyObject), pointer :: t => null()

    ! Set filename for binary state point
    path_state_point = 'restart.' // trim(to_str(current_batch)) // '.binary'

    ! Write message
    message = "Creating state point " // trim(path_state_point) // "..."
    call write_message()

    ! Open binary state point file for writing
    open(UNIT=UNIT_STATE, FILE=path_state_point, STATUS='replace', &
         ACCESS='stream')

    ! Write revision number for state point file
    write(UNIT_STATE) REVISION_STATEPOINT
    
    ! Write OpenMC version
    write(UNIT_STATE) VERSION_MAJOR, VERSION_MINOR, VERSION_RELEASE

    ! Write out random number seed
    write(UNIT_STATE) seed

    ! Write run information
    write(UNIT_STATE) run_mode, n_particles, n_batches, &
         n_inactive, gen_per_batch

    ! Write out current batch number
    write(UNIT_STATE) current_batch

    ! Write out keff and entropy for each batch
    if (run_mode == MODE_CRITICALITY) then
       write(UNIT_STATE) k_batch(1:current_batch)
       if (entropy_on) write(UNIT_STATE) entropy(1:current_batch)
    end if

    ! Write out global tallies sum and sum_sq
    write(UNIT_STATE) N_GLOBAL_TALLIES
    write(UNIT_STATE) global_tallies(:) % sum
    write(UNIT_STATE) global_tallies(:) % sum_sq

    ! Write number of meshes
    write(UNIT_STATE) n_meshes

    ! Write information for meshes
    do i = 1, n_meshes
       write(UNIT_STATE) meshes(i) % type
       write(UNIT_STATE) meshes(i) % n_dimension
       write(UNIT_STATE) meshes(i) % dimension
       write(UNIT_STATE) meshes(i) % lower_left
       write(UNIT_STATE) meshes(i) % upper_right
       write(UNIT_STATE) meshes(i) % width
    end do

    ! Write number of tallies
    write(UNIT_STATE) n_tallies

    ! Write size of each tally
    do i = 1, n_tallies
       ! Get pointer to tally
       t => tallies(i)

       write(UNIT_STATE) size(t % scores, 1)
       write(UNIT_STATE) size(t % scores, 2)

       ! Write number of filters
       write(UNIT_STATE) t % n_filters
       do j = 1, t % n_filters
          write(UNIT_STATE) t % filters(j)
          write(UNIT_STATE) t % n_filter_bins(t % filters(j))

          select case (t % filters(j))
          case(FILTER_UNIVERSE)
             write(UNIT_STATE) t % universe_bins(:) % scalar
          case(FILTER_MATERIAL)
             write(UNIT_STATE) t % material_bins(:) % scalar
          case(FILTER_CELL)
             write(UNIT_STATE) t % cell_bins(:) % scalar
          case(FILTER_CELLBORN)
             write(UNIT_STATE) t % cellborn_bins(:) % scalar
          case(FILTER_SURFACE)
             write(UNIT_STATE) t % surface_bins(:) % scalar
          case(FILTER_MESH)
             write(UNIT_STATE) t % mesh
          case(FILTER_ENERGYIN)
             write(UNIT_STATE) t % energy_in
          case(FILTER_ENERGYOUT)
             write(UNIT_STATE) t % energy_out
          end select
       end do

       ! Write number of nuclide bins
       write(UNIT_STATE) t % n_nuclide_bins

       ! Write nuclide bins
       do j = 1, t % n_nuclide_bins
          if (t % nuclide_bins(j) % scalar > 0) then
             write(UNIT_STATE) nuclides(t % nuclide_bins(j) % scalar) % zaid
          else
             write(UNIT_STATE) t % nuclide_bins(j) % scalar
          end if
       end do

       ! Write number of score bins
       write(UNIT_STATE) t % n_score_bins
       write(UNIT_STATE) t % score_bins(:) % scalar
    end do

    if (tallies_on) then
       ! Write tally sum and sum_sq
       do i = 1, n_tallies
          do k = 1, size(t % scores, 2)
             do j = 1, size(t % scores, 1)
                write(UNIT_STATE) t % scores(j,k) % sum
                write(UNIT_STATE) t % scores(j,k) % sum_sq
             end do
          end do
       end do
    end if

    ! Close binary state point file
    close(UNIT_STATE)

  end subroutine create_state_point

!===============================================================================
! LOAD_STATE_POINT loads data from a state point file to either continue a run
! or to print intermediate tally results
!===============================================================================

  subroutine load_state_point()

    integer :: i, j, k ! loop indices
    integer :: mode    ! specified run mode
    integer :: temp(3) ! temporary variable
    integer, allocatable :: int_array(:)
    real(8), allocatable :: real_array(:)

    ! Write message
    message = "Loading state point " // trim(path_state_point) // "..."
    call write_message(1)

    ! Open binary state point file for writing
    open(UNIT=UNIT_STATE, FILE=path_state_point, STATUS='old', &
         ACCESS='stream')

    ! Raad revision number for state point file and make sure it matches with
    ! current version
    read(UNIT_STATE) temp(1)
    if (temp(1) /= REVISION_STATEPOINT) then
       message = "State point binary version does not match current version " &
            // "in OpenMC."
       call fatal_error()
    end if
    
    ! Read OpenMC version
    read(UNIT_STATE) temp(1:3)
    if (temp(1) /= VERSION_MAJOR .or. temp(2) /= VERSION_MINOR &
         .or. temp(3) /= VERSION_RELEASE) then
       message = "State point file was created with a different version " // &
            "of OpenMC."
       call warning()
    end if

    ! Read and overwrite random number seed
    read(UNIT_STATE) seed

    ! Read and overwrite run information
    read(UNIT_STATE) mode, n_particles, n_batches, &
         n_inactive, gen_per_batch

    ! Read batch number to restart at
    read(UNIT_STATE) restart_batch

    ! Read keff and entropy for each batch
    if (mode == MODE_CRITICALITY) then
       read(UNIT_STATE) k_batch(1:restart_batch)
       if (entropy_on) read(UNIT_STATE) entropy(1:restart_batch)
    end if

    if (master) then
       ! Read number of global tallies and make sure it matches
       read(UNIT_STATE) temp(1)
       if (temp(1) /= N_GLOBAL_TALLIES) then
          message = "Number of global tallies does not match in state point."
          call fatal_error()
       end if

       ! Read global tally data
       read(UNIT_STATE) global_tallies(:) % sum
       read(UNIT_STATE) global_tallies(:) % sum_sq

       ! Read number of meshes
       read(UNIT_STATE) temp(1)
       if (temp(1) /= n_meshes) then
          message = "Number of meshes does not match in state point."
          call fatal_error()
       end if

       MESH_LOOP: do i = 1, n_meshes
          ! Read type of mesh and dimension
          read(UNIT_STATE) temp(1:2)

          ! Allocate temporary arrays
          allocate(int_array(temp(2)))
          allocate(real_array(temp(2)))

          ! Read dimension, lower_left, upper_right, width
          read(UNIT_STATE) int_array
          read(UNIT_STATE) real_array
          read(UNIT_STATE) real_array
          read(UNIT_STATE) real_array

          ! Deallocate temporary arrays
          deallocate(int_array)
          deallocate(real_array)
       end do MESH_LOOP

       ! Read number of tallies and make sure it matches
       read(UNIT_STATE) temp(1)
       if (temp(1) /= n_tallies) then
          message = "Number of tallies does not match in state point."
          call fatal_error()
       end if

       TALLY_METADATA: do i = 1, n_tallies
          ! Read dimensions of tally filters and scores and make sure they
          ! match
          read(UNIT_STATE) temp(1:2)
          if (temp(1) /= size(tallies(i) % scores, 1) .or. &
               temp(2) /= size(tallies(i) % scores, 2)) then
             message = "Tally dimensions do not match in state point."
             call fatal_error()
          end if

          ! Read number of filters
          read(UNIT_STATE) temp(1)

          FILTER_LOOP: do j = 1, temp(1)
             ! Read filter type and number of bins
             read(UNIT_STATE) temp(2:3)

             ! Read filter bins
             select case (temp(2))
             case (FILTER_MESH)
                allocate(int_array(1))
                read(UNIT_STATE) int_array
                deallocate(int_array)
             case (FILTER_ENERGYIN, FILTER_ENERGYOUT)
                allocate(real_array(temp(3) + 1))
                read(UNIT_STATE) real_array
                deallocate(real_array)
             case default
                allocate(int_array(temp(3)))
                read(UNIT_STATE) int_array
             end select
          end do FILTER_LOOP

          ! Read number of nuclides
          read(UNIT_STATE) temp(1)

          ! Read nuclide bins
          allocate(int_array(temp(1)))
          read(UNIT_STATE) int_array
          deallocate(int_array)

          ! Read number of scores
          read(UNIT_STATE) temp(1)

          ! Read nuclide bins
          allocate(int_array(temp(1)))
          read(UNIT_STATE) int_array
          deallocate(int_array)
       end do TALLY_METADATA

          ! Read sum and sum squared
       if (restart_batch > n_inactive) then
          TALLY_SCORES: do i = 1, n_tallies
             do k = 1, size(tallies(i) % scores, 2)
                do j = 1, size(tallies(i) % scores, 1)
                   read(UNIT_STATE) tallies(i) % scores(j,k) % sum
                   read(UNIT_STATE) tallies(i) % scores(j,k) % sum_sq
                end do
             end do
          end do TALLY_SCORES
       end if
    end if

    ! Close binary state point file
    close(UNIT_STATE)

  end subroutine load_state_point

!===============================================================================
! REPLAY_BATCH_HISTORY displays batch keff and entropy for each batch stored in
! a state point file
!===============================================================================

  subroutine replay_batch_history

    real(8), save :: temp(2) = ZERO

    ! Write message at beginning
    if (current_batch == 1) then
       message = "Replaying history from state point..."
       call write_message(1)
    end if

    ! For criticality calculations, turn on tallies if we've reached active
    ! batches
    if (current_batch == n_inactive) tallies_on = .true.

    ! Add to number of realizations
    if (current_batch > n_inactive) then
       n_realizations = n_realizations + 1

       temp(1) = temp(1) + k_batch(current_batch)
       temp(2) = temp(2) + k_batch(current_batch)*k_batch(current_batch)

       keff = temp(1) / n_realizations
       keff_std = sqrt((temp(2)/n_realizations - keff*keff) &
            / (n_realizations - 1))
    else
       keff = k_batch(current_batch)
    end if

    ! print out batch keff
    call print_batch_keff()

    ! Write message at end
    if (current_batch == restart_batch) then
       message = "Resuming simulation..."
       call write_message(1)
    end if

  end subroutine replay_batch_history

end module state_point
