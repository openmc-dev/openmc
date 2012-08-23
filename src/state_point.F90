module state_point

  use error,        only: warning, fatal_error
  use global
  use math,         only: t_percentile
  use output,       only: write_message, print_batch_keff
  use string,       only: to_str
  use tally_header, only: TallyObject

#ifdef MPI
  use mpi
#endif

  implicit none

contains

!===============================================================================
! CREATE_STATE_POINT creates a state point binary file that can be used for
! restarting a run or for getting intermediate tally results
!===============================================================================

  subroutine create_state_point()

    integer :: i ! loo pindex
    integer :: n ! temporary array length
    type(TallyObject), pointer :: t => null()

#ifdef MPI
    integer :: fh                      ! file handle
    integer :: temp
    integer :: size_offset_kind        ! size of MPI_OFFSET_KIND (bytes)
    integer :: size_bank               ! size of MPI_BANK type
    integer(MPI_OFFSET_KIND) :: offset ! offset in memory (0=beginning of file)
#else
    integer :: j, k ! loop indices
#endif

    ! Set filename for binary state point
    path_state_point = 'statepoint.' // trim(to_str(current_batch)) // '.binary'

    ! Write message
    message = "Creating state point " // trim(path_state_point) // "..."
    call write_message()

#ifdef MPI
    ! ==========================================================================
    ! PARALLEL I/O USING MPI-2 ROUTINES

    ! Open binary source file for reading
    call MPI_FILE_OPEN(MPI_COMM_WORLD, path_state_point, MPI_MODE_CREATE + &
         MPI_MODE_WRONLY, MPI_INFO_NULL, fh, mpi_err)

    if (master) then
       ! =======================================================================
       ! RUN INFORMATION AND TALLY METADATA

       call state_point_header(fh)

       ! =======================================================================
       ! TALLY RESULTS

       if (tallies_on) then
          ! Indicate that tallies are on
          temp = 1
          call MPI_FILE_WRITE(fh, temp, 1, MPI_INTEGER, &
               MPI_STATUS_IGNORE, mpi_err)

          ! Write all tally scores
          TALLY_SCORES: do i = 1, n_tallies
             t => tallies(i)
             
             n = size(t % scores, 1) * size(t % scores, 2)
             call MPI_FILE_WRITE(fh, t % scores, n, MPI_TALLYSCORE, &
                  MPI_STATUS_IGNORE, mpi_err)
          end do TALLY_SCORES
       else
          ! Indicate that tallies are off
          temp = 0
          call MPI_FILE_WRITE(fh, temp, 1, MPI_INTEGER, &
               MPI_STATUS_IGNORE, mpi_err)
       end if
    end if

    ! =======================================================================
    ! SOURCE BANK

    if (run_mode == MODE_CRITICALITY) then
       ! Get current offset for master
       if (master) call MPI_FILE_GET_POSITION(fh, offset, mpi_err)

       ! Determine offset on master process and broadcast to all processors
       call MPI_SIZEOF(offset, size_offset_kind, mpi_err)
       select case (size_offset_kind)
       case (4)
          call MPI_BCAST(offset, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)
       case (8)
          call MPI_BCAST(offset, 1, MPI_INTEGER8, 0, MPI_COMM_WORLD, mpi_err)
       end select

       ! Set proper offset for source data on this processor
       call MPI_TYPE_SIZE(MPI_BANK, size_bank, mpi_err)
       offset = offset + size_bank*maxwork*rank

       ! Write all source sites
       call MPI_FILE_WRITE_AT(fh, offset, source_bank(1), work, MPI_BANK, &
            MPI_STATUS_IGNORE, mpi_err)
    end if
       
    ! Close binary source file
    call MPI_FILE_CLOSE(fh, mpi_err)

#else
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
    write(UNIT_STATE) run_mode, n_particles, n_batches

    ! Write out current batch number
    write(UNIT_STATE) current_batch

    ! Write out information for criticality run
    if (run_mode == MODE_CRITICALITY) then
       write(UNIT_STATE) n_inactive, gen_per_batch
       write(UNIT_STATE) k_batch(1:current_batch)
       write(UNIT_STATE) entropy(1:current_batch)
    end if

    ! Write out global tallies sum and sum_sq
    write(UNIT_STATE) N_GLOBAL_TALLIES
    GLOBAL_TALLIES_LOOP: do i = 1, N_GLOBAL_TALLIES
       write(UNIT_STATE) global_tallies(i) % sum
       write(UNIT_STATE) global_tallies(i) % sum_sq
    end do GLOBAL_TALLIES_LOOP

    ! Write number of meshes
    write(UNIT_STATE) n_meshes

    ! Write information for meshes
    MESH_LOOP: do i = 1, n_meshes
       write(UNIT_STATE) meshes(i) % type
       write(UNIT_STATE) meshes(i) % n_dimension
       write(UNIT_STATE) meshes(i) % dimension
       write(UNIT_STATE) meshes(i) % lower_left
       write(UNIT_STATE) meshes(i) % upper_right
       write(UNIT_STATE) meshes(i) % width
    end do MESH_LOOP

    ! Write number of tallies
    write(UNIT_STATE) n_tallies

    TALLY_METADATA: do i = 1, n_tallies
       ! Get pointer to tally
       t => tallies(i)

       ! Write size of each tally
       write(UNIT_STATE) size(t % scores, 1)
       write(UNIT_STATE) size(t % scores, 2)

       ! Write number of filters
       write(UNIT_STATE) t % n_filters

       FILTER_LOOP: do j = 1, t % n_filters
          ! Write type of filter
          write(UNIT_STATE) t % filters(j)

          ! Write number of bins for this filter
          write(UNIT_STATE) t % n_filter_bins(t % filters(j))

          ! Write filter bins
          select case (t % filters(j))
          case(FILTER_UNIVERSE)
             write(UNIT_STATE) t % universe_bins
          case(FILTER_MATERIAL)
             write(UNIT_STATE) t % material_bins
          case(FILTER_CELL)
             write(UNIT_STATE) t % cell_bins
          case(FILTER_CELLBORN)
             write(UNIT_STATE) t % cellborn_bins
          case(FILTER_SURFACE)
             write(UNIT_STATE) t % surface_bins
          case(FILTER_MESH)
             write(UNIT_STATE) t % mesh
          case(FILTER_ENERGYIN)
             write(UNIT_STATE) t % energy_in
          case(FILTER_ENERGYOUT)
             write(UNIT_STATE) t % energy_out
          end select
       end do FILTER_LOOP

       ! Write number of nuclide bins
       write(UNIT_STATE) t % n_nuclide_bins

       ! Write nuclide bins
       NUCLIDE_LOOP: do j = 1, t % n_nuclide_bins
          if (t % nuclide_bins(j) > 0) then
             write(UNIT_STATE) nuclides(t % nuclide_bins(j)) % zaid
          else
             write(UNIT_STATE) t % nuclide_bins(j)
          end if
       end do NUCLIDE_LOOP

       ! Write number of score bins
       write(UNIT_STATE) t % n_score_bins
       write(UNIT_STATE) t % score_bins
    end do TALLY_METADATA

    if (tallies_on) then
       ! Indicate that tallies are on
       write(UNIT_STATE) 1

       ! Write tally sum and sum_sq
       TALLY_SCORES: do i = 1, n_tallies
          do k = 1, size(t % scores, 2)
             do j = 1, size(t % scores, 1)
                write(UNIT_STATE) t % scores(j,k) % sum
                write(UNIT_STATE) t % scores(j,k) % sum_sq
             end do
          end do
       end do TALLY_SCORES
    else
       ! Indicate that tallies are off
       write(UNIT_STATE) 0
    end if

    ! Write out source bank 
    if (run_mode == MODE_CRITICALITY) then
       write(UNIT_STATE) source_bank
    end if

    ! Close binary state point file
    close(UNIT_STATE)
#endif

  end subroutine create_state_point

#ifdef MPI
!===============================================================================
! STATE_POINT_HEADER
!===============================================================================

  subroutine state_point_header(fh)

    integer, intent(inout) :: fh ! file handle

    integer :: i  ! loop index
    integer :: j  ! loop index
    integer :: n  ! temporary array length
    type(TallyObject), pointer :: t => null()

    ! Write revision number for state point file
    call MPI_FILE_WRITE(fh, REVISION_STATEPOINT, 1, MPI_INTEGER, &
         MPI_STATUS_IGNORE, mpi_err)
    
    ! Write OpenMC version
    call MPI_FILE_WRITE(fh, VERSION_MAJOR, 1, MPI_INTEGER, &
         MPI_STATUS_IGNORE, mpi_err)
    call MPI_FILE_WRITE(fh, VERSION_MINOR, 1, MPI_INTEGER, &
         MPI_STATUS_IGNORE, mpi_err)
    call MPI_FILE_WRITE(fh, VERSION_RELEASE, 1, MPI_INTEGER, &
         MPI_STATUS_IGNORE, mpi_err)

    ! Write out random number seed
    call MPI_FILE_WRITE(fh, seed, 1, MPI_INTEGER8, &
         MPI_STATUS_IGNORE, mpi_err)

    ! Write run information
    call MPI_FILE_WRITE(fh, run_mode, 1, MPI_INTEGER, &
         MPI_STATUS_IGNORE, mpi_err)
    call MPI_FILE_WRITE(fh, n_particles, 1, MPI_INTEGER8, &
         MPI_STATUS_IGNORE, mpi_err)
    call MPI_FILE_WRITE(fh, n_batches, 1, MPI_INTEGER, &
         MPI_STATUS_IGNORE, mpi_err)

    ! Write out current batch number
    call MPI_FILE_WRITE(fh, current_batch, 1, MPI_INTEGER, &
         MPI_STATUS_IGNORE, mpi_err)

    ! Write out information for criticality run
    if (run_mode == MODE_CRITICALITY) then
       call MPI_FILE_WRITE(fh, n_inactive, 1, MPI_INTEGER, &
            MPI_STATUS_IGNORE, mpi_err)
       call MPI_FILE_WRITE(fh, gen_per_batch, 1, MPI_INTEGER, &
            MPI_STATUS_IGNORE, mpi_err)
       call MPI_FILE_WRITE(fh, k_batch, current_batch, MPI_REAL8, &
            MPI_STATUS_IGNORE, mpi_err)
       call MPI_FILE_WRITE(fh, entropy, current_batch, MPI_REAL8, &
            MPI_STATUS_IGNORE, mpi_err)
    end if

    ! Write out global tallies sum and sum_sq
    call MPI_FILE_WRITE(fh, N_GLOBAL_TALLIES, 1, MPI_INTEGER, &
         MPI_STATUS_IGNORE, mpi_err)
    call MPI_FILE_WRITE(fh, global_tallies, N_GLOBAL_TALLIES, &
         MPI_TALLYSCORE, MPI_STATUS_IGNORE, mpi_err)

    ! Write number of meshes
    call MPI_FILE_WRITE(fh, n_meshes, 1, MPI_INTEGER, &
         MPI_STATUS_IGNORE, mpi_err)

    ! Write information for meshes
    MESH_LOOP: do i = 1, n_meshes
       call MPI_FILE_WRITE(fh, meshes(i) % type, 1, MPI_INTEGER, &
            MPI_STATUS_IGNORE, mpi_err)
       call MPI_FILE_WRITE(fh, meshes(i) % n_dimension, 1, MPI_INTEGER, &
            MPI_STATUS_IGNORE, mpi_err)

       n = meshes(i) % n_dimension
       call MPI_FILE_WRITE(fh, meshes(i) % dimension, n, MPI_INTEGER, &
            MPI_STATUS_IGNORE, mpi_err)
       call MPI_FILE_WRITE(fh, meshes(i) % lower_left, n, MPI_REAL8, &
            MPI_STATUS_IGNORE, mpi_err)
       call MPI_FILE_WRITE(fh, meshes(i) % upper_right, n, MPI_REAL8, &
            MPI_STATUS_IGNORE, mpi_err)
       call MPI_FILE_WRITE(fh, meshes(i) % width, n, MPI_REAL8, &
            MPI_STATUS_IGNORE, mpi_err)
    end do MESH_LOOP

    ! Write number of tallies
    call MPI_FILE_WRITE(fh, n_tallies, 1, MPI_INTEGER, &
         MPI_STATUS_IGNORE, mpi_err)

    TALLY_METADATA: do i = 1, n_tallies
       ! Get pointer to tally
       t => tallies(i)

       ! Write size of each tally
       n = t % n_score_bins * t % n_nuclide_bins
       call MPI_FILE_WRITE(fh, n, 1, MPI_INTEGER, MPI_STATUS_IGNORE, mpi_err)
       call MPI_FILE_WRITE(fh, t % n_total_bins, 1, MPI_INTEGER, &
            MPI_STATUS_IGNORE, mpi_err)

       ! Write number of filters
       call MPI_FILE_WRITE(fh, t % n_filters, 1, MPI_INTEGER, &
            MPI_STATUS_IGNORE, mpi_err)

       FILTER_LOOP: do j = 1, t % n_filters
          ! Write type of filter
          call MPI_FILE_WRITE(fh, t % filters(j), 1, MPI_INTEGER, &
               MPI_STATUS_IGNORE, mpi_err)

          ! Write number of bins for this filter
          n = t % n_filter_bins(t % filters(j))
          call MPI_FILE_WRITE(fh, t % n_filter_bins(t % filters(j)), &
               1, MPI_INTEGER, MPI_STATUS_IGNORE, mpi_err)

          select case (t % filters(j))
          case(FILTER_UNIVERSE)
             call MPI_FILE_WRITE(fh, t % universe_bins, n, &
                  MPI_INTEGER, MPI_STATUS_IGNORE, mpi_err)
          case(FILTER_MATERIAL)
             call MPI_FILE_WRITE(fh, t % material_bins, n, &
                  MPI_INTEGER, MPI_STATUS_IGNORE, mpi_err)
          case(FILTER_CELL)
             call MPI_FILE_WRITE(fh, t % cell_bins, n, &
                  MPI_INTEGER, MPI_STATUS_IGNORE, mpi_err)
          case(FILTER_CELLBORN)
             call MPI_FILE_WRITE(fh, t % cellborn_bins, n, &
                  MPI_INTEGER, MPI_STATUS_IGNORE, mpi_err)
          case(FILTER_SURFACE)
             call MPI_FILE_WRITE(fh, t % surface_bins, n, &
                  MPI_INTEGER, MPI_STATUS_IGNORE, mpi_err)
          case(FILTER_MESH)
             call MPI_FILE_WRITE(fh, t % mesh, 1, &
                  MPI_INTEGER, MPI_STATUS_IGNORE, mpi_err)
          case(FILTER_ENERGYIN)
             call MPI_FILE_WRITE(fh, t % energy_in, n+1, &
                  MPI_REAL8, MPI_STATUS_IGNORE, mpi_err)
          case(FILTER_ENERGYOUT)
             call MPI_FILE_WRITE(fh, t % energy_out, n+1, &
                  MPI_REAL8, MPI_STATUS_IGNORE, mpi_err)
          end select
       end do FILTER_LOOP

       ! Write number of nuclide bins
       call MPI_FILE_WRITE(fh, t % n_nuclide_bins, 1, MPI_INTEGER, &
            MPI_STATUS_IGNORE, mpi_err)

       ! Write nuclide bins
       NUCLIDE_LOOP: do j = 1, t % n_nuclide_bins
          if (t % nuclide_bins(j) > 0) then
             call MPI_FILE_WRITE(fh, nuclides(t % nuclide_bins(j)) % zaid, &
                  1, MPI_INTEGER, MPI_STATUS_IGNORE, mpi_err)
          else
             call MPI_FILE_WRITE(fh, t % nuclide_bins(j), 1, &
                  MPI_INTEGER, MPI_STATUS_IGNORE, mpi_err)
          end if
       end do NUCLIDE_LOOP

       ! Write number of score bins
       call MPI_FILE_WRITE(fh, t % n_score_bins, 1, MPI_INTEGER, &
            MPI_STATUS_IGNORE, mpi_err)
       call MPI_FILE_WRITE(fh, t % score_bins, t % n_score_bins, &
            MPI_INTEGER, MPI_STATUS_IGNORE, mpi_err)
    end do TALLY_METADATA

  end subroutine state_point_header
#endif

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
    read(UNIT_STATE) mode, n_particles, n_batches

    ! Read batch number to restart at
    read(UNIT_STATE) restart_batch

    ! Read information specific to criticality run
    if (mode == MODE_CRITICALITY) then
       read(UNIT_STATE) n_inactive, gen_per_batch
       read(UNIT_STATE) k_batch(1:restart_batch)
       read(UNIT_STATE) entropy(1:restart_batch)
    end if

    if (master) then
       ! Read number of global tallies and make sure it matches
       read(UNIT_STATE) temp(1)
       if (temp(1) /= N_GLOBAL_TALLIES) then
          message = "Number of global tallies does not match in state point."
          call fatal_error()
       end if

       ! Read global tally data
       do i = 1, N_GLOBAL_TALLIES
          read(UNIT_STATE) global_tallies(i) % sum
          read(UNIT_STATE) global_tallies(i) % sum_sq
       end do

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
                deallocate(int_array)
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
       read(UNIT_STATE) temp(1)
       if (temp(1) == 1) then
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

    ! Read source bank for criticality run
    if (mode == MODE_CRITICALITY) then
       read(UNIT_STATE) source_bank
    end if

    ! Close binary state point file
    close(UNIT_STATE)

  end subroutine load_state_point

!===============================================================================
! REPLAY_BATCH_HISTORY displays batch keff and entropy for each batch stored in
! a state point file
!===============================================================================

  subroutine replay_batch_history

    real(8), save :: temp(2) = ZERO ! temporary values for keff
    real(8)       :: alpha          ! significance level for CI
    real(8)       :: t_value        ! t-value for confidence intervals

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

       ! calculate mean keff
       keff = temp(1) / n_realizations

       if (n_realizations > 1) then
          if (confidence_intervals) then
             ! Calculate t-value for confidence intervals
             alpha = ONE - CONFIDENCE_LEVEL
             t_value = t_percentile(ONE - alpha/TWO, n_realizations - 1)
          else
             t_value = ONE
          end if

          keff_std = t_value * sqrt((temp(2)/n_realizations - keff*keff) &
               / (n_realizations - 1))
       end if
    else
       keff = k_batch(current_batch)
    end if

    ! print out batch keff
    if (master) call print_batch_keff()

    ! Write message at end
    if (current_batch == restart_batch) then
       message = "Resuming simulation..."
       call write_message(1)
    end if

  end subroutine replay_batch_history

end module state_point
