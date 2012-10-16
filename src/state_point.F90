module state_point

!===============================================================================
! STATE_POINT -- This module handles writing and reading binary state point
! files. State points are contain complete tally results, source sites, and
! various other data. They can be used to restart a run or to reconstruct
! confidence intervals for tallies (this requires post-processing via Python
! scripts).
!
! Modifications to this module should be made with care. There are essentially
! three different ways to write or read state points: 1) normal Fortran file
! I/O, 2) MPI-IO, and 3) HDF5. The HDF5 functionality is contained in the
! hdf5_interface module. If you plan to change the state point, you will need to
! change all methods. You should also increment REVISION_STATEPOINT in the
! constants module.
!
! State points can be written at any batch during a simulation, or at specified
! intervals, using the <state_point ... /> tag.
!===============================================================================

  use error,        only: warning, fatal_error
  use global
  use math,         only: t_percentile
  use output,       only: write_message, print_batch_keff, time_stamp
  use source,       only: write_source_binary
  use string,       only: to_str
  use tally_header, only: TallyObject
  use tally,        only: setup_active_usertallies

#ifdef MPI
  use mpi
#endif

  implicit none

contains

!===============================================================================
! WRITE_STATE_POINT creates a state point binary file that can be used for
! restarting a run or for getting intermediate tally results
!===============================================================================

  subroutine write_state_point()

    integer :: i                       ! loop index
#ifdef MPI
    integer :: fh                      ! file handle
    integer :: n                       ! temporary array length
    integer :: temp                    ! temporary variable
    integer :: size_offset_kind        ! size of MPI_OFFSET_KIND (bytes)
    integer :: size_bank               ! size of MPI_BANK type
    integer(MPI_OFFSET_KIND) :: offset ! offset in memory (0=beginning of file)
#else
    integer :: j, k                    ! loop indices
#endif
    type(TallyObject), pointer :: t => null()

    ! Set filename for binary state point
    path_state_point = 'statepoint.' // trim(to_str(current_batch)) // '.binary'

    ! Write message
    message = "Creating state point " // trim(path_state_point) // "..."
    call write_message(1)

#ifdef MPI
    ! ==========================================================================
    ! PARALLEL I/O USING MPI-2 ROUTINES

    ! Open binary source file for reading
    call MPI_FILE_OPEN(MPI_COMM_WORLD, path_state_point, MPI_MODE_CREATE + &
         MPI_MODE_WRONLY, MPI_INFO_NULL, fh, mpi_err)

    ! ==========================================================================
    ! RUN INFORMATION AND TALLY METADATA

    if (master) call write_state_point_header(fh)

    ! ==========================================================================
    ! TALLY RESULTS

    if (.not. reduce_tallies) then
       ! If using the no-tally-reduction method, we need to collect tally
       ! results before writing them to the state point file.

       call write_tally_scores_nr(fh)

    elseif (master) then
       ! Write number of realizations
       call MPI_FILE_WRITE(fh, n_realizations, 1, MPI_INTEGER, &
            MPI_STATUS_IGNORE, mpi_err)

       ! Write global tallies
       call MPI_FILE_WRITE(fh, N_GLOBAL_TALLIES, 1, MPI_INTEGER, &
            MPI_STATUS_IGNORE, mpi_err)
       call MPI_FILE_WRITE(fh, global_tallies, N_GLOBAL_TALLIES, &
            MPI_TALLYSCORE, MPI_STATUS_IGNORE, mpi_err)

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

    ! ==========================================================================
    ! SOURCE BANK

    if (run_mode == MODE_CRITICALITY) then
       if (source_separate) then
          ! If the user has specified that the source sites should be written in
          ! a separate file, we make a call to the appropriate subroutine to
          ! write it separately

          path_source = "source." // trim(to_str(current_batch)) // ".binary"
          call write_source_binary()
       else
          ! Otherwise, write the source sites in the state point file

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

    ! Write current date and time
    write(UNIT_STATE) time_stamp()

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

       ! Number of realizations
       write(UNIT_STATE) t % n_realizations

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

    ! Number of realizations for global tallies
    write(UNIT_STATE) n_realizations

    ! Write out global tallies sum and sum_sq
    write(UNIT_STATE) N_GLOBAL_TALLIES
    GLOBAL_TALLIES_LOOP: do i = 1, N_GLOBAL_TALLIES
       write(UNIT_STATE) global_tallies(i) % sum
       write(UNIT_STATE) global_tallies(i) % sum_sq
    end do GLOBAL_TALLIES_LOOP

    if (tallies_on) then
       ! Indicate that tallies are on
       write(UNIT_STATE) 1

       TALLY_SCORES: do i = 1, n_tallies
          ! Get pointer to tally
          t => tallies(i)

          ! Write tally sum and sum_sq for each bin
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
       if (source_separate) then
          ! If the user has specified that the source sites should be written in
          ! a separate file, we make a call to the appropriate subroutine to
          ! write it separately

          path_source = "source." // trim(to_str(current_batch)) // ".binary"
          call write_source_binary()
       else
          ! Otherwise, write the source sites in the state point file

          write(UNIT_STATE) source_bank
       end if
    end if

    ! Close binary state point file
    close(UNIT_STATE)
#endif

  end subroutine write_state_point

#ifdef MPI
!===============================================================================
! WRITE_STATE_POINT_HEADER uses MPI-IO routines to write basic run information
! and tally metadata
!===============================================================================

  subroutine write_state_point_header(fh)

    integer, intent(inout) :: fh ! file handle

    integer       :: i            ! loop index
    integer       :: j            ! loop index
    integer       :: n            ! temporary array length
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

    ! Write current date and time
    call MPI_FILE_WRITE(fh, time_stamp(), 19, MPI_CHARACTER, &
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

       ! Write number of realizations
       call MPI_FILE_WRITE(fh, t % n_realizations, 1, MPI_INTEGER, &
            MPI_STATUS_IGNORE, mpi_err)

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

  end subroutine write_state_point_header
#endif

#ifdef MPI
!===============================================================================
! WRITE_TALLY_SCORES_NR
!===============================================================================

  subroutine write_tally_scores_nr(fh)

    integer, intent(in) :: fh ! file handle

    integer :: i      ! loop index
    integer :: n      ! number of filter bins
    integer :: m      ! number of score bins
    integer :: temp   ! temporary variable
    integer :: n_bins ! total number of bins
    real(8), allocatable :: tally_temp(:,:,:) ! contiguous array of scores
    real(8) :: global_temp(2,N_GLOBAL_TALLIES)
    real(8) :: dummy  ! temporary receive buffer for non-root reduces
    type(TallyObject), pointer :: t => null()

    ! ==========================================================================
    ! COLLECT AND WRITE GLOBAL TALLIES

    if (master) then
       ! Write number of realizations
       call MPI_FILE_WRITE(fh, n_realizations, 1, MPI_INTEGER, &
            MPI_STATUS_IGNORE, mpi_err)

       ! Write number of global tallies
       call MPI_FILE_WRITE(fh, N_GLOBAL_TALLIES, 1, MPI_INTEGER, &
            MPI_STATUS_IGNORE, mpi_err)
    end if

    ! Copy global tallies into temporary array for reducing
    n_bins = 2 * N_GLOBAL_TALLIES
    global_temp(1,:) = global_tallies(:) % sum
    global_temp(2,:) = global_tallies(:) % sum_sq

    if (master) then
       ! The MPI_IN_PLACE specifier allows the master to copy values into a
       ! receive buffer without having a temporary variable
       call MPI_REDUCE(MPI_IN_PLACE, global_temp, n_bins, MPI_REAL8, MPI_SUM, &
            0, MPI_COMM_WORLD, mpi_err)
       
       ! Write out global tallies sum and sum_sq
       call MPI_FILE_WRITE(fh, global_temp, n_bins, MPI_REAL8, &
            MPI_STATUS_IGNORE, mpi_err)

       ! Transfer values to value on master
       if (current_batch == n_batches) then
          global_tallies(:) % sum    = global_temp(1,:)
          global_tallies(:) % sum_sq = global_temp(2,:)
       end if
    else
       ! Receive buffer not significant at other processors
       call MPI_REDUCE(global_temp, dummy, n_bins, MPI_REAL8, MPI_SUM, &
            0, MPI_COMM_WORLD, mpi_err)
    end if

    if (tallies_on) then
       ! Indicate that tallies are on
       if (master) then
          temp = 1
          call MPI_FILE_WRITE(fh, temp, 1, MPI_INTEGER, &
               MPI_STATUS_IGNORE, mpi_err)
       end if

       ! Write all tally scores
       TALLY_SCORES: do i = 1, n_tallies
          t => tallies(i)

          ! Determine size of tally scores array
          m = size(t % scores, 1)
          n = size(t % scores, 2)
          n_bins = m*n*2

          ! Allocate array for storing sums and sums of squares, but
          ! contiguously in memory for each
          allocate(tally_temp(2,m,n))
          tally_temp(1,:,:) = t % scores(:,:) % sum
          tally_temp(2,:,:) = t % scores(:,:) % sum_sq

          if (master) then
             ! The MPI_IN_PLACE specifier allows the master to copy values into
             ! a receive buffer without having a temporary variable
             call MPI_REDUCE(MPI_IN_PLACE, tally_temp, n_bins, MPI_REAL8, &
                  MPI_SUM, 0, MPI_COMM_WORLD, mpi_err)

             ! Write reduced tally results to file
             call MPI_FILE_WRITE(fh, tally_temp, n_bins, MPI_REAL8, &
                  MPI_STATUS_IGNORE, mpi_err)

             ! At the end of the simulation, store the results back in the
             ! regular TallyScores array
             if (current_batch == n_batches) then
                t % scores(:,:) % sum = tally_temp(1,:,:)
                t % scores(:,:) % sum_sq = tally_temp(2,:,:)
             end if
          else
             ! Receive buffer not significant at other processors
             call MPI_REDUCE(tally_temp, dummy, n_bins, MPI_REAL8, MPI_SUM, &
                  0, MPI_COMM_WORLD, mpi_err)
          end if

          ! Deallocate temporary copy of tally results
          deallocate(tally_temp)
       end do TALLY_SCORES
    else
       if (master) then
       ! Indicate that tallies are off
          temp = 0
          call MPI_FILE_WRITE(fh, temp, 1, MPI_INTEGER, &
               MPI_STATUS_IGNORE, mpi_err)
       end if
    end if

  end subroutine write_tally_scores_nr
#endif

!===============================================================================
! LOAD_STATE_POINT loads data from a state point file to either continue a run
! or to print intermediate tally results
!===============================================================================

  subroutine load_state_point()

    integer :: i, j    ! loop indices
    integer :: mode    ! specified run mode
    integer :: temp(3) ! temporary variable
    integer, allocatable :: int_array(:)
    real(8), allocatable :: real_array(:)
    character(19)        :: current_time  ! current date and time

#ifdef MPI
    integer :: fh                      ! file handle
    integer :: n                       ! temporary array size
    integer :: size_offset_kind        ! size of MPI_OFFSET_KIND (bytes)
    integer :: size_bank               ! size of MPI_BANK type
    integer(MPI_OFFSET_KIND) :: offset ! offset in memory (0=beginning of file)
#else
    integer :: k ! loop index
#endif

    ! Write message
    message = "Loading state point " // trim(path_state_point) // "..."
    call write_message(1)

#ifdef MPI
    ! Open binary state point file for reading
    call MPI_FILE_OPEN(MPI_COMM_WORLD, path_state_point, MPI_MODE_RDONLY, &
         MPI_INFO_NULL, fh, mpi_err)

    ! ==========================================================================
    ! RUN INFORMATION AND TALLY METADATA

    ! Raad revision number for state point file and make sure it matches with
    ! current version
    call MPI_FILE_READ_ALL(fh, temp, 1, MPI_INTEGER, MPI_STATUS_IGNORE, mpi_err)
    if (temp(1) /= REVISION_STATEPOINT) then
       message = "State point binary version does not match current version " &
            // "in OpenMC."
       call fatal_error()
    end if
    
    ! Read OpenMC version
    call MPI_FILE_READ_ALL(fh, temp, 3, MPI_INTEGER, MPI_STATUS_IGNORE, mpi_err)
    if (temp(1) /= VERSION_MAJOR .or. temp(2) /= VERSION_MINOR &
         .or. temp(3) /= VERSION_RELEASE) then
       message = "State point file was created with a different version " // &
            "of OpenMC."
       call warning()
    end if

    ! Read date and time
    call MPI_FILE_READ_ALL(fh, current_time, 19, MPI_CHARACTER, &
         MPI_STATUS_IGNORE, mpi_err)

    ! Read and overwrite random number seed
    call MPI_FILE_READ_ALL(fh, seed, 1, MPI_INTEGER8, &
         MPI_STATUS_IGNORE, mpi_err)

    ! Read and overwrite run information
    call MPI_FILE_READ_ALL(fh, mode, 1, MPI_INTEGER, &
         MPI_STATUS_IGNORE, mpi_err)
    call MPI_FILE_READ_ALL(fh, n_particles, 1, MPI_INTEGER8, &
         MPI_STATUS_IGNORE, mpi_err)
    call MPI_FILE_READ_ALL(fh, n_batches, 1, MPI_INTEGER, &
         MPI_STATUS_IGNORE, mpi_err)

    ! Read batch number to restart at
    call MPI_FILE_READ_ALL(fh, restart_batch, 1, MPI_INTEGER, &
         MPI_STATUS_IGNORE, mpi_err)

    ! Read information specific to criticality run
    if (mode == MODE_CRITICALITY) then
       call MPI_FILE_READ_ALL(fh, n_inactive, 1, MPI_INTEGER, &
            MPI_STATUS_IGNORE, mpi_err)
       call MPI_FILE_READ_ALL(fh, gen_per_batch, 1, MPI_INTEGER, &
            MPI_STATUS_IGNORE, mpi_err)
       call MPI_FILE_READ_ALL(fh, k_batch, restart_batch, MPI_REAL8, &
            MPI_STATUS_IGNORE, mpi_err)
       call MPI_FILE_READ_ALL(fh, entropy, restart_batch, MPI_REAL8, &
            MPI_STATUS_IGNORE, mpi_err)
    end if

    if (master) then
       ! Read number of meshes
       call MPI_FILE_READ(fh, temp, 1, MPI_INTEGER, MPI_STATUS_IGNORE, mpi_err)
       if (temp(1) /= n_meshes) then
          message = "Number of meshes does not match in state point."
          call fatal_error()
       end if

       MESH_LOOP: do i = 1, n_meshes
          ! Read type of mesh and dimension
          call MPI_FILE_READ(fh, temp, 2, MPI_INTEGER, MPI_STATUS_IGNORE, &
               mpi_err)

          ! Skip mesh data
          call MPI_FILE_GET_POSITION(fh, offset, mpi_err)
          offset = offset + temp(2)*(4 + 3*8)
          call MPI_FILE_SEEK(fh, offset, MPI_SEEK_SET, mpi_err)
       end do MESH_LOOP

       ! Read number of tallies and make sure it matches
       call MPI_FILE_READ(fh, temp, 1, MPI_INTEGER, MPI_STATUS_IGNORE, mpi_err)
       if (temp(1) /= n_tallies) then
          message = "Number of tallies does not match in state point."
          call fatal_error()
       end if

       TALLY_METADATA: do i = 1, n_tallies
          ! Read number of realizations for global tallies
          call MPI_FILE_READ(fh, tallies(i) % n_realizations, 1, &
               MPI_INTEGER, MPI_STATUS_IGNORE, mpi_err)

          ! Read dimensions of tally filters and scores and make sure they
          ! match
          call MPI_FILE_READ(fh, temp, 2, MPI_INTEGER, &
               MPI_STATUS_IGNORE, mpi_err)
          if (temp(1) /= size(tallies(i) % scores, 1) .or. &
               temp(2) /= size(tallies(i) % scores, 2)) then
             message = "Tally dimensions do not match in state point."
             call fatal_error()
          end if

          ! Read number of filters
          call MPI_FILE_READ(fh, temp, 1, MPI_INTEGER, &
               MPI_STATUS_IGNORE, mpi_err)

          FILTER_LOOP: do j = 1, temp(1)
             ! Read filter type and number of bins
             call MPI_FILE_READ(fh, temp(2), 2, MPI_INTEGER, &
                  MPI_STATUS_IGNORE, mpi_err)

             ! Read filter bins
             select case (temp(2))
             case (FILTER_MESH)
                allocate(int_array(1))
                call MPI_FILE_READ(fh, int_array, 1, MPI_INTEGER, &
                     MPI_STATUS_IGNORE, mpi_err)
                deallocate(int_array)
             case (FILTER_ENERGYIN, FILTER_ENERGYOUT)
                allocate(real_array(temp(3) + 1))
                call MPI_FILE_READ(fh, real_array, temp(3) + 1, MPI_REAL8, &
                     MPI_STATUS_IGNORE, mpi_err)
                deallocate(real_array)
             case default
                allocate(int_array(temp(3)))
                call MPI_FILE_READ(fh, int_array, temp(3), MPI_INTEGER, &
                     MPI_STATUS_IGNORE, mpi_err)
                deallocate(int_array)
             end select
          end do FILTER_LOOP

          ! Read number of nuclides
          call MPI_FILE_READ(fh, temp, 1, MPI_INTEGER, &
               MPI_STATUS_IGNORE, mpi_err)

          ! Read nuclide bins
          allocate(int_array(temp(1)))
          call MPI_FILE_READ(fh, int_array, temp(1), MPI_INTEGER, &
               MPI_STATUS_IGNORE, mpi_err)
          deallocate(int_array)

          ! Read number of scores
          call MPI_FILE_READ(fh, temp, 1, MPI_INTEGER, &
               MPI_STATUS_IGNORE, mpi_err)

          ! Read nuclide bins
          allocate(int_array(temp(1)))
          call MPI_FILE_READ(fh, int_array, temp(1), MPI_INTEGER, &
               MPI_STATUS_IGNORE, mpi_err)
          deallocate(int_array)
       end do TALLY_METADATA

       ! Read number of realizations for global tallies
       call MPI_FILE_READ(fh, n_realizations, 1, MPI_INTEGER, &
            MPI_STATUS_IGNORE, mpi_err)

       ! Read number of global tallies and make sure it matches
       call MPI_FILE_READ(fh, temp, 1, MPI_INTEGER, MPI_STATUS_IGNORE, mpi_err)
       if (temp(1) /= N_GLOBAL_TALLIES) then
          message = "Number of global tallies does not match in state point."
          call fatal_error()
       end if

       ! Read global tally data
       call MPI_FILE_READ(fh, global_tallies, N_GLOBAL_TALLIES, &
            MPI_TALLYSCORE, MPI_STATUS_IGNORE, mpi_err)

       ! Check if tally results are present
       call MPI_FILE_READ(fh, temp, 1, MPI_INTEGER, MPI_STATUS_IGNORE, mpi_err)

       ! =======================================================================
       ! TALLY RESULTS

       ! Read sum and sum squared
       if (temp(1) == 1) then
          TALLY_SCORES: do i = 1, n_tallies
             n = size(tallies(i) % scores, 1) * size(tallies(i) % scores, 2)
             call MPI_FILE_READ(fh, tallies(i) % scores, n, MPI_TALLYSCORE, &
                  MPI_STATUS_IGNORE, mpi_err)
          end do TALLY_SCORES
       end if
    end if

    ! ==========================================================================
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
       call MPI_FILE_READ_AT(fh, offset, source_bank(1), work, MPI_BANK, &
            MPI_STATUS_IGNORE, mpi_err)
    end if

    ! Close binary state point file
    call MPI_FILE_CLOSE(fh, mpi_err)

#else
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

    ! Read date and time
    read(UNIT_STATE) current_time

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
          ! Read number of realizations
          read(UNIT_STATE) tallies(i) % n_realizations

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

       ! Read number of realizations for global tallies
       read(UNIT_STATE) n_realizations

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
    if (mode == MODE_CRITICALITY .and. run_mode /= MODE_TALLIES) then
       read(UNIT_STATE) source_bank
    end if

    ! Close binary state point file
    close(UNIT_STATE)
#endif

  end subroutine load_state_point

!===============================================================================
! REPLAY_BATCH_HISTORY displays batch keff and entropy for each batch stored in
! a state point file
!===============================================================================

  subroutine replay_batch_history

    integer       :: n = 0          ! number of realizations
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
    if (current_batch == n_inactive) then
       tallies_on = .true.
       call setup_active_usertallies()
    end if

    ! Add to number of realizations
    if (current_batch > n_inactive) then
       n = n + 1

       temp(1) = temp(1) + k_batch(current_batch)
       temp(2) = temp(2) + k_batch(current_batch)*k_batch(current_batch)

       ! calculate mean keff
       keff = temp(1) / n

       if (n > 1) then
          if (confidence_intervals) then
             ! Calculate t-value for confidence intervals
             alpha = ONE - CONFIDENCE_LEVEL
             t_value = t_percentile(ONE - alpha/TWO, n - 1)
          else
             t_value = ONE
          end if

          keff_std = t_value * sqrt((temp(2)/n - keff*keff)/(n - 1))
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
