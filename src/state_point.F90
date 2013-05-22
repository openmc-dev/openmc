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


  use constants
  use error,              only: fatal_error, warning
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
    integer                 :: i
    integer                 :: j
    integer, allocatable    :: temp_array(:)
    type(TallyObject), pointer :: t => null()

    ! Set filename for state point
    filename = trim(path_output) // 'statepoint.' // &
               trim(to_str(current_batch))

    ! Append appropriate extension
#ifdef HDF5
    filename = trim(filename) // '.h5'
#else
    filename = trim(filename) // '.binary'
#endif

    ! Write message
    message = "Creating state point " // trim(filename) // "..."
    call write_message(1)

    ! Create statepoint file 
    call file_create(filename, 'serial')

    if (master) then
      ! Write file type
      call write_data(FILETYPE_STATEPOINT, "filetype")

      ! Write revision number for state point file
      call write_data(REVISION_STATEPOINT, "revision")

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
        call write_data(k_generation, "k_generation", &
             length=current_batch*gen_per_batch)
        call write_data(entropy, "entropy", length=current_batch*gen_per_batch)
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

        ! Write number of user score bins
        call write_data(t % n_user_score_bins, "n_user_score_bins", &
             group="tallies/tally" // to_str(i))

      end do TALLY_METADATA

    end if

    ! Check for the no-tally-reduction method
    if (.not. reduce_tallies) then
      ! If using the no-tally-reduction method, we need to collect tally
      ! results before writing them to the state point file.

      call write_tally_results_nr()

    elseif (master) then

      ! Write number of global realizations
      call write_data(n_realizations, "n_realizations")

      ! Write global tallies
      call write_data(N_GLOBAL_TALLIES, "n_global_tallies")
      call write_tally_result(global_tallies, "global_tallies", &
           n1=N_GLOBAL_TALLIES, n2=1)

      ! Write tallies
      if (tallies_on) then

        ! Indicate that tallies are on
        call write_data(1, "tallies_present", group="tallies")

        ! Write all tally results
        TALLY_RESULTS: do i = 1, n_tallies

          ! Set point to current tally
          t => tallies(i)

          ! Write sum and sum_sq for each bin
          call write_tally_result(t % results, "results", &
               group="tallies/tally" // to_str(i), &
               n1=size(t % results, 1), n2=size(t % results, 2))

        end do TALLY_RESULTS

      else

        ! Indicate tallies are off
        call write_data(0, "tallies_present", group="tallies")

      end if

    end if

    ! Check for eigenvalue calculation
    if (run_mode == MODE_EIGENVALUE .and. source_write) then

      ! Check for writing source out separately
      if (source_separate) then

        ! Close statepoint file 
        call file_close('serial')

        ! Set filename for source
        filename = trim(path_output) // 'source.' // &
                   trim(to_str(current_batch))
#ifdef HDF5
        filename = trim(filename) // '.h5'
#else
        filename = trim(filename) // '.binary'
#endif

        ! Write message
        message = "Creating source file " // trim(filename) // "..."
        call write_message(1)

        ! Create statepoint file 
        call file_create(filename, 'parallel')

#ifdef HDF5
# ifdef MPI
      else
        ! Close HDF5 serial file and reopen in parallel
        call file_close('serial')
        call file_open(filename, 'parallel', 'w') 
# endif
#endif

      end if

      ! Write out source
      call write_source_bank()

      ! Close file, all files in parallel mode
      call file_close('parallel') ! even if no MPI, this will work for HDF5

    else

      ! Close file if not in eigenvalue mode or no source writing
      call file_close('serial')

    end if

  end subroutine write_state_point

!===============================================================================
! WRITE_TALLY_RESULTS_NR
!===============================================================================

  subroutine write_tally_results_nr()

    integer :: i      ! loop index
    integer :: n      ! number of filter bins
    integer :: m      ! number of score bins
    integer :: n_bins ! total number of bins
    real(8), allocatable :: tally_temp(:,:,:) ! contiguous array of results
    real(8), target :: global_temp(2,N_GLOBAL_TALLIES)
    real(8) :: dummy  ! temporary receive buffer for non-root reduces
    type(TallyObject), pointer :: t => null()
    type(TallyResult), allocatable :: tallyresult_temp(:,:)

    ! ==========================================================================
    ! COLLECT AND WRITE GLOBAL TALLIES

    if (master) then
      ! Write number of realizations
      call write_data(n_realizations, "n_realizations")

      ! Write number of global tallies
      call write_data(N_GLOBAL_TALLIES, "n_global_tallies")
    end if

    ! Copy global tallies into temporary array for reducing
    n_bins = 2 * N_GLOBAL_TALLIES
    global_temp(1,:) = global_tallies(:) % sum
    global_temp(2,:) = global_tallies(:) % sum_sq 

    if (master) then
      ! The MPI_IN_PLACE specifier allows the master to copy values into a
      ! receive buffer without having a temporary variable
#ifdef MPI
      call MPI_REDUCE(MPI_IN_PLACE, global_temp, n_bins, MPI_REAL8, MPI_SUM, &
           0, MPI_COMM_WORLD, mpi_err)
#endif

      ! Transfer values to value on master
      if (current_batch == n_batches) then
        global_tallies(:) % sum    = global_temp(1,:)
        global_tallies(:) % sum_sq = global_temp(2,:)
      end if

      ! Put reduced value in temporary tally result
      allocate(tallyresult_temp(N_GLOBAL_TALLIES, 1))
      tallyresult_temp(:,1) % sum    = global_temp(1,:)
      tallyresult_temp(:,1) % sum_sq = global_temp(2,:)

   
      ! Write out global tallies sum and sum_sq
      call write_tally_result(tallyresult_temp, "global_tallies", &
           n1=N_GLOBAL_TALLIES, n2=1)

      ! Deallocate temporary tally result
      deallocate(tallyresult_temp)
    else
      ! Receive buffer not significant at other processors
#ifdef MPI
      call MPI_REDUCE(global_temp, dummy, n_bins, MPI_REAL8, MPI_SUM, &
           0, MPI_COMM_WORLD, mpi_err)
#endif
    end if

    if (tallies_on) then
      ! Indicate that tallies are on
      if (master) then
        call write_data(1, "tallies_present", group="tallies")
      end if

      ! Write all tally results
      TALLY_RESULTS: do i = 1, n_tallies
        t => tallies(i)

        ! Determine size of tally results array
        m = size(t % results, 1)
        n = size(t % results, 2)
        n_bins = m*n*2

        ! Allocate array for storing sums and sums of squares, but
        ! contiguously in memory for each
        allocate(tally_temp(2,m,n))
        tally_temp(1,:,:) = t % results(:,:) % sum
        tally_temp(2,:,:) = t % results(:,:) % sum_sq

        if (master) then
          ! The MPI_IN_PLACE specifier allows the master to copy values into
          ! a receive buffer without having a temporary variable
#ifdef MPI
          call MPI_REDUCE(MPI_IN_PLACE, tally_temp, n_bins, MPI_REAL8, &
               MPI_SUM, 0, MPI_COMM_WORLD, mpi_err)
#endif
          ! At the end of the simulation, store the results back in the
          ! regular TallyResults array
          if (current_batch == n_batches) then
            t % results(:,:) % sum = tally_temp(1,:,:)
            t % results(:,:) % sum_sq = tally_temp(2,:,:)
          end if

         ! Put in temporary tally result
         allocate(tallyresult_temp(m,n))
         tallyresult_temp(:,:) % sum    = tally_temp(1,:,:)
         tallyresult_temp(:,:) % sum_sq = tally_temp(2,:,:)
 
         ! Write reduced tally results to file
          call write_tally_result(t % results, "results", &
               group="tallies/tally" // to_str(i), n1=m, n2=n)

          ! Deallocate temporary tally result
          deallocate(tallyresult_temp)
        else
          ! Receive buffer not significant at other processors
#ifdef MPI
          call MPI_REDUCE(tally_temp, dummy, n_bins, MPI_REAL8, MPI_SUM, &
               0, MPI_COMM_WORLD, mpi_err)
#endif
        end if

        ! Deallocate temporary copy of tally results
        deallocate(tally_temp)
      end do TALLY_RESULTS
    else
      if (master) then
        ! Indicate that tallies are off
        call write_data(0, "tallies_present", group="tallies")
      end if
    end if

  end subroutine write_tally_results_nr

!===============================================================================
! LOAD_STATE_POINT
!===============================================================================

  subroutine load_state_point()

    character(MAX_FILE_LEN) :: filename
    character(MAX_FILE_LEN) :: path_temp
    character(19)           :: current_time
    integer                 :: i
    integer                 :: j
    integer                 :: int_array(3)
    integer, allocatable    :: temp_array(:)
    real(8)                 :: real_array(3) 
    type(TallyObject), pointer :: t => null()

    ! Write message
    message = "Loading state point " // trim(path_state_point) // "..."
    call write_message(1)

    ! Open file for reading
    call file_open(path_state_point, 'parallel', 'r')

    ! Read filetype
    call read_data(int_array(1), "filetype", option="collective")

    ! Read revision number for state point file and make sure it matches with
    ! current version
    call read_data(int_array(1), "revision", option="collective")
    if (int_array(1) /= REVISION_STATEPOINT) then
      message = "State point version does not match current version " &
                // "in OpenMC."
      call fatal_error()
    end if

    ! Read OpenMC version
    call read_data(int_array(1), "version_major", option="collective")
    call read_data(int_array(2), "version_minor", option="collective")
    call read_data(int_array(3), "version_release", option="collective")
    if (int_array(1) /= VERSION_MAJOR .or. int_array(2) /= VERSION_MINOR &
        .or. int_array(3) /= VERSION_RELEASE) then
      message = "State point file was created with a different version " &
                // "of OpenMC."
      call warning()
    end if

    ! Read date and time
    call read_data(current_time, "date_and_time", option="collective")

    ! Read path to input
    call read_data(path_temp, "path", option="collective")

    ! Read and overwrite random number seed
    call read_data(seed, "seed", option="collective")

    ! Read and overwrite run information except number of batches
    call read_data(run_mode, "run_mode", option="collective")
    call read_data(n_particles, "n_particles", option="collective")
    call read_data(int_array(1), "n_batches", option="collective")

    ! Take maximum of statepoint n_batches and input n_batches
    n_batches = max(n_batches, int_array(1))

    ! Read batch number to restart at
    call read_data(restart_batch, "current_batch", option="collective")

    ! Read information specific to eigenvalue run
    if (run_mode == MODE_EIGENVALUE) then
      call read_data(int_array(1), "n_inactive", option="collective")
      call read_data(gen_per_batch, "gen_per_batch", option="collective")
      call read_data(k_generation, "k_generation", &
           length=restart_batch*gen_per_batch, option="collective")
      call read_data(entropy, "entropy", length=restart_batch*gen_per_batch, &
           option="collective")
      call read_data(k_col_abs, "k_col_abs", option="collective")
      call read_data(k_col_tra, "k_col_tra", option="collective")
      call read_data(k_abs_tra, "k_abs_tra", option="collective")
      call read_data(real_array(1:2), "k_combined", length=2, &
            option="collective")

      ! Take maximum of statepoint n_inactive and input n_inactive
      n_inactive = max(n_inactive, int_array(1))
    end if

    ! Read number of meshes
    call read_data(n_meshes, "n_meshes", group="tallies", option="collective")

    ! Read and overwrite mesh information
    MESH_LOOP: do i = 1, n_meshes
      call read_data(meshes(i) % id, "id", &
           group="tallies/mesh" // to_str(i), option="collective")
      call read_data(meshes(i) % type, "type", &
           group="tallies/mesh" // to_str(i), option="collective")
      call read_data(meshes(i) % n_dimension, "n_dimension", &
           group="tallies/mesh" // to_str(i), option="collective")
      call read_data(meshes(i) % dimension, "dimension", &
           group="tallies/mesh" // to_str(i), &
           length=meshes(i) % n_dimension, option="collective")
      call read_data(meshes(i) % lower_left, "lower_left", &
           group="tallies/mesh" // to_str(i), &
           length=meshes(i) % n_dimension, option="collective")
      call read_data(meshes(i) % upper_right, "upper_right", &
           group="tallies/mesh" // to_str(i), &
           length=meshes(i) % n_dimension, option="collective")
      call read_data(meshes(i) % width, "width", &
           group="tallies/mesh" // to_str(i), &
           length=meshes(i) % n_dimension, option="collective")
    end do MESH_LOOP

    ! Read and overwrite number of tallies
    call read_data(n_tallies, "n_tallies", group="tallies", option="collective")

    ! Read in tally metadata
    TALLY_METADATA: do i = 1, n_tallies

      ! Get pointer to tally
      t => tallies(i)

      ! Read tally id
      call read_data(t % id, "id", group="tallies/tally" // to_str(i), &
           option="collective")

      ! Read number of realizations
      call read_data(t % n_realizations, "n_realizations", &
           group="tallies/tally" // to_str(i), option="collective")

      ! Read size of tally results
      call read_data(int_array(1), "total_score_bins", &
           group="tallies/tally" // to_str(i), option="collective")
      call read_data(int_array(2), "total_filter_bins", &
           group="tallies/tally" // to_str(i), option="collective")

      ! Check size of tally results array
      if (int_array(1) /= t % total_score_bins .and. &
          int_array(2) /= t % total_filter_bins) then
        message = "Input file tally structure is different from restart."
        call fatal_error()
      end if

      ! Read number of filters
      call read_data(t % n_filters, "n_filters", &
           group="tallies/tally" // to_str(i), option="collective")

      ! Read filter information
      FILTER_LOOP: do j = 1, t % n_filters

        ! Read type of filter
        call read_data(t % filters(j) % type, "type", &
             group="tallies/tally" // trim(to_str(i)) // "/filter" // to_str(j), &
             option="collective")

        ! Read number of bins for this filter
        call read_data(t % filters(j) % n_bins, "n_bins", &
             group="tallies/tally" // trim(to_str(i)) // "/filter" // to_str(j), &
             option="collective")

        ! Read bins
        if (t % filters(j) % type == FILTER_ENERGYIN .or. &
            t % filters(j) % type == FILTER_ENERGYOUT) then
          call read_data(t % filters(j) % real_bins, "bins", &
               group="tallies/tally" // trim(to_str(i)) // "/filter" // to_str(j), &
               length=size(t % filters(j) % real_bins), option="collective")
        else
          call read_data(t % filters(j) % int_bins, "bins", &
               group="tallies/tally" // trim(to_str(i)) // "/filter" // to_str(j), &
               length=size(t % filters(j) % int_bins), option="collective")
        end if

      end do FILTER_LOOP

      ! Read number of nuclide bins
      call read_data(t % n_nuclide_bins, "n_nuclide_bins", &
           group="tallies/tally" // to_str(i), option="collective")

      ! Set up nuclide bin array and then write
      allocate(temp_array(t % n_nuclide_bins))
      call read_data(temp_array, "nuclide_bins", &
           group="tallies/tally" // to_str(i), length=t % n_nuclide_bins, &
           option="collective")
      NUCLIDE_LOOP: do j = 1, t % n_nuclide_bins
        if (temp_array(j) > 0) then
          nuclides(t % nuclide_bins(j)) % zaid = temp_array(j)
        else
          t % nuclide_bins(j) = temp_array(j)
        end if
      end do NUCLIDE_LOOP
      deallocate(temp_array)

      ! Write number of score bins, score bins, and scatt order
      call read_data(t % n_score_bins, "n_score_bins", &
           group="tallies/tally" // to_str(i), option="collective")
      call read_data(t % score_bins, "score_bins", &
           group="tallies/tally" // to_str(i), length=t % n_score_bins, &
           option="collective")
      call read_data(t % scatt_order, "scatt_order", &
           group="tallies/tally" // to_str(i), length=t % n_score_bins, &
           option="collective")

      ! Write number of user score bins
      call read_data(t % n_user_score_bins, "n_user_score_bins", &
           group="tallies/tally" // to_str(i), option="collective")

    end do TALLY_METADATA

    ! Read tallies to master
    if (master) then

      ! Read number of realizations for global tallies
      call read_data(n_realizations, "n_realizations", option="independent")

      ! Read number of global tallies
      call read_data(int_array(1), "n_global_tallies", option="independent")
      if (int_array(1) /= N_GLOBAL_TALLIES) then
        message = "Number of global tallies does not match in state point."
        call fatal_error()
      end if

      ! Read global tally data
      call read_tally_result(global_tallies, "global_tallies", &
           n1=N_GLOBAL_TALLIES, n2=1)

      ! Check if tally results are present
      call read_data(int_array(1), "tallies_present", group="tallies", &
           option="independent")

      ! Read in sum and sum squared
      if (int_array(1) == 1) then
        TALLY_RESULTS: do i = 1, n_tallies

          ! Set pointer to tally
          t => tallies(i)

          ! Read sum and sum_sq for each bin
          call read_tally_result(t % results, "results", &
               group="tallies/tally" // to_str(i), &
               n1=size(t % results, 1), n2=size(t % results, 2))

        end do TALLY_RESULTS
      end if
    end if

    ! Read source if in eigenvalue mode 
    if (run_mode == MODE_EIGENVALUE .and. run_mode /= MODE_TALLIES) then

      ! Check if source was written out separately
      if (source_separate) then

        ! Close statepoint file 
        call file_close('parallel')

        ! Set filename for source
        filename = trim(path_output) // 'source.' // &
                   trim(to_str(restart_batch))
#ifdef HDF5
        filename = trim(filename) // '.h5'
#else
        filename = trim(filename) // '.binary'
#endif

        ! Write message
        message = "Loading source file " // trim(filename) // "..."
        call write_message(1)

        ! Create statepoint file
        call file_open(filename, 'parallel', 'r')

      end if

      ! Write out source
      call read_source_bank()

      ! Close file
      if (source_separate) then
        call file_close('parallel')
      else
        call file_close('parallel')
      end if

    else

      ! Close file if not in eigenvalue mode
      call file_close('parallel')

    end if

  end subroutine load_state_point

  subroutine read_source
! TODO write this routine
! TODO what if n_particles does not match source bank
  end subroutine read_source
end module state_point
