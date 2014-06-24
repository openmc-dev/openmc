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

#ifdef MPI
  use mpi
#endif

  implicit none

  type(BinaryOutput) :: sp ! statepoint/source output file

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

    if (master) then
      ! Create statepoint file 
      call sp % file_create(filename)

      ! Write file type
      call sp % write_data(FILETYPE_STATEPOINT, "filetype")

      ! Write revision number for state point file
      call sp % write_data(REVISION_STATEPOINT, "revision")

      ! Write OpenMC version
      call sp % write_data(VERSION_MAJOR, "version_major")
      call sp % write_data(VERSION_MINOR, "version_minor")
      call sp % write_data(VERSION_RELEASE, "version_release")

      ! Write current date and time
      call sp % write_data(time_stamp(), "date_and_time")

      ! Write path to input
      call sp % write_data(path_input, "path")

      ! Write out random number seed
      call sp % write_data(seed, "seed")

      ! Write run information
      call sp % write_data(run_mode, "run_mode")
      call sp % write_data(n_particles, "n_particles")
      call sp % write_data(n_batches, "n_batches")

      ! Write out current batch number
      call sp % write_data(current_batch, "current_batch")

      ! Write out information for eigenvalue run
      if (run_mode == MODE_EIGENVALUE) then
        call sp % write_data(n_inactive, "n_inactive")
        call sp % write_data(gen_per_batch, "gen_per_batch")
        call sp % write_data(k_generation, "k_generation", &
             length=current_batch*gen_per_batch)
        call sp % write_data(entropy, "entropy", length=current_batch*gen_per_batch)
        call sp % write_data(k_col_abs, "k_col_abs")
        call sp % write_data(k_col_tra, "k_col_tra")
        call sp % write_data(k_abs_tra, "k_abs_tra")
        call sp % write_data(k_combined, "k_combined", length=2)

        ! Write out CMFD info
        if (cmfd_on) then
          call sp % write_data(1, "cmfd_on")
          call sp % write_data(cmfd % indices, "indices", length=4, group="cmfd")
          call sp % write_data(cmfd % k_cmfd, "k_cmfd", length=current_batch, &
               group="cmfd")
          call sp % write_data(cmfd % cmfd_src, "cmfd_src", &
               length=(/cmfd % indices(4), cmfd % indices(1), &
               cmfd % indices(2), cmfd % indices(3)/), &
               group="cmfd")
          call sp % write_data(cmfd % entropy, "cmfd_entropy", &
                          length=current_batch, group="cmfd")
          call sp % write_data(cmfd % balance, "cmfd_balance", &
               length=current_batch, group="cmfd")
          call sp % write_data(cmfd % dom, "cmfd_dominance", &
               length = current_batch, group="cmfd")
          call sp % write_data(cmfd % src_cmp, "cmfd_srccmp", &
               length = current_batch, group="cmfd")
        else
          call sp % write_data(0, "cmfd_on")
        end if
      end if

      ! Write number of meshes
      call sp % write_data(n_meshes, "n_meshes", group="tallies")

      ! Write information for meshes
      MESH_LOOP: do i = 1, n_meshes
        call sp % write_data(meshes(i) % id, "id", &
             group="tallies/mesh" // to_str(i))
        call sp % write_data(meshes(i) % type, "type", &
             group="tallies/mesh" // to_str(i))
        call sp % write_data(meshes(i) % n_dimension, "n_dimension", &
             group="tallies/mesh" // to_str(i))
        call sp % write_data(meshes(i) % dimension, "dimension", &
             group="tallies/mesh" // to_str(i), &
             length=meshes(i) % n_dimension)
        call sp % write_data(meshes(i) % lower_left, "lower_left", &
             group="tallies/mesh" // to_str(i), &
             length=meshes(i) % n_dimension)
        call sp % write_data(meshes(i) % upper_right, "upper_right", &
             group="tallies/mesh" // to_str(i), &
             length=meshes(i) % n_dimension)
        call sp % write_data(meshes(i) % width, "width", &
             group="tallies/mesh" // to_str(i), &
             length=meshes(i) % n_dimension)
      end do MESH_LOOP

      ! Write number of tallies
      call sp % write_data(n_tallies, "n_tallies", group="tallies")

      ! Write all tally information except results
      TALLY_METADATA: do i = 1, n_tallies
        !Get pointer to tally
        t => tallies(i)

        ! Write id
        call sp % write_data(t % id, "id", group="tallies/tally" // to_str(i))

        ! Write number of realizations
        call sp % write_data(t % n_realizations, "n_realizations", &
             group="tallies/tally" // to_str(i))

        ! Write size of each tally
        call sp % write_data(t % total_score_bins, "total_score_bins", &
             group="tallies/tally" // to_str(i))
        call sp % write_data(t % total_filter_bins, "total_filter_bins", &
             group="tallies/tally" // to_str(i))

        ! Write number of filters
        call sp % write_data(t % n_filters, "n_filters", &
             group="tallies/tally" // to_str(i))

        ! Write filter information
        FILTER_LOOP: do j = 1, t % n_filters

          ! Write type of filter
          call sp % write_data(t % filters(j) % type, "type", &
               group="tallies/tally" // trim(to_str(i)) // "/filter" // to_str(j))

          ! Write number of bins for this filter
          call sp % write_data(t % filters(j) % n_bins, "n_bins", &
               group="tallies/tally" // trim(to_str(i)) // "/filter" // to_str(j))

          ! Write bins
          if (t % filters(j) % type == FILTER_ENERGYIN .or. &
              t % filters(j) % type == FILTER_ENERGYOUT) then
            call sp % write_data(t % filters(j) % real_bins, "bins", &
                 group="tallies/tally" // trim(to_str(i)) // "/filter" // to_str(j), &
                 length=size(t % filters(j) % real_bins))
          else
            call sp % write_data(t % filters(j) % int_bins, "bins", &
                 group="tallies/tally" // trim(to_str(i)) // "/filter" // to_str(j), &
                 length=size(t % filters(j) % int_bins))
          end if

        end do FILTER_LOOP

        ! Write number of nuclide bins
        call sp % write_data(t % n_nuclide_bins, "n_nuclide_bins", &
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
        call sp % write_data(temp_array, "nuclide_bins", &
             group="tallies/tally" // to_str(i), length=t % n_nuclide_bins)
        deallocate(temp_array)

        ! Write number of score bins, score bins, and scatt order
        call sp % write_data(t % n_score_bins, "n_score_bins", &
             group="tallies/tally" // to_str(i))
        call sp % write_data(t % score_bins, "score_bins", &
             group="tallies/tally" // to_str(i), length=t % n_score_bins)
        call sp % write_data(t % scatt_order, "scatt_order", &
             group="tallies/tally" // to_str(i), length=t % n_score_bins)

        ! Write number of user score bins
        call sp % write_data(t % n_user_score_bins, "n_user_score_bins", &
             group="tallies/tally" // to_str(i))

      end do TALLY_METADATA

      ! Indicate where source bank is stored in statepoint
      if (source_separate) then
        call sp % write_data(0, "source_present")
      else
        call sp % write_data(1, "source_present")
      end if

    end if

    ! Check for the no-tally-reduction method
    if (.not. reduce_tallies) then
      ! If using the no-tally-reduction method, we need to collect tally
      ! results before writing them to the state point file.

      call write_tally_results_nr()

    elseif (master) then

      ! Write number of global realizations
      call sp % write_data(n_realizations, "n_realizations")

      ! Write global tallies
      call sp % write_data(N_GLOBAL_TALLIES, "n_global_tallies")
      call sp % write_tally_result(global_tallies, "global_tallies", &
           n1=N_GLOBAL_TALLIES, n2=1)

      ! Write tallies
      if (tallies_on) then

        ! Indicate that tallies are on
        call sp % write_data(1, "tallies_present", group="tallies")

        ! Write all tally results
        TALLY_RESULTS: do i = 1, n_tallies

          ! Set point to current tally
          t => tallies(i)

          ! Write sum and sum_sq for each bin
          call sp % write_tally_result(t % results, "results", &
               group="tallies/tally" // to_str(i), &
               n1=size(t % results, 1), n2=size(t % results, 2))

        end do TALLY_RESULTS

      else

        ! Indicate tallies are off
        call sp % write_data(0, "tallies_present", group="tallies")

      end if

      ! Close the file for serial writing
      call sp % file_close()

    end if

  end subroutine write_state_point

!===============================================================================
! WRITE_SOURCE_POINT
!===============================================================================

  subroutine write_source_point()

    type(BinaryOutput) :: sp
    character(MAX_FILE_LEN) :: filename

    ! Check to write out source for a specified batch
    if (sourcepoint_batch % contains(current_batch)) then

      ! Create or open up file
      if (source_separate) then

        ! Set filename
        filename = trim(path_output) // 'source.' // trim(to_str(current_batch))
#ifdef HDF5
        filename = trim(filename) // '.h5'
#else
        filename = trim(filename) // '.binary'
#endif

        ! Write message for new file creation
        message = "Creating source file " // trim(filename) // "..."
        call write_message(1)

        ! Create separate source file
        call sp % file_create(filename, serial = .false.)

        ! Write file type
        call sp % write_data(FILETYPE_SOURCE, "filetype")

      else

        ! Set filename for state point
        filename = trim(path_output) // 'statepoint.' // &
                   trim(to_str(current_batch))
#ifdef HDF5
        filename = trim(filename) // '.h5'
#else
        filename = trim(filename) // '.binary'
#endif

        ! Reopen statepoint file in parallel
        call sp % file_open(filename, 'w', serial = .false.)

      end if

      ! Write out source
      call sp % write_source_bank()

      ! Close file
      call sp % file_close()

    end if

    ! Also check to write source separately in overwritten file
    if (source_latest) then

      ! Set filename
      filename = trim(path_output) // 'source'
#ifdef HDF5
      filename = trim(filename) // '.h5'
#else
      filename = trim(filename) // '.binary'
#endif

      ! Write message for new file creation
      message = "Creating source file " // trim(filename) // "..."
      call write_message(1)

      ! Always create this file because it will be overwritten
      call sp % file_create(filename, serial = .false.)

      ! Write file type
      call sp % write_data(FILETYPE_SOURCE, "filetype")

      ! Write out source
      call sp % write_source_bank()

      ! Close file
      call sp % file_close()

    end if

  end subroutine write_source_point

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
#ifdef MPI
    real(8) :: dummy  ! temporary receive buffer for non-root reduces
#endif
    type(TallyObject), pointer :: t => null()
    type(TallyResult), allocatable :: tallyresult_temp(:,:)

    ! ==========================================================================
    ! COLLECT AND WRITE GLOBAL TALLIES

    if (master) then
      ! Write number of realizations
      call sp % write_data(n_realizations, "n_realizations")

      ! Write number of global tallies
      call sp % write_data(N_GLOBAL_TALLIES, "n_global_tallies")
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
      call sp % write_tally_result(tallyresult_temp, "global_tallies", &
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
        call sp % write_data(1, "tallies_present", group="tallies")
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
          call sp % write_tally_result(t % results, "results", &
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
        call sp % write_data(0, "tallies_present", group="tallies")
      end if
    end if

  end subroutine write_tally_results_nr

!===============================================================================
! LOAD_STATE_POINT
!===============================================================================

  subroutine load_state_point()

    character(MAX_FILE_LEN) :: path_temp
    character(19)           :: current_time
    integer                 :: i
    integer                 :: j
    integer                 :: length(4)
    integer                 :: int_array(3)
    integer, allocatable    :: temp_array(:)
    logical                 :: source_present
    real(8)                 :: real_array(3) 
    type(TallyObject), pointer :: t => null()

    ! Write message
    message = "Loading state point " // trim(path_state_point) // "..."
    call write_message(1)

    ! Open file for reading
    call sp % file_open(path_state_point, 'r', serial = .false.)

    ! Read filetype
    call sp % read_data(int_array(1), "filetype")

    ! Read revision number for state point file and make sure it matches with
    ! current version
    call sp % read_data(int_array(1), "revision")
    if (int_array(1) /= REVISION_STATEPOINT) then
      message = "State point version does not match current version " &
                // "in OpenMC."
      call fatal_error()
    end if

    ! Read OpenMC version
    call sp % read_data(int_array(1), "version_major")
    call sp % read_data(int_array(2), "version_minor")
    call sp % read_data(int_array(3), "version_release")
    if (int_array(1) /= VERSION_MAJOR .or. int_array(2) /= VERSION_MINOR &
        .or. int_array(3) /= VERSION_RELEASE) then
      message = "State point file was created with a different version " &
                // "of OpenMC."
      call warning()
    end if

    ! Read date and time
    call sp % read_data(current_time, "date_and_time")

    ! Read path to input
    call sp % read_data(path_temp, "path")

    ! Read and overwrite random number seed
    call sp % read_data(seed, "seed")

    ! Read and overwrite run information except number of batches
    call sp % read_data(run_mode, "run_mode")
    call sp % read_data(n_particles, "n_particles")
    call sp % read_data(int_array(1), "n_batches")

    ! Take maximum of statepoint n_batches and input n_batches
    n_batches = max(n_batches, int_array(1))

    ! Read batch number to restart at
    call sp % read_data(restart_batch, "current_batch")

    ! Read information specific to eigenvalue run
    if (run_mode == MODE_EIGENVALUE) then
      call sp % read_data(int_array(1), "n_inactive")
      call sp % read_data(gen_per_batch, "gen_per_batch")
      call sp % read_data(k_generation, "k_generation", &
           length=restart_batch*gen_per_batch)
      call sp % read_data(entropy, "entropy", length=restart_batch*gen_per_batch)
      call sp % read_data(k_col_abs, "k_col_abs")
      call sp % read_data(k_col_tra, "k_col_tra")
      call sp % read_data(k_abs_tra, "k_abs_tra")
      call sp % read_data(real_array(1:2), "k_combined", length=2)

      ! Take maximum of statepoint n_inactive and input n_inactive
      n_inactive = max(n_inactive, int_array(1))

      ! Read in to see if CMFD was on
      call sp % read_data(int_array(1), "cmfd_on")

      ! Write out CMFD info
      if (int_array(1) == 1) then
        call sp % read_data(cmfd % indices, "indices", length=4, group="cmfd")
        call sp % read_data(cmfd % k_cmfd, "k_cmfd", length=restart_batch, &
             group="cmfd")
        length = cmfd % indices([4,1,2,3])
        call sp % read_data(cmfd % cmfd_src, "cmfd_src", &
             length=length, group="cmfd")
        call sp % read_data(cmfd % entropy, "cmfd_entropy", &
                       length=restart_batch, group="cmfd")
        call sp % read_data(cmfd % balance, "cmfd_balance", &
             length=restart_batch, group="cmfd")
        call sp % read_data(cmfd % dom, "cmfd_dominance", &
             length = restart_batch, group="cmfd")
        call sp % read_data(cmfd % src_cmp, "cmfd_srccmp", &
             length = restart_batch, group="cmfd")
      end if
    end if

    ! Read number of meshes
    call sp % read_data(n_meshes, "n_meshes", group="tallies")

    ! Read and overwrite mesh information
    MESH_LOOP: do i = 1, n_meshes
      call sp % read_data(meshes(i) % id, "id", &
           group="tallies/mesh" // to_str(i))
      call sp % read_data(meshes(i) % type, "type", &
           group="tallies/mesh" // to_str(i))
      call sp % read_data(meshes(i) % n_dimension, "n_dimension", &
           group="tallies/mesh" // to_str(i))
      call sp % read_data(meshes(i) % dimension, "dimension", &
           group="tallies/mesh" // to_str(i), &
           length=meshes(i) % n_dimension)
      call sp % read_data(meshes(i) % lower_left, "lower_left", &
           group="tallies/mesh" // to_str(i), &
           length=meshes(i) % n_dimension)
      call sp % read_data(meshes(i) % upper_right, "upper_right", &
           group="tallies/mesh" // to_str(i), &
           length=meshes(i) % n_dimension)
      call sp % read_data(meshes(i) % width, "width", &
           group="tallies/mesh" // to_str(i), &
           length=meshes(i) % n_dimension)
    end do MESH_LOOP

    ! Read and overwrite number of tallies
    call sp % read_data(n_tallies, "n_tallies", group="tallies")

    ! Read in tally metadata
    TALLY_METADATA: do i = 1, n_tallies

      ! Get pointer to tally
      t => tallies(i)

      ! Read tally id
      call sp % read_data(t % id, "id", group="tallies/tally" // to_str(i))

      ! Read number of realizations
      call sp % read_data(t % n_realizations, "n_realizations", &
           group="tallies/tally" // to_str(i))

      ! Read size of tally results
      call sp % read_data(int_array(1), "total_score_bins", &
           group="tallies/tally" // to_str(i))
      call sp % read_data(int_array(2), "total_filter_bins", &
           group="tallies/tally" // to_str(i))

      ! Check size of tally results array
      if (int_array(1) /= t % total_score_bins .and. &
          int_array(2) /= t % total_filter_bins) then
        message = "Input file tally structure is different from restart."
        call fatal_error()
      end if

      ! Read number of filters
      call sp % read_data(t % n_filters, "n_filters", &
           group="tallies/tally" // to_str(i))

      ! Read filter information
      FILTER_LOOP: do j = 1, t % n_filters

        ! Read type of filter
        call sp % read_data(t % filters(j) % type, "type", &
             group="tallies/tally" // trim(to_str(i)) // "/filter" // to_str(j))

        ! Read number of bins for this filter
        call sp % read_data(t % filters(j) % n_bins, "n_bins", &
             group="tallies/tally" // trim(to_str(i)) // "/filter" // to_str(j))

        ! Read bins
        if (t % filters(j) % type == FILTER_ENERGYIN .or. &
            t % filters(j) % type == FILTER_ENERGYOUT) then
          call sp % read_data(t % filters(j) % real_bins, "bins", &
               group="tallies/tally" // trim(to_str(i)) // "/filter" // to_str(j), &
               length=size(t % filters(j) % real_bins))
        else
          call sp % read_data(t % filters(j) % int_bins, "bins", &
               group="tallies/tally" // trim(to_str(i)) // "/filter" // to_str(j), &
               length=size(t % filters(j) % int_bins))
        end if

      end do FILTER_LOOP

      ! Read number of nuclide bins
      call sp % read_data(t % n_nuclide_bins, "n_nuclide_bins", &
           group="tallies/tally" // to_str(i))

      ! Set up nuclide bin array and then write
      allocate(temp_array(t % n_nuclide_bins))
      call sp % read_data(temp_array, "nuclide_bins", &
           group="tallies/tally" // to_str(i), length=t % n_nuclide_bins)
      NUCLIDE_LOOP: do j = 1, t % n_nuclide_bins
        if (temp_array(j) > 0) then
          nuclides(t % nuclide_bins(j)) % zaid = temp_array(j)
        else
          t % nuclide_bins(j) = temp_array(j)
        end if
      end do NUCLIDE_LOOP
      deallocate(temp_array)

      ! Write number of score bins, score bins, and scatt order
      call sp % read_data(t % n_score_bins, "n_score_bins", &
           group="tallies/tally" // to_str(i))
      call sp % read_data(t % score_bins, "score_bins", &
           group="tallies/tally" // to_str(i), length=t % n_score_bins)
      call sp % read_data(t % scatt_order, "scatt_order", &
           group="tallies/tally" // to_str(i), length=t % n_score_bins)

      ! Write number of user score bins
      call sp % read_data(t % n_user_score_bins, "n_user_score_bins", &
           group="tallies/tally" // to_str(i))

    end do TALLY_METADATA

    ! Check for source in statepoint if needed
    call sp % read_data(int_array(1), "source_present")
    if (int_array(1) == 1) then
      source_present = .true.
    else
      source_present = .false.
    end if

    ! Check to make sure source bank is present
    if (path_source_point == path_state_point .and. .not. source_present) then
      message = "Source bank must be contained in statepoint restart file"
      call fatal_error()
    end if 

    ! Read tallies to master
    if (master) then

      ! Read number of realizations for global tallies
      call sp % read_data(n_realizations, "n_realizations", collect=.false.)

      ! Read number of global tallies
      call sp % read_data(int_array(1), "n_global_tallies", collect=.false.)
      if (int_array(1) /= N_GLOBAL_TALLIES) then
        message = "Number of global tallies does not match in state point."
        call fatal_error()
      end if

      ! Read global tally data
      call sp % read_tally_result(global_tallies, "global_tallies", &
           n1=N_GLOBAL_TALLIES, n2=1)

      ! Check if tally results are present
      call sp % read_data(int_array(1), "tallies_present", group="tallies", collect=.false.)

      ! Read in sum and sum squared
      if (int_array(1) == 1) then
        TALLY_RESULTS: do i = 1, n_tallies

          ! Set pointer to tally
          t => tallies(i)

          ! Read sum and sum_sq for each bin
          call sp % read_tally_result(t % results, "results", &
               group="tallies/tally" // to_str(i), &
               n1=size(t % results, 1), n2=size(t % results, 2))

        end do TALLY_RESULTS
      end if
    end if

    ! Read source if in eigenvalue mode 
    if (run_mode == MODE_EIGENVALUE) then

      ! Check if source was written out separately
      if (.not. source_present) then

        ! Close statepoint file 
        call sp % file_close()

        ! Write message
        message = "Loading source file " // trim(path_source_point) // "..."
        call write_message(1)

        ! Open source file 
        call sp % file_open(path_source_point, 'r', serial = .false.)

        ! Read file type
        call sp % read_data(int_array(1), "filetype")

      end if

      ! Write out source
      call sp % read_source_bank()

    end if

    ! Close file
    call sp % file_close()

  end subroutine load_state_point

  subroutine read_source
! TODO write this routine
! TODO what if n_particles does not match source bank
  end subroutine read_source
end module state_point
