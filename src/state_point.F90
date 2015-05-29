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
  use string,             only: to_str, zero_padded, count_digits
  use output_interface
  use tally_header,       only: TallyObject
  use mesh_header,        only: StructuredMesh
  use dict_header,        only: ElemKeyValueII, ElemKeyValueCI

#ifdef MPI
  use mpi
#endif

  implicit none

  type(BinaryOutput)        :: sp      ! Statepoint/source output file

contains

!===============================================================================
! WRITE_STATE_POINT
!===============================================================================

  subroutine write_state_point()

    character(MAX_FILE_LEN)       :: filename
    integer                       :: i, j, k
    integer, allocatable          :: id_array(:)
    integer, allocatable          :: key_array(:)
    type(StructuredMesh), pointer :: mesh
    type(TallyObject), pointer    :: tally
    type(ElemKeyValueII), pointer :: current
    type(ElemKeyValueII), pointer :: next
    character(8)                  :: moment_name  ! name of moment (e.g, P3)
    integer                       :: n_order      ! loop index for moment orders
    integer                       :: nm_order     ! loop index for Ynm moment orders

    ! Set filename for state point
    filename = trim(path_output) // 'statepoint.' // &
        & zero_padded(current_batch, count_digits(n_max_batches))

    ! Append appropriate extension
#ifdef HDF5
    filename = trim(filename) // '.h5'
#else
    filename = trim(filename) // '.binary'
#endif

    ! Write message
    call write_message("Creating state point " // trim(filename) // "...", 1)

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

      ! Indicate whether source bank is stored in statepoint
      if (source_separate) then
        call sp % write_data(0, "source_present")
      else
        call sp % write_data(1, "source_present")
      end if

      ! Write out information for eigenvalue run
      if (run_mode == MODE_EIGENVALUE) then
        call sp % write_data(n_inactive, "n_inactive")
        call sp % write_data(gen_per_batch, "gen_per_batch")
        call sp % write_data(k_generation, "k_generation", &
             length=current_batch*gen_per_batch)
        call sp % write_data(entropy, "entropy", &
             length=current_batch*gen_per_batch)
        call sp % write_data(k_col_abs, "k_col_abs")
        call sp % write_data(k_col_tra, "k_col_tra")
        call sp % write_data(k_abs_tra, "k_abs_tra")
        call sp % write_data(k_combined, "k_combined", length=2)

        ! Write out CMFD info
        if (cmfd_on) then
#ifdef HDF5
          call sp % open_group("cmfd")
          call sp % close_group()
#endif
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

#ifdef HDF5
      call sp % open_group("tallies")
      call sp % close_group()
#endif

      ! Write number of meshes
      call sp % write_data(n_meshes, "n_meshes", group="tallies/meshes")

      if (n_meshes > 0) then

        ! Print list of mesh IDs
        current => mesh_dict % keys()

        allocate(id_array(n_meshes))
        allocate(key_array(n_meshes))
        i = 1

        do while (associated(current))
          key_array(i) = current % key
          id_array(i) = current % value

          ! Move to next mesh
          next => current % next
          deallocate(current)
          current => next
          i = i + 1
        end do

        call sp % write_data(id_array, "ids", &
             group="tallies/meshes", length=n_meshes)
        call sp % write_data(key_array, "keys", &
             group="tallies/meshes", length=n_meshes)

        deallocate(key_array)

        ! Write information for meshes
        MESH_LOOP: do i = 1, n_meshes

          mesh => meshes(id_array(i))

          call sp % write_data(mesh % id, "id", &
               group="tallies/meshes/mesh " // trim(to_str(mesh % id)))
          call sp % write_data(mesh % type, "type", &
               group="tallies/meshes/mesh " // trim(to_str(mesh % id)))
          call sp % write_data(mesh % n_dimension, "n_dimension", &
               group="tallies/meshes/mesh " // trim(to_str(mesh % id)))
          call sp % write_data(mesh % dimension, "dimension", &
               group="tallies/meshes/mesh " // trim(to_str(mesh % id)), &
               length=mesh % n_dimension)
          call sp % write_data(mesh % lower_left, "lower_left", &
               group="tallies/meshes/mesh " // trim(to_str(mesh % id)), &
               length=mesh % n_dimension)
          call sp % write_data(mesh % upper_right, "upper_right", &
               group="tallies/meshes/mesh " // trim(to_str(mesh % id)), &
               length=mesh % n_dimension)
          call sp % write_data(mesh % width, "width", &
               group="tallies/meshes/mesh " // trim(to_str(mesh % id)), &
               length=mesh % n_dimension)
        end do MESH_LOOP

        deallocate(id_array)

      end if

      ! Write number of tallies
      call sp % write_data(n_tallies, "n_tallies", group="tallies")

      if (n_tallies > 0) then

        ! Print list of tally IDs
        allocate(id_array(n_tallies))
        allocate(key_array(n_tallies))

        ! Write all tally information except results
        do i = 1, n_tallies
          tally => tallies(i)
          key_array(i) = tally % id
          id_array(i) = i
        end do

        call sp % write_data(id_array, "ids", &
             group="tallies", length=n_tallies)
        call sp % write_data(key_array, "keys", &
             group="tallies", length=n_tallies)

        deallocate(key_array)

        ! Write all tally information except results
        TALLY_METADATA: do i = 1, n_tallies

          ! Get pointer to tally
          tally => tallies(i)

          call sp % write_data(tally % estimator, "estimator", &
               group="tallies/tally " // trim(to_str(tally % id)))
          call sp % write_data(tally % n_realizations, "n_realizations", &
               group="tallies/tally " // trim(to_str(tally % id)))
          call sp % write_data(tally % n_filters, "n_filters", &
               group="tallies/tally " // trim(to_str(tally % id)))

          ! Write filter information
          FILTER_LOOP: do j = 1, tally % n_filters

            call sp % write_data(tally % filters(j) % type, "type", &
                 group="tallies/tally " // trim(to_str(tally % id)) // &
                 "/filter " // to_str(j))
            call sp % write_data(tally % filters(j) % offset, "offset", &
                 group="tallies/tally " // trim(to_str(tally % id)) // &
                 "/filter " // to_str(j))
            call sp % write_data(tally % filters(j) % n_bins, "n_bins", &
                 group="tallies/tally " // trim(to_str(tally % id)) // &
                 "/filter " // to_str(j))
            if (tally % filters(j) % type == FILTER_ENERGYIN .or. &
                tally % filters(j) % type == FILTER_ENERGYOUT) then
              call sp % write_data(tally % filters(j) % real_bins, "bins", &
                   group="tallies/tally " // trim(to_str(tally % id)) // &
                   "/filter " // to_str(j), &
                   length=size(tally % filters(j) % real_bins))
            else
              call sp % write_data(tally % filters(j) % int_bins, "bins", &
                   group="tallies/tally " // trim(to_str(tally % id)) // &
                   "/filter " // to_str(j), &
                   length=size(tally % filters(j) % int_bins))
            end if

          end do FILTER_LOOP

          call sp % write_data(tally % n_nuclide_bins, "n_nuclides", &
               group="tallies/tally " // trim(to_str(tally % id)))

          ! Set up nuclide bin array and then write
          allocate(key_array(tally % n_nuclide_bins))
          NUCLIDE_LOOP: do j = 1, tally % n_nuclide_bins
            if (tally % nuclide_bins(j) > 0) then
              key_array(j) = nuclides(tally % nuclide_bins(j)) % zaid
            else
              key_array(j) = tally % nuclide_bins(j)
            end if
          end do NUCLIDE_LOOP
          call sp % write_data(key_array, "nuclides", &
               group="tallies/tally " // trim(to_str(tally % id)), &
               length=tally % n_nuclide_bins)
          deallocate(key_array)

          call sp % write_data(tally % n_score_bins, "n_score_bins", &
               group="tallies/tally " // trim(to_str(tally % id)))
          call sp % write_data(tally % score_bins, "score_bins", &
               group="tallies/tally " // trim(to_str(tally % id)), &
               length=tally % n_score_bins)
          call sp % write_data(tally % n_user_score_bins, "n_user_score_bins", &
               group="tallies/tally " // to_str(tally % id))

          ! Write explicit moment order strings for each score bin
          k = 1
          MOMENT_LOOP: do j = 1, tally % n_user_score_bins
            select case(tally % score_bins(k))
            case (SCORE_SCATTER_N, SCORE_NU_SCATTER_N)
              moment_name = 'P' // trim(to_str(tally % moment_order(k)))
              call sp % write_data(moment_name, "order" // trim(to_str(k)), &
                   group="tallies/tally " // trim(to_str(tally % id)) // &
                         "/moments")
              k = k + 1
            case (SCORE_SCATTER_PN, SCORE_NU_SCATTER_PN)
              do n_order = 0, tally % moment_order(k)
                moment_name = 'P' // trim(to_str(n_order))
                call sp % write_data(moment_name, "order" // trim(to_str(k)), &
                   group="tallies/tally " // trim(to_str(tally % id)) // &
                         "/moments")
                k = k + 1
              end do
            case (SCORE_SCATTER_YN, SCORE_NU_SCATTER_YN, SCORE_FLUX_YN, &
                  SCORE_TOTAL_YN)
              do n_order = 0, tally % moment_order(k)
                do nm_order = -n_order, n_order
                  moment_name = 'Y' // trim(to_str(n_order)) // ',' // &
                    trim(to_str(nm_order))
                  call sp % write_data(moment_name, "order" // &
                       trim(to_str(k)), &
                       group="tallies/tally " // trim(to_str(tally % id)) // &
                             "/moments")
                    k = k + 1
                end do
              end do
            case default
              moment_name = ''
              call sp % write_data(moment_name, "order" // trim(to_str(k)), &
                   group="tallies/tally " // trim(to_str(tally % id)) // &
                         "/moments")
              k = k + 1
            end select

          end do MOMENT_LOOP

        end do TALLY_METADATA

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
          tally => tallies(i)

          ! Write sum and sum_sq for each bin
          call sp % write_tally_result(tally % results, "results", &
               group="tallies/tally " // trim(to_str(tally % id)), &
               n1=size(tally % results, 1), n2=size(tally % results, 2))

        end do TALLY_RESULTS

      else

        ! Indicate tallies are off
        call sp % write_data(0, "tallies_present", group="tallies")

      end if

      ! Close the file for serial writing
      call sp % file_close()

    end if

    if (master .and. n_tallies > 0) then
      deallocate(id_array)
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
        filename = trim(path_output) // 'source.' // &
            & zero_padded(current_batch, count_digits(n_max_batches))

#ifdef HDF5
        filename = trim(filename) // '.h5'
#else
        filename = trim(filename) // '.binary'
#endif

        ! Write message for new file creation
        call write_message("Creating source file " // trim(filename) // "...", &
             &1)

        ! Create separate source file
        call sp % file_create(filename, serial = .false.)

        ! Write file type
        call sp % write_data(FILETYPE_SOURCE, "filetype")

      else

        ! Set filename for state point
        filename = trim(path_output) // 'statepoint.' // &
            & zero_padded(current_batch, count_digits(n_max_batches))
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
      call write_message("Creating source file " // trim(filename) // "...", 1)

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
    integer, allocatable       :: id_array(:)
    type(ElemKeyValueII), pointer :: current
    type(ElemKeyValueII), pointer :: next
    type(TallyObject), pointer :: tally
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
      if (current_batch == n_max_batches .or. satisfy_triggers) then
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

        ! Build list of tally IDs
        current => tally_dict % keys()
        allocate(id_array(n_tallies))
        i = 1

        do while (associated(current))
          id_array(i) = current % value
          ! Move to next tally
          next => current % next
          deallocate(current)
          current => next
          i = i + 1
        end do

      end if

      ! Write all tally results
      TALLY_RESULTS: do i = 1, n_tallies

        tally => tallies(i)

        ! Determine size of tally results array
        m = size(tally % results, 1)
        n = size(tally % results, 2)
        n_bins = m*n*2

        ! Allocate array for storing sums and sums of squares, but
        ! contiguously in memory for each
        allocate(tally_temp(2,m,n))
        tally_temp(1,:,:) = tally % results(:,:) % sum
        tally_temp(2,:,:) = tally % results(:,:) % sum_sq

        if (master) then
          ! The MPI_IN_PLACE specifier allows the master to copy values into
          ! a receive buffer without having a temporary variable
#ifdef MPI
          call MPI_REDUCE(MPI_IN_PLACE, tally_temp, n_bins, MPI_REAL8, &
               MPI_SUM, 0, MPI_COMM_WORLD, mpi_err)
#endif

          ! At the end of the simulation, store the results back in the
          ! regular TallyResults array
          if (current_batch == n_max_batches .or. satisfy_triggers) then
            tally % results(:,:) % sum = tally_temp(1,:,:)
            tally % results(:,:) % sum_sq = tally_temp(2,:,:)
          end if

         ! Put in temporary tally result
         allocate(tallyresult_temp(m,n))
         tallyresult_temp(:,:) % sum    = tally_temp(1,:,:)
         tallyresult_temp(:,:) % sum_sq = tally_temp(2,:,:)

         ! Write reduced tally results to file
          call sp % write_tally_result(tally % results, "results", &
               group="tallies/tally " // trim(to_str(tally % id)), n1=m, n2=n)

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

      deallocate(id_array)

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

    character(MAX_FILE_LEN)    :: path_temp
    character(19)              :: current_time
    integer                    :: i, j, k
    integer                    :: length(4)
    integer                    :: int_array(3)
    integer, allocatable       :: id_array(:)
    integer, allocatable       :: key_array(:)
    integer                    :: curr_key
    integer, allocatable       :: temp_array(:)
    logical                    :: source_present
    real(8)                    :: real_array(3)
    type(StructuredMesh), pointer :: mesh
    type(TallyObject), pointer :: tally
    integer                    :: n_order      ! loop index for moment orders
    integer                    :: nm_order     ! loop index for Ynm moment orders
    character(8)               :: moment_name  ! name of moment (e.g, P3, Y-1,1)

    ! Write message
    call write_message("Loading state point " // trim(path_state_point) &
         &// "...", 1)

    ! Open file for reading
    call sp % file_open(path_state_point, 'r', serial = .false.)

    ! Read filetype
    call sp % read_data(int_array(1), "filetype")

    ! Read revision number for state point file and make sure it matches with
    ! current version
    call sp % read_data(int_array(1), "revision")
    if (int_array(1) /= REVISION_STATEPOINT) then
      call fatal_error("State point version does not match current version &
           &in OpenMC.")
    end if

    ! Read OpenMC version
    call sp % read_data(int_array(1), "version_major")
    call sp % read_data(int_array(2), "version_minor")
    call sp % read_data(int_array(3), "version_release")
    if (int_array(1) /= VERSION_MAJOR .or. int_array(2) /= VERSION_MINOR &
        .or. int_array(3) /= VERSION_RELEASE) then
      if (master) call warning("State point file was created with a different &
           &version of OpenMC.")
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

    ! Check for source in statepoint if needed
    call sp % read_data(int_array(1), "source_present")
    if (int_array(1) == 1) then
      source_present = .true.
    else
      source_present = .false.
    end if

    if (restart_batch > n_batches) then
      call fatal_error("The number batches specified in settings.xml is fewer &
           & than the number of batches in the given statepoint file.")
    end if

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

      ! Read in CMFD info
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
    call sp % read_data(n_meshes, "n_meshes", group="tallies/meshes")

    if (n_meshes > 0) then

      ! Read list of mesh keys-> IDs
      allocate(id_array(n_meshes))
      allocate(key_array(n_meshes))

      call sp % read_data(id_array, "ids", &
           group="tallies/meshes", length=n_meshes)
      call sp % read_data(key_array, "keys", &
           group="tallies/meshes", length=n_meshes)

      ! Read and overwrite mesh information
      MESH_LOOP: do i = 1, n_meshes

        mesh => meshes(id_array(i))
        curr_key = key_array(id_array(i))

        call sp % read_data(mesh % id, "id", &
             group="tallies/meshes/mesh " // trim(to_str(curr_key)))
        call sp % read_data(mesh % type, "type", &
             group="tallies/meshes/mesh " // trim(to_str(curr_key)))
        call sp % read_data(mesh % n_dimension, "n_dimension", &
             group="tallies/meshes/mesh " // trim(to_str(meshes(i) % id)))
        call sp % read_data(mesh % dimension, "dimension", &
             group="tallies/meshes/mesh " // trim(to_str(curr_key)), &
             length=mesh % n_dimension)
        call sp % read_data(mesh % lower_left, "lower_left", &
             group="tallies/meshes/mesh " // trim(to_str(curr_key)), &
             length=mesh % n_dimension)
        call sp % read_data(mesh % upper_right, "upper_right", &
             group="tallies/meshes/mesh " // trim(to_str(curr_key)), &
             length=mesh % n_dimension)
        call sp % read_data(mesh % width, "width", &
             group="tallies/meshes/mesh " // trim(to_str(curr_key)), &
             length=meshes(i) % n_dimension)

      end do MESH_LOOP

      deallocate(id_array)
      deallocate(key_array)

    end if

    ! Read and overwrite number of tallies
    call sp % read_data(n_tallies, "n_tallies", group="tallies")

    ! Read list of tally keys-> IDs
    allocate(id_array(n_tallies))
    allocate(key_array(n_tallies))

    call sp % read_data(id_array, "ids", group="tallies", length=n_tallies)
    call sp % read_data(key_array, "keys", group="tallies", length=n_tallies)

    ! Read in tally metadata
    TALLY_METADATA: do i = 1, n_tallies

      ! Get pointer to tally
      tally => tallies(i)
      curr_key = key_array(id_array(i))

      call sp % read_data(tally % estimator, "estimator", &
           group="tallies/tally " // trim(to_str(curr_key)))
      call sp % read_data(tally % n_realizations, "n_realizations", &
           group="tallies/tally " // trim(to_str(curr_key)))
      call sp % read_data(tally % n_filters, "n_filters", &
           group="tallies/tally " // trim(to_str(curr_key)))

      FILTER_LOOP: do j = 1, tally % n_filters
        call sp % read_data(tally % filters(j) % type, "type", &
             group="tallies/tally " // trim(to_str(curr_key)) // &
             "/filter " // to_str(j))
        call sp % read_data(tally % filters(j) % offset, "offset", &
              group="tallies/tally " // trim(to_str(curr_key)) // &
               "/filter " // to_str(j))
        call sp % read_data(tally % filters(j) % n_bins, "n_bins", &
             group="tallies/tally " // trim(to_str(curr_key)) // &
             "/filter " // to_str(j))
        if (tally % filters(j) % type == FILTER_ENERGYIN .or. &
            tally % filters(j) % type == FILTER_ENERGYOUT) then
          call sp % read_data(tally % filters(j) % real_bins, "bins", &
               group="tallies/tally " // trim(to_str(curr_key)) // &
               "/filter " // to_str(j), &
               length=size(tally % filters(j) % real_bins))
        else
          call sp % read_data(tally % filters(j) % int_bins, "bins", &
               group="tallies/tally " // trim(to_str(curr_key)) // &
               "/filter " // to_str(j), &
               length=size(tally % filters(j) % int_bins))
        end if

      end do FILTER_LOOP

      call sp % read_data(tally % n_nuclide_bins, "n_nuclides", &
           group="tallies/tally " // trim(to_str(curr_key)))

      ! Set up nuclide bin array and then read
      allocate(temp_array(tally % n_nuclide_bins))
      call sp % read_data(temp_array, "nuclides", &
           group="tallies/tally " // trim(to_str(curr_key)), &
           length=tally % n_nuclide_bins)

      NUCLIDE_LOOP: do j = 1, tally % n_nuclide_bins
        if (temp_array(j) > 0) then
          tally % nuclide_bins(j) = temp_array(j)
        else
          tally % nuclide_bins(j) = temp_array(j)
        end if
      end do NUCLIDE_LOOP

      deallocate(temp_array)

      ! Write number of score bins, score bins, user score bins
      call sp % read_data(tally % n_score_bins, "n_score_bins", &
           group="tallies/tally " // trim(to_str(curr_key)))
      call sp % read_data(tally % score_bins, "score_bins", &
           group="tallies/tally " // trim(to_str(curr_key)), &
           length=tally % n_score_bins)
      call sp % read_data(tally % n_user_score_bins, "n_user_score_bins", &
           group="tallies/tally " // trim(to_str(curr_key)))

      ! Read explicit moment order strings for each score bin
      k = 1
      MOMENT_LOOP: do j = 1, tally % n_user_score_bins
        select case(tally % score_bins(k))
        case (SCORE_SCATTER_N, SCORE_NU_SCATTER_N)
          call sp % read_data(moment_name, "order" // trim(to_str(k)), &
               group="tallies/tally " // trim(to_str(curr_key)) // "/moments")
          k = k + 1
        case (SCORE_SCATTER_PN, SCORE_NU_SCATTER_PN)
          do n_order = 0, tally % moment_order(k)
            call sp % read_data(moment_name, "order" // trim(to_str(k)), &
                 group="tallies/tally " // trim(to_str(curr_key)) // "/moments")
            k = k + 1
          end do
        case (SCORE_SCATTER_YN, SCORE_NU_SCATTER_YN, SCORE_FLUX_YN, &
              SCORE_TOTAL_YN)
          do n_order = 0, tally % moment_order(k)
            do nm_order = -n_order, n_order
              call sp % read_data(moment_name, "order" // trim(to_str(k)), &
                   group="tallies/tally " // trim(to_str(curr_key)) // &
                         "/moments")
              k = k + 1
            end do
          end do
        case default
          call sp % read_data(moment_name, "order" // trim(to_str(k)), &
               group="tallies/tally " // trim(to_str(curr_key)) // "/moments")
          k = k + 1
        end select

      end do MOMENT_LOOP

    end do TALLY_METADATA

    ! Check to make sure source bank is present
    if (path_source_point == path_state_point .and. .not. source_present) then
      call fatal_error("Source bank must be contained in statepoint restart &
           &file")
    end if

    ! Read tallies to master
    if (master) then

      ! Read number of realizations for global tallies
      call sp % read_data(n_realizations, "n_realizations", collect=.false.)

      ! Read number of global tallies
      call sp % read_data(int_array(1), "n_global_tallies", collect=.false.)
      if (int_array(1) /= N_GLOBAL_TALLIES) then
        call fatal_error("Number of global tallies does not match in state &
             &point.")
      end if

      ! Read global tally data
      call sp % read_tally_result(global_tallies, "global_tallies", &
           n1=N_GLOBAL_TALLIES, n2=1)

      ! Check if tally results are present
      call sp % read_data(int_array(1), "tallies_present", &
           group="tallies", collect=.false.)

      ! Read in sum and sum squared
      if (int_array(1) == 1) then
        TALLY_RESULTS: do i = 1, n_tallies

          ! Set pointer to tally
          tally => tallies(i)
          curr_key = key_array(id_array(i))

          ! Read sum and sum_sq for each bin
          call sp % read_tally_result(tally % results, "results", &
               group="tallies/tally " // trim(to_str(curr_key)), &
               n1=size(tally % results, 1), n2=size(tally % results, 2))

        end do TALLY_RESULTS

      end if
    end if

    deallocate(id_array)
    deallocate(key_array)

    ! Read source if in eigenvalue mode
    if (run_mode == MODE_EIGENVALUE) then

      ! Check if source was written out separately
      if (.not. source_present) then

        ! Close statepoint file
        call sp % file_close()

        ! Write message
        call write_message("Loading source file " // trim(path_source_point) &
             &// "...", 1)

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
