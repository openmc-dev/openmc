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
  use hdf5_interface
  use output,             only: write_message, time_stamp
  use string,             only: to_str, zero_padded, count_digits
  use tally_header,       only: TallyObject
  use mesh_header,        only: StructuredMesh
  use dict_header,        only: ElemKeyValueII, ElemKeyValueCI

#ifdef MPI
  use message_passing
#endif

  use hdf5

  use, intrinsic :: ISO_C_BINDING, only: c_loc, c_ptr

  implicit none

contains

!===============================================================================
! WRITE_STATE_POINT
!===============================================================================

  subroutine write_state_point()

    integer                       :: i, j, k
    integer                       :: n_order      ! loop index for moment orders
    integer                       :: nm_order     ! loop index for Ynm moment orders
    integer, allocatable          :: id_array(:)
    integer, allocatable          :: key_array(:)
    integer(HID_T) :: file_id
    integer(HID_T) :: cmfd_group
    integer(HID_T) :: tallies_group, tally_group
    integer(HID_T) :: meshes_group, mesh_group
    integer(HID_T) :: filter_group
    character(8), allocatable :: moment_names(:) ! names of moments (e.g, P3)
    character(MAX_FILE_LEN)       :: filename
    type(StructuredMesh), pointer :: meshp
    type(TallyObject), pointer    :: tally
    type(ElemKeyValueII), pointer :: current
    type(ElemKeyValueII), pointer :: next

    ! Set filename for state point
    filename = trim(path_output) // 'statepoint.' // &
         & zero_padded(current_batch, count_digits(n_max_batches))
    filename = trim(filename) // '.h5'

    ! Write message
    call write_message("Creating state point " // trim(filename) // "...", 1)

    if (master) then
      ! Create statepoint file
      file_id = file_create(filename)

      ! Write file type
      call write_dataset(file_id, "filetype", 'statepoint')

      ! Write revision number for state point file
      call write_dataset(file_id, "revision", REVISION_STATEPOINT)

      ! Write OpenMC version
      call write_dataset(file_id, "version_major", VERSION_MAJOR)
      call write_dataset(file_id, "version_minor", VERSION_MINOR)
      call write_dataset(file_id, "version_release", VERSION_RELEASE)

      ! Write current date and time
      call write_dataset(file_id, "date_and_time", time_stamp())

      ! Write path to input
      call write_dataset(file_id, "path", path_input)

      ! Write out random number seed
      call write_dataset(file_id, "seed", seed)

      ! Write run information
      call write_dataset(file_id, "run_mode", run_mode)
      call write_dataset(file_id, "n_particles", n_particles)
      call write_dataset(file_id, "n_batches", n_batches)

      ! Write out current batch number
      call write_dataset(file_id, "current_batch", current_batch)

      ! Indicate whether source bank is stored in statepoint
      if (source_separate) then
        call write_dataset(file_id, "source_present", 0)
      else
        call write_dataset(file_id, "source_present", 1)
      end if

      ! Write out information for eigenvalue run
      if (run_mode == MODE_EIGENVALUE) then
        call write_dataset(file_id, "n_inactive", n_inactive)
        call write_dataset(file_id, "gen_per_batch", gen_per_batch)
        call write_dataset(file_id, "k_generation", k_generation)
        call write_dataset(file_id, "entropy", entropy)
        call write_dataset(file_id, "k_col_abs", k_col_abs)
        call write_dataset(file_id, "k_col_tra", k_col_tra)
        call write_dataset(file_id, "k_abs_tra", k_abs_tra)
        call write_dataset(file_id, "k_combined", k_combined)

        ! Write out CMFD info
        if (cmfd_on) then
          call write_dataset(file_id, "cmfd_on", 1)

          cmfd_group = create_group(file_id, "cmfd")
          call write_dataset(cmfd_group, "indices", cmfd%indices)
          call write_dataset(cmfd_group, "k_cmfd", cmfd%k_cmfd)
          call write_dataset(cmfd_group, "cmfd_src", cmfd%cmfd_src)
          call write_dataset(cmfd_group, "cmfd_entropy", cmfd%entropy)
          call write_dataset(cmfd_group, "cmfd_balance", cmfd%balance)
          call write_dataset(cmfd_group, "cmfd_dominance", cmfd%dom)
          call write_dataset(cmfd_group, "cmfd_srccmp", cmfd%src_cmp)
          call close_group(cmfd_group)
        else
          call write_dataset(file_id, "cmfd_on", 0)
        end if
      end if

      tallies_group = create_group(file_id, "tallies")

      ! Write number of meshes
      meshes_group = create_group(tallies_group, "meshes")
      call write_dataset(meshes_group, "n_meshes", n_meshes)

      if (n_meshes > 0) then

        ! Print list of mesh IDs
        current => mesh_dict%keys()

        allocate(id_array(n_meshes))
        allocate(key_array(n_meshes))
        i = 1

        do while (associated(current))
          key_array(i) = current%key
          id_array(i) = current%value

          ! Move to next mesh
          next => current%next
          deallocate(current)
          current => next
          i = i + 1
        end do

        call write_dataset(meshes_group, "ids", id_array)
        call write_dataset(meshes_group, "keys", key_array)

        deallocate(key_array)

        ! Write information for meshes
        MESH_LOOP: do i = 1, n_meshes
          meshp => meshes(id_array(i))
          mesh_group = create_group(meshes_group, "mesh " // trim(to_str(meshp%id)))

          call write_dataset(mesh_group, "id", meshp%id)
          call write_dataset(mesh_group, "type", meshp%type)
          call write_dataset(mesh_group, "n_dimension", meshp%n_dimension)
          call write_dataset(mesh_group, "dimension", meshp%dimension)
          call write_dataset(mesh_group, "lower_left", meshp%lower_left)
          call write_dataset(mesh_group, "upper_right", meshp%upper_right)
          call write_dataset(mesh_group, "width", meshp%width)

          call close_group(mesh_group)
        end do MESH_LOOP

        deallocate(id_array)
      end if

      call close_group(meshes_group)

      ! Write number of tallies
      call write_dataset(tallies_group, "n_tallies", n_tallies)

      if (n_tallies > 0) then

        ! Print list of tally IDs
        allocate(id_array(n_tallies))
        allocate(key_array(n_tallies))

        ! Write all tally information except results
        do i = 1, n_tallies
          tally => tallies(i)
          key_array(i) = tally%id
          id_array(i) = i
        end do

        call write_dataset(tallies_group, "ids", id_array)
        call write_dataset(tallies_group, "keys", key_array)

        deallocate(key_array)

        ! Write all tally information except results
        TALLY_METADATA: do i = 1, n_tallies

          ! Get pointer to tally
          tally => tallies(i)
          tally_group = create_group(tallies_group, "tally " // &
               trim(to_str(tally%id)))

          call write_dataset(tally_group, "estimator", tally%estimator)
          call write_dataset(tally_group, "n_realizations", tally%n_realizations)
          call write_dataset(tally_group, "n_filters", tally%n_filters)

          ! Write filter information
          FILTER_LOOP: do j = 1, tally%n_filters
            filter_group = create_group(tally_group, "filter " // &
                 trim(to_str(j)))

            call write_dataset(filter_group, "type", tally%filters(j)%type)
            call write_dataset(filter_group, "offset", tally%filters(j)%offset)
            call write_dataset(filter_group, "n_bins", tally%filters(j)%n_bins)
            if (tally%filters(j)%type == FILTER_ENERGYIN .or. &
                 tally%filters(j)%type == FILTER_ENERGYOUT) then
              call write_dataset(filter_group, "bins", &
                   tally%filters(j)%real_bins)
            else
              call write_dataset(filter_group, "bins", &
                   tally%filters(j)%int_bins)
            end if

            call close_group(filter_group)
          end do FILTER_LOOP

          call write_dataset(tally_group, "n_nuclides", tally%n_nuclide_bins)

          ! Set up nuclide bin array and then write
          allocate(key_array(tally%n_nuclide_bins))
          NUCLIDE_LOOP: do j = 1, tally%n_nuclide_bins
            if (tally%nuclide_bins(j) > 0) then
              key_array(j) = nuclides(tally%nuclide_bins(j))%zaid
            else
              key_array(j) = tally%nuclide_bins(j)
            end if
          end do NUCLIDE_LOOP
          call write_dataset(tally_group, "nuclides", key_array)
          deallocate(key_array)

          call write_dataset(tally_group, "n_score_bins", tally%n_score_bins)
          call write_dataset(tally_group, "score_bins", tally%score_bins)
          call write_dataset(tally_group, "n_user_score_bins", tally%n_user_score_bins)

          ! Write explicit moment order strings for each score bin
          k = 1
          allocate(moment_names(tally%n_score_bins))
          MOMENT_LOOP: do j = 1, tally%n_user_score_bins
            select case(tally%score_bins(k))
            case (SCORE_SCATTER_N, SCORE_NU_SCATTER_N)
              moment_names(k) = 'P' // trim(to_str(tally%moment_order(k)))
              k = k + 1
            case (SCORE_SCATTER_PN, SCORE_NU_SCATTER_PN)
              do n_order = 0, tally%moment_order(k)
                moment_names(k) = 'P' // trim(to_str(n_order))
                k = k + 1
              end do
            case (SCORE_SCATTER_YN, SCORE_NU_SCATTER_YN, SCORE_FLUX_YN, &
                 SCORE_TOTAL_YN)
              do n_order = 0, tally%moment_order(k)
                do nm_order = -n_order, n_order
                  moment_names(k) = 'Y' // trim(to_str(n_order)) // ',' // &
                       trim(to_str(nm_order))
                  k = k + 1
                end do
              end do
            case default
              moment_names(k) = ''
              k = k + 1
            end select
          end do MOMENT_LOOP

          call write_dataset(tally_group, "moment_orders", moment_names)
          deallocate(moment_names)

          call close_group(tally_group)
        end do TALLY_METADATA

      end if

      call close_group(tallies_group)
    end if

    ! Check for the no-tally-reduction method
    if (.not. reduce_tallies) then
      ! If using the no-tally-reduction method, we need to collect tally
      ! results before writing them to the state point file.

      call write_tally_results_nr(file_id)

    elseif (master) then

      ! Write number of global realizations
      call write_dataset(file_id, "n_realizations", n_realizations)

      ! Write global tallies
      call write_dataset(file_id, "n_global_tallies", N_GLOBAL_TALLIES)
      call write_dataset(file_id, "global_tallies", global_tallies)

      ! Write tallies
      tallies_group = open_group(file_id, "tallies")
      if (tallies_on) then
        ! Indicate that tallies are on
        call write_dataset(tallies_group, "tallies_present", 1)

        ! Write all tally results
        TALLY_RESULTS: do i = 1, n_tallies
          ! Set point to current tally
          tally => tallies(i)

          ! Write sum and sum_sq for each bin
          tally_group = open_group(tallies_group, "tally " // to_str(tally%id))
          call write_dataset(tally_group, "results", tally%results)
          call close_group(tally_group)
        end do TALLY_RESULTS

      else
        ! Indicate tallies are off
        call write_dataset(tallies_group, "tallies_present", 0)
      end if

      call close_group(tallies_group)
      call file_close(file_id)
    end if

    if (master .and. n_tallies > 0) then
      deallocate(id_array)
    end if

  end subroutine write_state_point

!===============================================================================
! WRITE_SOURCE_POINT
!===============================================================================

  subroutine write_source_point()

    logical :: parallel
    integer(HID_T) :: file_id
    character(MAX_FILE_LEN) :: filename

    ! When using parallel HDF5, the file is written to collectively by all
    ! processes. With MPI-only, the file is opened and written by the master
    ! (note that the call to write_source_bank is by all processes since slave
    ! processes need to send source bank data to the master.
#ifdef PHDF5
    parallel = .true.
#else
    parallel = .false.
#endif

    ! Check to write out source for a specified batch
    if (sourcepoint_batch%contains(current_batch)) then
      if (source_separate) then
        filename = trim(path_output) // 'source.' // &
             & zero_padded(current_batch, count_digits(n_max_batches))
        filename = trim(filename) // '.h5'
        call write_message("Creating source file " // trim(filename) &
             // "...", 1)

        ! Create separate source file
        if (master .or. parallel) then
          file_id = file_create(filename, parallel=.true.)
          call write_dataset(file_id, "filetype", 'source')
        end if
      else
        filename = trim(path_output) // 'statepoint.' // &
             zero_padded(current_batch, count_digits(n_max_batches))
        filename = trim(filename) // '.h5'

        if (master .or. parallel) then
          file_id = file_open(filename, 'w', parallel=.true.)
        end if
      end if

      call write_source_bank(file_id)
      if (master .or. parallel) call file_close(file_id)
    end if

    ! Also check to write source separately in overwritten file
    if (source_latest) then
      filename = trim(path_output) // 'source' // '.h5'
      call write_message("Creating source file " // trim(filename) // "...", 1)
      if (master .or. parallel) then
        file_id = file_create(filename, parallel=.true.)
        call write_dataset(file_id, "filetype", 'source')
      end if

      call write_source_bank(file_id)

      if (master .or. parallel) call file_close(file_id)
    end if

  end subroutine write_source_point

!===============================================================================
! WRITE_TALLY_RESULTS_NR
!===============================================================================

  subroutine write_tally_results_nr(file_id)
    integer(HID_T), intent(in) :: file_id

    integer :: i      ! loop index
    integer :: n      ! number of filter bins
    integer :: m      ! number of score bins
    integer :: n_bins ! total number of bins
    integer(HID_T) :: tallies_group, tally_group
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
      call write_dataset(file_id, "n_realizations", n_realizations)

      ! Write number of global tallies
      call write_dataset(file_id, "n_global_tallies", N_GLOBAL_TALLIES)

      tallies_group = open_group(file_id, "tallies")
    end if

    ! Copy global tallies into temporary array for reducing
    n_bins = 2 * N_GLOBAL_TALLIES
    global_temp(1,:) = global_tallies(:)%sum
    global_temp(2,:) = global_tallies(:)%sum_sq

    if (master) then
      ! The MPI_IN_PLACE specifier allows the master to copy values into a
      ! receive buffer without having a temporary variable
#ifdef MPI
      call MPI_REDUCE(MPI_IN_PLACE, global_temp, n_bins, MPI_REAL8, MPI_SUM, &
           0, MPI_COMM_WORLD, mpi_err)
#endif

      ! Transfer values to value on master
      if (current_batch == n_max_batches .or. satisfy_triggers) then
        global_tallies(:)%sum    = global_temp(1,:)
        global_tallies(:)%sum_sq = global_temp(2,:)
      end if

      ! Put reduced value in temporary tally result
      allocate(tallyresult_temp(N_GLOBAL_TALLIES, 1))
      tallyresult_temp(:,1)%sum    = global_temp(1,:)
      tallyresult_temp(:,1)%sum_sq = global_temp(2,:)

      ! Write out global tallies sum and sum_sq
      call write_dataset(file_id, "global_tallies", tallyresult_temp)

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
        call write_dataset(tallies_group, "tallies_present", 1)

        ! Build list of tally IDs
        current => tally_dict%keys()
        allocate(id_array(n_tallies))
        i = 1

        do while (associated(current))
          id_array(i) = current%value
          ! Move to next tally
          next => current%next
          deallocate(current)
          current => next
          i = i + 1
        end do

      end if

      ! Write all tally results
      TALLY_RESULTS: do i = 1, n_tallies

        tally => tallies(i)

        ! Determine size of tally results array
        m = size(tally%results, 1)
        n = size(tally%results, 2)
        n_bins = m*n*2

        ! Allocate array for storing sums and sums of squares, but
        ! contiguously in memory for each
        allocate(tally_temp(2,m,n))
        tally_temp(1,:,:) = tally%results(:,:)%sum
        tally_temp(2,:,:) = tally%results(:,:)%sum_sq

        if (master) then
          tally_group = open_group(tallies_group, "tally " // &
               trim(to_str(tally%id)))

          ! The MPI_IN_PLACE specifier allows the master to copy values into
          ! a receive buffer without having a temporary variable
#ifdef MPI
          call MPI_REDUCE(MPI_IN_PLACE, tally_temp, n_bins, MPI_REAL8, &
               MPI_SUM, 0, MPI_COMM_WORLD, mpi_err)
#endif

          ! At the end of the simulation, store the results back in the
          ! regular TallyResults array
          if (current_batch == n_max_batches .or. satisfy_triggers) then
            tally%results(:,:)%sum = tally_temp(1,:,:)
            tally%results(:,:)%sum_sq = tally_temp(2,:,:)
          end if

          ! Put in temporary tally result
          allocate(tallyresult_temp(m,n))
          tallyresult_temp(:,:)%sum    = tally_temp(1,:,:)
          tallyresult_temp(:,:)%sum_sq = tally_temp(2,:,:)

          ! Write reduced tally results to file
          call write_dataset(tally_group, "results", tally%results)

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

        if (master) call close_group(tally_group)
      end do TALLY_RESULTS

      deallocate(id_array)

    else
      if (master) then
        ! Indicate that tallies are off
        call write_dataset(tallies_group, "tallies_present", 0)
      end if
    end if

    if (master) call close_group(tallies_group)

  end subroutine write_tally_results_nr

!===============================================================================
! LOAD_STATE_POINT
!===============================================================================

  subroutine load_state_point()

    integer                    :: i, j
    integer                    :: int_array(3)
    integer                    :: curr_key
    integer, allocatable       :: id_array(:)
    integer, allocatable       :: key_array(:)
    integer, allocatable       :: temp_array(:)
    integer(HID_T) :: file_id
    integer(HID_T) :: cmfd_group
    integer(HID_T) :: tallies_group, tally_group
    integer(HID_T) :: meshes_group, mesh_group
    integer(HID_T) :: filter_group
    real(8)                    :: real_array(3)
    logical                    :: source_present
    character(MAX_FILE_LEN)    :: path_temp
    character(19)              :: current_time
    type(StructuredMesh), pointer :: meshp
    type(TallyObject), pointer :: tally

    ! Write message
    call write_message("Loading state point " // trim(path_state_point) &
         // "...", 1)

    ! Open file for reading
    file_id = file_open(path_state_point, 'r', parallel=.true.)

    ! Read filetype
    call read_dataset(file_id, "filetype", int_array(1))

    ! Read revision number for state point file and make sure it matches with
    ! current version
    call read_dataset(file_id, "revision", int_array(1))
    if (int_array(1) /= REVISION_STATEPOINT) then
      call fatal_error("State point version does not match current version &
           &in OpenMC.")
    end if

    ! Read OpenMC version
    call read_dataset(file_id, "version_major", int_array(1))
    call read_dataset(file_id, "version_minor", int_array(2))
    call read_dataset(file_id, "version_release", int_array(3))
    if (int_array(1) /= VERSION_MAJOR .or. int_array(2) /= VERSION_MINOR &
         .or. int_array(3) /= VERSION_RELEASE) then
      if (master) call warning("State point file was created with a different &
           &version of OpenMC.")
    end if

    ! Read date and time
    call read_dataset(file_id, "date_and_time", current_time)

    ! Read path to input
    call read_dataset(file_id, "path", path_temp)

    ! Read and overwrite random number seed
    call read_dataset(file_id, "seed", seed)

    ! Read and overwrite run information except number of batches
    call read_dataset(file_id, "run_mode", run_mode)
    call read_dataset(file_id, "n_particles", n_particles)
    call read_dataset(file_id, "n_batches", int_array(1))

    ! Take maximum of statepoint n_batches and input n_batches
    n_batches = max(n_batches, int_array(1))

    ! Read batch number to restart at
    call read_dataset(file_id, "current_batch", restart_batch)

    ! Check for source in statepoint if needed
    call read_dataset(file_id, "source_present", int_array(1))
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
      call read_dataset(file_id, "n_inactive", int_array(1))
      call read_dataset(file_id, "gen_per_batch", gen_per_batch)
      call read_dataset(file_id, "k_generation", &
           k_generation(1:restart_batch*gen_per_batch))
      call read_dataset(file_id, "entropy", &
           entropy(1:restart_batch*gen_per_batch))
      call read_dataset(file_id, "k_col_abs", k_col_abs)
      call read_dataset(file_id, "k_col_tra", k_col_tra)
      call read_dataset(file_id, "k_abs_tra", k_abs_tra)
      call read_dataset(file_id, "k_combined", real_array(1:2))

      ! Take maximum of statepoint n_inactive and input n_inactive
      n_inactive = max(n_inactive, int_array(1))

      ! Read in to see if CMFD was on
      call read_dataset(file_id, "cmfd_on", int_array(1))

      ! Read in CMFD info
      if (int_array(1) == 1) then
        cmfd_group = open_group(file_id, "cmfd")
        call read_dataset(cmfd_group, "indices", cmfd%indices)
        call read_dataset(cmfd_group, "k_cmfd", cmfd%k_cmfd(1:restart_batch))
        call read_dataset(cmfd_group, "cmfd_src", cmfd%cmfd_src)
        call read_dataset(cmfd_group, "cmfd_entropy", &
             cmfd%entropy(1:restart_batch))
        call read_dataset(cmfd_group, "cmfd_balance", &
             cmfd%balance(1:restart_batch))
        call read_dataset(cmfd_group, "cmfd_dominance", &
             cmfd%dom(1:restart_batch))
        call read_dataset(cmfd_group, "cmfd_srccmp", &
             cmfd%src_cmp(1:restart_batch))
        call close_group(cmfd_group)
      end if
    end if

    ! Read number of meshes
    tallies_group = open_group(file_id, "tallies")
    meshes_group = open_group(tallies_group, "meshes")
    call read_dataset(meshes_group, "n_meshes", n_meshes)

    if (n_meshes > 0) then

      ! Read list of mesh keys-> IDs
      allocate(id_array(n_meshes))
      allocate(key_array(n_meshes))

      call read_dataset(meshes_group, "ids", id_array)
      call read_dataset(meshes_group, "keys", key_array)

      ! Read and overwrite mesh information
      MESH_LOOP: do i = 1, n_meshes

        meshp => meshes(id_array(i))
        curr_key = key_array(id_array(i))

        mesh_group = open_group(meshes_group, "mesh " // &
             trim(to_str(curr_key)))
        call read_dataset(mesh_group, "id", meshp%id)
        call read_dataset(mesh_group, "type", meshp%type)
        call read_dataset(mesh_group, "n_dimension", meshp%n_dimension)
        call read_dataset(mesh_group, "dimension", meshp%dimension)
        call read_dataset(mesh_group, "lower_left", meshp%lower_left)
        call read_dataset(mesh_group, "upper_right", meshp%upper_right)
        call read_dataset(mesh_group, "width", meshp%width)
        call close_group(mesh_group)
      end do MESH_LOOP

      deallocate(id_array)
      deallocate(key_array)

    end if

    call close_group(meshes_group)

    ! Read and overwrite number of tallies
    call read_dataset(tallies_group, "n_tallies", n_tallies)

    ! Read list of tally keys-> IDs
    allocate(id_array(n_tallies))
    allocate(key_array(n_tallies))

    call read_dataset(tallies_group, "ids", id_array)
    call read_dataset(tallies_group, "keys", key_array)

    ! Read in tally metadata
    TALLY_METADATA: do i = 1, n_tallies

      ! Get pointer to tally
      tally => tallies(i)
      curr_key = key_array(id_array(i))
      tally_group = open_group(tallies_group, "tally " // &
           trim(to_str(curr_key)))

      call read_dataset(tally_group, "estimator", tally%estimator)
      call read_dataset(tally_group, "n_realizations", tally%n_realizations)
      call read_dataset(tally_group, "n_filters", tally%n_filters)

      FILTER_LOOP: do j = 1, tally%n_filters
        filter_group = open_group(tally_group, "filter " // trim(to_str(j)))

        call read_dataset(filter_group, "type", tally%filters(j)%type)
        call read_dataset(filter_group, "offset", tally%filters(j)%offset)
        call read_dataset(filter_group, "n_bins", tally%filters(j)%n_bins)
        if (tally%filters(j)%type == FILTER_ENERGYIN .or. &
             tally%filters(j)%type == FILTER_ENERGYOUT) then
          call read_dataset(filter_group, "bins", tally%filters(j)%real_bins)
        else
          call read_dataset(filter_group, "bins", tally%filters(j)%int_bins)
        end if

        call close_group(filter_group)
      end do FILTER_LOOP

      call read_dataset(tally_group, "n_nuclides", tally%n_nuclide_bins)

      ! Set up nuclide bin array and then read
      allocate(temp_array(tally%n_nuclide_bins))
      call read_dataset(tally_group, "nuclides", temp_array)

      NUCLIDE_LOOP: do j = 1, tally%n_nuclide_bins
        if (temp_array(j) > 0) then
          tally%nuclide_bins(j) = temp_array(j)
        else
          tally%nuclide_bins(j) = temp_array(j)
        end if
      end do NUCLIDE_LOOP

      deallocate(temp_array)

      ! Write number of score bins, score bins, user score bins
      call read_dataset(tally_group, "n_score_bins", tally%n_score_bins)
      call read_dataset(tally_group, "score_bins", tally%score_bins)
      call read_dataset(tally_group, "n_user_score_bins", tally%n_user_score_bins)

      call close_group(tally_group)
    end do TALLY_METADATA

    ! Check to make sure source bank is present
    if (path_source_point == path_state_point .and. .not. source_present) then
      call fatal_error("Source bank must be contained in statepoint restart &
           &file")
    end if

    ! Read tallies to master
    if (master) then

      ! Read number of realizations for global tallies
      call read_dataset(file_id, "n_realizations", n_realizations, indep=.true.)

      ! Read number of global tallies
      call read_dataset(file_id, "n_global_tallies", int_array(1), indep=.false.)
      if (int_array(1) /= N_GLOBAL_TALLIES) then
        call fatal_error("Number of global tallies does not match in state &
             &point.")
      end if

      ! Read global tally data
      call read_dataset(file_id, "global_tallies", global_tallies)

      ! Check if tally results are present
      tallies_group = open_group(file_id, "tallies")
      call read_dataset(file_id, "tallies_present", int_array(1), indep=.true.)

      ! Read in sum and sum squared
      if (int_array(1) == 1) then
        TALLY_RESULTS: do i = 1, n_tallies

          ! Set pointer to tally
          tally => tallies(i)
          curr_key = key_array(id_array(i))

          ! Read sum and sum_sq for each bin
          tally_group = open_group(tallies_group, "tally " // &
               trim(to_str(curr_key)))
          call read_dataset(tally_group, "results", tally%results)
          call close_group(tally_group)
        end do TALLY_RESULTS
      end if

      call close_group(tallies_group)
    end if

    deallocate(id_array)
    deallocate(key_array)

    ! Read source if in eigenvalue mode
    if (run_mode == MODE_EIGENVALUE) then

      ! Check if source was written out separately
      if (.not. source_present) then

        ! Close statepoint file
        call file_close(file_id)

        ! Write message
        call write_message("Loading source file " // trim(path_source_point) &
             // "...", 1)

        ! Open source file
        file_id = file_open(path_source_point, 'r', parallel=.true.)

        ! Read file type
        call read_dataset(file_id, "filetype", int_array(1))

      end if

      ! Write out source
      call read_source_bank(file_id)

    end if

    ! Close file
    call file_close(file_id)

  end subroutine load_state_point

!===============================================================================
! WRITE_SOURCE_BANK writes OpenMC source_bank data
!===============================================================================

  subroutine write_source_bank(group_id)
    use bank_header, only: Bank

    integer(HID_T), intent(in) :: group_id

    integer :: hdf5_err
    integer(HID_T) :: dset     ! data set handle
    integer(HID_T) :: dspace   ! data or file space handle
    integer(HID_T) :: memspace ! memory space handle
    integer(HSIZE_T) :: offset(1) ! source data offset
    integer(HSIZE_T) :: dims(1)
    type(c_ptr) :: f_ptr
#ifdef PHDF5
    integer :: data_xfer_mode
    integer(HID_T) :: plist    ! property list
#else
    integer :: i
#ifdef MPI
    type(Bank), allocatable, target :: temp_source(:)
#endif
#endif

#ifdef PHDF5
    ! Set size of total dataspace for all procs and rank
    dims(1) = n_particles
    call h5screate_simple_f(1, dims, dspace, hdf5_err)
    call h5dcreate_f(group_id, "source_bank", hdf5_bank_t, dspace, dset, hdf5_err)

    ! Create another data space but for each proc individually
    dims(1) = work
    call h5screate_simple_f(1, dims, memspace, hdf5_err)

    ! Select hyperslab for this dataspace
    offset(1) = work_index(rank)
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, offset, dims, hdf5_err)

    ! Set up the property list for parallel writing
    call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)
    call h5pset_dxpl_mpio_f(plist, H5FD_MPIO_COLLECTIVE_F, hdf5_err)

    ! Set up pointer to data
    f_ptr = c_loc(source_bank)

    ! Write data to file in parallel
    call h5dwrite_f(dset, hdf5_bank_t, f_ptr, hdf5_err, &
         file_space_id=dspace, mem_space_id=memspace, &
         xfer_prp=plist)

    ! Close all ids
    call h5sclose_f(dspace, hdf5_err)
    call h5sclose_f(memspace, hdf5_err)
    call h5dclose_f(dset, hdf5_err)
    call h5pclose_f(plist, hdf5_err)

#else

    if (master) then
      ! Create dataset big enough to hold all source sites
      dims(1) = n_particles
      call h5screate_simple_f(1, dims, dspace, hdf5_err)
      call h5dcreate_f(group_id, "source_bank", hdf5_bank_t, &
           dspace, dset, hdf5_err)

      ! Save source bank sites since the souce_bank array is overwritten below
#ifdef MPI
      allocate(temp_source(work))
      temp_source(:) = source_bank(:)
#endif

      do i = 0, n_procs - 1
        ! Create memory space
        dims(1) = work_index(i+1) - work_index(i)
        call h5screate_simple_f(1, dims, memspace, hdf5_err)

#ifdef MPI
        ! Receive source sites from other processes
        if (i > 0) then
          call MPI_RECV(source_bank, int(dims(1)), MPI_BANK, i, i, &
               MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_err)
        end if
#endif

        ! Select hyperslab for this dataspace
        call h5dget_space_f(dset, dspace, hdf5_err)
        offset(1) = work_index(i)
        call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, offset, dims, hdf5_err)

        ! Set up pointer to data and write data to hyperslab
        f_ptr = c_loc(source_bank)
        call h5dwrite_f(dset, hdf5_bank_t, f_ptr, hdf5_err, &
             file_space_id=dspace, mem_space_id=memspace)

        call h5sclose_f(memspace, hdf5_err)
        call h5sclose_f(dspace, hdf5_err)
      end do

      ! Close all ids
      call h5dclose_f(dset, hdf5_err)

      ! Restore state of source bank
#ifdef MPI
      source_bank(:) = temp_source(:)
      deallocate(temp_source)
#endif
    else
#ifdef MPI
      call MPI_SEND(source_bank, int(work), MPI_BANK, 0, rank, &
           MPI_COMM_WORLD, mpi_err)
#endif
    end if

#endif

  end subroutine write_source_bank

!===============================================================================
! READ_SOURCE_BANK reads OpenMC source_bank data
!===============================================================================

  subroutine read_source_bank(group_id)
    use bank_header, only: Bank

    integer(HID_T), intent(in) :: group_id

    integer :: hdf5_err
    integer(HID_T) :: dset     ! data set handle
    integer(HID_T) :: dspace   ! data space handle
    integer(HID_T) :: memspace ! memory space handle
    integer(HSIZE_T) :: dims(1)
    integer(HSIZE_T) :: offset(1) ! offset of data
    type(c_ptr) :: f_ptr
#ifdef PHDF5
    integer :: data_xfer_mode
    integer(HID_T) :: plist    ! property list
#endif

    ! Open the dataset
    call h5dopen_f(group_id, "source_bank", dset, hdf5_err)

    ! Create another data space but for each proc individually
    dims(1) = work
    call h5screate_simple_f(1, dims, memspace, hdf5_err)

    ! Select hyperslab for each process
    call h5dget_space_f(dset, dspace, hdf5_err)
    offset(1) = work_index(rank)
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, offset, dims, hdf5_err)

    ! Set up pointer to data
    f_ptr = c_loc(source_bank)

#ifdef PHDF5
    ! Read data in parallel
    call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)
    call h5pset_dxpl_mpio_f(plist, H5FD_MPIO_COLLECTIVE_F, hdf5_err)
    call h5dread_f(dset, hdf5_bank_t, f_ptr, hdf5_err, &
         file_space_id=dspace, mem_space_id=memspace, &
         xfer_prp=plist)
    call h5pclose_f(plist, hdf5_err)
#else
    call h5dread_f(dset, hdf5_bank_t, f_ptr, hdf5_err, &
         file_space_id=dspace, mem_space_id=memspace)
#endif

    ! Close all ids
    call h5sclose_f(dspace, hdf5_err)
    call h5sclose_f(memspace, hdf5_err)
    call h5dclose_f(dset, hdf5_err)

  end subroutine read_source_bank

end module state_point
