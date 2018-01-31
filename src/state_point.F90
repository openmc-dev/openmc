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

  use hdf5

  use cmfd_header
  use constants
  use eigenvalue,         only: openmc_get_keff
  use endf,               only: reaction_name
  use error,              only: fatal_error, warning, write_message
  use hdf5_interface
  use mesh_header,        only: RegularMesh, meshes, n_meshes
  use message_passing
  use mgxs_header,        only: nuclides_MG
  use nuclide_header,     only: nuclides
  use output,             only: time_stamp
  use random_lcg,         only: openmc_get_seed, openmc_set_seed
  use settings
  use simulation_header
  use string,             only: to_str, count_digits, zero_padded, to_f_string
  use tally_header
  use tally_filter_header
  use tally_derivative_header, only: tally_derivs
  use timer_header

  implicit none

contains

!===============================================================================
! OPENMC_STATEPOINT_WRITE writes an HDF5 statepoint file to disk
!===============================================================================

  subroutine openmc_statepoint_write(filename) bind(C)
    type(C_PTR), intent(in), optional :: filename

    integer :: i, j, k
    integer :: i_xs
    integer :: n_order      ! loop index for moment orders
    integer :: nm_order     ! loop index for Ynm moment orders
    integer, allocatable :: id_array(:)
    integer(HID_T) :: file_id
    integer(HID_T) :: cmfd_group, tallies_group, tally_group, meshes_group, &
                      filters_group, filter_group, derivs_group, &
                      deriv_group, runtime_group
    integer(C_INT) :: err
    real(C_DOUBLE) :: k_combined(2)
    character(MAX_WORD_LEN), allocatable :: str_array(:)
    character(C_CHAR), pointer :: string(:)
    character(len=:, kind=C_CHAR), allocatable :: filename_

    if (present(filename)) then
      call c_f_pointer(filename, string, [MAX_FILE_LEN])
      filename_ = to_f_string(string)
    else
      ! Set filename for state point
      filename_ = trim(path_output) // 'statepoint.' // &
           & zero_padded(current_batch, count_digits(n_max_batches))
      filename_ = trim(filename_) // '.h5'
    end if

    ! Write message
    call write_message("Creating state point " // trim(filename_) // "...", 5)

    if (master) then
      ! Create statepoint file
      file_id = file_create(filename_)

      ! Write file type
      call write_attribute(file_id, "filetype", "statepoint")

      ! Write revision number for state point file
      call write_attribute(file_id, "version", VERSION_STATEPOINT)

      ! Write OpenMC version
      call write_attribute(file_id, "openmc_version", VERSION)
#ifdef GIT_SHA1
      call write_attribute(file_id, "git_sha1", GIT_SHA1)
#endif

      ! Write current date and time
      call write_attribute(file_id, "date_and_time", time_stamp())

      ! Write path to input
      call write_attribute(file_id, "path", path_input)

      ! Write out random number seed
      call write_dataset(file_id, "seed", openmc_get_seed())

      ! Write run information
      if (run_CE) then
        call write_dataset(file_id, "energy_mode", "continuous-energy")
      else
        call write_dataset(file_id, "energy_mode", "multi-group")
      end if
      select case(run_mode)
      case (MODE_FIXEDSOURCE)
        call write_dataset(file_id, "run_mode", "fixed source")
      case (MODE_EIGENVALUE)
        call write_dataset(file_id, "run_mode", "eigenvalue")
      end select
      call write_dataset(file_id, "n_particles", n_particles)
      call write_dataset(file_id, "n_batches", n_batches)

      ! Write out current batch number
      call write_dataset(file_id, "current_batch", current_batch)

      ! Indicate whether source bank is stored in statepoint
      if (source_separate) then
        call write_attribute(file_id, "source_present", 0)
      else
        call write_attribute(file_id, "source_present", 1)
      end if

      ! Write out information for eigenvalue run
      if (run_mode == MODE_EIGENVALUE) then
        call write_dataset(file_id, "n_inactive", n_inactive)
        call write_dataset(file_id, "generations_per_batch", gen_per_batch)
        k = k_generation % size()
        call write_dataset(file_id, "k_generation", k_generation % data(1:k))
        if (entropy_on) then
          call write_dataset(file_id, "entropy", entropy % data(1:k))
        end if
        call write_dataset(file_id, "k_col_abs", k_col_abs)
        call write_dataset(file_id, "k_col_tra", k_col_tra)
        call write_dataset(file_id, "k_abs_tra", k_abs_tra)
        err = openmc_get_keff(k_combined)
        call write_dataset(file_id, "k_combined", k_combined)

        ! Write out CMFD info
        if (cmfd_on) then
          call write_attribute(file_id, "cmfd_on", 1)

          cmfd_group = create_group(file_id, "cmfd")
          call write_dataset(cmfd_group, "indices", cmfd % indices)
          call write_dataset(cmfd_group, "k_cmfd", cmfd % k_cmfd)
          call write_dataset(cmfd_group, "cmfd_src", cmfd % cmfd_src)
          call write_dataset(cmfd_group, "cmfd_entropy", cmfd % entropy)
          call write_dataset(cmfd_group, "cmfd_balance", cmfd % balance)
          call write_dataset(cmfd_group, "cmfd_dominance", cmfd % dom)
          call write_dataset(cmfd_group, "cmfd_srccmp", cmfd % src_cmp)
          call close_group(cmfd_group)
        else
          call write_attribute(file_id, "cmfd_on", 0)
        end if
      end if

      tallies_group = create_group(file_id, "tallies")

      ! Write number of meshes
      meshes_group = create_group(tallies_group, "meshes")
      call write_attribute(meshes_group, "n_meshes", n_meshes)

      if (n_meshes > 0) then
        ! Write IDs of meshes
        allocate(id_array(n_meshes))
        do i = 1, n_meshes
          id_array(i) = meshes(i) % id
        end do
        call write_attribute(meshes_group, "ids", id_array)
        deallocate(id_array)

        ! Write information for meshes
        MESH_LOOP: do i = 1, n_meshes
          call meshes(i) % to_hdf5(meshes_group)
        end do MESH_LOOP
      end if

      call close_group(meshes_group)

      ! Write information for derivatives.
      if (size(tally_derivs) > 0) then
        derivs_group = create_group(tallies_group, "derivatives")
        do i = 1, size(tally_derivs)
          associate(deriv => tally_derivs(i))
            deriv_group = create_group(derivs_group, "derivative " &
                 // trim(to_str(deriv % id)))
            select case (deriv % variable)
            case (DIFF_DENSITY)
              call write_dataset(deriv_group, "independent variable", "density")
              call write_dataset(deriv_group, "material", deriv % diff_material)
            case (DIFF_NUCLIDE_DENSITY)
              call write_dataset(deriv_group, "independent variable", &
                   "nuclide_density")
              call write_dataset(deriv_group, "material", deriv % diff_material)
              call write_dataset(deriv_group, "nuclide", &
                   nuclides(deriv % diff_nuclide) % name)
            case (DIFF_TEMPERATURE)
              call write_dataset(deriv_group, "independent variable", &
                   "temperature")
              call write_dataset(deriv_group, "material", deriv % diff_material)
            case default
              call fatal_error("Independent variable for derivative " &
                   // trim(to_str(deriv % id)) // " not defined in &
                   &state_point.F90.")
            end select
            call close_group(deriv_group)
          end associate
        end do
        call close_group(derivs_group)
      end if

      ! Write number of filters
      filters_group = create_group(tallies_group, "filters")
      call write_attribute(filters_group, "n_filters", n_filters)

      if (n_filters > 0) then
        ! Write IDs of filters
        allocate(id_array(n_filters))
        do i = 1, n_filters
          id_array(i) = filters(i) % obj % id
        end do
        call write_attribute(filters_group, "ids", id_array)
        deallocate(id_array)

        ! Write filter information
        FILTER_LOOP: do i = 1, n_filters
          filter_group = create_group(filters_group, "filter " // &
               trim(to_str(filters(i) % obj % id)))
          call filters(i) % obj % to_statepoint(filter_group)
          call close_group(filter_group)
        end do FILTER_LOOP
      end if

      call close_group(filters_group)

      ! Write number of tallies
      call write_attribute(tallies_group, "n_tallies", n_tallies)

      if (n_tallies > 0) then
        ! Write array of tally IDs
        allocate(id_array(n_tallies))
        do i = 1, n_tallies
          id_array(i) = tallies(i) % obj % id
        end do
        call write_attribute(tallies_group, "ids", id_array)
        deallocate(id_array)

        ! Write all tally information except results
        TALLY_METADATA: do i = 1, n_tallies

          ! Get pointer to tally
          associate (tally => tallies(i) % obj)
          tally_group = create_group(tallies_group, "tally " // &
               trim(to_str(tally % id)))

          ! Write the name for this tally
          call write_dataset(tally_group, "name", tally % name)

          select case(tally % estimator)
          case (ESTIMATOR_ANALOG)
            call write_dataset(tally_group, "estimator", "analog")
          case (ESTIMATOR_TRACKLENGTH)
            call write_dataset(tally_group, "estimator", "tracklength")
          case (ESTIMATOR_COLLISION)
            call write_dataset(tally_group, "estimator", "collision")
          end select
          call write_dataset(tally_group, "n_realizations", &
               tally % n_realizations)

          call write_dataset(tally_group, "n_filters", size(tally % filter))
          if (size(tally % filter) > 0) then
            ! Write IDs of filters
            allocate(id_array(size(tally % filter)))
            do j = 1, size(tally % filter)
              id_array(j) = filters(tally % filter(j)) % obj % id
            end do
            call write_dataset(tally_group, "filters", id_array)
            deallocate(id_array)
          end if

          ! Set up nuclide bin array and then write
          allocate(str_array(tally % n_nuclide_bins))
          NUCLIDE_LOOP: do j = 1, tally % n_nuclide_bins
            if (tally % nuclide_bins(j) > 0) then
              if (run_CE) then
                i_xs = index(nuclides(tally % nuclide_bins(j)) % name, '.')
                if (i_xs > 0) then
                  str_array(j) = nuclides(tally % nuclide_bins(j)) % name(1 : i_xs-1)
                else
                  str_array(j) = nuclides(tally % nuclide_bins(j)) % name
                end if
              else
                i_xs = index(nuclides_MG(tally % nuclide_bins(j)) % obj % name, '.')
                if (i_xs > 0) then
                  str_array(j) = nuclides_MG(tally % nuclide_bins(j)) % obj % name(1 : i_xs-1)
                else
                  str_array(j) = nuclides_MG(tally % nuclide_bins(j)) % obj % name
                end if
              end if
            else
              str_array(j) = 'total'
            end if
          end do NUCLIDE_LOOP
          call write_dataset(tally_group, "nuclides", str_array)
          deallocate(str_array)

          ! Write derivative information.
          if (tally % deriv /= NONE) then
            call write_dataset(tally_group, "derivative", &
                 tally_derivs(tally % deriv) % id)
          end if

          ! Write scores.
          call write_dataset(tally_group, "n_score_bins", tally % n_score_bins)
          allocate(str_array(size(tally % score_bins)))
          do j = 1, size(tally % score_bins)
            str_array(j) = reaction_name(tally % score_bins(j))
          end do
          call write_dataset(tally_group, "score_bins", str_array)
          call write_dataset(tally_group, "n_user_score_bins", &
               tally % n_user_score_bins)

          deallocate(str_array)

          ! Write explicit moment order strings for each score bin
          k = 1
          allocate(str_array(tally % n_score_bins))
          MOMENT_LOOP: do j = 1, tally % n_user_score_bins
            select case(tally % score_bins(k))
            case (SCORE_SCATTER_N, SCORE_NU_SCATTER_N)
              str_array(k) = trim(to_str(tally % moment_order(k)))
              k = k + 1
            case (SCORE_SCATTER_PN, SCORE_NU_SCATTER_PN)
              do n_order = 0, tally % moment_order(k)
                str_array(k) = trim(to_str(n_order))
                k = k + 1
              end do
            case (SCORE_SCATTER_YN, SCORE_NU_SCATTER_YN, SCORE_FLUX_YN, &
                 SCORE_TOTAL_YN)
              do n_order = 0, tally % moment_order(k)
                do nm_order = -n_order, n_order
                  str_array(k) = 'Y' // trim(to_str(n_order)) // ',' // &
                       trim(to_str(nm_order))
                  k = k + 1
                end do
              end do
            case default
              str_array(k) = ''
              k = k + 1
            end select
          end do MOMENT_LOOP

          call write_dataset(tally_group, "moment_orders", str_array)
          deallocate(str_array)

          call close_group(tally_group)
          end associate
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
      call write_dataset(file_id, "global_tallies", global_tallies)

      ! Write tallies
      if (active_tallies % size() > 0) then
        ! Indicate that tallies are on
        call write_attribute(file_id, "tallies_present", 1)

        tallies_group = open_group(file_id, "tallies")

        ! Write all tally results
        TALLY_RESULTS: do i = 1, n_tallies
          associate (tally => tallies(i) % obj)
            ! Write sum and sum_sq for each bin
            tally_group = open_group(tallies_group, "tally " &
                 // to_str(tally % id))
            call tally % write_results_hdf5(tally_group)
            call close_group(tally_group)
          end associate
        end do TALLY_RESULTS

        call close_group(tallies_group)
      else
        ! Indicate tallies are off
        call write_attribute(file_id, "tallies_present", 0)
      end if


      ! Write out the runtime metrics.
      runtime_group = create_group(file_id, "runtime")
      call write_dataset(runtime_group, "total initialization", &
           time_initialize % get_value())
      call write_dataset(runtime_group, "reading cross sections", &
           time_read_xs % get_value())
      call write_dataset(runtime_group, "simulation", &
           time_inactive % get_value() + time_active % get_value())
      call write_dataset(runtime_group, "transport", &
           time_transport % get_value())
      if (run_mode == MODE_EIGENVALUE) then
        call write_dataset(runtime_group, "inactive batches", &
             time_inactive % get_value())
      end if
      call write_dataset(runtime_group, "active batches", &
           time_active % get_value())
      if (run_mode == MODE_EIGENVALUE) then
        call write_dataset(runtime_group, "synchronizing fission bank", &
             time_bank % get_value())
        call write_dataset(runtime_group, "sampling source sites", &
             time_bank_sample % get_value())
        call write_dataset(runtime_group, "SEND-RECV source sites", &
             time_bank_sendrecv % get_value())
      end if
      call write_dataset(runtime_group, "accumulating tallies", &
           time_tallies % get_value())
      if (cmfd_run) then
        call write_dataset(runtime_group, "CMFD", time_cmfd % get_value())
        call write_dataset(runtime_group, "CMFD building matrices", &
             time_cmfdbuild % get_value())
        call write_dataset(runtime_group, "CMFD solving matrices", &
             time_cmfdsolve % get_value())
      end if
      call write_dataset(runtime_group, "total", time_total % get_value())
      call close_group(runtime_group)

      call file_close(file_id)
    end if
  end subroutine openmc_statepoint_write

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
             // "...", 5)

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
      call write_message("Creating source file " // trim(filename) // "...", 5)
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
    real(8), target :: global_temp(3,N_GLOBAL_TALLIES)
#ifdef OPENMC_MPI
    integer :: mpi_err ! MPI error code
    real(8) :: dummy   ! temporary receive buffer for non-root reduces
#endif
    type(TallyObject) :: dummy_tally

    ! ==========================================================================
    ! COLLECT AND WRITE GLOBAL TALLIES

    if (master) then
      ! Write number of realizations
      call write_dataset(file_id, "n_realizations", n_realizations)

      ! Write number of global tallies
      call write_dataset(file_id, "n_global_tallies", N_GLOBAL_TALLIES)

      tallies_group = open_group(file_id, "tallies")
    end if


#ifdef OPENMC_MPI
    ! Reduce global tallies
    n_bins = size(global_tallies)
    call MPI_REDUCE(global_tallies, global_temp, n_bins, MPI_REAL8, MPI_SUM, &
         0, mpi_intracomm, mpi_err)
#endif

    if (master) then
      ! Transfer values to value on master
      if (current_batch == n_max_batches .or. satisfy_triggers) then
        global_tallies(:,:) = global_temp(:,:)
      end if

      ! Write out global tallies sum and sum_sq
      call write_dataset(file_id, "global_tallies", global_temp)
    end if

    if (active_tallies % size() > 0) then
      ! Indicate that tallies are on
      if (master) then
        call write_attribute(file_id, "tallies_present", 1)
      end if

      ! Write all tally results
      TALLY_RESULTS: do i = 1, n_tallies
        associate (t => tallies(i) % obj)
          ! Determine size of tally results array
          m = size(t % results, 2)
          n = size(t % results, 3)
          n_bins = m*n*2

          ! Allocate array for storing sums and sums of squares, but
          ! contiguously in memory for each
          allocate(tally_temp(2,m,n))
          tally_temp(1,:,:) = t % results(RESULT_SUM,:,:)
          tally_temp(2,:,:) = t % results(RESULT_SUM_SQ,:,:)

          if (master) then
            tally_group = open_group(tallies_group, "tally " // &
                 trim(to_str(t % id)))

            ! The MPI_IN_PLACE specifier allows the master to copy values into
            ! a receive buffer without having a temporary variable
#ifdef OPENMC_MPI
            call MPI_REDUCE(MPI_IN_PLACE, tally_temp, n_bins, MPI_REAL8, &
                 MPI_SUM, 0, mpi_intracomm, mpi_err)
#endif

            ! At the end of the simulation, store the results back in the
            ! regular TallyResults array
            if (current_batch == n_max_batches .or. satisfy_triggers) then
              t % results(RESULT_SUM,:,:) = tally_temp(1,:,:)
              t % results(RESULT_SUM_SQ,:,:) = tally_temp(2,:,:)
            end if

            ! Put in temporary tally result
            allocate(dummy_tally % results(3,m,n))
            dummy_tally % results(RESULT_SUM,:,:) = tally_temp(1,:,:)
            dummy_tally % results(RESULT_SUM_SQ,:,:) = tally_temp(2,:,:)

            ! Write reduced tally results to file
            call dummy_tally % write_results_hdf5(tally_group)

            ! Deallocate temporary tally result
            deallocate(dummy_tally % results)
          else
            ! Receive buffer not significant at other processors
#ifdef OPENMC_MPI
            call MPI_REDUCE(tally_temp, dummy, n_bins, MPI_REAL8, MPI_SUM, &
                 0, mpi_intracomm, mpi_err)
#endif
          end if

          ! Deallocate temporary copy of tally results
          deallocate(tally_temp)

          if (master) call close_group(tally_group)
        end associate
      end do TALLY_RESULTS

      if (master) call close_group(tallies_group)
    else
      ! Indicate that tallies are off
      if (master) call write_dataset(file_id, "tallies_present", 0)
    end if


  end subroutine write_tally_results_nr

!===============================================================================
! LOAD_STATE_POINT
!===============================================================================

  subroutine load_state_point()

    integer :: i
    integer :: n
    integer :: int_array(3)
    integer, allocatable :: array(:)
    integer(C_INT64_T) :: seed
    integer(HID_T) :: file_id
    integer(HID_T) :: cmfd_group
    integer(HID_T) :: tallies_group
    integer(HID_T) :: tally_group
    logical :: source_present
    character(MAX_WORD_LEN) :: word

    ! Write message
    call write_message("Loading state point " // trim(path_state_point) &
         // "...", 5)

    ! Open file for reading
    file_id = file_open(path_state_point, 'r', parallel=.true.)

    ! Read filetype
    call read_attribute(word, file_id, "filetype")
    if (word /= 'statepoint') then
      call fatal_error("OpenMC tried to restart from a non-statepoint file.")
    end if

    ! Read revision number for state point file and make sure it matches with
    ! current version
    call read_attribute(array, file_id, "version")
    if (any(array /= VERSION_STATEPOINT)) then
      call fatal_error("State point version does not match current version &
           &in OpenMC.")
    end if

    ! Read and overwrite random number seed
    call read_dataset(seed, file_id, "seed")
    call openmc_set_seed(seed)

    ! It is not impossible for a state point to be generated from a CE run but
    ! to be loaded in to an MG run (or vice versa), check to prevent that.
    call read_dataset(word, file_id, "energy_mode")
    if (word == "multi-group" .and. run_CE) then
      call fatal_error("State point file is from multi-group run but &
                       & current run is continous-energy!")
    else if (word == "continuous-energy" .and. .not. run_CE) then
      call fatal_error("State point file is from continuous-energy run but &
                       & current run is multi-group!")
    end if

    ! Read and overwrite run information except number of batches
    call read_dataset(word, file_id, "run_mode")
    select case(word)
    case ('fixed source')
      run_mode = MODE_FIXEDSOURCE
    case ('eigenvalue')
      run_mode = MODE_EIGENVALUE
    end select
    call read_dataset(n_particles, file_id, "n_particles")
    call read_dataset(int_array(1), file_id, "n_batches")

    ! Take maximum of statepoint n_batches and input n_batches
    n_batches = max(n_batches, int_array(1))

    ! Read batch number to restart at
    call read_dataset(restart_batch, file_id, "current_batch")

    ! Check for source in statepoint if needed
    call read_attribute(int_array(1), file_id, "source_present")
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
      call read_dataset(int_array(1), file_id, "n_inactive")
      call read_dataset(gen_per_batch, file_id, "generations_per_batch")

      n = restart_batch*gen_per_batch
      call k_generation % resize(n)
      call read_dataset(k_generation % data(1:n), file_id, "k_generation")

      if (entropy_on) then
        call entropy % resize(n)
        call read_dataset(entropy % data(1:n), file_id, "entropy")
      end if
      call read_dataset(k_col_abs, file_id, "k_col_abs")
      call read_dataset(k_col_tra, file_id, "k_col_tra")
      call read_dataset(k_abs_tra, file_id, "k_abs_tra")

      ! Take maximum of statepoint n_inactive and input n_inactive
      n_inactive = max(n_inactive, int_array(1))

      ! Read in to see if CMFD was on
      call read_attribute(int_array(1), file_id, "cmfd_on")

      ! Read in CMFD info
      if (int_array(1) == 1) then
        cmfd_group = open_group(file_id, "cmfd")
        call read_dataset(cmfd % indices, cmfd_group, "indices")
        call read_dataset(cmfd % k_cmfd(1:restart_batch), cmfd_group, "k_cmfd")
        call read_dataset(cmfd % cmfd_src, cmfd_group, "cmfd_src")
        call read_dataset(cmfd % entropy(1:restart_batch), cmfd_group, &
             "cmfd_entropy")
        call read_dataset(cmfd % balance(1:restart_batch), cmfd_group, &
             "cmfd_balance")
        call read_dataset(cmfd % dom(1:restart_batch), cmfd_group, &
             "cmfd_dominance")
        call read_dataset(cmfd % src_cmp(1:restart_batch), cmfd_group, &
             "cmfd_srccmp")
        call close_group(cmfd_group)
      end if
    end if

    ! Check to make sure source bank is present
    if (path_source_point == path_state_point .and. .not. source_present) then
      call fatal_error("Source bank must be contained in statepoint restart &
           &file")
    end if

    ! Read tallies to master. If we are using Parallel HDF5, all processes
    ! need to be included in the HDF5 calls.
#ifdef PHDF5
    if (.true.) then
#else
    if (master) then
#endif

      ! Read number of realizations for global tallies
      call read_dataset(n_realizations, file_id, "n_realizations", indep=.true.)

      ! Read global tally data
      call read_dataset(global_tallies, file_id, "global_tallies")

      ! Check if tally results are present
      call read_attribute(int_array(1), file_id, "tallies_present")

      ! Read in sum and sum squared
      if (int_array(1) == 1) then
        tallies_group = open_group(file_id, "tallies")

        TALLY_RESULTS: do i = 1, n_tallies
          associate (t => tallies(i) % obj)
            ! Read sum, sum_sq, and N for each bin
            tally_group = open_group(tallies_group, "tally " // &
                 trim(to_str(t % id)))
            call t % read_results_hdf5(tally_group)
            call read_dataset(t % n_realizations, tally_group, &
                 "n_realizations")
            call close_group(tally_group)
          end associate
        end do TALLY_RESULTS

        call close_group(tallies_group)
      end if
    end if

    ! Read source if in eigenvalue mode
    if (run_mode == MODE_EIGENVALUE) then

      ! Check if source was written out separately
      if (.not. source_present) then

        ! Close statepoint file
        call file_close(file_id)

        ! Write message
        call write_message("Loading source file " // trim(path_source_point) &
             // "...", 5)

        ! Open source file
        file_id = file_open(path_source_point, 'r', parallel=.true.)

        ! Read file type
        call read_dataset(int_array(1), file_id, "filetype")

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
    integer(HID_T) :: plist    ! property list
#else
    integer :: i
#ifdef OPENMC_MPI
    integer :: mpi_err ! MPI error code
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
#ifdef OPENMC_MPI
      allocate(temp_source(work))
      temp_source(:) = source_bank(:)
#endif

      do i = 0, n_procs - 1
        ! Create memory space
        dims(1) = work_index(i+1) - work_index(i)
        call h5screate_simple_f(1, dims, memspace, hdf5_err)

#ifdef OPENMC_MPI
        ! Receive source sites from other processes
        if (i > 0) then
          call MPI_RECV(source_bank, int(dims(1)), MPI_BANK, i, i, &
               mpi_intracomm, MPI_STATUS_IGNORE, mpi_err)
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
#ifdef OPENMC_MPI
      source_bank(:) = temp_source(:)
      deallocate(temp_source)
#endif
    else
#ifdef OPENMC_MPI
      call MPI_SEND(source_bank, int(work), MPI_BANK, 0, rank, &
           mpi_intracomm, mpi_err)
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
    integer(HSIZE_T) :: dims(1)     ! dimensions on one processor
    integer(HSIZE_T) :: dims_all(1) ! overall dimensions
    integer(HSIZE_T) :: maxdims(1)  ! maximum dimensions
    integer(HSIZE_T) :: offset(1) ! offset of data
    type(c_ptr) :: f_ptr
#ifdef PHDF5
    integer(HID_T) :: plist    ! property list
#endif

    ! Open the dataset
    call h5dopen_f(group_id, "source_bank", dset, hdf5_err)

    ! Create another data space but for each proc individually
    dims(1) = work
    call h5screate_simple_f(1, dims, memspace, hdf5_err)

    ! Make sure source bank is big enough
    call h5dget_space_f(dset, dspace, hdf5_err)
    call h5sget_simple_extent_dims_f(dspace, dims_all, maxdims, hdf5_err)
    if (size(source_bank, KIND=HSIZE_T) > dims_all(1)) then
      call fatal_error("Number of source sites in source file is less than &
           &number of source particles per generation.")
    end if

    ! Select hyperslab for each process
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
