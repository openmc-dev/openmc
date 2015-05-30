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
  use dict_header,        only: ElemKeyValueII, ElemKeyValueCI
  use error,              only: fatal_error, warning
  use geometry_header,    only: Cell, Universe, Lattice, Surface
  use global
  use output,             only: write_message, time_stamp
  use source,             only: write_source_bank, read_source_bank
  use string,             only: to_str, zero_padded, count_digits
  use material_header,    only: Material
  use mesh_header,        only: StructuredMesh
  use output_interface
  use tally_header,       only: TallyObject, TallyResult
  use tally,              only: write_tally_result, read_tally_result

#ifdef HDF5
  use hdf5
  use hdf5_interface,   only: dims1, hdf5_rank, dset, dspace, hdf5_err, &
                              hdf5_open_group, hdf5_close_group, &
                              hdf5_tallyresult_t, memspace, plist, start1, &
                              count1, block1, f_ptr
  use, intrinsic :: ISO_C_BINDING
#endif

#ifdef MPI
  use mpi
#endif

  implicit none

  type(BinaryOutput)        :: sp      ! Statepoint/source output file
  type(BinaryOutput)        :: fh      ! general file handle

contains

!===============================================================================
! WRITE_STATE_POINT
!===============================================================================

  subroutine write_state_point()

    character(MAX_FILE_LEN)       :: filename
    integer                       :: i, j, k
    integer, allocatable          :: id_array(:)
    integer, allocatable          :: key_array(:)
    integer, allocatable          :: filter_map_array(:)
    integer                       :: otf_size_results_filters
    type(StructuredMesh), pointer :: mesh
    type(TallyObject), pointer    :: tally
    type(ElemKeyValueII), pointer :: current
    type(ElemKeyValueII), pointer :: next
    character(8)                  :: moment_name  ! name of moment (e.g, P3)
    integer                       :: n_order      ! loop index for moment orders
    integer                       :: nm_order     ! loop index for Ynm moment orders

    ! Start statepoint timer
    call time_statepoint % start()

    ! Set filename for state point
    filename = trim(path_output) // 'statepoint.' // &
        & zero_padded(current_batch, count_digits(n_max_batches))

    if (dd_run) then
      filename = trim(filename) // '.domain_' // &
          & zero_padded(domain_decomp % meshbin, &
                        count_digits(domain_decomp % n_domains))
    end if

    ! Append appropriate extension
#ifdef HDF5
    filename = trim(filename) // '.h5'
#else
    filename = trim(filename) // '.binary'
#endif

    ! Write message
    call write_message("Creating state point " // trim(filename) // "...", 1)

    if (master .or. (dd_run .and. domain_decomp % local_master)) then

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

      ! Write domain decomposition information
      if (dd_run) then
        call sp % write_data(1, "domain_decomp")
        call sp % write_data(domain_decomp % n_domains, "n_domains")
        call sp % write_data(domain_decomp % meshbin, "domain_id")
      else
        call sp % write_data(0, "domain_decomp")
        call sp % write_data(NONE, "n_domains")
        call sp % write_data(NONE, "domain_id")
      end if

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

          ! Write on-the-fly allocation tally info
          if (tally % on_the_fly_allocation) then
            otf_size_results_filters = tally % next_filter_idx - 1
            call sp % write_data(otf_size_results_filters, &
                 "otf_size_results_filters", group="tallies/tally" // &
                  trim(to_str(tally % id)))
            ! Write otf filter bin mapping
            allocate(filter_map_array(otf_size_results_filters))
            do j = 1, otf_size_results_filters
              filter_map_array(j) = tally % reverse_filter_index_map % get_key(j)
            end do
            call sp % write_data(filter_map_array, "otf_filter_bin_map", &
                 group="tallies/tally" // trim(to_str(tally % id)), &
                 length=otf_size_results_filters)
            deallocate(filter_map_array)
          else
            call sp % write_data(NONE, "otf_size_results_filters", &
                 group="tallies/tally" // trim(to_str(tally % id)))
          end if

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

      if (dd_run) then
        call fatal_error('no_reduce not implemented with domain decomposition')
      end if

      call write_tally_results_nr()

    elseif (master .or. (dd_run .and. domain_decomp % local_master)) then

      ! Write number of global realizations
      call sp % write_data(n_realizations, "n_realizations")

      ! Write global tallies
      call sp % write_data(N_GLOBAL_TALLIES, "n_global_tallies")
      call write_tally_result(sp, global_tallies, "global_tallies", &
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
          if (tally % on_the_fly_allocation) then

            otf_size_results_filters = tally % next_filter_idx - 1
            call write_tally_result(sp, &
                tally % results(:,1:otf_size_results_filters), "results", &
                group="tallies/tally" // trim(to_str(tally % id)), &
                n1=size(tally % results, 1), n2=otf_size_results_filters)

          else

            call write_tally_result(sp, tally % results, "results", &
                group="tallies/tally" // trim(to_str(tally % id)), &
                n1=size(tally % results, 1), n2=size(tally % results, 2))

          endif

        end do TALLY_RESULTS

      else

        ! Indicate tallies are off
        call sp % write_data(0, "tallies_present", group="tallies")

      end if

      ! Close the file for serial writing
      call sp % file_close()

    end if

#ifdef HDF5
    ! Write OTF tallies from DD runs, which were deferred for HDF5
    ! TODO: this function unscrambles the results arrays and writes to one
    ! aggregated statepoint for all domains.  It's terribly inefficient as
    ! written, so for now we leave it commented out and just write separate
    ! statepoints for each domain, which will need to be unscrambled in
    ! post-processing
!    call write_state_point_otf_tally_data(filename)
#endif

    if (master .and. n_tallies > 0) then
      deallocate(id_array)
    end if

    ! Stop statepoint timer
    call time_statepoint % stop()

  end subroutine write_state_point

#ifdef HDF5
!===============================================================================
! WRITE_STATE_POINT_OTF_TALLY_DATA
!===============================================================================

  subroutine write_state_point_otf_tally_data(filename)

    character(MAX_FILE_LEN), intent(in) :: filename

    character(MAX_FILE_LEN) :: groupname
    integer :: i, j
    integer :: n, m
    integer :: idx
    type(TallyObject), pointer :: t => null()

    integer(HID_T) :: file_id
    integer(HID_T) :: group_id

    ! Set up OTF tally datasets
    if (master) then
      do i = 1, n_tallies

        ! Set point to current tally
        t => tallies(i)

        ! Write sum and sum_sq for each bin
        if (t % on_the_fly_allocation) then

          n = t % total_score_bins
          m = t % total_filter_bins

          hdf5_rank = 1
          dims1(1) = n*m
          groupname = "tallies/tally" // to_str(i)
          call h5fopen_f(trim(filename), H5F_ACC_RDWR_F, file_id, hdf5_err)
          call h5gopen_f(file_id, trim(groupname), group_id, hdf5_err)
          call h5screate_simple_f(hdf5_rank, dims1, dspace, hdf5_err)
          call h5dcreate_f(group_id, 'results', hdf5_tallyresult_t, dspace, &
              dset, hdf5_err)
          call h5dclose_f(dset, hdf5_err)
          call h5sclose_f(dspace, hdf5_err)
          call h5gclose_f(group_id, hdf5_err)
          call h5fclose_f(file_id, hdf5_err)

        end if

      end do
    end if

# ifdef MPI
    ! All other domains need to wait for the datasets to be created before they
    ! can write to it
    call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)
# endif

    ! Write tally data to the file
    if (master .or. (dd_run .and. domain_decomp % local_master)) then

      ! Setup file access property list with parallel I/O access
      call h5pcreate_f(H5P_FILE_ACCESS_F, plist, hdf5_err)
      call h5pset_fapl_mpio_f(plist, domain_decomp % comm_domain_masters, &
          MPI_INFO_NULL, hdf5_err)

      ! Open the file
      call h5fopen_f(trim(filename), H5F_ACC_RDWR_F, file_id, hdf5_err, &
          access_prp = plist)

      ! Close property list
      call h5pclose_f(plist, hdf5_err)

      ! Create the property list to describe independent parallel I/O
      call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)
      call h5pset_dxpl_mpio_f(plist, H5FD_MPIO_INDEPENDENT_F, hdf5_err)

      count1 = 1

      do i = 1, n_tallies

        ! Set point to current tally
        t => tallies(i)

        if (t % on_the_fly_allocation) then

          ! Skip this tally if no bins were scored to
          if (t % next_filter_idx == 1) cycle

          if (master) then
            call write_message("Writing distributed OTF tally " // &
                               trim(to_str(t % id)) // "...", 6)
          end if

          block1 = size(t % results(:,1))

          ! Open the group
          groupname = 'tallies/tally' // trim(to_str(i))
          call h5gopen_f(file_id, trim(groupname), group_id, hdf5_err)

          ! Open the dataset
          call h5dopen_f(group_id, 'results', dset, hdf5_err)

          ! Open the dataspace and memory space
          call h5dget_space_f(dset, dspace, hdf5_err)
          call h5screate_simple_f(1, block1 * (t % next_filter_idx - 1), &
              memspace, hdf5_err)

          ! For on-the-fly distributed tallies, we need to unscramble. We
          ! do this by selecting an irregular hyperslab in the dataset, and
          ! ensure that the filter bins are in the same order. This means we
          ! need to sort the results bins.
          call heapsort_results(t)

          ! Loop through OTF filter bins and write
          call h5sselect_none_f(dspace, hdf5_err)
          do j = 1, t % next_filter_idx - 1

            idx = t % reverse_filter_index_map % get_key(j)
            start1 = (idx - 1) * block1

            ! Select the hyperslab
            call h5sselect_hyperslab_f(dspace, H5S_SELECT_OR_F, start1, &
                count1, hdf5_err, block = block1)

          end do

          ! Write the data
          f_ptr = c_loc(t % results(1,1))
          call h5dwrite_f(dset, hdf5_tallyresult_t, f_ptr, hdf5_err, &
              file_space_id = dspace, mem_space_id = memspace, &
              xfer_prp = plist)

          ! Close the dataspace and memory space
          call h5sclose_f(dspace, hdf5_err)
          call h5sclose_f(memspace, hdf5_err)

          ! Close dataset and group
          call h5dclose_f(dset, hdf5_err)
          call h5gclose_f(group_id, hdf5_err)

        end if

      end do

      ! Close property list
      call h5pclose_f(plist, hdf5_err)

      ! Close the file
      call h5fclose_f(file_id, hdf5_err)

    end if

  end subroutine write_state_point_otf_tally_data
#endif

!===============================================================================
! WRITE_SOURCE_POINT
!===============================================================================

  subroutine write_source_point()

    type(BinaryOutput) :: sp
    character(MAX_FILE_LEN) :: filename

    if (dd_run) then
      if (master) call warning('Source bank writing not implemented for DD runs.')
      return
    end if

    ! Check to write out source for a specified batch
    if (sourcepoint_batch % contains(current_batch)) then

      ! Create or open up file
      if (source_separate) then

        ! Set filename
        filename = trim(path_output) // 'source.' // &
            & zero_padded(current_batch, count_digits(n_max_batches))

        if (dd_run) then
          filename = trim(filename) // '.domain_' // &
              & zero_padded(domain_decomp % meshbin, &
                            count_digits(domain_decomp % n_domains))
        end if

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
            & zero_padded(current_batch, count_digits(n_batches))

        if (dd_run) then
          filename = trim(filename) // '.domain_' // &
              & zero_padded(domain_decomp % meshbin, &
                            count_digits(domain_decomp % n_domains))
        end if

#ifdef HDF5
        filename = trim(filename) // '.h5'
#else
        filename = trim(filename) // '.binary'
#endif

        ! Reopen statepoint file in parallel
        call sp % file_open(filename, 'w', serial = .false.)

      end if

      ! Write out source
      call write_source_bank(sp)

      ! Close file
      call sp % file_close()

    end if

    ! Also check to write source separately in overwritten file
    if (source_latest) then

      ! Set filename
      filename = trim(path_output) // 'source'

      if (dd_run) then
        filename = trim(filename) // '.domain_' // &
            & zero_padded(domain_decomp % meshbin, &
                          count_digits(domain_decomp % n_domains))
      end if

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
      call write_source_bank(sp)

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
      call write_tally_result(sp, tallyresult_temp, "global_tallies", &
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
          call write_tally_result(sp, tally % results, "results", &
               group="tallies/tally" // trim(to_str(tally % id)), n1=m, n2=n)

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
    integer, allocatable       :: filter_map_array(:)
    integer                    :: otf_size_results_filters
    integer                    :: dummy_filter_index
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

      ! Read on-the-fly allocation tally info
      if (tally % on_the_fly_allocation) then
        call sp % read_data(otf_size_results_filters, "otf_size_results_filters", &
             group="tallies/tally" // trim(to_str(tally % id)))
        
        ! Read otf filter bin mapping
        allocate(filter_map_array(otf_size_results_filters))
        call sp % read_data(filter_map_array, "otf_filter_bin_map", &
             group="tallies/tally" // trim(to_str(tally % id)), &
             length=otf_size_results_filters)

        ! Reset the filter map on the tally object
        do j = 1, otf_size_results_filters
          dummy_filter_index = tally % otf_filter_index(filter_map_array(j))
        end do

        deallocate(filter_map_array)

      else
        call sp % read_data(otf_size_results_filters, "otf_size_results_filters", &
             group="tallies/tally" // trim(to_str(tally % id)))
      end if

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
      call read_tally_result(sp, global_tallies, "global_tallies", &
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
          call read_tally_result(sp, tally % results, "results", &
               group="tallies/tally" // trim(to_str(curr_key)), &
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
      call read_source_bank(sp)

    end if

    ! Close file
    call sp % file_close()

  end subroutine load_state_point

  subroutine read_source
! TODO write this routine
! TODO what if n_particles does not match source bank
  end subroutine read_source

!===============================================================================
! SYNCHRONIZE_OTF_MATERIALS ensures that all processes within each domain have
! loaded the same otf materials.  These will be sorted to the same order in
! write_distribmat_comps
!===============================================================================

  subroutine synchronize_otf_materials(mat, reduce_comm)

    type(Material), pointer, intent(inout) :: mat
    integer,                    intent(in) :: reduce_comm

    integer :: j
    integer :: p
    integer :: idx
    integer :: real_bin
    integer :: reduce_n_procs
    integer :: reduce_rank
    logical :: reduce_master
    integer :: n_comps
    integer :: max_comps ! max number of compositions on any proc
    integer, allocatable :: proc_n_comps(:) ! number of compositions on proc
    integer, allocatable :: proc_comp_map(:,:) ! index -> real bin on proc
    integer              :: n_request
    integer              :: request(10000) ! TODO: make this something real!

    !=======================================================================
    ! SEND MAPS TO LOCAL MASTER

    reduce_master = .false.
    call MPI_COMM_SIZE(reduce_comm, reduce_n_procs, mpi_err)
    call MPI_COMM_RANK(reduce_comm, reduce_rank, mpi_err)
    if (reduce_rank == 0) reduce_master = .true.
    
    ! Master needs to know the max number of compositions
    call MPI_REDUCE(mat % next_comp_idx - 1, max_comps, 1, MPI_INTEGER, &
         MPI_MAX, 0, reduce_comm, mpi_err) 

    n_request = 0
    if (reduce_master) then

      allocate(proc_n_comps(reduce_n_procs - 1))
      allocate(proc_comp_map(max_comps, reduce_n_procs - 1))

      ! Receive composition information from all other procs in domain
      do p = 1, reduce_n_procs - 1

        ! Receive number of compositions
        n_request = n_request + 1
        call MPI_IRECV(proc_n_comps(p), 1, MPI_INTEGER, p, p, &
            reduce_comm, request(n_request), mpi_err)

        ! Receive composition maps
        n_request = n_request + 1
        call MPI_IRECV(proc_comp_map(:,p), max_comps, MPI_INTEGER, p, &
            p, reduce_comm, request(n_request), mpi_err)

      end do

    else

      ! Build composition  map
      n_comps = mat % next_comp_idx - 1
      allocate(proc_comp_map(n_comps, 1))
      do j = 1, n_comps
        real_bin = mat % reverse_comp_index_map % get_key(j)
        proc_comp_map(j, 1) = real_bin
      end do

      ! Send number of compositions
      n_request = n_request + 1
      call MPI_ISEND(n_comps, 1, MPI_INTEGER, &
          0, reduce_rank, &
          reduce_comm, request(n_request), mpi_err)

      ! Send composition maps
      n_request = n_request + 1
      call MPI_ISEND(proc_comp_map(:,1), n_comps, MPI_INTEGER, &
          0, reduce_rank, &
          reduce_comm, request(n_request), mpi_err)

    end if

    ! Wait for composition information to synchronize
    call MPI_WAITALL(n_request, request, MPI_STATUSES_IGNORE, mpi_err)

    !=======================================================================
    ! SYNCHRONIZE MAPS

    if (reduce_master) then

      ! go through all received maps and load materials
      do p = 1, reduce_n_procs - 1
        do j = 1, proc_n_comps(p)
          real_bin = proc_comp_map(j, p)
          idx = mat % otf_comp_index(real_bin)
        end do
      end do

      ! Build composition map
      deallocate(proc_comp_map)
      n_comps = mat % next_comp_idx - 1
      allocate(proc_comp_map(n_comps, 1))
      do j = 1, n_comps
        real_bin = mat % reverse_comp_index_map % get_key(j)
        proc_comp_map(j, 1) = real_bin
      end do

      ! broadcast the updated master n_comp
      call MPI_BCAST(n_comps, 1, MPI_INTEGER, 0, reduce_comm, mpi_err)

      ! broadcast the updated master composition map
      call MPI_BCAST(proc_comp_map(:, 1), n_comps, MPI_INTEGER, 0, &
          reduce_comm, mpi_err)

    else

      ! receive master n_comps
      call MPI_BCAST(n_comps, 1, MPI_INTEGER, 0, reduce_comm, mpi_err)

      deallocate(proc_comp_map)
      allocate(proc_comp_map(n_comps, 1))

      ! recieve master composition map
      call MPI_BCAST(proc_comp_map(:, 1), n_comps, MPI_INTEGER, 0, &
          reduce_comm, mpi_err)

      ! load any materials present on master
      do j = 1, n_comps
        real_bin = proc_comp_map(j, 1)
        idx = mat % otf_comp_index(real_bin)
      end do

    end if

  end subroutine synchronize_otf_materials

!===============================================================================
! GET_COMP_ARRAY_OWNERSHIP_SLICE
!===============================================================================

  subroutine get_comp_array_ownership_boundaries(mat, slice_start, slice_end)

    type(Material), pointer, intent(inout) :: mat
    integer,                 intent(out)   :: slice_start, slice_end

    integer :: arr_size
    integer :: minwork, remainder
    integer :: n_slice, slice_idx, slice_size

    ! Determine slicing params based on OTF or DD
    arr_size = mat % n_comp
    slice_idx = rank
    n_slice = n_procs
    if (mat % otf_compositions) arr_size = mat % next_comp_idx - 1
    if (dd_run) then
      slice_idx = domain_decomp % rank
      n_slice = domain_decomp % n_domain_procs
    end if

    ! Determine equal slicing - if there are more processes than array positions,
    ! slice_start should be higher than slice_end
    minwork = (arr_size)/n_slice
    remainder = mod(arr_size, n_slice)
    if (slice_idx < remainder) then
      slice_size = minwork + 1
      slice_start = slice_idx * (minwork + 1) + 1
    else
      slice_size = minwork
      slice_start = slice_idx * (minwork + 1) - slice_idx + remainder + 1
    endif
    slice_end = slice_start + slice_size - 1

  end subroutine get_comp_array_ownership_boundaries

!===============================================================================
! WRITE_DISTRIBMAT_COMPS
!===============================================================================

  subroutine write_distribmat_comps(filename)

    character(MAX_FILE_LEN), intent(in)    :: filename

    character(MAX_FILE_LEN) :: groupname
    integer :: i, j
    integer :: idx
    integer :: slice_start, slice_end
    type(BinaryOutput) :: fh
    type(Material), pointer :: mat => null()

#ifdef HDF5
    integer(HID_T) :: file_id
    integer(HID_T) :: group_id
    integer(HSIZE_T) :: chunk(1)
    integer(HID_T) :: chunk_plist
#endif

#ifndef HDF5
    call warning('Distribmat material writing only implemented for HDF5.')
    return
#else

    ! Start matdump timer
    call time_matdump % start()

    ! Create files and write headers (master only)
    if (master) then

      call fh % file_create(filename)

      do i = 1, n_materials
        mat => materials(i)
        if (mat % n_comp > 1) then

          groupname = 'mat-' // trim(to_str(mat % id))

          ! Create file and write header
          call fh % file_open(filename, 'w')
          call fh % write_data(mat % n_nuclides, 'n_nuclides', &
              group=trim(groupname), record=1)
          call fh % write_data(mat % n_comp, 'n_instances', &
              group=trim(groupname), record=2)
          call fh % file_close()

          ! Create the full dataset initially so all other procs can write to it
          hdf5_rank = 1
          dims1(1) = mat % n_comp * mat % n_nuclides
          call h5fopen_f(trim(filename), H5F_ACC_RDWR_F, file_id, hdf5_err)
          call h5gopen_f(file_id, trim(groupname), group_id, hdf5_err)
          call h5screate_simple_f(hdf5_rank, dims1, dspace, hdf5_err)
          call h5pcreate_f(H5P_DATASET_CREATE_F, chunk_plist, hdf5_err)
          ! Tune chunking and chunk caching to the filesystem if performance is needed
          chunk(1) = mat % n_nuclides *  800
!          call h5pset_chunk_f(chunk_plist, 1, chunk, hdf5_err)
!          call h5pset_chunk_cache_f(chunk_plist, 0_8, 0_8, 1.0_4, hdf5_err) ! Turn chunk caching off
!          call h5pset_chunk_cache_f(chunk_plist, 211_8, 16777216_8, 1.0_4, hdf5_err) ! OR: tune chunk cache to filesystem
          ! Set the fill value if needed for debugging (slow for large datasets)
!          call h5pset_fill_value_f(chunk_plist, H5T_NATIVE_DOUBLE, -1.0_8, hdf5_err)
          call h5dcreate_f(group_id, 'comps', H5T_NATIVE_DOUBLE, dspace, dset, hdf5_err, &
              dcpl_id = chunk_plist)
          call h5pclose_f(chunk_plist, hdf5_err)
          call h5dclose_f(dset, hdf5_err)
          call h5sclose_f(dspace, hdf5_err)
          call h5gclose_f(group_id, hdf5_err)
          call h5fclose_f(file_id, hdf5_err)

        end if
      end do
    end if

# ifdef MPI
    ! For parallel IO we need to wait for master to create the files
    call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)
# endif

    do i = 1, n_materials
      mat => materials(i)
      if (mat % n_comp > 1) then
        if (mat % otf_compositions) then

          call write_message("Preparing distributed material " // &
                             trim(to_str(mat % id)) // " for writing...", 6)

          ! All processes need to participate in IO for efficiency on HPC
          ! machines. We want this anyways for large arrays, where each
          ! processor dumps a slice. As a result, we need to make sure all
          ! processes have the same array before setting slice boundaries. If
          ! enough particles were run the OTF arrays would have the same
          ! elements anyways, but for large problems this isn't always
          ! guaranteed.
          if (dd_run) then
            call synchronize_otf_materials(mat, domain_decomp % comm)
          else
            call synchronize_otf_materials(mat, MPI_COMM_WORLD)
          end if

          ! For on-the-fly distributes materials, we need to unscramble. We
          ! do this by selecting an irregular hyperslab in the dataset, and
          ! ensure that the materials are in the same order.  This means we
          ! need to sort the compositions.
          call heapsort_matcomps(mat)

        end if
      end if
    end do

    ! Setup file access property list with parallel I/O access
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist, hdf5_err)
    call h5pset_fapl_mpio_f(plist, MPI_COMM_WORLD, MPI_INFO_NULL, hdf5_err)

    ! Open the file
    call h5fopen_f(trim(filename), H5F_ACC_RDWR_F, file_id, hdf5_err, &
        access_prp = plist)

    ! Close property list
    call h5pclose_f(plist, hdf5_err)

    ! Create the property list to describe collective parallel I/O
    call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)
    call h5pset_dxpl_mpio_f(plist, H5FD_MPIO_COLLECTIVE_F, hdf5_err)

    count1 = 1

    do i = 1, n_materials
      mat => materials(i)
      if (mat % n_comp > 1) then
        
        groupname = 'mat-' // trim(to_str(mat % id))

        if (master) then
          call write_message("Writing distributed material " // &
                             trim(to_str(mat % id)) // "...", 6)
        end if

        block1 = mat % n_nuclides

        ! Calculate boundaries of compositions array slice for this process
        call get_comp_array_ownership_boundaries(mat, slice_start, slice_end)

        ! Open the group
        call h5gopen_f(file_id, trim(groupname), group_id, hdf5_err)

        ! Open the dataset
        call h5dopen_f(group_id, 'comps', dset, hdf5_err)

        if (mat % otf_compositions) then
          ! OTF mats are spread across all ranks

          ! Open the dataspace and memory space
          call h5dget_space_f(dset, dspace, hdf5_err)
          call h5screate_simple_f(1, block1 * (slice_end - slice_start + 1), &
              memspace, hdf5_err)

          ! Select the irregular hyperslab
          call h5sselect_none_f(dspace, hdf5_err)
          do j = slice_start, slice_end

            idx = mat % reverse_comp_index_map % get_key(j)
            start1 = (idx - 1) * block1

            ! Select the hyperslab
            call h5sselect_hyperslab_f(dspace, H5S_SELECT_OR_F, start1, &
                count1, hdf5_err, block = block1)

          end do

          ! Write the data
          f_ptr = c_loc(mat % otf_comp(1, slice_start))
          call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, f_ptr, hdf5_err, &
              file_space_id = dspace, mem_space_id = memspace, &
              xfer_prp = plist)

          ! Close the dataspace and memory space
          call h5sclose_f(dspace, hdf5_err)
          call h5sclose_f(memspace, hdf5_err)

        else

          ! TODO: change non-otf distrib comps to use 2d array instead of
          ! comp datastructure - would be more efficient, and allow this
          ! section of code to be combined with the previous

          ! For normal distribmats, the compositions are in order
          do j = 1, mat % n_comp

            ! Open the dataspace and memory space
            call h5dget_space_f(dset, dspace, hdf5_err)
            call h5screate_simple_f(1, block1, memspace, hdf5_err)

            start1 = (j - 1) * block1

            ! Select the hyperslab
            call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, start1, &
                count1, hdf5_err, block = block1)

            ! Write the data
            f_ptr = c_loc(mat % comp(j) % atom_density)
            call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, f_ptr, hdf5_err, &
                file_space_id = dspace, mem_space_id = memspace, &
                xfer_prp = plist)

            ! Close the dataspace and memory space
            call h5sclose_f(dspace, hdf5_err)
            call h5sclose_f(memspace, hdf5_err)

          end do

        end if

        ! Close dataset and group
        call h5dclose_f(dset, hdf5_err)
        call h5gclose_f(group_id, hdf5_err)

      end if
    end do

    ! Close property list
    call h5pclose_f(plist, hdf5_err)

    ! Close the file
    call h5fclose_f(file_id, hdf5_err)

#endif

#ifdef MPI
    ! Everyone should wait here - master needs a good matdump time
    call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)
#endif

    ! Start matdump timer
    call time_matdump % stop()

  end subroutine write_distribmat_comps

!===============================================================================
! HEAPSORT_MATCOMPS performs a heapsort on OTF compositions within a material
! to order the array by real distribcell index
!===============================================================================

  subroutine heapsort_matcomps(mat)
   
    type(Material), pointer, intent(inout) :: mat

    integer :: start, n, bottom

    ! Build the heap - O(log(n))
    n = mat % next_comp_idx - 1
    do start = (n - 2) / 2, 0, -1
      call siftdown_matcomps(mat, start, n);
    end do
   
    ! Do the sort - O(n)
    do bottom = n - 1, 1, -1
      call swap_matcomps(mat, 1, bottom + 1)
      call siftdown_matcomps(mat, 0, bottom)
    end do
   
  end subroutine heapsort_matcomps

!===============================================================================
! SWAP_MATCOMPS swaps two sections of the otf_comp array in a material, and
! updates the otf mapping dictionaries
!===============================================================================
  
  subroutine swap_matcomps(mat, a, b)
    type(Material), pointer, intent(inout) :: mat
    integer, intent(in) :: a, b ! actual composition indices in otf_comp

    real(8), allocatable :: tmp(:)
    integer :: real_inst_a, real_inst_b

    allocate(tmp(mat % n_nuclides))

    ! Swap the maps
    real_inst_a = mat % reverse_comp_index_map % get_key(a)
    real_inst_b = mat % reverse_comp_index_map % get_key(b)
    call mat % reverse_comp_index_map % add_key(a, real_inst_b)
    call mat % reverse_comp_index_map % add_key(b, real_inst_a)
    call mat % comp_index_map % add_key(real_inst_a, b)
    call mat % comp_index_map % add_key(real_inst_b, a)

    ! Swap the compositions
    tmp = mat % otf_comp(:, b)
    mat % otf_comp(:, b) = mat % otf_comp(:, a)
    mat % otf_comp(:, a) = tmp

    deallocate(tmp)

  end subroutine swap_matcomps

!===============================================================================
! SIFTDOWN_MATCOMPS
!===============================================================================

  subroutine siftdown_matcomps(mat, start, bottom)
   
    type(Material), pointer, intent(inout) :: mat
    integer, intent(in) :: start, bottom

    integer :: child, root
    integer :: real_inst_a, real_inst_b

    root = start
    do while(root*2 + 1 < bottom)
      child = root * 2 + 1
   
      if (child + 1 < bottom) then
        real_inst_a = mat % reverse_comp_index_map % get_key(child + 1)
        real_inst_b = mat % reverse_comp_index_map % get_key(child + 2)
        if (real_inst_a < real_inst_b) child = child + 1
      end if
   
      real_inst_a = mat % reverse_comp_index_map % get_key(root + 1)
      real_inst_b = mat % reverse_comp_index_map % get_key(child + 1)
      if (real_inst_a < real_inst_b) then
        call swap_matcomps(mat, root + 1, child + 1)
        root = child
      else
        return
      end if

    end do      
   
  end subroutine siftdown_matcomps

!===============================================================================
! HEAPSORT_RESULTS performs a heapsort on OTF results arrays within a tally
! to order the array by real filter indices
!===============================================================================

  subroutine heapsort_results(t)
   
    type(TallyObject), pointer, intent(inout) :: t

    integer :: start, n, bottom

    ! Build the heap - O(log(n))
    n = t % next_filter_idx - 1
    do start = (n - 2) / 2, 0, -1
      call siftdown_results(t, start, n);
    end do
   
    ! Do the sort - O(n)
    do bottom = n - 1, 1, -1
      call swap_results(t, 1, bottom + 1)
      call siftdown_results(t, 0, bottom)
    end do
   
  end subroutine heapsort_results

!===============================================================================
! SWAP_RESULTS swaps two sections of the results array in a tally, and
! updates the otf mapping dictionaries
!===============================================================================
  
  subroutine swap_results(t, a, b)
    type(TallyObject), pointer, intent(inout) :: t
    integer, intent(in) :: a, b ! actual composition indices in otf_comp

    type(TallyResult), allocatable :: tmp(:)
    integer :: real_inst_a, real_inst_b

    allocate(tmp(t % total_score_bins))

    ! Swap the maps
    real_inst_a = t % reverse_filter_index_map % get_key(a)
    real_inst_b = t % reverse_filter_index_map % get_key(b)
    call t % reverse_filter_index_map % add_key(a, real_inst_b)
    call t % reverse_filter_index_map % add_key(b, real_inst_a)
    call t % filter_index_map % add_key(real_inst_a, b)
    call t % filter_index_map % add_key(real_inst_b, a)

    ! Swap the results
    tmp = t % results(:, b)
    t % results(:, b) = t % results(:, a)
    t % results(:, a) = tmp

    deallocate(tmp)

  end subroutine swap_results

!===============================================================================
! SIFTDOWN_RESULTS
!===============================================================================

  subroutine siftdown_results(t, start, bottom)
   
    type(TallyObject), pointer, intent(inout) :: t
    integer, intent(in) :: start, bottom

    integer :: child, root
    integer :: real_inst_a, real_inst_b

    root = start
    do while(root*2 + 1 < bottom)
      child = root * 2 + 1
   
      if (child + 1 < bottom) then
        real_inst_a = t % reverse_filter_index_map % get_key(child + 1)
        real_inst_b = t % reverse_filter_index_map % get_key(child + 2)
        if (real_inst_a < real_inst_b) child = child + 1
      end if
   
      real_inst_a = t % reverse_filter_index_map % get_key(root + 1)
      real_inst_b = t % reverse_filter_index_map % get_key(child + 1)
      if (real_inst_a < real_inst_b) then
        call swap_results(t, root + 1, child + 1)
        root = child
      else
        return
      end if

    end do      
   
  end subroutine siftdown_results

end module state_point
