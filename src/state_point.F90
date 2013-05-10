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

  use error,         only: warning, fatal_error
  use global
  use hdf5_interface
  use math,          only: t_percentile
  use output,        only: write_message, print_batch_keff, time_stamp
  use string,        only: to_str
  use tally_header,  only: TallyObject

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
    character(MAX_FILE_LEN) :: filename
    type(TallyObject), pointer :: t => null()

# ifdef HDF5
    integer(HSIZE_T)     :: dims(1)          ! dimensions of 1D arrays
    integer(HSIZE_T)     :: dims2(2)         ! dimensions of 2D arrays
    integer(HID_T)       :: dspace           ! identifier for dataspace
    integer(HID_T)       :: dset             ! identifier for dataset
    integer(HID_T)       :: tallies_group    ! "tallies" group
    integer(HID_T)       :: temp_group       ! group for i-th tally or mesh
    type(c_ptr)          :: f_ptr            ! Pointer for h5dwrite
#endif

    ! Set filename for state point
#ifdef HDF5
    filename = trim(path_output) // 'statepoint.' // &
         trim(to_str(current_batch)) // '.h5'
#else
    filename = trim(path_output) // 'statepoint.' // &
         trim(to_str(current_batch)) // '.binary'
#endif

    ! Write message
    message = "Creating state point " // trim(filename) // "..."
    call write_message(1)

    ! Create statepoint file 
#   ifdef HDF5
      if (master) call hdf5_file_create(filename, hdf5_state_point) 
#   elif MPI
      call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, MPI_MODE_CREATE + &
           MPI_MODE_WRONLY, MPI_INFO_NULL, fh, mpi_err)
#   else
      open(UNIT=UNIT_STATE, FILE=filename, STATUS='replace', &
           ACCESS='stream')
#   endif

    ! ==========================================================================
    ! RUN INFORMATION AND TALLY METADATA

    if (master) call write_state_point_header(fh)

    ! ==========================================================================
    ! TALLY RESULTS

    if (.not. reduce_tallies) then
      ! If using the no-tally-reduction method, we need to collect tally
      ! results before writing them to the state point file.

      call write_tally_results_nr(fh)

    elseif (master) then
      ! Write number of realizations
#     ifdef HDF5
        call hdf5_write_integer(hdf5_state_point, "n_realizations", &
             n_realizations)
#     elif MPI
        call MPI_FILE_WRITE(fh, n_realizations, 1, MPI_INTEGER, &
             MPI_STATUS_IGNORE, mpi_err)
#     else
        write(UNIT_STATE) n_realizations
#     endif

      ! Write global tallies
#     ifdef HDF5
        call hdf5_write_integer(hdf5_state_point, "n_global_tallies", &
             N_GLOBAL_TALLIES)
        dims(1) = N_GLOBAL_TALLIES
        call h5screate_simple_f(1, dims, dspace, hdf5_err)
        call h5dcreate_f(hdf5_state_point, "global_tallies", &
             hdf5_tallyresult_t, dspace, dset, hdf5_err)
        f_ptr = c_loc(global_tallies(1))
        call h5dwrite_f(dset, hdf5_tallyresult_t, f_ptr, hdf5_err)
        call h5dclose_f(dset, hdf5_err)
        call h5sclose_f(dspace, hdf5_err)
#     elif MPI
        call MPI_FILE_WRITE(fh, N_GLOBAL_TALLIES, 1, MPI_INTEGER, &
             MPI_STATUS_IGNORE, mpi_err)
        call MPI_FILE_WRITE(fh, global_tallies, N_GLOBAL_TALLIES, &
             MPI_TALLYRESULT, MPI_STATUS_IGNORE, mpi_err)
#     else
        write(UNIT_STATE) N_GLOBAL_TALLIES
        GLOBAL_TALLIES_LOOP: do i = 1, N_GLOBAL_TALLIES
          write(UNIT_STATE) global_tallies(i) % sum
          write(UNIT_STATE) global_tallies(i) % sum_sq
        end do GLOBAL_TALLIES_LOOP
#     endif

      ! Write tallies
      if (tallies_on) then

#ifdef HDF5
        ! Open group "tallies"
        call h5gopen_f(hdf5_state_point, "tallies", tallies_group, hdf5_err)
#endif

        ! Indicate that tallies are on
#       ifdef HDF5
          call hdf5_write_integer(tallies_group, "tallies_present", 1)
#       elif MPI
          temp = 1
          call MPI_FILE_WRITE(fh, temp, 1, MPI_INTEGER, &
               MPI_STATUS_IGNORE, mpi_err)
#       else
          write(UNIT_STATE) 1
#       endif

        ! Write all tally results
        TALLY_RESULTS: do i = 1, n_tallies
          t => tallies(i)
#ifdef HDF5
          ! Open group for the i-th tally
          call h5gopen_f(tallies_group, "tally" // to_str(i), &
               temp_group, hdf5_err)
#endif

          ! Write sum and sum_sq for each bin
#         ifdef HDF5
            dims2 = shape(t % results)
            call h5screate_simple_f(2, dims2, dspace, hdf5_err)
            call h5dcreate_f(temp_group, "results", hdf5_tallyresult_t, &
                 dspace, dset, hdf5_err)
            f_ptr = c_loc(t % results(1, 1))
            call h5dwrite_f(dset, hdf5_tallyresult_t, f_ptr, hdf5_err)
            call h5dclose_f(dset, hdf5_err)
            call h5sclose_f(dspace, hdf5_err)
#         elif MPI
            n = size(t % results, 1) * size(t % results, 2)
            call MPI_FILE_WRITE(fh, t % results, n, MPI_TALLYRESULT, &
                 MPI_STATUS_IGNORE, mpi_err)
#         else
            do k = 1, size(t % results, 2)
              do j = 1, size(t % results, 1)
                write(UNIT_STATE) t % results(j,k) % sum
                write(UNIT_STATE) t % results(j,k) % sum_sq
              end do
            end do
#         endif

#ifdef HDF5
          ! Close group for the i-th tally
          call h5gclose_f(temp_group, hdf5_err)
#endif
        end do TALLY_RESULTS

#ifdef HDF5
        ! Close the tallies group
        call h5gclose_f(tallies_group, hdf5_err)
#endif

      else
        ! Indicate that tallies are off
#       ifdef HDF5
          call h5gopen_f(hdf5_state_point, "tallies", tallies_group, hdf5_err)
          call hdf5_write_integer(tallies_group, "tallies_present", 0)
          call h5gclose_f(tallies_group, hdf5_err)
#       elif MPI
          temp = 0
          call MPI_FILE_WRITE(fh, temp, 1, MPI_INTEGER, &
               MPI_STATUS_IGNORE, mpi_err)
#       else
          write(UNIT_STATE) 0
#       endif
      end if
    end if

    ! Append more output data in statepoint file here, source last

    ! Write out source bank
    if (run_mode == MODE_EIGENVALUE) then
      if (source_write) then
#       ifdef HDF5
          path_source = "source." // trim(to_str(current_batch)) // ".h5"
#       else
          path_source = "source." // trim(to_str(current_batch)) // ".binary"
#       endif
        call write_source()
      end if
    endif    

    ! Close statepoint file if it hasn't already been closed
    if (.not. source_separate) then 
#     ifdef HDF5
        if (master) call h5fclose_f(hdf5_state_point, hdf5_err)
#     elif MPI
        call MPI_FILE_CLOSE(fh, mpi_err)
#     else 
        close(UNIT_STATE)
#     endif
    end if

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

#ifdef HDF5
    integer(HSIZE_T)     :: dims(1)          ! dimensions of 1D arrays
    integer(HID_T)       :: tallies_group    ! "tallies" group
    integer(HID_T)       :: temp_group       ! group for i-th tally or mesh
    integer(HID_T)       :: filter_group     ! group for i-th filter
    integer, allocatable :: temp_array(:)    ! nuclide bin array
#endif

    ! Write revision number for state point file
#   ifdef HDF5
      call hdf5_write_integer(hdf5_state_point, "revision_statepoint", &
           REVISION_STATEPOINT)
#   elif MPI
      call MPI_FILE_WRITE(fh, REVISION_STATEPOINT, 1, MPI_INTEGER, &
           MPI_STATUS_IGNORE, mpi_err)
#   else
      write(UNIT_STATE) REVISION_STATEPOINT
#   endif

    ! Write OpenMC version
#   ifdef HDF5
      call hdf5_write_integer(hdf5_state_point, "version_major", &
           VERSION_MAJOR)
      call hdf5_write_integer(hdf5_state_point, "version_minor", &
           VERSION_MINOR)
      call hdf5_write_integer(hdf5_state_point, "version_release", &
           VERSION_RELEASE)
#   elif MPI
      call MPI_FILE_WRITE(fh, VERSION_MAJOR, 1, MPI_INTEGER, &
           MPI_STATUS_IGNORE, mpi_err)
      call MPI_FILE_WRITE(fh, VERSION_MINOR, 1, MPI_INTEGER, &
           MPI_STATUS_IGNORE, mpi_err)
      call MPI_FILE_WRITE(fh, VERSION_RELEASE, 1, MPI_INTEGER, &
           MPI_STATUS_IGNORE, mpi_err)
#   else
      write(UNIT_STATE) VERSION_MAJOR, VERSION_MINOR, VERSION_RELEASE
#   endif

    ! Write current date and time
#   ifdef HDF5
      call h5ltmake_dataset_string_f(hdf5_state_point, "date_and_time", &
           time_stamp(), hdf5_err)
#   elif MPI
      call MPI_FILE_WRITE(fh, time_stamp(), 19, MPI_CHARACTER, &
           MPI_STATUS_IGNORE, mpi_err)
#   else
      write(UNIT_STATE) time_stamp()
#   endif

    ! Write path to input
#   ifdef HDF5
      call h5ltmake_dataset_string_f(hdf5_state_point, "path", &
           path_input, hdf5_err)
#   elif MPI
      call MPI_FILE_WRITE(fh, path_input, MAX_FILE_LEN, MPI_CHARACTER, &
           MPI_STATUS_IGNORE, mpi_err)
#   else
      write(UNIT_STATE) path_input
#   endif

    ! Write out random number seed
#   ifdef HDF5
      call hdf5_write_long(hdf5_state_point, "seed", seed)
#   elif MPI 
      call MPI_FILE_WRITE(fh, seed, 1, MPI_INTEGER8, &
           MPI_STATUS_IGNORE, mpi_err)
#   else
      write(UNIT_STATE) seed
#   endif

    ! Write run information
#   ifdef HDF5
      call hdf5_write_integer(hdf5_state_point, "run_mode", run_mode)
      call hdf5_write_long(hdf5_state_point, "n_particles", n_particles)
      call hdf5_write_integer(hdf5_state_point, "n_batches", n_batches)
#   elif MPI
      call MPI_FILE_WRITE(fh, run_mode, 1, MPI_INTEGER, &
           MPI_STATUS_IGNORE, mpi_err)
      call MPI_FILE_WRITE(fh, n_particles, 1, MPI_INTEGER8, &
           MPI_STATUS_IGNORE, mpi_err)
      call MPI_FILE_WRITE(fh, n_batches, 1, MPI_INTEGER, &
           MPI_STATUS_IGNORE, mpi_err)
#   else
      write(UNIT_STATE) run_mode, n_particles, n_batches
#   endif

    ! Write out current batch number
#   ifdef HDF5
      call hdf5_write_integer(hdf5_state_point, "current_batch", &
           current_batch)
#   elif MPI
      call MPI_FILE_WRITE(fh, current_batch, 1, MPI_INTEGER, &
           MPI_STATUS_IGNORE, mpi_err)
#   else
      write(UNIT_STATE) current_batch
#   endif

    ! Write out information for eigenvalue run
    if (run_mode == MODE_EIGENVALUE) then
#     ifdef HDF5
        call hdf5_write_integer(hdf5_state_point, "n_inactive", &
             n_inactive)
        call hdf5_write_integer(hdf5_state_point, "gen_per_batch", &
             gen_per_batch)
        dims(1) = current_batch
        call h5ltmake_dataset_double_f(hdf5_state_point, "k_batch", 1, &
             dims, k_batch, hdf5_err)
        dims(1) = current_batch*gen_per_batch
        call h5ltmake_dataset_double_f(hdf5_state_point, "entropy", 1, &
             dims, entropy, hdf5_err)
        call hdf5_write_double(hdf5_state_point, "k_col_abs", k_col_abs)
        call hdf5_write_double(hdf5_state_point, "k_col_tra", k_col_tra)
        call hdf5_write_double(hdf5_state_point, "k_abs_tra", k_abs_tra)
        dims(1) = 2
        call h5ltmake_dataset_double_f(hdf5_state_point, "k_combined", 1, &
             dims, k_combined, hdf5_err)
#     elif MPI
        call MPI_FILE_WRITE(fh, n_inactive, 1, MPI_INTEGER, &
             MPI_STATUS_IGNORE, mpi_err)
        call MPI_FILE_WRITE(fh, gen_per_batch, 1, MPI_INTEGER, &
             MPI_STATUS_IGNORE, mpi_err)
        call MPI_FILE_WRITE(fh, k_batch, current_batch, MPI_REAL8, &
             MPI_STATUS_IGNORE, mpi_err)
        call MPI_FILE_WRITE(fh, entropy, current_batch*gen_per_batch, &
             MPI_REAL8, MPI_STATUS_IGNORE, mpi_err)
        call MPI_FILE_WRITE(fh, k_col_abs, 1, MPI_REAL8, &
             MPI_STATUS_IGNORE, mpi_err)
        call MPI_FILE_WRITE(fh, k_col_tra, 1, MPI_REAL8, &
             MPI_STATUS_IGNORE, mpi_err)
        call MPI_FILE_WRITE(fh, k_abs_tra, 1, MPI_REAL8, &
             MPI_STATUS_IGNORE, mpi_err)
        call MPI_FILE_WRITE(fh, k_combined, 2, MPI_REAL8, &
             MPI_STATUS_IGNORE, mpi_err)
#     else
        write(UNIT_STATE) n_inactive, gen_per_batch
        write(UNIT_STATE) k_batch(1:current_batch)
        write(UNIT_STATE) entropy(1:current_batch*gen_per_batch)
        write(UNIT_STATE) k_col_abs
        write(UNIT_STATE) k_col_tra
        write(UNIT_STATE) k_abs_tra
        write(UNIT_STATE) k_combined
#     endif
    end if

#ifdef HDF5
    ! Create group "tallies"
    call h5gcreate_f(hdf5_state_point, "tallies", tallies_group, hdf5_err)
#endif

    ! Write number of meshes
#   ifdef HDF5
      call hdf5_write_integer(tallies_group, "n_meshes", n_meshes)
#   elif MPI
      call MPI_FILE_WRITE(fh, n_meshes, 1, MPI_INTEGER, &
           MPI_STATUS_IGNORE, mpi_err)
#   else
      write(UNIT_STATE) n_meshes
#   endif

    ! Write information for meshes
    MESH_LOOP: do i = 1, n_meshes
#     ifdef HDF5
        ! Create temporary group for each mesh
        call h5gcreate_f(tallies_group, "mesh" // to_str(i), &
             temp_group, hdf5_err)

        call hdf5_write_integer(temp_group, "id", meshes(i) % id)
        call hdf5_write_integer(temp_group, "type", meshes(i) % type)
        call hdf5_write_integer(temp_group, "n_dimension", &
             meshes(i) % n_dimension)
        dims(1) = meshes(i) % n_dimension
        call hdf5_write_array(temp_group, "dimension", &
             meshes(i) % dimension, 1, dims)
        call hdf5_write_array(temp_group, "lower_left", &
             meshes(i) % lower_left, 1, dims)
        call hdf5_write_array(temp_group, "upper_right", &
             meshes(i) % upper_right, 1, dims)
        call hdf5_write_array(temp_group, "width", &
             meshes(i) % width, 1, dims)

        ! Close temporary group for mesh
        call h5gclose_f(temp_group, hdf5_err)
#     elif MPI
        call MPI_FILE_WRITE(fh, meshes(i) % id, 1, MPI_INTEGER, &
             MPI_STATUS_IGNORE, mpi_err)
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
#     else
        write(UNIT_STATE) meshes(i) % id
        write(UNIT_STATE) meshes(i) % type
        write(UNIT_STATE) meshes(i) % n_dimension
        write(UNIT_STATE) meshes(i) % dimension
        write(UNIT_STATE) meshes(i) % lower_left
        write(UNIT_STATE) meshes(i) % upper_right
        write(UNIT_STATE) meshes(i) % width
#     endif
    end do MESH_LOOP

    ! Write number of tallies
#   ifdef HDF5
      call hdf5_write_integer(tallies_group, "n_tallies", n_tallies)
#   elif MPI
      call MPI_FILE_WRITE(fh, n_tallies, 1, MPI_INTEGER, &
           MPI_STATUS_IGNORE, mpi_err)
#   else
      write(UNIT_STATE) n_tallies
#   endif

    TALLY_METADATA: do i = 1, n_tallies
      ! Get pointer to tally
      t => tallies(i)

#ifdef HDF5
      ! Create group for this tally
      call h5gcreate_f(tallies_group, "tally" // to_str(i), &
           temp_group, hdf5_err)
#endif

      ! Write id
#     ifdef HDF5
        call hdf5_write_integer(temp_group, "id", t % id)
#     elif MPI
        call MPI_FILE_WRITE(fh, t % id, 1, MPI_INTEGER, &
             MPI_STATUS_IGNORE, mpi_err)
#     else
        write(UNIT_STATE) t % id
#     endif

      ! Write number of realizations
#     ifdef HDF5
        call hdf5_write_integer(temp_group, "n_realizations", &
             t % n_realizations)
#     elif HDf5
        call MPI_FILE_WRITE(fh, t % n_realizations, 1, MPI_INTEGER, &
             MPI_STATUS_IGNORE, mpi_err)
#     else
        write(UNIT_STATE) t % n_realizations
#     endif

      ! Write size of each tally
#     ifdef HDF5
        call hdf5_write_integer(temp_group, "total_score_bins", &
             t % total_score_bins)
        call hdf5_write_integer(temp_group, "total_filter_bins", &
             t % total_filter_bins)
#     elif MPI
        call MPI_FILE_WRITE(fh, t % total_score_bins, 1, MPI_INTEGER, &
             MPI_STATUS_IGNORE, mpi_err)
        call MPI_FILE_WRITE(fh, t % total_filter_bins, 1, MPI_INTEGER, &
             MPI_STATUS_IGNORE, mpi_err)
#     else
        write(UNIT_STATE) t % total_score_bins
        write(UNIT_STATE) t % total_filter_bins
#     endif

      ! Write number of filters
#     ifdef HDF5
        call hdf5_write_integer(temp_group, "n_filters", t % n_filters)
#     elif MPI
        call MPI_FILE_WRITE(fh, t % n_filters, 1, MPI_INTEGER, &
             MPI_STATUS_IGNORE, mpi_err)
#     else
        write(UNIT_STATE) t % n_filters
#     endif

      FILTER_LOOP: do j = 1, t % n_filters
#ifdef HDF5
        ! Create filter group
        call h5gcreate_f(temp_group, "filter" // to_str(j), filter_group, &
             hdf5_err)
#endif

        ! Write type of filter
#       ifdef HDF5
          call hdf5_write_integer(filter_group, "type", &
               t % filters(j) % type)
#       elif MPI
          call MPI_FILE_WRITE(fh, t % filters(j) % type, 1, MPI_INTEGER, &
               MPI_STATUS_IGNORE, mpi_err)
#       else
          write(UNIT_STATE) t % filters(j) % type
#       endif

        ! Write number of bins for this filter
#       ifdef HDF5
          call hdf5_write_integer(filter_group, "n_bins", &
               t % filters(j) % n_bins)
#       elif MPI
          call MPI_FILE_WRITE(fh, t % filters(j) % n_bins, &
               1, MPI_INTEGER, MPI_STATUS_IGNORE, mpi_err)
#       else
          write(UNIT_STATE) t % filters(j) % n_bins
#       endif

        ! Write bins
        if (t % filters(j) % type == FILTER_ENERGYIN .or. &
             t % filters(j) % type == FILTER_ENERGYOUT) then
#         ifdef HDF5
            dims(1) = size(t % filters(j) % real_bins)
            call h5ltmake_dataset_double_f(filter_group, "bins", 1, &
                 dims, t % filters(j) % real_bins, hdf5_err)
#         elif MPI
            n = size(t % filters(j) % real_bins)
            call MPI_FILE_WRITE(fh, t % filters(j) % real_bins, n, &
                 MPI_REAL8, MPI_STATUS_IGNORE, mpi_err)
#         else
            write(UNIT_STATE) t % filters(j) % real_bins
#         endif
        else
#         ifdef HDF5
            dims(1) = size(t % filters(j) % int_bins)
            call h5ltmake_dataset_int_f(filter_group, "bins", 1, &
                 dims, t % filters(j) % int_bins, hdf5_err)
#         elif MPI
            n = size(t % filters(j) % int_bins)
            call MPI_FILE_WRITE(fh, t % filters(j) % int_bins, n, &
                 MPI_INTEGER, MPI_STATUS_IGNORE, mpi_err)
#         else
            write(UNIT_STATE) t % filters(j) % int_bins
#         endif
        end if

#ifdef HDF5
        ! Write name of type
        select case (t % filters(j) % type)
        case(FILTER_UNIVERSE)
          call h5ltmake_dataset_string_f(filter_group, "type_name", &
               "universe", hdf5_err)
        case(FILTER_MATERIAL)
          call h5ltmake_dataset_string_f(filter_group, "type_name", &
               "material", hdf5_err)
        case(FILTER_CELL)
          call h5ltmake_dataset_string_f(filter_group, "type_name", &
               "cell", hdf5_err)
        case(FILTER_CELLBORN)
          call h5ltmake_dataset_string_f(filter_group, "type_name", &
               "cellborn", hdf5_err)
        case(FILTER_SURFACE)
          call h5ltmake_dataset_string_f(filter_group, "type_name", &
               "surface", hdf5_err)
        case(FILTER_MESH)
          call h5ltmake_dataset_string_f(filter_group, "type_name", &
               "mesh", hdf5_err)
        case(FILTER_ENERGYIN)
          call h5ltmake_dataset_string_f(filter_group, "type_name", &
               "energy", hdf5_err)
        case(FILTER_ENERGYOUT)
          call h5ltmake_dataset_string_f(filter_group, "type_name", &
               "energyout", hdf5_err)
        end select

        ! Close group for this filter
        call h5gclose_f(filter_group, hdf5_err)
#endif
      end do FILTER_LOOP

      ! Write number of nuclide bins
#     ifdef HDF5
        call hdf5_write_integer(temp_group, "n_nuclide_bins", &
             t % n_nuclide_bins)
#     elif MPI
        call MPI_FILE_WRITE(fh, t % n_nuclide_bins, 1, MPI_INTEGER, &
             MPI_STATUS_IGNORE, mpi_err)
#     else
        write(UNIT_STATE) t % n_nuclide_bins
#     endif

#ifdef HDF5
      ! Allocate array for HDf5
      allocate(temp_array(t % n_nuclide_bins))
#endif

      ! Write nuclide bins
      NUCLIDE_LOOP: do j = 1, t % n_nuclide_bins
        if (t % nuclide_bins(j) > 0) then
#         ifdef HDF5
            temp_array(j) = nuclides(t % nuclide_bins(j)) % zaid
#         elif MPI
            call MPI_FILE_WRITE(fh, nuclides(t % nuclide_bins(j)) % zaid, &
                 1, MPI_INTEGER, MPI_STATUS_IGNORE, mpi_err)
#         else
            write(UNIT_STATE) nuclides(t % nuclide_bins(j)) % zaid
#         endif
        else
#         ifdef HDF5
            temp_array(j) = t % nuclide_bins(j)
#         elif MPI
            call MPI_FILE_WRITE(fh, t % nuclide_bins(j), 1, &
                 MPI_INTEGER, MPI_STATUS_IGNORE, mpi_err)
#         else
            write(UNIT_STATE) t % nuclide_bins(j)
#         endif
        end if
      end do NUCLIDE_LOOP

#ifdef HDF5
      ! Write and deallocate nuclide bins
      dims(1) = t % n_nuclide_bins
      call h5ltmake_dataset_int_f(temp_group, "nuclide_bins", 1, &
           dims, temp_array, hdf5_err)
      deallocate(temp_array)
#endif

      ! Write number of score bins, score bins, and scatt order
#     ifdef HDF5
        call hdf5_write_integer(temp_group, "n_score_bins", &
             t % n_score_bins)
        dims(1) = t % n_score_bins
        call h5ltmake_dataset_int_f(temp_group, "score_bins", 1, &
             dims, t % score_bins, hdf5_err)
        call h5ltmake_dataset_int_f(temp_group, "scatt_order", 1, &
             dims, t % scatt_order, hdf5_err)
#     elif MPI
        call MPI_FILE_WRITE(fh, t % n_score_bins, 1, MPI_INTEGER, &
             MPI_STATUS_IGNORE, mpi_err)
        call MPI_FILE_WRITE(fh, t % score_bins, t % n_score_bins, &
             MPI_INTEGER, MPI_STATUS_IGNORE, mpi_err)
        call MPI_FILE_WRITE(fh, t % scatt_order, t % n_score_bins, &
             MPI_INTEGER, MPI_STATUS_IGNORE, mpi_err)
#     else
        write(UNIT_STATE) t % n_score_bins
        write(UNIT_STATE) t % score_bins
        write(UNIT_STATE) t % scatt_order
#     endif

      ! Write number of user score bins
#     ifdef HDF5
        call hdf5_write_integer(temp_group, "n_user_score_bins", &
             t % n_user_score_bins)
#     elif MPI
        call MPI_FILE_WRITE(fh, t % n_user_score_bins, 1, MPI_INTEGER, &
             MPI_STATUS_IGNORE, mpi_err)
#     else
        write(UNIT_STATE) t % n_user_score_bins
#     endif

#ifdef HDF5
      ! Close tally group
      call h5gclose_f(temp_group, hdf5_err)
#endif
    end do TALLY_METADATA

# ifdef HDF5
    ! Close tallies group
    call h5gclose_f(tallies_group, hdf5_err)
#endif

  end subroutine write_state_point_header
#endif

#ifdef MPI
!===============================================================================
! WRITE_TALLY_RESULTS_NR
!===============================================================================

  subroutine write_tally_results_nr(fh)

    integer, intent(in) :: fh ! file handle

    integer :: i      ! loop index
    integer :: n      ! number of filter bins
    integer :: m      ! number of score bins
    integer :: temp   ! temporary variable
    integer :: n_bins ! total number of bins
    real(8), allocatable :: tally_temp(:,:,:) ! contiguous array of results
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
          call MPI_REDUCE(MPI_IN_PLACE, tally_temp, n_bins, MPI_REAL8, &
               MPI_SUM, 0, MPI_COMM_WORLD, mpi_err)

          ! Write reduced tally results to file
          call MPI_FILE_WRITE(fh, tally_temp, n_bins, MPI_REAL8, &
               MPI_STATUS_IGNORE, mpi_err)

          ! At the end of the simulation, store the results back in the
          ! regular TallyResults array
          if (current_batch == n_batches) then
            t % results(:,:) % sum = tally_temp(1,:,:)
            t % results(:,:) % sum_sq = tally_temp(2,:,:)
          end if
        else
          ! Receive buffer not significant at other processors
          call MPI_REDUCE(tally_temp, dummy, n_bins, MPI_REAL8, MPI_SUM, &
               0, MPI_COMM_WORLD, mpi_err)
        end if

        ! Deallocate temporary copy of tally results
        deallocate(tally_temp)
      end do TALLY_RESULTS
    else
      if (master) then
        ! Indicate that tallies are off
        temp = 0
        call MPI_FILE_WRITE(fh, temp, 1, MPI_INTEGER, &
             MPI_STATUS_IGNORE, mpi_err)
      end if
    end if

  end subroutine write_tally_results_nr
#endif

!===============================================================================
! LOAD_STATE_POINT loads data from a state point file to either continue a run
! or to print intermediate tally results
!===============================================================================

  subroutine load_state_point()

    integer :: i, j           ! loop indices
    integer :: mode           ! specified run mode
    integer :: temp(3)        ! temporary variable
    integer, allocatable :: int_array(:)
    real(8), allocatable :: real_array(:)
    character(19)        :: current_time  ! current date and time
    character(MAX_FILE_LEN) :: path_temp

#ifdef MPI
    integer :: fh                      ! file handle
    integer :: n                       ! temporary array size
    integer :: size_offset_kind        ! size of MPI_OFFSET_KIND (bytes)
    integer :: size_bank               ! size of MPI_BANK type
    integer(MPI_OFFSET_KIND) :: offset ! offset in memory (0=beginning of file)
#else
    integer :: k ! loop index
#endif

#ifdef HDF5
    integer(HID_T)   :: hdf5_state_point ! identifier for state point file
    integer(HID_T)   :: tally_group      ! identifier for tally group
    integer(HID_T)   :: tallies_group
    integer(HID_T)   :: temp_group
    integer(HID_T)   :: filter_group
    integer(HID_T)   :: dset             ! identifier for dataset
    integer(HSIZE_T) :: dims(1)          ! dimensions of 1D arrays
    integer(HID_T) :: hdf5_source
    integer(HID_T) :: filespace
    integer(HID_T) :: memspace
    character(6) :: dsetname = 'source'
    integer(HID_T) :: dset_id
    integer(HID_T) :: plist_id
    type(c_ptr) :: f_ptr
    integer(HSIZE_T) :: maxdims(1)
#endif

    ! Write message
    message = "Loading state point " // trim(path_state_point) // "..."
    call write_message(1)

    ! Open binary state point file for reading
#   ifdef HDF5
      call h5fopen_f(path_state_point, H5F_ACC_RDONLY_F, &
           hdf5_state_point, hdf5_err)
#   elif MPI
      call MPI_FILE_OPEN(MPI_COMM_WORLD, path_state_point, MPI_MODE_RDONLY, &
           MPI_INFO_NULL, fh, mpi_err)
#   else
      open(UNIT=UNIT_STATE, FILE=path_state_point, STATUS='old', &
           ACCESS='stream')
#   endif

    ! ==========================================================================
    ! RUN INFORMATION AND TALLY METADATA

    ! Read revision number for state point file and make sure it matches with
    ! current version
#   ifdef HDF5
      call hdf5_read_integer(hdf5_state_point, "revision_statepoint", temp(1))
#   elif MPI
      call MPI_FILE_READ_ALL(fh, temp, 1, MPI_INTEGER, MPI_STATUS_IGNORE, &
           mpi_err)
#   else
      read(UNIT_STATE) temp(1)
#   endif
    if (temp(1) /= REVISION_STATEPOINT) then
      message = "State point binary version does not match current version " &
           // "in OpenMC."
      call fatal_error()
    end if

    ! Read OpenMC version
#   ifdef HDF5
      call hdf5_read_integer(hdf5_state_point, "version_major", temp(1))
      call hdf5_read_integer(hdf5_state_point, "version_minor", temp(2))
      call hdf5_read_integer(hdf5_state_point, "version_release", temp(3))
#   elif MPI
      call MPI_FILE_READ_ALL(fh, temp, 3, MPI_INTEGER, MPI_STATUS_IGNORE, &
           mpi_err)
#   else
      read(UNIT_STATE) temp(1:3)
#   endif
    if (temp(1) /= VERSION_MAJOR .or. temp(2) /= VERSION_MINOR &
         .or. temp(3) /= VERSION_RELEASE) then
      message = "State point file was created with a different version " // &
           "of OpenMC."
      call warning()
    end if

    ! Read date and time
#   ifdef HDF5
      call h5ltread_dataset_string_f(hdf5_state_point, "date_and_time", &
           current_time, hdf5_err) ! TODO check this
#   elif MPI
      call MPI_FILE_READ_ALL(fh, current_time, 19, MPI_CHARACTER, &
           MPI_STATUS_IGNORE, mpi_err)
#   else
      read(UNIT_STATE) current_time
#   endif

    ! Read path to input
#   ifdef HDF5
      call h5ltread_dataset_string_f(hdf5_state_point, "path", &
           path_temp, hdf5_err)
#   elif MPI
      call MPI_FILE_READ_ALL(fh, path_temp, MAX_FILE_LEN, MPI_CHARACTER, &
           MPI_STATUS_IGNORE, mpi_err)
#   else
      read(UNIT_STATE) path_temp
#   endif

    ! Read and overwrite random number seed
#   ifdef HDF5
      call hdf5_read_long(hdf5_state_point, "seed", seed)
#   elif MPI
      call MPI_FILE_READ_ALL(fh, seed, 1, MPI_INTEGER8, &
           MPI_STATUS_IGNORE, mpi_err)
#   else
      read(UNIT_STATE) seed
#   endif

    ! Read and overwrite run information except for n_batches
#   ifdef HDF5
      call hdf5_read_integer(hdf5_state_point, "run_mode", mode)
      call hdf5_read_long(hdf5_state_point, "n_particles", n_particles)
      call hdf5_read_integer(hdf5_state_point, "n_batches", temp(1))
#   elif MPI
      call MPI_FILE_READ_ALL(fh, mode, 1, MPI_INTEGER, &
           MPI_STATUS_IGNORE, mpi_err)
      call MPI_FILE_READ_ALL(fh, n_particles, 1, MPI_INTEGER8, &
           MPI_STATUS_IGNORE, mpi_err)
      call MPI_FILE_READ_ALL(fh, temp(1), 1, MPI_INTEGER, &
           MPI_STATUS_IGNORE, mpi_err)
#   else
      read(UNIT_STATE) mode, n_particles, temp(1)
#   endif

    ! Allow user to specify more than n_batches
    n_batches = max(n_batches, temp(1))

    ! Read batch number to restart at
#   ifdef HDF5
      call hdf5_read_integer(hdf5_state_point, "current_batch", restart_batch)
#   elif MPI
      call MPI_FILE_READ_ALL(fh, restart_batch, 1, MPI_INTEGER, &
           MPI_STATUS_IGNORE, mpi_err)
#   else
      read(UNIT_STATE) restart_batch
#   endif

    ! Read information specific to eigenvalue run
    if (mode == MODE_EIGENVALUE) then
#     ifdef HDF5
        call hdf5_read_integer(hdf5_state_point, "n_inactive", temp(1))
        call hdf5_read_integer(hdf5_state_point, "gen_per_batch", &
             gen_per_batch)
        dims(1) = restart_batch
        call h5ltread_dataset_double_f(hdf5_state_point, "k_batch", &
             k_batch(1:restart_batch), dims, hdf5_err)
        dims(1) = restart_batch*gen_per_batch
        call h5ltread_dataset_double_f(hdf5_state_point, "entropy", &
             entropy(1:restart_batch*gen_per_batch), dims, hdf5_err)
        call hdf5_read_double(hdf5_state_point, "k_col_abs", k_col_abs)
        call hdf5_read_double(hdf5_state_point, "k_col_tra", k_col_tra)
        call hdf5_read_double(hdf5_state_point, "k_abs_tra", k_abs_tra)
        ! not reading in k-combined because below code doesnt
#     elif MPI
        call MPI_FILE_READ_ALL(fh, temp(1), 1, MPI_INTEGER, &
             MPI_STATUS_IGNORE, mpi_err)
        call MPI_FILE_READ_ALL(fh, gen_per_batch, 1, MPI_INTEGER, &
             MPI_STATUS_IGNORE, mpi_err)
        call MPI_FILE_READ_ALL(fh, k_batch, restart_batch, MPI_REAL8, &
             MPI_STATUS_IGNORE, mpi_err)
        call MPI_FILE_READ_ALL(fh, entropy, restart_batch*gen_per_batch, &
             MPI_REAL8, MPI_STATUS_IGNORE, mpi_err)
        call MPI_FILE_READ_ALL(fh, k_col_abs, 1, MPI_REAL8, &
             MPI_STATUS_IGNORE, mpi_err)
        call MPI_FILE_READ_ALL(fh, k_col_tra, 1, MPI_REAL8, &
             MPI_STATUS_IGNORE, mpi_err)
        call MPI_FILE_READ_ALL(fh, k_abs_tra, 1, MPI_REAL8, &
             MPI_STATUS_IGNORE, mpi_err)
        allocate(real_array(2))
        call MPI_FILE_READ_ALL(fh, real_array, 2, MPI_REAL8, &
             MPI_STATUS_IGNORE, mpi_err)
        deallocate(real_array)
#     else
        read(UNIT_STATE) temp(1), gen_per_batch
        read(UNIT_STATE) k_batch(1:restart_batch)
        read(UNIT_STATE) entropy(1:restart_batch*gen_per_batch)
        read(UNIT_STATE) k_col_abs
        read(UNIT_STATE) k_col_tra
        read(UNIT_STATE) k_abs_tra
        allocate(real_array(2))
        read(UNIT_STATE) real_array
        deallocate(real_array)
#     endif

      ! Allow user to modify n_inactive
      n_inactive = max(n_inactive, temp(1))

    end if

    if (master) then
      ! Read number of meshes
#     ifdef HDF5
        call h5gopen_f(hdf5_state_point, "tallies", tallies_group, hdf5_err) ! TODO Close group
        call hdf5_read_integer(tallies_group, "n_meshes", temp(1))
#     elif MPI
        call MPI_FILE_READ(fh, temp, 1, MPI_INTEGER, MPI_STATUS_IGNORE, &
             mpi_err)
#     else
        read(UNIT_STATE) temp(1)
#     endif
      if (temp(1) /= n_meshes) then
        message = "Number of meshes does not match in state point."
        call fatal_error()
      end if

      MESH_LOOP: do i = 1, n_meshes
        ! Read id, mesh type, and dimension
#       ifdef HDF5
          ! Nothing performed for HDF5, skip reading of mesh
#       elif MPI
          call MPI_FILE_READ(fh, temp, 3, MPI_INTEGER, MPI_STATUS_IGNORE, &
               mpi_err)

          ! Skip mesh data
          call MPI_FILE_GET_POSITION(fh, offset, mpi_err)
          offset = offset + temp(3)*(4 + 3*8)
          call MPI_FILE_SEEK(fh, offset, MPI_SEEK_SET, mpi_err)
#       else
          ! Read id, mesh type, and dimension
          read(UNIT_STATE) temp(1:3)

          ! Allocate temporary arrays
          allocate(int_array(temp(3)))
          allocate(real_array(temp(3)))

          ! Read dimension, lower_left, upper_right, width
          read(UNIT_STATE) int_array
          read(UNIT_STATE) real_array
          read(UNIT_STATE) real_array
          read(UNIT_STATE) real_array

          ! Deallocate temporary arrays
          deallocate(int_array)
          deallocate(real_array)
#       endif
      end do MESH_LOOP

      ! Read number of tallies and make sure it matches
#     ifdef HDF5
        call hdf5_read_integer(tallies_group, "n_tallies", temp(1))
#     elif MPI
        call MPI_FILE_READ(fh, temp, 1, MPI_INTEGER, MPI_STATUS_IGNORE, &
             mpi_err)
#     else
        read(UNIT_STATE) temp(1)
#     endif
      if (temp(1) /= n_tallies) then
        message = "Number of tallies does not match in state point."
        call fatal_error()
      end if
!call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)
!call MPI_ABORT(MPI_COMM_WORLD, -8, mpi_err)

      TALLY_METADATA: do i = 1, n_tallies

        ! Read tally id
#       ifdef HDF5
          call h5gopen_f(tallies_group, "tally" // to_str(i), &
                         temp_group, hdf5_err) 
          call hdf5_read_integer(temp_group, "id", tallies(i) % id)
#       elif MPI
          call MPI_FILE_READ(fh, tallies(i) % id, 1, MPI_INTEGER, &
               MPI_STATUS_IGNORE, mpi_err)
#       else
          read(UNIT_STATE) tallies(i) % id
#       endif

        ! Read number of realizations for tallies
#       ifdef HDF5
          call hdf5_read_integer(temp_group, "n_realizations", &
               tallies(i) % n_realizations)
#       elif MPI 
          call MPI_FILE_READ(fh, tallies(i) % n_realizations, 1, &
               MPI_INTEGER, MPI_STATUS_IGNORE, mpi_err)
#       else
          read(UNIT_STATE) tallies(i) % n_realizations
#       endif

        ! Read dimensions of tally filters and results and make sure they
        ! match
#       ifdef HDF5
          call hdf5_read_integer(temp_group, "total_score_bins", &
               temp(1))
          call hdf5_read_integer(temp_group, "total_filter_bins", &
               temp(2))
#       elif MPI
          call MPI_FILE_READ(fh, temp, 2, MPI_INTEGER, &
               MPI_STATUS_IGNORE, mpi_err)
#       else
          read(UNIT_STATE) temp(1:2)
#       endif
        if (temp(1) /= size(tallies(i) % results, 1) .or. &
             temp(2) /= size(tallies(i) % results, 2)) then
          message = "Tally dimensions do not match in state point."
          call fatal_error()
        end if

        ! Read number of filters
#       ifdef HDF5
          call hdf5_read_integer(temp_group, "n_filters", temp(1))
#       elif MPI
          call MPI_FILE_READ(fh, temp, 1, MPI_INTEGER, &
               MPI_STATUS_IGNORE, mpi_err)
#       else
          read(UNIT_STATE) temp(1)
#       endif
!call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)
!call MPI_ABORT(MPI_COMM_WORLD, -8, mpi_err)

        FILTER_LOOP: do j = 1, temp(1)
#ifdef HDF5
          ! Open up the filter group
          call h5gopen_f(temp_group, "filter" // to_str(j), filter_group, &
               hdf5_err)
#endif

          ! Read filter type and number of bins
#         ifdef HDF5
            call hdf5_read_integer(filter_group, "type", temp(2))
            call hdf5_read_integer(filter_group, "n_bins", temp(3))
#         elif MPI
            call MPI_FILE_READ(fh, temp(2), 2, MPI_INTEGER, &
                 MPI_STATUS_IGNORE, mpi_err)
#         else
            read(UNIT_STATE) temp(2:3)
#         endif

          ! Read filter bins
          select case (temp(2))
          case (FILTER_MESH)
            allocate(int_array(1))
#           ifdef HDF5
              ! Skip HDF5 reading of this
#           elif MPI
              call MPI_FILE_READ(fh, int_array, 1, MPI_INTEGER, &
                   MPI_STATUS_IGNORE, mpi_err)
#           else
              read(UNIT_STATE) int_array
#           endif
            deallocate(int_array)
          case (FILTER_ENERGYIN, FILTER_ENERGYOUT)
            allocate(real_array(temp(3) + 1))
#           ifdef HDF5
              ! Skip HDF5 reading of this
#           elif MPI
              call MPI_FILE_READ(fh, real_array, temp(3) + 1, MPI_REAL8, &
                   MPI_STATUS_IGNORE, mpi_err)
#           else
              read(UNIT_STATE) real_array
#           endif
            deallocate(real_array)
          case default
            allocate(int_array(temp(3)))
#           ifdef HDF5
              ! Skip HDF5 reading of this
#           elif MPI
              call MPI_FILE_READ(fh, int_array, temp(3), MPI_INTEGER, &
                   MPI_STATUS_IGNORE, mpi_err)
#           else
              read(UNIT_STATE) int_array
#           endif
            deallocate(int_array)
          end select

#ifdef HDF5
            ! Close the filter group
            call h5gclose_f(filter_group, hdf5_err)
#endif

        end do FILTER_LOOP

        ! Read number of nuclides
#       ifdef HDF5
          call hdf5_read_integer(temp_group, "n_nuclide_bins", &
               temp(1))
#       elif MPI
          call MPI_FILE_READ(fh, temp, 1, MPI_INTEGER, &
               MPI_STATUS_IGNORE, mpi_err)
#       else
          read(UNIT_STATE) temp(1)
#       endif

        ! Read nuclide bins
        allocate(int_array(temp(1)))
#       ifdef HDF5
          ! Skip HDF5 for reading this
#       elif MPI
          call MPI_FILE_READ(fh, int_array, temp(1), MPI_INTEGER, &
               MPI_STATUS_IGNORE, mpi_err)
#       else
          read(UNIT_STATE) int_array
#       endif
        deallocate(int_array)

        ! Read number of score bins, score bins, and scatt_order
#       ifdef HDF5
          ! Skip HDF5 for reading this
#       elif MPI
          call MPI_FILE_READ(fh, temp, 1, MPI_INTEGER, &
               MPI_STATUS_IGNORE, mpi_err)
#       else
          read(UNIT_STATE) temp(1)
#       endif
        allocate(int_array(temp(1)))
#       ifdef HDF5
          ! Skip HDF5 for reading this
#       elif MPI
          call MPI_FILE_READ(fh, int_array, temp(1), MPI_INTEGER, &
               MPI_STATUS_IGNORE, mpi_err)
          call MPI_FILE_READ(fh, int_array, temp(1), MPI_INTEGER, &
               MPI_STATUS_IGNORE, mpi_err)
#       else
          read(UNIT_STATE) int_array
          read(UNIT_STATE) int_array
#       endif
        deallocate(int_array)
        
        ! Read number of user score bins
#       ifdef HDF5
          call hdf5_read_integer(temp_group, "n_user_score_bins", temp(1))
#       elif MPI
          call MPI_FILE_READ(fh, temp, 1, MPI_INTEGER, &
               MPI_STATUS_IGNORE, mpi_err)
#       else
          read(UNIT_STATE) temp(1)
#       endif

#ifdef HDF5
        ! Close HDF5 temp group
        call h5gclose_f(temp_group, hdf5_err)
#endif
      end do TALLY_METADATA

      ! Read number of realizations for global tallies
#     ifdef HDF5
        call hdf5_read_integer(hdf5_state_point, "n_realizations", &
             n_realizations)
#     elif MPI
        call MPI_FILE_READ(fh, n_realizations, 1, MPI_INTEGER, &
             MPI_STATUS_IGNORE, mpi_err)
#     else
        read(UNIT_STATE) n_realizations
#     endif

      ! Read number of global tallies and make sure it matches
#     ifdef HDF5
        call hdf5_read_integer(hdf5_state_point, "n_global_tallies", &
             temp(1))
#     elif MPI
        call MPI_FILE_READ(fh, temp, 1, MPI_INTEGER, MPI_STATUS_IGNORE, &
        mpi_err)
#     else
        read(UNIT_STATE) temp(1)
#     endif
      if (temp(1) /= N_GLOBAL_TALLIES) then
        message = "Number of global tallies does not match in state point."
        call fatal_error()
      end if
!call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)
!call MPI_ABORT(MPI_COMM_WORLD, -5, mpi_err)

      ! Read global tally data
#     ifdef HDF5
        ! Open global tallies dataset
        call h5dopen_f(hdf5_state_point, "global_tallies", dset, hdf5_err)

        ! Read global tallies
        f_ptr = c_loc(global_tallies(1))
        call h5dread_f(dset, hdf5_tallyresult_t, f_ptr, hdf5_err)

        ! Close global tallies dataset
        call h5dclose_f(dset, hdf5_err)
#     elif MPI
        call MPI_FILE_READ(fh, global_tallies, N_GLOBAL_TALLIES, &
             MPI_TALLYRESULT, MPI_STATUS_IGNORE, mpi_err)
#     else
        do i = 1, N_GLOBAL_TALLIES
          read(UNIT_STATE) global_tallies(i) % sum
          read(UNIT_STATE) global_tallies(i) % sum_sq
        end do
#     endif

      ! Check if tally results are present
#     ifdef HDF5
        call hdf5_read_integer(tallies_group, "tallies_present", temp(1))
#     elif MPI
        call MPI_FILE_READ(fh, temp, 1, MPI_INTEGER, MPI_STATUS_IGNORE, &
             mpi_err)
#     else
        read(UNIT_STATE) temp(1)
#     endif

      ! =======================================================================
      ! TALLY RESULTS

      ! Read sum and sum squared
      if (temp(1) == 1) then
        TALLY_RESULTS: do i = 1, n_tallies
#         ifdef HDF5
            ! Open tally group
            call h5gopen_f(tallies_group, "tally" // to_str(i), tally_group, &
                 hdf5_err)

            ! Read number of realizations
            call hdf5_read_integer(tally_group, "n_realizations", &
                 tallies(i) % n_realizations)

            ! Open dataset for tally results
            call h5dopen_f(tally_group, "results", dset, hdf5_err)

            ! Read sum and sum_sq for each tally bin
            f_ptr = c_loc(tallies(i) % results(1,1))
            call h5dread_f(dset, hdf5_tallyresult_t, f_ptr, hdf5_err)

            ! Close dataset for tally results
            call h5dclose_f(dset, hdf5_err)

            ! Close tally group
            call h5gclose_f(tally_group, hdf5_err)
#         elif MPI
            n = size(tallies(i) % results, 1) * size(tallies(i) % results, 2)
            call MPI_FILE_READ(fh, tallies(i) % results, n, MPI_TALLYRESULT, &
                 MPI_STATUS_IGNORE, mpi_err)
#         else
            do k = 1, size(tallies(i) % results, 2)
              do j = 1, size(tallies(i) % results, 1)
                read(UNIT_STATE) tallies(i) % results(j,k) % sum
                read(UNIT_STATE) tallies(i) % results(j,k) % sum_sq
              end do
            end do
#         endif
        end do TALLY_RESULTS
      end if
    end if

    ! Read source bank
    if (run_mode == MODE_EIGENVALUE) then
#     ifdef HDF5
         call read_source()
#     elif MPI
         call read_source(fh)
#     else
         call read_source()
#     endif
    endif

    ! Close statepoint file
    if (.not. source_separate) then
#     ifdef HDF5
        if (master) call h5fclose_f(hdf5_state_point, hdf5_err)
#     elif MPI
        call MPI_FILE_CLOSE(fh, mpi_err)
#     else 
        close(UNIT_STATE)
#     endif
    end if 

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

!===============================================================================
! WRITE_SOURCE writes out the final source distribution to a binary or HDF5
! file that can be used as a starting source in a new simulation
!===============================================================================

#ifdef HDF5
  subroutine write_source()
#elif MPI
  subroutine write_source(fh)
#else
  subroutine write_source()
#endif

#ifdef MPI
    integer :: fh                      ! file handle
    integer(MPI_OFFSET_KIND) :: offset ! offset in memory
#ifdef HDF5
    integer(HID_T) :: hdf5_source
    integer(HSIZE_T) :: dims(1)
    integer(HID_T) :: filespace
    integer(HID_T) :: memspace
    character(6) :: dsetname = 'source'
    integer(HID_T) :: dset_id
    integer(HID_T) :: plist_id
    type(c_ptr) :: f_ptr
#endif
#endif

#ifdef MPI
# ifdef HDF5
    ! Parallel HDF5 must be written out separately
    source_separate = .true.
# endif
#endif

    ! Check if source separate
    if (source_separate) then
#     ifdef MPI
#       ifdef HDF5
          if (master) call h5fclose_f(hdf5_state_point, hdf5_err)
          call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf5_err)
          call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, &
               hdf5_err)
          call h5fcreate_f(trim(path_source), H5F_ACC_TRUNC_F, &
               hdf5_source, hdf5_err, access_prp=plist_id) 
          call h5pclose_f(plist_id, hdf5_err)
          dims(1) = n_particles
          call h5screate_simple_f(1, dims, filespace, hdf5_err)
          call h5dcreate_f(hdf5_source, dsetname, hdf5_bank_t, filespace, &
               dset_id, hdf5_err) 
          call h5sclose_f(filespace, hdf5_err)
          call h5screate_simple_f(1, (/work/), memspace, hdf5_err)
          call h5dget_space_f(dset_id, filespace, hdf5_err)
          call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
               (/bank_first-1/), (/work/), hdf5_err)
          call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf5_err)
          call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdf5_err)
          f_ptr = c_loc(source_bank(1))
          call h5dwrite_f(dset_id, hdf5_bank_t, f_ptr, hdf5_err, &
          file_space_id = filespace, mem_space_id = memspace, &
          xfer_prp = plist_id)
          call h5sclose_f(filespace, hdf5_err)
          call h5sclose_f(memspace, hdf5_err)
          call h5dclose_f(dset_id, hdf5_err)
          call h5pclose_f(plist_id, hdf5_err)
          call h5fclose_f(hdf5_source, hdf5_err)
#       else
          call MPI_FILE_CLOSE(fh, mpi_err)
          call MPI_FILE_OPEN(MPI_COMM_WORLD, path_source, MPI_MODE_CREATE + &
               MPI_MODE_WRONLY, MPI_INFO_NULL, fh, mpi_err)
          if (master) then
            offset = 0
            call MPI_FILE_WRITE_AT(fh, offset, n_particles, 1, MPI_INTEGER8, &
                 MPI_STATUS_IGNORE, mpi_err)
          end if

          ! Set proper offset for source data on this processor
          offset = 8*(1 + rank*maxwork*8)

          ! Write all source sites
          call MPI_FILE_WRITE_AT(fh, offset, source_bank(1), work, MPI_BANK, &
               MPI_STATUS_IGNORE, mpi_err)
          call MPI_FILE_CLOSE(fh, mpi_err)
#       endif
#     else
        ! Open binary source file for writing
        open(UNIT=UNIT_SOURCE, FILE=path_source, STATUS='replace', &
             ACCESS='stream')

        ! Write the number of particles
        write(UNIT=UNIT_SOURCE) n_particles

        ! Write information from the source bank
        write(UNIT=UNIT_SOURCE) source_bank(1:work)

        ! Close binary source file
        close(UNIT=UNIT_SOURCE)
#     endif

    ! append source to statepoint file
    else
#     ifdef HDF5
        dims(1) = work
        call h5screate_simple_f(1, dims, filespace, hdf5_err)
        call h5dcreate_f(hdf5_state_point, "source_bank", hdf5_bank_t, &
             filespace, dset_id, hdf5_err)
        f_ptr = c_loc(source_bank(1))
        call h5dwrite_f(dset_id, hdf5_bank_t, f_ptr, hdf5_err)
        call h5dclose_f(dset_id, hdf5_err)
        call h5sclose_f(filespace, hdf5_err)
#     elif MPI
        if (master) then
          offset = 0
          call MPI_FILE_WRITE_AT(fh, offset, n_particles, 1, MPI_INTEGER8, &
               MPI_STATUS_IGNORE, mpi_err)
        end if

        ! Set proper offset for source data on this processor
        offset = 8*(1 + rank*maxwork*8)

        ! Write all source sites
        call MPI_FILE_WRITE_AT(fh, offset, source_bank(1), work, MPI_BANK, &
             MPI_STATUS_IGNORE, mpi_err)
#     else
        write(UNIT_STATE) source_bank
#     endif

    end if

  end subroutine write_source

!===============================================================================
! READ_SOURCE reads a source distribution
!===============================================================================

#ifdef HDF5
  subroutine read_source()
#elif MPI
  subroutine read_source(fh)
#else
  subroutine read_source()
#endif

    integer    :: i        ! loop over repeating sites
    integer(8) :: n_sites  ! number of sites in binary file
    integer    :: n_repeat ! number of times to repeat a site
#ifdef MPI
    integer    :: fh       ! file handle
    integer(MPI_OFFSET_KIND) :: offset ! offset in memory (0=beginning of file)
    integer    :: n_read   ! number of sites to read on a single process
#endif
#ifdef HDF5
    integer(HID_T) :: hdf5_source
    integer(HSIZE_T) :: dims(1)
    integer(HID_T) :: filespace
    integer(HID_T) :: memspace
    character(6) :: dsetname = 'source'
    integer(HID_T) :: dset_id
    integer(HID_T) :: plist_id
    type(c_ptr) :: f_ptr
#endif

#ifdef MPI
# ifdef HDF5
    ! Parallel HDF5 must be written out separately
    source_separate = .true.
# endif
#endif
!call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)
!call    MPI_ABORT(MPI_COMM_WORLD,-2,mpi_err)
    ! Check if source separate
    if (source_separate) then
#     ifdef MPI
#       ifdef HDF5
          path_source = "source." // trim(to_str(restart_batch)) // ".h5"
          call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf5_err)
          call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, hdf5_err)
          call h5fopen_f(trim(path_source), H5F_ACC_RDONLY_F, hdf5_source, hdf5_err, access_prp=plist_id) 
          call h5pclose_f(plist_id, hdf5_err)
          call h5dopen_f(hdf5_source, dsetname, dset_id, hdf5_err)
          call h5dget_space_f(dset_id, filespace, hdf5_err)
          call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, (/bank_first-1/), (/work/), hdf5_err)
          call h5screate_simple_f(1, (/work/), memspace, hdf5_err)
          call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, hdf5_err)
          call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, hdf5_err)
          f_ptr = c_loc(source_bank(1))
          call h5dread_f(dset_id, hdf5_bank_t, f_ptr, hdf5_err, file_space_id = filespace, mem_space_id = memspace, &
               xfer_prp = plist_id)
          call h5sclose_f(filespace, hdf5_err)
          call h5sclose_f(memspace, hdf5_err)
          call h5dclose_f(dset_id, hdf5_err)
          call h5pclose_f(plist_id, hdf5_err)
          call h5fclose_f(hdf5_source, hdf5_err)
#       else
          ! Close statepoint file
          call MPI_FILE_CLOSE(fh, mpi_err)

          ! Open binary source file for reading
          call MPI_FILE_OPEN(MPI_COMM_WORLD, path_source, MPI_MODE_RDONLY, &
               MPI_INFO_NULL, fh, mpi_err)

          ! Read number of source sites in file
          offset = 0
          call MPI_FILE_READ_AT(fh, offset, n_sites, 1, MPI_INTEGER8, &
               MPI_STATUS_IGNORE, mpi_err)

          ! Set proper offset for source data on this processor
          offset = 8*(1 + rank*maxwork*8)

          ! Read all source sites
          call MPI_FILE_READ_AT(fh, offset, source_bank(1), work, MPI_BANK, &
               MPI_STATUS_IGNORE, mpi_err)

          ! Close binary source file
          call MPI_FILE_CLOSE(fh, mpi_err)
#       endif
#     else
        ! Open binary source file for reading
        open(UNIT=UNIT_SOURCE, FILE=path_source, STATUS='old', &
             ACCESS='stream')

        ! Read number of source sites in file
        read(UNIT=UNIT_SOURCE) n_sites

        ! Read in the source bank
        read(UNIT=UNIT_SOURCE) source_bank(1:n_particles)

        ! Close binary source file
        close(UNIT=UNIT_SOURCE)
#     endif

    ! Read from statepoint file
    else
#     ifdef HDF5
        ! Open dataset for source bank
        call h5dopen_f(hdf5_state_point, "source_bank", dset_id, hdf5_err)

        ! Read source bank
        f_ptr = c_loc(source_bank(1))
        call h5dread_f(dset_id, hdf5_bank_t, f_ptr, hdf5_err)

        ! Close dataset for source bank
        call h5dclose_f(dset_id, hdf5_err)

        ! Close HDF5 state point file
        call h5fclose_f(hdf5_state_point, hdf5_err)
#     elif MPI
        ! Read number of source sites in file
        offset = 0
        call MPI_FILE_READ_AT(fh, offset, n_sites, 1, MPI_INTEGER8, &
             MPI_STATUS_IGNORE, mpi_err)

        ! Set proper offset for source data on this processor
        offset = 8*(1 + rank*maxwork*8)

        ! Read all source sites
        call MPI_FILE_READ_AT(fh, offset, source_bank(1), work, MPI_BANK, &
             MPI_STATUS_IGNORE, mpi_err)
#     else
        read(UNIT=UNIT_STATE) source_bank(1:n_particles)
#     endif

    end if

  end subroutine read_source

end module state_point
