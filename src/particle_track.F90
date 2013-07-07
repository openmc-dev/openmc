module particle_track

  use constants
  use global
  !use output_interface  only: file_create, file_close, write_double1Darray
  use string,           only: to_str

#ifdef HDF5
  use hdf5
#endif

  implicit none

#ifdef HDF5
  integer, private         :: hdf5_err        ! HDF5 error code
  integer(HID_T), private  :: track_dset      ! dataset handle
  integer(HID_T), private  :: track_dspace    ! master dataspace handle
  integer(HID_T), private  :: track_selspace  ! selected dataspace handle
  integer(HID_T), private  :: track_fh        ! HDF5 file handle
  integer(HID_T), private  :: cparms          ! chunk parameters
  integer(HSIZE_T), private  :: n_tracks        ! total number of tracks
#endif

contains

!===============================================================================
! INITIALIZE_PARTICLE_TRACK opens a particle track output file.
! 
! TODO: This subroutine needs to be modified to work with HDF5 files.  It
! should also probably write a header that identifies the file as a particle
! track and maybe adds particle identifying information (batch #, etc.).
!===============================================================================

  subroutine initialize_particle_track()
    character(MAX_FILE_LEN) :: filename
    ! Replace the following with '#ifdef HDF5'.
#if 0
    integer(HSIZE_T)        :: dims(2), max_dims(2), i
    integer(HID_T)          :: cparms
    integer(HID_T)          :: dapl
#endif

    filename = trim(path_output) // 'track_' // trim(to_str(current_batch)) &
         // '_' // trim(to_str(current_gen)) // '_' // trim(to_str(p % id))

    ! Replace the following with '#ifdef HDF5'.
#if 0
    filename = filename // '.h5'
    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, track_fh, hdf5_err)
    dims = (/3, 1/)
    i = 3
    max_dims = (/i, H5S_UNLIMITED_F/)
    call h5screate_simple_f(2, dims, track_dspace, hdf5_err, max_dims)
    call h5pcreate_f(H5P_DATASET_CREATE_F, cparms, hdf5_err)
    call h5pset_chunk_f(cparms, 2, dims, hdf5_err)
    call h5pset_fill_value_f(cparms, H5T_NATIVE_DOUBLE, 0, hdf5_err)
    call h5dcreate_f(track_fh, 'coordinates', H5T_NATIVE_DOUBLE, track_dspace, &
         track_dset, hdf5_err, cparms)
    n_tracks = 1
#else
    filename = filename // '.binary'
    open(UNIT=UNIT_TRACK, FILE=filename, ACTION="write", &
         STATUS='replace', ACCESS='stream')
#endif

  end subroutine initialize_particle_track

!===============================================================================
! WRITE_PARTICLE_TRACK outputs particle position to a binary file.
! 
! TODO: This subroutine needs to be modified to work with HDF5 files.  Perhaps
! it should also be made somehow more general so that it can output
! information other than just particle position.
!===============================================================================

  subroutine write_particle_track()
    ! Replace the following with '#ifdef HDF5'.
#if 0
    integer(HSIZE_T)             :: dims(2), max_dims(2), offset(2), i
    integer(HSIZE_T), parameter  :: cnt(2) = (/3, 1/), write_dims(2) = (/3, 1/)
    real(8)                      :: coords(3, 1)
    integer(HSIZE_T)                      :: test

    i = 3
    dims = (/i, n_tracks/)
    call h5dget_space_f(track_dset, track_dspace, hdf5_err)
    call h5dset_extent_f(track_dset, dims, hdf5_err)
    call h5dget_space_f(track_dset, track_dspace, hdf5_err)
    i = 1
    offset = (/i, n_tracks/)
    call h5sselect_hyperslab_f(track_dspace, H5S_SELECT_SET_F, offset, cnt, &
         hdf5_err)
    !call h5sget_select_npoints_f(track_dspace, test, hdf5_err)
    print *, test
    coords(:, 1) = p % coord0 % xyz
    print *, coords
    if (n_tracks<20) then
    call h5dwrite_f(track_dset, H5T_NATIVE_DOUBLE, coords, &
         write_dims, hdf5_err, file_space_id=track_dspace)
    end if
    !call h5sselect_hyperslab_f(track_dspace, H5S_SELECT_SET_F, offset, dims, &
    !     hdf5_err)
    !i = 3
    n_tracks = n_tracks + 1
    !dims = (/i, n_tracks/)
    !call h5sset_extent_simple_f(track_dspace, 2, dims, max_dims, hdf5_err)
    !call h5dset_extent_f(track_dset, dims, hdf5_err)
#else
    write(UNIT_TRACK) p % coord0 % xyz
#endif
  end subroutine write_particle_track

!===============================================================================
! FINALIZE_PARTICLE_TRACK closes the particle track file.
! 
! TODO: This subroutine needs to be modified to work with HDF5 files.
!===============================================================================

  subroutine finalize_particle_track()
    ! Replace the following with '#ifdef HDF5'.
#if 0
    integer(HSIZE_T)             :: dims(2), max_dims(2), i

    call h5sget_simple_extent_dims_f(track_dspace, dims, max_dims, hdf5_err)
    !i = 3
    !dims = (/i, dims(2)-1/)
    !call h5dextend_f(track_dset, dims, hdf5_err)

    call h5dclose_f(track_dset, hdf5_err)
    call h5sclose_f(track_dspace, hdf5_err)
    call h5pclose_f(cparms, hdf5_err)
    call h5fclose_f(track_fh, hdf5_err)
#else
    close(UNIT=UNIT_TRACK)
#endif
  end subroutine finalize_particle_track

end module particle_track
