!===============================================================================
! PARTICLE_TRACK handles output of particle tracks to disk.
! 
! TODO: This module writes binary output files via a stream access i.e. it
! writes the particle's location to disk everytime write_particle_track is
! called.  But HDF5 does not allow for a simple implimentation of stream access.
! An attempt was made to write HDF5 files in a stream fashion using an
! extendable dataset, but it always failed at the write call because of an out-
! of-bounds error.  The extendable dataset code is commented out by enclosing it
! within "#if 0"/"#endif" blocks, and has been replaced by code that writes
! particle coordinates to an ever-expanding array which is written to disk in
! its entirety when finalize_particle_track is called.
!
! In the interest of consistency, maybe this module should be changed so that it
! uses the extendable HDF5 dataset, or it writes binary output all in one call
! using the output_interface (like the current HDF5 implimentation).
!===============================================================================

module particle_track

  use constants
  use global
  use string,            only: to_str

#if 0
  use hdf5
#endif
#ifdef HDF5
  use output_interface,  only: file_create, file_close, write_data
#endif

  implicit none

#if 0
  integer, private           :: hdf5_err      ! HDF error code
  integer(HID_T), private    :: track_dset    ! dataset handle
  integer(HID_T), private    :: track_dspace  ! dataspace handle
  integer(HID_T), private    :: track_fh      ! HDF file handle
  integer(HID_T), private    :: cparms        ! chunk parameters
  integer(HSIZE_T), private  :: n_tracks      ! total number of tracks
#endif
#ifdef HDF5
  character(MAX_FILE_LEN), private  :: fname        ! file name
  integer, private                  :: n_tracks     ! total number of tracks
  real(8), private, allocatable     :: coords(:,:)  ! track coordinates
#endif

contains

!===============================================================================
! INITIALIZE_PARTICLE_TRACK opens a particle track output file.
!===============================================================================

  subroutine initialize_particle_track()
    character(MAX_FILE_LEN) :: filename
#if 0
    integer(HSIZE_T)        :: dims(2), max_dims(2), i, j(1)
    integer(HID_T)          :: cparms
    filename = trim(path_output) // 'track_' // trim(to_str(current_batch)) &
         // '_' // trim(to_str(current_gen)) // '_' // trim(to_str(p % id)) &
         // '.h5'
    n_tracks = 0
    ! Create file.
    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, track_fh, hdf5_err)
    ! Create a dataspace with an unlimited max dimension.
    dims = (/3, 1/)
    i = 3
    max_dims = (/i, H5S_UNLIMITED_F/)
    call h5screate_simple_f(2, dims, track_dspace, hdf5_err, max_dims)
    ! Set dataspace chunking.
    call h5pcreate_f(H5P_DATASET_CREATE_F, cparms, hdf5_err)
    call h5pset_chunk_f(cparms, 2, dims, hdf5_err)
    ! Set fill value.
    call h5pset_fill_value_f(cparms, H5T_NATIVE_DOUBLE, 0, hdf5_err)
    ! Create dataset.
    call h5dcreate_f(track_fh, 'coordinates', H5T_NATIVE_DOUBLE, track_dspace, &
         track_dset, hdf5_err, cparms)
#endif
#ifdef HDF5
    filename = trim(path_output) // 'track_' // trim(to_str(current_batch)) &
         // '_' // trim(to_str(current_gen)) // '_' // trim(to_str(p % id)) &
         // '.h5'
    n_tracks = 0
    fname = filename
    allocate(coords(3,1))
#else
    filename = trim(path_output) // 'track_' // trim(to_str(current_batch)) &
         // '_' // trim(to_str(current_gen)) // '_' // trim(to_str(p % id)) &
         // '.binary'
    open(UNIT=UNIT_TRACK, FILE=filename, ACTION="write", &
         STATUS='replace', ACCESS='stream')
#endif
  end subroutine initialize_particle_track

!===============================================================================
! WRITE_PARTICLE_TRACK outputs particle position to a binary file.
!===============================================================================

  subroutine write_particle_track()
#if 0
    integer(HSIZE_T)             :: dims(2), max_dims(2), offset(2), i
    integer(HSIZE_T), parameter  :: cnt(2) = (/3, 1/), write_dims(2) = (/3, 1/)
    real(8)                      :: coords(3)

    ! There is another set of track coordinates.  Incriment the counter.
    n_tracks = n_tracks + 1
    ! Make the dataset the right size to fit the track coordinates.
    i = 3
    dims = (/i, n_tracks/)
    call h5dset_extent_f(track_dset, dims, hdf5_err)
    ! Get the dataspace with the updated dimensions.
    call h5dget_space_f(track_dset, track_dspace, hdf5_err)
    ! Select the hyperslab where the latest coordinates will be written.
    i = 1
    offset = (/i, n_tracks/)
    if (n_tracks < 20) then
    call h5sselect_hyperslab_f(track_dspace, H5S_SELECT_SET_F, offset, cnt, &
         hdf5_err)
    endif
    ! Write the coordinates to the dataset.
    coords = p % coord0 % xyz
    call h5dwrite_f(track_dset, H5T_NATIVE_DOUBLE, coords, &
         write_dims, hdf5_err, file_space_id=track_dspace)
#endif
#ifdef HDF5
    real(8)  :: old_coords(3, n_tracks)

    old_coords = coords
    n_tracks = n_tracks + 1
    deallocate(coords)
    allocate(coords(3, n_tracks))
    coords(:, 1:n_tracks-1) = old_coords
    coords(:, n_tracks) = p % coord0 % xyz
#else
    write(UNIT_TRACK) p % coord0 % xyz
#endif
  end subroutine write_particle_track

!===============================================================================
! FINALIZE_PARTICLE_TRACK closes the particle track file.
!===============================================================================

  subroutine finalize_particle_track()
#if 0
    call h5dclose_f(track_dset, hdf5_err)
    call h5sclose_f(track_dspace, hdf5_err)
    call h5pclose_f(cparms, hdf5_err)
    call h5fclose_f(track_fh, hdf5_err)
#endif
#ifdef HDF5
    call file_create(fname, 'serial')
    call write_data(coords, 'coordinates', length=(/3, n_tracks/))
    call file_close('serial')
    deallocate(coords)
#else
    close(UNIT=UNIT_TRACK)
#endif
  end subroutine finalize_particle_track

end module particle_track
