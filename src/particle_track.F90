!===============================================================================
! PARTICLE_TRACK handles output of particle tracks (the paths taken by particles
! as they are transported through the geometry).
!===============================================================================

module particle_track

  use global
  use output_interface,  only: BinaryOutput
  use particle_header,   only: Particle
  use string,            only: to_str

  implicit none

  integer, private                  :: n_tracks     ! total number of tracks
  real(8), private, allocatable     :: coords(:,:)  ! track coordinates

contains

!===============================================================================
! INITIALIZE_PARTICLE_TRACK
!===============================================================================

  subroutine initialize_particle_track()
    n_tracks = 0
    allocate(coords(1,1))
  end subroutine initialize_particle_track

!===============================================================================
! WRITE_PARTICLE_TRACK copies particle position to an array.
!===============================================================================

  subroutine write_particle_track(p)
    type(Particle), intent(in)  :: p

    real(8)  :: old_coords(3, n_tracks)

    ! Save the coordinates gathered in previous calls.
    old_coords = coords

    ! Add another column to coords.
    n_tracks = n_tracks + 1
    deallocate(coords)
    allocate(coords(3, n_tracks))

    ! Put the old coordinates back into the array.
    coords(:, 1:n_tracks-1) = old_coords

    ! Write current coordinates into the newest column.
    coords(:, n_tracks) = p % coord0 % xyz
  end subroutine write_particle_track

!===============================================================================
! FINALIZE_PARTICLE_TRACK writes the particle track array to disk.
!
! output_interface currently does not support writing 2D binary arrays so there
! are two different versions of this subroutine; an HDF version which simply
! writes the array and a binary version which flattens the array into a 1D shape
! before writing.
!===============================================================================

#ifdef HDF5

  subroutine finalize_particle_track(p)
    type(Particle), intent(in)  :: p

    character(MAX_FILE_LEN)  :: fname
    type(BinaryOutput)       :: binout

    fname = trim(path_output) // 'track_' // trim(to_str(current_batch)) &
         // '_' // trim(to_str(current_gen)) // '_' // trim(to_str(p % id)) &
         // '.h5'
    call file_create(fname)
    call write_data(coords, 'coordinates', length=(/3, n_tracks/))
    call file_close()
    deallocate(coords)
  end subroutine finalize_particle_track

#else

  subroutine finalize_particle_track(p)
    type(Particle), intent(in)  :: p

    character(MAX_FILE_LEN)  :: fname
    type(BinaryOutput)       :: binout
    integer                  :: i
    real(8)                  :: flat_coords(3*n_tracks)

    fname = trim(path_output) // 'track_' // trim(to_str(current_batch)) &
         // '_' // trim(to_str(current_gen)) // '_' // trim(to_str(p % id)) &
         // '.binary'
    call binout % file_create(fname)
    do i=1, n_tracks
      flat_coords(3*i-2 : 3*i) = coords(:,i)
    end do
    call binout % write_data(flat_coords, 'coordinates', length=3*n_tracks)
    call binout % file_close()
    deallocate(coords)
  end subroutine finalize_particle_track

#endif

end module particle_track
