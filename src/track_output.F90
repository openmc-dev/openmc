!===============================================================================
! TRACK_OUTPUT handles output of particle tracks (the paths taken by particles
! as they are transported through the geometry).
!===============================================================================

module track_output

  use global
  use output_interface,  only: BinaryOutput
  use particle_header,   only: Particle
  use string,            only: to_str

  implicit none

  integer, private                  :: n_tracks     ! total number of tracks
  real(8), private, allocatable     :: coords(:,:)  ! track coordinates
!$omp threadprivate(n_tracks, coords)

contains

!===============================================================================
! INITIALIZE_PARTICLE_TRACK
!===============================================================================

  subroutine initialize_particle_track()
    n_tracks = 0
  end subroutine initialize_particle_track

!===============================================================================
! WRITE_PARTICLE_TRACK copies particle position to an array.
!===============================================================================

  subroutine write_particle_track(p)
    type(Particle), intent(in)  :: p

    real(8), allocatable :: new_coords(:, :)

    ! Add another column to coords
    n_tracks = n_tracks + 1
    if (allocated(coords)) then
      allocate(new_coords(3, n_tracks))
      new_coords(:, 1:n_tracks-1) = coords
      call move_alloc(FROM=new_coords, TO=coords)
    else
      allocate(coords(3,1))
    end if

    ! Write current coordinates into the newest column.
    coords(:, n_tracks) = p % coord0 % xyz
  end subroutine write_particle_track

!===============================================================================
! FINALIZE_PARTICLE_TRACK writes the particle track array to disk.
!===============================================================================

  subroutine finalize_particle_track(p)
    type(Particle), intent(in)  :: p

    integer                  :: length(2)
    character(MAX_FILE_LEN)  :: fname
    type(BinaryOutput)       :: binout

#ifdef HDF5
    fname = trim(path_output) // 'track_' // trim(to_str(current_batch)) &
         // '_' // trim(to_str(current_gen)) // '_' // trim(to_str(p % id)) &
         // '.h5'
#else
    fname = trim(path_output) // 'track_' // trim(to_str(current_batch)) &
         // '_' // trim(to_str(current_gen)) // '_' // trim(to_str(p % id)) &
         // '.binary'
#endif
    call binout % file_create(fname)
    length = [3, n_tracks]
    call binout % write_data(coords, 'coordinates', length=length)
    call binout % file_close()
    deallocate(coords)
  end subroutine finalize_particle_track

end module track_output
