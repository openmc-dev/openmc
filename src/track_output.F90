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

  type, private :: TrackCoordinates
    real(8), allocatable :: coords(:,:)
  end type TrackCoordinates

  type(TrackCoordinates), private, allocatable :: tracks(:)
!$omp threadprivate(tracks)

contains

!===============================================================================
! INITIALIZE_PARTICLE_TRACK allocates the array to store particle track
! information
!===============================================================================

  subroutine initialize_particle_track()
    allocate(tracks(1))
  end subroutine initialize_particle_track

!===============================================================================
! WRITE_PARTICLE_TRACK copies particle position to an array.
!===============================================================================

  subroutine write_particle_track(p)
    type(Particle), intent(in)  :: p
    real(8), allocatable :: new_coords(:, :)

    integer :: i
    integer :: n_tracks

    ! Add another column to coords
    i = size(tracks)
    if (allocated(tracks(i) % coords)) then
      n_tracks = size(tracks(i) % coords, 2)
      allocate(new_coords(3, n_tracks + 1))
      new_coords(:, 1:n_tracks) = tracks(i) % coords
      call move_alloc(FROM=new_coords, TO=tracks(i) % coords)
    else
      n_tracks = 0
      allocate(tracks(i) % coords(3, 1))
    end if

    ! Write current coordinates into the newest column.
    n_tracks = n_tracks + 1
    tracks(i) % coords(:, n_tracks) = p % coord(1) % xyz
  end subroutine write_particle_track

!===============================================================================
! ADD_PARTICLE_TRACK creates a new entry in the track coordinates for a
! secondary particle
!===============================================================================

  subroutine add_particle_track()
    type(TrackCoordinates), allocatable :: new_tracks(:)

    integer :: i

    ! Determine current number of particle tracks
    i = size(tracks)

    ! Create array one larger than current
    allocate(new_tracks(i + 1))

    ! Copy memory and move allocation
    new_tracks(1:i) = tracks(i)
    call move_alloc(FROM=new_tracks, TO=tracks)

  end subroutine add_particle_track

!===============================================================================
! FINALIZE_PARTICLE_TRACK writes the particle track array to disk.
!===============================================================================

  subroutine finalize_particle_track(p)
    type(Particle), intent(in)  :: p

    integer                  :: length(2)
    character(MAX_FILE_LEN)  :: fname
    type(BinaryOutput)       :: binout

    integer :: i
    integer, allocatable :: n_coords(:)
    integer :: n_particle_tracks

#ifdef HDF5
    fname = trim(path_output) // 'track_' // trim(to_str(current_batch)) &
         // '_' // trim(to_str(current_gen)) // '_' // trim(to_str(p % id)) &
         // '.h5'
#else
    fname = trim(path_output) // 'track_' // trim(to_str(current_batch)) &
         // '_' // trim(to_str(current_gen)) // '_' // trim(to_str(p % id)) &
         // '.binary'
#endif

    ! Determine total number of particles and number of coordinates for each
    n_particle_tracks = size(tracks)
    allocate(n_coords(n_particle_tracks))
    do i = 1, n_particle_tracks
      n_coords(i) = size(tracks(i) % coords, 2)
    end do

!$omp critical (FinalizeParticleTrack)
    call binout % file_create(fname)
    call binout % write_data(n_particle_tracks, 'n_particles')
    call binout % write_data(n_coords, 'n_coords', length=n_particle_tracks)
    do i = 1, n_particle_tracks
      length(:) = [3, n_coords(i)]
      call binout % write_data(tracks(i) % coords, 'coordinates_' // &
           trim(to_str(i)), length=length)
    end do
    call binout % file_close()
!$omp end critical (FinalizeParticleTrack)
    deallocate(tracks)
  end subroutine finalize_particle_track

end module track_output
