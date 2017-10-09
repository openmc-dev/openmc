!===============================================================================
! TRACK_OUTPUT handles output of particle tracks (the paths taken by particles
! as they are transported through the geometry).
!===============================================================================

module track_output

  use hdf5

  use constants
  use hdf5_interface
  use particle_header, only: Particle
  use settings,        only: path_output
  use simulation_header
  use string,          only: to_str

  implicit none
  private

  type TrackCoordinates
    real(8), allocatable :: coords(:,:)
  end type TrackCoordinates

  type(TrackCoordinates), allocatable :: tracks(:)
!$omp threadprivate(tracks)

  public :: initialize_particle_track
  public :: write_particle_track
  public :: add_particle_track
  public :: finalize_particle_track

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
    if (allocated(tracks(i)%coords)) then
      n_tracks = size(tracks(i)%coords, 2)
      allocate(new_coords(3, n_tracks + 1))
      new_coords(:, 1:n_tracks) = tracks(i)%coords
      call move_alloc(FROM=new_coords, TO=tracks(i)%coords)
    else
      n_tracks = 0
      allocate(tracks(i)%coords(3, 1))
    end if

    ! Write current coordinates into the newest column.
    n_tracks = n_tracks + 1
    tracks(i)%coords(:, n_tracks) = p%coord(1)%xyz
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

    integer :: i
    integer :: n_particle_tracks
    integer(HID_T) :: file_id
    character(MAX_FILE_LEN) :: fname
    integer, allocatable :: n_coords(:)

    fname = trim(path_output) // 'track_' // trim(to_str(current_batch)) &
         // '_' // trim(to_str(current_gen)) // '_' // trim(to_str(p%id)) &
         // '.h5'

    ! Determine total number of particles and number of coordinates for each
    n_particle_tracks = size(tracks)
    allocate(n_coords(n_particle_tracks))
    do i = 1, n_particle_tracks
      n_coords(i) = size(tracks(i)%coords, 2)
    end do

!$omp critical (FinalizeParticleTrack)
    file_id = file_create(fname)
    call write_attribute(file_id, 'filetype', 'track')
    call write_attribute(file_id, 'version', VERSION_TRACK)
    call write_attribute(file_id, 'n_particles', n_particle_tracks)
    call write_attribute(file_id, 'n_coords', n_coords)
    do i = 1, n_particle_tracks
      call write_dataset(file_id, 'coordinates_' // trim(to_str(i)), &
           tracks(i)%coords)
    end do
    call file_close(file_id)
!$omp end critical (FinalizeParticleTrack)
    deallocate(tracks)
  end subroutine finalize_particle_track

end module track_output
