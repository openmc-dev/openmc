module particle_restart

  use, intrinsic :: ISO_FORTRAN_ENV

  use bank_header,      only: Bank
  use constants
  use geometry_header,  only: BASE_UNIVERSE
  use global
  use hdf5_interface,   only: file_open, file_close, read_dataset
  use output,           only: write_message, print_particle
  use particle_header,  only: Particle
  use random_lcg,       only: set_particle_seed
  use tracking,         only: transport

  use hdf5, only: HID_T

  implicit none
  private
  public ::  run_particle_restart

contains

!===============================================================================
! RUN_PARTICLE_RESTART is the main routine that runs the particle restart
!===============================================================================

  subroutine run_particle_restart()

    integer(8) :: particle_seed
    integer :: previous_run_mode
    type(Particle) :: p

    ! Set verbosity high
    verbosity = 10

    ! Initialize the particle to be tracked
    call p % initialize()

    ! Read in the restart information
    call read_particle_restart(p, previous_run_mode)

    ! Set all tallies to 0 for now (just tracking errors)
    n_tallies = 0

    ! Compute random number seed
    select case (previous_run_mode)
    case (MODE_EIGENVALUE)
      particle_seed = ((current_batch - 1)*gen_per_batch + &
           current_gen - 1)*n_particles + p % id
    case (MODE_FIXEDSOURCE)
      particle_seed = p % id
    end select

    call set_particle_seed(particle_seed)

    ! Transport neutron
    call transport(p)

    ! Write output if particle made it
    call print_particle(p)

  end subroutine run_particle_restart

!===============================================================================
! READ_PARTICLE_RESTART reads the particle restart file
!===============================================================================

  subroutine read_particle_restart(p, previous_run_mode)
    type(Particle), intent(inout) :: p
    integer, intent(inout) :: previous_run_mode

    integer :: int_scalar
    integer(HID_T) :: file_id
    character(MAX_WORD_LEN) :: tempstr

    ! Write meessage
    call write_message("Loading particle restart file " &
         &// trim(path_particle_restart) // "...", 1)

    ! Open file
    file_id = file_open(path_particle_restart, 'r')

    ! Read data from file
    call read_dataset(tempstr, file_id, 'filetype')
    call read_dataset(int_scalar, file_id, 'revision')
    call read_dataset(current_batch, file_id, 'current_batch')
    call read_dataset(gen_per_batch, file_id, 'gen_per_batch')
    call read_dataset(current_gen, file_id, 'current_gen')
    call read_dataset(n_particles, file_id, 'n_particles')
    call read_dataset(tempstr, file_id, 'run_mode')
    select case (tempstr)
    case ('k-eigenvalue')
      previous_run_mode = MODE_EIGENVALUE
    case ('fixed source')
      previous_run_mode = MODE_FIXEDSOURCE
    end select
    call read_dataset(p % id, file_id, 'id')
    call read_dataset(p % wgt, file_id, 'weight')
    call read_dataset(p % E, file_id, 'energy')
    call read_dataset(p % g, file_id, 'energy_group')
    call read_dataset(p % coord(1) % xyz, file_id, 'xyz')
    call read_dataset(p % coord(1) % uvw, file_id, 'uvw')

    ! Set particle last attributes
    p % last_wgt = p % wgt
    p % last_xyz = p % coord(1)%xyz
    p % last_uvw = p % coord(1)%uvw
    p % last_E   = p % E
    p % last_g   = p % g

    ! Close hdf5 file
    call file_close(file_id)

  end subroutine read_particle_restart

end module particle_restart
