module particle_restart

  use, intrinsic :: ISO_FORTRAN_ENV

  use bank_header,      only: Bank
  use constants
  use geometry_header,  only: BASE_UNIVERSE
  use global
  use output,           only: write_message, print_particle
  use output_interface, only: BinaryOutput
  use particle_header,  only: Particle
  use random_lcg,       only: set_particle_seed
  use tracking,         only: transport

  implicit none
  private
  public ::  run_particle_restart

  ! Binary file
  type(BinaryOutput) :: pr

contains

!===============================================================================
! RUN_PARTICLE_RESTART is the main routine that runs the particle restart
!===============================================================================

  subroutine run_particle_restart()

    integer(8) :: particle_seed
    type(Particle) :: p

    ! Set verbosity high
    verbosity = 10

    ! Initialize the particle to be tracked
    call p % initialize()

    ! Read in the restart information
    call read_particle_restart(p)

    ! Set all tallies to 0 for now (just tracking errors)
    n_tallies = 0

    ! Compute random number seed
    particle_seed = ((current_batch - 1)*gen_per_batch + &
         current_gen - 1)*n_particles + p % id
    call set_particle_seed(particle_seed)

    ! Transport neutron
    call transport(p)

    ! Write output if particle made it
    call print_particle(p)

  end subroutine run_particle_restart

!===============================================================================
! READ_PARTICLE_RESTART reads the particle restart file
!===============================================================================

  subroutine read_particle_restart(p)

    integer :: int_scalar
    type(Particle), intent(inout) :: p

    ! Write meessage
    message = "Loading particle restart file " // trim(path_particle_restart) &
              // "..."
    call write_message(1)

    ! Open file
    call pr % file_open(path_particle_restart, 'r')

    ! Read data from file
    call pr % read_data(int_scalar, 'filetype')
    call pr % read_data(int_scalar, 'revision')
    call pr % read_data(current_batch, 'current_batch')
    call pr % read_data(gen_per_batch, 'gen_per_batch')
    call pr % read_data(current_gen, 'current_gen')
    call pr % read_data(n_particles, 'n_particles')
    call pr % read_data(p % id, 'id')
    call pr % read_data(p % wgt, 'weight')
    call pr % read_data(p % E, 'energy')
    call pr % read_data(p % coord % xyz, 'xyz', length=3)
    call pr % read_data(p % coord % uvw, 'uvw', length=3)

    ! Set particle last attributes
    p % last_wgt = p % wgt
    p % last_xyz = p % coord % xyz
    p % last_uvw = p % coord % uvw
    p % last_E   = p % E

    ! Close hdf5 file
    call pr % file_close()

  end subroutine read_particle_restart

end module particle_restart
