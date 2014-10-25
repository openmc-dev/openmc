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

  type(BinaryOutput)        :: pr      ! Binary file

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

    integer :: int_scalar
    integer, intent(inout) :: previous_run_mode
    type(Particle), intent(inout) :: p

    ! Write meessage
    call write_message("Loading particle restart file " &
         &// trim(path_particle_restart) // "...", 1)

    ! Open file
    call pr % file_open(path_particle_restart, 'r')

    ! Read data from file
    call pr % read_data(int_scalar, 'filetype')
    call pr % read_data(int_scalar, 'revision')
    call pr % read_data(current_batch, 'current_batch')
    call pr % read_data(gen_per_batch, 'gen_per_batch')
    call pr % read_data(current_gen, 'current_gen')
    call pr % read_data(n_particles, 'n_particles')
    call pr % read_data(previous_run_mode, 'run_mode')
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
