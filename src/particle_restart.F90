module particle_restart 

  use, intrinsic :: ISO_FORTRAN_ENV

  use bank_header,     only: Bank
  use constants
  use geometry_header, only: BASE_UNIVERSE
  use global
  use particle_header, only: deallocate_coord
  use output,          only: write_message, print_particle
  use output_interface
  use physics,         only: transport
  use random_lcg,      only: set_particle_seed
  use source,          only: initialize_particle

  implicit none
  private
  public ::  run_particle_restart

  ! Short names for output and error units
  integer :: ou = OUTPUT_UNIT
  integer :: eu = ERROR_UNIT

contains

!===============================================================================
! READ_PARTICLE_RESTART reads the particle restart file
!===============================================================================

  subroutine read_particle_restart()

    integer :: filetype ! type of restart file
    integer :: revision ! revision of restart file

    ! Write message to user that we are loading a particle restart file
    message = "Loading particle restart file " // trim(path_particle_restart) &
              // "..."
    call write_message(1)

    ! Open up the file for reading
    call file_open(path_particle_restart, 'serial', 'r')

    ! Read data from file
    call read_data(filetype, 'filetype')
    call read_data(revision, 'revision')
    call read_data(current_batch, 'current_batch')
    call read_data(gen_per_batch, 'gen_per_batch')
    call read_data(current_gen, 'current_gen')
    call read_data(n_particles, 'n_particles')
    call read_data(p % id, 'id')
    call read_data(p % wgt, 'weight')
    call read_data(p % E, 'energy')
    call read_data(p % coord % xyz, 'xyz', length = 3)
    call read_data(p % coord % uvw, 'uvw', length = 3)

    ! Set particle last attributes
    p % last_wgt = p % wgt
    p % last_xyz = p % coord % xyz
    p % last_E   = p % E

    ! Close hdf5 file
    call file_close('serial')

  end subroutine read_particle_restart

!===============================================================================
! RUN_PARTICLE_RESTART is the main routine that runs the particle restart
!===============================================================================

  subroutine run_particle_restart()

    integer(8) :: particle_seed

    ! Initialize the particle to be tracked
    allocate(p)
    call initialize_particle()

    ! Read in the restart information
    call read_particle_restart()

    ! Set all tallies to 0 for now (just tracking errors)
    n_tallies = 0

    ! Compute random number seed
    particle_seed = ((current_batch - 1)*gen_per_batch + &
         current_gen - 1)*n_particles + p % id
    call set_particle_seed(particle_seed)

    ! Transport neutron
    call transport()

    ! Write output if particle made it
    call print_particle()

  end subroutine run_particle_restart

end module particle_restart
