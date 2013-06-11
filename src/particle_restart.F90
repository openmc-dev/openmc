module particle_restart 

  use, intrinsic :: ISO_FORTRAN_ENV

  use bank_header,     only: Bank
  use constants
  use geometry_header, only: BASE_UNIVERSE
  use global
  use particle_header, only: deallocate_coord
  use output,          only: write_message, print_particle
  use physics,         only: transport
  use random_lcg,      only: set_particle_seed
  use source,          only: initialize_particle

#ifdef HDF5
  use hdf5_interface 
#endif

  implicit none
  private
  public ::  run_particle_restart

#ifdef HDF5
  integer(HID_T) :: hdf5_particle_file
#endif

  ! Short names for output and error units
  integer :: ou = OUTPUT_UNIT
  integer :: eu = ERROR_UNIT

contains

#ifdef HDF5

!===============================================================================
! READ_HDF5_PARTICLE_RESTART
!===============================================================================

  subroutine read_hdf5_particle_restart()

    ! write meessage
    message = "Loading particle restart file " // trim(path_particle_restart) &
              // "..."
    call write_message(1)

    ! open hdf5 file
    call h5fopen_f(path_particle_restart, H5F_ACC_RDONLY_F, hdf5_particle_file,&
                   hdf5_err)

    ! read data from file
    call hdf5_read_integer(hdf5_particle_file, 'current_batch', current_batch)
    call hdf5_read_integer(hdf5_particle_file, 'gen_per_batch', gen_per_batch)
    call hdf5_read_integer(hdf5_particle_file, 'current_gen', current_gen)
    call hdf5_read_long(hdf5_particle_file, 'n_particles', n_particles, hdf5_integer8_t)
    call hdf5_read_long(hdf5_particle_file, 'id', p % id, hdf5_integer8_t)
    call hdf5_read_double(hdf5_particle_file, 'weight', p % wgt)
    call hdf5_read_double(hdf5_particle_file, 'energy', p % E)
    dims1 = (/3/)
    call h5ltread_dataset_double_f(hdf5_particle_file, 'xyz', p % coord % xyz, &
         dims1, hdf5_err)
    call h5ltread_dataset_double_f(hdf5_particle_file, 'uvw', p % coord % uvw, &
         dims1, hdf5_err)

    ! set particle last attributes
    p % last_wgt = p % wgt
    p % last_xyz = p % coord % xyz
    p % last_E   = p % E

    ! close hdf5 file
    call h5fclose_f(hdf5_particle_file, hdf5_err)

  end subroutine read_hdf5_particle_restart

#endif

!===============================================================================
! READ_BINARY_PARTICLE_RESTART
!===============================================================================

  subroutine read_binary_particle_restart()

    integer :: filetype
    integer :: revision

    ! write meessage
    message = "Loading particle restart file " // trim(path_particle_restart) &
              // "..."
    call write_message(1)

    ! open file
    open(UNIT=UNIT_PARTICLE, FILE=path_particle_restart, STATUS='old', &
         ACCESS='stream')

    ! read data from file
    read(UNIT_PARTICLE) filetype
    read(UNIT_PARTICLE) revision
    read(UNIT_PARTICLE) current_batch
    read(UNIT_PARTICLE) gen_per_batch
    read(UNIT_PARTICLE) current_gen
    read(UNIT_PARTICLE) n_particles
    read(UNIT_PARTICLE) p % id
    read(UNIT_PARTICLE) p % wgt
    read(UNIT_PARTICLE) p % E
    read(UNIT_PARTICLE) p % coord % xyz
    read(UNIT_PARTICLE) p % coord % uvw

    ! set particle last attributes
    p % last_wgt = p % wgt
    p % last_xyz = p % coord % xyz
    p % last_E   = p % E

    ! close hdf5 file
    close(UNIT_PARTICLE)

  end subroutine read_binary_particle_restart

!===============================================================================
! RUN_PARTICLE_RESTART
!===============================================================================

  subroutine run_particle_restart()

    integer(8) :: particle_seed

    ! initialize the particle to be tracked
    allocate(p)
    call initialize_particle()

    ! read in the restart information
#ifdef HDF5
    call read_hdf5_particle_restart()
#else
    call read_binary_particle_restart()
#endif

    ! set all tallies to 0 for now (just tracking errors)
    n_tallies = 0

    ! compute random number seed
    particle_seed = ((current_batch - 1)*gen_per_batch + &
         current_gen - 1)*n_particles + p % id
    call set_particle_seed(particle_seed)

    ! transport neutron
    call transport()

    ! write output if particle made it
    call print_particle()

  end subroutine run_particle_restart

end module particle_restart
