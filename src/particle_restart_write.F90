module particle_restart_write

  use, intrinsic :: ISO_FORTRAN_ENV

  use bank_header,  only: Bank
  use global
  use string,       only: to_str

#ifdef HDF5
  use hdf5_interface
#endif

  implicit none
  private
  public ::  write_particle_restart

#ifdef HDF5
  integer(HID_T) :: hdf5_particle_file
#endif

contains

!===============================================================================
! WRITE_PARTICLE_RESTART
!===============================================================================

  subroutine write_particle_restart()

    ! Dont write another restart file if in particle restart mode
    if (run_mode == MODE_PARTICLE) return

    ! write out binary or HDF5 file
#ifdef HDF5
    call write_particle_restart_hdf5()
#else
    call write_particle_restart_binary()
#endif

  end subroutine write_particle_restart

#ifdef HDF5

!===============================================================================
! WRITE_PARTICLE_RESTART_HDF5
!===============================================================================

  subroutine write_particle_restart_hdf5()

    character(MAX_FILE_LEN) :: filename
    type(Bank), pointer     :: src => null()

    ! set up file name
    filename = trim(path_output) // 'particle_' // trim(to_str(current_batch)) &
         // '_' // trim(to_str(current_work)) // '.h5'

    ! create hdf5 file
    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, hdf5_particle_file, hdf5_err)

    ! get information about source particle
    src => source_bank(current_work)

    ! write data to file
    call hdf5_write_integer(hdf5_particle_file, 'current_batch', current_batch)
    call hdf5_write_integer(hdf5_particle_file, 'gen_per_batch', gen_per_batch)
    call hdf5_write_integer(hdf5_particle_file, 'current_gen', current_gen)
    call hdf5_write_long(hdf5_particle_file, 'n_particles', n_particles, hdf5_integer8_t)
    call hdf5_write_long(hdf5_particle_file, 'id', p % id, hdf5_integer8_t)
    call hdf5_write_double(hdf5_particle_file, 'weight', src % wgt)
    call hdf5_write_double(hdf5_particle_file, 'energy', src % E)
    dims1 = (/3/)
    call h5ltmake_dataset_double_f(hdf5_particle_file, 'xyz', 1, dims1, &
         src % xyz, hdf5_err)
    call h5ltmake_dataset_double_f(hdf5_particle_file, 'uvw', 1, dims1, &
         src % uvw, hdf5_err)

    ! close hdf5 file
    call h5fclose_f(hdf5_particle_file, hdf5_err)

  end subroutine write_particle_restart_hdf5

#endif

!===============================================================================
! WRITE_PARTICLE_RESTART_BINARY
!===============================================================================

  subroutine write_particle_restart_binary()

    character(MAX_FILE_LEN) :: filename
    type(Bank), pointer     :: src => null()

    ! set up file name
    filename = trim(path_output) // 'particle_' // trim(to_str(current_batch)) &
         // '_' // trim(to_str(current_work)) // '.binary'

    ! create hdf5 file
    open(UNIT=UNIT_PARTICLE, FILE=filename, STATUS='replace', &
         ACCESS='stream')

    ! get information about source particle
    src => source_bank(current_work)

    ! write data to file
    write(UNIT_PARTICLE) current_batch
    write(UNIT_PARTICLE) gen_per_batch
    write(UNIT_PARTICLE) current_gen
    write(UNIT_PARTICLE) n_particles
    write(UNIT_PARTICLE) p % id
    write(UNIT_PARTICLE) src % wgt
    write(UNIT_PARTICLE) src % E
    write(UNIT_PARTICLE) src % xyz
    write(UNIT_PARTICLE) src % uvw

    ! close hdf5 file
    close(UNIT_PARTICLE)

  end subroutine write_particle_restart_binary

end module particle_restart_write
