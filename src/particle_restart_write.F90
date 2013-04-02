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
#endif HDF5

contains

!===============================================================================
! WRITE_PARTICLE_RESTART
!===============================================================================

  subroutine write_particle_restart()
#ifdef HDF5
    character(MAX_FILE_LEN) :: filename
    integer(HSIZE_T)        :: dims1(1)
    type(Bank), pointer     :: src => null()

    ! set up file name
    filename = 'particle_'//trim(to_str(rank))//'.h5'

    ! create hdf5 file
    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, hdf5_particle_file, hdf5_err)

    ! get information about source particle
    src => source_bank(current_work)

    ! write data to file
    call hdf5_write_integer(hdf5_particle_file, 'current_batch', current_batch)
    call hdf5_write_integer(hdf5_particle_file, 'gen_per_batch', gen_per_batch)
    call hdf5_write_integer(hdf5_particle_file, 'current_gen', current_gen)
    call hdf5_write_long(hdf5_particle_file, 'n_particles', n_particles)
    call hdf5_write_long(hdf5_particle_file, 'id', p % id)
    call hdf5_write_double(hdf5_particle_file, 'weight', src % wgt)
    call hdf5_write_double(hdf5_particle_file, 'energy', src % E)
    dims1 = (/3/)
    call h5ltmake_dataset_double_f(hdf5_particle_file, 'xyz', 1, dims1, &
         src % xyz, hdf5_err)
    call h5ltmake_dataset_double_f(hdf5_particle_file, 'uvw', 1, dims1, &
         src % uvw, hdf5_err)

    ! close hdf5 file
    call h5fclose_f(hdf5_particle_file, hdf5_err)
#endif
  end subroutine write_particle_restart

end module particle_restart_write
