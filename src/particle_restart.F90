module particle_restart 

  use constants
  use global
  use string

#ifdef HDF5
  use hdf5_interface
#endif

  private
  public ::  write_particle_restart

  integer(HID_T) :: hdf5_particle_file

contains

!===============================================================================
! WRITE_PARTICLE_RESTART
!===============================================================================

  subroutine write_particle_restart()

    character(MAX_FILE_LEN) :: filename

    ! set up file name
    filename = 'particle_'//to_str(p % id)//'.h5'

    ! create hdf5 file
    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, hdf5_particle_file, hdf5_err)

    ! write data to file


    ! close hdf5 file
    call h5fclose_f(hdf5_particle_file, hdf5_err)

  end subroutine write_particle_restart

!===============================================================================
! READ_PARTICLE_RESTART
!===============================================================================

  subroutine read_particle_restart

  end subroutine read_particle_restart

!===============================================================================
! RUN_PARTICLE_RESTART
!===============================================================================

  subroutine run_particle_restart

  end subroutine run_particle_restart

end module particle_restart
