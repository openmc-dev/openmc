module hdf5_interface

  use constants
  use global

#ifdef HDF5
  use hdf5, only: HSIZE_T, h5open_f, h5fcreate_f, h5fclose_f, h5close_f
  use h5lt, only: h5ltmake_dataset_int_f
#endif

  implicit none

contains

#ifdef HDF5

!===============================================================================
! OPEN_HDF5_OUTPUT
!===============================================================================

  subroutine open_hdf5_output()

    character(9), parameter :: filename = "output.h5" ! File name
    integer :: error  ! Error flag

    ! Initialize FORTRAN interface.
    call h5open_f (error)

    ! Create a new file using default properties.
    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, hdf5_output_file, error)

  end subroutine open_hdf5_output

!===============================================================================
! WRITE_HDF5_SUMMARY
!===============================================================================

  subroutine write_hdf5_summary()

    integer          :: error ! Error flag
    integer          :: rank = 1
    integer(HSIZE_T) :: dims(1) = (/1/)
     
    ! Write version information
    call h5ltmake_dataset_int_f(hdf5_output_file, "/version_major", rank, &
         dims, (/ VERSION_MAJOR /), error)
    call h5ltmake_dataset_int_f(hdf5_output_file, "/version_minor", rank, &
         dims, (/ VERSION_MINOR /), error)
    call h5ltmake_dataset_int_f(hdf5_output_file, "/version_release", rank, & 
         dims, (/ VERSION_RELEASE /), error)

  end subroutine write_hdf5_summary

!===============================================================================
! CLOSE_HDF5_OUTPUT
!===============================================================================

  subroutine close_hdf5_output()

    integer :: error  ! Error flag

    ! Terminate access to the file.
    call h5fclose_f(hdf5_output_file, error)

    ! Close FORTRAN interface.
    call h5close_f(error)

  end subroutine close_hdf5_output

#endif

end module hdf5_interface
