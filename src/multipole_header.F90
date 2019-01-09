module multipole_header

  use constants
  use hdf5_interface
  use error,            only: fatal_error
  use string,           only: to_str

  implicit none

contains

!===============================================================================
! CHECK_WMP_VERSION checks for the right version of WMP data within HDF5
! files
!===============================================================================

  subroutine check_wmp_version(file_id)
    integer(HID_T), intent(in) :: file_id

    integer, allocatable :: version(:)

    if (attribute_exists(file_id, 'version')) then
      call read_attribute(version, file_id, 'version')
      if (version(1) /= WMP_VERSION(1)) then
        call fatal_error("WMP data format uses version " // trim(to_str(&
             version(1))) // "." // trim(to_str(version(2))) // " whereas &
             &your installation of OpenMC expects version " // trim(to_str(&
             WMP_VERSION(1))) // ".x data.")
      end if
    else
      call fatal_error("WMP data does not indicate a version. Your &
           &installation of OpenMC expects version " // trim(to_str(&
           WMP_VERSION(1))) // ".x data.")
    end if
  end subroutine check_wmp_version

end module multipole_header
