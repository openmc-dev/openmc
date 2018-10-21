module multipole_header

  use constants
  use hdf5_interface
  use error,            only: fatal_error
  use string,           only: to_str

  implicit none

  !========================================================================
  ! Multipole related constants

  ! Constants that determine which value to access
  integer, parameter :: MP_EA = 1       ! Pole

  ! Residue indices
  integer, parameter :: MP_RS = 2, &    ! Residue scattering
                        MP_RA = 3, &    ! Residue absorption
                        MP_RF = 4       ! Residue fission

  ! Polynomial fit indices
  integer, parameter :: FIT_S = 1, &    ! Scattering
                        FIT_A = 2, &    ! Absorption
                        FIT_F = 3       ! Fission

!===============================================================================
! MULTIPOLE contains all the components needed for the windowed multipole
! temperature dependent cross section libraries for the resolved resonance
! region.
!===============================================================================

  type MultipoleArray

    character(20) :: name                   ! name of nuclide, e.g. U235

    ! Isotope Properties
    logical                 :: fissionable  ! Is this isotope fissionable?
    complex(8), allocatable :: data(:,:)    ! Poles and residues
    real(8)                 :: sqrtAWR      ! Square root of the atomic
                                            ! weight ratio

    ! Windows
    integer :: fit_order                    ! Order of the fit. 1 linear,
                                            !   2 quadratic, etc.
    real(8) :: E_min                        ! Start energy for the windows
    real(8) :: E_max                        ! End energy for the windows
    real(8) :: spacing                      ! The actual spacing in sqrt(E)
                                            !   space.
    integer, allocatable :: windows(:, :)   ! Contains the indexes of the poles
                                            !   at the start and end of the
                                            !   window
    real(8), allocatable :: curvefit(:,:,:) ! Contains the fitting function.
                                            !   (reaction type, coeff index,
                                            !   window index)
    integer, allocatable :: broaden_poly(:) ! if 1, broaden, if 0, don't.

  contains
    procedure :: from_hdf5 => multipole_from_hdf5

  end type MultipoleArray

contains

!===============================================================================
! FROM_HDF5 loads multipole data from a HDF5 file.
!===============================================================================

  subroutine multipole_from_hdf5(this, group_id)
    class(MultipoleArray), intent(inout) :: this
    integer(HID_T),        intent(in)    :: group_id

    integer :: n_poles, n_residues, n_windows
    integer(HSIZE_T) :: dims_1d(1), dims_2d(2), dims_3d(3)
    integer(HID_T) :: dset

    ! Get name of nuclide from group
    this % name = get_name(group_id)

    ! Get rid of leading '/'
    this % name = trim(this % name(2:))

    ! Read scalar values.
    call read_dataset(this % spacing, group_id, "spacing")
    call read_dataset(this % sqrtAWR, group_id, "sqrtAWR")
    call read_dataset(this % E_min, group_id, "E_min")
    call read_dataset(this % E_max, group_id, "E_max")

    ! Read the "data" array.  Use its shape to figure out the number of poles
    ! and residue types in this data.
    dset = open_dataset(group_id, "data")
    call get_shape(dset, dims_2d)
    n_residues = int(dims_2d(1), 4) - 1
    n_poles = int(dims_2d(2), 4)
    allocate(this % data(n_residues+1, n_poles))
    if (n_poles > 0) call read_dataset(this % data, dset)
    call close_dataset(dset)

    ! Check to see if this data includes fission residues.
    this % fissionable = (n_residues == 3)

    ! Read the "windows" array and use its shape to figure out the number of
    ! windows.
    dset = open_dataset(group_id, "windows")
    call get_shape(dset, dims_2d)
    n_windows = int(dims_2d(2), 4)
    allocate(this % windows(dims_2d(1), n_windows))
    call read_dataset(this % windows, dset)
    call close_dataset(dset)

    ! Read the "broaden_poly" arrays.
    dset = open_dataset(group_id, "broaden_poly")
    call get_shape(dset, dims_1d)
    if (dims_1d(1) /= n_windows) call fatal_error("broaden_poly array shape is&
         &not consistent with the windows array shape in WMP library for"&
         // trim(this % name) // ".")
    allocate(this % broaden_poly(n_windows))
    call read_dataset(this % broaden_poly, dset)
    call close_dataset(dset)

    ! Read the "curvefit" array.
    dset = open_dataset(group_id, "curvefit")
    call get_shape(dset, dims_3d)
    if (dims_3d(3) /= n_windows) call fatal_error("curvefit array shape is not&
         &consistent with the windows array shape in WMP library for"&
         // trim(this % name) // ".")
    allocate(this % curvefit(dims_3d(1), dims_3d(2), dims_3d(3)))
    call read_dataset(this % curvefit, dset)
    call close_dataset(dset)
    this % fit_order = int(dims_3d(2), 4) - 1

  end subroutine multipole_from_hdf5

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
