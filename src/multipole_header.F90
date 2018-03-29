module multipole_header

  use hdf5

  use constants
  use dict_header,      only: DictIntInt
  use error,            only: fatal_error
  use hdf5_interface

  implicit none

  !========================================================================
  ! Multipole related constants

  ! Formalisms
  integer, parameter :: FORM_MLBW = 2, &
                        FORM_RM   = 3, &
                        FORM_RML  = 7

  ! Constants that determine which value to access
  integer, parameter :: MP_EA = 1       ! Pole

  ! Reich-Moore indices
  integer, parameter :: RM_RT = 2, &    ! Residue total
                        RM_RA = 3, &    ! Residue absorption
                        RM_RF = 4       ! Residue fission

  ! Multi-level Breit Wigner indices
  integer, parameter :: MLBW_RT = 2, &  ! Residue total
                        MLBW_RX = 3, &  ! Residue compettitive
                        MLBW_RA = 4, &  ! Residue absorption
                        MLBW_RF = 5     ! Residue fission

  ! Polynomial fit indices
  integer, parameter :: FIT_T = 1, &    ! Total
                        FIT_A = 2, &    ! Absorption
                        FIT_F = 3       ! Fission

!===============================================================================
! MULTIPOLE contains all the components needed for the windowed multipole
! temperature dependent cross section libraries for the resolved resonance
! region.
!===============================================================================

  type MultipoleArray

    !=========================================================================
    ! Isotope Properties

    logical                 :: fissionable     ! Is this isotope fissionable?
    integer, allocatable    :: l_value(:)      ! The l index of the pole
    integer                 :: num_l           ! Number of unique l values
    real(8), allocatable    :: pseudo_k0RS(:)  ! The value (sqrt(2*mass neutron
                                               !   /reduced planck constant)
                                               !   * AWR/(AWR + 1)
                                               !   * scattering radius for
                                               !    each l
    complex(8), allocatable :: data(:,:)       ! Poles and residues
    real(8)                 :: sqrtAWR         ! Square root of the atomic
                                               ! weight ratio
    integer                 :: formalism       ! R-matrix formalism

    !=========================================================================
    ! Windows

    integer :: fit_order                    ! Order of the fit. 1 linear,
                                            !   2 quadratic, etc.
    real(8) :: start_E                      ! Start energy for the windows
    real(8) :: end_E                        ! End energy for the windows
    real(8) :: spacing                      ! The actual spacing in sqrt(E)
                                            !   space.
    ! spacing = sqrt(multipole_w % endE - multipole_w % startE)
    !           / multipole_w % windows
    integer, allocatable :: w_start(:)      ! Contains the index of the pole at
                                            !   the start of the window
    integer, allocatable :: w_end(:)        ! Contains the index of the pole at
                                            !   the end of the window
    real(8), allocatable :: curvefit(:,:,:) ! Contains the fitting function.
                                            !   (reaction type, coeff index,
                                            !   window index)

    integer, allocatable :: broaden_poly(:) ! if 1, broaden, if 0, don't.

  contains

    procedure :: from_hdf5 => multipole_from_hdf5

  end type MultipoleArray

contains

!===============================================================================
! FROM_HDF5 loads multipole data from an HDF5 file.
!===============================================================================

  subroutine multipole_from_hdf5(this, filename)
    class(MultipoleArray), intent(inout) :: this
    character(len=*),      intent(in)    :: filename

    character(len=10) :: version
    integer :: i, n_poles, n_residue_types, n_windows
    integer(HSIZE_T) :: dims_1d(1), dims_2d(2), dims_3d(3)
    integer(HID_T) :: file_id
    integer(HID_T) :: group_id
    integer(HID_T) :: dset
    type(DictIntInt) :: l_val_dict

    ! Open file for reading and move into the /isotope group
    file_id = file_open(filename, 'r', parallel=.true.)
    group_id = open_group(file_id, "/nuclide")

    ! Check the file version number.
    call read_dataset(version, file_id, "version")
    if (version /= VERSION_MULTIPOLE) call fatal_error("The current multipole&
         & format version is " // trim(VERSION_MULTIPOLE) // " but the file "&
         // trim(filename) // " uses version " // trim(version) // ".")

    ! Read scalar values.
    call read_dataset(this % formalism, group_id, "formalism")
    call read_dataset(this % spacing, group_id, "spacing")
    call read_dataset(this % sqrtAWR, group_id, "sqrtAWR")
    call read_dataset(this % start_E, group_id, "start_E")
    call read_dataset(this % end_E, group_id, "end_E")

    ! Read the "data" array.  Use its shape to figure out the number of poles
    ! and residue types in this data.
    dset = open_dataset(group_id, "data")
    call get_shape(dset, dims_2d)
    n_residue_types = int(dims_2d(1), 4) - 1
    n_poles = int(dims_2d(2), 4)
    allocate(this % data(n_residue_types+1, n_poles))
    call read_dataset(this % data, dset)
    call close_dataset(dset)

    ! Check to see if this data includes fission residues.
    if (this % formalism == FORM_RM) then
      this % fissionable = (n_residue_types == 3)
    else
      ! Assume FORM_MLBW.
      this % fissionable = (n_residue_types == 4)
    end if

    ! Read the "l_value" array.
    allocate(this % l_value(n_poles))
    call read_dataset(this % l_value, group_id, "l_value")

    ! Figure out the number of unique l values in the l_value array.
    do i = 1, n_poles
      if (.not. l_val_dict % has(this % l_value(i))) then
        call l_val_dict % set(this % l_value(i), 0)
      end if
    end do
    this % num_l = l_val_dict % size()
    call l_val_dict % clear()

    ! Read the "pseudo_K0RS" array.
    allocate(this % pseudo_k0RS(this % num_l))
    call read_dataset(this % pseudo_k0RS, group_id, "pseudo_K0RS")

    ! Read the "w_start" array and use its shape to figure out the number of
    ! windows.
    dset = open_dataset(group_id, "w_start")
    call get_shape(dset, dims_1d)
    n_windows = int(dims_1d(1), 4)
    allocate(this % w_start(n_windows))
    call read_dataset(this % w_start, dset)
    call close_dataset(dset)

    ! Read the "w_end" and "broaden_poly" arrays.
    allocate(this % w_end(n_windows))
    call read_dataset(this % w_end, group_id, "w_end")
    allocate(this % broaden_poly(n_windows))
    call read_dataset(this % broaden_poly, group_id, "broaden_poly")

    ! Read the "curvefit" array.
    dset = open_dataset(group_id, "curvefit")
    call get_shape(dset, dims_3d)
    allocate(this % curvefit(dims_3d(1), dims_3d(2), dims_3d(3)))
    call read_dataset(this % curvefit, dset)
    call close_dataset(dset)
    this % fit_order = int(dims_3d(2), 4) - 1

    ! Close the group and file.
    call close_group(group_id)
    call file_close(file_id)
  end subroutine multipole_from_hdf5
end module multipole_header
