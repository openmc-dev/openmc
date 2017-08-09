module tally_header

  use, intrinsic :: ISO_C_BINDING

  use hdf5

  use constants,           only: NONE, N_FILTER_TYPES
  use tally_filter_header, only: TallyFilterContainer
  use trigger_header,      only: TriggerObject

  implicit none

!===============================================================================
! TALLYDERIVATIVE describes a first-order derivative that can be applied to
! tallies.
!===============================================================================

  type TallyDerivative
    integer :: id
    integer :: variable
    integer :: diff_material
    integer :: diff_nuclide
    real(8) :: flux_deriv
  end type TallyDerivative

!===============================================================================
! TALLYOBJECT describes a user-specified tally. The region of phase space to
! tally in is given by the TallyFilters and the results are stored in a
! TallyResult array.
!===============================================================================

  type TallyObject
    ! Basic data

    integer :: id                   ! user-defined identifier
    character(len=104) :: name = "" ! user-defined name
    integer :: type                 ! volume, surface current
    integer :: estimator            ! collision, track-length
    real(8) :: volume               ! volume of region
    integer, allocatable :: filter(:) ! index in filters array

    ! The stride attribute is used for determining the index in the results
    ! array for a matching_bin combination. Since multiple dimensions are
    ! mapped onto one dimension in the results array, the stride attribute gives
    ! the stride for a given filter type within the results array

    integer, allocatable :: stride(:)

    ! This array provides a way to lookup what index in the filters array a
    ! certain filter is. For example, if find_filter(FILTER_CELL) > 0, then the
    ! value is the index in filters(:).

    integer :: find_filter(N_FILTER_TYPES) = 0

    ! Individual nuclides to tally
    integer              :: n_nuclide_bins = 0
    integer, allocatable :: nuclide_bins(:)
    logical              :: all_nuclides = .false.

    ! Values to score, e.g. flux, absorption, etc.
    ! scat_order is the scattering order for each score.
    ! It is to be 0 if the scattering order is 0, or if the score is not a
    ! scattering response.
    integer              :: n_score_bins = 0
    integer, allocatable :: score_bins(:)
    integer, allocatable :: moment_order(:)
    integer              :: n_user_score_bins = 0

    ! Results for each bin -- the first dimension of the array is for scores
    ! (e.g. flux, total reaction rate, fission reaction rate, etc.) and the
    ! second dimension of the array is for the combination of filters
    ! (e.g. specific cell, specific energy group, etc.)

    integer :: total_filter_bins
    integer :: total_score_bins
    real(C_DOUBLE), allocatable :: results(:,:,:)

    ! reset property - allows a tally to be reset after every batch
    logical :: reset = .false.

    ! Number of realizations of tally random variables
    integer :: n_realizations = 0

    ! Tally precision triggers
    integer                           :: n_triggers = 0  ! # of triggers
    type(TriggerObject),  allocatable :: triggers(:)     ! Array of triggers

    ! Index for the TallyDerivative for differential tallies.
    integer :: deriv = NONE

  contains
    procedure :: write_results_hdf5
    procedure :: read_results_hdf5
  end type TallyObject

contains

  subroutine write_results_hdf5(this, group_id)
    class(TallyObject), intent(in) :: this
    integer(HID_T),     intent(in) :: group_id

    integer :: hdf5_err
    integer(HID_T) :: dset, dspace
    integer(HID_T) :: memspace
    integer(HSIZE_T) :: dims(3)
    integer(HSIZE_T) :: dims_slab(3)
    integer(HSIZE_T) :: offset(3) = [1,0,0]

    ! Create file dataspace
    dims_slab(:) = shape(this % results)
    dims_slab(1) = 2
    call h5screate_simple_f(3, dims_slab, dspace, hdf5_err)

    ! Create memory dataspace that contains only SUM and SUM_SQ values
    dims(:) = shape(this % results)
    call h5screate_simple_f(3, dims, memspace, hdf5_err)
    call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, offset, dims_slab, &
         hdf5_err)

    ! Create and write to dataset
    call h5dcreate_f(group_id, "results", H5T_NATIVE_DOUBLE, dspace, dset, &
         hdf5_err)
    call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, this % results, dims_slab, &
         hdf5_err, mem_space_id=memspace)

    ! Close identifiers
    call h5dclose_f(dset, hdf5_err)
    call h5sclose_f(memspace, hdf5_err)
    call h5sclose_f(dspace, hdf5_err)
  end subroutine write_results_hdf5

  subroutine read_results_hdf5(this, group_id)
    class(TallyObject), intent(inout) :: this
    integer(HID_T),     intent(in) :: group_id

    integer :: hdf5_err
    integer(HID_T) :: dset, dspace
    integer(HID_T) :: memspace
    integer(HSIZE_T) :: dims(3)
    integer(HSIZE_T) :: dims_slab(3)
    integer(HSIZE_T) :: offset(3) = [1,0,0]

    ! Create file dataspace
    dims_slab(:) = shape(this % results)
    dims_slab(1) = 2
    call h5screate_simple_f(3, dims_slab, dspace, hdf5_err)

    ! Create memory dataspace that contains only SUM and SUM_SQ values
    dims(:) = shape(this % results)
    call h5screate_simple_f(3, dims, memspace, hdf5_err)
    call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, offset, dims_slab, &
         hdf5_err)

    ! Create and write to dataset
    call h5dopen_f(group_id, "results", dset, hdf5_err)
    call h5dread_f(dset, H5T_NATIVE_DOUBLE, this % results, dims_slab, &
         hdf5_err, mem_space_id=memspace)

    ! Close identifiers
    call h5dclose_f(dset, hdf5_err)
    call h5sclose_f(memspace, hdf5_err)
    call h5sclose_f(dspace, hdf5_err)
  end subroutine read_results_hdf5

end module tally_header
