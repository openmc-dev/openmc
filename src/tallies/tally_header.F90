module tally_header

  use, intrinsic :: ISO_C_BINDING

  use hdf5

  use constants,           only: NONE, N_FILTER_TYPES, ZERO
  use error
  use dict_header,         only: DictIntInt
  use nuclide_header,      only: nuclide_dict
  use string,              only: to_lower, to_f_string
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

  integer(C_INT32_T), bind(C) :: n_tallies = 0 ! # of tallies

  type(TallyObject),     allocatable, target :: tallies(:)
  type(TallyDerivative), allocatable :: tally_derivs(:)
!$omp threadprivate(tally_derivs)

  ! Dictionary that maps user IDs to indices in 'tallies'
  type(DictIntInt) :: tally_dict

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

!===============================================================================
!                               C API FUNCTIONS
!===============================================================================

  function openmc_extend_tallies(n, index_start, index_end) result(err) bind(C)
    ! Extend the tallies array by n elements
    integer(C_INT32_T), value, intent(in) :: n
    integer(C_INT32_T), intent(out) :: index_start
    integer(C_INT32_T), intent(out) :: index_end
    integer(C_INT) :: err

    type(TallyObject), allocatable :: temp(:) ! temporary tallies array

    if (n_tallies == 0) then
      ! Allocate tallies array
      allocate(tallies(n))
    else
      ! Allocate tallies array with increased size
      allocate(temp(n_tallies + n))

      ! Copy original tallies to temporary array
      temp(1:n_tallies) = tallies

      ! Move allocation from temporary array
      call move_alloc(FROM=temp, TO=tallies)
    end if

    ! Return indices in tallies array
    index_start = n_tallies + 1
    index_end = n_tallies + n
    n_tallies = index_end

    err = 0
  end function openmc_extend_tallies


  function openmc_get_tally(id, index) result(err) bind(C)
    ! Returns the index in the tallies array of a tally with a given ID
    integer(C_INT32_T), value :: id
    integer(C_INT32_T), intent(out) :: index
    integer(C_INT) :: err

    if (allocated(tallies)) then
      if (tally_dict % has_key(id)) then
        index = tally_dict % get_key(id)
        err = 0
      else
        err = E_TALLY_INVALID_ID
      end if
    else
      err = E_TALLY_NOT_ALLOCATED
    end if
  end function openmc_get_tally


  function openmc_tally_get_id(index, id) result(err) bind(C)
    ! Return the ID of a tally
    integer(C_INT32_T), value       :: index
    integer(C_INT32_T), intent(out) :: id
    integer(C_INT) :: err

    if (index >= 1 .and. index <= size(tallies)) then
      id = tallies(index) % id
      err = 0
    else
      err = E_OUT_OF_BOUNDS
    end if
  end function openmc_tally_get_id


  function openmc_tally_get_nuclides(index, nuclides, n) result(err) bind(C)
    ! Return the list of nuclides assigned to a tally
    integer(C_INT32_T), value :: index
    type(C_PTR), intent(out) :: nuclides
    integer(C_INT), intent(out) :: n
    integer(C_INT) :: err

    err = E_UNASSIGNED
    if (index >= 1 .and. index <= size(tallies)) then
      associate (t => tallies(index))
        if (allocated(t % nuclide_bins)) then
          nuclides = C_LOC(t % nuclide_bins(1))
          n = size(t % nuclide_bins)
          err = 0
        end if
      end associate
    else
      err = E_OUT_OF_BOUNDS
    end if
  end function openmc_tally_get_nuclides


  function openmc_tally_results(index, ptr, shape_) result(err) bind(C)
    ! Returns a pointer to a tally results array along with its shape. This
    ! allows a user to obtain in-memory tally results from Python directly.
    integer(C_INT32_T), intent(in), value :: index
    type(C_PTR),        intent(out) :: ptr
    integer(C_INT),     intent(out) :: shape_(3)
    integer(C_INT) :: err

    err = E_UNASSIGNED
    if (index >= 1 .and. index <= size(tallies)) then
      if (allocated(tallies(index) % results)) then
        ptr = C_LOC(tallies(index) % results(1,1,1))
        shape_(:) = shape(tallies(index) % results)
        err = 0
      end if
    else
      err = E_OUT_OF_BOUNDS
    end if
  end function openmc_tally_results


  function openmc_tally_set_nuclides(index, n, nuclides) result(err) bind(C)
    ! Sets the nuclides in the tally which results should be scored for
    integer(C_INT32_T), value  :: index
    integer(C_INT), value      :: n
    type(C_PTR),    intent(in) :: nuclides(n)
    integer(C_INT) :: err

    integer :: i
    character(C_CHAR), pointer :: string(:)
    character(len=:, kind=C_CHAR), allocatable :: nuclide_

    err = E_UNASSIGNED
    if (index >= 1 .and. index <= size(tallies)) then
      associate (t => tallies(index))
        if (allocated(t % nuclide_bins)) deallocate(t % nuclide_bins)
        allocate(t % nuclide_bins(n))
        t % n_nuclide_bins = n

        do i = 1, n
          ! Convert C string to Fortran string
          call c_f_pointer(nuclides(i), string, [10])
          nuclide_ = to_lower(to_f_string(string))

          select case (nuclide_)
          case ('total')
            t % nuclide_bins(i) = -1
          case default
            if (nuclide_dict % has_key(nuclide_)) then
              t % nuclide_bins(i) = nuclide_dict % get_key(nuclide_)
            else
              err = E_NUCLIDE_NOT_LOADED
              return
            end if
          end select
        end do

        ! Recalculate total number of scoring bins
        t % total_score_bins = t % n_score_bins * t % n_nuclide_bins

        ! (Re)allocate results array
        if (allocated(t % results)) deallocate(t % results)
        allocate(t % results(3, t % total_score_bins, t % total_filter_bins))
        t % results(:,:,:) = ZERO

        err = 0
      end associate
    else
      err = E_OUT_OF_BOUNDS
    end if
  end function openmc_tally_set_nuclides

end module tally_header
