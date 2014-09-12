module tally_header

  use constants,    only: NONE, N_FILTER_TYPES
  use dict_header,  only: DictIntInt

  implicit none

!===============================================================================
! TALLYMAPELEMENT gives an index to a tally which is to be scored and the
! corresponding bin for the filter variable
!===============================================================================

  type TallyMapElement
    integer :: index_tally
    integer :: index_bin
  end type TallyMapElement

!===============================================================================
! TALLYMAPITEM contains a list of tally/bin combinations for each mappable
! filter bin specified.
!===============================================================================

  type TallyMapItem
    type(TallyMapElement), allocatable :: elements(:)
  end type TallyMapItem

!===============================================================================
! TALLYMAP contains a list of pairs of indices to tallies and the corresponding
! bin for a given filter. There is one TallyMap for each mappable filter
! type. The items array is as long as the corresponding array for that filter,
! e.g. for tally_maps(FILTER_CELL), items is n_cells long.
!===============================================================================

  type TallyMap
    type(TallyMapItem), allocatable :: items(:)
  end type TallyMap

!===============================================================================
! TALLYRESULT provides accumulation of results in a particular tally bin
!===============================================================================

  type TallyResult
    real(8) :: value    = 0.
    real(8) :: sum      = 0.
    real(8) :: sum_sq   = 0.
  end type TallyResult

!===============================================================================
! TALLYFILTER describes a filter that limits what events score to a tally. For
! example, a cell filter indicates that only particles in a specified cell
! should score to the tally.
!===============================================================================

  type TallyFilter
    integer :: type = NONE
    integer :: n_bins = 0
    integer, allocatable :: int_bins(:)
    real(8), allocatable :: real_bins(:) ! Only used for energy filters

    ! Type-Bound procedures
    contains
      procedure :: clear => tallyfilter_clear ! Deallocates TallyFilter
  end type TallyFilter

!===============================================================================
! TALLYOBJECT describes a user-specified tally. The region of phase space to
! tally in is given by the TallyFilters and the results are stored in a
! TallyResult array.
!===============================================================================

  type TallyObject
    ! Basic data

    integer :: id                   ! user-defined identifier
    character(len=52) :: label = "" ! user-defined label
    integer :: type                 ! volume, surface current
    integer :: estimator            ! collision, track-length
    real(8) :: volume               ! volume of region

    ! Information about what filters should be used

    integer                        :: n_filters    ! Number of filters
    type(TallyFilter), allocatable :: filters(:)   ! Filter data (type/bins)

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
    type(TallyResult), allocatable :: results(:,:)

    ! For domain-decomposed runs, the results array is allocated on the fly
    ! as filter bins need to be accessed.  In this case, the filter indices
    ! calculated by matching bins will not actually correspond to the index
    ! in the results array that stores the desired values.  The following
    ! dictionary stores a mapping to keep track of these differences.

    type(DictIntInt) :: filter_index_map

    ! Whether or not this tally will do on-the-fly memory allocation
    logical :: on_the_fly_allocation = .false.

    ! Size of results array when variable
    integer :: size_results_filters

    ! reset property - allows a tally to be reset after every batch
    logical :: reset = .false.

    ! Number of realizations of tally random variables
    integer :: n_realizations = 0

    ! Type-Bound procedures
    contains
      procedure :: get_filter_index => filter_index
      procedure :: clear => tallyobject_clear ! Deallocates TallyObject
  end type TallyObject

  contains

!===============================================================================
! FILTER_INDEX returns the filter index in the results array for a given set of
! matching bins
!===============================================================================

    function filter_index(this, matching_bins) result(idx)

      use global, only: rank

      class(TallyObject), intent(inout) :: this 
      integer, allocatable, intent(in) :: matching_bins(:)
      
      integer :: idx
      type(TallyResult), allocatable :: temp(:,:)

      ! Get index in total array
      idx = sum((matching_bins(1:this % n_filters) - 1) * this % stride) + 1

      ! If the results array is fully allocated, this index is valid
      if (.not. this % on_the_fly_allocation) return

      ! If we're doing on-the-fly memory allocation, we must use the map
      if (this % filter_index_map % has_key(idx)) then

        idx = this % filter_index_map % get_key(idx)

      else
        
        print *, rank, 'allocating bin', idx

        ! This is the first time this filter index has been needed, so we must 
        ! allocate TallyResult object for it
  
        ! TODO: this should be done with a linked list

        ! Allocate results array with increased size
        allocate(temp(this % total_score_bins, this % size_results_filters + 1))

        ! Copy original results to temporary array
        temp(:, 1:this % size_results_filters) = &
             this % results(:, 1:this % size_results_filters)

        ! Move allocation from temporary array
        call move_alloc(FROM=temp, TO=this % results)

        ! Increment counter and add to mapping        
        this % size_results_filters = this % size_results_filters + 1
        call this % filter_index_map % add_key(idx, this % size_results_filters)
        idx = this % size_results_filters

      end if

    end function filter_index

!===============================================================================
! TALLYFILTER_CLEAR deallocates a TallyFilter element and sets it to its as
! initialized state.
!===============================================================================

    subroutine tallyfilter_clear(this)
      class(TallyFilter), intent(inout) :: this ! The TallyFilter to be cleared

      this % type = NONE
      this % n_bins = 0
      if (allocated(this % int_bins)) &
           deallocate(this % int_bins)
      if (allocated(this % real_bins)) &
           deallocate(this % real_bins)

    end subroutine tallyfilter_clear

!===============================================================================
! TALLYOBJECT_CLEAR deallocates a TallyObject element and sets it to its as
! initialized state.
!===============================================================================

    subroutine tallyobject_clear(this)
      class(TallyObject), intent(inout) :: this ! The TallyObject to be cleared

      integer :: i  ! Loop Index

      ! This routine will go through each item in TallyObject and set the value
      ! to its default, as-initialized values, including deallocations.
      this % label = ""

      if (allocated(this % filters)) then
        do i = 1, size(this % filters)
          call this % filters(i) % clear()
        end do
        deallocate(this % filters)
      end if

      if (allocated(this % stride)) &
           deallocate(this % stride)

      this % find_filter = 0

      this % n_nuclide_bins = 0
      if (allocated(this % nuclide_bins)) &
           deallocate(this % nuclide_bins)
      this % all_nuclides = .false.

      this % n_score_bins = 0
      if (allocated(this % score_bins)) &
           deallocate(this % score_bins)
      if (allocated(this % moment_order)) &
           deallocate(this % moment_order)
      this % n_user_score_bins = 0

      if (allocated(this % results)) &
           deallocate(this % results)

      this % reset = .false.

      this % n_realizations = 0

    end subroutine tallyobject_clear

end module tally_header
