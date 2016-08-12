module tally_header

  use constants,           only: NONE, N_FILTER_TYPES, OTF_HEADROOM
  use tally_filter_header, only: TallyFilterContainer
  use dict_header,        only: DictIntInt
  use trigger_header,      only: TriggerObject

  use, intrinsic :: ISO_C_BINDING

  implicit none

!===============================================================================
! TALLYRESULT provides accumulation of results in a particular tally bin
!===============================================================================

  type, bind(C) :: TallyResult
    real(C_DOUBLE) :: value    = 0.
    real(C_DOUBLE) :: sum      = 0.
    real(C_DOUBLE) :: sum_sq   = 0.
  end type TallyResult

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
    type(TallyFilterContainer), allocatable :: filters(:)

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

    type(DictIntInt) :: filter_index_map          ! real_bin --> local_bin
    type(DictIntInt) :: reverse_filter_index_map  ! local_bin --> real_bin

    ! Whether or not this tally will do on-the-fly memory allocation
    logical :: on_the_fly_allocation = .false.

    ! Size of results array when variable
    integer :: size_results_filters

    ! Next index in the results array when variable
    integer :: next_filter_idx = 1

    ! reset property - allows a tally to be reset after every batch
    logical :: reset = .false.

    ! Number of realizations of tally random variables
    integer :: n_realizations = 0

    ! Tally precision triggers
    integer                           :: n_triggers = 0  ! # of triggers
    type(TriggerObject),  allocatable :: triggers(:)     ! Array of triggers

    ! Type-Bound procedures
    contains
      procedure :: get_filter_index => filter_index
      procedure :: otf_filter_index
      procedure :: grow_results_array
  end type TallyObject

  contains

!===============================================================================
! FILTER_INDEX returns the filter index in the results array for a given set of
! matching bins
!===============================================================================

    function filter_index(this, matching_bins) result(idx)

      class(TallyObject),   intent(inout) :: this
      integer, allocatable, intent(in)    :: matching_bins(:)

      integer :: idx

      ! Get index in total array
      idx = sum((matching_bins(1:size(this % filters)) - 1) * this % stride) + 1

      ! If the results array is fully allocated, this index is valid
      if (this % on_the_fly_allocation) then
        idx = this % otf_filter_index(idx)
      end if

    end function filter_index

!===============================================================================
! OTF_FILTER_INDEX returns the filter index in the results array when OTF tally
! allocation is active
!===============================================================================

    function otf_filter_index(this, real_bin) result(idx)

      class(TallyObject), intent(inout) :: this
      integer,            intent(in)    :: real_bin

      integer :: idx

      ! If we're doing on-the-fly memory allocation, we must use the map
      if (this % filter_index_map % has_key(real_bin)) then

        idx = this % filter_index_map % get_key(real_bin)

      else

!$omp critical (otf_tally_allocation)

        ! This is the first time this filter index has been needed

        ! Grow the results array if we've used it all
        if (this % next_filter_idx > this % size_results_filters) then
          call this % grow_results_array()
        end if

        ! Update the map
        call this % filter_index_map % add_key(real_bin, this % next_filter_idx)
        call this % reverse_filter_index_map % add_key( &
             this % next_filter_idx, real_bin)

        ! Set the return index
        idx = this % next_filter_idx

        ! Increment the next index
        this % next_filter_idx = this % next_filter_idx + 1
!$omp end critical (otf_tally_allocation)

      end if

    end function otf_filter_index

!===============================================================================
! GROW_RESULTS_ARRAY
!===============================================================================

    subroutine grow_results_array(this)

      class(TallyObject), intent(inout) :: this

      integer :: newsize
      type(TallyResult), allocatable :: temp(:,:)

      newsize = ceiling(real(this % size_results_filters, 8) * OTF_HEADROOM)

      ! Allocate results array with increased size
      allocate(temp(this % total_score_bins, newsize))

      ! Copy original results to temporary array
      temp(:, 1:this % size_results_filters) = &
           this % results(:, 1:this % size_results_filters)

      ! Move allocation from temporary array
      call move_alloc(FROM=temp, TO=this % results)

      ! Update size
      this % size_results_filters = newsize

    end subroutine grow_results_array

end module tally_header
