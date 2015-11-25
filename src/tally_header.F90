module tally_header

  use constants,          only: NONE, N_FILTER_TYPES
  use trigger_header,     only: TriggerObject
  use, intrinsic :: ISO_C_BINDING

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

  type, bind(C) :: TallyResult
    real(C_DOUBLE) :: value    = 0.
    real(C_DOUBLE) :: sum      = 0.
    real(C_DOUBLE) :: sum_sq   = 0.
  end type TallyResult

!===============================================================================
! TALLYFILTER describes a filter that limits what events score to a tally. For
! example, a cell filter indicates that only particles in a specified cell
! should score to the tally.
!===============================================================================

  type TallyFilter
    integer :: type = NONE
    integer :: n_bins = 0
    integer :: offset = 0 ! Only used for distribcell filters
    integer, allocatable :: int_bins(:)
    real(8), allocatable :: real_bins(:) ! Only used for energy filters
  end type TallyFilter

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

    ! reset property - allows a tally to be reset after every batch
    logical :: reset = .false.

    ! Number of realizations of tally random variables
    integer :: n_realizations = 0

    ! Tally precision triggers
    integer                           :: n_triggers = 0  ! # of triggers
    type(TriggerObject),  allocatable :: triggers(:)     ! Array of triggers
  end type TallyObject

end module tally_header
