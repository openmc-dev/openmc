module tally_header

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
! TALLYMAPITEM
!===============================================================================

  type TallyMapItem
     type(TallyMapElement), allocatable :: elements(:)
  end type TallyMapItem

!===============================================================================
! TALLYMAP contains a list of pairs of indices to tallies and the corresponding
! bin for a given filter.
!===============================================================================

  type TallyMap
     type(TallyMapItem), allocatable :: items(:)
  end type TallyMap

!===============================================================================
! TALLYSCORE provides accumulation of scores in a particular tally bin
!===============================================================================

  type TallyScore
     integer :: n_events = 0
     real(8) :: value    = 0.
     real(8) :: sum      = 0.
     real(8) :: sum_sq   = 0.
  end type TallyScore

!===============================================================================
! TALLYOBJECT describes a user-specified tally. The region of phase space to
! tally in is given by the TallyBins and the scores are stored in a TallyScore
! array.
!===============================================================================

  type TallyObject
     ! Basic data

     integer :: id               ! user-defined identifier
     character(len=52) :: label  ! user-defined label
     integer :: type             ! volume, surface current
     integer :: estimator        ! collision, track-length
     real(8) :: volume           ! volume of region

     ! Information about what filters should be used

     integer              :: n_filters
     integer, allocatable :: filters(:)

     ! Filter bin specifications

     integer, allocatable :: universe_bins(:)
     integer, allocatable :: material_bins(:)
     integer, allocatable :: cell_bins(:)
     integer, allocatable :: cellborn_bins(:)
     integer, allocatable :: surface_bins(:)
     integer              :: mesh = 0
     real(8), allocatable :: energy_in(:)
     real(8), allocatable :: energy_out(:)

     ! Total number of filter bins
     integer :: n_total_bins = 0

     ! The following attributes do not necessarily need to be stored but they
     ! greatly simplify logic in many places. n_bins gives the number of bins
     ! for each filter type, e.g. n_filter_bins(FILTER_CELL) would be the size
     ! of cell_bins. The stride attribute is used for determining the index in
     ! the scores array for a bin combination. Since multiple dimensions are
     ! mapped onto one dimension in the scores array, the stride attribute gives
     ! the stride for a given filter type within the scores array

     integer, allocatable :: n_filter_bins(:)
     integer, allocatable :: stride(:)

     ! Individual nuclides to tally
     integer, allocatable :: nuclide_bins(:)
     integer :: n_nuclide_bins = 0
     logical :: all_nuclides = .false.

     ! Values to score, e.g. flux, absorption, etc.
     integer, allocatable :: score_bins(:)
     integer :: n_score_bins = 0
     
     ! Scores for each bin -- the first dimenion of the array is for scores
     ! (e.g. flux, total reaction rate, fission reaction rate, etc.) and the
     ! second dimension of the array is for the combination of filters
     ! (e.g. specific cell, specific energy group, etc.)

     type(TallyScore), allocatable :: scores(:,:)

  end type TallyObject

end module tally_header
