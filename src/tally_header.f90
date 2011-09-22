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
     integer :: n_events    = 0
     real(8) :: val_history = 0.
     real(8) :: val         = 0.
     real(8) :: val_sq      = 0.
  end type TallyScore

!===============================================================================
! TALLYFILTER represents a single bin filter for a tally to score in. This could
! be, for example, a single cell OR a collection of cells.
!===============================================================================

  type TallyFilter
     integer              :: scalar = 0
     integer, allocatable :: array(:)
  end type TallyFilter

!===============================================================================
! TALLY describes a user-specified tally. The region of phase space to tally in
! is given by the TallyBins and the scores are stored in a TallyScore array.
!===============================================================================

  type TallyObject
     ! Basic data

     integer :: uid
     integer :: type
     real(8) :: volume
     integer :: cell_type
     integer :: reaction_type
     integer :: material_type

     ! Tally bin specifications

     type(TallyFilter), pointer :: cell_bins(:) => null()
     type(TallyFilter), pointer :: surface_bins(:) => null()
     type(TallyFilter), pointer :: universe_bins(:) => null()
     type(TallyFilter), pointer :: material_bins(:) => null()
     type(TallyFilter), pointer :: mesh_bins(:) => null()
     type(TallyFilter), pointer :: bornin_bins(:) => null()
     real(8), allocatable       :: energy_in(:)
     real(8), allocatable       :: energy_out(:)

     ! Number of bins for each filter

     integer :: n_cell_bins     = 0
     integer :: n_surface_bins  = 0
     integer :: n_universe_bins = 0
     integer :: n_material_bins = 0
     integer :: n_mesh_bins     = 0
     integer :: n_bornin_bins   = 0
     integer :: n_energy_in     = 0
     integer :: n_energy_out    = 0

     ! Macroscopic properties to score

     type(TallyFilter), pointer :: macro_bins(:) => null()
     integer :: n_macro_bins = 0
     
     ! Scores for each bin -- the most natural way to have scores would be to
     ! have a dimension for each different type of bin, but older Fortran
     ! standards are limited to 7 dimensions or less, so instead we map a
     ! combination of the filter bins into one integer and have that as the
     ! index into a one-dimensional array

     type(TallyScore), allocatable :: scores(:,:)

  end type TallyObject

end module tally_header
