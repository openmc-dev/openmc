module tally_header

  implicit none

!===============================================================================
! TALLYSCORE provides accumulation of scores in a particular tally bin
!===============================================================================

  type TallyScore
     integer :: n_events
     real(8) :: val
     real(8) :: val_sq
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

     ! Macroscopic properties to score

     type(TallyFilter), pointer :: macro_bins(:) => null()
     
     ! Scores for each bin -- the most natural way to have scores would be to
     ! have a dimension for each different type of bin, but older Fortran
     ! standards are limited to 7 dimensions or less, so instead we map a
     ! combination of the bins into one integer and have that as the index into
     ! a one-dimensional array

     type(TallyScore), allocatable :: score(:)

  end type TallyObject

end module tally_header
