module tally_filter_header

  use constants,       only: MAX_LINE_LEN
  use particle_header, only: Particle
  use stl_vector,      only: VectorInt, VectorReal

  use hdf5

  implicit none

!===============================================================================
! TALLYFILTERMATCH stores every valid bin and weight for a filter
!===============================================================================

  type TallyFilterMatch
    ! Index of the bin and weight being used in the current filter combination
    integer          :: i_bin
    type(VectorInt)  :: bins
    type(VectorReal) :: weights

    ! Indicates whether all valid bins for this filter have been found
    logical          :: bins_present = .false.
  end type TallyFilterMatch

!===============================================================================
! TALLYFILTER describes a filter that limits what events score to a tally. For
! example, a cell filter indicates that only particles in a specified cell
! should score to the tally.
!===============================================================================

  type, abstract :: TallyFilter
    integer :: id
    integer :: n_bins = 0
  contains
    procedure(get_all_bins_),  deferred :: get_all_bins
    procedure(to_statepoint_), deferred :: to_statepoint
    procedure(text_label_),    deferred :: text_label
    procedure                           :: initialize => filter_initialize
  end type TallyFilter

  abstract interface

!===============================================================================
! GET_NEXT_BIN gives the index for the next valid filter bin and a weight that
! will be applied to the flux.
!
! In principle, a filter can have multiple valid bins.  If current_bin =
! NO_BIN_FOUND, then this method should give the first valid bin.  Providing the
! first valid bin should then give the second valid bin, and so on.  When there
! are no valid bins left, the next_bin should be NO_VALID_BIN.

    subroutine get_all_bins_(this, p, estimator, match)
      import TallyFilter
      import Particle
      import TallyFilterMatch
      class(TallyFilter), intent(in)  :: this
      type(Particle),     intent(in)  :: p
      integer,            intent(in)  :: estimator
      type(TallyFilterMatch), intent(inout) :: match
    end subroutine get_all_bins_

!===============================================================================
! TO_STATEPOINT writes all the information needed to reconstruct the filter to
! the given filter_group.

    subroutine to_statepoint_(this, filter_group)
      import TallyFilter
      import HID_T
      class(TallyFilter), intent(in) :: this
      integer(HID_T),     intent(in) :: filter_group
    end subroutine to_statepoint_

!===============================================================================
! TEXT_LABEL returns a string describing the given filter bin.  For example, an
! energy filter might return the string "Incoming Energy [0.625E-6, 20.0)".
! This is used to write the tallies.out file.

    function text_label_(this, bin) result(label)
      import TallyFilter
      import MAX_LINE_LEN
      class(TallyFilter), intent(in) :: this
      integer,            intent(in) :: bin
      character(MAX_LINE_LEN)        :: label
    end function text_label_

  end interface

!===============================================================================
! TALLYFILTERCONTAINER contains an allocatable TallyFilter object for arrays of
! TallyFilters
!===============================================================================

  type TallyFilterContainer
    class(TallyFilter), allocatable :: obj
  end type TallyFilterContainer

contains

!===============================================================================
! INITIALIZE sets up any internal data, as necessary.  If this procedure is not
! overriden by the derived class, then it will do nothing by default.

  subroutine filter_initialize(this)
    class(TallyFilter), intent(inout) :: this
  end subroutine filter_initialize

end module tally_filter_header
