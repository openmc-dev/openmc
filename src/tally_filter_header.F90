module tally_filter_header

  use constants,       only: MAX_LINE_LEN
  use particle_header, only: Particle

  use hdf5

  implicit none

!===============================================================================
! TALLYFILTER describes a filter that limits what events score to a tally. For
! example, a cell filter indicates that only particles in a specified cell
! should score to the tally.
!===============================================================================

  type, abstract :: TallyFilter
    integer :: n_bins = 0
  contains
    procedure(get_next_bin_),  deferred :: get_next_bin
    procedure(to_statepoint_), deferred :: to_statepoint
    procedure                           :: to_summary => to_summary_
    procedure(text_label_),    deferred :: text_label
    procedure                           :: initialize => initialize_
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

    subroutine get_next_bin_(this, p, estimator, current_bin, next_bin, weight)
      import TallyFilter
      import Particle
      class(TallyFilter), intent(in)  :: this
      type(Particle),     intent(in)  :: p
      integer,            intent(in)  :: estimator
      integer,            intent(in)  :: current_bin
      integer,            intent(out) :: next_bin
      real(8),            intent(out) :: weight
    end subroutine get_next_bin_

!===============================================================================
! TO_STATPEOINT writes all the information needed to reconstruct the filter to
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
! TO_SUMMARY writes all the information needed to reconstruct the filter to the
! given filter_group.  If this procedure is not overridden by the derived class,
! then it will call to_statepoint by default.

  subroutine to_summary_(this, filter_group)
    class(TallyFilter), intent(in) :: this
    integer(HID_T),     intent(in) :: filter_group

    call this % to_statepoint(filter_group)
  end subroutine to_summary_

!===============================================================================
! INITIALIZE sets up any internal data, as necessary.  If this procedure is not
! overriden by the derived class, then it will do nothing by default.

  subroutine initialize_(this)
    class(TallyFilter), intent(inout) :: this
  end subroutine initialize_

end module tally_filter_header
