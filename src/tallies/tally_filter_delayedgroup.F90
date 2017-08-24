module tally_filter_delayedgroup

  use, intrinsic :: ISO_C_BINDING

  use hdf5

  use constants,          only: ONE, MAX_LINE_LEN
  use error,              only: fatal_error
  use hdf5_interface
  use particle_header,    only: Particle
  use string,             only: to_str
  use tally_filter_header

  implicit none
  private

!===============================================================================
! DELAYEDGROUPFILTER bins outgoing fission neutrons in their delayed groups.
! The get_all_bins functionality is not actually used.  The bins are manually
! iterated over in the scoring subroutines.
!===============================================================================

  type, public, extends(TallyFilter) :: DelayedGroupFilter
    integer, allocatable :: groups(:)
  contains
    procedure :: get_all_bins => get_all_bins_dg
    procedure :: to_statepoint => to_statepoint_dg
    procedure :: text_label => text_label_dg
  end type DelayedGroupFilter

contains

  subroutine get_all_bins_dg(this, p, estimator, match)
    class(DelayedGroupFilter), intent(in)  :: this
    type(Particle),            intent(in)  :: p
    integer,                   intent(in)  :: estimator
    type(TallyFilterMatch),         intent(inout) :: match

    call match % bins % push_back(1)
    call match % weights % push_back(ONE)
  end subroutine get_all_bins_dg

  subroutine to_statepoint_dg(this, filter_group)
    class(DelayedGroupFilter), intent(in) :: this
    integer(HID_T),            intent(in) :: filter_group

    call write_dataset(filter_group, "type", "delayedgroup")
    call write_dataset(filter_group, "n_bins", this % n_bins)
    call write_dataset(filter_group, "bins", this % groups)
  end subroutine to_statepoint_dg

  function text_label_dg(this, bin) result(label)
    class(DelayedGroupFilter), intent(in) :: this
    integer,                   intent(in) :: bin
    character(MAX_LINE_LEN)               :: label

    label = "Delayed Group " // to_str(this % groups(bin))
  end function text_label_dg

end module tally_filter_delayedgroup
