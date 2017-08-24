module tally_filter_delayedgroup

  use, intrinsic :: ISO_C_BINDING

  use hdf5

  use constants,          only: ONE, MAX_LINE_LEN, MAX_DELAYED_GROUPS
  use error,              only: fatal_error
  use hdf5_interface
  use particle_header,    only: Particle
  use string,             only: to_str
  use tally_filter_header
  use xml_interface

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
    procedure :: from_xml
    procedure :: get_all_bins => get_all_bins_dg
    procedure :: to_statepoint => to_statepoint_dg
    procedure :: text_label => text_label_dg
  end type DelayedGroupFilter

contains

  subroutine from_xml(this, node)
    class(DelayedGroupFilter), intent(inout) :: this
    type(XMLNode), intent(in) :: node

    integer :: i
    integer :: n

    n = node_word_count(node, "bins")

    ! Allocate and store bins
    this % n_bins = n
    allocate(this % groups(n))
    call get_node_array(node, "bins", this % groups)

    ! Check that bins are all are between 1 and MAX_DELAYED_GROUPS
    do i = 1, n
      if (this % groups(i) < 1 .or. &
           this % groups(i) > MAX_DELAYED_GROUPS) then
        call fatal_error("Encountered delayedgroup bin with index " &
             // trim(to_str(this % groups(i))) // " that is outside &
             &the range of 1 to MAX_DELAYED_GROUPS ( " &
             // trim(to_str(MAX_DELAYED_GROUPS)) // ")")
      end if
    end do
  end subroutine from_xml

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
