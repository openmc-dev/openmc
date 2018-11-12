module tally_filter_delayedgroup

  use, intrinsic :: ISO_C_BINDING

  use tally_filter_header

  implicit none
  private

!===============================================================================
! DELAYEDGROUPFILTER bins outgoing fission neutrons in their delayed groups.
! The get_all_bins functionality is not actually used.  The bins are manually
! iterated over in the scoring subroutines.
!===============================================================================

  type, public, extends(TallyFilter) :: DelayedGroupFilter
  contains
    procedure :: groups
  end type DelayedGroupFilter

contains

  function groups(this, i) result(group)
    class(DelayedGroupFilter), intent(in) :: this
    integer,                   intent(in) :: i
    integer                               :: group
    interface
      function delayedgroup_filter_groups(filt, i) result(group) bind(C)
        import C_PTR, C_INT
        type(C_PTR),    value :: filt
        integer(C_INT), value :: i
        integer(C_INT)        :: group
      end function
    end interface
    group = delayedgroup_filter_groups(this % ptr, i)
  end function groups

end module tally_filter_delayedgroup
