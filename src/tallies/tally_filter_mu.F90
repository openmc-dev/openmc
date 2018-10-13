module tally_filter_mu

  use tally_filter_header

  implicit none

!===============================================================================
! MUFILTER bins the incoming-outgoing direction cosine.  This is only used for
! scatter reactions.
!===============================================================================

  type, extends(CppTallyFilter) :: MuFilter
  end type MuFilter

end module tally_filter_mu
