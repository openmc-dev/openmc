module tally_filter_universe

  use tally_filter_header

  implicit none

!===============================================================================
! UNIVERSEFILTER specifies which geometric universes tally events reside in.
!===============================================================================

  type, extends(CppTallyFilter) :: UniverseFilter
  end type UniverseFilter

end module tally_filter_universe
