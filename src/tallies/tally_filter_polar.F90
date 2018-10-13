module tally_filter_polar

  use tally_filter_header

  implicit none

!===============================================================================
! POLARFILTER bins the incident neutron polar angle (relative to the global
! z-axis).
!===============================================================================

  type, extends(CppTallyFilter) :: PolarFilter
  end type PolarFilter

end module tally_filter_polar
