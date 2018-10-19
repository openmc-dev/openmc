module tally_filter_sptl_legendre

  use tally_filter_header

  implicit none

!===============================================================================
! SPATIALLEGENDREFILTER gives Legendre moments of the particle's normalized
! position along an axis
!===============================================================================

  type, extends(CppTallyFilter) :: SpatialLegendreFilter
  end type SpatialLegendreFilter

end module tally_filter_sptl_legendre
