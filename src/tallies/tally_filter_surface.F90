module tally_filter_surface

  use tally_filter_header

  implicit none

!===============================================================================
! SURFACEFILTER specifies which surface particles are crossing
!===============================================================================

  type, extends(CppTallyFilter) :: SurfaceFilter
    ! True if this filter is used for surface currents
    logical              :: current = .false.
  end type SurfaceFilter

end module tally_filter_surface
