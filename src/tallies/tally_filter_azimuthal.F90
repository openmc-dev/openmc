module tally_filter_azimuthal

  use tally_filter_header

  implicit none

!===============================================================================
! AZIMUTHALFILTER bins the incident neutron azimuthal angle (relative to the
! global xy-plane).
!===============================================================================

  type, extends(CppTallyFilter) :: AzimuthalFilter
  end type AzimuthalFilter

end module tally_filter_azimuthal
