module tally_filter_zernike

  use tally_filter_header

  implicit none

!===============================================================================
! ZERNIKEFILTER gives Zernike polynomial moments of a particle's position
!===============================================================================

  type, extends(CppTallyFilter) :: ZernikeFilter
  end type ZernikeFilter

!===============================================================================
! ZERNIKERADIALFILTER gives even order radial Zernike polynomial moments of a
! particle's position
!===============================================================================

  type, extends(ZernikeFilter) :: ZernikeRadialFilter
  end type ZernikeRadialFilter

end module tally_filter_zernike
