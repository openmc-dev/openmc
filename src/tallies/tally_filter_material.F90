module tally_filter_material

  use tally_filter_header

  implicit none

!===============================================================================
! MATERIAL specifies which material tally events reside in.
!===============================================================================

  type, extends(CppTallyFilter) :: MaterialFilter
  end type MaterialFilter

end module tally_filter_material
