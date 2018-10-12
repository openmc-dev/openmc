module tally_filter_cellborn

  use tally_filter_header

  implicit none

!===============================================================================
! CELLBORNFILTER specifies which cell the particle was born in.
!===============================================================================

  type, public, extends(CppTallyFilter) :: CellbornFilter
  end type CellbornFilter

end module tally_filter_cellborn
