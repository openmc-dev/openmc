module tally_filter_cell

  use tally_filter_header

  implicit none

!===============================================================================
! CELLFILTER specifies which geometric cells tally events reside in.
!===============================================================================

  type, extends(CppTallyFilter) :: CellFilter
  end type CellFilter

end module tally_filter_cell
