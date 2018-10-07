module tally_filter_cellfrom

  use, intrinsic :: ISO_C_BINDING

  use tally_filter_cell

  implicit none
  private

!===============================================================================
! CELLFROMFILTER specifies which geometric cells particles exit when crossing a
! surface.
!===============================================================================

  type, public, extends(CellFilter) :: CellFromFilter
  end type CellFromFilter

end module tally_filter_cellfrom
