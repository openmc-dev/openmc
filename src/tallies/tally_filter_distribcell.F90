module tally_filter_distribcell

  use tally_filter_header

  implicit none
  private

!===============================================================================
! DISTRIBCELLFILTER specifies which distributed geometric cells tally events
! reside in.
!===============================================================================

  type, public, extends(CppTallyFilter) :: DistribcellFilter
  contains
    procedure :: initialize => initialize_distribcell
  end type DistribcellFilter

contains

  subroutine initialize_distribcell(this)
    class(DistribcellFilter), intent(inout) :: this
    call this % initialize_cpp_inner()
    this % n_bins = this % n_bins_cpp()
  end subroutine initialize_distribcell

end module tally_filter_distribcell
