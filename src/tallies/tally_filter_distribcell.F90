module tally_filter_distribcell

  use tally_filter_header

  implicit none
  private

  type, public, extends(TallyFilter) :: DistribcellFilter
  contains
    procedure :: initialize => initialize_distribcell
  end type DistribcellFilter

contains

  subroutine initialize_distribcell(this)
    class(DistribcellFilter), intent(inout) :: this
    call this % initialize_cpp()
    this % n_bins = this % n_bins_cpp()
  end subroutine initialize_distribcell

end module tally_filter_distribcell
