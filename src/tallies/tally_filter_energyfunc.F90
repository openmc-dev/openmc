module tally_filter_energyfunc

  use tally_filter_header

  implicit none

!===============================================================================
! EnergyFunctionFilter multiplies tally scores by an arbitrary function of
! incident energy described by a piecewise linear-linear interpolation.
!===============================================================================

  type, extends(CppTallyFilter) :: EnergyFunctionFilter
  end type EnergyFunctionFilter

end module tally_filter_energyfunc
