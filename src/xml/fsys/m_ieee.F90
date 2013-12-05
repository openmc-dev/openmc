module m_ieee

  implicit none
  private

  public :: generate_nan

contains

  function generate_nan() result(nan)
    real :: nan
    real :: zero
    zero = 0.0
    nan = 0.0/zero
  end function generate_nan 

end module m_ieee 
