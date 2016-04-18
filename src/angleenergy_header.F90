module angleenergy_header

!===============================================================================
! ANGLEENERGY (abstract) defines a correlated or uncorrelated angle-energy
! distribution that is a function of incoming energy. Each derived type must
! implement a sample() subroutine that returns an outgoing energy and scattering
! cosine given an incoming energy.
!===============================================================================

  type, abstract :: AngleEnergy
  contains
    procedure(angleenergy_sample_), deferred :: sample
  end type AngleEnergy

  abstract interface
    subroutine angleenergy_sample_(this, E_in, E_out, mu)
      import AngleEnergy
      class(AngleEnergy), intent(in) :: this
      real(8), intent(in) :: E_in
      real(8), intent(out) :: E_out
      real(8), intent(out) :: mu
    end subroutine angleenergy_sample_
  end interface

  type :: AngleEnergyContainer
    class(AngleEnergy), allocatable :: obj
  end type AngleEnergyContainer

end module angleenergy_header
