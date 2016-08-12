module angleenergy_header

  use hdf5, only: HID_T

!===============================================================================
! ANGLEENERGY (abstract) defines a correlated or uncorrelated angle-energy
! distribution that is a function of incoming energy. Each derived type must
! implement a sample() subroutine that returns an outgoing energy and scattering
! cosine given an incoming energy.
!===============================================================================

  type, abstract :: AngleEnergy
  contains
    procedure(angleenergy_sample_), deferred :: sample
    procedure(angleenergy_from_hdf5_), deferred :: from_hdf5
  end type AngleEnergy

  abstract interface
    subroutine angleenergy_sample_(this, E_in, E_out, mu)
      import AngleEnergy
      class(AngleEnergy), intent(in) :: this
      real(8), intent(in) :: E_in
      real(8), intent(out) :: E_out
      real(8), intent(out) :: mu
    end subroutine angleenergy_sample_

    subroutine angleenergy_from_hdf5_(this, group_id)
      import AngleEnergy, HID_T
      class(AngleEnergy), intent(inout) :: this
      integer(HID_T), intent(in) :: group_id
    end subroutine angleenergy_from_hdf5_
  end interface

  type :: AngleEnergyContainer
    class(AngleEnergy), allocatable :: obj
  end type AngleEnergyContainer

end module angleenergy_header
