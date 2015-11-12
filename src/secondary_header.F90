module secondary_header

  use endf_header, only: Tab1
  use interpolation, only: interpolate_tab1
  use random_lcg, only: prn

  type, abstract :: AngleEnergy
  contains
    procedure(iSampleAngleEnergy), deferred :: sample
  end type AngleEnergy

  type :: AngleEnergyContainer
    class(AngleEnergy), allocatable :: obj
  end type AngleEnergyContainer

  type :: SecondaryDistribution
    type(Tab1), allocatable :: applicability(:)
    type(AngleEnergyContainer), allocatable :: distribution(:)
  contains
    procedure :: sample => secondary_sample
  end type SecondaryDistribution

  abstract interface
    subroutine iSampleAngleEnergy(this, E_in, E_out, mu)
      import AngleEnergy
      class(AngleEnergy), intent(in) :: this
      real(8), intent(in) :: E_in
      real(8), intent(out) :: E_out
      real(8), intent(out) :: mu
    end subroutine iSampleAngleEnergy
  end interface

contains

  subroutine secondary_sample(this, E_in, E_out, mu)
    class(SecondaryDistribution), intent(in) :: this
    real(8), intent(in) :: E_in
    real(8), intent(out) :: E_out
    real(8), intent(out) :: mu

    real(8) :: p_valid

    do i = 1, size(this%applicability)
      ! Determine probability that i-th energy distribution is sampled
      p_valid = interpolate_tab1(this%applicability(i), E_in)

      ! If i-th distribution is sampled, sample energy from the distribution
      if (prn() <= p_valid) then
        call this%distribution(i)%obj%sample(E_in, E_out, mu)
        exit
      end if
    end do
  end subroutine secondary_sample

end module secondary_header
