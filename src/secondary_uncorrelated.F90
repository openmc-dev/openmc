module secondary_uncorrelated

  use angle_distribution, only: AngleDistribution
  use constants, only: ONE, TWO
  use energy_distribution, only: EnergyDistribution
  use secondary_header, only: AngleEnergy
  use random_lcg, only: prn

  type, extends(AngleEnergy) :: UncorrelatedAngleEnergy
    type(AngleDistribution) :: angle
    class(EnergyDistribution), allocatable :: energy
  contains
    procedure :: sample => uncorrelated_sample
  end type UncorrelatedAngleEnergy

contains

  subroutine uncorrelated_sample(this, E_in, E_out, mu)
    class(UncorrelatedAngleEnergy), intent(in) :: this
    real(8), intent(in) :: E_in
    real(8), intent(out) :: E_out
    real(8), intent(out) :: mu

    ! Sample cosine of scattering angle
    if (allocated(this%angle%energy)) then
      mu = this%angle%sample(E_in)
    else
      ! no angle distribution given => assume isotropic for all energies
      mu = TWO*prn() - ONE
    end if

    ! Sample outgoing energy
    E_out = this%energy%sample(E_in)
  end subroutine uncorrelated_sample

end module secondary_uncorrelated
