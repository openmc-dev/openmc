module secondary_uncorrelated

  use angle_distribution, only: AngleDistribution
  use angleenergy_header, only: AngleEnergy
  use constants, only: ONE, TWO
  use energy_distribution, only: EnergyDistribution
  use random_lcg, only: prn

!===============================================================================
! UNCORRELATEDANGLEENERGY represents an uncorrelated angle-energy
! distribution. This corresponds to when an energy distribution is given in ENDF
! File 5/6 and an angular distribution is given in ENDF File 4.
!===============================================================================

  type, extends(AngleEnergy) :: UncorrelatedAngleEnergy
    logical :: fission = .false.
    type(AngleDistribution) :: angle
    class(EnergyDistribution), allocatable :: energy
  contains
    procedure :: sample => uncorrelated_sample
  end type UncorrelatedAngleEnergy

contains

  subroutine uncorrelated_sample(this, E_in, E_out, mu)
    class(UncorrelatedAngleEnergy), intent(in) :: this
    real(8), intent(in)  :: E_in  ! incoming energy
    real(8), intent(out) :: E_out ! sampled outgoing energy
    real(8), intent(out) :: mu    ! sampled scattering cosine

    ! Sample cosine of scattering angle
    if (this%fission) then
      ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< REMOVE THIS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      ! For fission, the angle is not used, so just assign a dummy value
      mu = ONE
      ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< REMOVE THIS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    elseif (allocated(this%angle%energy)) then
      mu = this%angle%sample(E_in)
    else
      ! no angle distribution given => assume isotropic for all energies
      mu = TWO*prn() - ONE
    end if

    ! Sample outgoing energy
    E_out = this%energy%sample(E_in)
  end subroutine uncorrelated_sample

end module secondary_uncorrelated
