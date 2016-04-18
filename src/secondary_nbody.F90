module secondary_nbody

  use angleenergy_header, only: AngleEnergy
  use constants,          only: ONE, TWO, PI
  use math,               only: maxwell_spectrum
  use random_lcg,         only: prn

!===============================================================================
! NBODYPHASESPACE gives the energy distribution for particles emitted from
! neutron and charged-particle reactions. This corresponds to ACE law 66 and
! ENDF File 6, LAW=6.
!===============================================================================

  type, extends(AngleEnergy) :: NBodyPhaseSpace
    integer :: n_bodies
    real(8) :: mass_ratio
    real(8) :: A
    real(8) :: Q
  contains
    procedure :: sample => nbody_sample
  end type NBodyPhaseSpace

contains

  subroutine nbody_sample(this, E_in, E_out, mu)
    class(NBodyPhaseSpace), intent(in) :: this
    real(8), intent(in)  :: E_in  ! incoming energy
    real(8), intent(out) :: E_out ! sampled outgoing energy
    real(8), intent(out) :: mu    ! sampled outgoing energy

    real(8) :: Ap      ! total mass of particles in neutron masses
    real(8) :: E_max   ! maximum possible COM energy
    real(8) :: x, y, v
    real(8) :: r1, r2, r3, r4, r5, r6

    ! By definition, the distribution of the angle is isotropic for an N-body
    ! phase space distribution
    mu = TWO*prn() - ONE

    ! Determine E_max parameter
    Ap = this%mass_ratio
    E_max = (Ap - ONE)/Ap * (this%A/(this%A + ONE)*E_in + this%Q)

    ! x is essentially a Maxwellian distribution
    x = maxwell_spectrum(ONE)

    select case (this%n_bodies)
    case (3)
      y = maxwell_spectrum(ONE)
    case (4)
      r1 = prn()
      r2 = prn()
      r3 = prn()
      y = -log(r1*r2*r3)
    case (5)
      r1 = prn()
      r2 = prn()
      r3 = prn()
      r4 = prn()
      r5 = prn()
      r6 = prn()
      y = -log(r1*r2*r3*r4) - log(r5) * cos(PI/TWO*r6)**2
    end select

    ! Now determine v and E_out
    v = x/(x+y)
    E_out = E_max * v
  end subroutine nbody_sample

end module secondary_nbody
