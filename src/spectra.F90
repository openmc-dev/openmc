module spectra

use constants,  only: ZERO, ONE, TWO, PI
use random_lcg, only: prn

implicit none

contains

!===============================================================================
! MAXWELL_SPECTRUM samples an energy from the Maxwell fission distribution based
! on a direct sampling scheme. The probability distribution function for a
! Maxwellian is given as p(x) = 2/(T*sqrt(pi))*sqrt(x/T)*exp(-x/T). This PDF can
! be sampled using rule C64 in the Monte Carlo Sampler LA-9721-MS.
!===============================================================================

  function maxwell_spectrum(T) result(E_out)

    real(8), intent(in)  :: T     ! tabulated function of incoming E
    real(8)              :: E_out ! sampled energy

    real(8) :: r1, r2, r3  ! random numbers
    real(8) :: c           ! cosine of pi/2*r3

    r1 = prn()
    r2 = prn()
    r3 = prn()

    ! determine cosine of pi/2*r
    c = cos(PI/TWO*r3)

    ! determine outgoing energy
    E_out = -T*(log(r1) + log(r2)*c*c)

  end function maxwell_spectrum

!===============================================================================
! WATT_SPECTRUM samples the outgoing energy from a Watt energy-dependent fission
! spectrum. Although fitted parameters exist for many nuclides, generally the
! continuous tabular distributions (LAW 4) should be used in lieu of the Watt
! spectrum. This direct sampling scheme is an unpublished scheme based on the
! original Watt spectrum derivation (See F. Brown's MC lectures).
!===============================================================================

  function watt_spectrum(a, b) result(E_out)

    real(8), intent(in) :: a     ! Watt parameter a
    real(8), intent(in) :: b     ! Watt parameter b
    real(8)             :: E_out ! energy of emitted neutron

    real(8) :: w ! sampled from Maxwellian

    w     = maxwell_spectrum(a)
    E_out = w + a*a*b/4.0_8 + (TWO*prn() - ONE)*sqrt(a*a*b*w)

  end function watt_spectrum

!===============================================================================
! ROTATE_ANGLE rotates direction cosines through a polar angle whose cosine is
! mu and through an azimuthal angle sampled uniformly. Note that this is done
! with direct sampling rather than rejection as is done in MCNP and SERPENT.
!===============================================================================

  function rotate_angle(uvw0, mu, phi) result(uvw)
    real(8), intent(in) :: uvw0(3) ! directional cosine
    real(8), intent(in) :: mu      ! cosine of angle in lab or CM
    real(8), optional   :: phi     ! azimuthal angle
    real(8)             :: uvw(3)  ! rotated directional cosine

    real(8) :: phi_   ! azimuthal angle
    real(8) :: sinphi ! sine of azimuthal angle
    real(8) :: cosphi ! cosine of azimuthal angle
    real(8) :: a      ! sqrt(1 - mu^2)
    real(8) :: b      ! sqrt(1 - w^2)
    real(8) :: u0     ! original cosine in x direction
    real(8) :: v0     ! original cosine in y direction
    real(8) :: w0     ! original cosine in z direction

    ! Copy original directional cosines
    u0 = uvw0(1)
    v0 = uvw0(2)
    w0 = uvw0(3)

    ! Sample azimuthal angle in [0,2pi) if none provided
    if (present(phi)) then
      phi_ = phi
    else
      phi_ = TWO * PI * prn()
    end if

    ! Precompute factors to save flops
    sinphi = sin(phi_)
    cosphi = cos(phi_)
    a = sqrt(max(ZERO, ONE - mu*mu))
    b = sqrt(max(ZERO, ONE - w0*w0))

    ! Need to treat special case where sqrt(1 - w**2) is close to zero by
    ! expanding about the v component rather than the w component
    if (b > 1e-10) then
      uvw(1) = mu*u0 + a*(u0*w0*cosphi - v0*sinphi)/b
      uvw(2) = mu*v0 + a*(v0*w0*cosphi + u0*sinphi)/b
      uvw(3) = mu*w0 - a*b*cosphi
    else
      b = sqrt(ONE - v0*v0)
      uvw(1) = mu*u0 + a*(u0*v0*cosphi + w0*sinphi)/b
      uvw(2) = mu*v0 - a*b*cosphi
      uvw(3) = mu*w0 + a*(v0*w0*cosphi - u0*sinphi)/b
    end if

  end function rotate_angle

end module spectra