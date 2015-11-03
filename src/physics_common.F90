module physics_common

  use constants
  use global,                 only: weight_cutoff, weight_survive
  use particle_header,        only: Particle
  use random_lcg,             only: prn

  implicit none

contains

!===============================================================================
! RUSSIAN_ROULETTE
!===============================================================================

  subroutine russian_roulette(p)

    type(Particle), intent(inout) :: p

    if (p % wgt < weight_cutoff) then
      if (prn() < p % wgt / weight_survive) then
        p % wgt = weight_survive
        p % last_wgt = p % wgt
      else
        p % wgt = ZERO
        p % last_wgt = ZERO
        p % alive = .false.
      end if
    end if

  end subroutine russian_roulette

!===============================================================================
! ROTATE_ANGLE rotates direction cosines through a polar angle whose cosine is
! mu and through an azimuthal angle sampled uniformly. Note that this is done
! with direct sampling rather than rejection as is done in MCNP and SERPENT.
!===============================================================================

  function rotate_angle(uvw0, mu) result(uvw)

    real(8), intent(in) :: uvw0(3) ! directional cosine
    real(8), intent(in) :: mu      ! cosine of angle in lab or CM
    real(8)             :: uvw(3)  ! rotated directional cosine

    real(8) :: phi    ! azimuthal angle
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

    ! Sample azimuthal angle in [0,2pi)
    phi = TWO * PI * prn()

    ! Precompute factors to save flops
    sinphi = sin(phi)
    cosphi = cos(phi)
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
end module physics_common