module distribution_multivariate

  use constants, only: ONE, TWO, PI
  use distribution_univariate, only: Distribution
  use math, only: rotate_angle
  use random_lcg, only: prn

!===============================================================================
! UNITSPHEREDISTRIBUTION type defines a probability density function for points
! on the unit sphere. Extensions of this type are used to sample angular
! distributions for starting soures
!===============================================================================

  type, abstract :: UnitSphereDistribution
    real(8) :: reference_uvw(3)
  contains
    procedure(iSample), deferred :: sample
  end type UnitSphereDistribution

  abstract interface
    function iSample(this) result(uvw)
      import UnitSphereDistribution
      class(UnitSphereDistribution), intent(in) :: this
      real(8) :: uvw(3)
    end function iSample
  end interface

!===============================================================================
! Derived classes of UnitSphereDistribution
!===============================================================================

  ! Explicit distribution of polar and azimuthal angles
  type, extends(UnitSphereDistribution) :: PolarAzimuthal
    class(Distribution), allocatable :: mu
    class(Distribution), allocatable :: phi
  contains
    procedure :: sample => polar_azimuthal_sample
  end type PolarAzimuthal

  ! Uniform distribution on the unit sphere
  type, extends(UnitSphereDistribution) :: Isotropic
  contains
    procedure :: sample => isotropic_sample
  end type Isotropic

  ! Monodirectional distribution
  type, extends(UnitSphereDistribution) :: Monodirectional
  contains
    procedure :: sample => monodirectional_sample
  end type Monodirectional

contains

  function polar_azimuthal_sample(this) result(uvw)
    class(PolarAzimuthal), intent(in) :: this
    real(8) :: uvw(3)

    real(8) :: mu     ! cosine of polar angle
    real(8) :: phi    ! azimuthal angle

    ! Sample cosine of polar angle
    mu = this%mu%sample()
    if (mu == ONE) then
      uvw(:) = this%reference_uvw
    else
      ! Sample azimuthal angle
      phi = this%phi%sample()
      uvw(:) = rotate_angle(this%reference_uvw, mu, phi)
    end if
  end function polar_azimuthal_sample

  function isotropic_sample(this) result(uvw)
    class(Isotropic), intent(in) :: this
    real(8) :: uvw(3)

    real(8) :: phi
    real(8) :: mu

    phi = TWO*PI*prn()
    mu = TWO*prn() - ONE
    uvw(1) = mu
    uvw(2) = sqrt(ONE - mu*mu) * cos(phi)
    uvw(3) = sqrt(ONE - mu*mu) * sin(phi)
  end function isotropic_sample

  function monodirectional_sample(this) result(uvw)
    class(Monodirectional), intent(in) :: this
    real(8) :: uvw(3)

    uvw(:) = this%reference_uvw
  end function monodirectional_sample

end module distribution_multivariate
