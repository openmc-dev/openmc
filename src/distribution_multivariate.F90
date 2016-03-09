module distribution_multivariate

  use constants,               only: ONE, TWO, PI
  use distribution_univariate, only: Distribution
  use random_lcg,              only: prn
  use math,                    only: rotate_angle

  implicit none

!===============================================================================
! UNITSPHEREDISTRIBUTION type defines a probability density function for points
! on the unit sphere. Extensions of this type are used to sample angular
! distributions for starting sources
!===============================================================================

  type, abstract :: UnitSphereDistribution
    real(8) :: reference_uvw(3)
  contains
    procedure(unitsphere_distribution_sample_), deferred :: sample
  end type UnitSphereDistribution

  abstract interface
    function unitsphere_distribution_sample_(this) result(uvw)
      import UnitSphereDistribution
      class(UnitSphereDistribution), intent(in) :: this
      real(8) :: uvw(3)
    end function unitsphere_distribution_sample_
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

!===============================================================================
! SPATIALDISTRIBUTION type defines a probability density function for arbitrary
! points in Euclidean space.
!===============================================================================

  type, abstract :: SpatialDistribution
  contains
    procedure(spatial_distribution_sample_), deferred :: sample
  end type SpatialDistribution

  abstract interface
    function spatial_distribution_sample_(this) result(xyz)
      import SpatialDistribution
      class(SpatialDistribution), intent(in) :: this
      real(8) :: xyz(3)
    end function spatial_distribution_sample_
  end interface

  type, extends(SpatialDistribution) :: CartesianIndependent
    class(Distribution), allocatable :: x
    class(Distribution), allocatable :: y
    class(Distribution), allocatable :: z
  contains
    procedure :: sample => cartesian_independent_sample
  end type CartesianIndependent

  type, extends(SpatialDistribution) :: SpatialBox
    real(8) :: lower_left(3)
    real(8) :: upper_right(3)
    logical :: only_fissionable = .false.
  contains
    procedure :: sample => spatial_box_sample
  end type SpatialBox

  type, extends(SpatialDistribution) :: SpatialPoint
    real(8) :: xyz(3)
  contains
    procedure :: sample => spatial_point_sample
  end type SpatialPoint

contains

  function polar_azimuthal_sample(this) result(uvw)
    class(PolarAzimuthal), intent(in) :: this
    real(8) :: uvw(3)

    real(8) :: mu     ! cosine of polar angle
    real(8) :: phi    ! azimuthal angle

    ! Sample cosine of polar angle
    mu = this % mu % sample()
    if (mu == ONE) then
      uvw(:) = this % reference_uvw
    else
      ! Sample azimuthal angle
      phi = this % phi % sample()
      uvw(:) = rotate_angle(this % reference_uvw, mu, phi)
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

    uvw(:) = this % reference_uvw
  end function monodirectional_sample

  function cartesian_independent_sample(this) result(xyz)
    class(CartesianIndependent), intent(in) :: this
    real(8) :: xyz(3)

    xyz(1) = this % x % sample()
    xyz(2) = this % y % sample()
    xyz(3) = this % z % sample()
  end function cartesian_independent_sample

  function spatial_box_sample(this) result(xyz)
    class(SpatialBox), intent(in) :: this
    real(8) :: xyz(3)

    integer :: i
    real(8) :: r(3)

    r = [ (prn(), i = 1,3) ]
    xyz(:) = this % lower_left + r*(this % upper_right - this % lower_left)
  end function spatial_box_sample

  function spatial_point_sample(this) result(xyz)
    class(SpatialPoint), intent(in) :: this
    real(8) :: xyz(3)

    xyz(:) = this % xyz
  end function spatial_point_sample

end module distribution_multivariate
