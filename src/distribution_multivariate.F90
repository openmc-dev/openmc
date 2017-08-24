module distribution_multivariate

  use constants,               only: ONE, TWO, PI
  use distribution_univariate
  use error,                   only: fatal_error
  use random_lcg,              only: prn
  use math,                    only: rotate_angle
  use xml_interface

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
    procedure(spatial_distribution_from_xml_), deferred :: from_xml
    procedure(spatial_distribution_sample_), deferred :: sample
  end type SpatialDistribution

  abstract interface
    subroutine spatial_distribution_from_xml_(this, node)
      import SpatialDistribution, XMLNode
      class(SpatialDistribution), intent(inout) :: this
      type(XMLNode), intent(in) :: node
    end subroutine spatial_distribution_from_xml_

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
    procedure :: from_xml => cartesian_independent_from_xml
    procedure :: sample => cartesian_independent_sample
  end type CartesianIndependent

  type, extends(SpatialDistribution) :: SpatialBox
    real(8) :: lower_left(3)
    real(8) :: upper_right(3)
    logical :: only_fissionable = .false.
  contains
    procedure :: from_xml => spatial_box_from_xml
    procedure :: sample => spatial_box_sample
  end type SpatialBox

  type, extends(SpatialDistribution) :: SpatialPoint
    real(8) :: xyz(3)
  contains
    procedure :: from_xml => spatial_point_from_xml
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

  subroutine cartesian_independent_from_xml(this, node)
    class(CartesianIndependent), intent(inout) :: this
    type(XMLNode), intent(in) :: node

    type(XMLNode) :: node_dist

    ! Read distribution for x coordinate
    if (check_for_node(node, "x")) then
      node_dist = node % child("x")
      call distribution_from_xml(this % x, node_dist)
    else
      allocate(Discrete :: this % x)
      select type (dist => this % x)
      type is (Discrete)
        allocate(dist % x(1), dist % p(1))
        dist % x(1) = ZERO
        dist % p(1) = ONE
      end select
    end if

    ! Read distribution for y coordinate
    if (check_for_node(node, "y")) then
      node_dist = node % child("y")
      call distribution_from_xml(this % y, node_dist)
    else
      allocate(Discrete :: this % y)
      select type (dist => this % y)
      type is (Discrete)
        allocate(dist % x(1), dist % p(1))
        dist % x(1) = ZERO
        dist % p(1) = ONE
      end select
    end if

    if (check_for_node(node, "z")) then
      node_dist = node % child("z")
      call distribution_from_xml(this % z, node_dist)
    else
      allocate(Discrete :: this % z)
      select type (dist => this % z)
      type is (Discrete)
        allocate(dist % x(1), dist % p(1))
        dist % x(1) = ZERO
        dist % p(1) = ONE
      end select
    end if

  end subroutine cartesian_independent_from_xml

  function cartesian_independent_sample(this) result(xyz)
    class(CartesianIndependent), intent(in) :: this
    real(8) :: xyz(3)

    xyz(1) = this % x % sample()
    xyz(2) = this % y % sample()
    xyz(3) = this % z % sample()
  end function cartesian_independent_sample

  subroutine spatial_box_from_xml(this, node)
    class(SpatialBox), intent(inout) :: this
    type(XMLNode), intent(in) :: node

    real(8), allocatable :: temp_real(:)

    ! Make sure correct number of parameters are given
    if (node_word_count(node, "parameters") /= 6) then
      call fatal_error('Box/fission spatial source must have &
           &six parameters specified.')
    end if

    ! Read lower-right/upper-left coordinates
    allocate(temp_real(6))
    call get_node_array(node, "parameters", temp_real)
    this % lower_left(:) = temp_real(1:3)
    this % upper_right(:) = temp_real(4:6)
    deallocate(temp_real)
  end subroutine spatial_box_from_xml

  function spatial_box_sample(this) result(xyz)
    class(SpatialBox), intent(in) :: this
    real(8) :: xyz(3)

    integer :: i
    real(8) :: r(3)

    r = [ (prn(), i = 1,3) ]
    xyz(:) = this % lower_left + r*(this % upper_right - this % lower_left)
  end function spatial_box_sample

  subroutine spatial_point_from_xml(this, node)
    class(SpatialPoint), intent(inout) :: this
    type(XMLNode), intent(in) :: node

    ! Make sure correct number of parameters are given
    if (node_word_count(node, "parameters") /= 3) then
      call fatal_error('Point spatial source must have &
           &three parameters specified.')
    end if

    ! Read location of point source
    call get_node_array(node, "parameters", this % xyz)
  end subroutine spatial_point_from_xml

  function spatial_point_sample(this) result(xyz)
    class(SpatialPoint), intent(in) :: this
    real(8) :: xyz(3)

    xyz(:) = this % xyz
  end function spatial_point_sample

end module distribution_multivariate
