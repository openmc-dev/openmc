module surface_header

  use constants, only: ONE, TWO, ZERO, INFINITY, FP_COINCIDENT

  implicit none

!===============================================================================
! SURFACE type defines a first- or second-order surface that can be used to
! construct closed volumes (cells)
!===============================================================================

  type, abstract :: Surface
    integer :: id                     ! Unique ID
    integer, allocatable :: &
         neighbor_pos(:), &           ! List of cells on positive side
         neighbor_neg(:)              ! List of cells on negative side
    integer :: bc                     ! Boundary condition
    character(len=52) :: name = ""    ! User-defined name
  contains
    procedure(iEvaluate), deferred :: evaluate
    procedure(iDistance), deferred :: distance
    procedure(iReflect),  deferred :: reflect
    procedure(iNormal),   deferred :: normal
  end type Surface

  type :: SurfaceContainer
    class(Surface), allocatable :: obj
  end type SurfaceContainer

  type, extends(Surface) :: SurfaceXPlane
    ! x = x0
    real(8) :: x0
  contains
    procedure :: evaluate => x_plane_evaluate
    procedure :: reflect  => x_plane_reflect
    procedure :: distance => x_plane_distance
    procedure :: normal => x_plane_normal
  end type SurfaceXPlane

  type, extends(Surface) :: SurfaceYPlane
    ! y = y0
    real(8) :: y0
  contains
    procedure :: evaluate => y_plane_evaluate
    procedure :: reflect  => y_plane_reflect
    procedure :: distance => y_plane_distance
    procedure :: normal => y_plane_normal
  end type SurfaceYPlane

  type, extends(Surface) :: SurfaceZPlane
    ! z = z0
    real(8) :: z0
  contains
    procedure :: evaluate => z_plane_evaluate
    procedure :: reflect  => z_plane_reflect
    procedure :: distance => z_plane_distance
    procedure :: normal => z_plane_normal
  end type SurfaceZPlane

  type, extends(Surface) :: SurfacePlane
    ! Ax + By + Cz = D
    real(8) :: A
    real(8) :: B
    real(8) :: C
    real(8) :: D
  contains
    procedure :: evaluate => plane_evaluate
    procedure :: reflect  => plane_reflect
    procedure :: distance => plane_distance
    procedure :: normal => plane_normal
  end type SurfacePlane

  type, extends(Surface) :: SurfaceXCylinder
    ! (y - y0)^2 + (z - z0)^2 = R^2
    real(8) :: y0
    real(8) :: z0
    real(8) :: r
  contains
    procedure :: evaluate => x_cylinder_evaluate
    procedure :: reflect  => x_cylinder_reflect
    procedure :: distance => x_cylinder_distance
    procedure :: normal => x_cylinder_normal
  end type SurfaceXCylinder

  type, extends(Surface) :: SurfaceYCylinder
    ! (x - x0)^2 + (z - z0)^2 = R^2
    real(8) :: x0
    real(8) :: z0
    real(8) :: r
  contains
    procedure :: evaluate => y_cylinder_evaluate
    procedure :: reflect  => y_cylinder_reflect
    procedure :: distance => y_cylinder_distance
    procedure :: normal => y_cylinder_normal
  end type SurfaceYCylinder

  type, extends(Surface) :: SurfaceZCylinder
    ! (x - x0)^2 + (y - y0)^2 = R^2
    real(8) :: x0
    real(8) :: y0
    real(8) :: r
  contains
    procedure :: evaluate => z_cylinder_evaluate
    procedure :: reflect  => z_cylinder_reflect
    procedure :: distance => z_cylinder_distance
    procedure :: normal => z_cylinder_normal
  end type SurfaceZCylinder

  type, extends(Surface) :: SurfaceSphere
    ! (x - x0)^2 + (y - y0)^2 + (z - z0)^2 = R^2
    real(8) :: x0
    real(8) :: y0
    real(8) :: z0
    real(8) :: r
  contains
    procedure :: evaluate => sphere_evaluate
    procedure :: reflect  => sphere_reflect
    procedure :: distance => sphere_distance
    procedure :: normal => sphere_normal
  end type SurfaceSphere

  type, extends(Surface) :: SurfaceXCone
    ! (y - y0)^2 + (z - z0)^2 = R^2*(x - x0)^2
    real(8) :: x0
    real(8) :: y0
    real(8) :: z0
    real(8) :: r2
  contains
    procedure :: evaluate => x_cone_evaluate
    procedure :: reflect  => x_cone_reflect
    procedure :: distance => x_cone_distance
    procedure :: normal => x_cone_normal
  end type SurfaceXCone

  type, extends(Surface) :: SurfaceYCone
    ! (x - x0)^2 + (z - z0)^2 = R^2*(y - y0)^2
    real(8) :: x0
    real(8) :: y0
    real(8) :: z0
    real(8) :: r2
  contains
    procedure :: evaluate => y_cone_evaluate
    procedure :: reflect  => y_cone_reflect
    procedure :: distance => y_cone_distance
    procedure :: normal => y_cone_normal
  end type SurfaceYCone

  type, extends(Surface) :: SurfaceZCone
    ! (x - x0)^2 + (y - y0)^2 = R^2*(z - z0)^2
    real(8) :: x0
    real(8) :: y0
    real(8) :: z0
    real(8) :: r2
  contains
    procedure :: evaluate => z_cone_evaluate
    procedure :: reflect  => z_cone_reflect
    procedure :: distance => z_cone_distance
    procedure :: normal => z_cone_normal
  end type SurfaceZCone

  abstract interface
    pure function iEvaluate(this, xyz) result(f)
      import Surface
      class(Surface), intent(in) :: this
      real(8), intent(in) :: xyz(3)
      real(8) :: f
    end function iEvaluate

    pure function iDistance(this, xyz, uvw, coincident) result(d)
      import Surface
      class(Surface), intent(in) :: this
      real(8), intent(in) :: xyz(3)
      real(8), intent(in) :: uvw(3)
      logical, intent(in) :: coincident
      real(8) :: d
    end function iDistance

    subroutine iReflect(this, xyz, uvw)
      import Surface
      class(Surface), intent(in) :: this
      real(8), intent(in) :: xyz(3)
      real(8), intent(inout) :: uvw(3)
    end subroutine iReflect

    pure function iNormal(this, xyz) result(uvw)
      import Surface
      class(Surface), intent(in) :: this
      real(8), intent(in) :: xyz(3)
      real(8) :: uvw(3)
    end function iNormal
  end interface

contains

!===============================================================================
! SurfaceXPlane Implementation
!===============================================================================

  pure function x_plane_evaluate(this, xyz) result(f)
    class(SurfaceXPlane), intent(in) :: this
    real(8), intent(in) :: xyz(3)
    real(8) :: f

    f = xyz(1) - this%x0
  end function x_plane_evaluate

  pure function x_plane_distance(this, xyz, uvw, coincident) result(d)
    class(SurfaceXPlane), intent(in) :: this
    real(8), intent(in) :: xyz(3)
    real(8), intent(in) :: uvw(3)
    logical, intent(in) :: coincident
    real(8) :: d

    real(8) :: f

    f = this%x0 - xyz(1)
    if (coincident .or. abs(f) < FP_COINCIDENT .or. uvw(1) == ZERO) then
      d = INFINITY
    else
      d = f/uvw(1)
      if (d < ZERO) d = INFINITY
    end if
  end function x_plane_distance

  subroutine x_plane_reflect(this, xyz, uvw)
    class(SurfaceXPlane), intent(in) :: this
    real(8), intent(in)    :: xyz(3)
    real(8), intent(inout) :: uvw(3)

    uvw(1) = -uvw(1)
  end subroutine x_plane_reflect

  pure function x_plane_normal(this, xyz) result(uvw)
    class(SurfaceXPlane), intent(in) :: this
    real(8), intent(in) :: xyz(3)
    real(8) :: uvw(3)

    uvw(:) = [1, 0, 0]
  end function x_plane_normal

!===============================================================================
! SurfaceYPlane Implementation
!===============================================================================

  pure function y_plane_evaluate(this, xyz) result(f)
    class(SurfaceYPlane), intent(in) :: this
    real(8), intent(in) :: xyz(3)
    real(8) :: f

    f = xyz(2) - this%y0
  end function y_plane_evaluate

  pure function y_plane_distance(this, xyz, uvw, coincident) result(d)
    class(SurfaceYPlane), intent(in) :: this
    real(8), intent(in) :: xyz(3)
    real(8), intent(in) :: uvw(3)
    logical, intent(in) :: coincident
    real(8) :: d

    real(8) :: f

    f = this%y0 - xyz(2)
    if (coincident .or. abs(f) < FP_COINCIDENT .or. uvw(2) == ZERO) then
      d = INFINITY
    else
      d = f/uvw(2)
      if (d < ZERO) d = INFINITY
    end if
  end function y_plane_distance

  subroutine y_plane_reflect(this, xyz, uvw)
    class(SurfaceYPlane), intent(in)    :: this
    real(8), intent(in)    :: xyz(3)
    real(8), intent(inout) :: uvw(3)

    uvw(2) = -uvw(2)
  end subroutine y_plane_reflect

  pure function y_plane_normal(this, xyz) result(uvw)
    class(SurfaceYPlane), intent(in) :: this
    real(8), intent(in) :: xyz(3)
    real(8) :: uvw(3)

    uvw(:) = [0, 1, 0]
  end function y_plane_normal

!===============================================================================
! SurfaceZPlane Implementation
!===============================================================================

  pure function z_plane_evaluate(this, xyz) result(f)

    class(SurfaceZPlane), intent(in) :: this
    real(8),        intent(in) :: xyz(3)
    real(8)             :: f

    f = xyz(3) - this%z0

  end function z_plane_evaluate

  pure function z_plane_distance(this, xyz, uvw, coincident) result(d)
    class(SurfaceZPlane), intent(in) :: this
    real(8), intent(in) :: xyz(3)
    real(8), intent(in) :: uvw(3)
    logical, intent(in) :: coincident
    real(8) :: d

    real(8) :: f

    f = this%z0 - xyz(3)
    if (coincident .or. abs(f) < FP_COINCIDENT .or. uvw(3) == ZERO) then
      d = INFINITY
    else
      d = f/uvw(3)
      if (d < ZERO) d = INFINITY
    end if
  end function z_plane_distance

  subroutine z_plane_reflect(this, xyz, uvw)
    class(SurfaceZPlane), intent(in) :: this
    real(8), intent(in)    :: xyz(3)
    real(8), intent(inout) :: uvw(3)

    uvw(3) = -uvw(3)
  end subroutine z_plane_reflect

  pure function z_plane_normal(this, xyz) result(uvw)
    class(SurfaceZPlane), intent(in) :: this
    real(8), intent(in) :: xyz(3)
    real(8) :: uvw(3)

    uvw(:) = [0, 0, 1]
  end function z_plane_normal

!===============================================================================
! SurfacePlane Implementation
!===============================================================================

  pure function plane_evaluate(this, xyz) result(f)
    class(SurfacePlane), intent(in) :: this
    real(8), intent(in) :: xyz(3)
    real(8) :: f

    f = this%A*xyz(1) + this%B*xyz(2) + this%C*xyz(3) - this%D
  end function plane_evaluate

  pure function plane_distance(this, xyz, uvw, coincident) result(d)
    class(SurfacePlane), intent(in) :: this
    real(8), intent(in) :: xyz(3)
    real(8), intent(in) :: uvw(3)
    logical, intent(in) :: coincident
    real(8) :: d

    real(8) :: f
    real(8) :: tmp

    f = this%A*xyz(1) + this%B*xyz(2) + this%C*xyz(3) - this%D
    tmp = this%A*uvw(1) + this%B*uvw(2) + this%C*uvw(3)
    if (coincident .or. abs(f) < FP_COINCIDENT .or. tmp == ZERO) then
      d = INFINITY
    else
      d = -f/tmp
      if (d < ZERO) d = INFINITY
    end if
  end function plane_distance

  subroutine plane_reflect(this, xyz, uvw)
    class(SurfacePlane), intent(in)    :: this
    real(8),        intent(in)    :: xyz(3)
    real(8),        intent(inout) :: uvw(3)

    real(8) :: n(3)

    ! Construct normal vector
    n(:) = [this%A, this%B, this%C]

    ! Reflect direction according to normal
    uvw = uvw - TWO*dot_product(n, uvw)/dot_product(n, n) * n
  end subroutine plane_reflect

  pure function plane_normal(this, xyz) result(uvw)
    class(SurfacePlane), intent(in) :: this
    real(8), intent(in) :: xyz(3)
    real(8) :: uvw(3)

    uvw(:) = [this%A, this%B, this%C]
  end function plane_normal

!===============================================================================
! SurfaceXCylinder Implementation
!===============================================================================

  pure function x_cylinder_evaluate(this, xyz) result(f)
    class(SurfaceXCylinder), intent(in) :: this
    real(8), intent(in) :: xyz(3)
    real(8) :: f

    real(8) :: y, z

    y = xyz(2) - this%y0
    z = xyz(3) - this%z0
    f = y*y + z*z - this%r*this%r
  end function x_cylinder_evaluate

  pure function x_cylinder_distance(this, xyz, uvw, coincident) result(d)
    class(SurfaceXCylinder), intent(in) :: this
    real(8), intent(in) :: xyz(3)
    real(8), intent(in) :: uvw(3)
    logical, intent(in) :: coincident
    real(8) :: d

    real(8) :: y, z, k, a, c, quad

    a = ONE - uvw(1)*uvw(1)  ! v^2 + w^2
    if (a == ZERO) then
      d = INFINITY
    else
      y = xyz(2) - this%y0
      z = xyz(3) - this%z0
      k = y*uvw(2) + z*uvw(3)
      c = y*y + z*z - this%r*this%r
      quad = k*k - a*c

      if (quad < ZERO) then
        ! no intersection with cylinder

        d = INFINITY

      elseif (coincident .or. abs(c) < FP_COINCIDENT) then
        ! particle is on the cylinder, thus one distance is positive/negative
        ! and the other is zero. The sign of k determines if we are facing in or
        ! out

        if (k >= ZERO) then
          d = INFINITY
        else
          d = (-k + sqrt(quad))/a
        end if

      elseif (c < ZERO) then
        ! particle is inside the cylinder, thus one distance must be negative
        ! and one must be positive. The positive distance will be the one with
        ! negative sign on sqrt(quad)

        d = (-k + sqrt(quad))/a

      else
        ! particle is outside the cylinder, thus both distances are either
        ! positive or negative. If positive, the smaller distance is the one
        ! with positive sign on sqrt(quad)

        d = (-k - sqrt(quad))/a
        if (d < ZERO) d = INFINITY

      end if
    end if
  end function x_cylinder_distance

  subroutine x_cylinder_reflect(this, xyz, uvw)
    class(SurfaceXCylinder), intent(in)    :: this
    real(8),        intent(in)    :: xyz(3)
    real(8),        intent(inout) :: uvw(3)

    real(8) :: y, z, dot_prod

    ! Find y-y0, z-z0 and dot product of direction and surface normal
    y = xyz(2) - this%y0
    z = xyz(3) - this%z0
    dot_prod = uvw(2)*y + uvw(3)*z

    ! Reflect direction according to normal
    uvw(2) = uvw(2) - TWO*dot_prod*y/(this%r*this%r)
    uvw(3) = uvw(3) - TWO*dot_prod*z/(this%r*this%r)
  end subroutine x_cylinder_reflect

  pure function x_cylinder_normal(this, xyz) result(uvw)
    class(SurfaceXCylinder), intent(in) :: this
    real(8), intent(in) :: xyz(3)
    real(8) :: uvw(3)

    uvw(1) = ZERO
    uvw(2) = TWO*(xyz(2) - this%y0)
    uvw(3) = TWO*(xyz(3) - this%z0)
  end function x_cylinder_normal

!===============================================================================
! SurfaceYCylinder Implementation
!===============================================================================

  pure function y_cylinder_evaluate(this, xyz) result(f)
    class(SurfaceYCylinder), intent(in) :: this
    real(8), intent(in) :: xyz(3)
    real(8) :: f

    real(8) :: x, z

    x = xyz(1) - this%x0
    z = xyz(3) - this%z0
    f = x*x + z*z - this%r*this%r
  end function y_cylinder_evaluate

  pure function y_cylinder_distance(this, xyz, uvw, coincident) result(d)
    class(SurfaceYCylinder), intent(in) :: this
    real(8), intent(in) :: xyz(3)
    real(8), intent(in) :: uvw(3)
    logical, intent(in) :: coincident
    real(8) :: d

    real(8) :: x, z, k, a, c, quad

    a = ONE - uvw(2)*uvw(2)  ! u^2 + w^2
    if (a == ZERO) then
      d = INFINITY
    else
      x = xyz(1) - this%x0
      z = xyz(3) - this%z0
      k = x*uvw(1) + z*uvw(3)
      c = x*x + z*z - this%r*this%r
      quad = k*k - a*c

      if (quad < ZERO) then
        ! no intersection with cylinder

        d = INFINITY

      elseif (coincident .or. abs(c) < FP_COINCIDENT) then
        ! particle is on the cylinder, thus one distance is positive/negative
        ! and the other is zero. The sign of k determines if we are facing in or
        ! out

        if (k >= ZERO) then
          d = INFINITY
        else
          d = (-k + sqrt(quad))/a
        end if

      elseif (c < ZERO) then
        ! particle is inside the cylinder, thus one distance must be negative
        ! and one must be positive. The positive distance will be the one with
        ! negative sign on sqrt(quad)

        d = (-k + sqrt(quad))/a

      else
        ! particle is outside the cylinder, thus both distances are either
        ! positive or negative. If positive, the smaller distance is the one
        ! with positive sign on sqrt(quad)

        d = (-k - sqrt(quad))/a
        if (d < ZERO) d = INFINITY

      end if
    end if
  end function y_cylinder_distance

  subroutine y_cylinder_reflect(this, xyz, uvw)
    class(SurfaceYCylinder), intent(in) :: this
    real(8), intent(in)    :: xyz(3)
    real(8), intent(inout) :: uvw(3)

    real(8) :: x, z, dot_prod

    ! Find x-x0, z-z0 and dot product of direction and surface normal
    x = xyz(1) - this%x0
    z = xyz(3) - this%z0
    dot_prod = uvw(1)*x + uvw(3)*z

    ! Reflect direction according to normal
    uvw(1) = uvw(1) - TWO*dot_prod*x/(this%r*this%r)
    uvw(3) = uvw(3) - TWO*dot_prod*z/(this%r*this%r)
  end subroutine y_cylinder_reflect

  pure function y_cylinder_normal(this, xyz) result(uvw)
    class(SurfaceYCylinder), intent(in) :: this
    real(8), intent(in) :: xyz(3)
    real(8) :: uvw(3)

    uvw(1) = TWO*(xyz(1) - this%x0)
    uvw(2) = ZERO
    uvw(3) = TWO*(xyz(3) - this%z0)
  end function y_cylinder_normal

!===============================================================================
! SurfaceZCylinder Implementation
!===============================================================================

  pure function z_cylinder_evaluate(this, xyz) result(f)
    class(SurfaceZCylinder), intent(in) :: this
    real(8), intent(in) :: xyz(3)
    real(8) :: f

    real(8) :: x, y

    x = xyz(1) - this%x0
    y = xyz(2) - this%y0
    f = x*x + y*y - this%r*this%r
  end function z_cylinder_evaluate

  pure function z_cylinder_distance(this, xyz, uvw, coincident) result(d)
    class(SurfaceZCylinder), intent(in) :: this
    real(8), intent(in) :: xyz(3)
    real(8), intent(in) :: uvw(3)
    logical, intent(in) :: coincident
    real(8) :: d

    real(8) :: x, y, k, a, c, quad

    a = ONE - uvw(3)*uvw(3)  ! u^2 + v^2
    if (a == ZERO) then
      d = INFINITY
    else
      x = xyz(1) - this%x0
      y = xyz(2) - this%y0
      k = x*uvw(1) + y*uvw(2)
      c = x*x + y*y - this%r*this%r
      quad = k*k - a*c

      if (quad < ZERO) then
        ! no intersection with cylinder

        d = INFINITY

      elseif (coincident .or. abs(c) < FP_COINCIDENT) then
        ! particle is on the cylinder, thus one distance is positive/negative
        ! and the other is zero. The sign of k determines if we are facing in or
        ! out

        if (k >= ZERO) then
          d = INFINITY
        else
          d = (-k + sqrt(quad))/a
        end if

      elseif (c < ZERO) then
        ! particle is inside the cylinder, thus one distance must be negative
        ! and one must be positive. The positive distance will be the one with
        ! negative sign on sqrt(quad)

        d = (-k + sqrt(quad))/a

      else
        ! particle is outside the cylinder, thus both distances are either
        ! positive or negative. If positive, the smaller distance is the one
        ! with positive sign on sqrt(quad)

        d = (-k - sqrt(quad))/a
        if (d < ZERO) d = INFINITY

      end if
    end if
  end function z_cylinder_distance

  subroutine z_cylinder_reflect(this, xyz, uvw)
    class(SurfaceZCylinder), intent(in) :: this
    real(8), intent(in)    :: xyz(3)
    real(8), intent(inout) :: uvw(3)

    real(8) :: x, y, dot_prod

    ! Find x-x0, y-y0 and dot product of direction and surface normal
    x = xyz(1) - this%x0
    y = xyz(2) - this%y0
    dot_prod = uvw(1)*x + uvw(2)*y

    ! Reflect direction according to normal
    uvw(1) = uvw(1) - TWO*dot_prod*x/(this%r*this%r)
    uvw(2) = uvw(2) - TWO*dot_prod*y/(this%r*this%r)
  end subroutine z_cylinder_reflect

  pure function z_cylinder_normal(this, xyz) result(uvw)
    class(SurfaceZCylinder), intent(in) :: this
    real(8), intent(in) :: xyz(3)
    real(8) :: uvw(3)

    uvw(1) = TWO*(xyz(1) - this%x0)
    uvw(2) = TWO*(xyz(2) - this%y0)
    uvw(3) = ZERO
  end function z_cylinder_normal

!===============================================================================
! SphereImplementation
!===============================================================================

  pure function sphere_evaluate(this, xyz) result(f)
    class(SurfaceSphere), intent(in) :: this
    real(8), intent(in) :: xyz(3)
    real(8) :: f

    real(8) :: x, y, z

    x = xyz(1) - this%x0
    y = xyz(2) - this%y0
    z = xyz(3) - this%z0
    f = x*x + y*y + z*z - this%r*this%r
  end function sphere_evaluate

  pure function sphere_distance(this, xyz, uvw, coincident) result(d)
    class(SurfaceSphere), intent(in) :: this
    real(8), intent(in) :: xyz(3)
    real(8), intent(in) :: uvw(3)
    logical, intent(in) :: coincident
    real(8) :: d

    real(8) :: x, y, z, k, c, quad

    x = xyz(1) - this%x0
    y = xyz(2) - this%y0
    z = xyz(3) - this%z0
    k = x*uvw(1) + y*uvw(2) + z*uvw(3)
    c = x*x + y*y + z*z - this%r*this%r
    quad = k*k - c

    if (quad < ZERO) then
      ! no intersection with sphere

      d = INFINITY

    elseif (coincident .or. abs(c) < FP_COINCIDENT) then
      ! particle is on the sphere, thus one distance is positive/negative and
      ! the other is zero. The sign of k determines if we are facing in or out

      if (k >= ZERO) then
        d = INFINITY
      else
        d = -k + sqrt(quad)
      end if

    elseif (c < ZERO) then
      ! particle is inside the sphere, thus one distance must be negative and
      ! one must be positive. The positive distance will be the one with
      ! negative sign on sqrt(quad)

      d = -k + sqrt(quad)

    else
      ! particle is outside the sphere, thus both distances are either positive
      ! or negative. If positive, the smaller distance is the one with positive
      ! sign on sqrt(quad)

      d = -k - sqrt(quad)
      if (d < ZERO) d = INFINITY

    end if
  end function sphere_distance

  subroutine sphere_reflect(this, xyz, uvw)
    class(SurfaceSphere), intent(in)    :: this
    real(8),        intent(in)    :: xyz(3)
    real(8),        intent(inout) :: uvw(3)

    real(8) :: x, y, z, dot_prod

    x = xyz(1) - this%x0
    y = xyz(2) - this%y0
    z = xyz(3) - this%z0
    dot_prod = uvw(1)*x + uvw(2)*y + uvw(3)*z

    uvw(1) = uvw(1) - TWO*dot_prod*x/(this%r*this%r)
    uvw(2) = uvw(2) - TWO*dot_prod*y/(this%r*this%r)
    uvw(3) = uvw(3) - TWO*dot_prod*z/(this%r*this%r)
  end subroutine sphere_reflect

  pure function sphere_normal(this, xyz) result(uvw)
    class(SurfaceSphere), intent(in) :: this
    real(8), intent(in) :: xyz(3)
    real(8) :: uvw(3)

    uvw(:) = TWO*(xyz - [this%x0, this%y0, this%z0])
  end function sphere_normal

!===============================================================================
! SurfaceXCone Implementation
!===============================================================================

  pure function x_cone_evaluate(this, xyz) result(f)
    class(SurfaceXCone), intent(in) :: this
    real(8), intent(in) :: xyz(3)
    real(8) :: f

    real(8) :: x, y, z

    x = xyz(1) - this%x0
    y = xyz(2) - this%y0
    z = xyz(3) - this%z0
    f = y*y + z*z - this%r2*x*x
  end function x_cone_evaluate

  pure function x_cone_distance(this, xyz, uvw, coincident) result(d)
    class(SurfaceXCone), intent(in) :: this
    real(8), intent(in) :: xyz(3)
    real(8), intent(in) :: uvw(3)
    logical, intent(in) :: coincident
    real(8) :: d

    real(8) :: x, y, z, k, a, b, c, quad

    x = xyz(1) - this%x0
    y = xyz(2) - this%y0
    z = xyz(3) - this%z0
    a = uvw(2)*uvw(2) + uvw(3)*uvw(3) - this%r2*uvw(1)*uvw(1)
    k = y*uvw(2) + z*uvw(3) - this%r2*x*uvw(1)
    c = y*y + z*z - this%r2*x*x
    quad = k*k - a*c

    if (quad < ZERO) then
      ! no intersection with cone

      d = INFINITY

    elseif (coincident .or. abs(c) < FP_COINCIDENT) then
      ! particle is on the cone, thus one distance is positive/negative and the
      ! other is zero. The sign of k determines which distance is zero and which
      ! is not.

      if (k >= ZERO) then
        d = (-k - sqrt(quad))/a
      else
        d = (-k + sqrt(quad))/a
      end if

    else
      ! calculate both solutions to the quadratic
      quad = sqrt(quad)
      d = (-k - quad)/a
      b = (-k + quad)/a

      ! determine the smallest positive solution
      if (d < ZERO) then
        if (b > ZERO) then
          d = b
        end if
      else
        if (b > ZERO) d = min(d, b)
      end if
    end if

    ! If the distance was negative, set boundary distance to infinity
    if (d <= ZERO) d = INFINITY
  end function x_cone_distance

  subroutine x_cone_reflect(this, xyz, uvw)
    class(SurfaceXCone), intent(in)    :: this
    real(8),        intent(in)    :: xyz(3)
    real(8),        intent(inout) :: uvw(3)

    real(8) :: x, y, z, r, dot_prod

    ! Find x-x0, y-y0, z-z0 and dot product of direction and surface normal
    x = xyz(1) - this%x0
    y = xyz(2) - this%y0
    z = xyz(3) - this%z0
    r = this%r2
    dot_prod = (uvw(2)*y + uvw(3)*z - r*uvw(1)*x)/((r + ONE)*r*x*x)

    ! Reflect direction according to normal
    uvw(1) = uvw(1) + TWO*dot_prod*r*x
    uvw(2) = uvw(2) - TWO*dot_prod*y
    uvw(3) = uvw(3) - TWO*dot_prod*z
  end subroutine x_cone_reflect

  pure function x_cone_normal(this, xyz) result(uvw)
    class(SurfaceXCone), intent(in) :: this
    real(8), intent(in) :: xyz(3)
    real(8) :: uvw(3)

    uvw(1) = -TWO*this%r2*(xyz(1) - this%x0)
    uvw(2) = TWO*(xyz(2) - this%y0)
    uvw(3) = TWO*(xyz(3) - this%z0)
  end function x_cone_normal

!===============================================================================
! SurfaceYCone Implementation
!===============================================================================

  pure function y_cone_evaluate(this, xyz) result(f)
    class(SurfaceYCone), intent(in) :: this
    real(8), intent(in) :: xyz(3)
    real(8) :: f

    real(8) :: x, y, z

    x = xyz(1) - this%x0
    y = xyz(2) - this%y0
    z = xyz(3) - this%z0
    f = x*x + z*z - this%r2*y*y
  end function y_cone_evaluate

  pure function y_cone_distance(this, xyz, uvw, coincident) result(d)
    class(SurfaceYCone), intent(in) :: this
    real(8), intent(in) :: xyz(3)
    real(8), intent(in) :: uvw(3)
    logical, intent(in) :: coincident
    real(8) :: d

    real(8) :: x, y, z, k, a, b, c, quad

    x = xyz(1) - this%x0
    y = xyz(2) - this%y0
    z = xyz(3) - this%z0
    a = uvw(1)*uvw(1) + uvw(3)*uvw(3) - this%r2*uvw(2)*uvw(2)
    k = x*uvw(1) + z*uvw(3) - this%r2*y*uvw(2)
    c = x*x + z*z - this%r2*y*y
    quad = k*k - a*c

    if (quad < ZERO) then
      ! no intersection with cone

      d = INFINITY

    elseif (coincident .or. abs(c) < FP_COINCIDENT) then
      ! particle is on the cone, thus one distance is positive/negative and the
      ! other is zero. The sign of k determines which distance is zero and which
      ! is not.

      if (k >= ZERO) then
        d = (-k - sqrt(quad))/a
      else
        d = (-k + sqrt(quad))/a
      end if

    else
      ! calculate both solutions to the quadratic
      quad = sqrt(quad)
      d = (-k - quad)/a
      b = (-k + quad)/a

      ! determine the smallest positive solution
      if (d < ZERO) then
        if (b > ZERO) then
          d = b
        end if
      else
        if (b > ZERO) d = min(d, b)
      end if
    end if

    ! If the distance was negative, set boundary distance to infinity
    if (d <= ZERO) d = INFINITY
  end function y_cone_distance

  subroutine y_cone_reflect(this, xyz, uvw)
    class(SurfaceYCone), intent(in) :: this
    real(8), intent(in)    :: xyz(3)
    real(8), intent(inout) :: uvw(3)

    real(8) :: x, y, z, r, dot_prod

    ! Find x-x0, y-y0, z-z0 and dot product of direction and surface normal
    x = xyz(1) - this%x0
    y = xyz(2) - this%y0
    z = xyz(3) - this%z0
    r = this%r2
    dot_prod = (uvw(1)*x + uvw(3)*z - r*uvw(2)*y)/((r + ONE)*r*y*y)

    ! Reflect direction according to normal
    uvw(1) = uvw(1) - TWO*dot_prod*x
    uvw(2) = uvw(2) + TWO*dot_prod*r*y
    uvw(3) = uvw(3) - TWO*dot_prod*z
  end subroutine y_cone_reflect

  pure function y_cone_normal(this, xyz) result(uvw)
    class(SurfaceYCone), intent(in) :: this
    real(8), intent(in) :: xyz(3)
    real(8) :: uvw(3)

    uvw(1) = TWO*(xyz(1) - this%x0)
    uvw(2) = -TWO*this%r2*(xyz(2) - this%y0)
    uvw(3) = TWO*(xyz(3) - this%z0)
  end function y_cone_normal

!===============================================================================
! SurfaceZConeImplementation
!===============================================================================

  pure function z_cone_evaluate(this, xyz) result(f)
    class(SurfaceZCone), intent(in) :: this
    real(8), intent(in) :: xyz(3)
    real(8) :: f

    real(8) :: x, y, z

    x = xyz(1) - this%x0
    y = xyz(2) - this%y0
    z = xyz(3) - this%z0
    f = x*x + y*y - this%r2*z*z
  end function z_cone_evaluate

  pure function z_cone_distance(this, xyz, uvw, coincident) result(d)
    class(SurfaceZCone), intent(in) :: this
    real(8), intent(in) :: xyz(3)
    real(8), intent(in) :: uvw(3)
    logical, intent(in) :: coincident
    real(8) :: d

    real(8) :: x, y, z, k, a, b, c, quad

    x = xyz(1) - this%x0
    y = xyz(2) - this%y0
    z = xyz(3) - this%z0
    a = uvw(1)*uvw(1) + uvw(2)*uvw(2) - this%r2*uvw(3)*uvw(3)
    k = x*uvw(1) + y*uvw(2) - this%r2*z*uvw(3)
    c = x*x + y*y - this%r2*z*z
    quad = k*k - a*c

    if (quad < ZERO) then
      ! no intersection with cone

      d = INFINITY

    elseif (coincident .or. abs(c) < FP_COINCIDENT) then
      ! particle is on the cone, thus one distance is positive/negative and the
      ! other is zero. The sign of k determines which distance is zero and which
      ! is not.

      if (k >= ZERO) then
        d = (-k - sqrt(quad))/a
      else
        d = (-k + sqrt(quad))/a
      end if

    else
      ! calculate both solutions to the quadratic
      quad = sqrt(quad)
      d = (-k - quad)/a
      b = (-k + quad)/a

      ! determine the smallest positive solution
      if (d < ZERO) then
        if (b > ZERO) then
          d = b
        end if
      else
        if (b > ZERO) d = min(d, b)
      end if
    end if

    ! If the distance was negative, set boundary distance to infinity
    if (d <= ZERO) d = INFINITY
  end function z_cone_distance

  subroutine z_cone_reflect(this, xyz, uvw)
    class(SurfaceZCone), intent(in) :: this
    real(8), intent(in)    :: xyz(3)
    real(8), intent(inout) :: uvw(3)

    real(8) :: x, y, z, r, dot_prod

    ! Find x-x0, y-y0, z-z0 and dot product of direction and surface normal
    x = xyz(1) - this%x0
    y = xyz(2) - this%y0
    z = xyz(3) - this%z0
    r = this%r2
    dot_prod = (uvw(1)*x + uvw(2)*y - r*uvw(3)*z)/((r + ONE)*r*z*z)

    ! Reflect direction according to normal
    uvw(1) = uvw(1) - TWO*dot_prod*x
    uvw(2) = uvw(2) - TWO*dot_prod*y
    uvw(3) = uvw(3) + TWO*dot_prod*r*z
  end subroutine z_cone_reflect

  pure function z_cone_normal(this, xyz) result(uvw)
    class(SurfaceZCone), intent(in) :: this
    real(8), intent(in) :: xyz(3)
    real(8) :: uvw(3)

    uvw(1) = TWO*(xyz(1) - this%x0)
    uvw(2) = TWO*(xyz(2) - this%y0)
    uvw(3) = -TWO*this%r2*(xyz(3) - this%z0)
  end function z_cone_normal

end module surface_header
