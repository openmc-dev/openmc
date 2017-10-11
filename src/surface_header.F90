module surface_header

  use, intrinsic :: ISO_C_BINDING

  use constants, only: NONE, ONE, TWO, ZERO, HALF, INFINITY, FP_COINCIDENT
  use dict_header, only: DictIntInt

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
    integer :: i_periodic = NONE      ! Index of corresponding periodic surface
    character(len=104) :: name = ""   ! User-defined name
  contains
    procedure :: reflect
    procedure(surface_evaluate_), deferred :: evaluate
    procedure(surface_normal_),   deferred :: normal
  end type Surface

  abstract interface
    pure function surface_evaluate_(this, xyz) result(f)
      import Surface
      class(Surface), intent(in) :: this
      real(8), intent(in) :: xyz(3)
      real(8) :: f
    end function surface_evaluate_

    pure function surface_normal_(this, xyz) result(uvw)
      import Surface
      class(Surface), intent(in) :: this
      real(8), intent(in) :: xyz(3)
      real(8) :: uvw(3)
    end function surface_normal_
  end interface

!===============================================================================
! SURFACECONTAINER allows us to store an array of different types of surfaces
!===============================================================================

  type :: SurfaceContainer
    class(Surface), allocatable :: obj
  end type SurfaceContainer

!===============================================================================
! All the derived types below are extensions of the abstract Surface type. They
! inherent the reflect() type-bound procedure and must implement evaluate(),
! and normal()
!===============================================================================

  type, extends(Surface) :: SurfaceXPlane
    ! x = x0
    real(8) :: x0
  contains
    procedure :: evaluate => x_plane_evaluate
    procedure :: normal => x_plane_normal
  end type SurfaceXPlane

  type, extends(Surface) :: SurfaceYPlane
    ! y = y0
    real(8) :: y0
  contains
    procedure :: evaluate => y_plane_evaluate
    procedure :: normal => y_plane_normal
  end type SurfaceYPlane

  type, extends(Surface) :: SurfaceZPlane
    ! z = z0
    real(8) :: z0
  contains
    procedure :: evaluate => z_plane_evaluate
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
    procedure :: normal => plane_normal
  end type SurfacePlane

  type, extends(Surface) :: SurfaceXCylinder
    ! (y - y0)^2 + (z - z0)^2 = R^2
    real(8) :: y0
    real(8) :: z0
    real(8) :: r
  contains
    procedure :: evaluate => x_cylinder_evaluate
    procedure :: normal => x_cylinder_normal
  end type SurfaceXCylinder

  type, extends(Surface) :: SurfaceYCylinder
    ! (x - x0)^2 + (z - z0)^2 = R^2
    real(8) :: x0
    real(8) :: z0
    real(8) :: r
  contains
    procedure :: evaluate => y_cylinder_evaluate
    procedure :: normal => y_cylinder_normal
  end type SurfaceYCylinder

  type, extends(Surface) :: SurfaceZCylinder
    ! (x - x0)^2 + (y - y0)^2 = R^2
    real(8) :: x0
    real(8) :: y0
    real(8) :: r
  contains
    procedure :: evaluate => z_cylinder_evaluate
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
    procedure :: normal => z_cone_normal
  end type SurfaceZCone

  type, extends(Surface) :: SurfaceQuadric
    ! Ax^2 + By^2 + Cz^2 + Dxy + Eyz + Fxz + Gx + Hy + Jz + K = 0
    real(8) :: A, B, C, D, E, F, G, H, J, K
  contains
    procedure :: evaluate => quadric_evaluate
    procedure :: normal => quadric_normal
  end type SurfaceQuadric

  integer(C_INT32_T), bind(C) :: n_surfaces  ! # of surfaces

  type(SurfaceContainer), allocatable, target :: surfaces(:)

  ! Dictionary that maps user IDs to indices in 'surfaces'
  type(DictIntInt) :: surface_dict

contains

!===============================================================================
! REFLECT determines the direction a particle will travel if it is specularly
! reflected from the surface at a given position and direction. The position is
! needed because the reflection is performed using the surface normal, which
! depends on the position for second-order surfaces.
!===============================================================================

  pure subroutine reflect(this, xyz, uvw)
    class(Surface), intent(in) :: this
    real(8), intent(in) :: xyz(3)
    real(8), intent(inout) :: uvw(3)

    real(8) :: projection
    real(8) :: magnitude
    real(8) :: n(3)

    ! Construct normal vector
    n(:) = this%normal(xyz)

    ! Determine projection of direction onto normal and squared magnitude of
    ! normal
    projection = n(1)*uvw(1) + n(2)*uvw(2) + n(3)*uvw(3)
    magnitude = n(1)*n(1) + n(2)*n(2) + n(3)*n(3)

    ! Reflect direction according to normal
    uvw(:) = uvw - TWO*projection/magnitude * n
  end subroutine reflect

!===============================================================================
! SurfaceXPlane Implementation
!===============================================================================

  pure function x_plane_evaluate(this, xyz) result(f)
    class(SurfaceXPlane), intent(in) :: this
    real(8), intent(in) :: xyz(3)
    real(8) :: f

    f = xyz(1) - this%x0
  end function x_plane_evaluate

  pure function x_plane_normal(this, xyz) result(uvw)
    class(SurfaceXPlane), intent(in) :: this
    real(8), intent(in) :: xyz(3)
    real(8) :: uvw(3)

    uvw(:) = [ONE, ZERO, ZERO]
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

  pure function y_plane_normal(this, xyz) result(uvw)
    class(SurfaceYPlane), intent(in) :: this
    real(8), intent(in) :: xyz(3)
    real(8) :: uvw(3)

    uvw(:) = [ZERO, ONE, ZERO]
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

  pure function z_plane_normal(this, xyz) result(uvw)
    class(SurfaceZPlane), intent(in) :: this
    real(8), intent(in) :: xyz(3)
    real(8) :: uvw(3)

    uvw(:) = [ZERO, ZERO, ONE]
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

  pure function z_cylinder_normal(this, xyz) result(uvw)
    class(SurfaceZCylinder), intent(in) :: this
    real(8), intent(in) :: xyz(3)
    real(8) :: uvw(3)

    uvw(1) = TWO*(xyz(1) - this%x0)
    uvw(2) = TWO*(xyz(2) - this%y0)
    uvw(3) = ZERO
  end function z_cylinder_normal

!===============================================================================
! SurfaceSphere Implementation
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

  pure function y_cone_normal(this, xyz) result(uvw)
    class(SurfaceYCone), intent(in) :: this
    real(8), intent(in) :: xyz(3)
    real(8) :: uvw(3)

    uvw(1) = TWO*(xyz(1) - this%x0)
    uvw(2) = -TWO*this%r2*(xyz(2) - this%y0)
    uvw(3) = TWO*(xyz(3) - this%z0)
  end function y_cone_normal

!===============================================================================
! SurfaceZCone Implementation
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

  pure function z_cone_normal(this, xyz) result(uvw)
    class(SurfaceZCone), intent(in) :: this
    real(8), intent(in) :: xyz(3)
    real(8) :: uvw(3)

    uvw(1) = TWO*(xyz(1) - this%x0)
    uvw(2) = TWO*(xyz(2) - this%y0)
    uvw(3) = -TWO*this%r2*(xyz(3) - this%z0)
  end function z_cone_normal

!===============================================================================
! SurfaceQuadric Implementation
!===============================================================================

  pure function quadric_evaluate(this, xyz) result(f)
    class(SurfaceQuadric), intent(in) :: this
    real(8), intent(in) :: xyz(3)
    real(8) :: f

    associate (x => xyz(1), y => xyz(2), z => xyz(3))
      f = x*(this%A*x + this%D*y + this%G) + &
           y*(this%B*y + this%E*z + this%H) + &
           z*(this%C*z + this%F*x + this%J) + this%K
    end associate
  end function quadric_evaluate

  pure function quadric_normal(this, xyz) result(uvw)
    class(SurfaceQuadric), intent(in) :: this
    real(8), intent(in) :: xyz(3)
    real(8) :: uvw(3)

    associate (x => xyz(1), y => xyz(2), z => xyz(3))
      uvw(1) = TWO*this%A*x + this%D*y + this%F*z + this%G
      uvw(2) = TWO*this%B*y + this%D*x + this%E*z + this%H
      uvw(3) = TWO*this%C*z + this%E*y + this%F*x + this%J
    end associate
  end function quadric_normal

!===============================================================================
! FREE_MEMORY_SURFACES deallocates global arrays defined in this module
!===============================================================================

  subroutine free_memory_surfaces()
    n_surfaces = 0
    if (allocated(surfaces)) deallocate(surfaces)
    call surface_dict % clear()
  end subroutine free_memory_surfaces

end module surface_header
