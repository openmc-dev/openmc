module math

  use, intrinsic :: ISO_C_BINDING

  use constants
  use random_lcg, only: prn

  implicit none
  private
  public :: t_percentile
  public :: calc_pn
  public :: calc_rn
  public :: rotate_angle
  public :: spline
  public :: spline_interpolate
  public :: spline_integrate

  interface

    pure function t_percentile(p, df) bind(C, name='t_percentile_c') &
         result(t)
      use ISO_C_BINDING
      implicit none
      real(C_DOUBLE), value, intent(in) :: p
      integer(C_INT), value, intent(in) :: df
      real(C_DOUBLE) :: t
    end function t_percentile

    pure subroutine calc_pn(n, x, pnx) bind(C, name='calc_pn_c')
      use ISO_C_BINDING
      implicit none
      integer(C_INT), value, intent(in) :: n
      real(C_DOUBLE), value, intent(in) :: x
      real(C_DOUBLE), intent(out) :: pnx(n + 1)
    end subroutine calc_pn

    pure subroutine calc_rn(n, uvw, rn) bind(C, name='calc_rn_c')
      use ISO_C_BINDING
      implicit none
      integer(C_INT), value, intent(in) :: n
      real(C_DOUBLE), intent(in)  :: uvw(3)
      real(C_DOUBLE), intent(out) :: rn(2 * n + 1)
    end subroutine calc_rn

    subroutine rotate_angle_c_intfc(uvw, mu, phi) bind(C, name='rotate_angle_c')
      use ISO_C_BINDING
      implicit none
      real(C_DOUBLE), intent(inout) :: uvw(3)
      real(C_DOUBLE), value, intent(in)    :: mu
      real(C_DOUBLE), optional, intent(in) :: phi
    end subroutine rotate_angle_c_intfc

    subroutine spline(n, x, y, z) bind(C, name='spline_c')
      use ISO_C_BINDING
      implicit none
      integer(C_INT), value, intent(in) :: n
      real(C_DOUBLE), intent(in)        :: x(n)
      real(C_DOUBLE), intent(in)        :: y(n)
      real(C_DOUBLE), intent(in)        :: z(n)
    end subroutine spline

    function spline_interpolate(n, x, y, z, xint) &
         bind(C, name='spline_interpolate_c') result(yint)
      use ISO_C_BINDING
      implicit none
      integer(C_INT), value, intent(in) :: n
      real(C_DOUBLE), intent(in)        :: x(n)
      real(C_DOUBLE), intent(in)        :: y(n)
      real(C_DOUBLE), intent(in)        :: z(n)
      real(C_DOUBLE), value, intent(in) :: xint
      real(C_DOUBLE)                    :: yint
    end function spline_interpolate

    function spline_integrate(n, x, y, z, xa, xb) &
         bind(C, name='spline_integrate_c') result(s)
      use ISO_C_BINDING
      implicit none
      integer(C_INT), value, intent(in) :: n
      real(C_DOUBLE), intent(in)        :: x(n)
      real(C_DOUBLE), intent(in)        :: y(n)
      real(C_DOUBLE), intent(in)        :: z(n)
      real(C_DOUBLE), value, intent(in) :: xa
      real(C_DOUBLE), value, intent(in) :: xb
      real(C_DOUBLE)                    :: s
    end function spline_integrate
  end interface

contains

!===============================================================================
! ROTATE_ANGLE rotates direction cosines through a polar angle whose cosine is
! mu and through an azimuthal angle sampled uniformly. Note that this is done
! with direct sampling rather than rejection as is done in MCNP and SERPENT.
!===============================================================================

  function rotate_angle(uvw0, mu, phi) result(uvw)
    real(C_DOUBLE), intent(in) :: uvw0(3)       ! directional cosine
    real(C_DOUBLE), intent(in) :: mu            ! cosine of angle in lab or CM
    real(C_DOUBLE), intent(in), optional :: phi ! azimuthal angle

    real(C_DOUBLE) :: uvw(3)  ! rotated directional cosine

    uvw = uvw0
    call rotate_angle_c_intfc(uvw, mu, phi)

  end function rotate_angle

end module math
