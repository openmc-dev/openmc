module math

  use, intrinsic :: ISO_C_BINDING

  use constants
  use random_lcg, only: prn

  implicit none
  private
  public :: t_percentile
  public :: calc_pn
  public :: calc_rn
  public :: calc_zn
  public :: evaluate_legendre
  public :: rotate_angle
  public :: maxwell_spectrum
  public :: watt_spectrum
  public :: faddeeva
  public :: w_derivative
  public :: broaden_wmp_polynomials

!===============================================================================
! FADDEEVA_W evaluates the scaled complementary error function.  This
! interfaces with the MIT C library
!===============================================================================

  interface

    pure function t_percentile_cc(p, df) bind(C, name='t_percentile_c') &
         result(t)
      use ISO_C_BINDING
      implicit none
      real(C_DOUBLE), value, intent(in) :: p
      integer(C_INT), value, intent(in) :: df
      real(C_DOUBLE) :: t
    end function t_percentile_cc

    pure function calc_pn_cc(n, x) bind(C, name='calc_pn_c') result(pnx)
      use ISO_C_BINDING
      implicit none
      integer(C_INT), value, intent(in) :: n
      real(C_DOUBLE), value, intent(in) :: x
      real(C_DOUBLE) :: pnx
    end function calc_pn_cc

    subroutine calc_rn_cc(n, uvw, rn) bind(C, name='calc_rn_c')
      use ISO_C_BINDING
      implicit none
      integer(C_INT), value, intent(in) :: n
      real(C_DOUBLE), intent(in) :: uvw(3)
      real(C_DOUBLE), intent(in) :: rn(2 * n + 1)
    end subroutine calc_rn_cc

    pure function evaluate_legendre_cc(n, data, x) &
         bind(C, name='evaluate_legendre_c') result(val)
      use ISO_C_BINDING
      implicit none
      integer(C_INT), value, intent(in) :: n
      real(C_DOUBLE), intent(in) :: data(n)
      real(C_DOUBLE), value, intent(in) :: x
      real(C_DOUBLE) :: val
    end function evaluate_legendre_cc

    subroutine rotate_angle_cc(uvw, mu, phi) bind(C, name='rotate_angle_c')
      use ISO_C_BINDING
      implicit none
      real(C_DOUBLE), intent(inout) :: uvw(3)
      real(C_DOUBLE), value, intent(in)    :: mu
      real(C_DOUBLE), value, intent(in) :: phi
    end subroutine rotate_angle_cc

    function maxwell_spectrum_cc(T) bind(C, name='maxwell_spectrum_c') &
         result(E_out)
      use ISO_C_BINDING
      implicit none
      real(C_DOUBLE), value, intent(in) :: T
      real(C_DOUBLE) :: E_out
    end function maxwell_spectrum_cc

    function watt_spectrum_cc(a, b) bind(C, name='watt_spectrum_c') &
         result(E_out)
      use ISO_C_BINDING
      implicit none
      real(C_DOUBLE), value, intent(in) :: a
      real(C_DOUBLE), value, intent(in) :: b
      real(C_DOUBLE) :: E_out
    end function watt_spectrum_cc

    function faddeeva_cc(z) bind(C, name='faddeeva_c') result(wv)
      use ISO_C_BINDING
      implicit none
      complex(C_DOUBLE_COMPLEX), value, intent(in) :: z
      complex(C_DOUBLE_COMPLEX) :: wv
    end function faddeeva_cc

    function w_derivative_cc(z, order) bind(C, name='w_derivative_c') &
           result(wv)
      use ISO_C_BINDING
      implicit none
      complex(C_DOUBLE_COMPLEX), value, intent(in) :: z
      integer(C_INT), value,            intent(in) :: order
      complex(C_DOUBLE_COMPLEX) :: wv
    end function w_derivative_cc

    subroutine broaden_wmp_polynomials_cc(E, dopp, n, factors) &
           bind(C, name='broaden_wmp_polynomials_c')
      use ISO_C_BINDING
      implicit none
      real(C_DOUBLE), value, intent(in) :: E
      real(C_DOUBLE), value, intent(in) :: dopp
      integer(C_INT), value, intent(in) :: n
      real(C_DOUBLE), intent(inout) :: factors(n)
    end subroutine broaden_wmp_polynomials_cc

    function faddeeva_w(z, relerr) bind(C, name='Faddeeva_w') result(w)
      use ISO_C_BINDING
      implicit none
      complex(C_DOUBLE_COMPLEX), value :: z
      real(C_DOUBLE),            value :: relerr
      complex(C_DOUBLE_COMPLEX)        :: w
    end function faddeeva_w
  end interface

contains

!===============================================================================
! T_PERCENTILE calculates the percentile of the Student's t distribution with a
! specified probability level and number of degrees of freedom
!===============================================================================

  pure function t_percentile(p, df) result(t) bind(C)

    real(C_DOUBLE), intent(in)  :: p  ! probability level
    integer(C_INT), intent(in)  :: df ! degrees of freedom
    real(C_DOUBLE)              :: t  ! corresponding t-value

    t = t_percentile_cc(p, df)

  end function t_percentile

!===============================================================================
! CALC_PN calculates the n-th order Legendre polynomial at the value of x.
! Since this function is called repeatedly during the neutron transport process,
! neither n or x is checked to see if they are in the applicable range.
! This is left to the client developer to use where applicable. x is to be in
! the domain of [-1,1], and 0<=n<=5. If x is outside of the range, the return
! value will be outside the expected range; if n is outside the stated range,
! the return value will be 1.0.
!===============================================================================

  pure function calc_pn(n, x) result(pnx) bind(C)

    integer(C_INT), intent(in) :: n   ! Legendre order requested
    real(C_DOUBLE), intent(in) :: x   ! Independent variable the Legendre is to
                                      ! be evaluated at; x must be in the
                                      ! domain [-1,1]
    real(C_DOUBLE)             :: pnx ! The Legendre poly of order n evaluated
                                      ! at x

    pnx = calc_pn_cc(n, x)

  end function calc_pn

!===============================================================================
! EVALUATE_LEGENDRE Find the value of f(x) given a set of Legendre coefficients
! and the value of x
!===============================================================================

  pure function evaluate_legendre(n, data, x) result(val) bind(C)
    integer(C_INT), intent(in) :: n
    real(C_DOUBLE), intent(in) :: data(n)
    real(C_DOUBLE), intent(in) :: x
    real(C_DOUBLE)             :: val

    val = evaluate_legendre_cc(size(data), data, x)

  end function evaluate_legendre

!===============================================================================
! CALC_RN calculates the n-th order real spherical harmonics for a given angle
! (in terms of (u,v,w)).  All Rn,m values are provided (where -n<=m<=n)
!===============================================================================

  subroutine calc_rn(n, uvw, rn) bind(C)

    integer(C_INT), intent(in) :: n      ! Order requested
    real(C_DOUBLE), intent(in) :: uvw(3) ! Direction of travel;
                                         ! assumed to be on unit sphere
    real(C_DOUBLE)             :: rn(2*n + 1) ! The resultant R_n(uvw)

    call calc_rn_cc(n, uvw, rn)

  end subroutine calc_rn

!===============================================================================
! CALC_ZN calculates the n-th order modified Zernike polynomial moment for a
! given angle (rho, theta) location in the unit disk. The normlization of the
! polynomials is such that the integral of Z_pq*Z_pq over the unit disk is
! exactly pi
!===============================================================================

  subroutine calc_zn(n, rho, phi, zn) bind(C)
    ! This procedure uses the modified Kintner's method for calculating Zernike
    ! polynomials as outlined in Chong, C. W., Raveendran, P., & Mukundan,
    ! R. (2003). A comparative analysis of algorithms for fast computation of
    ! Zernike moments. Pattern Recognition, 36(3), 731-742.

    integer(C_INT), intent(in) :: n      ! Maximum order
    real(C_DOUBLE), intent(in) :: rho    ! Radial location in the unit disk
    real(C_DOUBLE), intent(in) :: phi    ! Theta (radians) location in the unit disk
    real(C_DOUBLE), intent(out) :: zn(:) ! The resulting list of coefficients

    real(C_DOUBLE) :: sin_phi, cos_phi   ! Sine and Cosine of phi
    real(C_DOUBLE) :: sin_phi_vec(n+1)   ! Contains sin(n*phi)
    real(C_DOUBLE) :: cos_phi_vec(n+1)   ! Contains cos(n*phi)
    real(C_DOUBLE) :: zn_mat(n+1, n+1)   ! Matrix form of the coefficients which is
                                         ! easier to work with
    real(C_DOUBLE) :: k1, k2, k3, k4     ! Variables for R_m_n calculation
    integer(C_INT) :: i,p,q              ! Loop counters

    real(C_DOUBLE), parameter :: SQRT_N_1(0:10) = [&
         sqrt(1.0_8), sqrt(2.0_8), sqrt(3.0_8), sqrt(4.0_8), &
         sqrt(5.0_8), sqrt(6.0_8), sqrt(7.0_8), sqrt(8.0_8), &
         sqrt(9.0_8), sqrt(10.0_8), sqrt(11.0_8)]
    real(C_DOUBLE), parameter :: SQRT_2N_2(0:10) = SQRT_N_1*sqrt(2.0_8)

    ! n == radial degree
    ! m == azimuthal frequency

    ! ==========================================================================
    ! Determine vector of sin(n*phi) and cos(n*phi). This takes advantage of the
    ! following recurrence relations so that only a single sin/cos have to be
    ! evaluated (http://mathworld.wolfram.com/Multiple-AngleFormulas.html)
    !
    ! sin(nx) = 2 cos(x) sin((n-1)x) - sin((n-2)x)
    ! cos(nx) = 2 cos(x) cos((n-1)x) - cos((n-2)x)

    sin_phi = sin(phi)
    cos_phi = cos(phi)

    sin_phi_vec(1) = 1.0_8
    cos_phi_vec(1) = 1.0_8

    sin_phi_vec(2) = 2.0_8 * cos_phi
    cos_phi_vec(2) = cos_phi

    do i = 3, n+1
      sin_phi_vec(i) = 2.0_8 * cos_phi * sin_phi_vec(i-1) - sin_phi_vec(i-2)
      cos_phi_vec(i) = 2.0_8 * cos_phi * cos_phi_vec(i-1) - cos_phi_vec(i-2)
    end do

    do i = 1, n+1
      sin_phi_vec(i) = sin_phi_vec(i) * sin_phi
    end do

    ! ==========================================================================
    ! Calculate R_pq(rho)

    ! Fill the main diagonal first (Eq. 3.9 in Chong)
    do p = 0, n
      zn_mat(p+1, p+1) = rho**p
    end do

    ! Fill in the second diagonal (Eq. 3.10 in Chong)
    do q = 0, n-2
      zn_mat(q+2+1, q+1) = (q+2) * zn_mat(q+2+1, q+2+1) - (q+1) * zn_mat(q+1, q+1)
    end do

    ! Fill in the rest of the values using the original results (Eq. 3.8 in Chong)
    do p = 4, n
      k2 = 2 * p * (p - 1) * (p - 2)
      do q = p-4, 0, -2
        k1 = (p + q) * (p - q) * (p - 2) / 2
        k3 = -q**2*(p - 1) - p * (p - 1) * (p - 2)
        k4 = -p * (p + q - 2) * (p - q - 2) / 2
        zn_mat(p+1, q+1) = ((k2 * rho**2 + k3) * zn_mat(p-2+1, q+1) + k4 * zn_mat(p-4+1, q+1)) / k1
      end do
    end do

    ! Roll into a single vector for easier computation later
    ! The vector is ordered (0,0), (1,-1), (1,1), (2,-2), (2,0),
    ! (2, 2), ....   in (n,m) indices
    ! Note that the cos and sin vectors are offset by one
    ! sin_phi_vec = [sin(x), sin(2x), sin(3x) ...]
    ! cos_phi_vec = [1.0, cos(x), cos(2x)... ]
    i = 1
    do p = 0, n
      do q = -p, p, 2
        if (q < 0) then
          zn(i) = zn_mat(p+1, abs(q)+1) * sin_phi_vec(abs(q)) * SQRT_2N_2(p)
        else if (q == 0) then
          zn(i) = zn_mat(p+1, q+1) * SQRT_N_1(p)
        else
          zn(i) = zn_mat(p+1, q+1) * cos_phi_vec(abs(q)+1) * SQRT_2N_2(p)
        end if
        i = i + 1
      end do
    end do
  end subroutine calc_zn

!===============================================================================
! ROTATE_ANGLE rotates direction cosines through a polar angle whose cosine is
! mu and through an azimuthal angle sampled uniformly. Note that this is done
! with direct sampling rather than rejection as is done in MCNP and SERPENT.
!===============================================================================

  subroutine rotate_angle(uvw0, mu, uvw, phi) bind(C)
    real(C_DOUBLE), intent(in)  :: uvw0(3) ! directional cosine
    real(C_DOUBLE), intent(in)  :: mu      ! cosine of angle in lab or CM
    real(C_DOUBLE), intent(out) :: uvw(3)  ! rotated directional cosine
    real(C_DOUBLE), optional    :: phi     ! azimuthal angle

    real(C_DOUBLE) :: phi_   ! azimuthal angle
    real(C_DOUBLE) :: sinphi ! sine of azimuthal angle
    real(C_DOUBLE) :: cosphi ! cosine of azimuthal angle
    real(C_DOUBLE) :: a      ! sqrt(1 - mu^2)
    real(C_DOUBLE) :: b      ! sqrt(1 - w^2)
    real(C_DOUBLE) :: u0     ! original cosine in x direction
    real(C_DOUBLE) :: v0     ! original cosine in y direction
    real(C_DOUBLE) :: w0     ! original cosine in z direction

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

  end subroutine rotate_angle

!===============================================================================
! MAXWELL_SPECTRUM samples an energy from the Maxwell fission distribution based
! on a direct sampling scheme. The probability distribution function for a
! Maxwellian is given as p(x) = 2/(T*sqrt(pi))*sqrt(x/T)*exp(-x/T). This PDF can
! be sampled using rule C64 in the Monte Carlo Sampler LA-9721-MS.
!===============================================================================

  function maxwell_spectrum(T) result(E_out) bind(C)

    real(C_DOUBLE), intent(in)  :: T     ! tabulated function of incoming E
    real(C_DOUBLE)              :: E_out ! sampled energy

    real(C_DOUBLE) :: r1, r2, r3  ! random numbers
    real(C_DOUBLE) :: c           ! cosine of pi/2*r3

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

  function watt_spectrum(a, b) result(E_out) bind(C)

    real(C_DOUBLE), intent(in) :: a     ! Watt parameter a
    real(C_DOUBLE), intent(in) :: b     ! Watt parameter b
    real(C_DOUBLE)             :: E_out ! energy of emitted neutron

    real(C_DOUBLE) :: w ! sampled from Maxwellian

    w     = maxwell_spectrum(a)
    E_out = w + a*a*b/4. + (TWO*prn() - ONE)*sqrt(a*a*b*w)

  end function watt_spectrum

!===============================================================================
! FADDEEVA the Faddeeva function, using Stephen Johnson's implementation
!===============================================================================

  function faddeeva(z) result(wv) bind(C)
    complex(C_DOUBLE_COMPLEX), intent(in) :: z  ! The point to evaluate Z at
    complex(C_DOUBLE_COMPLEX)             :: wv ! The resulting w(z) value
    real(C_DOUBLE) :: relerr ! Target relative error in inner loop of MIT
                             !  Faddeeva

    ! Technically, the value we want is given by the equation:
    ! w(z) = I/Pi * Integrate[Exp[-t^2]/(z-t), {t, -Infinity, Infinity}]
    ! as shown in Equation 63 from Hwang, R. N. "A rigorous pole
    ! representation of multilevel cross sections and its practical
    ! applications." Nuclear Science and Engineering 96.3 (1987): 192-209.
    !
    ! The MIT Faddeeva function evaluates w(z) = exp(-z^2)erfc(-iz). These
    ! two forms of the Faddeeva function are related by a transformation.
    !
    ! If we call the integral form w_int, and the function form w_fun:
    ! For imag(z) > 0, w_int(z) = w_fun(z)
    ! For imag(z) < 0, w_int(z) = -conjg(w_fun(conjg(z)))

    ! Note that faddeeva_w will interpret zero as machine epsilon

    relerr = ZERO
    if (aimag(z) > ZERO) then
      wv = faddeeva_w(z, relerr)
    else
      wv = -conjg(faddeeva_w(conjg(z), relerr))
    end if

  end function faddeeva

  recursive function w_derivative(z, order) result(wv) bind(C)
    complex(C_DOUBLE_COMPLEX), intent(in) :: z ! The point to evaluate Z at
    integer(C_INT),            intent(in) :: order
    complex(C_DOUBLE_COMPLEX)     :: wv     ! The resulting w(z) value

    select case(order)
    case (0)
      wv = faddeeva(z)
    case (1)
      wv = -TWO * z * faddeeva(z) + TWO * ONEI / SQRT_PI
    case default
      wv = -TWO * z * w_derivative(z, order-1) &
           - TWO * (order-1) * w_derivative(z, order-2)
    end select
  end function w_derivative

!===============================================================================
! BROADEN_WMP_POLYNOMIALS Doppler broadens the windowed multipole curvefit.  The
! curvefit is a polynomial of the form
! a/E + b/sqrt(E) + c + d sqrt(E) ...
!===============================================================================

  subroutine broaden_wmp_polynomials(E, dopp, n, factors) bind(C)
    real(C_DOUBLE), intent(in) :: E          ! Energy to evaluate at
    real(C_DOUBLE), intent(in) :: dopp       ! sqrt(atomic weight ratio / kT),
                                      !  kT given in eV.
    integer(C_INT), intent(in) :: n          ! number of components to polynomial
    real(C_DOUBLE), intent(out):: factors(n) ! output leading coefficient

    integer :: i

    real(8) :: sqrtE               ! sqrt(energy)
    real(8) :: beta                ! sqrt(atomic weight ratio * E / kT)
    real(8) :: half_inv_dopp2      ! 0.5 / dopp**2
    real(8) :: quarter_inv_dopp4   ! 0.25 / dopp**4
    real(8) :: erf_beta            ! error function of beta
    real(8) :: exp_m_beta2         ! exp(-beta**2)

    sqrtE = sqrt(E)
    beta = sqrtE * dopp
    half_inv_dopp2 = HALF / dopp**2
    quarter_inv_dopp4 = half_inv_dopp2**2

    if (beta > 6.0_8) then
      ! Save time, ERF(6) is 1 to machine precision.
      ! beta/sqrtpi*exp(-beta**2) is also approximately 1 machine epsilon.
      erf_beta = ONE
      exp_m_beta2 = ZERO
    else
      erf_beta = erf(beta)
      exp_m_beta2 = exp(-beta**2)
    end if

    ! Assume that, for sure, we'll use a second order (1/E, 1/V, const)
    ! fit, and no less.

    factors(1) = erf_beta / E
    factors(2) = ONE / sqrtE
    factors(3) = factors(1) * (half_inv_dopp2 + E) &
         + exp_m_beta2 / (beta * SQRT_PI)

    ! Perform recursive broadening of high order components
    do i = 1, n-3
      if (i /= 1) then
        factors(i+3) = -factors(i-1) * (i - ONE) * i * quarter_inv_dopp4 &
             + factors(i+1) * (E + (ONE + TWO * i) * half_inv_dopp2)
      else
        ! Although it's mathematically identical, factors(0) will contain
        ! nothing, and we don't want to have to worry about memory.
        factors(i+3) = factors(i+1)*(E + (ONE + TWO * i) * half_inv_dopp2)
      end if
    end do
  end subroutine broaden_wmp_polynomials

end module math
