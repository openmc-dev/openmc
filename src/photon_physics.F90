module photon_physics

  use algorithm,       only: binary_search
  use constants
  use particle_header
  use photon_header,   only: PhotonInteraction, BremsstrahlungData, &
                             compton_profile_pz, ttb_e_grid, ttb
  use random_lcg,      only: prn
  use settings

contains

!===============================================================================
! KLEIN_NISHINA
!===============================================================================

  subroutine klein_nishina(alpha, alpha_out, mu)
    real(8), intent(in)  :: alpha
    real(8), intent(out) :: alpha_out
    real(8), intent(out) :: mu

    real(8) :: beta ! 1 + 2a
    real(8) :: t    ! (1 + 2a)/(9 + 2a)
    real(8) :: r, s, x
    real(8) :: gamma

    beta = ONE + TWO*alpha
    if (alpha < THREE) then
      ! Kahn's rejection method
      t = beta/(beta + 8.0_8)
      do
        if (prn() < t) then
          ! Left branch of flow chart
          r = TWO*prn()
          x = ONE + alpha*r
          if (prn() < FOUR/x*(ONE - ONE/x)) then
            mu = 1 - r
            exit
          end if
        else
          ! Right branch of flow chart
          x = beta/(ONE + TWO*alpha*prn())
          mu = ONE + (ONE - x)/alpha
          if (prn() < HALF*(mu**2 + ONE/x)) exit
        end if
      end do
      alpha_out = alpha/x

    else
      ! Koblinger's direct method
      gamma = ONE - beta**(-2)
      s = prn()*(FOUR/alpha + HALF*gamma + &
           (ONE - (ONE + beta)/alpha**2)*log(beta))
      if (s <= 2./alpha) then
        ! For first term, x = 1 + 2ar
        ! Therefore, a' = a/(1 + 2ar)
        alpha_out = alpha/(ONE + TWO*alpha*prn())
      elseif (s <= FOUR/alpha) then
        ! For third term, x = beta/(1 + 2ar)
        ! Therefore, a' = a(1 + 2ar)/beta
        alpha_out = alpha*(ONE + TWO*alpha*prn())/beta
      elseif (s <= FOUR/alpha + HALF*gamma) then
        ! For fourth term, x = 1/sqrt(1 - gamma*r)
        ! Therefore, a' = a*sqrt(1 - gamma*r)
        alpha_out = alpha*sqrt(ONE - gamma*prn())
      else
        ! For third term, x = beta^r
        ! Therefore, a' = a/beta^r
        alpha_out = alpha/beta**prn()
      end if

      ! Calculate cosine of scattering angle based on basic relation
      mu = ONE + ONE/alpha - ONE/alpha_out
    end if

  end subroutine klein_nishina

!===============================================================================
! COMPTON_SCATTER
!===============================================================================

  subroutine compton_scatter(el, alpha, alpha_out, mu, i_shell, use_doppler)
    type(PhotonInteraction), intent(in) :: el
    real(8), intent(in)  :: alpha
    real(8), intent(out) :: alpha_out
    real(8), intent(out) :: mu
    integer, intent(out) :: i_shell
    logical, intent(in), optional :: use_doppler

    real(8) :: x
    real(8) :: form_factor_xmax
    real(8) :: form_factor_x
    real(8) :: e_out
    logical :: use_doppler_

    if (present(use_doppler)) then
      use_doppler_ = use_doppler
    else
      use_doppler_ = .false.
    end if

    form_factor_xmax = ZERO
    do
      ! Sample Klein-Nishina distribution for trial energy and angle
      call klein_nishina(alpha, alpha_out, mu)

      ! Note that the parameter used here does not correspond exactly to the
      ! momentum transfer q in ENDF-102 Eq. (27.2). Rather, this is the
      ! parameter as defined by Hubbell, where the actual data comes from
      x = MASS_ELECTRON_EV/PLANCK_C*alpha*sqrt(HALF*(ONE - mu))

      ! Calculate S(x, Z) and S(x_max, Z)
      form_factor_x = el % incoherent_form_factor % evaluate(x)
      if (form_factor_xmax == ZERO) then
        form_factor_xmax = el % incoherent_form_factor % evaluate(&
             MASS_ELECTRON_EV/PLANCK_C*alpha)
      end if

      ! Perform rejection on form factor
      if (prn() < form_factor_x / form_factor_xmax) then
        if (use_doppler_) then
          call compton_doppler(el, alpha, mu, e_out, i_shell)
          alpha_out = e_out/MASS_ELECTRON_EV
        else
          i_shell = 0
        end if
        exit
      end if
    end do

  end subroutine compton_scatter

!===============================================================================
! COMPTON_DOPPLER
!===============================================================================

  subroutine compton_doppler(el, alpha, mu, e_out, i_shell)
    type(PhotonInteraction), intent(in) :: el
    real(8), intent(in)  :: alpha
    real(8), intent(in)  :: mu
    real(8), intent(out) :: e_out
    integer, intent(out) :: i_shell

    integer :: i
    integer :: n
    real(8) :: rn, m
    real(8) :: c, c_l, c_max
    real(8) :: pz_l, pz_r, pz, pz_max
    real(8) :: p_l, p_r
    real(8) :: e, e_b
    real(8) :: e_out1, e_out2
    real(8) :: a, b, quad
    real(8) :: f
    real(8) :: momentum_sq

    n = size(compton_profile_pz)

    do
      ! Sample electron shell
      rn = prn()
      c = ZERO
      do i_shell = 1, size(el % electron_pdf)
        c = c + el % electron_pdf(i_shell)
        if (rn < c) exit
      end do

      ! Determine binding energy of shell
      e_b = el % binding_energy(i_shell)

      ! Determine p_z,max
      e = alpha*MASS_ELECTRON_EV
      if (e < e_b) then
        e_out = alpha/(1 + alpha*(1 - mu))*MASS_ELECTRON_EV
        exit
      end if

      pz_max = -FINE_STRUCTURE*(e_b - (e - e_b)*alpha*(ONE - mu)) / &
           sqrt(TWO*e*(e - e_b)*(ONE - mu) + e_b**2)
      if (pz_max < ZERO) then
        e_out = alpha/(1 + alpha*(1 - mu))*MASS_ELECTRON_EV
        exit
      end if

      ! Determine profile cdf value corresponding to p_z,max
      if (pz_max > compton_profile_pz(n)) then
        c_max = el % profile_cdf(n, i_shell)
      else
        i = binary_search(compton_profile_pz, n, pz_max)
        pz_l = compton_profile_pz(i)
        pz_r = compton_profile_pz(i + 1)
        p_l = el % profile_pdf(i, i_shell)
        p_r = el % profile_pdf(i + 1, i_shell)
        c_l = el % profile_cdf(i, i_shell)
        if (pz_l == pz_r) then
          c_max = c_l
        elseif (p_l == p_r) then
          c_max = c_l + (pz_max - pz_l)*p_l
        else
          m = (p_l - p_r)/(pz_l - pz_r)
          c_max = c_l + ((m*(pz_max - pz_l) + p_l)**2 - p_l**2)/(TWO*m)
        end if
      end if

      ! Sample value on bounded cdf
      c = prn()*c_max

      ! Determine pz corresponding to sampled cdf value
      i = binary_search(el % profile_cdf(:, i_shell), n, c)
      pz_l = compton_profile_pz(i)
      pz_r = compton_profile_pz(i + 1)
      p_l = el % profile_pdf(i, i_shell)
      p_r = el % profile_pdf(i + 1, i_shell)
      c_l = el % profile_cdf(i, i_shell)
      if (pz_l == pz_r) then
        pz = pz_l
      elseif (p_l == p_r) then
        pz = pz_l + (c - c_l)/p_l
      else
        m = (p_l - p_r)/(pz_l - pz_r)
        pz = pz_l + (sqrt(p_l**2 + TWO*m*(c - c_l)) - p_l)/m
      end if

      ! Determine outgoing photon energy corresponding to electron momentum
      ! (solve Eq. 39 in LA-UR-04-0487 for E')
      momentum_sq = (pz/FINE_STRUCTURE)**2
      f = ONE + alpha*(ONE - mu)
      a = momentum_sq - f*f
      b = TWO*e*(f - momentum_sq*mu)
      c = e**2*(momentum_sq - ONE)

      quad = b**2 - FOUR*a*c
      if (quad < 0) then
        e_out = alpha/(1 + alpha*(1 - mu))*MASS_ELECTRON_EV
        exit
      end if
      quad = sqrt(quad)
      e_out1 = -(b + quad)/(TWO*a)
      e_out2 = -(b - quad)/(TWO*a)

      ! Determine solution to quadratic equation that is positive
      if (e_out1 > ZERO) then
        if (e_out2 > ZERO) then
          ! If both are positive, pick one at random
          if (prn() < HALF) then
            e_out = e_out1
          else
            e_out = e_out2
          end if
        else
          e_out = e_out1
        end if
      else
        if (e_out2 > ZERO) then
          e_out = e_out2
        else
          ! No positive solution -- resample
          cycle
        end if
      end if
      if (e_out < e - e_b) exit
    end do

  end subroutine compton_doppler

!===============================================================================
! RAYLEIGH_SCATTER
!===============================================================================

  subroutine rayleigh_scatter(el, alpha, mu)
    type(PhotonInteraction), intent(in) :: el
    real(8), intent(in)  :: alpha
    real(8), intent(out) :: mu

    integer :: i
    real(8) :: F
    real(8) :: F_max
    real(8) :: x2
    real(8) :: x2_max
    real(8) :: r

    do
      ! Determine maximum value of x^2
      x2_max = (MASS_ELECTRON_EV/PLANCK_C*alpha)**2

      ! Determine F(x^2_max, Z)
      F_max = el % coherent_int_form_factor % evaluate(x2_max)

      ! Sample cumulative distribution
      F = prn()*F_max

      ! Determine x^2 corresponding to F
      i = binary_search(el%coherent_int_form_factor%y, &
           size(el%coherent_int_form_factor%y), F)
      r = (F - el%coherent_int_form_factor%y(i)) / &
           (el%coherent_int_form_factor%y(i+1) - el%coherent_int_form_factor%y(i))
      x2 = el%coherent_int_form_factor%x(i) + r*(el%coherent_int_form_factor%x(i+1) - &
           el%coherent_int_form_factor%x(i))

      ! Calculate mu
      mu = ONE - TWO*x2/x2_max

      if (prn() < HALF*(ONE + mu**2)) exit
    end do

  end subroutine rayleigh_scatter

!===============================================================================
! ATOMIC_RELAXATION
!===============================================================================

  recursive subroutine atomic_relaxation(p, elm, i_shell)
    type(Particle), intent(inout) :: p
    type(PhotonInteraction), intent(in) :: elm
    integer, intent(in) :: i_shell

    integer :: i_hole
    integer :: i_transition
    integer :: primary
    integer :: secondary
    real(8) :: c
    real(8) :: rn
    real(8) :: E
    real(8) :: mu
    real(8) :: phi
    real(8) :: uvw(3)

    ! If no transitions, assume fluorescent photon from captured free electron
    if (elm % shells(i_shell) % n_transitions == 0) then
      mu = TWO*prn() - ONE
      phi = TWO*PI*prn()
      uvw(1) = mu
      uvw(2) = sqrt(ONE - mu*mu)*cos(phi)
      uvw(3) = sqrt(ONE - mu*mu)*sin(phi)
      E = elm % shells(i_shell) % binding_energy
      call particle_create_secondary(p, uvw, E, PHOTON, run_ce=.true._C_BOOL)
      return
    end if

    ! Sample transition
    rn = prn()
    c = ZERO
    do i_transition = 1, elm % shells(i_shell) % n_transitions
      c = c + elm % shells(i_shell) % &
           transition_probability(i_transition)
      if (rn < c) exit
    end do

    ! Get primary and secondary subshell designators
    primary = elm % shells(i_shell) % transition_subshells(1, i_transition)
    secondary = elm % shells(i_shell) % transition_subshells(2, i_transition)

    ! Sample angle isotropically
    mu = TWO*prn() - ONE
    phi = TWO*PI*prn()
    uvw(1) = mu
    uvw(2) = sqrt(ONE - mu*mu)*cos(phi)
    uvw(3) = sqrt(ONE - mu*mu)*sin(phi)

    ! Get the transition energy
    E = elm % shells(i_shell) % transition_energy(i_transition)

    if (secondary /= 0) then
      ! Non-radiative transition -- Auger/Coster-Kronig effect

      ! Create auger electron
      call particle_create_secondary(p, uvw, E, ELECTRON, run_ce=.true._C_BOOL)

      ! Fill hole left by emitted auger electron
      i_hole = elm % shell_dict % get(secondary)
      call atomic_relaxation(p, elm, i_hole)
    else
      ! Radiative transition -- get X-ray energy

      ! Create fluorescent photon
      call particle_create_secondary(p, uvw, E, PHOTON, run_ce=.true._C_BOOL)

    end if

    ! Fill hole created by electron transitioning to the photoelectron hole
    i_hole = elm % shell_dict % get(primary)
    call atomic_relaxation(p, elm, i_hole)

  end subroutine atomic_relaxation

!===============================================================================
! PAIR_PRODUCTION samples the kinetic energy and direction of the electron and
! positron created when a photon is absorbed near an atomic nucleus. The
! simulation procedure follows the semiempirical model outlined in F. Salvat, J.
! M. Fernández-Varea, and J. Sempau, "PENELOPE-2011: A Code System for Monte
! Carlo Simulation of Electron and Photon Transport," OECD-NEA,
! Issy-les-Moulineaux, France (2011).
!===============================================================================

  subroutine pair_production(elm, alpha, E_electron, E_positron, mu_electron, &
       mu_positron)
    type(PhotonInteraction), intent(in) :: elm
    real(8), intent(in)  :: alpha
    real(8), intent(out) :: E_electron
    real(8), intent(out) :: E_positron
    real(8), intent(out) :: mu_electron
    real(8), intent(out) :: mu_positron

    integer :: i
    real(8) :: f
    real(8) :: c
    real(8) :: a
    real(8) :: b
    real(8) :: q
    real(8) :: rn
    real(8) :: beta
    real(8) :: e, e_min, e_max
    real(8) :: t1, t2, t3, t4
    real(8) :: u1, u2
    real(8) :: phi1, phi2
    real(8) :: phi1_max, phi2_max
    real(8), parameter :: r(99) = (/ &
         122.81_8, 73.167_8, 69.228_8, 67.301_8, 64.696_8, 61.228_8, &
         57.524_8, 54.033_8, 50.787_8, 47.851_8, 46.373_8, 45.401_8, &
         44.503_8, 43.815_8, 43.074_8, 42.321_8, 41.586_8, 40.953_8, &
         40.524_8, 40.256_8, 39.756_8, 39.144_8, 38.462_8, 37.778_8, &
         37.174_8, 36.663_8, 35.986_8, 35.317_8, 34.688_8, 34.197_8, &
         33.786_8, 33.422_8, 33.068_8, 32.740_8, 32.438_8, 32.143_8, &
         31.884_8, 31.622_8, 31.438_8, 31.142_8, 30.950_8, 30.758_8, &
         30.561_8, 30.285_8, 30.097_8, 29.832_8, 29.581_8, 29.411_8, &
         29.247_8, 29.085_8, 28.930_8, 28.721_8, 28.580_8, 28.442_8, &
         28.312_8, 28.139_8, 27.973_8, 27.819_8, 27.675_8, 27.496_8, &
         27.285_8, 27.093_8, 26.911_8, 26.705_8, 26.516_8, 26.304_8, &
         26.108_8, 25.929_8, 25.730_8, 25.577_8, 25.403_8, 25.245_8, &
         25.100_8, 24.941_8, 24.790_8, 24.655_8, 24.506_8, 24.391_8, &
         24.262_8, 24.145_8, 24.039_8, 23.922_8, 23.813_8, 23.712_8, &
         23.621_8, 23.523_8, 23.430_8, 23.331_8, 23.238_8, 23.139_8, &
         23.048_8, 22.967_8, 22.833_8, 22.694_8, 22.624_8, 22.545_8, &
         22.446_8, 22.358_8, 22.264_8 /)

    ! The reduced screening radius r is the ratio of the screening radius to
    ! the Compton wavelength of the electron, where the screening radius is
    ! obtained under the assumption that the Coulomb field of the nucleus is
    ! exponentially screened by atomic electrons. This allows us to use a
    ! simplified atomic form factor and analytical approximations of the
    ! screening functions in the pair production DCS instead of computing the
    ! screening functions numerically. The reduced screening radii above for
    ! Z = 1-99 come from F. Salvat, J. M. Fernández-Varea, and J. Sempau,
    ! "PENELOPE-2011: A Code System for Monte Carlo Simulation of Electron and
    ! Photon Transport," OECD-NEA, Issy-les-Moulineaux, France (2011).

    ! Compute the minimum and maximum values of the electron reduced energy,
    ! i.e. the fraction of the photon energy that is given to the electron
    e_min = ONE/alpha
    e_max = ONE - ONE/alpha

    ! Compute the high-energy Coulomb correction
    a = elm % Z / FINE_STRUCTURE
    c = a**2*(ONE/(ONE + a**2) + 0.202059_8 - 0.03693_8*a**2 + 0.00835_8*a**4 &
         - 0.00201_8*a**6 + 0.00049_8*a**8 - 0.00012_8*a**10 + 0.00003_8*a**12)

    ! The analytical approximation of the DCS underestimates the cross section
    ! at low energies. The correction factor f compensates for this.
    q = sqrt(TWO/alpha)
    f = q*(-0.1774_8 - 12.10_8*a + 11.18_8*a**2) &
         + q**2*(8.523_8 + 73.26_8*a - 44.41_8*a**2) &
         + q**3*(-13.52_8 - 121.1_8*a + 96.41_8*a**2) &
         + q**4*(8.946_8 + 62.05_8*a - 63.41_8*a**2)

    ! Calculate phi_1(1/2) and phi_2(1/2). The unnormalized PDF for the reduced
    ! energy is given by p = 2*(1/2 - e)^2*phi_1(e) + phi_2(e), where phi_1 and
    ! phi_2 are non-negative and maximum at e = 1/2.
    b = TWO*r(elm % Z)/alpha
    t1 = TWO*log(ONE + b**2)
    t2 = b*atan(ONE/b)
    t3 = b**2*(FOUR - FOUR*t2 - THREE*log(ONE + ONE/b**2))
    t4 = FOUR*log(r(elm % Z)) - FOUR*c + f
    phi1_max = 7.0_8/THREE - t1 - 6.0_8*t2 - t3 + t4
    phi2_max = 11.0_8/6.0_8 - t1 - THREE*t2 + HALF*t3 + t4

    ! To aid sampling, the unnormalized PDF can be expressed as
    ! p = u_1*U_1(e)*pi_1(e) + u_2*U_2(e)*pi_2(e), where pi_1 and pi_2 are
    ! normalized PDFs on the interval (e_min, e_max) from which values of e can
    ! be sampled using the inverse transform method, and
    ! U_1 = phi_1(e)/phi_1(1/2) and U_2 = phi_2(e)/phi_2(1/2) are valid
    ! rejection functions. The reduced energy can now be sampled using a
    ! combination of the composition and rejection methods.
    u1 = TWO/THREE*(HALF - ONE/alpha)**2*phi1_max
    u2 = phi2_max
    do
      rn = prn()

      ! Sample the index i in (1, 2) using the point probabilities
      ! p(1) = u_1/(u_1 + u_2) and p(2) = u_2/(u_1 + u_2)
      if (prn() < u1/(u1 + u2)) then
        i = 1

        ! Sample e from pi_1 using the inverse transform method
        if (rn >= HALF) then
          e = HALF + (HALF - ONE/alpha)*(TWO*rn - ONE)**(ONE/THREE)
        else
          e = HALF - (HALF - ONE/alpha)*(ONE - TWO*rn)**(ONE/THREE)
        end if
      else
        i = 2

        ! Sample e from pi_2 using the inverse transform method
        e = ONE/alpha + (HALF - ONE/alpha)*TWO*rn
      end if

      ! Calculate phi_i(e) and deliver e if rn <= U_i(e)
      b = r(elm % Z)/(TWO*alpha*e*(ONE - e))
      t1 = TWO*log(ONE + b**2)
      t2 = b*atan(ONE/b)
      t3 = b**2*(FOUR - FOUR*t2 - THREE*log(ONE + ONE/b**2))
      if (i == 1) then
        phi1 = 7.0_8/THREE - t1 - 6.0_8*t2 - t3 + t4
        if (prn() <= phi1/phi1_max) exit
      else
        phi2 = 11.0_8/6.0_8 - t1 - THREE*t2 + HALF*t3 + t4
        if (prn() <= phi2/phi2_max) exit
      end if
    end do

    ! Compute the kinetic energy of the electron and the positron
    E_electron = (alpha*e - ONE)*MASS_ELECTRON_EV
    E_positron = (alpha*(ONE - e) - ONE)*MASS_ELECTRON_EV

    ! Sample the scattering angle of the electron. The cosine of the polar
    ! angle of the direction relative to the incident photon is sampled from
    ! p(mu) = C/(1 - beta*mu)^2 using the inverse transform method.
    beta = sqrt(E_electron*(E_electron + TWO*MASS_ELECTRON_EV)) &
         / (E_electron + MASS_ELECTRON_EV)
    rn = TWO*prn() - ONE
    mu_electron = (rn + beta)/(rn*beta + ONE)

    ! Sample the scattering angle of the positron
    beta = sqrt(E_positron*(E_positron + TWO*MASS_ELECTRON_EV)) &
         / (E_positron + MASS_ELECTRON_EV)
    rn = TWO*prn() - ONE
    mu_positron = (rn + beta)/(rn*beta + ONE)

  end subroutine pair_production

!===============================================================================
! THICK_TARGET_BREMSSTRAHLUNG
!===============================================================================

  subroutine thick_target_bremsstrahlung(p, E_lost) bind(C)
    type(Particle), intent(inout) :: p
    real(C_DOUBLE), intent(out) :: E_lost

    integer :: i, j
    integer :: i_e, i_w
    integer :: n
    integer :: n_e
    real(8) :: a
    real(8) :: f
    real(8) :: e, e_l, e_r
    real(8) :: y, y_l, y_r
    real(8) :: w, w_l, w_r
    real(8) :: p_l, p_r
    real(8) :: c, c_l, c_max
    type(BremsstrahlungData), pointer :: mat

    if (p % material == MATERIAL_VOID) return

    if (p % E < energy_cutoff(PHOTON)) return

    ! Get bremsstrahlung data for this material and particle type
    if (p % type == POSITRON) then
      mat => ttb(p % material) % positron
    else
      mat => ttb(p % material) % electron
    end if

    e = log(p % E)
    n_e = size(ttb_e_grid)

    ! Find the lower bounding index of the incident electron energy
    j = binary_search(ttb_e_grid, n_e, e)
    if (j == n_e) j = j - 1

    ! Get the interpolation bounds
    e_l = ttb_e_grid(j)
    e_r = ttb_e_grid(j+1)
    y_l = mat % yield(j)
    y_r = mat % yield(j+1)

    ! Calculate the interpolation weight w_j+1 of the bremsstrahlung energy PDF
    ! interpolated in log energy, which can be interpreted as the probability
    ! of index j+1
    f = (e - e_l)/(e_r - e_l)

    ! Get the photon number yield for the given energy using linear
    ! interpolation on a log-log scale
    y = exp(y_l + (y_r - y_l)*f)

    ! Sample number of secondary bremsstrahlung photons
    n = int(y + prn())

    E_lost = ZERO
    if (n == 0) return

    ! Sample index of the tabulated PDF in the energy grid, j or j+1
    if (prn() <= f .or. j == 1) then
      i_e = j + 1

      ! Interpolate the maximum value of the CDF at the incoming particle
      ! energy on a log-log scale
      p_l = mat % pdf(i_e-1, i_e)
      p_r = mat % pdf(i_e, i_e)
      c_l = mat % cdf(i_e-1, i_e)
      a = log(p_r/p_l)/(e_r - e_l) + ONE
      c_max = c_l + exp(e_l)*p_l/a*(exp(a*(e - e_l)) - ONE)
    else
      i_e = j

      ! Maximum value of the CDF
      c_max = mat % cdf(i_e, i_e)
    end if

    ! Sample the energies of the emitted photons
    do i = 1, n
      ! Generate a random number r and determine the index i for which
      ! cdf(i) <= r*cdf,max <= cdf(i+1)
      c = prn()*c_max
      i_w = binary_search(mat % cdf(:i_e,i_e), i_e, c)

      ! Sample the photon energy
      w_l = ttb_e_grid(i_w)
      w_r = ttb_e_grid(i_w+1)
      p_l = mat % pdf(i_w, i_e)
      p_r = mat % pdf(i_w+1, i_e)
      c_l = mat % cdf(i_w, i_e)
      a = log(p_r/p_l)/(w_r - w_l) + ONE
      w = exp(w_l)*(a*(c - c_l)/(exp(w_l)*p_l) + ONE)**(ONE/a)

      if (w > energy_cutoff(PHOTON)) then
        ! Create secondary photon
        call particle_create_secondary(p, p % coord(1) % uvw, w, PHOTON, &
             run_ce=.true._C_BOOL)
        E_lost = E_lost + w
      end if
    end do

  end subroutine thick_target_bremsstrahlung

end module photon_physics
