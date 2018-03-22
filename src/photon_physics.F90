module photon_physics

  use algorithm,       only: binary_search
  use constants
  use particle_header, only: Particle
  use photon_header,   only: PhotonInteraction, Bremsstrahlung, &
                             compton_profile_pz, ttb_e_grid, ttb_k_grid, ttb
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

  subroutine compton_scatter(el, alpha, alpha_out, mu, use_doppler)
    type(PhotonInteraction), intent(in) :: el
    real(8), intent(in)  :: alpha
    real(8), intent(out) :: alpha_out
    real(8), intent(out) :: mu
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
      x = MASS_ELECTRON/PLANCK_C*alpha*sqrt(HALF*(ONE - mu))

      ! Calculate S(x, Z) and S(x_max, Z)
      form_factor_x = el % incoherent_form_factor % evaluate(x)
      if (form_factor_xmax == ZERO) then
        form_factor_xmax = el % incoherent_form_factor % evaluate(&
             MASS_ELECTRON/PLANCK_C*alpha)
      end if

      ! Perform rejection on form factor
      if (prn() < form_factor_x / form_factor_xmax) then
        if (use_doppler_) then
          call compton_doppler(el, alpha, mu, e_out)
          alpha_out = e_out/MASS_ELECTRON
        end if
        exit
      end if
    end do

  end subroutine compton_scatter

!===============================================================================
! COMPTON_DOPPLER
!===============================================================================

  subroutine compton_doppler(el, alpha, mu, e_out)
    type(PhotonInteraction), intent(in) :: el
    real(8), intent(in)  :: alpha
    real(8), intent(in)  :: mu
    real(8), intent(out) :: e_out

    integer :: i, i_shell
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
      e = alpha*MASS_ELECTRON
      if (e < e_b) then
        e_out = alpha/(1 + alpha*(1 - mu))*MASS_ELECTRON
        exit
      end if

      pz_max = -FINE_STRUCTURE*(e_b - (e - e_b)*alpha*(ONE - mu)) / &
           sqrt(TWO*e*(e - e_b)*(ONE - mu) + e_b**2)
      if (pz_max < ZERO) then
        e_out = alpha/(1 + alpha*(1 - mu))*MASS_ELECTRON
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
        e_out = alpha/(1 + alpha*(1 - mu))*MASS_ELECTRON
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
      x2_max = (MASS_ELECTRON/PLANCK_C*alpha)**2

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
      call p % create_secondary(uvw, E, PHOTON, run_ce=.true.)
      return
    end if

    ! Sample transition
    rn = prn()
    c = ZERO
    do i_transition = 1, elm % shells(i_shell) % n_transitions - 1
      c = c + elm % shells(i_shell) % &
           transition_probability(i_transition + 1)
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
      call p % create_secondary(uvw, E, ELECTRON, run_ce=.true.)

      ! Fill hole left by emitted auger electron
      i_hole = elm % shell_dict % get(secondary)
      call atomic_relaxation(p, elm, i_hole)
    else
      ! Radiative transition -- get X-ray energy

      ! Create fluorescent photon
      call p % create_secondary(uvw, E, PHOTON, run_ce=.true.)

    end if

    ! Fill hole created by electron transitioning to the photoelectron hole
    i_hole = elm % shell_dict % get(primary)
    call atomic_relaxation(p, elm, i_hole)

  end subroutine atomic_relaxation

!===============================================================================
! THICK_TARGET_BREMSSTRAHLUNG
!===============================================================================

  subroutine thick_target_bremsstrahlung(p, E_lost)
    type(Particle), intent(inout) :: p
    real(8),        intent(inout) :: E_lost

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
    type(Bremsstrahlung), pointer :: mat

    if (p % E < energy_cutoff(PHOTON)) return

    ! Get bremsstrahlung data for this material
    mat => ttb(p % material)

    e = log(p % E)
    n_e = size(ttb_e_grid)

    ! Find the lower bounding index of the incident electron energy
    j = binary_search(ttb_e_grid, n_e, e)

    ! Get the interpolation bounds
    e_l = ttb_e_grid(j)
    e_r = ttb_e_grid(j+1)
    y_l = mat % yield(j)
    y_r = mat % yield(j+1)

    ! Calculate the interpolation weight w_j+1 of the bremsstrahlung energy PDF
    ! interpolated in log energy, which can be interpreted as the probability
    ! of index j+1
    f = (e - e_l) / (e_r - e_l)

    ! Get the photon number yield for the given energy using linear
    ! interpolation on a log-log scale
    y = exp(y_l + (y_r - y_l) * f)

    ! Sample number of secondary bremsstrahlung photons
    n = int(y + prn())

    E_lost = ZERO
    if (n == 0) return

    ! Sample index of the tabulated PDF in the energy grid, j or j+1
    if (prn() > f) then
      i_e = j

      ! Maximum value of the CDF
      c_max = mat % cdf(i_e, i_e)
    else
      i_e = j + 1

      ! Interpolate the maximum value of the CDF at the incoming particle
      ! energy on a log-log scale
      p_l = mat % pdf(i_e, i_e-1)
      p_r = mat % pdf(i_e, i_e)
      c_l = mat % cdf(i_e, i_e-1)
write(*,*) "p_r: ", p_r, "p_l: ", p_l, "p_r/p_l: ", p_r/p_l
write(*,*) "e_r: ", e_r, "e_l: ", e_l, "e_r/e_l: ", e_r/e_l
      a = (log(p_r/p_l)) / (e_r - e_l) + ONE
      c_max = c_l + (exp(e_l) * p_l)/a * (exp(a*(e - e_l)) - ONE)
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
      a = (log(p_r/p_l)) / (w_r - w_l) + ONE
      w = exp(w_l) * (a*(c - c_l)/(exp(w_l) * p_l) + ONE)**(ONE/a)


      if (w < energy_cutoff(PHOTON)) cycle

      ! Create secondary photon
      call p % create_secondary(p % coord(1) % uvw, w, PHOTON, run_ce=.true.)
      E_lost = E_lost + w
    end do

  end subroutine thick_target_bremsstrahlung

end module photon_physics
