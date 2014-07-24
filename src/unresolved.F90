module unresolved

  use ace_header,   only: Nuclide
  use constants,    only: ZERO, ONE, SQRT_PI, TWO, THREE, PI, FOUR, C_1, &
                          K_BOLTZMANN
  use error,        only: fatal_error, warning
  use faddeeva,     only: quickw, faddeeva_w
  use global
  use output,       only: write_message
  use random_lcg,   only: prn
  use search,       only: binary_search, near_neighb

  implicit none

  logical                   :: OTF_URR       ! are we treating the URR OTF?
  character(MAX_LINE_LEN)   :: URR_METHOD    ! method for URR data treatment
  character(MAX_LINE_LEN)   :: URR_FREQUENCY ! frequency of urr xs realization
  integer, allocatable      :: urr_zaids(:)  ! ZAID #'s for URR nuclides
  type ENDFFilename
    character(MAX_LINE_LEN) :: filename      ! ENDF filename for URR nuclide
  end type ENDFFilename
  type(ENDFFilename), allocatable :: urr_endf_filenames(:) ! ENDF filename

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! RESONANCE is an object containing information about a pseudo-resonance that is
! contributing to the cross section value at an energy grid point in the URR
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  type Resonance

! TODO: account for temperature correlation?
!        need to store last resonace ensemble?

    ! partial and total xs contributions
    real(8) :: dsig_t    ! contribution of resonance to E_0 total xs
    real(8) :: dsig_n    ! contribution of resonance to E_0 scatter xs
    real(8) :: dsig_gam  ! contribution of resonance to E_0 capture xs
    real(8) :: dsig_f    ! contribution of resonance to E_0 fission xs
    real(8) :: dsig_x    ! contribution of resonance to E_0 competitive xs

    ! sampled unresolved resonance parameters
    real(8) :: D_lJ      ! sampled nuclear level spacing
    real(8) :: Gam_t     ! sampled total width
    real(8) :: Gam_n     ! sampled neutron width
    real(8) :: Gam_gam   ! sampled radiative width
    real(8) :: Gam_f     ! sampled fission width
    real(8) :: Gam_x     ! sampled competitive width

    ! resonance energy variables
    real(8) :: E_lam     ! sampled resonance energy
    real(8) :: E_lam_up  ! highest-energy level that has been added
    real(8) :: E_lam_low ! lowest-energy level that has been added
    real(8) :: E_lam_tmp ! temporary storage of resonance energy just below ES

    ! counter for the number of resonances added for a given spin sequence
    integer :: i_res

    ! type-bound procedures
    contains

      ! reset the resonance object before starting the next spin sequence
      procedure :: reset => reset_resonance

      ! sample unresolved resonance parameters
      procedure :: parameters => sample_parameters 

      ! sample level spacing
      procedure :: level_spacing => level_spacing

      ! sample channel widths
      procedure :: channel_width => channel_width

      ! interface for calculation of partial cross sections at E_0
      procedure :: xs => calc_xs

      ! calculate SLBW partial cross sections
      procedure :: slbw => slbw_xs

  end type Resonance

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! CROSSSECTION is an object containing data for a partial (or total) cross
! section
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  type CrossSection

    ! cross section value accumulator at an energy
    real(8) :: val

    ! cross section value at an energy with one less contributing resonance
    real(8) :: val_last

    ! relative error between cross section values computed with N and N+1
    ! contributing resonances
    real(8) :: rel_err

    ! type-bound procedures
    contains

      ! clear ladder partial cross sections, relative errors
      procedure :: reset => reset_xs

      ! accumulate resonance contribution to ladder partial cross section
      procedure :: resonance => accum_resonance

      ! add contribution of potential scattering cross section
      procedure :: potential => pot_scatter_xs

  end type CrossSection

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! Tabulated Chi-Squared Distribution
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

! TODO: extend tabulated values beyond Mathematica's precision
  ! tabulated chi-square distribution for 1-4 degrees of freedom with values
  ! taken to preserve integral probabilities for equiprobable bins
  real(8), dimension(20,4), parameter :: chi2=reshape((/&
     & 1.31003e-3_8,  9.19501e-3_8, 0.0250905e0_8, 0.049254e0_8, &
     & 0.0820892e0_8, 0.124169e0_8, 0.176268e0_8,  0.239417e0_8, &
     & 0.314977e0_8,  0.404749e0_8, 0.511145e0_8,  0.637461e0_8, &
     & 0.788315e0_8,  0.970419e0_8, 1.194e0_8,     1.47573e0_8,  &
     & 1.84547e0_8,   2.36522e0_8,  3.20371e0_8,   5.58201e0_8,  &
     & 0.0508548e0_8, 0.156167e0_8, 0.267335e0_8,  0.38505e0_8,  &
     & 0.510131e0_8,  0.643564e0_8, 0.786543e0_8,  0.940541e0_8, &
     & 1.1074e0_8,    1.28947e0_8,  1.48981e0_8,   1.71249e0_8,  &
     & 1.96314e0_8,   2.24984e0_8,  2.58473e0_8,   2.98744e0_8,  &
     & 3.49278e0_8,   4.17238e0_8,  5.21888e0_8,   7.99146e0_8,  &
     & 0.206832e0_8,  0.470719e0_8, 0.691933e0_8,  0.901674e0_8, &
     & 1.10868e0_8,   1.31765e0_8,  1.53193e0_8,   1.75444e0_8,  &
     & 1.98812e0_8,   2.23621e0_8,  2.50257e0_8,   2.79213e0_8,  &
     & 3.11143e0_8,   3.46967e0_8,  3.88053e0_8,   4.36586e0_8,  &
     & 4.96417e0_8,   5.75423e0_8,  6.94646e0_8,   10.0048e0_8,  &
     & 0.459462e0_8,  0.893735e0_8, 1.21753e0_8,   1.50872e0_8,  &
     & 1.78605e0_8,   2.05854e0_8,  2.33194e0_8,   2.61069e0_8,  &
     & 2.89878e0_8,   3.20032e0_8,  3.51995e0_8,   3.86331e0_8,  &
     & 4.23776e0_8,   4.65345e0_8,  5.12533e0_8,   5.67712e0_8,  &
     & 6.35044e0_8,   7.22996e0_8,  8.541e0_8,     11.8359e0_8   &
                                                              &/),(/20,4/))

contains

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! CALCULATE_URR_XS_OTF calculates xs values in the URR on-the-fly
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine calculate_urr_xs_otf(i_nuc, E)

    type(Nuclide), pointer :: nuc => null() ! nuclide object pointer
    type(Resonance)    :: res     ! resonance object
    type(CrossSection) :: sig_t   ! total xs object
    type(CrossSection) :: sig_n   ! elastic scattering xs object
    type(CrossSection) :: sig_gam ! radiative capture xs object
    type(CrossSection) :: sig_f   ! fission xs object
    type(CrossSection) :: sig_x   ! competitive inelastic scattering xs object
    integer :: i_nuc  ! nuclide index
    real(8) :: E      ! neutron energy [eV]
    integer :: i_E    ! energy grid index
    integer :: i_near ! nearest-neighbor energy index
    integer :: i_l    ! orbital quantum #
    integer :: i_J    ! total angular momentum quantum #
    integer :: i_r    ! resonance index
    real(8) :: f      ! interpolation factor

!$omp threadprivate(nuc) 

    ! Set pointer to nuclide
    nuc => nuclides(i_nuc)

    nuc % E = E
    i_E     = binary_search(nuc % ES, nuc % NE, E)
    i_near  = near_neighb(i_E, nuc % ES(i_E), nuc % ES(i_E + 1), E)
    if (nuc % INT == LINEAR_LINEAR) then
      f = (E - nuc % ES(i_E)) &
        & / (nuc % ES(i_E + 1) - nuc % ES(i_E))
    else if (nuc % INT == LOG_LOG) then
      f = log(E / nuc % ES(i_E)) / log(nuc % ES(i_E + 1) / nuc % ES(i_E))
    else
      message = 'Interpolations other than lin-lin or log-log currently not &
        & supported in OTF URR treatments'
      call fatal_error()
    end if

    ! reset xs objects
    call sig_t   % reset()
    call sig_n   % reset()
    call sig_gam % reset()
    call sig_f   % reset()
    call sig_x   % reset()

    ! loop over orbital quantum #'s
    ORBITAL_ANG_MOM_LOOP: do i_l = 1, nuc % NLS

      ! set current orbital angular momentum quantum #
      nuc % L = i_l - 1

      ! loop over total angular momentum quantum #'s
      TOTAL_ANG_MOM_LOOP: do i_J = 1, nuc % NJS(i_l)

        ! set current total angular momentum quantum #
        nuc % J = nuc % J_grid(i_l) % vals(i_J)

        ! set current partial width degrees of freedom
        nuc % AMUX = int(nuc % AMUX_grid(i_l) % vals(i_J))
        nuc % AMUN = int(nuc % AMUN_grid(i_l) % vals(i_J))
        nuc % AMUG = int(nuc % AMUG_grid(i_l) % vals(i_J))
        nuc % AMUF = int(nuc % AMUF_grid(i_l) % vals(i_J))

        ! set current mean unresolved resonance parameters
        if (nuc % INT == LINEAR_LINEAR) then
          nuc % D   = nuc % D_means(i_E, i_l)           % vals(i_J) &
            & + f * (nuc % D_means(i_E + 1, i_l)        % vals(i_J) &
            & -      nuc % D_means(i_E, i_l)            % vals(i_J))
          nuc % GN0 = nuc % Gam_n_means(i_E, i_l)       % vals(i_J) &
            & + f * (nuc % Gam_n_means(i_E + 1, i_l)    % vals(i_J) &
            & -      nuc % Gam_n_means(i_E, i_l)        % vals(i_J))
          nuc % GG  = nuc % Gam_gam_means(i_E, i_l)     % vals(i_J) &
            & + f * (nuc % Gam_gam_means(i_E + 1, i_l)  % vals(i_J) &
            & -      nuc % Gam_gam_means(i_E, i_l)      % vals(i_J))
          nuc % GF  = nuc % Gam_f_means(i_E, i_l)       % vals(i_J) &
            & + f * (nuc % Gam_f_means(i_E + 1, i_l)    % vals(i_J) &
            & -      nuc % Gam_f_means(i_E, i_l)        % vals(i_J))
          nuc % GX  = nuc % Gam_x_means(i_E, i_l)       % vals(i_J) &
            & + f * (nuc % Gam_x_means(i_E + 1, i_l)    % vals(i_J) &
            & -      nuc % Gam_x_means(i_E, i_l)        % vals(i_J))      
        else if (nuc % INT == LOG_LOG) then
          nuc % D = exp((ONE - f) * log(nuc % D_means(i_E, i_l) % vals(i_J)) &
            & + f * log(nuc % D_means(i_E + 1, i_l) % vals(i_J)))
          nuc % GN0 = exp((ONE - f) * log(nuc % Gam_n_means(i_E, i_l) % vals(i_J)) &
            & + f * log(nuc % Gam_n_means(i_E + 1, i_l) % vals(i_J)))
          nuc % GG = exp((ONE - f) * log(nuc % Gam_gam_means(i_E, i_l) % vals(i_J)) &
            & + f * log(nuc % Gam_gam_means(i_E + 1, i_l) % vals(i_J)))
          nuc % GF = exp((ONE - f) * log(nuc % Gam_f_means(i_E, i_l) % vals(i_J)) &
            & + f * log(nuc % Gam_f_means(i_E + 1, i_l) % vals(i_J)))
          nuc % GX = exp((ONE - f) * log(nuc % Gam_x_means(i_E, i_l) % vals(i_J)) &
            & + f * log(nuc % Gam_x_means(i_E + 1, i_l) % vals(i_J)))
        else
          message = 'Interpolations other than lin-lin or log-log currently not &
            & supported in OTF URR treatments'
          call fatal_error()
        end if


        ! reset the resonance object for a new spin sequence
        call res % reset(i_nuc, E)
        
        ! loop over the addition of resonances to this ladder
        RESONANCES_LOOP: do i_r = 1, nuc % n_resonances

          ! sample unresolved resonance parameters for this spin
          ! sequence, at this energy
          call res % parameters(i_nuc)

          ! set current temperature
          nuc % T = nuc % kT / K_BOLTZMANN

          ! calculate the contribution to the partial cross sections,
          ! at this energy, from an additional resonance 
          call res % xs(i_nuc)

          ! add this contribution to the accumulated partial cross
          ! section values built up from all resonances
          call sig_t   % resonance(res % dsig_t)
          call sig_n   % resonance(res % dsig_n)
          call sig_gam % resonance(res % dsig_gam)
          call sig_f   % resonance(res % dsig_f)
          call sig_x   % resonance(res % dsig_x)

        end do RESONANCES_LOOP
      end do TOTAL_ANG_MOM_LOOP
    end do ORBITAL_ANG_MOM_LOOP

    ! add potential scattering contribution
    call sig_t % potential(i_nuc)
    call sig_n % potential(i_nuc)

    ! interpret MF3 data according to ENDF self-shielding flag (LSSF)
    if (nuc % LSSF == 0) then ! MF3 contains background xs values (add to MF2
                              ! contributions)
      micro_xs(i_nuc) % total      = sig_t % val + micro_xs(i_nuc) % total
      micro_xs(i_nuc) % elastic    = sig_n % val + micro_xs(i_nuc) % elastic
      micro_xs(i_nuc) % absorption = sig_f % val + sig_gam % val &
                                 & + micro_xs(i_nuc) % absorption
      micro_xs(i_nuc) % fission    = sig_f % val + micro_xs(i_nuc) % fission
    elseif (nuc % LSSF == 1) then ! MF3 contains infinite dilute xs values
                                  ! (so the calculated xs is correct, as is)
      micro_xs(i_nuc) % total      = sig_t % val
      micro_xs(i_nuc) % elastic    = sig_n % val
      micro_xs(i_nuc) % absorption = sig_f % val + sig_gam % val
      micro_xs(i_nuc) % fission    = sig_f % val
    else
      message = 'Self-shielding flag (LSSF) not allowed - must be 0 or 1.'
      call fatal_error()
    end if

    if (micro_xs(i_nuc) % elastic < ZERO) then
      micro_xs(i_nuc) % total   = micro_xs(i_nuc) % total &
                              & + abs(micro_xs(i_nuc) % elastic)
      micro_xs(i_nuc) % elastic = ZERO
      message = 'Encountered a negative elastic scattering xs in the URR'
      call warning()
    end if

    jt = jt + ONE
    xst = xst + micro_xs(i_nuc) % total
    jf = jf + ONE
    xsf = xsf + micro_xs(i_nuc) % fission
    jn = jn + ONE
    xsn = xsn + micro_xs(i_nuc) % elastic
    jg = jg + ONE
    xsg = xsg + micro_xs(i_nuc) % absorption - micro_xs(i_nuc) % fission
    jx = jx + ONE
    xsx = xsx + micro_xs(i_nuc) % total - micro_xs(i_nuc) % elastic - micro_xs(i_nuc) % absorption
!    write(*,'(ES10.3,ES10.3,ES10.3,ES10.3,ES10.3)') xst/jt,xsf/jf,xsn/jn,&
!      & xsg/jg,xsx/jx

  end subroutine calculate_urr_xs_otf

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! RESET_RESONANCE resets the resonance object before proceeding to add
! contributions from the next spin sequence
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine reset_resonance(this, i_nuc, E)

    class(Resonance), intent(inout) :: this ! pseudo-resonance object

    type(Nuclide), pointer :: nuc => null() ! nuclide pointer
    integer                :: i_nuc         ! nuclide index
    real(8)                :: E             ! neutron energy

    nuc => nuclides(i_nuc)

    ! prepare to construct a new ladder (reset resonance counter and energies)
    this % i_res     = 0
    this % E_lam_up  = E
    this % E_lam_low = E
    this % E_lam_tmp = E

  end subroutine reset_resonance

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! SAMPLE_PARAMETERS samples unresolved resonance parameters for the next
! pseudo-resonance added to the ladder
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine sample_parameters(this, i_nuc)

    class(Resonance), intent(inout) :: this ! pseudo-resonance object

    integer :: i_nuc ! nuclide index

    ! sample unresolved resonance parameters for this resonance
    call this % level_spacing(i_nuc)
    call this % channel_width(i_nuc)

  end subroutine sample_parameters

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! LEVEL_SPACING samples the energy spacing between adjacent resonances
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine level_spacing(this, i_nuc)

    class(Resonance), intent(inout) :: this ! pseudo-resonance object

    type(Nuclide), pointer :: nuc => null() ! nuclide pointer
    integer                :: i_nuc         ! nuclide index

    nuc => nuclides(i_nuc)

    ! sample a level spacing from the Wigner distribution
    this % D_lJ = wigner_dist(nuc % D)

    ! set lowest energy (i.e. the first) resonance for this ladder well below
    ! the energy grid point such that the ladder spans a sufficient energy range
    if (this % i_res == 0) then
      this % E_lam = (nuc % E - nuc % n_resonances/2 * nuc % D) &
        & + (ONE - TWO * prn()) * this % D_lJ

    ! add subsequent resonance energies at the sampled spacing above the last
    ! resonance
    else
      this % E_lam = this % E_lam + this % D_lJ
    end if

! TODO: check why this doesn't work
!    if (1 == 2) then
!      rn = prn()
!      if (this % i_res == 0) then
!        this % E_lam     = nuc % E + (ONE - rn) * this % D_lJ
!        this % E_lam_tmp = nuc % E - rn * this % D_lJ
!        this % E_lam_up  = this % E_lam
!        this % E_lam_low = this % E_lam_tmp
!      else if (this % i_res == 1) then
!        this % E_lam     = this % E_lam_tmp
!      else
!        if (this % E_lam_up - nuc % E > nuc % E - this % E_lam_low) then
!          this % E_lam     = this % E_lam_low - this % D_lJ
!          this % E_lam_low = this % E_lam
!        else
!          this % E_lam     = this % E_lam_up + this % D_lJ
!          this % E_lam_up  = this % E_lam
!        end if
!      end if
!    end if

  end subroutine level_spacing

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! WIGNER_DIST samples the Wigner distribution for level spacings
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  function wigner_dist(D_mean) result(D_samp)

    real(8) :: D_mean ! mean level spacing
    real(8) :: D_samp ! sampled level spacing

    ! sample a level spacing by directly inverting the Wigner distribution CDF
    D_samp = D_mean * sqrt(-FOUR * log(prn()) / PI)

  end function wigner_dist

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! CHANNEL_WIDTH samples the channel partial widths
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine channel_width(this, i_nuc)

    class(Resonance), intent(inout) :: this ! pseudo-resonance object

    type(Nuclide), pointer :: nuc => null() ! nuclide object pointer
    integer                :: i_nuc         ! nuclide index
    integer                :: i_tabn        ! elastic chi-squared table index
    integer                :: i_tabg        ! capture chi-squared table index
    integer                :: i_tabf        ! fission chi-squared table index
    integer                :: i_tabx        ! competitivechi-squared table index
    real(8)                :: rho           ! derived variable
    real(8)                :: nu            ! derived variable

    nuc => nuclides(i_nuc)

! TODO: Actually sample a chi-squared distribution rather than using
! tabulated values. Look at the third Monte Carlo Sampler?

    ! sample indices into table of equiprobable chi-squared function values
    i_tabn = ceiling(prn() * 20.0_8)
    i_tabg = ceiling(prn() * 20.0_8)
    i_tabf = ceiling(prn() * 20.0_8)
    i_tabx = ceiling(prn() * 20.0_8)

    ! compute factors needed to go from the mean reduced width amplitude that is
    ! provided by ENDF for scattering to what we want - a partial width
    rho = wavenumber(nuc % awr, this % E_lam) * nuc % ac
    nu  = penetration(nuc % L, rho) / rho

    ! use the sampled tabulated chi-squared values to calculate sample widths
    ! neutron width
    this % Gam_n   = nuc % GN0 * sqrt(this % E_lam) * nu &
      & * chi2(i_tabn, nuc % AMUN)
    if (nuc % AMUN > 0) this % Gam_n = this % Gam_n / dble(nuc % AMUN)

    ! constant radiative width
    ! (many channels -> many degrees of freedom -> Dirac delta)
    this % Gam_gam = nuc % GG

    ! fission width
    if (nuc % AMUF > 0) then
      this % Gam_f = nuc % GF  * chi2(i_tabf, nuc % AMUF)
      this % Gam_f = this % Gam_f / dble(nuc % AMUF)
    end if

    ! competitive width
    this % Gam_x = nuc % GX  * chi2(i_tabx, nuc % AMUX)
    if (nuc % AMUX > 0) this % Gam_x = this % Gam_x / dble(nuc % AMUX)

    ! total width (sum of partials)
    this % Gam_t = this % Gam_n + this % Gam_gam + this % Gam_f + this % Gam_x

  end subroutine channel_width

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! CALC_XS is an interface for the calculation of partial cross sections at E_0,
! the energy that the ladder is being generated about
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine calc_xs(this, i_nuc)

    class(Resonance), intent(inout) :: this ! pseudo-resonance object

    integer :: i_nuc ! nuclide index

    ! accumulate the resonance counter by 1 for this ladder realization
    this % i_res = this % i_res + 1

    ! calculate SLBW xs contributions from an additional resonance
    call this % slbw(i_nuc)

  end subroutine calc_xs

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! SLBW_XS calculates Single-Level Breit-Wigner cross sections at an energy point
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine slbw_xs(this, i_nuc)

    class(Resonance), intent(inout) :: this ! pseudo-resonance object

    type(Nuclide), pointer :: nuc => null() ! nuclide pointer
    integer                :: i_nuc         ! nuclide index

    ! nuclear variables
    real(8) :: A        ! target mass in units of neutron mass (AWR)
    real(8) :: AP       ! scattering radius
    real(8) :: ac       ! channel radius

    ! quantum variables
    integer :: L        ! orbital quantum #
    real(8) :: J        ! total angular momentum quantum #
    real(8) :: I        ! total spin
    real(8) :: g_J      ! spin statistical factor

    ! energy variables
    real(8) :: E_n      ! neutron energy in the lab system
    real(8) :: k_n      ! center-of-mass neutron wavenumber at E_n
    real(8) :: E_lam    ! resonance energy in the lab system
    real(8) :: k_lam    ! center-of-mass neutron wavenumber at E_lam
    real(8) :: E_shift  ! shifted resonance energy in the lab system

    ! broadening variables
    real(8) :: T        ! temperature
    real(8) :: theta    ! total width / Doppler width
    real(8) :: x        ! derived variable

    ! unresolved resonance parameters
    real(8) :: Gam_t    ! sampled total width
    real(8) :: Gam_t_n  ! sampled energy-dependent total width at E_n
    real(8) :: Gam_n    ! sampled neutron width
    real(8) :: Gam_n_n  ! sampled energy-dependent neutron width at E_n
    real(8) :: Gam_gam  ! sampled capture width
    real(8) :: Gam_f    ! sampled fission width
    real(8) :: Gam_x    ! sampled competitive width
    real(8) :: sig_lam  ! peak resonance cross section

    nuc => nuclides(i_nuc)

    ! set variables
    A       = nuc % awr
    AP      = nuc % AP
    ac      = nuc % ac
    L       = nuc % L
    J       = nuc % J
    I       = nuc % SPI
    g_J     = (TWO * J + ONE) / (FOUR * I + TWO)
    Gam_n   = this % Gam_n
    E_n     = nuc % E
    k_n     = wavenumber(A, E_n)
    E_lam   = this % E_lam
    k_lam   = wavenumber(A, E_lam)
    E_shift = E_lam + (Gam_n * (shift(L, k_lam*ac) - shift(L, k_n*ac))) &
      & / (TWO * penetration(L, k_lam*ac))
    Gam_n_n = Gam_n * penetration(L, k_n*ac) / penetration(L, k_lam*ac)
    Gam_t   = this % Gam_t
    Gam_t_n = Gam_t - Gam_n + Gam_n_n
    Gam_gam = this % Gam_gam
    Gam_f   = this % Gam_f
    Gam_x   = this % Gam_x
    T       = nuc % T
    theta   = Gam_t_n / (TWO * sqrt(K_BOLTZMANN * 1.0E6_8 * T * E_n / A))
    x       = (TWO * (E_n - E_shift)) / Gam_t_n
    sig_lam = FOUR * PI / k_lam**2 * g_J * Gam_n / Gam_t

! TODO: Correct negative scattering xs values to 0 b in the library version of
!       code for use in OpenMC

! TODO: Compute competitive xs contribution correctly

    ! this particular form comes from the NJOY2012 manual
    this % dsig_n   = sig_lam * &
      & ((cos(TWO * phase_shift(L, k_n*AP)) - (ONE - Gam_n_n / Gam_t_n)) &
      & * psi(theta, x) + sin(TWO * phase_shift(L, k_n*AP)) * chi(theta, x))
    this % dsig_gam = sig_lam * Gam_gam / Gam_t_n * psi(theta, x)
    this % dsig_f   = sig_lam * Gam_f   / Gam_t_n * psi(theta, x)
    this % dsig_x   = sig_lam * Gam_x   / Gam_t_n * psi(theta, x)
    this % dsig_t = this % dsig_n   &
                & + this % dsig_gam &
                & + this % dsig_f   &
                & + this % dsig_x

! TODO: Possibly add the option of using different forms

    ! this particular form comes from the ENDF 2012 manual
!    this % dsig_n   = FOUR * PI / k_n**2 &
!      & *((TWO * dble(L) + ONE) * (sin(phase_shift(L, k_n*AP)))**2 &
!      & + g_J * (Gam_n_n**2 - TWO * Gam_n_n * Gam_t_n &
!      & * (sin(phase_shift(L, k_n*AP)))**2 + TWO * Gam_n_n * (E_n - E_shift) &
!      & * sin(TWO * phase_shift(L, k_n*AP))) &
!      & / (FOUR * (E_n - E_shift)**2 + Gam_t_n**2))

    ! this particular form comes from a previous ENDF 102 manual via Kord
!    this % dsig_n  = 2603911.0_8 / E_n * ((A + ONE) / A)**2 &
!      & * (((Gam_n_n/Gam_t_n)**2 &
!      & - TWO*Gam_n_n/Gam_t_n*(sin(phase_shift(L, k_n*AP)))**2) &
!      & * psi(theta, x) + Gam_n_n/Gam_t_n * sin(TWO * phase_shift(L, k_n*AP)) &
!      & * chi(theta, x))

  end subroutine slbw_xs

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! PENETRATION calculates hard sphere penetrability factors
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  function penetration(L, rho) result(P)

    integer, intent(in) :: L   ! current orbital quantum #
    real(8)             :: rho ! derived variable, ka
    real(8)             :: P   ! penetration factor

    ! calculate penetrability for the current orbital quantum #
    select case(L)

      case(0)
        P = rho

      case(1)
        P = rho**3 / (ONE + rho**2)

      case(2)
        P = rho**5 / (9.0_8 + THREE * rho**2 + rho**4)

      case(3)
        P = rho**7 / (225.0_8 + 45.0_8 * rho**2 + 6.0_8 * rho**4 + rho**6)

      case(4)
        P = rho**9 / (11025.0_8 + 1575.0_8 * rho**2 + 135.0_8 * rho**4 &
          & + 10.0_8 * rho**6 + rho**8)

      case default
        message = 'Orbital quantum number not allowed'
        call fatal_error()
    end select

  end function penetration

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! PHASE_SHIFT calculates hard sphere phase shifts
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  function phase_shift(L, rho) result(phi)

    integer :: L   ! current orbital quantum #
    real(8) :: rho ! derived variable, ka
    real(8) :: phi ! hard sphere phase shift

    ! calculate phase shift for the current orbital quantum #
    select case(L)

      case(0)
        phi = rho

      case(1)
        phi = rho - atan(rho)

      case(2)
        phi = rho - atan(THREE * rho / (THREE - rho**2))

      case(3)
        phi = rho - atan((15.0_8 * rho - rho**3) / (15.0_8 - 6.0_8 * rho**2))

      case(4)
        phi = rho - atan((105.0_8 * rho - 10.0_8 * rho**3) &
          & / (105.0_8 - 45.0_8 * rho**2 + rho**4))

      case default
        message = 'Orbital quantum number not allowed'
        call fatal_error()
    end select

  end function phase_shift

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! SHIFT calculates resonance energy shift factors
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  function shift(L, rho) result(S)

    integer :: L   ! current orbital quantum #
    real(8) :: rho ! derived variable, ka
    real(8) :: S   ! shift factor (for shifting the resonance energy)

    ! calculate shift factor for current orbital quantum #
    select case(L)

      case(0)
        S = ZERO

      case(1)
        S = -ONE / (ONE + rho**2)

      case(2)
        S = -(18.0_8 + THREE * rho**2) / (9.0_8 + THREE * rho**2 + rho**4)

      case(3)
        S = -(675.0_8 + 90.0_8 * rho**2 + 6.0_8 * rho**4) &
          & / (225.0_8 + 45.0_8 * rho**2 + 6.0_8 * rho**4 + rho**6)

      case(4)
        S = -(44100.0_8 + 4725.0_8 * rho**2 + &
          & 270.0_8 * rho**4 + 10.0_8 * rho**6) &
          & / (11025.0_8 + 1575.0_8 * rho**2 + 135.0_8 * rho**4 &
          & +10.0_8 * rho**6 + rho**8)

      case default
        message = 'Orbital quantum number not allowed'
        call fatal_error()
    end select
    
  end function shift

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! PSI computes a value of the psi Doppler broadening function
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  function psi(theta, x) result(psi_val)

    real(8)    :: theta   ! 
    real(8)    :: x       ! 
    real(8)    :: psi_val ! calculated value of psi
    complex(8) :: w_val   ! complex return value of the Faddeeva evaluation
    real(8)    :: relerr  ! relative error of the Faddeeva evaluation

    relerr = ZERO

! TODO: Allow option of using different Faddeeva evaluations?

    ! call S.G. Johnson's Faddeeva evaluation
    w_val = faddeeva_w(cmplx(theta * x / TWO, theta / TWO, 8), relerr)

    ! compute psi
    psi_val = SQRT_PI / TWO * theta &
      & * real(real(w_val, 8), 8)

    ! QUICKW Faddeeva evaluation from Argonne (also used in NJOY - NJOY manual)
!    psi_val = SQRT_PI / TWO * theta &
!      & * real(real(quickw(cmplx(theta * x / TWO, theta / TWO, 8)), 8), 8)

  end function psi

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! CHI computes a value of the chi Doppler broadening function
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  function chi(theta, x) result(chi_val)

    real(8)    :: theta   ! 
    real(8)    :: x       ! 
    real(8)    :: chi_val ! calculated value of chi
    complex(8) :: w_val   ! complex return value of the Faddeeva evaluation
    real(8)    :: relerr  ! relative error of the Faddeeva evaluation

    relerr = ZERO

! TODO: Allow option of using different Faddeeva evaluations?

    ! call S.G. Johnson's Faddeeva evaluation
    w_val = faddeeva_w(cmplx(theta * x / TWO, theta / TWO, 8), relerr)

    ! compute chi
    chi_val = SQRT_PI / TWO * theta &
      & * real(dimag(w_val), 8)

    ! QUICKW Faddeeva evaluation from Argonne (also used in NJOY - NJOY manual)
!    chi_val = SQRT_PI / TWO * theta &
!      & * real(dimag(quickw(cmplx(theta * x / TWO, theta / TWO, 8))), 8)

  end function chi

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! WAVENUMBER computes a center-of-mass reference frame wavenumber
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  function wavenumber(A, E) result(k_val)

    real(8) :: A     ! atomic weight ratio
    real(8) :: E     ! evaluation energy
    real(8) :: k_val ! computed wavenumber

    ! compute center-of-mass neutron wavenumber evaluated at some energy
    k_val = C_1 * A / (A + ONE) * sqrt(E)

  end function wavenumber

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! ACCUM_RESONANCE accumulates the contribution to the ladder partial cross
! section due to the addition of a resonance
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine accum_resonance(this, sig)

    class(CrossSection), intent(inout) :: this ! cross section object

    real(8) :: sig ! contribution to xs from the new resonance

    ! set the current xs value as the last value
    this % val_last = this % val

    ! add xs contribution from a new resonance to the xs value at the current E
    this % val = this % val + sig

    ! compute the magnitude of the relative error between current and last xs
    if (this % val > ZERO) &
      this % rel_err  = abs((this % val_last - this % val) / this % val)

  end subroutine accum_resonance

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! POT_SCATTER_XS adds the contribution of potential scattering to the elastic
! and total cross sections
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine pot_scatter_xs(this, i_nuc)

    class(CrossSection), intent(inout) :: this ! cross section object

    type(Nuclide), pointer :: nuc => null() ! nuclide pointer
    integer :: i_nuc   ! index in nuclides
    real(8) :: sig_pot ! potential scattering cross section
    real(8) :: k_n     ! center-of-mass neutron wavenumber
    integer :: NLS     ! # of orbital quantum #'s
    real(8) :: E_n     ! neutron energy in the lab system
    real(8) :: A       ! target mass in units of neutron mass
    integer :: i_L     ! orbital quantum # index
    real(8) :: AP      ! scattering radius

    ! set nuclide variables
    nuc => nuclides(i_nuc)
    E_n = nuc % E
    A   = nuc % awr
    AP  = nuc % AP
    NLS = nuc % NLS

    ! compute neutron COM wavenumber
    k_n = wavenumber(A, E_n)

    ! compute potential scattering xs by adding contribution from each l-wave
    sig_pot = ZERO
    do i_L = 0, NLS - 1
      sig_pot = sig_pot &
        & + FOUR * PI / k_n**2 * (TWO * dble(i_L) + ONE) &
        & * (sin(phase_shift(i_L, k_n*AP)))**2
    end do

    ! add the potential scattering xs to this xs
    this % val = this % val + sig_pot

  end subroutine pot_scatter_xs

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! RESET_XS resets the xs values and relative error for a cross section object
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine reset_xs(this)

    class(CrossSection), intent(inout) :: this ! cross section object

    ! clear accumulated xs values for this realization and set error to INF
    this % val      = ZERO
    this % val_last = ZERO
    this % rel_err  = ONE

  end subroutine reset_xs

end module unresolved
