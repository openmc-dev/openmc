module unresolved

  use ace_header,    only: Nuclide, Reaction
  use error,         only: fatal_error, warning
  use faddeeva,      only: quickw, faddeeva_w
  use fission,       only: nu_total
  use global
  use output,        only: write_message
  use random_lcg,    only: prn
  use search,        only: binary_search
  use vector_header, only: Vector

  implicit none

  logical :: urr_method    ! new urr method?
  logical :: competitive   ! competitve inelastic scatter xs resonance structure?
  logical :: urr_pointwise ! pointwise cross section calculation?
  character(80), allocatable :: urr_endf_filenames(:) ! ENDF filename list
  character(80)              :: urr_formalism         ! URR formalism
  character(80)              :: urr_frequency         ! freq of realizations
  integer, allocatable :: urr_zaids(:)    ! ZAID's for URR nuclides
  integer, allocatable :: n_resonances(:) ! # URR resonances for each l-wave
  integer :: n_otf_urr_xs      ! number of nuclides to calc otf urr xs for
  integer :: n_avg_urr_xs      ! number of nuclides to calc average urr xs for
  integer :: n_urr_method      ! number of nuclides to treat with a new method
  integer :: n_s_wave          ! number of contributing s-wave resonances
  integer :: n_p_wave          ! number of contributing p-wave resonances
  integer :: n_d_wave          ! number of contributing d-wave resonances
  integer :: n_f_wave          ! number of contributing f-wave resonances
  integer :: urr_avg_batches   ! min number of batches for inf dil xs calc
  integer :: urr_avg_histories ! number of histories for inf dil xs calc
  real(8) :: urr_avg_tol       ! max rel error for inf dil xs calc termination

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
      procedure :: reset_resonance => resonance_reset

      ! sample unresolved resonance parameters
      procedure :: sample_parameters => parameters_sample

      ! sample level spacing
      procedure :: level_spacing => spacing_level

      ! sample channel widths
      procedure :: channel_width => width_channel

      ! interface for calculation of partial cross sections at E_0
      procedure :: calc_xs => xs_calc

      ! calculate SLBW partial cross sections
      procedure :: slbw_xs => xs_slbw

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

    ! history-based sum of cross section values
    real(8) :: xs_sum_tmp

    ! batch-based sum of cross section values
    real(8) :: xs_sum

    ! sum of squares of batch-based cross section value sums
    real(8) :: xs_sum2

    ! mean of batch-based cross section values
    real(8) :: xs_mean

    ! standard error of the mean cross section values
    real(8) :: xs_sem

    ! relative uncertainty
    real(8) :: rel_unc

    ! type-bound procedures
    contains

      ! accumulate resonance contribution to ladder partial cross section
      procedure :: accum_resonance => resonance_accum

      ! add contribution of potential scattering cross section
      procedure :: potential_xs => xs_potential

      ! accumulate values for a single history
      procedure :: accum_history => history_accum

      ! clear values for a single history
      procedure :: flush_history => history_flush

      ! accumulate values for a single batch
      procedure :: accum_batch => batch_accum

      ! clear values for a single batch
      procedure :: flush_batch => batch_flush

      ! calculate batch statistics
      procedure :: calc_stats => stats_calc

      ! clear batch statistics
      procedure :: flush_stats => stats_flush

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

  real(8), parameter :: Eid(18) = (/&
    2.000E+04_8,&
    2.300E+04_8,&
    2.600E+04_8,&
    3.000E+04_8,&
    3.500E+04_8,&
    4.000E+04_8,&
    4.500E+04_8,&
    4.509E+04_8,&
    5.000E+04_8,&
    5.500E+04_8,&
    6.000E+04_8,&
    7.000E+04_8,&
    8.000E+04_8,&
    9.000E+04_8,&
    1.000E+05_8,&
    1.200E+05_8,&
    1.400E+05_8,&
    1.490E+05_8/)

  real(8), parameter :: xsidn(18) = (/&
    1.385E+01_8,&
    1.369E+01_8,&
    1.357E+01_8,&
    1.342E+01_8,&
    1.328E+01_8,&
    1.315E+01_8,&
    1.306E+01_8,&
    1.305E+01_8,&
    1.293E+01_8,&
    1.279E+01_8,&
    1.266E+01_8,&
    1.244E+01_8,&
    1.224E+01_8,&
    1.207E+01_8,&
    1.191E+01_8,&
    1.163E+01_8,&
    1.139E+01_8,&
    1.129E+01_8/)

  real(8), parameter :: xsidf(18) = (/&
    ZERO,&
    ZERO,&
    ZERO,&
    ZERO,&
    ZERO,&
    ZERO,&
    ZERO,&
    ZERO,&
    ZERO,&
    ZERO,&
    ZERO,&
    ZERO,&
    ZERO,&
    ZERO,&
    ZERO,&
    ZERO,&
    ZERO,&
    ZERO/)

  real(8), parameter :: xsidg(18) = (/&
    5.296E-01_8,&
    4.957E-01_8,&
    4.681E-01_8,&
    4.373E-01_8,&
    4.065E-01_8,&
    3.826E-01_8,&
    3.631E-01_8,&
    3.623E-01_8,&
    3.194E-01_8,&
    2.888E-01_8,&
    2.641E-01_8,&
    2.286E-01_8,&
    2.043E-01_8,&
    1.871E-01_8,&
    1.739E-01_8,&
    1.566E-01_8,&
    1.459E-01_8,&
    1.423E-01_8/)

  real(8), parameter :: xsidx(18) = (/&
    0.000E+00_8,&
    0.000E+00_8,&
    0.000E+00_8,&
    0.000E+00_8,&
    0.000E+00_8,&
    0.000E+00_8,&
    0.000E+00_8,&
    0.000E+00_8,&
    6.315E-02_8,&
    1.304E-01_8,&
    1.936E-01_8,&
    3.016E-01_8,&
    3.890E-01_8,&
    4.615E-01_8,&
    5.229E-01_8,&
    6.207E-01_8,&
    6.959E-01_8,&
    7.244E-01_8/)

contains

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! RESONANCE_ENSEMBLE generates a single realization of a resonance ensemble in
! the URR
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine resonance_ensemble()

    type(Nuclide), pointer, save :: nuc => null() ! nuclide object pointer
    type(Resonance) :: res ! resonance object
    integer :: i_nuc       ! nuclide index
    integer :: i_l         ! orbital quantum number index
    integer :: i_J         ! total angular momentum quantum #
    integer :: i_E         ! tabulated URR parameters energy index
    integer :: i_res       ! resonance counter
    integer :: n_res       ! number of l-wave resonances to include
    integer :: n_above_urr ! number of resonances abover upper URR energy
    real(8) :: E_res ! current resonance (lab) energy (e.g. E_lam)
    real(8) :: m     ! energy interpolation factor

!$omp threadprivate(nuc) 

    ! loop over all nuclides
    NUCLIDE_LOOP: do i_nuc = 1, n_nuclides_total

      nuc => nuclides(i_nuc)

      ! skip nuclides without an OTF URR treatment
      if (.not. nuc % otf_urr_xs) cycle

      micro_xs(i_nuc) % use_ptable = .true.

      ! allocate a realization of a URR resonance ensemble
      call nuc % alloc_ensemble

      allocate(n_resonances(nuc % NLS(nuc % i_urr)))

      ! loop over orbital angular momenta
      ORBITAL_ANG_MOM_LOOP: do i_l = 1, nuc % NLS(nuc % i_urr)

        ! set current orbital angular momentum quantum #
        nuc % L = i_l - 1

        ! get the number of contributing l-wave resonances for this l
        n_res = l_wave_resonances(nuc % L)

        ! loop over total angular momenta
        TOTAL_ANG_MOM_LOOP: do i_J = 1, nuc % NJS(i_l)

          ! set current total angular momentum quantum #
          nuc % J = nuc % AJ(i_l) % data(i_J)

          ! set current partial width degrees of freedom
          nuc % AMUX = int(nuc % DOFX(i_l) % data(i_J))
          nuc % AMUN = int(nuc % DOFN(i_l) % data(i_J))
          nuc % AMUG = int(nuc % DOFG(i_l) % data(i_J))
          nuc % AMUF = int(nuc % DOFF(i_l) % data(i_J))

          ! set energy of the lowest-lying contributing resonance
          E_res = nuc % EL(nuc % i_urr) &
            & - (n_res/2 + ONE - TWO * prn()) * nuc % D_mean(i_l) % data(i_J) % data(1)

          ! resonance index
          i_res = 1

          n_above_urr = 0

          RESONANCE_LOOP: do while(n_above_urr < n_res/2)!E_res < nuc % EH(nuc % i_urr) &
!            & + n_res/2 * nuc % D_mean(i_l) % data(i_J) % data(nuc % NE))

            ! compute interpolation factor
            if (E_res < nuc % ES(1)) then
              i_E = 1
            else if (E_res > nuc % ES(nuc % NE)) then
              i_E = nuc % NE - 1
            else
              i_E = binary_search(nuc % ES, nuc % NE, E_res)
            end if
            m = interp_factor(E_res, nuc % ES(i_E), nuc % ES(i_E + 1), nuc % INT)

            ! set current mean unresolved resonance parameters
            nuc % D   = interpolator(m, &
              & nuc % D_mean(i_l) % data(i_J) % data(i_E), &
              & nuc % D_mean(i_l) % data(i_J) % data(i_E + 1), nuc % INT)
            nuc % GN0 = interpolator(m, &
              & nuc % GN0_mean(i_l) % data(i_J) % data(i_E), &
              & nuc % GN0_mean(i_l) % data(i_J) % data(i_E + 1), nuc % INT)
            nuc % GG  = interpolator(m, &
              & nuc % GG_mean(i_l) % data(i_J) % data(i_E), &
              & nuc % GG_mean(i_l) % data(i_J) % data(i_E + 1), nuc % INT)

            ! TODO: add in catch here for when threshold occurs between tabulated pts
            if (nuc % GF_mean(i_l) % data(i_J) % data(i_E) /= ZERO &
              & .and. nuc % GF_mean(i_l) % data(i_J) % data(i_E + 1) /= ZERO) then
              nuc % GF  = interpolator(m, &
                & nuc % GF_mean(i_l) % data(i_J) % data(i_E), &
                & nuc % GF_mean(i_l) % data(i_J) % data(i_E + 1), nuc % INT)
            else
              nuc % GF = ZERO
            end if

            ! TODO: add in catch here for when threshold occurs between tabulated pts
            if (nuc % GX_mean(i_l) % data(i_J) % data(i_E) /= ZERO &
              & .and. nuc % GX_mean(i_l) % data(i_J) % data(i_E + 1) /= ZERO) then
              nuc % GX  = interpolator(m, &
                & nuc % GX_mean(i_l) % data(i_J) % data(i_E), &
                & nuc % GX_mean(i_l) % data(i_J) % data(i_E + 1), nuc % INT)
            else
              nuc % GX = ZERO
            end if

            ! sample unresolved resonance parameters for this spin
            ! sequence, at this energy
            res % E_lam = E_res
            call res % channel_width(i_nuc)
            nuc % urr_resonances(i_res, i_l) % E_lam(i_J) = E_res
            nuc % urr_resonances(i_res, i_l) % GN(i_J)    = res % Gam_n
            nuc % urr_resonances(i_res, i_l) % GG(i_J)    = res % Gam_gam
            nuc % urr_resonances(i_res, i_l) % GF(i_J)    = res % Gam_f
            nuc % urr_resonances(i_res, i_l) % GX(i_J)    = res % Gam_x
            nuc % urr_resonances(i_res, i_l) % GT(i_J)    = res % Gam_t

            ! add an additional resonance
            i_res = i_res + 1
            E_res = E_res + wigner_dist(nuc % D)

            if (E_res > nuc % EH(nuc % i_urr)) n_above_urr = n_above_urr + 1

          end do RESONANCE_LOOP
        end do TOTAL_ANG_MOM_LOOP

        n_resonances(i_l) = i_res - 1

      end do ORBITAL_ANG_MOM_LOOP
    end do NUCLIDE_LOOP

  end subroutine resonance_ensemble

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! POINTWISE_URR generates pointwise energy-cross section data in the URR
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine pointwise_urr()

    type(Nuclide), pointer, save :: nuc => null() ! nuclide object pointer
    type(Vector), allocatable :: i_last(:) ! last i_low for (l,J)
    type(Resonance) :: res ! resonance object
    type(CrossSection) :: sig_n   ! elastic scattering xs object
    type(CrossSection) :: sig_gam ! radiative capture xs object
    type(CrossSection) :: sig_f   ! fission xs object
    type(CrossSection) :: sig_x   ! competitive inelastic scattering xs object
    type(CrossSection) :: sig_t   ! total xs object
    integer :: i_nuc  ! nuclide index
    integer :: i_E    ! URR resonance parameters energy grid index
    integer :: i_l    ! orbital quantum #
    integer :: i_J    ! total angular momentum quantum #
    integer :: n_res  ! number of resonances to include for a given l-wave
    integer :: i_low  ! index of lowest-lying resonance
    integer :: i_res  ! resonance counter
    integer :: i_grid ! cross section energy grid index
    integer :: n_pts  ! cross section energy grid point counter
    real(8) :: m        ! URR resonance parameters interpolation factor
    real(8) :: f        ! cross section energy grid interpolation factor
    real(8) :: xsidnval ! averaged elastic cross section
    real(8) :: xsidfval ! averaged fission cross section
    real(8) :: xsidgval ! averaged capture cross section
    real(8) :: xsidxval ! averaged competitive inelastic cross section

!$omp threadprivate(nuc)

    ! loop over all nuclides
    NUCLIDE_LOOP: do i_nuc = 1, n_nuclides_total

      nuc => nuclides(i_nuc)

      ! skip nuclides without an OTF URR treatment
      if (urr_pointwise .and. nuc % otf_urr_xs) then
        nuc % point_urr_xs = .true.
      else
        nuc % point_urr_xs = .false.
        cycle
      end if

      micro_xs(i_nuc) % use_ptable = .true.

      ! set current temperature
      nuc % T = nuc % kT / K_BOLTZMANN

      ! set current energy
      nuc % E = nuc % EL(nuc % i_urr)

      ! allocate pointwise URR cross sections
      call nuc % alloc_pointwise_tmp()

      n_pts = 0

      allocate(i_last(nuc % NLS(nuc % i_urr)))
      do i_l = 1, nuc % NLS(nuc % i_urr)
        allocate(i_last(i_l) % data(nuc % NJS(i_l)))
        i_last(i_l) % data = ONE
      end do

      ENERGY_LOOP: do

        n_pts = n_pts + 1

        if (nuc % E > nuc % ES(nuc % NE)) nuc % E = nuc % ES(nuc % NE)
        i_E = binary_search(nuc % ES, nuc % NE, nuc % E)
        m = interp_factor(nuc % E, nuc % ES(i_E), nuc % ES(i_E + 1), nuc % INT)

        ! reset xs objects
        call sig_n   % flush_history()
        call sig_gam % flush_history()
        call sig_f   % flush_history()
        call sig_x   % flush_history()
        call sig_t   % flush_history()

        ! loop over orbital quantum #'s
        ORBITAL_ANG_MOM_LOOP: do i_l = 1, nuc % NLS(nuc % i_urr)

          ! set current orbital angular momentum quantum #
          nuc % L = i_l - 1

          ! get the number of contributing l-wave resonances for this l
          n_res = l_wave_resonances(nuc % L)

          ! loop over total angular momentum quantum #'s
          TOTAL_ANG_MOM_LOOP: do i_J = 1, nuc % NJS(i_l)

            ! set current total angular momentum quantum #
            nuc % J = nuc % AJ(i_l) % data(i_J)

            ! compute statistical spin factor
            nuc % g_J = (TWO * nuc % J + ONE) &
              & / (FOUR * nuc % SPI(nuc % i_urr) + TWO)

            ! reset the resonance object for a new spin sequence
            call res % reset_resonance(i_nuc)

            ! find the nearest lower resonance
            i_low = int(i_last(i_l) % data(i_J))

            do while(nuc % urr_resonances(i_low, i_l) % E_lam(i_J) < nuc % E)
              i_low = i_low + 1
            end do
            i_low = i_low - 1
            i_last(i_l) % data(i_J) = dble(i_low)

            ! loop over the addition of resonances to this ladder
            RESONANCES_LOOP: do i_res=max(i_low-n_res/2,1),max(i_low-n_res/2,1)+20

              ! set URR resonance parameters
              res % E_lam   = nuc % urr_resonances(i_res, i_l) % E_lam(i_J)
              res % Gam_n   = nuc % urr_resonances(i_res, i_l) % GN(i_J)
              res % Gam_gam = nuc % urr_resonances(i_res, i_l) % GG(i_J)
              res % Gam_f   = nuc % urr_resonances(i_res, i_l) % GF(i_J)
              res % Gam_x   = nuc % urr_resonances(i_res, i_l) % GX(i_J)
              res % Gam_t   = nuc % urr_resonances(i_res, i_l) % GT(i_J)

              ! calculate the contribution to the partial cross sections,
              ! at this energy, from an additional resonance 
              call res % calc_xs(i_nuc)

              ! add this contribution to the accumulated partial cross
              ! section values built up from all resonances
! TODO: move sig_t outside of loop
              call sig_n   % accum_resonance(res % dsig_n)
              call sig_gam % accum_resonance(res % dsig_gam)
              call sig_f   % accum_resonance(res % dsig_f)
              call sig_x   % accum_resonance(res % dsig_x)
              call sig_t   % accum_resonance(res % dsig_t)

            end do RESONANCES_LOOP
          end do TOTAL_ANG_MOM_LOOP
        end do ORBITAL_ANG_MOM_LOOP

        ! add potential scattering contribution
        call sig_n % potential_xs(i_nuc)
        call sig_t % potential_xs(i_nuc)

        if (nuc % E < 1.0e6_8 * nuc % energy(1)) then
          i_grid = 1
        elseif (nuc % E > 1.0e6_8 * nuc % energy(nuc % n_grid)) then
          i_grid = nuc % n_grid - 1
        else
          i_grid = binary_search(1.0e6_8 * nuc % energy, nuc % n_grid, nuc % E)
        end if

        ! check for rare case where two energy points are the same                                                                                                         
        if (nuc % energy(i_grid) == nuc % energy(i_grid+1)) i_grid = i_grid + 1

        ! calculate xs energy grid interpolation factor
        f = interp_factor(nuc % E, 1.0e6_8 * nuc % energy(i_grid), &
          & 1.0e6_8 * nuc % energy(i_grid + 1), LINEAR_LINEAR)

        ! add energy point to grid
        nuc % urr_energy_tmp(n_pts) = nuc % E

        ! calculate evaluator-supplied backgrounds at the current energy
        ! elastic scattering xs
        nuc % urr_elastic_tmp(n_pts)&
          & = interpolator(f, nuc % elastic(i_grid),&
          & nuc % elastic(i_grid + 1),&
          & LINEAR_LINEAR)

        ! radiative capture xs
        nuc % urr_capture_tmp(n_pts)&
          & = interpolator(f, nuc % absorption(i_grid) - nuc % fission(i_grid),&
          & nuc % absorption(i_grid + 1) - nuc % fission(i_grid + 1),&
          & LINEAR_LINEAR)

        ! fission xs
        nuc % urr_fission_tmp(n_pts)&
          & = interpolator(f, nuc % fission(i_grid),&
          & nuc % fission(i_grid + 1),&
          & LINEAR_LINEAR)

        ! competitive first level inelastic scattering xs
        nuc % urr_inelastic_tmp(n_pts)&
          & = interpolator(f, nuc % total(i_grid)&
          &                 - nuc % absorption(i_grid)&
          &                 - nuc % elastic(i_grid),&
          &                   nuc % total(i_grid + 1)&
          &                 - nuc % absorption(i_grid + 1)&
          &                 - nuc % elastic(i_grid + 1), LINEAR_LINEAR)

        ! total xs
        nuc % urr_total_tmp(n_pts)&
          & = interpolator(f, nuc % total(i_grid),&
          & nuc % total(i_grid + 1),&
          & LINEAR_LINEAR)

        ! interpret MF3 data according to ENDF self-shielding factor flag (LSSF):
        ! MF3 contains background xs values (add to MF2 resonance contributions)
        if (nuc % LSSF == 0) then
          call fatal_error('LSSF = 0 not yet supported')

          ! elastic scattering xs
          nuc % urr_elastic_tmp(n_pts) = nuc % urr_elastic_tmp(n_pts) &
            & + sig_n % val

          ! radiative capture xs
          nuc % urr_capture_tmp(n_pts) = nuc % urr_capture_tmp(n_pts) &
            & + sig_gam % val

          ! fission xs
          nuc % urr_fission_tmp(n_pts) = nuc % urr_fission_tmp(n_pts) &
            & + sig_f % val

          ! competitive first level inelastic scattering xs
          nuc % urr_inelastic_tmp(n_pts) = nuc % urr_inelastic_tmp(n_pts) &
            & + sig_x % val

          ! total xs
          nuc % urr_total_tmp(n_pts) = nuc % urr_total_tmp(n_pts) &
            & + sig_t % val

        ! multipy the self-shielding factors by the average (infinite-dilute)
        ! cross sections
        elseif (nuc % LSSF == 1) then

! TODO: LOG_LOG?
          xsidnval = interpolator(m, xsidn(i_E), xsidn(i_E + 1), LINEAR_LINEAR)
          xsidfval = interpolator(m, xsidf(i_E), xsidf(i_E + 1), LINEAR_LINEAR)
          xsidgval = interpolator(m, xsidg(i_E), xsidg(i_E + 1), LINEAR_LINEAR)
          xsidxval = interpolator(m, xsidx(i_E), xsidx(i_E + 1), LINEAR_LINEAR)

          ! competitive xs
          if (xsidxval > ZERO) then
            if (competitive) then
              ! self-shielded treatment of competitive inelastic cross section
              nuc % urr_inelastic_tmp(n_pts) = sig_x % val / xsidxval &
                & * nuc % urr_inelastic_tmp(n_pts)
            else
              ! infinite-dilute treatment of competitive inelastic cross section
              nuc % urr_inelastic_tmp(n_pts) = nuc % urr_inelastic_tmp(n_pts)
            end if
          else
            ! use background competitive inelastic cross section, as is
            nuc % urr_inelastic_tmp(n_pts) = nuc % urr_inelastic_tmp(n_pts)
          end if

          ! elastic scattering xs
          nuc % urr_elastic_tmp(n_pts) = sig_n % val / xsidnval &
            & * nuc % urr_elastic_tmp(n_pts)

          ! set negative SLBW elastic xs to zero
          if (nuc % urr_elastic_tmp(n_pts) < ZERO) nuc % urr_elastic_tmp(n_pts) = ZERO

          ! radiative capture xs
          nuc % urr_capture_tmp(n_pts) = sig_gam % val / xsidgval &
            & * nuc % urr_capture_tmp(n_pts)

          ! fission xs
          if (xsidfval > ZERO) then
            nuc % urr_fission_tmp(n_pts) = sig_f % val / xsidfval &
              & * nuc % urr_fission_tmp(n_pts)
          else
            nuc % urr_fission_tmp(n_pts) = nuc % urr_fission_tmp(n_pts)
          end if

          nuc % urr_total_tmp(n_pts) = nuc % urr_elastic_tmp(n_pts)&
            &                        + nuc % urr_capture_tmp(n_pts)&
            &                        + nuc % urr_fission_tmp(n_pts)&
            &                        + nuc % urr_inelastic_tmp(n_pts)

        else
          call fatal_error('Self-shielding flag (LSSF) not allowed - must be 0 or 1.')
        end if

        nuc % E = nuc % E + nuc % urr_dE
        if (nuc % E >= nuc % EH(nuc % i_urr) + nuc % urr_dE) exit

      end do ENERGY_LOOP

      call nuc % alloc_pointwise(n_pts)
      nuc % urr_energy    = nuc % urr_energy_tmp(1:n_pts)
      nuc % urr_elastic   = nuc % urr_elastic_tmp(1:n_pts)
      nuc % urr_capture   = nuc % urr_capture_tmp(1:n_pts)
      nuc % urr_fission   = nuc % urr_fission_tmp(1:n_pts)
      nuc % urr_inelastic = nuc % urr_inelastic_tmp(1:n_pts)
      nuc % urr_total     = nuc % urr_total_tmp(1:n_pts)
      call nuc % dealloc_pointwise_tmp()

    end do NUCLIDE_LOOP

  end subroutine pointwise_urr

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! CALCULATE_AVG_URR_XS computes the infinite-dilute partial cross sections in
! the URR
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine calculate_avg_urr_xs()

    type(Nuclide), pointer, save :: nuc => null() ! nuclide object pointer
    type(Resonance)    :: res   ! resonance object
    type(CrossSection) :: sig_n ! elastic scattering xs object
    type(CrossSection) :: sig_g ! radiative capture xs object
    type(CrossSection) :: sig_f ! fission xs object
    type(CrossSection) :: sig_x ! competitive inelastic scattering xs object
    integer :: i_nuc ! nuclide index
    integer :: i_E   ! energy grid index
    integer :: i_b   ! batch index
    integer :: i_h   ! history index
    integer :: i_l   ! orbital quantum #
    integer :: i_J   ! total angular momentum quantum #
    integer :: i_r   ! resonance index
    integer :: n_res ! number of resonances to include for a given l-wave
    real(8) :: E ! neutron lab energy [eV]

    ! loop over all nuclides
    NUCLIDE_LOOP: do i_nuc = 1, n_nuclides_total

      nuc => nuclides(i_nuc)

      ! skip nuclides without an OTF URR treatment
      if (.not. nuc % avg_urr_xs) cycle

      micro_xs(i_nuc) % use_ptable = .true.

      ! allocate infinite-dilute cross section objects
      call nuc % alloc_avg_urr()

      ! loop over tabulated ENDF energies
      ENERGY_LOOP: do i_E = 1, nuc % NE

        nuc % E = nuc % ES(i_E)
        E = nuc % E

        ! reset accumulator of statistics
        call sig_n % flush_stats()
        call sig_f % flush_stats()
        call sig_g % flush_stats()
        call sig_x % flush_stats()

        i_b = 0

        ! loop over batches until convergence
        BATCH_LOOP: do

          i_b = i_b + 1

          call sig_n % flush_batch
          call sig_f % flush_batch
          call sig_g % flush_batch
          call sig_x % flush_batch

          ! loop over realizations
          HISTORY_LOOP: do i_h = 1, urr_avg_histories

            ! reset accumulator of histories
            call sig_n % flush_history()
            call sig_f % flush_history()
            call sig_g % flush_history()
            call sig_x % flush_history()

            ! loop over orbital quantum #'s
            ORBITAL_ANG_MOM_LOOP: do i_l = 1, nuc % NLS(nuc % i_urr)

              ! set current orbital angular momentum quantum #
              nuc % L = i_l - 1

              ! get the number of contributing l-wave resonances for this l
              n_res = l_wave_resonances(nuc % L)

              ! loop over total angular momentum quantum #'s
              TOTAL_ANG_MOM_LOOP: do i_J = 1, nuc % NJS(i_l)

                ! set current total angular momentum quantum #
                nuc % J = nuc % AJ(i_l) % data(i_J)

                ! compute statistical spin factor
                nuc % g_J = (TWO * nuc % J + ONE) / (FOUR * nuc % SPI(nuc % i_urr) + TWO)

                ! set current partial width degrees of freedom
                nuc % AMUX = int(nuc % DOFX(i_l) % data(i_J))
                nuc % AMUN = int(nuc % DOFN(i_l) % data(i_J))
                nuc % AMUG = int(nuc % DOFG(i_l) % data(i_J))
                nuc % AMUF = int(nuc % DOFF(i_l) % data(i_J))

                ! set current mean unresolved resonance parameters
                nuc % D   = nuc % D_mean(i_l)   % data(i_J) % data(i_E)
                nuc % GN0 = nuc % GN0_mean(i_l) % data(i_J) % data(i_E)
                nuc % GG  = nuc % GG_mean(i_l)  % data(i_J) % data(i_E)
                nuc % GF  = nuc % GF_mean(i_l)  % data(i_J) % data(i_E)
                nuc % GX  = nuc % GX_mean(i_l)  % data(i_J) % data(i_E)

                ! reset the resonance object for a new spin sequence
                call res % reset_resonance(i_nuc)
        
                ! loop over the addition of resonances to this ladder
                RESONANCES_LOOP: do i_r = 1, n_res

                  ! sample unresolved resonance parameters for this spin
                  ! sequence, at this energy
                  call res % sample_parameters(i_nuc)

                  ! set current temperature
                  nuc % T = nuc % kT / K_BOLTZMANN

                  ! calculate the contribution to the partial cross sections,
                  ! at this energy, from an additional resonance 
                  call res % calc_xs(i_nuc)

                  ! add this contribution to the accumulated partial cross
                  ! section values built up from all resonances
                  call sig_n   % accum_resonance(res % dsig_n)
                  call sig_g   % accum_resonance(res % dsig_gam)
                  call sig_f   % accum_resonance(res % dsig_f)
                  call sig_x   % accum_resonance(res % dsig_x)

                end do RESONANCES_LOOP
              end do TOTAL_ANG_MOM_LOOP
            end do ORBITAL_ANG_MOM_LOOP

            ! add potential scattering contribution
            call sig_n % potential_xs(i_nuc)

            ! set negative SLBW elastic xs to zero
            if (sig_n % val < ZERO) then
              sig_n % val = ZERO
            end if

            ! accumulate the result of this history
            call sig_n % accum_history()
            call sig_f % accum_history()
            call sig_g % accum_history()
            call sig_x % accum_history()

          end do HISTORY_LOOP

          ! accumulate the result of this batch
          call sig_n % accum_batch()
          call sig_f % accum_batch()
          call sig_g % accum_batch()
          call sig_x % accum_batch()

          ! calculate statistics for this batch
          call sig_n % calc_stats(i_b)
          call sig_f % calc_stats(i_b)
          call sig_g % calc_stats(i_b)
          call sig_x % calc_stats(i_b)

! TODO: format avg urr xs output
          if (i_b > urr_avg_batches .and. max(sig_n % rel_unc, sig_f % rel_unc, &
            & sig_g % rel_unc, sig_x % rel_unc) < urr_avg_tol) then
            if (1==1) then
              write(*, '(I5, ES10.3, ES10.3, ES10.3, ES10.3, ES10.3, ES10.3, ES10.3, ES10.3, ES10.3)') &
                & i_b, E, sig_n % xs_mean, sig_n % xs_sem, &
                &         sig_f % xs_mean, sig_f % xs_sem, &
                &         sig_g % xs_mean, sig_g % xs_sem, &
                &         sig_x % xs_mean, sig_x % xs_sem
            end if
            exit
          end if

        end do BATCH_LOOP

        ! set infinite-dilute xs values at this energy to converged means
        nuc % avg_urr_n(i_E) = sig_n % xs_mean
        nuc % avg_urr_f(i_E) = sig_f % xs_mean
        nuc % avg_urr_g(i_E) = sig_g % xs_mean
        nuc % avg_urr_x(i_E) = sig_x % xs_mean

      end do ENERGY_LOOP
    end do NUCLIDE_LOOP

  end subroutine calculate_avg_urr_xs

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! CALCULATE_URR_XS_OTF calculates xs values in the URR on-the-fly
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine calculate_urr_xs_otf(i_nuc, E)

    type(Nuclide), pointer, save :: nuc => null() ! nuclide object pointer
    type(Resonance) :: res ! resonance object
    type(CrossSection) :: sig_n   ! elastic scattering xs object
    type(CrossSection) :: sig_gam ! radiative capture xs object
    type(CrossSection) :: sig_f   ! fission xs object
    type(CrossSection) :: sig_x   ! competitive inelastic scattering xs object
    type(CrossSection) :: sig_t   ! total xs object
    integer :: i_nuc  ! nuclide index
    integer :: i_E    ! energy grid index
    integer :: i_l    ! orbital quantum #
    integer :: i_J    ! total angular momentum quantum #
    integer :: i_r    ! resonance index
    integer :: n_res ! number of resonances to include for a given l-wave
    integer :: i_energy
    real(8) :: xsidnval ! infinite-dilute n xs value from NJOY's MC^2 quadrature
    real(8) :: xsidfval ! infinite-dilute f xs value from NJOY's MC^2 quadrature
    real(8) :: xsidgval ! infinite-dilute g xs value from NJOY's MC^2 quadrature
    real(8) :: xsidxval ! infinite-dilute x xs value from NJOY's MC^2 quadrature
    real(8) :: inelastic_val ! competitive inelastic scattering cross section
    real(8) :: capture_val   ! radiative capture cross section
    real(8) :: E      ! neutron energy [eV]
    real(8) :: m      ! energy interpolation factor

!$omp threadprivate(nuc) 

    ! Set pointer to nuclide
    nuc => nuclides(i_nuc)

    micro_xs(i_nuc) % use_ptable = .true.

    if (nuc % point_urr_xs) then
      i_E = binary_search(nuc % urr_energy, size(nuc % urr_energy), E)
      m = interp_factor(E, nuc % urr_energy(i_E), nuc % urr_energy(i_E + 1), &
        & nuc % INT)
      micro_xs(i_nuc) % elastic    &
        & = interpolator(m, nuc % urr_elastic(i_E), nuc % urr_elastic(i_E + 1),&
        & nuc % INT)
      micro_xs(i_nuc) % fission    &
        & = interpolator(m, nuc % urr_fission(i_E), nuc % urr_fission(i_E + 1),&
        & nuc % INT)
      micro_xs(i_nuc) % absorption &
        & = interpolator(m, nuc % urr_capture(i_E) + nuc % urr_fission(i_E),&
        & nuc % urr_capture(i_E + 1) + nuc % urr_fission(i_E + 1), nuc % INT)
      if ((nuc % urr_inelastic(i_E) >= ZERO) .and. (nuc % urr_inelastic(i_E + 1) >= ZERO)) then
        inelastic_val &
          & = interpolator(m, nuc % urr_inelastic(i_E), nuc % urr_inelastic(i_E + 1),&
          & nuc % INT)
      else
        inelastic_val = ZERO
      end if
      micro_xs(i_nuc) % total = micro_xs(i_nuc) % elastic &
                            & + micro_xs(i_nuc) % absorption &
                            & + inelastic_val

      ! Determine nu-fission cross section
      if (nuc % fissionable) then
        micro_xs(i_nuc) % nu_fission = nu_total(nuc, E/1.0e6_8) * &
          micro_xs(i_nuc) % fission
      end if
      return
    end if

    ! set current temperature
    nuc % T = nuc % kT / K_BOLTZMANN

    ! set current energy and interpolation factor
    nuc % E = E
    i_E = binary_search(nuc % ES, nuc % NE, nuc % E)
    m = interp_factor(nuc % E, nuc % ES(i_E), nuc % ES(i_E + 1), nuc % INT)

    ! reset xs objects
    call sig_t   % flush_history()
    call sig_n   % flush_history()
    call sig_gam % flush_history()
    call sig_f   % flush_history()
    call sig_x   % flush_history()

    ! loop over orbital quantum #'s
    ORBITAL_ANG_MOM_LOOP: do i_l = 1, nuc % NLS(nuc % i_urr)

      ! set current orbital angular momentum quantum #
      nuc % L = i_l - 1

      ! get the number of contributing l-wave resonances for this l
      n_res = l_wave_resonances(nuc % L)

      ! loop over total angular momentum quantum #'s
      TOTAL_ANG_MOM_LOOP: do i_J = 1, nuc % NJS(i_l)

        ! set current total angular momentum quantum #
        nuc % J = nuc % AJ(i_l) % data(i_J)

        ! compute statistical spin factor
        nuc % g_J = (TWO * nuc % J + ONE) / (FOUR * nuc % SPI(nuc % i_urr) + TWO)

        ! set current partial width degrees of freedom
        nuc % AMUX = int(nuc % DOFX(i_l) % data(i_J))
        nuc % AMUN = int(nuc % DOFN(i_l) % data(i_J))
        nuc % AMUG = int(nuc % DOFG(i_l) % data(i_J))
        nuc % AMUF = int(nuc % DOFF(i_l) % data(i_J))

        ! set current mean unresolved resonance parameters
        nuc % D   = interpolator(m, &
          & nuc % D_mean(i_l) % data(i_J) % data(i_E), &
          & nuc % D_mean(i_l) % data(i_J) % data(i_E + 1), nuc % INT)

        nuc % GN0 = interpolator(m, &
          & nuc % GN0_mean(i_l) % data(i_J) % data(i_E), &
          & nuc % GN0_mean(i_l) % data(i_J) % data(i_E + 1), nuc % INT)

        nuc % GG  = interpolator(m, &
          & nuc % GG_mean(i_l) % data(i_J) % data(i_E), &
          & nuc % GG_mean(i_l) % data(i_J) % data(i_E + 1), nuc % INT)

        ! TODO: add in catch here for when threshold occurs between tabulated pts
        if (nuc % GF_mean(i_l) % data(i_J) % data(i_E) /= ZERO &
          & .and. nuc % GF_mean(i_l) % data(i_J) % data(i_E + 1) /= ZERO) then
          nuc % GF  = interpolator(m, &
            & nuc % GF_mean(i_l) % data(i_J) % data(i_E), &
            & nuc % GF_mean(i_l) % data(i_J) % data(i_E + 1), nuc % INT)
        else
          nuc % GF = ZERO
        end if

        ! TODO: add in catch here for when threshold occurs between tabulated pts
        if (nuc % GX_mean(i_l) % data(i_J) % data(i_E) /= ZERO &
          & .and. nuc % GX_mean(i_l) % data(i_J) % data(i_E + 1) /= ZERO) then
          nuc % GX  = interpolator(m, &
            & nuc % GX_mean(i_l) % data(i_J) % data(i_E), &
            & nuc % GX_mean(i_l) % data(i_J) % data(i_E + 1), nuc % INT)
        else
          nuc % GX = ZERO
        end if

        ! reset the resonance object for a new spin sequence
        call res % reset_resonance(i_nuc)

        ! loop over the addition of resonances to this ladder
        RESONANCES_LOOP: do i_r = 1, n_res

          ! sample unresolved resonance parameters for this spin
          ! sequence, at this energy
          call res % sample_parameters(i_nuc)

          ! calculate the contribution to the partial cross sections,
          ! at this energy, from an additional resonance 
          call res % calc_xs(i_nuc)

          ! add this contribution to the accumulated partial cross
          ! section values built up from all resonances
! TODO: move sig_t outside of loop
          call sig_t   % accum_resonance(res % dsig_t)
          call sig_n   % accum_resonance(res % dsig_n)
          call sig_gam % accum_resonance(res % dsig_gam)
          call sig_f   % accum_resonance(res % dsig_f)
          call sig_x   % accum_resonance(res % dsig_x)

        end do RESONANCES_LOOP
      end do TOTAL_ANG_MOM_LOOP
    end do ORBITAL_ANG_MOM_LOOP

    ! add potential scattering contribution
    call sig_t % potential_xs(i_nuc)
    call sig_n % potential_xs(i_nuc)

    ! determine energy table
    i_energy = 1
    do
      if (E <= Eid(i_energy + 1)) exit
      i_energy = i_energy + 1
      if (i_energy >= size(Eid)) then
        i_energy = size(Eid) - 1
        exit
      end if
    end do

    ! interpret MF3 data according to ENDF self-shielding factor flag (LSSF):
    ! MF3 contains background xs values (add to MF2 resonance contributions)
    if (nuc % LSSF == 0) then
      call fatal_error('LSSF = 0 not yet supported')
      micro_xs(i_nuc) % total      = sig_t % val + micro_xs(i_nuc) % total
      micro_xs(i_nuc) % elastic    = sig_n % val + micro_xs(i_nuc) % elastic
      micro_xs(i_nuc) % absorption = sig_f % val + sig_gam % val &
                                 & + micro_xs(i_nuc) % absorption
      micro_xs(i_nuc) % fission    = sig_f % val + micro_xs(i_nuc) % fission
    
    ! MF3 contains evaluator-supplied infinite dilute xs values that we multipy
    ! the self-shielding factors computed from MF2 by
    elseif (nuc % LSSF == 1) then

      ! TODO: LOG_LOG?
      xsidnval = xsidn(i_energy) &
        & + (E - Eid(i_energy)) / (Eid(i_energy+1) - Eid(i_energy)) &
        & * (xsidn(i_energy+1) - xsidn(i_energy))
      xsidfval = xsidf(i_energy) &
        & + (E - Eid(i_energy)) / (Eid(i_energy+1) - Eid(i_energy)) &
        & * (xsidf(i_energy+1) - xsidf(i_energy))
      xsidgval = xsidg(i_energy) &
        & + (E - Eid(i_energy)) / (Eid(i_energy+1) - Eid(i_energy)) &
        & * (xsidg(i_energy+1) - xsidg(i_energy))
      xsidxval = xsidx(i_energy) &
        & + (E - Eid(i_energy)) / (Eid(i_energy+1) - Eid(i_energy)) &
        & * (xsidx(i_energy+1) - xsidx(i_energy))

      if (xsidxval > ZERO) then

        if (competitive) then
          ! self-shielded treatment of competitive inelastic cross section
          inelastic_val = sig_x % val / xsidxval &
            & * (micro_xs(i_nuc) % total &
            & - micro_xs(i_nuc) % absorption &
            & - micro_xs(i_nuc) % elastic)
        else
          ! infinite-dilute treatment of competitive inelastic cross section
          inelastic_val = micro_xs(i_nuc) % total &
            & - micro_xs(i_nuc) % absorption &
            & - micro_xs(i_nuc) % elastic
        end if

      else
        inelastic_val = micro_xs(i_nuc) % total &
          & - micro_xs(i_nuc) % absorption &
          & - micro_xs(i_nuc) % elastic
      end if

      micro_xs(i_nuc) % elastic = sig_n % val / xsidnval &
        & * micro_xs(i_nuc) % elastic

      ! set negative SLBW elastic xs to zero
      if (micro_xs(i_nuc) % elastic < ZERO) then
        micro_xs(i_nuc) % elastic = ZERO
      end if

      capture_val = sig_gam % val / xsidgval &
        & * (micro_xs(i_nuc) % absorption - micro_xs(i_nuc) % fission)

      if (xsidfval > ZERO) then
        micro_xs(i_nuc) % fission = sig_f % val / xsidfval &
          & * micro_xs(i_nuc) % fission
      else
        micro_xs(i_nuc) % fission = micro_xs(i_nuc) % fission
      end if

      micro_xs(i_nuc) % absorption = micro_xs(i_nuc) % fission + capture_val

      micro_xs(i_nuc) % total = micro_xs(i_nuc) % elastic &
        & + micro_xs(i_nuc) % absorption &
        & + inelastic_val

    else
      call fatal_error('Self-shielding flag (LSSF) not allowed - must be 0 or 1.')
    end if

    ! Determine nu-fission cross section
    if (nuc % fissionable) then
      micro_xs(i_nuc) % nu_fission = nu_total(nuc, E/1.0e6_8) * &
           micro_xs(i_nuc) % fission
    end if

  end subroutine calculate_urr_xs_otf

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! RESONANCE_RESET resets the resonance object before proceeding to add
! contributions from the next spin sequence
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine resonance_reset(this, i_nuc)

    class(Resonance), intent(inout) :: this ! pseudo-resonance object

    type(Nuclide), pointer :: nuc => null() ! nuclide pointer
    integer                :: i_nuc         ! nuclide index

    nuc => nuclides(i_nuc)

    ! prepare to construct a new ladder (reset resonance counter and energies)
    this % i_res     = 0
    this % E_lam_up  = nuc % E
    this % E_lam_low = nuc % E
    this % E_lam_tmp = nuc % E

  end subroutine resonance_reset

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! PARAMETERS_SAMPLE samples unresolved resonance parameters for the next
! pseudo-resonance added to the ladder
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine parameters_sample(this, i_nuc)

    class(Resonance), intent(inout) :: this ! pseudo-resonance object

    integer :: i_nuc ! nuclide index

    ! sample unresolved resonance parameters for this resonance
    call this % level_spacing(i_nuc)
    call this % channel_width(i_nuc)

  end subroutine parameters_sample

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! SPACING_LEVEL samples the energy spacing between adjacent resonances
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine spacing_level(this, i_nuc)

    class(Resonance), intent(inout) :: this ! pseudo-resonance object

    type(Nuclide), pointer :: nuc => null() ! nuclide pointer
    integer                :: i_nuc         ! nuclide index
    integer :: n_res ! number of resonances to include for a given l-wave

    nuc => nuclides(i_nuc)

    n_res = l_wave_resonances(nuc % L)

    ! sample a level spacing from the Wigner distribution
    this % D_lJ = wigner_dist(nuc % D)

    ! set lowest energy (i.e. the first) resonance for this ladder well below
    ! the energy grid point such that the ladder spans a sufficient energy range
    if (this % i_res == 0) then
      this % E_lam = (nuc % E - n_res/2 * nuc % D) &
        & + (ONE - TWO * prn()) * this % D_lJ

    ! add subsequent resonance energies at the sampled spacing above the last
    ! resonance
    else
      this % E_lam = this % E_lam + this % D_lJ
    end if

  end subroutine spacing_level

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! L_WAVE_RESONANCES determines the number of resonances to include from the lth
! wave when computing a URR cross section on-the-fly
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  function l_wave_resonances(L) result(n_res)

    integer :: L     ! orbital quantum number
    integer :: n_res ! number of resonances to include for this l

    select case(L)
    case(0)
      n_res = n_s_wave
    case(1)
      n_res = n_p_wave
    case(2)
      n_res = n_d_wave
    case(3)
      n_res = n_f_wave
    case default
      call fatal_error('Only s-, p-, d-, and f-wave resonances are supported in ENDF-6')
    end select

  end function l_wave_resonances

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! WIGNER_DIST samples the Wigner distribution for level spacings
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  function wigner_dist(D_avg) result(D_samp)

    real(8) :: D_avg  ! mean level spacing
    real(8) :: D_samp ! sampled level spacing

    ! sample a level spacing by directly inverting the Wigner distribution CDF
    D_samp = D_avg * sqrt(-FOUR * log(prn()) / PI)

  end function wigner_dist

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! WIDTH_CHANNEL samples the channel partial widths
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine width_channel(this, i_nuc)

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
    i_tabf = ceiling(prn() * 20.0_8)
    i_tabg = ceiling(prn() * 20.0_8)
    i_tabx = ceiling(prn() * 20.0_8)

    ! compute factors needed to go from the mean reduced width amplitude that is
    ! provided by ENDF for scattering to what we want - a partial width
    rho = wavenumber(nuc % awr, this % E_lam) * nuc % ac(nuc % i_urr)
    nu  = penetration(nuc % L, rho) / rho

    ! use the sampled tabulated chi-squared values to calculate sample widths
    ! neutron width
    if (nuc % AMUN > 0) then
      this % Gam_n = nuc % GN0 * sqrt(this % E_lam) * nu &
        & * chi2(i_tabn, nuc % AMUN)
      this % Gam_n = this % Gam_n / dble(nuc % AMUN)
    else
      this % Gam_n = ZERO
    end if

    ! fission width
    if (nuc % AMUF > 0) then
      this % Gam_f = nuc % GF  * chi2(i_tabf, nuc % AMUF)
      this % Gam_f = this % Gam_f / dble(nuc % AMUF)
    else 
      this % Gam_f = ZERO
    end if

    ! constant radiative width
    ! (many channels -> many degrees of freedom -> Dirac delta)
    this % Gam_gam = nuc % GG

    ! competitive width
    if (nuc % AMUX > 0) then
      this % Gam_x = nuc % GX  * chi2(i_tabx, nuc % AMUX)
      this % Gam_x = this % Gam_x / dble(nuc % AMUX)
    else
      this % Gam_x = ZERO
    end if

    ! total width (sum of partials)
    this % Gam_t = this % Gam_n + this % Gam_f + this % Gam_gam + this % Gam_x

  end subroutine width_channel

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! XS_CALC is an interface for the calculation of partial cross sections at E_0,
! the energy that the ladder is being generated about
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine xs_calc(this, i_nuc)

    class(Resonance), intent(inout) :: this ! pseudo-resonance object

    integer :: i_nuc ! nuclide index

    ! accumulate the resonance counter by 1 for this ladder realization
    this % i_res = this % i_res + 1

    ! calculate SLBW xs contributions from an additional resonance
    call this % slbw_xs(i_nuc)

  end subroutine xs_calc

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! XS_SLBW calculates Single-Level Breit-Wigner cross sections at an energy point
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine xs_slbw(this, i_nuc)

    class(Resonance), intent(inout) :: this ! pseudo-resonance object

    type(Nuclide), pointer :: nuc => null() ! nuclide pointer
    integer                :: i_nuc         ! nuclide index

    ! energy variables
    real(8) :: k_n      ! center-of-mass neutron wavenumber at E_n
    real(8) :: k_lam    ! center-of-mass neutron wavenumber at E_lam
    real(8) :: E_shift  ! shifted resonance energy in the lab system

    ! broadening variables
    real(8) :: theta    ! total width / Doppler width
    real(8) :: x        ! derived variable

    ! unresolved resonance parameters
    real(8) :: Gam_t_n  ! sampled energy-dependent total width at E_n
    real(8) :: Gam_n_n  ! sampled energy-dependent neutron width at E_n
    real(8) :: sig_lam  ! peak resonance cross section
    real(8) :: sig_lam_Gam_t_n_psi

    nuc => nuclides(i_nuc)

    ! set variables
    k_n     = wavenumber(nuc % awr, nuc % E)
    k_lam   = wavenumber(nuc % awr, this % E_lam)
    E_shift = this % E_lam &
      & + (this % Gam_n * (shift(nuc % L, k_lam*nuc % ac(nuc % i_urr)) &
      & - shift(nuc % L, k_n*nuc % ac(nuc % i_urr)))) &
      & / (TWO * penetration(nuc % L, k_lam*nuc % ac(nuc % i_urr)))
    Gam_n_n = this % Gam_n &
      & * penetration(nuc % L, k_n*nuc % ac(nuc % i_urr)) &
      & / penetration(nuc % L, k_lam*nuc % ac(nuc % i_urr))
    Gam_t_n = this % Gam_t - this % Gam_n + Gam_n_n
    theta   = Gam_t_n &
      & / (TWO * sqrt(K_BOLTZMANN * 1.0E6_8 * nuc % T * nuc % E / nuc % awr))
    x       = (TWO * (nuc % E - E_shift)) / Gam_t_n
    sig_lam = FOUR * PI / (k_lam * k_lam) * nuc % g_J * this % Gam_n / this % Gam_t

! TODO: Correct negative scattering xs values to 0 b in the library version of
!       code for use in OpenMC

! TODO: Compute competitive xs contribution correctly

    ! this particular form comes from the NJOY2012 manual
    if (Gam_n_n > ZERO) then
      this % dsig_n = sig_lam * &
        & ((cos(TWO * phase_shift(nuc % L, k_n*nuc % AP(nuc % i_urr))) - (ONE - Gam_n_n / Gam_t_n)) &
        & * psi(theta, x) + sin(TWO * phase_shift(nuc % L, k_n*nuc % AP(nuc % i_urr))) * chi(theta, x))
    else
      call fatal_error('Encountered a non-positive elastic scattering width in the URR')
    end if

    sig_lam_Gam_t_n_psi = sig_lam * psi(theta, x) / Gam_t_n

    if (this % Gam_gam > ZERO) then
      this % dsig_gam = sig_lam_Gam_t_n_psi * this % Gam_gam
    else
      this % dsig_gam = ZERO
    end if

    if (this % Gam_f > ZERO) then
      this % dsig_f   = sig_lam_Gam_t_n_psi * this % Gam_f
    else
      this % dsig_f   = ZERO
    end if

    if (this % Gam_x > ZERO) then
      this % dsig_x   = sig_lam_Gam_t_n_psi * this % Gam_x
    else
      this % dsig_x   = ZERO
    end if

    this % dsig_t = this % dsig_n   &
                & + this % dsig_gam &
                & + this % dsig_f   &
                & + this % dsig_x

  end subroutine xs_slbw

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! PENETRATION calculates hard sphere penetrability factors
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  function penetration(L, rho) result(P)

    integer, intent(in) :: L    ! current orbital quantum #
    real(8)             :: rho  ! derived variable, ka
    real(8)             :: rho2 ! rho**2
    real(8)             :: P    ! penetration factor

    ! pre-compute exponentiation
    rho2 = rho * rho

    ! calculate penetrability for the current orbital quantum #
    select case(L)

      case(0)
        P = rho

      case(1)
        P = rho * rho2 / (ONE + rho2)

      case(2)
        P = rho * rho2 * rho2 / (9.0_8 + rho2 * (THREE + rho2))

      case(3)
        P = rho * rho2 * rho2 * rho2 &
          & / (225.0_8 + rho2 * (45.0_8 + rho2 * (6.0_8 + rho2)))

      case(4)
        P = rho * rho2 * rho2 * rho2 * rho2 &
          & / (11025.0_8 + rho2 * (1575.0_8 + rho2 * (135.0_8 + rho2 * (10.0_8 + rho2))))

      case default
        call fatal_error('Orbital quantum number not allowed')
    end select

  end function penetration

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! PHASE_SHIFT calculates hard sphere phase shifts
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  function phase_shift(L, rho) result(phi)

    integer :: L    ! current orbital quantum #
    real(8) :: rho  ! derived variable, ka
    real(8) :: rho2 ! rho**2
    real(8) :: rho3 ! rho**3
    real(8) :: rho4 ! rho**4
    real(8) :: phi  ! hard sphere phase shift

    ! pre-compute exponentiations
    rho2 = rho * rho
    rho3 = rho * rho * rho
    rho4 = rho * rho * rho * rho

    ! calculate phase shift for the current orbital quantum #
    select case(L)

      case(0)
        phi = rho

      case(1)
        phi = rho - atan(rho)

      case(2)
        phi = rho - atan(THREE * rho / (THREE - rho2))

      case(3)
        phi = rho - atan((15.0_8 * rho - rho3) / (15.0_8 - 6.0_8 * rho2))

      case(4)
        phi = rho - atan((105.0_8 * rho - 10.0_8 * rho3) &
          & / (105.0_8 - 45.0_8 * rho2 + rho4))

      case default
        call fatal_error('Orbital quantum number not allowed')
    end select

  end function phase_shift

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! SHIFT calculates resonance energy shift factors
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  function shift(L, rho) result(S)

    integer :: L    ! current orbital quantum #
    real(8) :: rho  ! derived variable, ka
    real(8) :: rho2 ! rho**2
    real(8) :: S    ! shift factor (for shifting the resonance energy)

    ! pre-compute exponentiation
    rho2 = rho * rho

    ! calculate shift factor for current orbital quantum #
    select case(L)

      case(0)
        S = ZERO

      case(1)
        S = -ONE / (ONE + rho2)

      case(2)
        S = -(18.0_8 + THREE * rho2) / (9.0_8 + rho2 * (THREE + rho2))

      case(3)
        S = -(675.0_8 + rho2 * (90.0_8 + 6.0_8 * rho2)) &
          & / (225.0_8 + rho2 * (45.0_8 + rho2 * (6.0_8 + rho2)))

      case(4)
        S = -(44100.0_8 + rho2 * (4725.0_8 + rho2 * (270.0_8 + 10.0_8 * rho2))) &
          & / (11025.0_8 + rho2 * (1575.0_8 + rho2 * (135.0_8 + rho2 * (10.0_8 + rho2))))

      case default
        call fatal_error('Orbital quantum number not allowed')

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

    ! evaluate the W (Faddeeva) function
    select case (w_eval)

    ! call S.G. Johnson's Faddeeva evaluation
    case (MIT_W)
      relerr = 1.0e-1
      w_val = faddeeva_w(cmplx(theta * x / TWO, theta / TWO, 8), relerr)
      psi_val = SQRT_PI / TWO * theta &
        & * real(real(w_val, 8), 8)

    ! QUICKW Faddeeva evaluation from Argonne (also used in NJOY - NJOY manual)
    case (QUICK_W)
      psi_val = SQRT_PI / TWO * theta &
        & * real(real(quickw(cmplx(theta * x / TWO, theta / TWO, 8)), 8), 8)

    case default
      call fatal_error('Unrecognized W function evaluation method')
    end select

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

    ! evaluate the W (Faddeeva) function
    select case (w_eval)

    ! call S.G. Johnson's Faddeeva evaluation
    case (MIT_W)

      relerr = 1.0e-1
      w_val = faddeeva_w(cmplx(theta * x / TWO, theta / TWO, 8), relerr)
      chi_val = SQRT_PI / TWO * theta &
        & * real(aimag(w_val), 8)

    ! QUICKW Faddeeva evaluation from Argonne (also used in NJOY - NJOY manual)
    case (QUICK_W)

      chi_val = SQRT_PI / TWO * theta &
        & * real(aimag(quickw(cmplx(theta * x / TWO, theta / TWO, 8))), 8)

    case default

      call fatal_error('Unrecognized W function evaluation method')

    end select

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
! RESONANCE_ACCUM accumulates the contribution to the ladder partial cross
! section due to the addition of a resonance
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine resonance_accum(this, sig)

    class(CrossSection), intent(inout) :: this ! cross section object

    real(8) :: sig ! contribution to xs from the new resonance

    ! set the current xs value as the last value
    this % val_last = this % val

    ! add xs contribution from a new resonance to the xs value at the current E
    this % val = this % val + sig

  end subroutine resonance_accum

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! XS_POTENTIAL adds the contribution of potential scattering to the elastic
! and total cross sections
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine xs_potential(this, i_nuc)

    class(CrossSection), intent(inout) :: this ! cross section object

    type(Nuclide), pointer :: nuc => null() ! nuclide pointer
    integer :: i_nuc   ! index in nuclides
    real(8) :: sig_pot ! potential scattering cross section
    real(8) :: k_n     ! center-of-mass neutron wavenumber
    integer :: i_l     ! orbital quantum # index

    ! set nuclide variables
    nuc => nuclides(i_nuc)

    ! compute neutron COM wavenumber
    k_n = wavenumber(nuc % awr, nuc % E)

    ! compute potential scattering xs by adding contribution from each l-wave
    sig_pot = ZERO
    do i_l = 0, nuc % NLS(nuc % i_urr) - 1
      sig_pot = sig_pot &
        & + FOUR * PI / (k_n * k_n) * (TWO * dble(i_l) + ONE) &
        & * (sin(phase_shift(i_l, k_n*nuc % AP(nuc % i_urr)))) &
        & * (sin(phase_shift(i_l, k_n*nuc % AP(nuc % i_urr))))
    end do

    ! add the potential scattering xs to this xs
    this % val = this % val + sig_pot

  end subroutine xs_potential

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! HISTORY_ACCUM adds the single-history xs realization to the single-batch
! accumulator
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine history_accum(this)

    class(CrossSection), intent(inout) :: this ! cross section object

    ! accumulate history xs value for this realization
    this % xs_sum_tmp = this % xs_sum_tmp + this % val

  end subroutine history_accum

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! HISTORY_FLUSH flushes the single-history xs realization
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine history_flush(this)

    class(CrossSection), intent(inout) :: this ! cross section object

    ! clear accumulated xs values for this realization
    this % val      = ZERO
    this % val_last = ZERO

  end subroutine history_flush

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! BATCH_ACCUM adds the single-batch xs realization to the overall accumulator
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine batch_accum(this)

    class(CrossSection), intent(inout) :: this ! cross section object

    ! accumulate batch-wise mean
    this % xs_sum = this % xs_sum + this % xs_sum_tmp / dble(urr_avg_histories)

    ! accumulate squared batch-wise mean
    this % xs_sum2 = this % xs_sum2 &
      & + (this % xs_sum_tmp / dble(urr_avg_histories)) &
      & * (this % xs_sum_tmp / dble(urr_avg_histories))

  end subroutine batch_accum

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! BATCH_FLUSH clears the single-batch sum of xs realizations for the batch
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine batch_flush(this)

    class(CrossSection), intent(inout) :: this ! cross section object

    ! clear accumulated xs values for this batch
    this % xs_sum_tmp = ZERO

  end subroutine batch_flush

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! STATS_CALC computes batch-based means and standard errors of those means
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine stats_calc(this, i_bat_int)

    class(CrossSection), intent(inout) :: this ! cross section object
    integer :: i_bat_int ! current batch index
    real(8) :: i_bat ! current batch index for calcs

    i_bat = dble(i_bat_int)

    ! compute mean xs value
    this % xs_mean  = this % xs_sum / i_bat

    if (this % xs_mean /= ZERO .and. i_bat_int > 1) then
      ! compute standard error of mean xs value
      this % xs_sem   = sqrt((ONE / (i_bat - ONE)) &
        & * (this % xs_sum2 / i_bat &
        & - (this % xs_mean) * (this % xs_mean)))

      ! compute relative uncertainty in xs value
      this % rel_unc  = this % xs_sem / this % xs_mean
    end if

  end subroutine stats_calc

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! STATS_FLUSH clears the batch statistics and accumulators
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine stats_flush(this)

    class(CrossSection), intent(inout) :: this ! cross section object

    ! clear mean, SEM, and overall accumulators
    this % xs_mean = ZERO
    this % xs_sem  = ZERO
    this % rel_unc = ZERO
    this % xs_sum  = ZERO
    this % xs_sum2 = ZERO

  end subroutine stats_flush

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! INTERP_FACTOR computes an interpolation factor
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  function interp_factor(val, val_low, val_up, scheme) result(factor)

    real(8) :: val     ! value we're interpolating to
    real(8) :: val_low ! lower bounding value
    real(8) :: val_up  ! upper bounding value
    real(8) :: factor  ! interpolation factor
    integer :: scheme ! interpolation scheme

    select case(scheme)
    case(LINEAR_LINEAR)

      factor = (val - val_low) / (val_up - val_low)

    case(LOG_LOG)

      factor = log(val / val_low) / log(val_up / val_low)

    case default
 
      call fatal_error('Interpolations other than lin-lin or log-log currently not &
        & supported in OTF URR treatments')

    end select

  end function interp_factor

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! INTERPOLATOR computes an interpolation (and an extrapolation, currently, too)
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  function interpolator(factor, val_low, val_up, scheme) result(val)

    real(8) :: factor  ! interpolation factor
    real(8) :: val_low ! lower bounding value
    real(8) :: val_up  ! upper bounding value
    real(8) :: val     ! interpolated value
    integer :: scheme ! interpolation scheme

    select case(scheme)
    case (LINEAR_LINEAR)

      val = val_low + factor * (val_up - val_low)

    case (LOG_LOG)

      val = exp((ONE - factor) * log(val_low) + factor * log(val_up))

    case default

      call fatal_error('Interpolations other than lin-lin or log-log currently not &
        & supported in OTF URR treatments')

    end select

  end function interpolator

end module unresolved
