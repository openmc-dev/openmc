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
  logical :: competitive   ! competitve reaction xs resonance structure?
  logical :: urr_pointwise ! pointwise cross section calculation?
  character(80), allocatable :: urr_endf_filenames(:) ! ENDF filename list
  character(80)              :: urr_formalism         ! URR formalism
  character(80)              :: urr_frequency         ! freq of realizations
  integer, allocatable :: urr_zaids(:)    ! ZAID's for URR nuclides
  integer, allocatable :: n_resonances(:) ! # URR resonances for each l-wave
  integer :: n_otf_urr_xs       ! number of nuclides to calc otf urr xs for
  integer :: n_avg_urr_nuclides ! number of nuclides to calc average urr xs for
  integer :: n_urr_method       ! number of nuclides to treat with a new method
  integer :: l_waves(4)         ! number of contributing l-wave
  integer :: urr_avg_batches    ! min number of batches for inf dil xs calc
  integer :: urr_avg_histories  ! number of histories for inf dil xs calc
  real(8) :: urr_avg_tol ! max rel error for inf dil xs calc termination

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
      procedure :: reset_resonance => reset_resonance

      ! sample unresolved resonance parameters
      procedure :: sample_parameters => sample_parameters

      ! sample level spacing
      procedure :: level_spacing => level_spacing

      ! sample channel widths
      procedure :: channel_width => channel_width

      ! interface for calculation of partial cross sections at E_0
      procedure :: calc_xs => calc_xs

      ! calculate SLBW partial cross sections
      procedure :: slbw_xs => slbw_xs

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
      procedure :: accum_resonance => accum_resonance

      ! add contribution of potential scattering cross section
      procedure :: potential_xs => potential_xs

      ! accumulate values for a single history
      procedure :: accum_history => accum_history

      ! clear values for a single history
      procedure :: flush_history => flush_history

      ! accumulate values for a single batch
      procedure :: accum_batch => accum_batch

      ! clear values for a single batch
      procedure :: flush_batch => flush_batch

      ! calculate batch statistics
      procedure :: calc_stats => calc_stats

      ! clear batch statistics
      procedure :: flush_stats => flush_stats

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
    integer :: i_res       ! resonance counter
    integer :: n_res       ! number of l-wave resonances to include
    integer :: n_above_urr ! number of resonances abover upper URR energy
    real(8) :: E_res ! current resonance (lab) energy (e.g. E_lam)
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

          ! set energy of the lowest-lying contributing URR resonance
          nuc % D = nuc % D_mean(i_l) % data(i_J) % data(1)
          if (i_l > nuc % NLS(nuc % i_urr - 1)) then
            ! the URR has more l-states than the RRR; place resonance energy
            ! randomly about lower URR energy bound
            E_res = nuc % EL(nuc % i_urr) &
                & + (ONE - TWO * prn()) * wigner_dist(nuc % D)
          else
            ! offset first URR resonance energy from the highest-energy RRR
            ! resonance with the same (l,J) spin sequence
            E_res = E_last_rrr(i_nuc, nuc % L, nuc % J) + wigner_dist(nuc % D)
          end if

          ! resonance index
          i_res = 1

          n_above_urr = 0

          RESONANCE_LOOP: do while(n_above_urr < n_res/2)

            ! sample unresolved resonance parameters for this spin
            ! sequence, at this energy
            res % E_lam = E_res
            call res % channel_width(i_nuc)
            call add_parameters(res, i_nuc, i_res, i_l, i_J)

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
    integer :: i_nuc     ! nuclide index
    integer :: i_E       ! URR resonance parameters energy grid index
    integer :: i_l       ! orbital quantum #
    integer :: i_J       ! total angular momentum quantum #
    integer :: n_res     ! number of contributing l-state resonances
    integer :: n_rrr_res ! number of RRR resonances we need to grab
    integer :: i_low     ! index of lowest-lying resonance
    integer :: i_res     ! resonance counter
    integer :: i_rrr_res ! RRR resonance index
    integer :: i_grid    ! xs energy grid index
    integer :: n_pts     ! xs energy grid point counter
    integer :: i_ES      ! index of current URR tabulated energy
    real(8) :: m            ! URR resonance parameters interpolation factor
    real(8) :: f            ! cross section energy grid interpolation factor
    real(8) :: avg_urr_n_xs ! averaged elastic cross section
    real(8) :: avg_urr_f_xs ! averaged fission cross section
    real(8) :: avg_urr_g_xs ! averaged capture cross section
    real(8) :: avg_urr_x_xs ! averaged competitive inelastic cross section
    real(8) :: dE_trial     ! trial energy spacing for xs grid
    real(8) :: xs_trial     ! trial xs via interpolation between gridpoints
    real(8) :: rel_err      ! relative error between interpolated and exact xs
    logical :: enhance ! refine energy-xs grid?
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

      n_pts = 1

      ! enforce xs continuity at RRR-URR energy crossover
      i_grid = binary_search(1.0e6_8 * nuc % energy, nuc % n_grid, nuc % E)
      nuc % urr_energy_tmp(n_pts)  = 1.0e6_8 * nuc % energy(i_grid)
      nuc % urr_elastic_tmp(n_pts) = nuc % elastic(i_grid)
      nuc % urr_capture_tmp(n_pts) = nuc % absorption(i_grid) &
                                 & - nuc % fission(i_grid)
      nuc % urr_fission_tmp(n_pts) = nuc % fission(i_grid)
      nuc % urr_total_tmp(n_pts)   = nuc % total(i_grid)

      nuc % E = 1.0e6_8 * nuc % energy(i_grid + 1)

      allocate(i_last(nuc % NLS(nuc % i_urr)))
      do i_l = 1, nuc % NLS(nuc % i_urr)
        allocate(i_last(i_l) % data(nuc % NJS(i_l)))
        i_last(i_l) % data = ONE
      end do

      i_ES = 1

      ENERGY_LOOP: do

        dE_trial = nuc % urr_dE
        nuc % E = nuc % E + dE_trial
        if (nuc % E > nuc % EH(nuc % i_urr) + dE_trial) exit
        if (i_ES < nuc % NE) then
          if (nuc % E > nuc % ES(i_ES + 1)) then
            i_ES = i_ES + 1
            write(*,'(A40,ES23.16,A12)') &
              & 'Reconstructing URR xs in', nuc % ES(i_ES), ' eV interval'
          end if
        end if
        n_pts = n_pts + 1
        enhance  = .true.

        do while(enhance)

          if (nuc % E > nuc % ES(nuc % NE)) then
            i_E = nuc % NE - 1
          else
            i_E = binary_search(nuc % ES, nuc % NE, nuc % E)
          end if

          m = interp_factor(nuc % E, nuc % ES(i_E), nuc % ES(i_E+1), nuc % INT)

          ! reset xs objects
          call flush_histories(sig_t, sig_n, sig_gam, sig_f, sig_x)

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
              if (i_low == 0) i_low = 1
              i_last(i_l) % data(i_J) = dble(i_low)

              ! loop over the addition of resonances to this ladder
              if (i_low - n_res/2 < 1) then
                ! if we're near the lower end of the URR, need to incorporate
                ! resolved resonance region resonances in order to fix-up
                ! (i.e. smooth out) cross sections at the RRR-URR crossover
                ! energy

                ! if the RRR has resonances with this l-state
                if (i_l <= nuc % NLS(nuc % i_urr - 1)) then

                  ! how many RRR resonances are contributing
                  n_rrr_res = abs(i_low - n_res/2) + 1

                  ! loop over contributing resolved resonance region resonances
                  RRR_RESONANCES_LOOP: do i_res = n_rrr_res, 1, -1

                    i_rrr_res = rrr_res(i_nuc, i_res, nuc % L, nuc % J)

                    call add_resonance(res, i_nuc, i_rrr_res, i_l, i_J, &
                      & nuc % LRF(nuc % i_urr - 1), &
                      & sig_t, sig_n, sig_gam, sig_f, sig_x)

                  end do RRR_RESONANCES_LOOP
                end if

                ! loop over contributing unresolved resonance region resonances
                URR_RESONANCES_LOOP: do i_res = 1, i_low + n_res/2 - 1

                  call add_resonance(res, i_nuc, i_res, i_l, i_J, &
                    & nuc % LRF(nuc % i_urr), &
                    & sig_t, sig_n, sig_gam, sig_f, sig_x)

                end do URR_RESONANCES_LOOP

              else
                ! we're firmly in the URR and can ignore anything going on in
                ! the upper resolved resonance region energies
                URR_LOOP: do i_res = i_low - n_res/2, i_low + n_res/2 - 1

                  call add_resonance(res, i_nuc, i_res, i_l, i_J, &
                    & nuc % LRF(nuc % i_urr), &
                    & sig_t, sig_n, sig_gam, sig_f, sig_x)

                end do URR_LOOP
              end if
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
            i_grid = binary_search(1.0e6_8*nuc % energy, nuc % n_grid, nuc % E)
          end if

          ! check for rare case where two energy points are the same
          if (nuc % energy(i_grid) == nuc % energy(i_grid+1)) &
            & i_grid = i_grid + 1

          ! calculate xs energy grid interpolation factor
          f = interp_factor(nuc % E, 1.0e6_8 * nuc % energy(i_grid), &
            & 1.0e6_8 * nuc % energy(i_grid + 1), LINEAR_LINEAR)

          ! calculate evaluator-supplied backgrounds at the current energy
          call calc_urr_bckgrnd(i_nuc, n_pts, f, i_grid, LINEAR_LINEAR)

          ! interpret MF3 data according to ENDF-6 LSSF flag:
          ! MF3 contains background xs, add to MF2 resonance contributions
          if (nuc % LSSF == 0) then
            call fatal_error('LSSF = 0 not yet supported')

            ! add resonance xs component to background
            call add2bckgrnd(i_nuc, n_pts, sig_t, sig_n, sig_gam, sig_f, sig_x)

          elseif (nuc % LSSF == 1) then
            ! multipy the self-shielding factors by the infinite-dilute xs

            call set_avg_urr_xs(m, i_nuc, i_E, &
              & avg_urr_n_xs, avg_urr_f_xs, avg_urr_g_xs, avg_urr_x_xs)

            ! competitive xs
            if (avg_urr_x_xs > ZERO) then
              if (competitive) then
                ! self-shielded treatment of competitive inelastic xs
                nuc % urr_inelastic_tmp(n_pts) = sig_x % val / avg_urr_x_xs &
                  & * nuc % urr_inelastic_tmp(n_pts)
              else
                ! infinite-dilute treatment of competitive inelastic xs
                nuc % urr_inelastic_tmp(n_pts) = nuc % urr_inelastic_tmp(n_pts)
              end if
            else
              ! use background competitive inelastic cross section, as is
              nuc % urr_inelastic_tmp(n_pts) = nuc % urr_inelastic_tmp(n_pts)
            end if

            ! elastic scattering xs
            nuc % urr_elastic_tmp(n_pts) = sig_n % val / avg_urr_n_xs &
              & * nuc % urr_elastic_tmp(n_pts)

            ! set negative SLBW elastic xs to zero
            if (nuc % urr_elastic_tmp(n_pts) < ZERO) &
              & nuc % urr_elastic_tmp(n_pts) = ZERO

            ! radiative capture xs
            nuc % urr_capture_tmp(n_pts) = sig_gam % val / avg_urr_g_xs &
              & * nuc % urr_capture_tmp(n_pts)

            ! fission xs
            if (avg_urr_f_xs > ZERO) then
              nuc % urr_fission_tmp(n_pts) = sig_f % val / avg_urr_f_xs &
                & * nuc % urr_fission_tmp(n_pts)
            else
              nuc % urr_fission_tmp(n_pts) = nuc % urr_fission_tmp(n_pts)
            end if

            nuc % urr_total_tmp(n_pts) &
              & = nuc % urr_elastic_tmp(n_pts)&
              & + nuc % urr_capture_tmp(n_pts)&
              & + nuc % urr_fission_tmp(n_pts)&
              & + nuc % urr_inelastic_tmp(n_pts)

          else
            call fatal_error('ENDF-6 LSSF not allowed - must be 0 or 1.')
          end if

          xs_trial = HALF &
            & * (nuc % urr_total_tmp(n_pts) + nuc % urr_total_tmp(n_pts - 1))
          dE_trial = HALF * dE_trial
          nuc % E  = nuc % E - dE_trial

          if (nuc % E > nuc % ES(nuc % NE)) then
            i_E = nuc % NE - 1
          else
            i_E = binary_search(nuc % ES, nuc % NE, nuc % E)
          end if

          m = interp_factor(nuc % E, nuc % ES(i_E), nuc % ES(i_E+1), nuc % INT)

          ! reset xs objects
          call flush_histories(sig_t, sig_n, sig_gam, sig_f, sig_x)

          ! loop over orbital quantum #'s
          ORBITAL_ANG_MOM_LOOP1: do i_l = 1, nuc % NLS(nuc % i_urr)

            ! set current orbital angular momentum quantum #
            nuc % L = i_l - 1

            ! get the number of contributing l-wave resonances for this l
            n_res = l_wave_resonances(nuc % L)

            ! loop over total angular momentum quantum #'s
            TOTAL_ANG_MOM_LOOP1: do i_J = 1, nuc % NJS(i_l)

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
              if (i_low == 0) i_low = 1
              i_last(i_l) % data(i_J) = dble(i_low)

              ! loop over the addition of resonances to this ladder
              if (i_low - n_res/2 < 1) then
                ! if we're near the lower end of the URR, need to incorporate
                ! resolved resonance region resonances in order to fix-up
                ! (i.e. smooth out) cross sections at the RRR-URR crossover
                ! energy

                ! if the RRR has resonances with this l-state
                if (i_l <= nuc % NLS(nuc % i_urr - 1)) then

                  ! how many RRR resonances are contributing
                  n_rrr_res = abs(i_low - n_res/2) + 1

                  ! loop over contributing resolved resonance region resonances
                  RRR_RESONANCES_LOOP1: do i_res = n_rrr_res, 1, -1

                    i_rrr_res = rrr_res(i_nuc, i_res, nuc % L, nuc % J)

                    call add_resonance(res, i_nuc, i_rrr_res, i_l, i_J, &
                      & nuc % LRF(nuc % i_urr - 1), &
                      & sig_t, sig_n, sig_gam, sig_f, sig_x)

                  end do RRR_RESONANCES_LOOP1
                end if

                ! loop over contributing unresolved resonance region resonances
                URR_RESONANCES_LOOP1: do i_res = 1, i_low + n_res/2 - 1

                  call add_resonance(res, i_nuc, i_res, i_l, i_J, &
                    & nuc % LRF(nuc % i_urr), &
                    & sig_t, sig_n, sig_gam, sig_f, sig_x)

                end do URR_RESONANCES_LOOP1

              else
                ! we're deep in the URR and can ignore anything going on in the
                ! upper resolved resonance region energies
! TODO: move sig_t outside of loop
                RESONANCES_LOOP1: do i_res &
                  & = i_low - n_res/2, i_low + n_res/2 - 1

                  call add_resonance(res, i_nuc, i_res, i_l, i_J, &
                    & nuc % LRF(nuc % i_urr), &
                    & sig_t, sig_n, sig_gam, sig_f, sig_x)

                end do RESONANCES_LOOP1
              end if
            end do TOTAL_ANG_MOM_LOOP1
          end do ORBITAL_ANG_MOM_LOOP1

          ! add potential scattering contribution
          call sig_n % potential_xs(i_nuc)
          call sig_t % potential_xs(i_nuc)

          if (nuc % E < 1.0e6_8 * nuc % energy(1)) then
            i_grid = 1
          elseif (nuc % E > 1.0e6_8 * nuc % energy(nuc % n_grid)) then
            i_grid = nuc % n_grid - 1
          else
            i_grid = binary_search(1.0e6_8*nuc % energy, nuc % n_grid, nuc % E)
          end if

          ! check for rare case where two energy points are the same
          if (nuc % energy(i_grid) == nuc % energy(i_grid+1)) &
            & i_grid = i_grid + 1

          ! calculate xs energy grid interpolation factor
          f = interp_factor(nuc % E, 1.0e6_8 * nuc % energy(i_grid), &
            & 1.0e6_8 * nuc % energy(i_grid + 1), LINEAR_LINEAR)

          ! calculate evaluator-supplied backgrounds at the current energy
          call calc_urr_bckgrnd(i_nuc, n_pts, f, i_grid, LINEAR_LINEAR)

          ! interpret MF3 data according to ENDF-6 LSSF flag:
          ! MF3 contains background xs, add to MF2 resonance contributions
          if (nuc % LSSF == 0) then
            call fatal_error('LSSF = 0 not yet supported')

            ! add resonance xs component to background
            call add2bckgrnd(i_nuc, n_pts, sig_t, sig_n, sig_gam, sig_f, sig_x)

          elseif (nuc % LSSF == 1) then
            ! multipy the self-shielding factors by the infinite-dilute xs

            call set_avg_urr_xs(m, i_nuc, i_E, &
              & avg_urr_n_xs, avg_urr_f_xs, avg_urr_g_xs, avg_urr_x_xs)

            ! competitive xs
            if (avg_urr_x_xs > ZERO) then
              if (competitive) then
                ! self-shielded treatment of competitive inelastic xs
                nuc % urr_inelastic_tmp(n_pts) = sig_x % val / avg_urr_x_xs &
                  & * nuc % urr_inelastic_tmp(n_pts)
              else
                ! infinite-dilute treatment of competitive inelastic xs
                nuc % urr_inelastic_tmp(n_pts) = nuc % urr_inelastic_tmp(n_pts)
              end if
            else
              ! use background competitive inelastic cross section, as is
              nuc % urr_inelastic_tmp(n_pts) = nuc % urr_inelastic_tmp(n_pts)
            end if

            ! elastic scattering xs
            nuc % urr_elastic_tmp(n_pts) = sig_n % val / avg_urr_n_xs &
              & * nuc % urr_elastic_tmp(n_pts)

            ! set negative SLBW elastic xs to zero
            if (nuc % urr_elastic_tmp(n_pts) < ZERO) &
              & nuc % urr_elastic_tmp(n_pts) = ZERO

            ! radiative capture xs
            nuc % urr_capture_tmp(n_pts) = sig_gam % val / avg_urr_g_xs &
              & * nuc % urr_capture_tmp(n_pts)

            ! fission xs
            if (avg_urr_f_xs > ZERO) then
              nuc % urr_fission_tmp(n_pts) = sig_f % val / avg_urr_f_xs &
                & * nuc % urr_fission_tmp(n_pts)
            else
              nuc % urr_fission_tmp(n_pts) = nuc % urr_fission_tmp(n_pts)
            end if

            nuc % urr_total_tmp(n_pts) = nuc % urr_elastic_tmp(n_pts)&
              &                        + nuc % urr_capture_tmp(n_pts)&
              &                        + nuc % urr_fission_tmp(n_pts)&
              &                        + nuc % urr_inelastic_tmp(n_pts)

          else
            call fatal_error('ENDF-6 LSSF not allowed - must be 0 or 1.')
          end if

          rel_err = abs(xs_trial - nuc % urr_total_tmp(n_pts)) &
            & / nuc % urr_total_tmp(n_pts)
          if (rel_err < nuc % urr_tol .or. dE_trial < 1.0e-11_8) then
            enhance = .false.
          end if

        end do

        ! add energy point to grid
        nuc % E = nuc % E + dE_trial
        nuc % urr_energy_tmp(n_pts) = nuc % E

        if (nuc % E > nuc % ES(nuc % NE)) then
          i_E = nuc % NE - 1
        else
          i_E = binary_search(nuc % ES, nuc % NE, nuc % E)
        end if

        m = interp_factor(nuc % E, nuc % ES(i_E), nuc % ES(i_E + 1), nuc % INT)

        ! reset xs objects
        call flush_histories(sig_t, sig_n, sig_gam, sig_f, sig_x)

        ! loop over orbital quantum #'s
        ORBITAL_ANG_MOM_LOOP2: do i_l = 1, nuc % NLS(nuc % i_urr)

          ! set current orbital angular momentum quantum #
          nuc % L = i_l - 1

          ! get the number of contributing l-wave resonances for this l
          n_res = l_wave_resonances(nuc % L)

          ! loop over total angular momentum quantum #'s
          TOTAL_ANG_MOM_LOOP2: do i_J = 1, nuc % NJS(i_l)

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
            if (i_low == 0) i_low = 1
            i_last(i_l) % data(i_J) = dble(i_low)

            ! loop over the addition of resonances to this ladder
            if (i_low - n_res/2 < 1) then
              ! if we're near the lower end of the URR, need to incorporate
              ! resolved resonance region resonances in order to fix-up
              ! (i.e. smooth out) cross sections at the RRR-URR crossover
              ! energy

              ! if the RRR has resonances with this l-state
              if (i_l <= nuc % NLS(nuc % i_urr - 1)) then

                ! how many RRR resonances are contributing
                n_rrr_res = abs(i_low - n_res/2) + 1

                ! loop over contributing resolved resonance region resonances
                RRR_RESONANCES_LOOP2: do i_res = n_rrr_res, 1, -1

                  i_rrr_res = rrr_res(i_nuc, i_res, nuc % L, nuc % J)

                  call add_resonance(res, i_nuc, i_rrr_res, i_l, i_J, &
                    & nuc % LRF(nuc % i_urr - 1), &
                    & sig_t, sig_n, sig_gam, sig_f, sig_x)

                end do RRR_RESONANCES_LOOP2
              end if

              ! loop over contributing unresolved resonance region resonances
              URR_RESONANCES_LOOP2: do i_res = 1, i_low + n_res/2 - 1

                call add_resonance(res, i_nuc, i_res, i_l, i_J, &
                  & nuc % LRF(nuc % i_urr), &
                  & sig_t, sig_n, sig_gam, sig_f, sig_x)

              end do URR_RESONANCES_LOOP2

            else
              ! we're firmly in the URR and can ignore anything going on in the
              ! upper resolved resonance region energies
              RESONANCES_LOOP2: do i_res = i_low - n_res/2, i_low + n_res/2 - 1

                call add_resonance(res, i_nuc, i_res, i_l, i_J, &
                  & nuc % LRF(nuc % i_urr), &
                  & sig_t, sig_n, sig_gam, sig_f, sig_x)

              end do RESONANCES_LOOP2
            end if
          end do TOTAL_ANG_MOM_LOOP2
        end do ORBITAL_ANG_MOM_LOOP2

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

        ! calculate evaluator-supplied backgrounds at the current energy
        call calc_urr_bckgrnd(i_nuc, n_pts, f, i_grid, LINEAR_LINEAR)

        ! interpret MF3 data according to ENDF-6 LSSF flag:
        ! MF3 contains background xs values, add to MF2 resonance contributions
        if (nuc % LSSF == 0) then
          call fatal_error('LSSF = 0 not yet supported')

          ! add resonance xs component to background
          call add2bckgrnd(i_nuc, n_pts, sig_t, sig_n, sig_gam, sig_f, sig_x)

        elseif (nuc % LSSF == 1) then
          ! multipy the self-shielding factors by the average (infinite-dilute)
          ! cross sections

! TODO: LOG_LOG?
          avg_urr_n_xs = interpolator(m, &
            & nuc % avg_urr_n(i_E), nuc % avg_urr_n(i_E + 1), LINEAR_LINEAR)
          avg_urr_f_xs = interpolator(m, &
            & nuc % avg_urr_f(i_E), nuc % avg_urr_f(i_E + 1), LINEAR_LINEAR)
          avg_urr_g_xs = interpolator(m, &
            & nuc % avg_urr_g(i_E), nuc % avg_urr_g(i_E + 1), LINEAR_LINEAR)
          avg_urr_x_xs = interpolator(m, &
            & nuc % avg_urr_x(i_E), nuc % avg_urr_x(i_E + 1), LINEAR_LINEAR)

          ! competitive xs
          if (avg_urr_x_xs > ZERO) then
            if (competitive) then
              ! self-shielded treatment of competitive inelastic cross section
              nuc % urr_inelastic_tmp(n_pts) = sig_x % val / avg_urr_x_xs &
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
          nuc % urr_elastic_tmp(n_pts) = sig_n % val / avg_urr_n_xs &
            & * nuc % urr_elastic_tmp(n_pts)

          ! set negative SLBW elastic xs to zero
          if (nuc % urr_elastic_tmp(n_pts) < ZERO) &
            & nuc % urr_elastic_tmp(n_pts) = ZERO

          ! radiative capture xs
          nuc % urr_capture_tmp(n_pts) = sig_gam % val / avg_urr_g_xs &
            & * nuc % urr_capture_tmp(n_pts)

          ! fission xs
          if (avg_urr_f_xs > ZERO) then
            nuc % urr_fission_tmp(n_pts) = sig_f % val / avg_urr_f_xs &
              & * nuc % urr_fission_tmp(n_pts)
          else
            nuc % urr_fission_tmp(n_pts) = nuc % urr_fission_tmp(n_pts)
          end if

          nuc % urr_total_tmp(n_pts) = nuc % urr_elastic_tmp(n_pts)&
            &                        + nuc % urr_capture_tmp(n_pts)&
            &                        + nuc % urr_fission_tmp(n_pts)&
            &                        + nuc % urr_inelastic_tmp(n_pts)

        else
          call fatal_error('Self-shielding flag (LSSF) not allowed -&
            & must be 0 or 1.')
        end if

      end do ENERGY_LOOP

      ! pass temporary, pre-allocated energy-xs vectors to dynamic vectors of
      ! the proper length
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
    type(Resonance)    :: res     ! resonance object
    type(CrossSection) :: sig_t   ! total xs object
    type(CrossSection) :: sig_n   ! elastic scattering xs object
    type(CrossSection) :: sig_gam ! radiative capture xs object
    type(CrossSection) :: sig_f   ! fission xs object
    type(CrossSection) :: sig_x   ! competitive inelastic scattering xs object
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
        call flush_statistics(sig_t, sig_n, sig_gam, sig_f, sig_x)

        i_b = 0

        ! loop over batches until convergence
        BATCH_LOOP: do

          i_b = i_b + 1

          ! reset batch accumulators
          call flush_batches(sig_t, sig_n, sig_gam, sig_f, sig_x)

          ! loop over realizations
          HISTORY_LOOP: do i_h = 1, urr_avg_histories

            ! reset accumulator of histories
            call flush_histories(sig_t, sig_n, sig_gam, sig_f, sig_x)

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
                  call accum_resonances(res,sig_t,sig_n,sig_gam,sig_f,sig_x)

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
            call sig_gam % accum_history()
            call sig_x % accum_history()

          end do HISTORY_LOOP

          ! accumulate the result of this batch
          call sig_n % accum_batch()
          call sig_f % accum_batch()
          call sig_gam % accum_batch()
          call sig_x % accum_batch()

          ! calculate statistics for this batch
          call sig_n % calc_stats(i_b)
          call sig_f % calc_stats(i_b)
          call sig_gam % calc_stats(i_b)
          call sig_x % calc_stats(i_b)

! TODO: format avg urr xs output
          if (i_b > urr_avg_batches &
            & .and. max(sig_n % rel_unc, sig_f % rel_unc, sig_gam % rel_unc, &
            & sig_x % rel_unc) < urr_avg_tol) then
            if (1==1) then
              write(*,'(I5,ES10.3,ES10.3,ES10.3,ES10.3,ES10.3,ES10.3,ES10.3,ES10.3,ES10.3)') &
                & i_b, E, sig_n % xs_mean, sig_n % xs_sem, &
                &         sig_f % xs_mean, sig_f % xs_sem, &
                &         sig_gam % xs_mean, sig_gam % xs_sem, &
                &         sig_x % xs_mean, sig_x % xs_sem
            end if
            exit
          end if

        end do BATCH_LOOP

        ! set infinite-dilute xs values at this energy to converged means
        nuc % avg_urr_n(i_E) = sig_n % xs_mean
        nuc % avg_urr_f(i_E) = sig_f % xs_mean
        nuc % avg_urr_g(i_E) = sig_gam % xs_mean
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
    real(8) :: avg_urr_n_xs ! infinite-dilute n xs from NJOY's MC^2 quadrature
    real(8) :: avg_urr_f_xs ! infinite-dilute f xs from NJOY's MC^2 quadrature
    real(8) :: avg_urr_g_xs ! infinite-dilute g xs from NJOY's MC^2 quadrature
    real(8) :: avg_urr_x_xs ! infinite-dilute x xs from NJOY's MC^2 quadrature
    real(8) :: inelastic_xs ! competitive inelastic scattering cross section
    real(8) :: capture_xs   ! radiative capture cross section
    real(8) :: E            ! neutron energy [eV]
    real(8) :: m            ! pointwise xs energy interpolation factor
    real(8) :: f            ! resonance parameters energy interpolation factor
!$omp threadprivate(nuc) 

    ! Set pointer to nuclide
    nuc => nuclides(i_nuc)

    micro_xs(i_nuc) % use_ptable = .true.

    if (nuc % point_urr_xs) then
      i_E = binary_search(nuc % urr_energy, size(nuc % urr_energy), E)
      m = interp_factor(E, nuc % urr_energy(i_E), nuc % urr_energy(i_E + 1), &
        & nuc % INT)
      micro_xs(i_nuc) % elastic = interpolator(m, &
        & nuc % urr_elastic(i_E), nuc % urr_elastic(i_E + 1), nuc % INT)
      micro_xs(i_nuc) % fission = interpolator(m, &
        & nuc % urr_fission(i_E), nuc % urr_fission(i_E + 1), nuc % INT)
      micro_xs(i_nuc) % absorption = interpolator(m, &
        & nuc % urr_capture(i_E) + nuc % urr_fission(i_E), &
        & nuc % urr_capture(i_E + 1) + nuc % urr_fission(i_E + 1), nuc % INT)
      if ((nuc % urr_inelastic(i_E) >= ZERO) &
        & .and. (nuc % urr_inelastic(i_E + 1) >= ZERO)) then
        inelastic_xs = interpolator(m, &
          & nuc % urr_inelastic(i_E), nuc % urr_inelastic(i_E + 1), nuc % INT)
      else
        inelastic_xs = ZERO
      end if
      micro_xs(i_nuc) % total = micro_xs(i_nuc) % elastic &
                            & + micro_xs(i_nuc) % absorption &
                            & + inelastic_xs

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
    call flush_histories(sig_t, sig_n, sig_gam, sig_f, sig_x)

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
          call accum_resonances(res, sig_t, sig_n, sig_gam, sig_f, sig_x)

        end do RESONANCES_LOOP
      end do TOTAL_ANG_MOM_LOOP
    end do ORBITAL_ANG_MOM_LOOP

    ! add potential scattering contribution
    call sig_t % potential_xs(i_nuc)
    call sig_n % potential_xs(i_nuc)

    ! determine energy table
    i_energy = 1
    do
      if (E <= nuc % ES(i_energy + 1)) exit
      i_energy = i_energy + 1
      if (i_energy >= nuc % NE) then
        i_energy = nuc % NE - 1
        exit
      end if
    end do

    ! interpret MF3 data according to ENDF-6 LSSF flag:
    ! MF3 contains background xs values, add to MF2 resonance contributions
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
      ! TODO: LOG_LOG or maybe just nuc % INT?

      ! tabulated unresolved resonance parameters interpolation factor
      f = interp_factor(E, nuc % ES(i_energy), nuc % ES(i_energy+1), nuc % INT)

      ! interpolate evaluator-supplied URR background xs
      call set_avg_urr_xs(f, i_nuc, i_energy, &
        avg_urr_n_xs, avg_urr_f_xs, avg_urr_g_xs, avg_urr_x_xs)

      if (avg_urr_x_xs > ZERO) then

        if (competitive) then
          ! self-shielded treatment of competitive inelastic cross section
          inelastic_xs = sig_x % val / avg_urr_x_xs &
            & * (micro_xs(i_nuc) % total &
            & - micro_xs(i_nuc) % absorption &
            & - micro_xs(i_nuc) % elastic)
        else
          ! infinite-dilute treatment of competitive inelastic cross section
          inelastic_xs = micro_xs(i_nuc) % total &
            & - micro_xs(i_nuc) % absorption &
            & - micro_xs(i_nuc) % elastic
        end if

      else
        inelastic_xs = micro_xs(i_nuc) % total &
          & - micro_xs(i_nuc) % absorption &
          & - micro_xs(i_nuc) % elastic
      end if

      micro_xs(i_nuc) % elastic = sig_n % val / avg_urr_n_xs &
        & * micro_xs(i_nuc) % elastic

      ! set negative SLBW elastic xs to zero
      if (micro_xs(i_nuc) % elastic < ZERO) then
        micro_xs(i_nuc) % elastic = ZERO
      end if

      capture_xs = sig_gam % val / avg_urr_g_xs &
        & * (micro_xs(i_nuc) % absorption - micro_xs(i_nuc) % fission)

      if (avg_urr_f_xs > ZERO) then
        micro_xs(i_nuc) % fission = sig_f % val / avg_urr_f_xs &
          & * micro_xs(i_nuc) % fission
      else
        micro_xs(i_nuc) % fission = micro_xs(i_nuc) % fission
      end if

      micro_xs(i_nuc) % absorption = micro_xs(i_nuc) % fission + capture_xs

      micro_xs(i_nuc) % total = micro_xs(i_nuc) % elastic &
        & + micro_xs(i_nuc) % absorption &
        & + inelastic_xs

    else
      call fatal_error('Self-shielding flag (LSSF) not allowed - must be 0 or 1.')
    end if

    ! Determine nu-fission cross section
    if (nuc % fissionable) then
      micro_xs(i_nuc) % nu_fission = nu_total(nuc, E / 1.0e6_8) &
        & * micro_xs(i_nuc) % fission
    end if

  end subroutine calculate_urr_xs_otf

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! RESET_RESONANCE resets the resonance object before proceeding to add
! contributions from the next spin sequence
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine reset_resonance(this, i_nuc)

    class(Resonance), intent(inout) :: this ! pseudo-resonance object

    type(Nuclide), pointer :: nuc => null() ! nuclide pointer
    integer :: i_nuc ! nuclide index

    nuc => nuclides(i_nuc)

    ! prepare to construct a new ladder (reset resonance counter and energies)
    this % i_res     = 0
    this % E_lam_up  = nuc % E
    this % E_lam_low = nuc % E
    this % E_lam_tmp = nuc % E

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
    integer :: i_nuc ! nuclide index
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

  end subroutine level_spacing

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
      n_res = l_waves(1)
    case(1)
      n_res = l_waves(2)
    case(2)
      n_res = l_waves(3)
    case(3)
      n_res = l_waves(4)
    case default
      call fatal_error('Only s-, p-, d-, and f-wave resonances are supported &
        & in ENDF-6')
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
! CHANNEL_WIDTH samples the channel partial widths
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine channel_width(this, i_nuc)

    class(Resonance), intent(inout) :: this ! pseudo-resonance object

    type(Nuclide), pointer :: nuc => null() ! nuclide object pointer
    integer :: i_nuc  ! nuclide index
    integer :: i_tabn ! elastic chi-squared table index
    integer :: i_tabg ! capture chi-squared table index
    integer :: i_tabf ! fission chi-squared table index
    integer :: i_tabx ! competitivechi-squared table index
    real(8) :: rho    ! derived variable
    real(8) :: nu     ! derived variable

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
    call this % slbw_xs(i_nuc)

  end subroutine calc_xs

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! SLBW_XS calculates Single-Level Breit-Wigner cross sections at an energy point
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine slbw_xs(this, i_nuc)

    class(Resonance), intent(inout) :: this ! pseudo-resonance object

    type(Nuclide), pointer :: nuc => null() ! nuclide pointer
    integer :: i_nuc ! nuclide index

    ! energy variables
    real(8) :: k_n     ! center-of-mass neutron wavenumber at E_n
    real(8) :: k_lam   ! center-of-mass neutron wavenumber at E_lam
    real(8) :: E_shift ! shifted resonance energy in the lab system

    ! broadening variables
    real(8) :: theta   ! total width / Doppler width
    real(8) :: x       ! derived variable

    ! unresolved resonance parameters
    real(8) :: Gam_t_n ! sampled energy-dependent total width at E_n
    real(8) :: Gam_n_n ! sampled energy-dependent neutron width at E_n
    real(8) :: sig_lam ! peak resonance cross section
    real(8) :: sig_lam_Gam_t_n_psi

    nuc => nuclides(i_nuc)

    ! set variables
    k_n = wavenumber(nuc % awr, nuc % E)

    k_lam = wavenumber(nuc % awr, this % E_lam)

    E_shift = this % E_lam &
      & + (this % Gam_n * (shift(nuc % L, k_lam*nuc % ac(nuc % i_urr)) &
      & - shift(nuc % L, k_n*nuc % ac(nuc % i_urr)))) &
      & / (TWO * penetration(nuc % L, k_lam*nuc % ac(nuc % i_urr)))

    Gam_n_n = this % Gam_n &
      & * penetration(nuc % L, k_n*nuc % ac(nuc % i_urr)) &
      & / penetration(nuc % L, k_lam*nuc % ac(nuc % i_urr))

    Gam_t_n = this % Gam_t - this % Gam_n + Gam_n_n

    theta = Gam_t_n &
      & / (TWO * sqrt(K_BOLTZMANN * 1.0E6_8 * nuc % T * nuc % E / nuc % awr))

    x = (TWO * (nuc % E - E_shift)) / Gam_t_n

    sig_lam = FOUR * PI / (k_lam * k_lam) * nuc % g_J &
          & * this % Gam_n / this % Gam_t

! TODO: Correct negative scattering xs values to 0 b in the library version of
!       code for use in OpenMC

! TODO: Compute competitive xs contribution correctly

    ! this particular form comes from the NJOY2012 manual
    if (Gam_n_n > ZERO) then
      this % dsig_n = sig_lam * &
        & ((cos(TWO * phase_shift(nuc % L, k_n*nuc % AP(nuc % i_urr))) &
        & - (ONE - Gam_n_n / Gam_t_n)) * psi(theta, x) &
        & + sin(TWO * phase_shift(nuc % L, k_n*nuc % AP(nuc % i_urr))) &
        & * chi(theta, x))
    else
      call fatal_error('Encountered a non-positive elastic scattering width &
        &in the URR')
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

  end subroutine slbw_xs

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
          & / (11025.0_8 + rho2 * (1575.0_8 + rho2 &
          & * (135.0_8 + rho2 * (10.0_8 + rho2))))

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
        S = -(44100.0_8 + rho2 * (4725.0_8 + rho2 * (270.0_8 + 10.0_8 * rho2)))&
          & / (11025.0_8 + rho2 * (1575.0_8 + rho2 * (135.0_8 &
          & + rho2 * (10.0_8 + rho2))))

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
      w_val = faddeeva_w(cmplx(theta * x * HALF, theta * HALF, 8), relerr)
      psi_val = SQRT_PI * HALF * theta &
        & * real(real(w_val, 8), 8)

    ! QUICKW Faddeeva evaluation from Argonne (also used in NJOY - NJOY manual)
    case (QUICK_W)
      psi_val = SQRT_PI * HALF * theta &
        & * real(real(quickw(cmplx(theta * x * HALF, theta * HALF, 8)), 8), 8)

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
      w_val = faddeeva_w(cmplx(theta * x * HALF, theta * HALF, 8), relerr)
      chi_val = SQRT_PI * HALF * theta &
        & * real(aimag(w_val), 8)

    ! QUICKW Faddeeva evaluation from Argonne (also used in NJOY - NJOY manual)
    case (QUICK_W)

      chi_val = SQRT_PI * HALF * theta &
        & * real(aimag(quickw(cmplx(theta * x * HALF, theta * HALF, 8))), 8)

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

  end subroutine accum_resonance

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! POTENTIAL_XS adds the contribution of potential scattering to the elastic
! and total cross sections
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine potential_xs(this, i_nuc)

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

  end subroutine potential_xs

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! ACCUM_HISTORY adds the single-history xs realization to the single-batch
! accumulator
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine accum_history(this)

    class(CrossSection), intent(inout) :: this ! cross section object

    ! accumulate history xs value for this realization
    this % xs_sum_tmp = this % xs_sum_tmp + this % val

  end subroutine accum_history

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! FLUSH_HISTORY flushes the single-history xs realization
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine flush_history(this)

    class(CrossSection), intent(inout) :: this ! cross section object

    ! clear accumulated xs values for this realization
    this % val      = ZERO
    this % val_last = ZERO

  end subroutine flush_history

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! ACCUM_BATCH adds the single-batch xs realization to the overall accumulator
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine accum_batch(this)

    class(CrossSection), intent(inout) :: this ! cross section object

    ! accumulate batch-wise mean
    this % xs_sum = this % xs_sum + this % xs_sum_tmp / dble(urr_avg_histories)

    ! accumulate squared batch-wise mean
    this % xs_sum2 = this % xs_sum2 &
      & + (this % xs_sum_tmp / dble(urr_avg_histories)) &
      & * (this % xs_sum_tmp / dble(urr_avg_histories))

  end subroutine accum_batch

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! FLUSH_BATCH clears the single-batch sum of xs realizations for the batch
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine flush_batch(this)

    class(CrossSection), intent(inout) :: this ! cross section object

    ! clear accumulated xs values for this batch
    this % xs_sum_tmp = ZERO

  end subroutine flush_batch

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! CALC_STATS computes batch-based means and standard errors of those means
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine calc_stats(this, i_bat_int)

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

  end subroutine calc_stats

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! FLUSH_STATS clears the batch statistics and accumulators
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine flush_stats(this)

    class(CrossSection), intent(inout) :: this ! cross section object

    ! clear mean, SEM, and overall accumulators
    this % xs_mean = ZERO
    this % xs_sem  = ZERO
    this % rel_unc = ZERO
    this % xs_sum  = ZERO
    this % xs_sum2 = ZERO

  end subroutine flush_stats

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
 
      call fatal_error('Interpolations other than lin-lin or log-log currently &
        &not supported in OTF URR treatments')

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

      call fatal_error('Interpolations other than lin-lin or log-log currently &
        &not supported in OTF URR treatments')

    end select

  end function interpolator

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! E_LAST_RRR determines the energy of the highest-energy resolved resonance
! region resonance for a given (l,J) spin sequence
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  function E_last_rrr(i_nuc, l_val, J_val) result(E_val)

    integer :: i_nuc ! nuclide index
    integer :: i_res ! RRR resonance index for a given l
    integer :: l_val ! orbital quantum number
    real(8) :: J_val ! total angular momentum quantum number
    real(8) :: E_val ! highest-energy RRR resonance energy for (l,J)
    type(Nuclide), pointer :: nuc => null() ! nuclide pointer
!$omp threadprivate(nuc)

    nuc => nuclides(i_nuc)

    do i_res = size(nuc % rm_resonances(l_val + 1) % E_lam), 1, -1
      if (nuc % rm_resonances(l_val + 1) % AJ(i_res) == J_val &
        .and. nuc%rm_resonances(l_val+1)%E_lam(i_res) < nuc%EL(nuc%i_urr)) then
        E_val = nuc % rm_resonances(l_val + 1) % E_lam(i_res)
        exit
      end if
    end do

  end function E_last_rrr

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! RRR_RES finds the index of the RRR resonance which we need to add the
! contribution of to a URR xs
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  function rrr_res(i_nuc, n_rrr_res, l_val, J_val) result(i_res)

    integer :: i_nuc     ! nuclide index
    integer :: n_rrr_res ! how many RRR resonances to go back
    integer :: cnt_res   ! how many RRR resonances have we gone back
    integer :: i_res     ! index of the RRR resonance cnt_res resonances back
    integer :: l_val ! orbital quantum number
    real(8) :: J_val ! total angular momentum quantum number
    type(Nuclide), pointer :: nuc => null() ! nuclide pointer
!$omp threadprivate(nuc)

    nuc => nuclides(i_nuc)

    cnt_res   = 0

    do i_res = size(nuc % rm_resonances(l_val + 1) % E_lam), 1, -1
      if (nuc % rm_resonances(l_val + 1) % AJ(i_res) == J_val &
        .and. nuc%rm_resonances(l_val+1)%E_lam(i_res) < nuc%EL(nuc%i_urr)) then
        cnt_res = cnt_res + 1
      end if
      if (cnt_res == n_rrr_res) exit
    end do

  end function rrr_res

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! ACCUM_RESONANCES accumulates contribution from an additional resonance
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine accum_resonances(res, sig_t, sig_n, sig_gam, sig_f, sig_x)

    type(Resonance) :: res ! resonance object
    type(CrossSection) :: sig_t   ! total xs object
    type(CrossSection) :: sig_n   ! elastic scattering xs object
    type(CrossSection) :: sig_gam ! radiative capture xs object
    type(CrossSection) :: sig_f   ! fission xs object
    type(CrossSection) :: sig_x   ! competitive inelastic scattering xs object

    call sig_t   % accum_resonance(res % dsig_t)
    call sig_n   % accum_resonance(res % dsig_n)
    call sig_gam % accum_resonance(res % dsig_gam)
    call sig_f   % accum_resonance(res % dsig_f)
    call sig_x   % accum_resonance(res % dsig_x)

  end subroutine accum_resonances

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! FLUSH_STATISTICS zeroes out statistics accumulators
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine flush_statistics(sig_t, sig_n, sig_gam, sig_f, sig_x)

    type(CrossSection) :: sig_t   ! total xs object
    type(CrossSection) :: sig_n   ! elastic scattering xs object
    type(CrossSection) :: sig_gam ! radiative capture xs object
    type(CrossSection) :: sig_f   ! fission xs object
    type(CrossSection) :: sig_x   ! competitive inelastic scattering xs object

    call sig_t   % flush_stats()
    call sig_n   % flush_stats()
    call sig_gam % flush_stats()
    call sig_f   % flush_stats()
    call sig_x   % flush_stats()    

  end subroutine flush_statistics

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! FLUSH_BATCHES zeroes out batch accumulators
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine flush_batches(sig_t, sig_n, sig_gam, sig_f, sig_x)

    type(CrossSection) :: sig_t   ! total xs object
    type(CrossSection) :: sig_n   ! elastic scattering xs object
    type(CrossSection) :: sig_gam ! radiative capture xs object
    type(CrossSection) :: sig_f   ! fission xs object
    type(CrossSection) :: sig_x   ! competitive inelastic scattering xs object

    call sig_t   % flush_batch()
    call sig_n   % flush_batch()
    call sig_gam % flush_batch()
    call sig_f   % flush_batch()
    call sig_x   % flush_batch()    

  end subroutine flush_batches

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! FLUSH_HISTORIES flushes cross section object histories
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine flush_histories(sig_t, sig_n, sig_gam, sig_f, sig_x)

    type(CrossSection) :: sig_n   ! elastic scattering xs object
    type(CrossSection) :: sig_gam ! radiative capture xs object
    type(CrossSection) :: sig_f   ! fission xs object
    type(CrossSection) :: sig_x   ! competitive inelastic scattering xs object
    type(CrossSection) :: sig_t   ! total xs object

    call sig_t   % flush_history()
    call sig_n   % flush_history()
    call sig_gam % flush_history()
    call sig_f   % flush_history()
    call sig_x   % flush_history()

  end subroutine flush_histories

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! ADD_PARAMETERS adds the URR resonance parameters for a single URR resonance to
! the realization
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine add_parameters(res, i_nuc, i_res, i_l, i_J)

    type(Nuclide), pointer, save :: nuc => null() ! nuclide object pointer
    type(Resonance) :: res ! resonance object
    integer :: i_nuc ! nuclide index
    integer :: i_res ! resonance counter
    integer :: i_l   ! orbital quantum number index
    integer :: i_J   ! total angular momentum quantum #
!$omp threadprivate(nuc)

    nuc => nuclides(i_nuc)

    nuc % urr_resonances(i_res, i_l) % E_lam(i_J) = res % E_lam
    nuc % urr_resonances(i_res, i_l) % GN(i_J)    = res % Gam_n
    nuc % urr_resonances(i_res, i_l) % GG(i_J)    = res % Gam_gam
    nuc % urr_resonances(i_res, i_l) % GF(i_J)    = res % Gam_f
    nuc % urr_resonances(i_res, i_l) % GX(i_J)    = res % Gam_x
    nuc % urr_resonances(i_res, i_l) % GT(i_J)    = res % Gam_t

  end subroutine add_parameters

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! SET_PARAMETERS sets the URR resonance parameters for a single URR resonance
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine set_parameters(res, i_nuc, i_res, i_l, i_J, LRF_val)

    type(Nuclide), pointer, save :: nuc => null() ! nuclide object pointer
    type(Resonance) :: res ! resonance object
    integer :: i_nuc   ! nuclide index
    integer :: i_res   ! resonance counter
    integer :: i_l     ! orbital quantum number index
    integer :: i_J     ! total angular momentum quantum #
    integer :: LRF_val ! ENDF-6 LRF resonance parameter representation flag
!$omp threadprivate(nuc)

    nuc => nuclides(i_nuc)

    select case(LRF_val)
    case (1)
      res % E_lam   = nuc % urr_resonances(i_res, i_l) % E_lam(i_J)
      res % Gam_n   = nuc % urr_resonances(i_res, i_l) % GN(i_J)
      res % Gam_gam = nuc % urr_resonances(i_res, i_l) % GG(i_J)
      res % Gam_f   = nuc % urr_resonances(i_res, i_l) % GF(i_J)
      res % Gam_x   = nuc % urr_resonances(i_res, i_l) % GX(i_J)
      res % Gam_t   = nuc % urr_resonances(i_res, i_l) % GT(i_J)
    
    case (3)
      res % E_lam   = nuc % rm_resonances(i_l) % E_lam(i_res)
      res % Gam_n   = nuc % rm_resonances(i_l) % GN(i_res)
      res % Gam_gam = nuc % rm_resonances(i_l) % GG(i_res)
      res % Gam_f   = nuc % rm_resonances(i_l) % GFA(i_res) &
                  & + nuc % rm_resonances(i_l) % GFB(i_res)
      res % Gam_x   = ZERO
      res % Gam_t   = res % Gam_n &
                  & + res % Gam_gam &
                  & + res % Gam_f &
                  & + res % Gam_x

    case default
      write(*,'(A80)') 'Not a supported resonance parameter representation'

    end select

  end subroutine set_parameters

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! SET_MEAN_PARAMETERS sets the URR mean resonance parameters at an energy
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine set_mean_parameters(i_nuc, i_E, E_res, i_l, i_J)

    type(Nuclide), pointer, save :: nuc => null() ! nuclide object pointer
    integer :: i_nuc ! nuclide index
    integer :: i_E   ! tabulated URR parameters energy index
    integer :: i_l   ! orbital quantum number index
    integer :: i_J   ! total angular momentum quantum #
    real(8) :: E_res ! current resonance (lab) energy (e.g. E_lam)
    real(8) :: m     ! energy interpolation factor
!$omp threadprivate(nuc)

    nuc => nuclides(i_nuc)

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

  end subroutine set_mean_parameters

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! ADD_RESONANCE add an additional contributing resonance
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine add_resonance(res, i_nuc, i_res, i_l, i_J, &
                    & LRF_val, sig_t, sig_n, sig_gam, sig_f, sig_x)

    type(Resonance) :: res ! resonance object
    type(CrossSection) :: sig_t   ! total xs object
    type(CrossSection) :: sig_n   ! elastic scattering xs object
    type(CrossSection) :: sig_gam ! radiative capture xs object
    type(CrossSection) :: sig_f   ! fission xs object
    type(CrossSection) :: sig_x   ! competitive inelastic scattering xs object
    integer :: i_nuc   ! nuclide index
    integer :: i_res   ! resonance index
    integer :: i_l     ! orbital quantum number index
    integer :: i_J     ! total angular momentum quantum #
    integer :: LRF_val ! ENDF-6 LRF resonance parameter representation flag

    ! set resonance parameters
    call set_parameters(res, i_nuc, i_res, i_l, i_J, LRF_val)

    ! calculate the contribution to the partial cross sections,
    ! at this energy, from an additional resonance 
    call res % calc_xs(i_nuc)

    ! add this contribution to the accumulated partial cross
    ! section values built up from all resonances
! TODO: move sig_t outside of loop
    call accum_resonances(res, sig_t, sig_n, sig_gam, sig_f, sig_x)

  end subroutine add_resonance

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! ADD_RESONANCE add an additional contributing resonance
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine set_avg_urr_xs(m, i_nuc, i_E, n_xs, f_xs, g_xs, x_xs)

    type(Nuclide), pointer, save :: nuc => null() ! nuclide object pointer
    integer :: i_nuc
    integer :: i_E
    real(8) :: m
    real(8) :: n_xs
    real(8) :: f_xs
    real(8) :: g_xs
    real(8) :: x_xs
!$omp threadprivate(nuc)

    nuc => nuclides(i_nuc)

! TODO: LOG_LOG?
    ! infinite-dilute elastic scattering
    n_xs = interpolator(m, &
      & nuc % avg_urr_n(i_E), nuc % avg_urr_n(i_E + 1), LINEAR_LINEAR)

    ! infinite-dilute fission
    f_xs = interpolator(m, &
      & nuc % avg_urr_f(i_E), nuc % avg_urr_f(i_E + 1), LINEAR_LINEAR)

    ! infinite-dilute capture
    g_xs = interpolator(m, &
      & nuc % avg_urr_g(i_E), nuc % avg_urr_g(i_E + 1), LINEAR_LINEAR)

    ! infinite-dilute competitive reaction xs
    x_xs = interpolator(m, &
      & nuc % avg_urr_x(i_E), nuc % avg_urr_x(i_E + 1), LINEAR_LINEAR)

  end subroutine set_avg_urr_xs

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! CALC_URR_BCKGRND interpolates the evaluator-supplied background xs at the
! current energy
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine calc_urr_bckgrnd(i_nuc, n_pts, f, i_grid, int_style)

    type(Nuclide), pointer, save :: nuc => null() ! nuclide object pointer
    integer :: i_nuc
    integer :: n_pts
    integer :: i_grid
    integer :: int_style ! interpolation scheme
    real(8) :: f
!$omp threadprivate(nuc)

    nuc => nuclides(i_nuc)

    select case(int_style)
    case(LINEAR_LINEAR)
      ! elastic scattering xs
      nuc % urr_elastic_tmp(n_pts) = interpolator(f, &
        & nuc % elastic(i_grid), nuc % elastic(i_grid + 1), &
        & int_style)

      ! radiative capture xs
      nuc % urr_capture_tmp(n_pts) = interpolator(f, &
        & nuc % absorption(i_grid) - nuc % fission(i_grid), &
        & nuc % absorption(i_grid + 1) - nuc % fission(i_grid + 1), &
        & int_style)

      ! fission xs
      nuc % urr_fission_tmp(n_pts) = interpolator(f, &
        & nuc % fission(i_grid), nuc % fission(i_grid + 1), &
        & int_style)

      ! competitive first level inelastic scattering xs
      nuc % urr_inelastic_tmp(n_pts) = interpolator(f, &
        &   nuc % total(i_grid) &
        & - nuc % absorption(i_grid) &
        & - nuc % elastic(i_grid), &
        &   nuc % total(i_grid + 1) &
        & - nuc % absorption(i_grid + 1) &
        & - nuc % elastic(i_grid + 1), &
        & int_style)

      ! total xs
      nuc % urr_total_tmp(n_pts) = interpolator(f, &
        & nuc % total(i_grid), nuc % total(i_grid + 1), &
        & int_style)

    case default
      call fatal_error('Only lin-lin interpolation on the evaluator-supplied &
        &URR background grid is supported.')

    end select

  end subroutine calc_urr_bckgrnd

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! ADD2BCKGRND adds the resonance xs component to the evaluator-supplied
! background
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine add2bckgrnd(i_nuc, n_pts, sig_t, sig_n, sig_gam, sig_f, sig_x)

    type(Nuclide), pointer, save :: nuc => null() ! nuclide object pointer
    type(CrossSection) :: sig_t   ! total xs object
    type(CrossSection) :: sig_n   ! elastic scattering xs object
    type(CrossSection) :: sig_gam ! radiative capture xs object
    type(CrossSection) :: sig_f   ! fission xs object
    type(CrossSection) :: sig_x   ! competitive inelastic scattering xs object
    integer :: i_nuc
    integer :: n_pts
!$omp threadprivate(nuc)

    nuc => nuclides(i_nuc)

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

  end subroutine add2bckgrnd

end module unresolved
