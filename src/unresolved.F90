module unresolved

  use error,         only: fatal_error
  use faddeeva,      only: quickw, faddeeva_w
  use fission,       only: nu_total
  use global
  use list_header,   only: ListInt, ListReal
  use random_lcg,    only: prn
  use search,        only: binary_search
  use vector_header, only: Vector, JaggedArray

  implicit none

  logical :: competitive      ! use competitve reaction xs resonance structure?
  logical :: run_fasturr      ! use special treatment for Fast/URR data?
  logical :: write_avg_urr_xs ! write averaged xs values to a file?
  logical :: write_urr_tables ! write probability tables to a file?
  integer :: background          ! where to get background cross sections
  integer :: band_spacing        ! cross section band spacing scheme
  integer :: E_spacing           ! probability table energy spacing scheme
  integer :: formalism           ! URR resonance formalism
  integer :: histories_avg_urr   ! histories for averaged cross section calc
  integer :: i_real              ! index of URR realization used for calc
  integer :: i_real_user         ! user-specified realization index
  integer :: l_waves(4)          ! number of contributing l-wave
  integer :: max_batches_avg_urr ! max batches for averaged cross section calc
  integer :: min_batches_avg_urr ! min batches for averaged cross section calc
  integer :: n_bands             ! number of probability table xs bands
  integer :: n_fasturr           ! number of URR isotopes being processed
  integer :: n_tables            ! number of probability tables (energies)
  integer :: n_reals             ! number of independent URR realizations
  integer :: real_freq           ! frequency of URR realizations
  integer :: represent_params    ! representation of URR parameters
  integer :: represent_urr       ! representation of URR xs
  real(8) :: first_bound         ! xs boundary between first two bands
  real(8) :: last_bound          ! xs boundary between last two bands
  real(8) :: max_dE_point_urr    ! max diff between reconstructed energies [eV]
  real(8) :: max_E_point_urr     ! max energy for pointwise xs recon [eV]
  real(8) :: min_dE_point_urr    ! min diff between reconstructed energies [eV]
  real(8) :: tol_avg_urr         ! max rel err for inf dil xs calc termination
  real(8) :: tol_point_urr       ! max pointwise xs reconstruction rel err
  real(8), allocatable :: Etables(:)  ! probability table energies
  real(8), allocatable :: xs_bands(:) ! probability table xs band boundaries
  character(MAX_FILE_LEN) :: path_avg_urr_xs ! path to averaged URR xs files
  character(80), allocatable :: endf_files(:) ! list of ENDF-6 filenames
!$omp threadprivate(i_real)

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! RESONANCE is an object containing information about a pseudo-resonance that
! is contributing to the cross section value at an energy grid point in the URR
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  type Resonance

    ! sampled unresolved resonance parameters
    real(8) :: D_lJ      ! sampled nuclear level spacing
    real(8) :: Gam_t     ! sampled total width
    real(8) :: Gam_n     ! sampled neutron width
    real(8) :: Gam_g     ! sampled radiative width
    real(8) :: Gam_f     ! sampled fission width
    real(8) :: Gam_x     ! sampled competitive width

    ! resonance energy variables
    real(8) :: E_lam     ! sampled resonance energy

    ! counter for the number of resonances added for a given spin sequence
    integer :: i_res

    ! partial and total xs contributions
    real(8) :: dxs_t    ! contribution of resonance to E_n total xs
    real(8) :: dxs_n    ! contribution of resonance to E_n scatter xs
    real(8) :: dxs_g    ! contribution of resonance to E_n capture xs
    real(8) :: dxs_f    ! contribution of resonance to E_n fission xs
    real(8) :: dxs_x    ! contribution of resonance to E_n competitive xs

  ! type-bound procedures
  contains

    ! sample unresolved resonance parameters
    procedure :: sample_parameters => sample_parameters

    ! sample level spacing
    procedure :: level_spacing => level_spacing

    ! sample channel widths
    procedure :: channel_width => channel_width

    ! interface for calculation of partial cross sections at E_n
    procedure :: calc_xs => calc_xs

    ! calculate SLBW partial cross sections
    procedure :: slbw_xs => slbw_xs

    ! calculate MLBW partial cross sections
    procedure :: mlbw_xs => mlbw_xs

  end type Resonance

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! SLBWRESONANCES is an object containing a vector of SLBW resonances'
! information
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  type SLBWResonances

    real(8), allocatable :: E_lam(:)
    real(8), allocatable :: AJ(:)
    real(8), allocatable :: GN(:)
    real(8), allocatable :: GG(:)
    real(8), allocatable :: GF(:)
    real(8), allocatable :: GX(:)
    real(8), allocatable :: GT(:)

    ! type-bound procedures
    contains

      ! allocate vector of SLBW resonances for a given l
      procedure :: alloc_slbw_resonances => alloc_slbw_resonances

      ! deallocate vector of SLBW resonances for a given l
      procedure :: dealloc_slbw_resonances => dealloc_slbw_resonances

  end type SLBWResonances

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! MLBWRESONANCES is an object containing a vector of MLBW resonances'
! information
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  type MLBWResonances

    real(8), allocatable :: E_lam(:)
    real(8), allocatable :: AJ(:)
    real(8), allocatable :: GN(:)
    real(8), allocatable :: GG(:)
    real(8), allocatable :: GF(:)
    real(8), allocatable :: GX(:)
    real(8), allocatable :: GT(:)

    ! type-bound procedures
    contains

      ! allocate vector of MLBW resonances for a given l
      procedure :: alloc_mlbw_resonances => alloc_mlbw_resonances

      ! deallocate vector of MLBW resonances for a given l
      procedure :: dealloc_mlbw_resonances => dealloc_mlbw_resonances

  end type MLBWResonances

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! RMRESONANCES is an object containing a vector of Reich-Moore resonance
! data
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  type RMResonances

    real(8), allocatable :: E_lam(:)
    real(8), allocatable :: AJ(:)
    real(8), allocatable :: GN(:)
    real(8), allocatable :: GG(:)
    real(8), allocatable :: GFA(:)
    real(8), allocatable :: GFB(:)

    ! type-bound procedures
    contains

      ! allocate vector of R-M resonances for a given l
      procedure :: alloc_rm_resonances => alloc_rm_resonances

      ! deallocate vector of R-M resonances for a given l
      procedure :: dealloc_rm_resonances => dealloc_rm_resonances

  end type RMResonances

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! CROSSSECTION is an object containing data for a partial (or total) cross
! section
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  type CrossSection

    real(8) :: xs          ! cross section value accumulator
    real(8) :: xs_tmp_sum  ! single-batch sum of xs history values
    integer :: cnt_tmp_sum ! single-batch band count
    real(8) :: xs_sum      ! sum of single-batch xs sums
    real(8) :: xs_sum2     ! sum of squares of single-batch xs sums
    integer :: cnt         ! sum of single-batch band counts
    integer :: cnt2        ! sum of squares of single-batch band counts
    real(8) :: xs_mean     ! overall mean band xs value per band count
    real(8) :: cnt_mean    ! overall mean band count per history
    real(8) :: xs_sem      ! standard error of the mean band xs
    real(8) :: cnt_sem     ! standard error of the mean band count
    real(8) :: rel_unc     ! relative uncertainty (SEM / mean)

  ! type-bound procedures
  contains

    ! accumulate resonance contribution to ladder partial cross section
    procedure :: accum_resonance => accum_resonance

    ! add contribution of potential scattering cross section
    procedure :: potential_xs => potential_xs

    ! clear batch statistics
    procedure :: flush_xs_stats => flush_xs_stats

  end type CrossSection

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! PROBABILITYTABLE is an object containing data for a single table
! (i.e. one isotope, one temperature, one energy)
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  type ProbabilityTable

    type(CrossSection), allocatable :: t(:) ! total xs object
    type(CrossSection), allocatable :: n(:) ! elastic scattering xs object
    type(CrossSection), allocatable :: g(:) ! radiative capture xs object
    type(CrossSection), allocatable :: f(:) ! fission xs object
    type(CrossSection), allocatable :: x(:) ! competitive xs object
    type(CrossSection) :: avg_t   ! infinite-dilute total xs object
    type(CrossSection) :: avg_n   ! infinite-dilute elastic xs object
    type(CrossSection) :: avg_g   ! infinite-dilute capture xs object
    type(CrossSection) :: avg_f   ! infinite-dilute fission xs object
    type(CrossSection) :: avg_x   ! infinite-dilute competitive xs object

  end type ProbabilityTable

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! LOCALSEQUENCE is an object containing vectors of resonance parameters for
! one spin sequence of a local realization of resonances about E_n - the vectors
! are the lengths of the number of contributing resonances for the corresponding
! spin sequence
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  type LocalSequence

     real(8), allocatable :: E_lam(:)
     real(8), allocatable :: Gam_n(:)
     real(8), allocatable :: Gam_f(:)
     real(8), allocatable :: Gam_g(:)
     real(8), allocatable :: Gam_x(:)
     real(8), allocatable :: Gam_t(:)

   ! type-bound procedures
   contains

     ! allocate vector of local resonances for a given l
     procedure :: alloc_local_sequence => alloc_local_sequence

     ! deallocate vector of local resonances for a given l
     procedure :: dealloc_local_sequence => dealloc_local_sequence

  end type LocalSequence

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! ISOTOPE is an object containing data for a single isotope with a URR that is
! to be processed
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  type Isotope

    real(8) :: AWR ! weight of nucleus in neutron masses
    real(8) :: T   ! current isotope temperature [K]
    type(ListReal) :: Tlist ! list of temperatures isotope has ACE data at
    type(ListInt) :: inuclist ! list of indices for different temperatures
    integer :: nT = 0 ! number of temperatures isotope has ACE data at
    logical :: fissionable ! is isotope fissionable?

    ! ENDF-6 nuclear data
    logical :: been_read = .false. ! has ENDF-6 data already been read in?
    integer, allocatable :: NLS(:)  ! number of l values in each energy region
    integer, allocatable :: NJS(:)  ! number of J values for each URR l
    integer, allocatable :: LRU(:)  ! resolved (1) or unresolved (2) parameters
    integer, allocatable :: LRF(:)  ! ENDF-6 resonance formalism number
    integer, allocatable :: NRO(:)  ! AP energy-dependence flag
    integer, allocatable :: NAPS(:) ! channel radius handling flag
    integer :: MAT     ! isotope MAT number
    integer :: ZAI     ! isotope ZAID number
    integer :: LRP     ! resonance parameter flag
    integer :: NER     ! number of resonance energy regions
    integer :: LSSF    ! self-shielding factor flag
    integer :: INT     ! ENDF-6 interpolation law for parameters and xs
    integer :: MF3_INT ! ENDF-6 interpolation law for File 3 xs
! TODO: check for interpolation law consistency
! TODO: use LRX if ZERO's aren't given for inelastic width in ENDF when LRX = 0
    integer :: LRX     ! competitive inelastic width flag
    integer :: NE      ! number of URR tabulated data energies
    real(8) :: E_ex1   ! first level inelastic scattering excitation energy
    real(8) :: E_ex2   ! second level inelastic scattering excitation energy
    real(8), allocatable :: EL(:)  ! lower energy bound of energy region
    real(8), allocatable :: EH(:)  ! upper energy bound of energy region
    real(8), allocatable :: AP(:)  ! scattering radius
    real(8), allocatable :: ac(:)  ! channel radius
    real(8), allocatable :: SPI(:) ! total spin
    real(8), allocatable :: ES(:)  ! URR tabulated data energies

    ! mean values
    type(JaggedArray), allocatable :: D_mean(:)   ! level spacing
    type(JaggedArray), allocatable :: GN0_mean(:) ! reduced neutron width
    type(JaggedArray), allocatable :: GG_mean(:)  ! radiative capture width
    type(JaggedArray), allocatable :: GF_mean(:)  ! fission width
    type(JaggedArray), allocatable :: GX_mean(:)  ! competitive inelastic width
    type(Vector), allocatable :: AJ(:)   ! total angular momentum
    type(Vector), allocatable :: DOFN(:) ! # neutron channels
    type(Vector), allocatable :: DOFG(:) ! # capture channels
    type(Vector), allocatable :: DOFF(:) ! # fission channels
    type(Vector), allocatable :: DOFX(:) ! # competitive channels

    ! current values
    real(8) :: E_last = ZERO ! last neutron energy
    real(8) :: E   ! neutron energy
    real(8) :: J   ! total angular momentum
    real(8) :: g_J ! statistical spin factor
    real(8) :: D   ! level spacing
    real(8) :: GN0 ! reduced neutron width
    real(8) :: GG  ! radiative capture width
    real(8) :: GF  ! fission width
    real(8) :: GX  ! competitive inelastic width
    integer :: L    ! orbital quantum number
    integer :: AMUN ! number of neutron channels (degrees of freedom)
    integer :: AMUG ! number of capture channels (degrees of freedom)
    integer :: AMUF ! number of fission channels (degrees of freedom)
    integer :: AMUX ! number of competitive channels (degrees of freedom)

    ! computed averaged, infinite-dilute URR cross section values and energies
    integer :: nEavg
    real(8), allocatable :: Eavg(:)
    real(8), allocatable :: avg_urr_t(:)
    real(8), allocatable :: avg_urr_n(:)
    real(8), allocatable :: avg_urr_f(:)
    real(8), allocatable :: avg_urr_g(:)
    real(8), allocatable :: avg_urr_x(:)

    ! ENDF-6 File 3 evaluator-supplied background cross sections and energies
    real(8), allocatable :: MF3_n_e(:)
    real(8), allocatable :: MF3_f_e(:)
    real(8), allocatable :: MF3_g_e(:)
    real(8), allocatable :: MF3_x_e(:)
    real(8), allocatable :: MF3_n(:)
    real(8), allocatable :: MF3_f(:)
    real(8), allocatable :: MF3_g(:)
    real(8), allocatable :: MF3_x(:)

    ! URR treatment parameters and indices
    logical :: otf_urr_xs   = .false. ! calculate URR xs on-the-fly?
    logical :: prob_bands   = .false. ! calculate probability tables?
    logical :: point_urr_xs = .false. ! calculate pointwise URR cross sections?
    integer :: i_urr ! index of URR energy range

    ! probability tables for given (energy, temperature) pairs
    integer :: ntabs                 ! number of probability tables (energies)
    real(8), allocatable :: Etabs(:) ! probability table energies
    type(ProbabilityTable), allocatable :: prob_tables(:,:)

    ! vectors of URR resonances for a given (realization, L, J)
    type(SLBWResonances), allocatable :: urr_resonances_tmp(:,:,:)
    type(SLBWResonances), allocatable :: urr_resonances(:,:,:)

    ! max resonances for a given (l,J)
    integer :: n_lam_tmp = 1000000
    integer, allocatable :: n_lam(:,:,:)

    ! vector of Reich-Moore resonances for each l
    type(RMResonances), allocatable :: rm_resonances(:)

    ! vector of SLBW resonances for each l
    type(SLBWResonances), allocatable :: slbw_resonances(:)

    ! vector of MLBW resonances for each l
    type(MLBWResonances), allocatable :: mlbw_resonances(:)

    ! local contributing resonances for (i_l, i_J)
    type(LocalSequence), allocatable :: local_realization(:,:)

    ! pointwise URR cross section data
    real(8), allocatable :: E_tmp(:) ! Scratch energy grid values
    real(8), allocatable :: n_tmp(:) ! scratch elastic xs
    real(8), allocatable :: g_tmp(:) ! scratch capture xs
    real(8), allocatable :: f_tmp(:) ! scratch fission xs
    real(8), allocatable :: x_tmp(:) ! scratch competitive xs
    real(8), allocatable :: t_tmp(:) ! scratch total xs
    real(8), allocatable :: urr_E(:) ! energy grid values
    real(8), allocatable :: urr_n(:) ! elastic xs
    real(8), allocatable :: urr_g(:) ! capture xs
    real(8), allocatable :: urr_f(:) ! fission xs
    real(8), allocatable :: urr_x(:) ! first level inelastic xs
    real(8), allocatable :: urr_t(:) ! total xs

    ! pointwise URR cross section parameters
    integer :: n_urr_gridpoints  = 100000000 ! max URR energy-xs gridpoints

  ! Type-Bound procedures
  contains

    ! allocate resonance energy range variables
    procedure :: alloc_energy_range => alloc_energy_range

    ! deallocate resonance energy range variables
    procedure :: dealloc_energy_ranges => dealloc_energy_ranges

    ! deallocate File 3 average (infinite-dilute) cross sections
    procedure :: dealloc_MF3 => dealloc_MF3

    ! allocate temporary URR resonance ensemble realization
    procedure :: alloc_ensemble_tmp => alloc_ensemble_tmp

    ! deallocate temporary URR resonance ensemble realization
    procedure :: dealloc_ensemble_tmp => dealloc_ensemble_tmp

    ! allocate URR resonance ensemble realization
    procedure :: alloc_ensemble => alloc_ensemble

    ! deallocate URR resonance ensemble realization
    procedure :: dealloc_ensemble => dealloc_ensemble

    ! allocate local URR resonance realization
    procedure :: alloc_local_realization => alloc_local_realization

    ! deallocate local URR resonance realization
    procedure :: dealloc_local_realization => dealloc_local_realization

    ! allocate temporary pointwise URR cross sections
    procedure :: alloc_pointwise_tmp => alloc_pointwise_tmp

    ! deallocate temporary pointwise URR cross sections
    procedure :: dealloc_pointwise_tmp => dealloc_pointwise_tmp

    ! allocate pointwise URR cross sections
    procedure :: alloc_pointwise => alloc_pointwise

    ! deallocate pointwise URR cross sections
    procedure :: dealloc_pointwise => dealloc_pointwise

    ! allocate probability tables
    procedure :: alloc_prob_tables => alloc_prob_tables

    ! deallocate probability tables
    procedure :: dealloc_prob_tables => dealloc_prob_tables

    ! zeros out statistics accumulators
    procedure :: flush_ptable_stats => flush_ptable_stats

    ! zeros out batch accumulators
    procedure :: flush_batches => flush_batches

    ! zeros out xs value accumulator
    procedure :: flush_histories => flush_histories

    ! calculate xs means and standard errors
    procedure :: calc_stats => calc_stats

    ! add the contribution of an additional reference
    procedure :: res_contrib => res_contrib

    ! update batch accumulators with history result
    procedure :: accum_history => accum_history

    ! accumulate values for a single batch
    procedure :: accum_batch => accum_batch

    ! set channel radius
    procedure :: channel_radius => channel_radius

    ! deallocate isotope
    procedure :: dealloc_isotope => dealloc_isotope

  end type Isotope

  type(Isotope), allocatable, target :: isotopes(:)

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! Tabulated Chi-Squared Distribution
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

! TODO: extend tabulated values beyond Mathematica's precision
  ! tabulated chi-square distribution for 1-4 degrees of freedom with values
  ! taken to preserve integral probabilities for equiprobable bins
  real(8), dimension(20,4), parameter :: chi2 = reshape((/&
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

  subroutine resonance_ensemble(iso)

    type(Isotope), pointer :: tope => null() ! isotope object pointer
    type(Resonance) :: res ! resonance object
    integer :: iso         ! isotope index
    integer :: i_l         ! orbital quantum number index
    integer :: i_J         ! total angular momentum quantum number
    integer :: i_E         ! tabulated URR parameters energy index
    integer :: i_ens       ! ensemble index
    integer :: i_res       ! (l, J) resonance counter
    integer :: n_res       ! number of l-wave resonances to include
    integer :: n_above_urr ! number of resonances abover upper URR energy
    real(8) :: E_res ! current resonance (lab) energy (e.g. E_lam)
    real(8) :: m     ! energy interpolation factor
!$omp threadprivate(tope)

    tope => isotopes(iso)

    ! allocate temporary URR resonance ensemble realizations
    call tope % alloc_ensemble_tmp(n_reals)

    allocate(tope % n_lam(n_reals, tope % NLS(tope % i_urr),&
      & maxval(tope % NJS(:))))

    res % i_res = 0

    ! loop over independent realizations
    ENSEMBLE_LOOP: do i_ens = 1, n_reals

      ! loop over orbital angular momenta
      ORBITAL_ANG_MOM_LOOP: do i_l = 1, tope % NLS(tope % i_urr)

        ! set current orbital angular momentum quantum number
        tope % L = i_l - 1

        ! get the number of contributing l-wave resonances for this l
        n_res = n_res_contrib(tope % L)

        ! loop over total angular momenta
        TOTAL_ANG_MOM_LOOP: do i_J = 1, tope % NJS(i_l)

          ! set current total angular momentum quantum number
          tope % J = tope % AJ(i_l) % data(i_J)

          ! set current partial width degrees of freedom
          tope % AMUX = int(tope % DOFX(i_l) % data(i_J))
          tope % AMUN = int(tope % DOFN(i_l) % data(i_J))
          tope % AMUG = int(tope % DOFG(i_l) % data(i_J))
          tope % AMUF = int(tope % DOFF(i_l) % data(i_J))

          ! set energy of the lowest-lying contributing URR resonance
          tope % D = tope % D_mean(i_l) % data(i_J) % data(1)
          if (i_l > tope % NLS(tope % i_urr - 1)) then
            ! the URR has more l-states than the RRR; place resonance energy
            ! randomly about lower URR energy bound
            E_res = tope % EL(tope % i_urr) &
              & + (ONE - TWO * prn()) * wigner_dist(tope % D)
          else
            ! offset first URR resonance energy from the highest-energy RRR
            ! resonance with the same (l,J) spin sequence
            E_res = E_last_rrr(iso, tope % L, tope % J) + wigner_dist(tope % D)
          end if

          i_res = 0
          n_above_urr = 0
          RESONANCE_LOOP: do while(n_above_urr < n_res + 1)

            i_res = i_res + 1

            ! compute interpolation factor
            if (E_res < tope % ES(1)) then
              i_E = 1
            else if (E_res > tope % ES(tope % NE)) then
              i_E = tope % NE - 1
            else
              i_E = binary_search(tope % ES, tope % NE, E_res)
            end if

            m = interp_factor(E_res, tope % ES(i_E), tope % ES(i_E + 1),&
              & tope % INT)

            ! set current mean unresolved resonance parameters
            tope % D = interpolator(m, &
              & tope % D_mean(i_l) % data(i_J) % data(i_E), &
              & tope % D_mean(i_l) % data(i_J) % data(i_E + 1), tope % INT)
            tope % GN0 = interpolator(m, &
              & tope % GN0_mean(i_l) % data(i_J) % data(i_E), &
              & tope % GN0_mean(i_l) % data(i_J) % data(i_E + 1), tope % INT)
            tope % GG = interpolator(m, &
              & tope % GG_mean(i_l) % data(i_J) % data(i_E), &
              & tope % GG_mean(i_l) % data(i_J) % data(i_E + 1), tope % INT)
            if (tope % INT == LINEAR_LINEAR &
              & .or. tope % GF_mean(i_l) % data(i_J) % data(i_E) > ZERO) then
              tope % GF = interpolator(m, &
                & tope % GF_mean(i_l) % data(i_J) % data(i_E), &
                & tope % GF_mean(i_l) % data(i_J) % data(i_E + 1), tope % INT)
            else
              tope % GF = ZERO
            end if
            if (tope % INT == LINEAR_LINEAR &
              .or. tope % GX_mean(i_l) % data(i_J) % data(i_E) > ZERO) then
              tope % GX = interpolator(m, &
                & tope % GX_mean(i_l) % data(i_J) % data(i_E), &
                & tope % GX_mean(i_l) % data(i_J) % data(i_E + 1), tope % INT)
            else
              tope % GX = ZERO
            end if

            ! sample unresolved resonance parameters for this spin
            ! sequence, at this energy
            res % E_lam = E_res
            tope % E = E_res
            call res % channel_width(iso, i_l, i_J)
            call add_parameters(res, iso, i_ens, i_res, i_l, i_J)

            ! add an additional resonance
            E_res = E_res + wigner_dist(tope % D)

            if (E_res > tope % EH(tope % i_urr)) n_above_urr = n_above_urr + 1

          end do RESONANCE_LOOP

          tope % n_lam(i_ens, i_l, i_J) = i_res

        end do TOTAL_ANG_MOM_LOOP
      end do ORBITAL_ANG_MOM_LOOP
    end do ENSEMBLE_LOOP

    ! allocate URR resonance ensemble realizations
    call tope % alloc_ensemble(n_reals)

    ! transfer scratch URR ensembles to permanent array
    do i_ens = 1, n_reals
      do i_l = 1, tope % NLS(tope % i_urr)
        do i_J = 1, tope % NJS(i_l)
          do i_res = 1, tope % n_lam(i_ens, i_l, i_J)
            tope % urr_resonances(i_ens, i_l, i_J) % E_lam(i_res) &
              & = tope % urr_resonances_tmp(i_ens, i_l, i_J) % E_lam(i_res)
            tope % urr_resonances(i_ens, i_l, i_J) % GN(i_res) &
              & = tope % urr_resonances_tmp(i_ens, i_l, i_J) % GN(i_res)
            tope % urr_resonances(i_ens, i_l, i_J) % GG(i_res) &
              & = tope % urr_resonances_tmp(i_ens, i_l, i_J) % GG(i_res)
            tope % urr_resonances(i_ens, i_l, i_J) % GF(i_res) &
              & = tope % urr_resonances_tmp(i_ens, i_l, i_J) % GF(i_res)
            tope % urr_resonances(i_ens, i_l, i_J) % GX(i_res) &
              & = tope % urr_resonances_tmp(i_ens, i_l, i_J) % GX(i_res)
            tope % urr_resonances(i_ens, i_l, i_J) % GT(i_res) &
              & = tope % urr_resonances_tmp(i_ens, i_l, i_J) % GT(i_res)
          end do
        end do
      end do
    end do

    ! deallocate temporary URR resonance ensemble realizations
    call tope % dealloc_ensemble_tmp(n_reals)

  end subroutine resonance_ensemble

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! POINTWISE_URR generates pointwise energy-cross section data in the URR
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine pointwise_urr(iso, i_nuc, T_K)

    type(Isotope), pointer :: tope => null() ! isotope object pointer
    type(Nuclide), pointer :: nuc => null() ! nuclide object pointer
    type(Resonance) :: res ! resonance object
    type(CrossSection) :: n ! elastic scattering xs object
    type(CrossSection) :: g ! radiative capture xs object
    type(CrossSection) :: f ! fission xs object
    type(CrossSection) :: x ! competitive inelastic scattering xs object
    type(CrossSection) :: t ! total xs object
    integer :: iso       ! isotope index
    integer :: i_nuc     ! nuclide index
    integer :: i_l       ! orbital quantum number
    integer :: i_J       ! total angular momentum quantum number
    integer :: n_res     ! number of contributing l-state resonances
    integer :: n_rrr_res ! number of RRR resonances we need to grab
    integer :: i_low     ! index of lowest-lying resonance
    integer :: i_res     ! resonance counter
    integer :: i_rrr_res ! RRR resonance index
    integer :: i_grid    ! xs energy grid index
    integer :: n_pts     ! xs energy grid point counter
    integer :: i_ES      ! index of current URR tabulated energy
    integer :: iavg      ! index in average cross section
    real(8) :: T_K          ! isotope temperature [K]
    real(8) :: favg         ! average cross section interpolation factor
    real(8) :: fact         ! cross section energy grid interpolation factor
    real(8) :: avg_urr_n_xs ! averaged elastic cross section
    real(8) :: avg_urr_f_xs ! averaged fission cross section
    real(8) :: avg_urr_g_xs ! averaged capture cross section
    real(8) :: avg_urr_x_xs ! averaged competitive inelastic cross section
    real(8) :: dE_trial     ! trial energy spacing for xs grid
    real(8) :: xs_trial     ! trial xs via interpolation between gridpoints
    real(8) :: rel_err      ! relative error between interpolated and exact xs
    logical :: enhance ! refine energy-xs grid?
!$omp threadprivate(tope)

    ! only one realization allowed when using a pointwise representation
    i_real = 1

    tope => isotopes(iso)
    nuc  => nuclides(i_nuc)

    tope % T = T_K

    ! set current energy
    tope % E = tope % EL(tope % i_urr)

    ! allocate pointwise URR cross sections
    call tope % alloc_pointwise_tmp()

    ! enforce xs continuity at RRR-URR energy crossover
    n_pts = 1
    i_grid = binary_search(1.0e6_8 * nuc % energy, nuc % n_grid, tope % E)
    tope % E_tmp(n_pts)  = 1.0e6_8 * nuc % energy(i_grid)
    tope % n_tmp(n_pts) = nuc % elastic(i_grid)
    tope % g_tmp(n_pts) = nuc % absorption(i_grid) &
      & - nuc % fission(i_grid)
    tope % f_tmp(n_pts) = nuc % fission(i_grid)
    tope % t_tmp(n_pts) = nuc % total(i_grid)
    tope % x_tmp(n_pts) = tope % t_tmp(n_pts) &
      & - tope % n_tmp(n_pts) &
      & - tope % g_tmp(n_pts) &
      & - tope % f_tmp(n_pts)
    tope % t_tmp(n_pts) = tope % n_tmp(n_pts) &
      & + tope % g_tmp(n_pts) + tope % f_tmp(n_pts)&
      & + tope % x_tmp(n_pts)

    tope % E = 1.0e6_8 * nuc % energy(i_grid + 1)

    i_ES = 1
    ENERGY_LOOP: do

      dE_trial = max_dE_point_urr
      tope % E = tope % E + dE_trial
      if (tope % E &
           > min(tope % EH(tope % i_urr), max_E_point_urr) + dE_trial)&
           exit
      if (i_ES < tope % NE) then
        if (tope % E > tope % ES(i_ES + 1)) then
          i_ES = i_ES + 1
          write(*,'(A40,ES23.16,A12)') &
            & 'Reconstructing URR xs in', tope % ES(i_ES), ' eV interval'
        end if
      end if
      n_pts = n_pts + 1
      enhance  = .true.

      do while(enhance)

        ! reset xs accumulators
        call flush_sigmas(t, n, g, f, x)

        ! loop over orbital quantum numbers
        LOC_ORBITAL_ANG_MOM_LOOP: do i_l = 1, tope % NLS(tope % i_urr)

          ! set current orbital angular momentum quantum number
          tope % L = i_l - 1

          ! get the number of contributing l-wave resonances for this l
          n_res = n_res_contrib(tope % L)

          ! loop over total angular momentum quantum numbers
          LOC_TOTAL_ANG_MOM_LOOP: do i_J = 1, tope % NJS(i_l)

            ! set current total angular momentum quantum number
            tope % J = tope % AJ(i_l) % data(i_J)

            ! zero the resonance counter
            res % i_res = 0

            ! find the nearest lower resonance
            if (tope % E &
              & < tope % urr_resonances(i_real, i_l, i_J) % E_lam(1)) then
              i_low = 1
            else
              i_low = binary_search(&
                & tope % urr_resonances(i_real, i_l, i_J) % E_lam(:),&
                & tope % n_lam(i_real, i_l, i_J), tope % E)
            end if

            ! loop over the addition of resonances to this ladder
            if (i_low - n_res/2 < 1) then
              ! if we're near the lower end of the URR, need to incorporate
              ! resolved resonance region resonances in order to fix-up
              ! (i.e. smooth out) cross sections at the RRR-URR crossover
              ! energy

              ! if the RRR has resonances with this l-state
              if (i_l <= tope % NLS(tope % i_urr - 1)) then

                ! how many RRR resonances are contributing
                n_rrr_res = abs(i_low - n_res/2) + 1

                ! loop over contributing resolved resonance region resonances
                LOC_RRR_RESONANCES_LOOP: do i_res = n_rrr_res, 1, -1
                  i_rrr_res = rrr_res(iso, i_res, tope % L, tope % J)
                  res % i_res = res % i_res + 1
                  call set_parameters(res, iso, i_rrr_res, i_l, i_J,&
                       tope % i_urr - 1)
                end do LOC_RRR_RESONANCES_LOOP
              end if

              ! loop over contributing unresolved resonance region resonances
              LOC_URR_RESONANCES_LOOP: do i_res = 1, n_res - res % i_res
                res % i_res = res % i_res + 1
                call set_parameters(res, iso, i_res, i_l, i_J, tope % i_urr)
              end do LOC_URR_RESONANCES_LOOP

            else
              ! we're firmly in the URR and can ignore anything going on in
              ! the upper resolved resonance region energies
              LOC_URR_LOOP: do i_res = i_low - n_res/2, i_low + n_res/2 - 1
                res % i_res = res % i_res + 1
                call set_parameters(res, iso, i_res, i_l, i_J, tope % i_urr)
              end do LOC_URR_LOOP
            end if
          end do LOC_TOTAL_ANG_MOM_LOOP
        end do LOC_ORBITAL_ANG_MOM_LOOP

        ! loop over orbital quantum numbers
        ORBITAL_ANG_MOM_LOOP: do i_l = 1, tope % NLS(tope % i_urr)

          ! set current orbital angular momentum quantum number
          tope % L = i_l - 1

          ! get the number of contributing l-wave resonances for this l
          n_res = n_res_contrib(tope % L)

          ! loop over total angular momentum quantum numbers
          TOTAL_ANG_MOM_LOOP: do i_J = 1, tope % NJS(i_l)

            ! set current total angular momentum quantum number
            tope % J = tope % AJ(i_l) % data(i_J)

            ! compute statistical spin factor
            tope % g_J = (TWO * tope % J + ONE) &
                 / (FOUR * tope % SPI(tope % i_urr) + TWO)

            ! loop over resonances localized about e_n
            RESONANCES_LOOP: do i_res = 1, n_res

              res % i_res = i_res
              res % E_lam = tope % local_realization(i_l, i_J) % E_lam(i_res)
              res % Gam_n = tope % local_realization(i_l, i_J) % Gam_n(i_res)
              res % Gam_g = tope % local_realization(i_l, i_J) % Gam_g(i_res)
              res % Gam_f = tope % local_realization(i_l, i_J) % Gam_f(i_res)
              res % Gam_x = tope % local_realization(i_l, i_J) % Gam_x(i_res)
              res % Gam_t = tope % local_realization(i_l, i_J) % Gam_t(i_res)

              call res % calc_xs(iso)
              call accum_resonances(res, t, n, g, f, x)

            end do RESONANCES_LOOP
          end do TOTAL_ANG_MOM_LOOP
        end do ORBITAL_ANG_MOM_LOOP

        ! add potential scattering contribution
        call n % potential_xs(iso)
        call t % potential_xs(iso)

        if (tope % E < 1.0e6_8 * nuc % energy(1)) then
          call fatal_error('URR energy grid extends below available pointwise&
            & data')
        elseif (tope % E > 1.0e6_8 * nuc % energy(nuc % n_grid)) then
          call fatal_error('URR energy grid extends above available pointwise&
            & data')
        else
          i_grid = binary_search(1.0e6_8 * nuc % energy, nuc % n_grid, &
            & tope % E)
        end if

        ! check for rare case where two energy points are the same
        if (nuc % energy(i_grid) == nuc % energy(i_grid + 1)) &
          & i_grid = i_grid + 1

        ! calculate xs energy grid interpolation factor
        fact = interp_factor(tope % E, 1.0e6_8 * nuc % energy(i_grid), &
          & 1.0e6_8 * nuc % energy(i_grid + 1), tope % INT)

        ! calculate evaluator-supplied backgrounds at the current energy
        call interp_ace_background(iso, i_nuc, n_pts, fact, i_grid)

        ! interpret MF3 data according to ENDF-6 LSSF flag:
        ! MF3 contains background xs, add to MF2 resonance contributions
        if (tope % LSSF == 0) then

          ! add resonance xs component to background
          call add_mf3_background(iso, i_nuc, n_pts, t, n, g, f, x)

        ! multipy the self-shielding factors by the infinite-dilute xs
        elseif (tope % LSSF == 1) then

          if (tope % E > tope % Eavg(tope % nEavg)) then
            iavg = tope % nEavg - 1
          else
            iavg = binary_search(tope % Eavg, tope % nEavg, tope % E)
          end if

          favg = interp_factor(tope % E, &
            & tope % Eavg(iavg), tope % Eavg(iavg + 1), tope % INT)

          call interp_avg_urr_xs(favg, iso, iavg, &
            & avg_urr_n_xs, avg_urr_f_xs, avg_urr_g_xs, avg_urr_x_xs)

          ! competitive xs
          if (avg_urr_x_xs > ZERO) then
            if (tope % E < tope % E_ex2 .and. competitive) then
              ! self-shielded treatment of competitive inelastic xs
              tope % x_tmp(n_pts) = x % xs / avg_urr_x_xs &
                & * tope % x_tmp(n_pts)
            else
              ! infinite-dilute treatment of competitive inelastic xs
              tope % x_tmp(n_pts) = tope % x_tmp(n_pts)
            end if
          else
            ! use background competitive inelastic cross section, as is
            tope % x_tmp(n_pts) = tope % x_tmp(n_pts)
          end if

          ! elastic scattering xs
          tope % n_tmp(n_pts) = n % xs / avg_urr_n_xs * tope % n_tmp(n_pts)

          ! set negative elastic xs and competitive xs to zero
          if (tope % n_tmp(n_pts) < ZERO) tope % n_tmp(n_pts) = ZERO
          if (tope % x_tmp(n_pts) < ZERO) tope % x_tmp(n_pts) = ZERO

          ! radiative capture xs
          if (avg_urr_g_xs > ZERO) then
            tope % g_tmp(n_pts) = g % xs / avg_urr_g_xs * tope % g_tmp(n_pts)
          else
            ! use background capture cross section, as is
            tope % g_tmp(n_pts) = tope % g_tmp(n_pts)
          end if

          ! fission xs
          if (avg_urr_f_xs > ZERO) then
            tope % f_tmp(n_pts) = f % xs / avg_urr_f_xs * tope % f_tmp(n_pts)
          else
            tope % f_tmp(n_pts) = tope % f_tmp(n_pts)
          end if

          tope % t_tmp(n_pts) = tope % n_tmp(n_pts) + tope % g_tmp(n_pts)&
            & + tope % f_tmp(n_pts) + tope % x_tmp(n_pts)

        else
          call fatal_error('ENDF-6 LSSF not allowed - must be 0 or 1.')
        end if

        xs_trial = HALF &
          & * (tope % t_tmp(n_pts) + tope % t_tmp(n_pts - 1))
        dE_trial = HALF * dE_trial
        tope % E  = tope % E - dE_trial

        ! reset xs accumulators
        call flush_sigmas(t, n, g, f, x)

        ! loop over orbital quantum numbers
        ORBITAL_ANG_MOM_LOOPb: do i_l = 1, tope % NLS(tope % i_urr)

          ! set current orbital angular momentum quantum number
          tope % L = i_l - 1

          ! get the number of contributing l-wave resonances for this l
          n_res = n_res_contrib(tope % L)

          ! loop over total angular momentum quantum numbers
          TOTAL_ANG_MOM_LOOPb: do i_J = 1, tope % NJS(i_l)

            ! set current total angular momentum quantum number
            tope % J = tope % AJ(i_l) % data(i_J)

            ! compute statistical spin factor
            tope % g_J = (TWO * tope % J + ONE) &
                 / (FOUR * tope % SPI(tope % i_urr) + TWO)

            ! loop over resonances localized about e_n
            RESONANCES_LOOPb: do i_res = 1, n_res

              res % i_res = i_res
              res % E_lam = tope % local_realization(i_l, i_J) % E_lam(i_res)
              res % Gam_n = tope % local_realization(i_l, i_J) % Gam_n(i_res)
              res % Gam_g = tope % local_realization(i_l, i_J) % Gam_g(i_res)
              res % Gam_f = tope % local_realization(i_l, i_J) % Gam_f(i_res)
              res % Gam_x = tope % local_realization(i_l, i_J) % Gam_x(i_res)
              res % Gam_t = tope % local_realization(i_l, i_J) % Gam_t(i_res)

              call res % calc_xs(iso)
              call accum_resonances(res, t, n, g, f, x)

            end do RESONANCES_LOOPb
          end do TOTAL_ANG_MOM_LOOPb
        end do ORBITAL_ANG_MOM_LOOPb

        ! add potential scattering contribution
        call n % potential_xs(iso)
        call t % potential_xs(iso)

        if (tope % E < 1.0e6_8 * nuc % energy(1)) then
          call fatal_error('URR energy grid extends below available pointwise&
            & data')
        elseif (tope % E > 1.0e6_8 * nuc % energy(nuc % n_grid)) then
          call fatal_error('URR energy grid extends above available pointwise&
            & data')
        else
          i_grid = binary_search(1.0e6_8 * nuc % energy, nuc % n_grid, &
            & tope % E)
        end if

        ! check for rare case where two energy points are the same
        if (nuc % energy(i_grid) == nuc % energy(i_grid + 1)) &
          & i_grid = i_grid + 1

        ! calculate xs energy grid interpolation factor
        fact = interp_factor(tope % E, 1.0e6_8 * nuc % energy(i_grid), &
          & 1.0e6_8 * nuc % energy(i_grid + 1), tope % INT)

        ! calculate evaluator-supplied backgrounds at the current energy
        call interp_ace_background(iso, i_nuc, n_pts, fact, i_grid)

        ! interpret MF3 data according to ENDF-6 LSSF flag:
        ! MF3 contains background xs, add to MF2 resonance contributions
        if (tope % LSSF == 0) then

          ! add resonance xs component to background
          call add_mf3_background(iso, i_nuc, n_pts, t, n, g, f, x)

        elseif (tope % LSSF == 1) then
          ! multipy the self-shielding factors by the infinite-dilute xs

          if (tope % E > tope % Eavg(tope % nEavg)) then
            iavg = tope % nEavg - 1
          else
            iavg = binary_search(tope % Eavg, tope % nEavg, tope % E)
          end if

          favg = interp_factor(tope % E, &
            & tope % Eavg(iavg), tope % Eavg(iavg + 1), tope % INT)

          call interp_avg_urr_xs(favg, iso, iavg, &
            & avg_urr_n_xs, avg_urr_f_xs, avg_urr_g_xs, avg_urr_x_xs)

          ! competitive xs
          if (avg_urr_x_xs > ZERO) then
            if (tope % E < tope % E_ex2 .and. competitive) then
              ! self-shielded treatment of competitive inelastic xs
              tope % x_tmp(n_pts) = x % xs / avg_urr_x_xs &
                & * tope % x_tmp(n_pts)
            else
              ! infinite-dilute treatment of competitive inelastic xs
              tope % x_tmp(n_pts) = tope % x_tmp(n_pts)
            end if
          else
            ! use background competitive inelastic cross section, as is
            tope % x_tmp(n_pts) = tope % x_tmp(n_pts)
          end if

          ! elastic scattering xs
          tope % n_tmp(n_pts) = n % xs / avg_urr_n_xs * tope % n_tmp(n_pts)

          ! set negative elastic xs and competitive xs to zero
          if (tope % n_tmp(n_pts) < ZERO) tope % n_tmp(n_pts) = ZERO
          if (tope % x_tmp(n_pts) < ZERO) tope % x_tmp(n_pts) = ZERO

          ! radiative capture xs
          if (avg_urr_g_xs > ZERO) then
            tope % g_tmp(n_pts) = g % xs / avg_urr_g_xs * tope % g_tmp(n_pts)
          else
            ! use background capture cross section, as is
            tope % g_tmp(n_pts) = tope % g_tmp(n_pts)
          end if

          ! fission xs
          if (avg_urr_f_xs > ZERO) then
            tope % f_tmp(n_pts) = f % xs / avg_urr_f_xs * tope % f_tmp(n_pts)
          else
            tope % f_tmp(n_pts) = tope % f_tmp(n_pts)
          end if

          tope % t_tmp(n_pts) = tope % n_tmp(n_pts) + tope % g_tmp(n_pts)&
            & + tope % f_tmp(n_pts) + tope % x_tmp(n_pts)

        else
          call fatal_error('ENDF-6 LSSF not allowed - must be 0 or 1.')
        end if

        rel_err = abs(xs_trial - tope % t_tmp(n_pts)) &
          & / tope % t_tmp(n_pts)
        if (rel_err < tol_point_urr .or. dE_trial < min_dE_point_urr) then
          enhance = .false.
        end if
      end do

      ! add energy point to grid
      tope % E = tope % E + dE_trial
      tope % E_tmp(n_pts) = tope % E

      ! reset xs accumulators
      call flush_sigmas(t, n, g, f, x)

      ! loop over orbital quantum numbers
      ORBITAL_ANG_MOM_LOOPc: do i_l = 1, tope % NLS(tope % i_urr)

        ! set current orbital angular momentum quantum number
        tope % L = i_l - 1

        ! get the number of contributing l-wave resonances for this l
        n_res = n_res_contrib(tope % L)

        ! loop over total angular momentum quantum numbers
        TOTAL_ANG_MOM_LOOPc: do i_J = 1, tope % NJS(i_l)

          ! set current total angular momentum quantum number
          tope % J = tope % AJ(i_l) % data(i_J)

          ! compute statistical spin factor
          tope % g_J = (TWO * tope % J + ONE) &
               / (FOUR * tope % SPI(tope % i_urr) + TWO)

          ! loop over resonances localized about e_n
          RESONANCES_LOOPc: do i_res = 1, n_res

            res % i_res = i_res
            res % E_lam = tope % local_realization(i_l, i_J) % E_lam(i_res)
            res % Gam_n = tope % local_realization(i_l, i_J) % Gam_n(i_res)
            res % Gam_g = tope % local_realization(i_l, i_J) % Gam_g(i_res)
            res % Gam_f = tope % local_realization(i_l, i_J) % Gam_f(i_res)
            res % Gam_x = tope % local_realization(i_l, i_J) % Gam_x(i_res)
            res % Gam_t = tope % local_realization(i_l, i_J) % Gam_t(i_res)

            call res % calc_xs(iso)
            call accum_resonances(res, t, n, g, f, x)

          end do RESONANCES_LOOPc
        end do TOTAL_ANG_MOM_LOOPc
      end do ORBITAL_ANG_MOM_LOOPc

      ! add potential scattering contribution
      call n % potential_xs(iso)
      call t % potential_xs(iso)

      if (tope % E < 1.0e6_8 * nuc % energy(1)) then
        call fatal_error('URR energy grid extends below available pointwise&
          & data')
      elseif (tope % E > 1.0e6_8 * nuc % energy(nuc % n_grid)) then
        call fatal_error('URR energy grid extends above available pointwise&
          & data')
      else
        i_grid = binary_search(1.0e6_8 * nuc % energy, nuc % n_grid, tope % E)
      end if

      ! check for rare case where two energy points are the same
      if (nuc % energy(i_grid) == nuc % energy(i_grid + 1)) i_grid = i_grid + 1

      ! calculate xs energy grid interpolation factor
      fact = interp_factor(tope % E, 1.0e6_8 * nuc % energy(i_grid), &
        & 1.0e6_8 * nuc % energy(i_grid + 1), tope % INT)

      ! calculate evaluator-supplied backgrounds at the current energy
      call interp_ace_background(iso, i_nuc, n_pts, fact, i_grid)

      ! interpret MF3 data according to ENDF-6 LSSF flag:
      ! MF3 contains background xs values, add to MF2 resonance contributions
      if (tope % LSSF == 0) then

        ! add resonance xs component to background
        call add_mf3_background(iso, i_nuc, n_pts, t, n, g, f, x)

      elseif (tope % LSSF == 1) then
        ! multipy the self-shielding factors by the average (infinite-dilute)
        ! cross sections

        if (tope % E > tope % Eavg(tope % nEavg)) then
          iavg = tope % nEavg - 1
        else
          iavg = binary_search(tope % Eavg, tope % nEavg, tope % E)
        end if

        favg = interp_factor(tope % E, &
          & tope % Eavg(iavg), tope % Eavg(iavg + 1), tope % INT)

        call interp_avg_urr_xs(favg, iso, iavg, &
          & avg_urr_n_xs, avg_urr_f_xs, avg_urr_g_xs, avg_urr_x_xs)

        ! competitive xs
        if (avg_urr_x_xs > ZERO) then
          if (tope % E < tope % E_ex2 .and. competitive) then
            ! self-shielded treatment of competitive inelastic xs
            tope % x_tmp(n_pts) = x % xs / avg_urr_x_xs &
              & * tope % x_tmp(n_pts)
          else
            ! infinite-dilute treatment of competitive inelastic xs
            tope % x_tmp(n_pts) = tope % x_tmp(n_pts)
          end if
        else
          ! use background competitive inelastic cross section, as is
          tope % x_tmp(n_pts) = tope % x_tmp(n_pts)
        end if

        ! elastic scattering xs
        tope % n_tmp(n_pts) = n % xs / avg_urr_n_xs * tope % n_tmp(n_pts)

        ! set negative elastic xs and competitive xs to zero
        if (tope % n_tmp(n_pts) < ZERO) tope % n_tmp(n_pts) = ZERO
        if (tope % x_tmp(n_pts) < ZERO) tope % x_tmp(n_pts) = ZERO

        ! radiative capture xs
        if (avg_urr_g_xs > ZERO) then
          tope % g_tmp(n_pts) = g % xs / avg_urr_g_xs * tope % g_tmp(n_pts)
        else
          ! use background capture cross section, as is
          tope % g_tmp(n_pts) = tope % g_tmp(n_pts)
        end if

        ! fission xs
        if (avg_urr_f_xs > ZERO) then
          tope % f_tmp(n_pts) = f % xs / avg_urr_f_xs * tope % f_tmp(n_pts)
        else
          tope % f_tmp(n_pts) = tope % f_tmp(n_pts)
        end if

        tope % t_tmp(n_pts) = tope % n_tmp(n_pts) + tope % g_tmp(n_pts)&
          & + tope % f_tmp(n_pts) + tope % x_tmp(n_pts)

      else
        call fatal_error('Self-shielding flag (LSSF) not allowed -&
          & must be 0 or 1.')
      end if

    end do ENERGY_LOOP

    ! pass temporary, pre-allocated energy-xs vectors to dynamic vectors of
    ! the proper length
    call tope % alloc_pointwise(n_pts)
    tope % urr_E = tope % E_tmp(1:n_pts)
    tope % urr_n = tope % n_tmp(1:n_pts)
    tope % urr_g = tope % g_tmp(1:n_pts)
    tope % urr_f = tope % f_tmp(1:n_pts)
    tope % urr_x = tope % x_tmp(1:n_pts)
    tope % urr_t = tope % t_tmp(1:n_pts)
    call tope % dealloc_pointwise_tmp()

  end subroutine pointwise_urr

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! CALC_URR_XS_OTF calculates unresolved resonance region cross sections, at a
! single energy, on-the-fly from resonance parameters for a single realization
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine calc_urr_xs_otf(iso, i_nuc, E, T_K)

    type(Isotope), pointer, save :: tope => null() ! isotope object pointer
    type(Nuclide), pointer :: nuc => null() ! nuclide object pointer
    type(Resonance) :: res ! resonance object
    type(CrossSection) :: n ! elastic scattering xs object
    type(CrossSection) :: g ! radiative capture xs object
    type(CrossSection) :: f ! fission xs object
    type(CrossSection) :: x ! competitive inelastic scattering xs object
    type(CrossSection) :: t ! total xs object
    integer :: iso       ! isotope index
    integer :: i_nuc     ! nuclide index
    integer :: iavg      ! average cross section index
    integer :: i_l       ! orbital quantum number
    integer :: i_J       ! total angular momentum quantum number
    integer :: n_res     ! number of contributing l-state resonances
    integer :: n_rrr_res ! number of RRR resonances we need to grab
    integer :: i_low     ! index of lowest-lying resonance
    integer :: i_res     ! resonance counter
    integer :: i_rrr_res ! RRR resonance index
    real(8) :: E            ! neutron energy
    real(8) :: T_K          ! isotope temperature [K]
    real(8) :: favg         ! average cross section interpolation factor
    real(8) :: avg_urr_n_xs ! averaged elastic cross section
    real(8) :: avg_urr_f_xs ! averaged fission cross section
    real(8) :: avg_urr_g_xs ! averaged capture cross section
    real(8) :: avg_urr_x_xs ! averaged competitive reaction cross section
    real(8) :: capture_xs   ! radiative capture cross section
    real(8) :: inelastic_xs ! competitive inelastic scattering cross section
!$omp threadprivate(tope, nuc)

    tope => isotopes(iso)
    nuc  => nuclides(i_nuc)

    ! set current temperature
    tope % T = T_K

    ! set current energy and interpolation factor
    tope % E = E

    ! reset xs accumulators
    call flush_sigmas(t, n, g, f, x)

    if (i_real_user == 0) then
      ! randomly select which realization to use
      i_real = ceiling(prn() * n_reals)
      if (i_real == 0) call fatal_error('i_real is sampled to be 0')
      if (i_real == n_reals + 1) &
        & call fatal_error('i_real is sampled to be > n_reals')
    else
      i_real = i_real_user
    end if

    ! loop over orbital quantum numbers
    LOC_ORBITAL_ANG_MOM_LOOP: do i_l = 1, tope % NLS(tope % i_urr)

      ! set current orbital angular momentum quantum number
      tope % L = i_l - 1

      ! get the number of contributing l-wave resonances for this l
      n_res = n_res_contrib(tope % L)

      ! loop over total angular momentum quantum numbers
      LOC_TOTAL_ANG_MOM_LOOP: do i_J = 1, tope % NJS(i_l)

        ! set current total angular momentum quantum number
        tope % J = tope % AJ(i_l) % data(i_J)

        ! zero the resonance counter
        res % i_res = 0

        ! find the nearest lower resonance
        if (tope % E &
          & < tope % urr_resonances(i_real, i_l, i_J) % E_lam(1)) then
          i_low = 1
        else
          i_low = binary_search(&
            & tope % urr_resonances(i_real, i_l, i_J) % E_lam(:),&
            & tope % n_lam(i_real, i_l, i_J), tope % E)
        end if

        ! loop over the addition of resonances to this ladder
        if (i_low - n_res/2 < 1) then
          ! if we're near the lower end of the URR, need to incorporate
          ! resolved resonance region resonances in order to fix-up
          ! (i.e. smooth out) cross sections at the RRR-URR crossover
          ! energy

          ! if the RRR has resonances with this l-state
          if (i_l <= tope % NLS(tope % i_urr - 1)) then

            ! how many RRR resonances are contributing
            n_rrr_res = abs(i_low - n_res/2) + 1

            ! loop over contributing resolved resonance region resonances
            LOC_RRR_RESONANCES_LOOP: do i_res = n_rrr_res, 1, -1
              i_rrr_res = rrr_res(iso, i_res, tope % L, tope % J)
              res % i_res = res % i_res + 1
              call set_parameters(res, iso, i_rrr_res, i_l, i_J,&
                   tope % i_urr - 1)
            end do LOC_RRR_RESONANCES_LOOP
          end if

          ! loop over contributing unresolved resonance region resonances
          LOC_URR_RESONANCES_LOOP: do i_res = 1, n_res - res % i_res
            res % i_res = res % i_res + 1
            call set_parameters(res, iso, i_res, i_l, i_J, tope % i_urr)
          end do LOC_URR_RESONANCES_LOOP

        else
          ! we're firmly in the URR and can ignore anything going on in
          ! the upper resolved resonance region energies
          LOC_URR_LOOP: do i_res = i_low - n_res/2, i_low + n_res/2 - 1
            res % i_res = res % i_res + 1
            call set_parameters(res, iso, i_res, i_l, i_J, tope % i_urr)
          end do LOC_URR_LOOP
        end if
      end do LOC_TOTAL_ANG_MOM_LOOP
    end do LOC_ORBITAL_ANG_MOM_LOOP

    ! loop over orbital quantum numbers
    ORBITAL_ANG_MOM_LOOP: do i_l = 1, tope % NLS(tope % i_urr)

      ! set current orbital angular momentum quantum number
      tope % L = i_l - 1

      ! get the number of contributing l-wave resonances for this l
      n_res = n_res_contrib(tope % L)

      ! loop over total angular momentum quantum numbers
      TOTAL_ANG_MOM_LOOP: do i_J = 1, tope % NJS(i_l)

        ! set current total angular momentum quantum number
        tope % J = tope % AJ(i_l) % data(i_J)

        ! compute statistical spin factor
        tope % g_J = (TWO * tope % J + ONE) &
             / (FOUR * tope % SPI(tope % i_urr) + TWO)

        ! loop over resonances localized about e_n
        RESONANCES_LOOP: do i_res = 1, n_res

          res % i_res = i_res
          res % E_lam = tope % local_realization(i_l, i_J) % E_lam(i_res)
          res % Gam_n = tope % local_realization(i_l, i_J) % Gam_n(i_res)
          res % Gam_g = tope % local_realization(i_l, i_J) % Gam_g(i_res)
          res % Gam_f = tope % local_realization(i_l, i_J) % Gam_f(i_res)
          res % Gam_x = tope % local_realization(i_l, i_J) % Gam_x(i_res)
          res % Gam_t = tope % local_realization(i_l, i_J) % Gam_t(i_res)

          call res % calc_xs(iso)
          call accum_resonances(res, t, n, g, f, x)

        end do RESONANCES_LOOP
      end do TOTAL_ANG_MOM_LOOP
    end do ORBITAL_ANG_MOM_LOOP

    ! add potential scattering contribution
    call n % potential_xs(iso)
    call t % potential_xs(iso)

    ! interpret MF3 data according to ENDF-6 LSSF flag:
    ! MF3 contains background xs, add to MF2 resonance contributions
    if (tope % LSSF == 0) then

      ! add resonance xs component to background
      call add_mf3_background(iso, i_nuc, NONE, t, n, g, f, x)

    else if (tope % LSSF == 1) then
      ! multipy the self-shielding factors by the infinite-dilute xs

      ! tabulated unresolved resonance parameters interpolation factor
      iavg = binary_search(tope % Eavg, tope % nEavg, tope % E)

      favg = interp_factor(tope % E, &
        & tope % Eavg(iavg), tope % Eavg(iavg + 1), tope % INT)

      ! interpolate averaged, infinite-dilute URR cross sections
      call interp_avg_urr_xs(favg, iso, iavg, &
        & avg_urr_n_xs, avg_urr_f_xs, avg_urr_g_xs, avg_urr_x_xs)

      ! competitive xs
      if (avg_urr_x_xs > ZERO) then
        if (tope % E < tope % E_ex2 .and. competitive) then
          ! self-shielded treatment of competitive inelastic xs
          inelastic_xs = x % xs / avg_urr_x_xs &
            & * (micro_xs(i_nuc) % total &
            & -  micro_xs(i_nuc) % absorption &
            & -  micro_xs(i_nuc) % elastic)
        else
          ! infinite-dilute treatment of competitive inelastic xs
          inelastic_xs = micro_xs(i_nuc) % total &
            & - micro_xs(i_nuc) % absorption &
            & - micro_xs(i_nuc) % elastic
        end if
      else
        ! use background competitive inelastic cross section, as is
        inelastic_xs = micro_xs(i_nuc) % total &
          & - micro_xs(i_nuc) % absorption &
          & - micro_xs(i_nuc) % elastic
      end if

      ! elastic scattering xs
      micro_xs(i_nuc) % elastic = n % xs / avg_urr_n_xs &
        & * micro_xs(i_nuc) % elastic

      ! set negative elastic xs and competitive xs to zero
      if (micro_xs(i_nuc) % elastic < ZERO) micro_xs(i_nuc) % elastic = ZERO
      if (inelastic_xs < ZERO) inelastic_xs = ZERO

      ! radiative capture xs
      if (avg_urr_g_xs > ZERO) then
        capture_xs = g % xs / avg_urr_g_xs &
          & * (micro_xs(i_nuc) % absorption - micro_xs(i_nuc) % fission)
      else
        ! use background capture cross section, as is
        capture_xs = micro_xs(i_nuc) % absorption - micro_xs(i_nuc) % fission
      end if

      ! fission xs
      if (avg_urr_f_xs > ZERO) then
        micro_xs(i_nuc) % fission = f % xs / avg_urr_f_xs &
          & * micro_xs(i_nuc) % fission
      else
        ! use background fission cross section, as is
        micro_xs(i_nuc) % fission = micro_xs(i_nuc) % fission
      end if

      ! absorption xs
      micro_xs(i_nuc) % absorption = micro_xs(i_nuc) % fission + capture_xs

      ! total xs
      micro_xs(i_nuc) % total = micro_xs(i_nuc) % elastic &
        & + micro_xs(i_nuc) % absorption + inelastic_xs

    else
      call fatal_error('ENDF-6 LSSF not allowed - must be 0 or 1.')
    end if

    ! Determine nu-fission cross section
    if (tope % fissionable) then
      micro_xs(i_nuc) % nu_fission = nu_total(nuc, E / 1.0e6_8) &
        & * micro_xs(i_nuc) % fission
    end if

  end subroutine calc_urr_xs_otf

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! PROB_TABLES computes probability tables for the URR
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine prob_tables(iso)

    type(Isotope), pointer :: tope => null() ! isotope object pointer
    type(Nuclide), pointer :: nuc => null() ! nuclide object pointer
    type(ProbabilityTable), pointer :: ptable => null() ! prob. table pointer
    type(Resonance) :: res ! resonance object
    character(6) :: zaid_str ! ENDF-6 MAT number as a string
    integer :: i_b    ! batch index
    integer :: i_band ! probability band index
    integer :: i_E    ! energy grid index
    integer :: i_grid ! File 3 energy grid index
    integer :: i_h    ! history index
    integer :: i_J    ! total angular momentum quantum number index
    integer :: i_l    ! orbital quantum number index
    integer :: i_nuc  ! nuclide index
    integer :: i_r    ! resonance index
    integer :: i_T    ! temperature index
    integer :: i_temp ! temperature index
    integer :: iso    ! isotope index
    integer :: n_res  ! number of resonances to include for a given l-wave
    integer :: avg_unit = 98 ! avg xs output file unit
    integer :: tab_unit = 99 ! tables output file unit
    real(8) :: E        ! neutron lab energy [eV]
    real(8) :: fact     ! File 3 energy grid interpolation factor
    real(8) :: xs_t_min ! min realized total xs
    real(8) :: xs_t_max ! max realized total xs

    xs_t_min = 1.0e6_8
    xs_t_max = ZERO

    tope => isotopes(iso)

    write(zaid_str, '(I6)') tope % ZAI
    if (master)&
      write(*, '(A80)') 'Generating probability tables for ZAID = '//&
      trim(adjustl(zaid_str))
    if (write_urr_tables)&
      open(unit = tab_unit, file = trim(adjustl(zaid_str))//'-urr-tables.dat')
    if (write_avg_urr_xs) then
      open(unit = avg_unit, file = trim(adjustl(zaid_str))//'-avg-urr-xs.dat')
      write(avg_unit, '(I6)') tope % ZAI
      write(avg_unit, '(I1)') represent_params
      write(avg_unit, '(ES13.6)') tol_avg_urr
      write(avg_unit, '(A13,A13,A13,A13,A13,A13)') &
        & 'energy [eV]', 'total', 'elastic', 'capture', 'fission', 'competitive'
    end if

    ! loop over energy mesh
    ENERGY_LOOP: do i_E = 1, tope % ntabs

      tope % E = tope % Etabs(i_E)
      E = tope % E

      ! reset accumulator of statistics
      call tope % flush_ptable_stats(i_E)

      i_b = 0

      ! loop over batches until convergence
      BATCH_LOOP: do

        i_b = i_b + 1

        ! reset batch accumulators
        call tope % flush_batches()

        ! loop over realizations
        HISTORY_LOOP: do i_h = 1, histories_avg_urr

          ! reset accumulator of histories
          call tope % flush_histories()

          ! loop over orbital quantum numbers
          LOC_ORBITAL_ANG_MOM_LOOP: do i_l = 1, tope % NLS(tope % i_urr)

            ! set current orbital angular momentum quantum number
            tope % L = i_l - 1

            ! get the number of contributing l-wave resonances for this l
            n_res = n_res_contrib(tope % L)

            ! loop over total angular momentum quantum numbers
            LOC_TOTAL_ANG_MOM_LOOP: do i_J = 1, tope % NJS(i_l)

              ! set current total angular momentum quantum number
              tope % J = tope % AJ(i_l) % data(i_J)

              ! set current partial width degrees of freedom
              tope % AMUX = int(tope % DOFX(i_l) % data(i_J))
              tope % AMUN = int(tope % DOFN(i_l) % data(i_J))
              tope % AMUG = int(tope % DOFG(i_l) % data(i_J))
              tope % AMUF = int(tope % DOFF(i_l) % data(i_J))

              ! zero the resonance counter
              res % i_res = 0

              ! set mean URR parameters to neutron energy
              call set_mean_parameters(iso, E, i_l, i_J)

              ! sample unresolved resonance parameters for this spin
              ! sequence, at this energy
              call res % sample_parameters(iso, i_l, i_J)

              ! loop over the addition of resonances to this ladder
              LOC_RESONANCES_LOOP: do i_r = 1, n_res

                if (represent_params == E_RESONANCE) then
                  ! interpolate mean URR parameters to current resonance energy
                  call set_mean_parameters(iso, res % E_lam, i_l, i_J)
                end if

                res % i_res = i_r

                ! sample unresolved resonance parameters for this spin
                ! sequence, at this energy
                call res % sample_parameters(iso, i_l, i_J)

              end do LOC_RESONANCES_LOOP
            end do LOC_TOTAL_ANG_MOM_LOOP
          end do LOC_ORBITAL_ANG_MOM_LOOP

          ! loop over orbital quantum numbers
          ORBITAL_ANG_MOM_LOOP: do i_l = 1, tope % NLS(tope % i_urr)

            ! set current orbital angular momentum quantum number
            tope % L = i_l - 1

            ! get the number of contributing l-wave resonances for this l
            n_res = n_res_contrib(tope % L)

            ! loop over total angular momentum quantum numbers
            TOTAL_ANG_MOM_LOOP: do i_J = 1, tope % NJS(i_l)

              ! set current total angular momentum quantum number
              tope % J = tope % AJ(i_l) % data(i_J)

              ! compute statistical spin factor
              tope % g_J = (TWO * tope % J + ONE) &
                   / (FOUR * tope % SPI(tope % i_urr) + TWO)

              ! loop over resonances localized about e_n
              RESONANCES_LOOP: do i_r = 1, n_res

                res % i_res = i_r
                res % E_lam = tope % local_realization(i_l, i_J) % E_lam(i_r)
                res % Gam_n = tope % local_realization(i_l, i_J) % Gam_n(i_r)
                res % Gam_g = tope % local_realization(i_l, i_J) % Gam_g(i_r)
                res % Gam_f = tope % local_realization(i_l, i_J) % Gam_f(i_r)
                res % Gam_x = tope % local_realization(i_l, i_J) % Gam_x(i_r)
                res % Gam_t = tope % local_realization(i_l, i_J) % Gam_t(i_r)

                TEMPERATURES_LOOP: do i_T = 1, tope % nT

                  ! set current temperature
                  tope % T = tope % Tlist % get_item(i_T)

                  ! calculate the contribution to the partial cross sections,
                  ! at this energy, from an additional resonance
                  call res % calc_xs(iso)

                  ! add this contribution to the accumulated partial cross
                  ! section values built up from all resonances
                  ! TODO: move t outside of loop
                  call tope % res_contrib(res, i_E, i_T)

                end do TEMPERATURES_LOOP
              end do RESONANCES_LOOP
            end do TOTAL_ANG_MOM_LOOP
          end do ORBITAL_ANG_MOM_LOOP

          do i_T = 1, tope % nT

            i_nuc = tope % inuclist % get_item(i_T)
            nuc  => nuclides(i_nuc)

            ! add potential scattering contribution
            call tope % prob_tables(i_E, i_T) % avg_n % potential_xs(iso)
            call tope % prob_tables(i_E, i_T) % avg_t % potential_xs(iso)

            ! set negative elastic xs to zero
            if (tope % prob_tables(i_E, i_T) % avg_n % xs < ZERO) then
              tope % prob_tables(i_E, i_T) % avg_t % xs &
                & = tope % prob_tables(i_E, i_T) % avg_t % xs &
                & + abs(tope % prob_tables(i_E, i_T) % avg_n % xs)
              tope % prob_tables(i_E, i_T) % avg_n % xs = ZERO
            end if

            ! add File 3 fission reaction contribution
            if (tope % prob_tables(i_E, i_T) % avg_f % xs > ZERO) then
              continue
            else
              tope % prob_tables(i_E, i_T) % avg_f % xs = ZERO
              if (background == ENDFFILE .and. allocated(tope % MF3_f_e)) then
                if (E < tope % MF3_f_e(1)) then
                  tope % prob_tables(i_E, i_T) % avg_f % xs = ZERO
                else
                  i_grid = binary_search(tope % MF3_f_e, size(tope % MF3_f_e),&
                    & E)
                  if (tope % INT == LINEAR_LINEAR &
                    & .or. tope % MF3_f(i_grid) > XS_CUTOFF) then
                    fact = interp_factor(E, tope % MF3_f_e(i_grid), &
                      & tope % MF3_f_e(i_grid + 1), tope % INT)
                    tope % prob_tables(i_E, i_T) % avg_f % xs &
                      & = interpolator(fact, tope % MF3_f(i_grid), &
                      & tope % MF3_f(i_grid + 1), tope % INT)
                    tope % prob_tables(i_E, i_T) % avg_t % xs &
                      & = tope % prob_tables(i_E, i_T) % avg_t % xs &
                      & + tope % prob_tables(i_E, i_T) % avg_f % xs
                  else
                    tope % prob_tables(i_E, i_T) % avg_f % xs = ZERO
                  end if
                end if
              else if (background == ACEFILE) then
                if (E < 1.0e6_8 * nuc % energy(1)) then
                  tope % prob_tables(i_E, i_T) % avg_f % xs = ZERO
                else
                  i_grid = binary_search(1.0e6_8 * nuc % energy, &
                    & size(nuc % energy), E)
                  if (tope % INT == LINEAR_LINEAR &
                    & .or. nuc % fission(i_grid) > XS_CUTOFF) then
                    fact = interp_factor(E, 1.0e6_8 * nuc % energy(i_grid), &
                      & 1.0e6_8 * nuc % energy(i_grid + 1), tope % INT)
                    tope % prob_tables(i_E, i_T) % avg_f % xs &
                      & = interpolator(fact, nuc % fission(i_grid), &
                      & nuc % fission(i_grid + 1), tope % INT)
                    tope % prob_tables(i_E, i_T) % avg_t % xs &
                      & = tope % prob_tables(i_E, i_T) % avg_t % xs &
                      & + tope % prob_tables(i_E, i_T) % avg_f % xs
                  else
                    tope % prob_tables(i_E, i_T) % avg_f % xs = ZERO
                  end if
                end if
              end if
            end if

            ! add File 3 infinite-dilute competitive reaction contribution
            if (tope % E < tope % E_ex2 .and. competitive) then
              continue
            else
              tope % prob_tables(i_E, i_T) % avg_x % xs = ZERO
              if (background == ENDFFILE .and. allocated(tope % MF3_x_e)) then
                if (E < tope % MF3_x_e(1)) then
                  tope % prob_tables(i_E, i_T) % avg_x % xs = ZERO
                else
                  i_grid = binary_search(tope % MF3_x_e, size(tope % MF3_x_e),&
                    & E)
                  if (tope % INT == LINEAR_LINEAR &
                    & .or. tope % MF3_x(i_grid) > XS_CUTOFF) then
                    fact = interp_factor(E, tope % MF3_x_e(i_grid), &
                      & tope % MF3_x_e(i_grid + 1), tope % INT)
                    tope % prob_tables(i_E, i_T) % avg_t % xs &
                      & = tope % prob_tables(i_E, i_T) % avg_t % xs &
                      & - tope % prob_tables(i_E, i_T) % avg_x % xs
                    tope % prob_tables(i_E, i_T) % avg_x % xs &
                      & = interpolator(fact, tope % MF3_x(i_grid), &
                      & tope % MF3_x(i_grid + 1), tope % INT)
                    tope % prob_tables(i_E, i_T) % avg_t % xs &
                      & = tope % prob_tables(i_E, i_T) % avg_t % xs &
                      & + tope % prob_tables(i_E, i_T) % avg_x % xs
                  else
                    tope % prob_tables(i_E, i_T) % avg_t % xs &
                      & = tope % prob_tables(i_E, i_T) % avg_t % xs &
                      & - tope % prob_tables(i_E, i_T) % avg_x % xs
                  end if
                end if
              else if (background == ACEFILE) then
                if (E < 1.0e6_8 * nuc % energy(1)) then
                  tope % prob_tables(i_E, i_T) % avg_x % xs = ZERO
                else
                  i_grid = binary_search(1.0e6_8 * nuc % energy, &
                    & size(nuc % energy), E)
                  if (tope % INT == LINEAR_LINEAR .or. nuc % total(i_grid) &
                    & - nuc % elastic(i_grid) &
                    & - nuc % absorption(i_grid) > XS_CUTOFF) then
                    fact = interp_factor(E, 1.0e6_8 * nuc % energy(i_grid), &
                      & 1.0e6_8 * nuc % energy(i_grid + 1), tope % INT)
                    tope % prob_tables(i_E, i_T) % avg_t % xs &
                      & = tope % prob_tables(i_E, i_T) % avg_t % xs &
                      & - tope % prob_tables(i_E, i_T) % avg_x % xs
                    tope % prob_tables(i_E, i_T) % avg_x % xs &
                      & = interpolator(fact, nuc % total(i_grid) &
                      & - nuc % elastic(i_grid) - nuc % absorption(i_grid), &
                      & nuc % total(i_grid + 1) - nuc % elastic(i_grid + 1) &
                      & - nuc % absorption(i_grid + 1), tope % INT)
                    tope % prob_tables(i_E, i_T) % avg_t % xs &
                      & = tope % prob_tables(i_E, i_T) % avg_t % xs &
                      & + tope % prob_tables(i_E, i_T) % avg_x % xs
                  else
                    tope % prob_tables(i_E, i_T) % avg_t % xs &
                      & = tope % prob_tables(i_E, i_T) % avg_t % xs &
                      & - tope % prob_tables(i_E, i_T) % avg_x % xs
                  end if
                end if
              end if
            end if

            ! set negative competitive xs to zero
            if (tope % prob_tables(i_E, i_T) % avg_x % xs < ZERO) then
              tope % prob_tables(i_E, i_T) % avg_t % xs &
                & = tope % prob_tables(i_E, i_T) % avg_t % xs &
                & + abs(tope % prob_tables(i_E, i_T) % avg_x % xs)
              tope % prob_tables(i_E, i_T) % avg_x % xs = ZERO
            end if

            if (tope % prob_tables(i_E, i_T) % avg_t % xs < xs_t_min)&
              & xs_t_min = tope % prob_tables(i_E, i_T) % avg_t % xs
            if (tope % prob_tables(i_E, i_T) % avg_t % xs > xs_t_max)&
              & xs_t_max = tope % prob_tables(i_E, i_T) % avg_t % xs

            ! accumulate the result of this history
            call tope % accum_history(i_E, i_T)
          end do

        end do HISTORY_LOOP

        ! accumulate the result of this batch
        call tope % accum_batch(i_E)

        ! calculate statistics for this batch
        call tope % calc_stats(i_b)

        ! TODO: pick the temperature used for cutoff criterion more thoughtfully
        if ((i_b > min_batches_avg_urr &
          & .and. max(tope % prob_tables(i_E, 1) % avg_t % rel_unc, &
          &           tope % prob_tables(i_E, 1) % avg_n % rel_unc, &
          &           tope % prob_tables(i_E, 1) % avg_g % rel_unc, &
          &           tope % prob_tables(i_E, 1) % avg_f % rel_unc, &
          &           tope % prob_tables(i_E, 1) % avg_x % rel_unc) &
          & < tol_avg_urr) &
          & .or. i_b == max_batches_avg_urr) exit
      end do BATCH_LOOP

      ! write probability tables out to a file
      if (write_urr_tables) then
        do i_temp = 1, tope % nT
          tope % T = tope % Tlist % get_item(i_temp)
          write(tab_unit, '(A13,ES13.6)') 'E [eV]', tope % Etabs(i_E)
          write(tab_unit, '(A13,ES13.6)') 'T [K]', tope % T
          write(tab_unit, '(A13,A13,A13,A13,A13,A13,A13,A13)') &
            & 'lower [b]', 'upper [b]', &
            & 'prob', 'total', 'elastic', 'capture', 'fission', 'competitive'

          if (master) then
            write(*, '(A32,ES10.3,A3)') 'Generated probability tables at ', &
              & tope % Etabs(i_E), ' eV'
          end if
          ptable => tope % prob_tables(i_E, i_temp)
          do i_band = 1, n_bands
            if (i_band == 1) then
              write(tab_unit, '(ES13.6,ES13.6,ES13.6,ES13.6,ES13.6,ES13.6,ES13.6,ES13.6)')&
                & xs_t_min, xs_bands(1), &
                & ptable % t(i_band) % cnt_mean, ptable % t(i_band) % xs_mean, &
                & ptable % n(i_band) % xs_mean, &
                & ptable % g(i_band) % xs_mean, &
                & ptable % f(i_band) % xs_mean, &
                & ptable % x(i_band) % xs_mean
            else if (i_band == n_bands) then
              write(tab_unit, '(ES13.6,ES13.6,ES13.6,ES13.6,ES13.6,ES13.6,ES13.6,ES13.6)')&
                & xs_bands(i_band - 1), xs_t_max, &
                & ptable % t(i_band) % cnt_mean, ptable % t(i_band) % xs_mean, &
                & ptable % n(i_band) % xs_mean, &
                & ptable % g(i_band) % xs_mean, &
                & ptable % f(i_band) % xs_mean, &
                & ptable % x(i_band) % xs_mean
            else
              write(tab_unit, '(ES13.6,ES13.6,ES13.6,ES13.6,ES13.6,ES13.6,ES13.6,ES13.6)')&
                & xs_bands(i_band - 1), xs_bands(i_band), &
                & ptable % t(i_band) % cnt_mean, ptable % t(i_band) % xs_mean, &
                & ptable % n(i_band) % xs_mean, &
                & ptable % g(i_band) % xs_mean, &
                & ptable % f(i_band) % xs_mean, &
                & ptable % x(i_band) % xs_mean
            end if
          end do
          write(tab_unit, '(A13,A13,ES13.6,ES13.6,ES13.6,ES13.6,ES13.6,ES13.6)') &
            & 'batch', 'averaged xs', &
            & ONE, ptable % avg_t % xs_mean, &
            & ptable % avg_n % xs_mean, &
            & ptable % avg_g % xs_mean, &
            & ptable % avg_f % xs_mean, &
            & ptable % avg_x % xs_mean
          write(tab_unit, '(I13,A13,ES13.6,ES13.6,ES13.6,ES13.6,ES13.6,ES13.6)') &
            & i_b, '1sigma', &
            & ZERO, ptable % avg_t % xs_sem, &
            & ptable % avg_n % xs_sem, &
            & ptable % avg_g % xs_sem, &
            & ptable % avg_f % xs_sem, &
            & ptable % avg_x % xs_sem
        end do
      end if

      ! write averaged URR cross sections out to a file
      if (write_avg_urr_xs) then
        ptable => tope % prob_tables(i_E, 1)
        write(avg_unit, '(ES13.6,ES13.6,ES13.6,ES13.6,ES13.6,ES13.6)') &
          & tope % Etabs(i_E), &
          & ptable % avg_t % xs_mean, &
          & ptable % avg_n % xs_mean, &
          & ptable % avg_g % xs_mean, &
          & ptable % avg_f % xs_mean, &
          & ptable % avg_x % xs_mean
      end if

    end do ENERGY_LOOP

    if (write_urr_tables) close(tab_unit)
    if (write_avg_urr_xs) close(avg_unit)

  end subroutine prob_tables

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! CALCULATE_URR_XS_OTF calculates URR cross section values on-the-fly,
! generating a new realization about each new E_n OR from pre-computed
! pointwise values reconstructed at simulation initialization
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

! TODO: fix-up the RRR-URR energy crossover as in calc_urr_xs_otf; probably need
! to utilize mlbw_resonances, slbw_resonances, in place of rm_resonances, where
! appropriate

  subroutine calculate_urr_xs_otf(iso, i_nuc, E, T_K)

    type(Isotope), pointer :: tope => null() ! isotope object pointer
    type(Nuclide), pointer :: nuc => null() ! nuclide object pointer
    type(Resonance) :: res ! resonance object
    type(CrossSection) :: n ! elastic scattering xs object
    type(CrossSection) :: g ! radiative capture xs object
    type(CrossSection) :: f ! fission xs object
    type(CrossSection) :: x ! competitive inelastic scattering xs object
    type(CrossSection) :: t ! total xs object
    integer :: iso      ! isotope index
    integer :: i_nuc    ! nuclide index
    integer :: i_E      ! first URR parameters energy mesh index
    integer :: iavg     ! average cross section index
    integer :: i_l      ! orbital quantum number index
    integer :: i_J      ! total angular momentum quantum number index
    integer :: i_r      ! resonance index
    integer :: n_res    ! number of resonances to include for a given l-wave
    real(8) :: avg_urr_n_xs ! infinite-dilute n xs
    real(8) :: avg_urr_f_xs ! infinite-dilute f xs
    real(8) :: avg_urr_g_xs ! infinite-dilute g xs
    real(8) :: avg_urr_x_xs ! infinite-dilute x xs
    real(8) :: inelastic_xs ! competitive inelastic scattering cross section
    real(8) :: capture_xs   ! radiative capture cross section
    real(8) :: E            ! neutron energy [eV]
    real(8) :: m            ! pointwise xs energy interpolation factor
    real(8) :: favg         ! average cross section interpolation factor
    real(8) :: T_K          ! isotope temperature [K]
!$omp threadprivate(tope, nuc)

    tope => isotopes(iso)
    nuc  => nuclides(i_nuc)

    if (tope % point_urr_xs) then
      i_E = binary_search(tope % urr_E, size(tope % urr_E), E)
      m = interp_factor(E, tope % urr_E(i_E), tope % urr_E(i_E + 1), &
        & LINEAR_LINEAR)
      micro_xs(i_nuc) % elastic = interpolator(m, &
        & tope % urr_n(i_E), tope % urr_n(i_E + 1), LINEAR_LINEAR)
      micro_xs(i_nuc) % fission = interpolator(m, &
        & tope % urr_f(i_E), tope % urr_f(i_E + 1), LINEAR_LINEAR)
      micro_xs(i_nuc) % absorption = interpolator(m, &
        & tope % urr_g(i_E) + tope % urr_f(i_E), &
        & tope % urr_g(i_E + 1) + tope % urr_f(i_E + 1), LINEAR_LINEAR)
      inelastic_xs = interpolator(m, &
        & tope % urr_x(i_E), tope % urr_x(i_E + 1), LINEAR_LINEAR)
      micro_xs(i_nuc) % total = micro_xs(i_nuc) % elastic &
                            & + micro_xs(i_nuc) % absorption &
                            & + inelastic_xs

      ! Determine nu-fission cross section
      if (tope % fissionable) then
        micro_xs(i_nuc) % nu_fission = nu_total(nuc, E/1.0e6_8) * &
          micro_xs(i_nuc) % fission
      end if
      return
    end if

    ! set current temperature and neutron energy
    tope % T = T_K
    tope % E = E

    ! reset xs accumulators
    call flush_sigmas(t, n, g, f, x)

    ! loop over orbital quantum numbers
    LOC_ORBITAL_ANG_MOM_LOOP: do i_l = 1, tope % NLS(tope % i_urr)

      ! set current orbital angular momentum quantum number
      tope % L = i_l - 1

      ! get the number of contributing l-wave resonances for this l
      n_res = n_res_contrib(tope % L)

      ! loop over total angular momentum quantum numbers
      LOC_TOTAL_ANG_MOM_LOOP: do i_J = 1, tope % NJS(i_l)

        ! set current total angular momentum quantum number
        tope % J = tope % AJ(i_l) % data(i_J)

        ! set current partial width degrees of freedom
        tope % AMUX = int(tope % DOFX(i_l) % data(i_J))
        tope % AMUN = int(tope % DOFN(i_l) % data(i_J))
        tope % AMUG = int(tope % DOFG(i_l) % data(i_J))
        tope % AMUF = int(tope % DOFF(i_l) % data(i_J))

        ! zero the resonance counter
        res % i_res = 0

        ! set mean URR parameters to neutron energy
        call set_mean_parameters(iso, E, i_l, i_J)

        ! sample unresolved resonance parameters for this spin
        ! sequence, at this energy
        call res % sample_parameters(iso, i_l, i_J)

        ! loop over the addition of resonances to this ladder
        LOC_RESONANCES_LOOP: do i_r = 1, n_res

          if (represent_params == E_RESONANCE) then
            ! interpolate mean URR parameters to current resonance energy
            call set_mean_parameters(iso, res % E_lam, i_l, i_J)
          end if

          res % i_res = i_r

          ! sample unresolved resonance parameters for this spin
          ! sequence, at this energy
          call res % sample_parameters(iso, i_l, i_J)

        end do LOC_RESONANCES_LOOP
      end do LOC_TOTAL_ANG_MOM_LOOP
    end do LOC_ORBITAL_ANG_MOM_LOOP

    ! loop over orbital quantum numbers
    ORBITAL_ANG_MOM_LOOP: do i_l = 1, tope % NLS(tope % i_urr)

      ! set current orbital angular momentum quantum number
      tope % L = i_l - 1

      ! get the number of contributing l-wave resonances for this l
      n_res = n_res_contrib(tope % L)

      ! loop over total angular momentum quantum numbers
      TOTAL_ANG_MOM_LOOP: do i_J = 1, tope % NJS(i_l)

        ! set current total angular momentum quantum number
        tope % J = tope % AJ(i_l) % data(i_J)

        ! compute statistical spin factor
        tope % g_J = (TWO * tope % J + ONE) &
             / (FOUR * tope % SPI(tope % i_urr) + TWO)

        ! loop over resonances localized about e_n
        RESONANCES_LOOP: do i_r = 1, n_res

          res % i_res = i_r
          res % E_lam = tope % local_realization(i_l, i_J) % E_lam(i_r)
          res % Gam_n = tope % local_realization(i_l, i_J) % Gam_n(i_r)
          res % Gam_g = tope % local_realization(i_l, i_J) % Gam_g(i_r)
          res % Gam_f = tope % local_realization(i_l, i_J) % Gam_f(i_r)
          res % Gam_x = tope % local_realization(i_l, i_J) % Gam_x(i_r)
          res % Gam_t = tope % local_realization(i_l, i_J) % Gam_t(i_r)

          ! calculate the contribution to the partial cross sections,
          ! at this energy, from an additional resonance
          call res % calc_xs(iso)

          ! add this contribution to the accumulated partial cross
          ! section values built up from all resonances
          ! TODO: move t outside of loop
          call accum_resonances(res, t, n, g, f, x)

        end do RESONANCES_LOOP
      end do TOTAL_ANG_MOM_LOOP
    end do ORBITAL_ANG_MOM_LOOP

    ! add potential scattering contribution
    call t % potential_xs(iso)
    call n % potential_xs(iso)

    ! interpret MF3 data according to ENDF-6 LSSF flag:
    ! MF3 contains background xs values, add to MF2 resonance contributions
    if (tope % LSSF == 0) then

      ! add resonance xs component to background
      call add_mf3_background(iso, i_nuc, NONE, t, n, g, f, x)

    ! MF3 contains evaluator-supplied background xs values that we multipy
    ! the self-shielding factors computed from MF2 by
    elseif (tope % LSSF == 1) then

      ! determine energy index
      iavg = binary_search(tope % Eavg, tope % nEavg, tope % E)

      ! tabulated unresolved resonance parameters interpolation factor
      favg = interp_factor(E, tope % Eavg(iavg), tope % Eavg(iavg + 1), &
        & tope % INT)

      ! interpolate infinite-dilute URR xs
      call interp_avg_urr_xs(favg, iso, iavg, &
        avg_urr_n_xs, avg_urr_f_xs, avg_urr_g_xs, avg_urr_x_xs)

      if (avg_urr_x_xs > ZERO) then
        if (tope % E < tope % E_ex2 .and. competitive) then
          ! self-shielded treatment of competitive inelastic cross section
          inelastic_xs = x % xs / avg_urr_x_xs &
            & * (micro_xs(i_nuc) % total &
            & -  micro_xs(i_nuc) % absorption &
            & -  micro_xs(i_nuc) % elastic)
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

      micro_xs(i_nuc) % elastic = n % xs / avg_urr_n_xs &
        & * micro_xs(i_nuc) % elastic

      ! set negative elastic xs and competitive xs to zero
      if (micro_xs(i_nuc) % elastic < ZERO) micro_xs(i_nuc) % elastic = ZERO
      if (inelastic_xs < ZERO) inelastic_xs = ZERO

      if (avg_urr_g_xs > ZERO) then
        capture_xs = g % xs / avg_urr_g_xs &
          & * (micro_xs(i_nuc) % absorption - micro_xs(i_nuc) % fission)
      else
        capture_xs = micro_xs(i_nuc) % absorption - micro_xs(i_nuc) % fission
      end if

      if (avg_urr_f_xs > ZERO) then
        micro_xs(i_nuc) % fission = f % xs / avg_urr_f_xs &
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
    if (tope % fissionable) then
      micro_xs(i_nuc) % nu_fission = nu_total(nuc, E / 1.0e6_8) &
        & * micro_xs(i_nuc) % fission
    end if

    ! set last neutron energy
    tope % E_last = tope % E

  end subroutine calculate_urr_xs_otf

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! CALCULATE_PROB_BAND_XS calculates a URR cross section from the
! OpenMC-computed probability tables
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine calculate_prob_band_xs(i_so, i_nuc, E)

    type(Isotope),  pointer, save :: tope => null()
    type(Nuclide),  pointer, save :: nuc  => null()
    integer, intent(in) :: i_so  ! index into isotopes array
    integer, intent(in) :: i_nuc ! index into nuclides array
    real(8), intent(in) :: E     ! energy
    integer :: i            ! loop index
    integer :: i_E          ! index for energy
    integer :: iT           ! index for temperature
    integer :: i_low        ! band index at lower bounding energy
    integer :: i_up         ! band index at upper bounding energy
    integer :: same_nuc_idx ! index of same nuclide
    real(8) :: f        ! interpolation factor
    real(8) :: r        ! pseudo-random number
    real(8) :: xs_n     ! elastic cross section
    real(8) :: xs_g     ! capture cross section
    real(8) :: capture  ! temporary capture cross section
    real(8) :: xs_f     ! fission cross section
    real(8) :: xs_x     ! inelastic cross section
    real(8) :: inelast  ! temporary inelastic cross section
    real(8) :: cum_prob ! cumulative band probability
    logical :: same_nuc ! do we know the xs for this nuclide at this energy?
!$omp threadprivate(tope, nuc)

    tope => isotopes(i_so)
    nuc  => nuclides(i_nuc)

    tope % T = nuc % kT / K_BOLTZMANN

    micro_xs(i_nuc) % use_ptable = .true.

    ! determine energy table
    i_E = 1
    do
      if (E < tope % Etabs(i_E + 1)) exit
      i_E = i_E + 1
    end do

    ! tabulated unresolved resonance parameters interpolation factor
    f = interp_factor(E, tope % Etabs(i_E), tope % Etabs(i_E + 1), tope % INT)

    ! if we're dealing with a nuclide that we've previously encountered at
    ! this energy but a different temperature, use the original random number
    ! to preserve correlation of temperature in probability tables
    same_nuc = .false.
    do i = 1, nuc % nuc_list % size()
      if (E /= ZERO &
        & .and. E / 1.0E6_8 &
        & == micro_xs(nuc % nuc_list % get_item(i)) % last_E) then
        same_nuc = .true.
        same_nuc_idx = i
        exit
      end if
    end do

    if (same_nuc) then
      r = micro_xs(nuc % nuc_list % get_item(same_nuc_idx)) % last_prn
    else
      r = prn()
      micro_xs(i_nuc) % last_prn = r
    end if

    iT = tope % Tlist % index(tope % T)
    i_low = 1
    cum_prob = tope % prob_tables(i_E, iT) % t(i_low) % cnt_mean
    do
      if (cum_prob > r) exit
      i_low = i_low + 1
      cum_prob = cum_prob &
        & + tope % prob_tables(i_E, iT) % t(i_low) % cnt_mean
    end do
    i_up = 1
    cum_prob = tope % prob_tables(i_E + 1, iT) % t(i_up) % cnt_mean
    do
      if (cum_prob > r) exit
      i_up = i_up + 1
      cum_prob = cum_prob &
        & + tope % prob_tables(i_E + 1, iT) % t(i_up) % cnt_mean
    end do

    ! elastic xs from probability bands
    xs_n = interpolator(f, &
      & tope % prob_tables(i_E, iT) % n(i_low) % xs_mean, &
      & tope % prob_tables(i_E + 1, iT) % n(i_up) % xs_mean, tope % INT)

    ! fission xs from probability bands
    if (tope % INT == LINEAR_LINEAR .or. &
      & tope % prob_tables(i_E, iT) % f(i_low) % xs_mean > ZERO) then
      xs_f = interpolator(f, &
        & tope % prob_tables(i_E, iT) % f(i_low) % xs_mean, &
        & tope % prob_tables(i_E + 1, iT) % f(i_up) % xs_mean, tope % INT)
    else
      xs_f = ZERO
    end if

    ! capture xs from probability bands
    if (tope % INT == LINEAR_LINEAR .or. &
      & tope % prob_tables(i_E, iT) % g(i_low) % xs_mean > ZERO) then
      xs_g = interpolator(f, &
        & tope % prob_tables(i_E, iT) % g(i_low) % xs_mean, &
        & tope % prob_tables(i_E + 1, iT) % g(i_up) % xs_mean, tope % INT)
    else
      xs_g = ZERO
    end if

    ! competitive xs from probability bands
    if (tope % INT == LINEAR_LINEAR .or. &
      & tope % prob_tables(i_E, iT) % x(i_low) % xs_mean > ZERO) then
      xs_x = interpolator(f, &
        & tope % prob_tables(i_E, iT) % x(i_low) % xs_mean, &
        & tope % prob_tables(i_E + 1, iT) % x(i_up) % xs_mean, tope % INT)
    else
      xs_x = ZERO
    end if

    if (tope % LSSF == 0) then
! TODO: move adding to ACE File 3 background to prob band calculation
      call fatal_error('LSSF=0 probability tables not yet supported for transport')
      micro_xs(i_nuc) % elastic = xs_n
      capture = xs_g
      micro_xs(i_nuc) % fission = xs_f
      inelast = xs_x
    else if (tope % LSSF == 1) then
! TODO: use stored background ACE cross sections for URR reactions
      if (xs_x > ZERO .and. competitive) then
        inelast = xs_x / interpolator(f, &
          & tope % prob_tables(i_E, iT) % avg_x % xs_mean, &
          & tope % prob_tables(i_E + 1, iT) % avg_x % xs_mean, tope % INT)&
          & * (micro_xs(i_nuc) % total - micro_xs(i_nuc) % absorption &
          & - micro_xs(i_nuc) % elastic)
      else
        inelast = xs_x
      end if
! TODO: move this division to prob table calculation
      micro_xs(i_nuc) % elastic = xs_n / interpolator(f, &
        & tope % prob_tables(i_E, iT) % avg_n % xs_mean, &
        & tope % prob_tables(i_E + 1, iT) % avg_n % xs_mean, tope % INT)&
        & * micro_xs(i_nuc) % elastic
      if (xs_g > ZERO) then
        capture = xs_g / interpolator(f, &
          & tope % prob_tables(i_E, iT) % avg_g % xs_mean, &
          & tope % prob_tables(i_E + 1, iT) % avg_g % xs_mean, tope % INT)&
          & * (micro_xs(i_nuc) % absorption - micro_xs(i_nuc) % fission)
      else
        capture = micro_xs(i_nuc) % absorption - micro_xs(i_nuc) % fission
      end if
      if (xs_f > ZERO) then
        micro_xs(i_nuc) % fission = xs_f / interpolator(f, &
          & tope % prob_tables(i_E, iT) % avg_f % xs_mean, &
          & tope % prob_tables(i_E + 1, iT) % avg_f % xs_mean, tope % INT)&
          & * micro_xs(i_nuc) % fission
      else
        micro_xs(i_nuc) % fission = micro_xs(i_nuc) % fission
      end if
    end if

    micro_xs(i_nuc) % absorption = capture + micro_xs(i_nuc) % fission
    micro_xs(i_nuc) % total = micro_xs(i_nuc) % elastic &
      & + micro_xs(i_nuc) % absorption + inelast

    ! Determine nu-fission cross section
    if (nuc % fissionable) then
      micro_xs(i_nuc) % nu_fission = nu_total(nuc, E / 1.0E6_8) * &
           micro_xs(i_nuc) % fission
    end if

    tope % E_last = E

  end subroutine calculate_prob_band_xs

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! SAMPLE_PARAMETERS samples unresolved resonance parameters for the next
! pseudo-resonance added to the ladder
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine sample_parameters(this, iso, i_l, i_J)

    class(Resonance), intent(inout) :: this ! pseudo-resonance object
    integer :: iso ! isotope index
    integer :: i_l ! orbital angular momentum quantum number index
    integer :: i_J ! total angular momentum quantum number index

    ! sample unresolved resonance parameters for this resonance
    call this % level_spacing(iso, i_l, i_J)
    call this % channel_width(iso, i_l, i_J)

  end subroutine sample_parameters

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! LEVEL_SPACING samples the energy spacing between adjacent resonances
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine level_spacing(this, iso, i_l, i_J)

    class(Resonance), intent(inout) :: this ! pseudo-resonance object
    type(Isotope), pointer :: tope => null() ! nuclide pointer
    integer :: iso   ! isotope index
    integer :: n_res ! number of resonances to include for a given l-wave
    integer :: i_l   ! orbital angular momentum quantum number index
    integer :: i_J   ! total angular momentum quantum number index

    tope => isotopes(iso)

    if (tope % E == tope % E_last) then
      ! Energy hasn't changed since last realization, so use the same one
      this % E_lam &
           = tope % local_realization(i_l, i_J) % E_lam(this % i_res)
    else
      ! sample a level spacing from the Wigner distribution
      this % D_lJ = wigner_dist(tope % D)
      if (this % i_res == 0) then
        ! set lowest energy (i.e. the first) resonance for this ladder well below
        ! the energy grid point such that the ladder spans a sufficient energy
        ! range
        n_res = n_res_contrib(tope % L)
        this % E_lam = (tope % E - n_res/2 * tope % D) &
             + (ONE - TWO * prn()) * this % D_lJ
      else
        ! add subsequent resonance energies at the sampled spacing above the
        ! last resonance
        this % E_lam = this % E_lam + this % D_lJ
        tope % local_realization(i_l, i_J) % E_lam(this % i_res) = this % E_lam
      end if
    end if

  end subroutine level_spacing

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! N_RES_CONTRIB determines the number of resonances to include from the lth
! wave that contribute to a URR cross section value
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  function n_res_contrib(L) result(n_res)

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
      call fatal_error('Only s, p, and d wave resonances are supported &
        & in ENDF-6')
      n_res = l_waves(4)
    case default
      call fatal_error('Only s, p, and d wave resonances are supported &
        & in ENDF-6')
    end select

  end function n_res_contrib

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
! CHANNEL_WIDTH samples the channel partial widths at E_lambda when generating
! a full resonance ensemble or at E_n when generating localized parameters for
! an on-the-fly cross section calculation
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine channel_width(this, iso, i_l, i_J)

    class(Resonance), intent(inout) :: this ! pseudo-resonance object
    type(Isotope), pointer :: tope => null() ! isotope object pointer
    integer :: iso    ! isotope index
    integer :: i_l    ! orbital angular momentum quantum number index
    integer :: i_J    ! total angular momentum quantum number index
    integer :: i_tabn ! elastic chi-squared table index
    integer :: i_tabg ! capture chi-squared table index
    integer :: i_tabf ! fission chi-squared table index
    integer :: i_tabx ! competitivechi-squared table index
    real(8) :: rho    ! derived variable
    real(8) :: nu     ! derived variable

    tope => isotopes(iso)

    ! TODO: Actually sample a chi-squared distribution rather than using
    ! tabulated values. Look at the third Monte Carlo Sampler?

    if (tope % E == tope % E_last) then
      ! Energy hasn't changed since last realization, so use the same one
      this % Gam_n = tope % local_realization(i_l, i_J) % Gam_n(this % i_res)
      this % Gam_f = tope % local_realization(i_l, i_J) % Gam_f(this % i_res)
      this % Gam_g = tope % local_realization(i_l, i_J) % Gam_g(this % i_res)
      this % Gam_x = tope % local_realization(i_l, i_J) % Gam_x(this % i_res)
      this % Gam_t = tope % local_realization(i_l, i_J) % Gam_t(this % i_res)

    else
      ! sample indices into table of equiprobable chi-squared function values
      i_tabn = ceiling(prn() * 20.0_8)
      i_tabf = ceiling(prn() * 20.0_8)
      i_tabg = ceiling(prn() * 20.0_8)
      i_tabx = ceiling(prn() * 20.0_8)

      ! use the sampled tabulated chi-squared values to calculate sample widths
      ! neutron width
      if (tope % AMUN > 0) then
        ! compute factors needed to go from the mean reduced width that is
        ! provided by ENDF for elastic scattering to what we want - a partial width
        ! (use absolute value of energy in order to handle bound levels which have
        ! negative resonance energies - this is an ENDF-6 convention, not theory)
        if (represent_params == E_NEUTRON) then
          rho = wavenumber(tope % AWR, abs(tope % E)) * tope % ac(tope % i_urr)
          nu  = penetration(tope % L, rho) / rho
          this % Gam_n = tope % GN0 * sqrt(abs(tope % E)) * nu &
            & * chi2(i_tabn, tope % AMUN)
        else if (represent_params == E_RESONANCE) then
          rho = wavenumber(tope % AWR, abs(this % E_lam)) * tope % ac(tope % i_urr)
          nu  = penetration(tope % L, rho) / rho
          this % Gam_n = tope % GN0 * sqrt(abs(this % E_lam)) * nu &
            & * chi2(i_tabn, tope % AMUN)          
        end if
      else
        call fatal_error('Non-positive neutron width sampled')
      end if

      ! fission width
      if (tope % AMUF > 0) then
        this % Gam_f = tope % GF * chi2(i_tabf, tope % AMUF) / dble(tope % AMUF)
      else
        this % Gam_f = ZERO
      end if

      ! constant radiative width
      ! (many channels --> many degrees of freedom --> Dirac delta)
      this % Gam_g = tope % GG

      ! competitive width
      if (tope % AMUX > 0) then
        this % Gam_x = tope % GX * chi2(i_tabx, tope % AMUX) / dble(tope % AMUX)
      else
        this % Gam_x = ZERO
      end if

      ! total width (sum of partials)
      this % Gam_t = this % Gam_n + this % Gam_f + this % Gam_g + this % Gam_x

      if (this % i_res /= 0) then
        tope % local_realization(i_l, i_J) % Gam_n(this % i_res) = this % Gam_n
        tope % local_realization(i_l, i_J) % Gam_f(this % i_res) = this % Gam_f
        tope % local_realization(i_l, i_J) % Gam_g(this % i_res) = this % Gam_g
        tope % local_realization(i_l, i_J) % Gam_x(this % i_res) = this % Gam_x
        tope % local_realization(i_l, i_J) % Gam_t(this % i_res) = this % Gam_t
      end if
    end if

  end subroutine channel_width

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! CALC_XS is an interface for the calculation of partial cross sections at E_n,
! the energy that the ladder is being generated about
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine calc_xs(this, iso)

    class(Resonance), intent(inout) :: this ! pseudo-resonance object
    integer :: iso ! isotope index

    ! calculate cross section contributions from an additional resonance
    select case(formalism)
    case (SLBW)
      call this % slbw_xs(iso)
    case (MLBW)
      call this % mlbw_xs(iso)
    case (REICH_MOORE)
      call fatal_error('Reich-Moore formalism not yet supported for the URR')
    case (MNBW)
      call fatal_error('MNBW formalism not yet supported for the URR')
    case default
      call fatal_error('Unrecognized URR resonance formalism')
    end select

  end subroutine calc_xs

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! SLBW_XS calculates single-level Breit-Wigner cross sections at E_n
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine slbw_xs(this, iso)

    class(Resonance), intent(inout) :: this ! pseudo-resonance object
    type(Isotope), pointer :: tope => null() ! nuclide pointer
    integer :: iso ! isotope index
    real(8) :: k_n     ! center-of-mass neutron wavenumber at E_n
    real(8) :: k_n_x   ! center-of-mass neutron wavenumber at E_n - QI
    real(8) :: k_lam   ! center-of-mass neutron wavenumber at E_lam
    real(8) :: k_lam_x ! center-of-mass neutron wavenumber at |E_lam - QI|
    real(8) :: E_shift ! shifted resonance energy in the lab system
    real(8) :: theta   ! total width / Doppler width
    real(8) :: x       ! derived variable
    real(8) :: Gam_t_n ! sampled energy-dependent total width at E_n
    real(8) :: Gam_n_n ! sampled energy-dependent neutron width at E_n
    real(8) :: Gam_x_n ! sampled energy-dependent competitive width at E_n
    real(8) :: sig_lam ! peak resonance cross section
    real(8) :: sig_lam_Gam_t_n_psi ! compound variable

    tope => isotopes(iso)

    k_n = wavenumber(tope % AWR, tope % E)

    ! (use absolute value of energy in order to handle bound levels which have
    ! negative resonance energies)
    k_lam = wavenumber(tope % AWR, abs(this % E_lam))

    ! if URR parameters have resonance energy dependence
    if (represent_params == E_RESONANCE) then

      E_shift = this % E_lam &
        & + this % Gam_n * (shift(tope % L, k_lam * tope % ac(tope % i_urr)) &
        & - shift(tope % L, k_n * tope % ac(tope % i_urr))) &
        & / (TWO * penetration(tope % L, k_lam * tope % ac(tope % i_urr)))

      Gam_n_n = this % Gam_n &
        & * penetration(tope % L, k_n   * tope % ac(tope % i_urr)) &
        & / penetration(tope % L, k_lam * tope % ac(tope % i_urr))

      if (tope % E >= tope % E_ex2) then
        ! two competitive reactions possible, can't calculate an energy-dependent
        ! width because it depends on the two (unprovided) reaction partial widths
        Gam_x_n = this % Gam_x
      else if (tope % E >= tope % E_ex1) then
        ! can compute an energy-dependent width for the one competitive reaction
        k_n_x = wavenumber(tope % AWR, tope % E - tope % E_ex1)
        k_lam_x = wavenumber(tope % AWR, abs(this % E_lam - tope % E_ex1))
        Gam_x_n = this % Gam_x &
          & * penetration(tope % L, k_n_x   * tope % ac(tope % i_urr)) &
          & / penetration(tope % L, k_lam_x * tope % ac(tope % i_urr))
      else
        Gam_x_n = ZERO
      end if

    else
      ! assume all URR parameters already have neutron energy dependence

      E_shift = this % E_lam
      Gam_n_n = this % Gam_n
      Gam_x_n = this % Gam_x

    end if

    Gam_t_n = this % Gam_t - this % Gam_n - this % Gam_x &
      & + Gam_n_n + Gam_x_n

    theta = HALF * Gam_t_n &
      & / sqrt(K_BOLTZMANN * 1.0E6_8 * tope % T * tope % E / tope % AWR)

    x = (TWO * (tope % E - E_shift)) / Gam_t_n

    sig_lam = FOUR * PI / (k_lam * k_lam) * tope % g_J &
      & * this % Gam_n / this % Gam_t

    ! this particular form comes from the NJOY2012 manual
    this % dxs_n = sig_lam * &
      & ((cos(TWO * phase_shift(tope % L, k_n * tope % AP(tope % i_urr))) &
      & - (ONE - Gam_n_n / Gam_t_n)) * psi(tope%T, theta, x) &
      & + sin(TWO * phase_shift(tope % L, k_n * tope % AP(tope % i_urr))) &
      & * chi(tope%T, theta, x))

    sig_lam_Gam_t_n_psi = sig_lam * psi(tope%T, theta, x) / Gam_t_n

    if (this % Gam_g > ZERO) then
      this % dxs_g = sig_lam_Gam_t_n_psi * this % Gam_g
    else
      this % dxs_g = ZERO
    end if

    if (this % Gam_f > ZERO) then
      this % dxs_f = sig_lam_Gam_t_n_psi * this % Gam_f
    else
      this % dxs_f = ZERO
    end if

    ! can only have a competitive resonance component w/ a single open channel
    if (Gam_x_n > ZERO .and. tope % E < tope % E_ex2) then
      this % dxs_x = sig_lam_Gam_t_n_psi * Gam_x_n
    else
      this % dxs_x = ZERO
    end if

    this % dxs_t = this % dxs_n &
               & + this % dxs_g &
               & + this % dxs_f &
               & + this % dxs_x

  end subroutine slbw_xs

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! MLBW_XS calculates a multi-level Breit-Wigner elastic scattering cross
! section and single-level Breit-Wigner cross sections for other reactions at
! E_n
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine mlbw_xs(this, iso)

    class(Resonance), intent(inout) :: this ! pseudo-resonance object
    type(Isotope), pointer :: tope => null() ! nuclide pointer
    integer :: iso ! isotope index
    real(8) :: k_n     ! center-of-mass neutron wavenumber at E_n
    real(8) :: k_n_x   ! center-of-mass neutron wavenumber at E_n - QI
    real(8) :: k_lam   ! center-of-mass neutron wavenumber at E_lam
    real(8) :: k_lam_x ! center-of-mass neutron wavenumber at |E_lam - QI|
    real(8) :: E_shift ! shifted resonance energy in the lab system
    real(8) :: theta   ! total width / Doppler width
    real(8) :: x       ! derived variable
    real(8) :: Gam_t_n ! sampled energy-dependent total width at E_n
    real(8) :: Gam_n_n ! sampled energy-dependent neutron width at E_n
    real(8) :: Gam_x_n ! sampled energy-dependent competitive width at E_n
    real(8) :: sig_lam ! peak resonance cross section
    real(8) :: sig_lam_Gam_t_n_psi ! compound variable

    tope => isotopes(iso)

    if (tope % J == tope % SPI(tope % i_urr)) then
      call fatal_error('Computing MLBW elastic scattering cross section&
           & for resonance with total orbital angular momentum quantum&
           & number, J, equal to the nuclear spin, I')
    end if

    k_n = wavenumber(tope % AWR, tope % E)

    ! (use absolute value of energy in order to handle bound levels which have
    ! negative resonance energies)
    k_lam = wavenumber(tope % AWR, abs(this % E_lam))

    ! if URR parameters have resonance energy dependence
    if (represent_params == E_RESONANCE) then

      E_shift = this % E_lam &
        & + this % Gam_n * (shift(tope % L, k_lam * tope % ac(tope % i_urr)) &
        & - shift(tope % L, k_n * tope % ac(tope % i_urr))) &
        & / (TWO * penetration(tope % L, k_lam * tope % ac(tope % i_urr)))

      Gam_n_n = this % Gam_n &
        & * penetration(tope % L, k_n   * tope % ac(tope % i_urr)) &
        & / penetration(tope % L, k_lam * tope % ac(tope % i_urr))

      if (tope % E >= tope % E_ex2) then
        ! two competitive reactions possible, can't calculate an energy-dependent
        ! width because it depends on the two (unprovided) reaction partial widths
        Gam_x_n = this % Gam_x
      else if (tope % E >= tope % E_ex1) then
        ! can compute an energy-dependent width for the one competitive reaction
        k_n_x = wavenumber(tope % AWR, tope % E - tope % E_ex1)
        k_lam_x = wavenumber(tope % AWR, abs(this % E_lam - tope % E_ex1))
        Gam_x_n = this % Gam_x &
          & * penetration(tope % L, k_n_x   * tope % ac(tope % i_urr)) &
          & / penetration(tope % L, k_lam_x * tope % ac(tope % i_urr))
      else
        Gam_x_n = ZERO
      end if

    else

      ! assume all URR parameters already have neutron energy dependence
      E_shift = this % E_lam
      Gam_n_n = this % Gam_n
      Gam_x_n = this % Gam_x

    end if

    Gam_t_n = this % Gam_t - this % Gam_n - this % Gam_x &
         + Gam_n_n + Gam_x_n

    theta = HALF * Gam_t_n &
         / sqrt(K_BOLTZMANN * 1.0E6_8 * tope % T * tope % E / tope % AWR)

    x = (TWO * (tope % E - E_shift)) / Gam_t_n

    sig_lam = FOUR * PI / (k_lam * k_lam) * tope % g_J &
         * this % Gam_n / this % Gam_t

    ! this particular form comes from the NJOY2012 manual
    this % dxs_n = sig_lam * &
         ((cos(TWO * phase_shift(tope % L, k_n * tope % AP(tope % i_urr))) &
         - (ONE - Gam_n_n / Gam_t_n)&
         + HALF * G_func(iso, E_shift, Gam_n_n, Gam_t_n, this % i_res)&
         / Gam_n_n) * psi(tope % T, theta, x) &
         + (sin(TWO * phase_shift(tope % L, k_n * tope % AP(tope % i_urr))) &
         + H_func(iso, E_shift, Gam_n_n, Gam_t_n, this % i_res)&
         / Gam_n_n) * chi(tope % T, theta, x))

    sig_lam_Gam_t_n_psi = sig_lam * psi(tope % T, theta, x) / Gam_t_n

    if (this % Gam_g > ZERO) then
      this % dxs_g = sig_lam_Gam_t_n_psi * this % Gam_g
    else
      this % dxs_g = ZERO
    end if

    if (this % Gam_f > ZERO) then
      this % dxs_f = sig_lam_Gam_t_n_psi * this % Gam_f
    else
      this % dxs_f = ZERO
    end if

    ! can only have a competitive resonance component w/ a single open channel
    if (Gam_x_n > ZERO .and. tope % E < tope % E_ex2) then
      this % dxs_x = sig_lam_Gam_t_n_psi * Gam_x_n
    else
      this % dxs_x = ZERO
    end if

    this % dxs_t = this % dxs_n &
         + this % dxs_g &
         + this % dxs_f &
         + this % dxs_x

  end subroutine mlbw_xs

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! G_FUNC calculates a value for the G-function appearing in the NJOY-2012 form
! of the MLBW resonance formalae
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  function G_func(iso, E_res, G_n, G_t, i_res) result(G_val)

    type(Isotope), pointer :: tope => null() ! isotope pointer
    integer :: iso   ! isotope index
    integer :: i_l   ! orbital angular momentum quantum number index
    integer :: i_J   ! total angular momentum quantum number index
    integer :: i_r   ! spin sequence resonance index
    integer :: i_res ! index of current nearest resonance
    real(8) :: G_val ! G-function value
    real(8) :: E_res ! energy of the resonance contributing to the xs at E_n
    real(8) :: G_n   ! neutron width of the resonance at E_res
    real(8) :: G_t   ! total width of the resonance at E_res

    tope => isotopes(iso)

    G_val = ZERO

    i_l = tope % L + 1

    TOTAL_ANG_MOM_LOOP: do i_J = 1, tope % NJS(i_l)
      RESONANCE_LOOP: do i_r = 1, n_res_contrib(tope % L)
! TODO: compute E-dependence of local realization's parameters
        if (i_r == i_res .and. tope % J == tope % AJ(i_l) % data(i_J)) cycle

        G_val = G_val&
             + G_n * tope % local_realization(i_l, i_J) % Gam_n(i_r)&
             * (G_t + tope % local_realization(i_l, i_J) % Gam_t(i_r))&
             / ((E_res - tope % local_realization(i_l, i_J) % E_lam(i_r))&
             * (E_res - tope % local_realization(i_l, i_J) % E_lam(i_r))&
             + (G_t + tope % local_realization(i_l, i_J) % Gam_t(i_r))&
             * (G_t + tope % local_realization(i_l, i_J) % Gam_t(i_r)) / 4.0_8)

      end do RESONANCE_LOOP
    end do TOTAL_ANG_MOM_LOOP

  end function G_func


!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! H_FUNC calculates a value for the H-function appearing in the NJOY-2012 form
! of the MLBW resonance formalae
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  function H_func(iso, E_res, G_n, G_t, i_res) result(H_val)

    type(Isotope), pointer :: tope => null() ! isotope pointer
    integer :: iso   ! isotope index
    integer :: i_l   ! orbital angular momentum quantum number index
    integer :: i_J   ! total angular momentum quantum number index
    integer :: i_r   ! spin sequence resonance index
    integer :: i_res ! index of current nearest resonance
    real(8) :: H_val ! G-function value
    real(8) :: E_res ! energy of the resonance contributing to the xs at E_n
    real(8) :: G_n   ! neutron width of the resonance at E_res
    real(8) :: G_t   ! total width of the resonance at E_res

    tope => isotopes(iso)

    H_val = ZERO

    i_l = tope % L + 1

    TOTAL_ANG_MOM_LOOP: do i_J = 1, tope % NJS(i_l)
      RESONANCE_LOOP: do i_r = 1, n_res_contrib(tope % L)
! TODO: compute E-dependence of local realization's parameters
        if (i_r == i_res .and. tope % J == tope % AJ(i_l) % data(i_J)) cycle

        H_val = H_val&
             + G_n * tope % local_realization(i_l, i_J) % Gam_n(i_r)&
             * (E_res - tope % local_realization(i_l, i_J) % E_lam(i_r))&
             / ((E_res - tope % local_realization(i_l, i_J) % E_lam(i_r))&
             * (E_res - tope % local_realization(i_l, i_J) % E_lam(i_r))&
             + (G_t + tope % local_realization(i_l, i_J) % Gam_t(i_r))&
             * (G_t + tope % local_realization(i_l, i_J) % Gam_t(i_r)) / 4.0_8)

      end do RESONANCE_LOOP
    end do TOTAL_ANG_MOM_LOOP

  end function H_func

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! PENETRATION calculates hard sphere penetrability factors
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  function penetration(L, rho) result(P)

    integer, intent(in) :: L    ! current orbital quantum number
    real(8) :: rho  ! derived variable, ka
    real(8) :: rho2 ! rho**2
    real(8) :: P    ! penetration factor

    ! pre-compute exponentiation
    rho2 = rho * rho

    ! calculate penetrability for the current orbital quantum number
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

    integer :: L    ! current orbital quantum number
    real(8) :: rho  ! derived variable, ka
    real(8) :: rho2 ! rho**2
    real(8) :: rho3 ! rho**3
    real(8) :: rho4 ! rho**4
    real(8) :: phi  ! hard sphere phase shift

    ! pre-compute exponentiations
    rho2 = rho * rho
    rho3 = rho * rho2
    rho4 = rho * rho3

    ! calculate phase shift for the current orbital quantum number
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

    integer :: L    ! current orbital quantum number
    real(8) :: rho  ! derived variable, ka
    real(8) :: rho2 ! rho**2
    real(8) :: S    ! shift factor (for shifting the resonance energy)

    ! pre-compute exponentiation
    rho2 = rho * rho

    ! calculate shift factor for current orbital quantum number
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
! PSI computes a value of the psi Doppler integral function
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  function psi(T, theta, x) result(psi_val)

    real(8)    :: T       ! temperature [K]
    real(8)    :: theta   ! psi argument
    real(8)    :: x       ! psi argument
    real(8)    :: psi_val ! calculated value of psi
    real(8)    :: relerr  ! relative error of the Faddeeva evaluation
    complex(8) :: w_val   ! complex return value of the Faddeeva evaluation

    if (T > ZERO) then

      ! evaluate the W (Faddeeva) function
      select case (w_eval)

        ! call S.G. Johnson's Faddeeva evaluation
      case (MIT_W)
        relerr = 1.0e-6
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

    else

      psi_val = ONE / (ONE + x*x)

    end if

  end function psi

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! CHI computes a value of the chi Doppler integral function
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  function chi(T, theta, x) result(chi_val)

    real(8)    :: T       ! temperature [K]
    real(8)    :: theta   !
    real(8)    :: x       !
    real(8)    :: chi_val ! calculated value of chi
    complex(8) :: w_val   ! complex return value of the Faddeeva evaluation
    real(8)    :: relerr  ! relative error of the Faddeeva evaluation

    if (T > ZERO) then

      ! evaluate the W (Faddeeva) function
      select case (w_eval)

        ! S.G. Johnson's Faddeeva evaluation
      case (MIT_W)

        relerr = 1.0e-6
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

    else

      chi_val = x / (ONE + x*x)

    end if

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

  subroutine accum_resonance(this, dxs)

    class(CrossSection), intent(inout) :: this ! cross section object
    real(8) :: dxs ! contribution to xs from the new resonance

    ! add xs contribution from a new resonance to the xs value at E_n
    this % xs = this % xs + dxs

  end subroutine accum_resonance

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! POTENTIAL_XS adds the contribution of potential scattering to the elastic
! and total cross sections
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine potential_xs(this, iso)

    class(CrossSection), intent(inout) :: this ! cross section object
    type(Isotope), pointer :: tope => null() ! nuclide pointer
    integer :: iso ! isotope index
    integer :: i_l ! orbital quantum number index
    real(8) :: sig_pot ! potential scattering cross section
    real(8) :: k_n     ! center-of-mass neutron wavenumber

    ! set nuclide variables
    tope => isotopes(iso)

    ! compute neutron COM wavenumber
    k_n = wavenumber(tope % AWR, tope % E)

    ! compute potential scattering xs by adding contribution from each l-wave
    sig_pot = ZERO
    do i_l = 0, tope % NLS(tope % i_urr) - 1
      sig_pot = sig_pot &
        & + FOUR * PI / (k_n * k_n) * (TWO * dble(i_l) + ONE) &
        & * (sin(phase_shift(i_l, k_n * tope % AP(tope % i_urr)))) &
        & * (sin(phase_shift(i_l, k_n * tope % AP(tope % i_urr))))
    end do

    ! add the potential scattering xs to this xs
    this % xs = this % xs + sig_pot

  end subroutine potential_xs

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! ACCUM_HISTORY adds the single-history xs realization to the single-batch
! accumulator
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine accum_history(this, i_E, i_T)

    class(Isotope), intent(inout) :: this ! cross section object
    integer :: i_E ! tabulated URR energy index
    integer :: i_T ! temperature index
    integer :: i_b ! total xs band that this realization falls in

    ! accumulate infinite-dilute xs values for this realization
    this % prob_tables(i_E, i_T) % avg_t % xs_tmp_sum &
      & = this % prob_tables(i_E, i_T) % avg_t % xs_tmp_sum &
      & + this % prob_tables(i_E, i_T) % avg_t % xs
    this % prob_tables(i_E, i_T) % avg_n % xs_tmp_sum &
      & = this % prob_tables(i_E, i_T) % avg_n % xs_tmp_sum &
      & + this % prob_tables(i_E, i_T) % avg_n % xs
    this % prob_tables(i_E, i_T) % avg_g % xs_tmp_sum &
      & = this % prob_tables(i_E, i_T) % avg_g % xs_tmp_sum &
      & + this % prob_tables(i_E, i_T) % avg_g % xs
    this % prob_tables(i_E, i_T) % avg_f % xs_tmp_sum &
      & = this % prob_tables(i_E, i_T) % avg_f % xs_tmp_sum &
      & + this % prob_tables(i_E, i_T) % avg_f % xs
    this % prob_tables(i_E, i_T) % avg_x % xs_tmp_sum &
      & = this % prob_tables(i_E, i_T) % avg_x % xs_tmp_sum &
      & + this % prob_tables(i_E, i_T) % avg_x % xs

    ! determine band number
    i_b = 1
    do
      if (this % prob_tables(i_E, i_T) % avg_t % xs < xs_bands(i_b)) exit
      i_b = i_b + 1
      if (i_b == n_bands) exit
    end do

    ! add partial xs values to the band's single-batch sums
    this % prob_tables(i_E, i_T) % t(i_b) % cnt_tmp_sum &
      & = this % prob_tables(i_E, i_T) % t(i_b) % cnt_tmp_sum + 1
    this % prob_tables(i_E, i_T) % t(i_b) % xs_tmp_sum &
      & = this % prob_tables(i_E, i_T) % t(i_b) % xs_tmp_sum &
      & + this % prob_tables(i_E, i_T) % avg_t % xs
    this % prob_tables(i_E, i_T) % n(i_b) % xs_tmp_sum &
      & = this % prob_tables(i_E, i_T) % n(i_b) % xs_tmp_sum &
      & + this % prob_tables(i_E, i_T) % avg_n % xs
    this % prob_tables(i_E, i_T) % g(i_b) % xs_tmp_sum &
      & = this % prob_tables(i_E, i_T) % g(i_b) % xs_tmp_sum &
      & + this % prob_tables(i_E, i_T) % avg_g % xs
    this % prob_tables(i_E, i_T) % f(i_b) % xs_tmp_sum &
      & = this % prob_tables(i_E, i_T) % f(i_b) % xs_tmp_sum &
      & + this % prob_tables(i_E, i_T) % avg_f % xs
    this % prob_tables(i_E, i_T) % x(i_b) % xs_tmp_sum &
      & = this % prob_tables(i_E, i_T) % x(i_b) % xs_tmp_sum &
      & + this % prob_tables(i_E, i_T) % avg_x % xs

  end subroutine accum_history

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! FLUSH_HISTORIES flushes the single-history probability table cross section
! values
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine flush_histories(this)

    class(Isotope), intent(inout) :: this ! isotope object
    integer :: i_E ! tabulated URR energy index
    integer :: i_T ! temperature index

    do i_E = 1, this % ntabs
      do i_T = 1, this % nT
        this % prob_tables(i_E, i_T) % avg_t % xs = ZERO
        this % prob_tables(i_E, i_T) % avg_n % xs = ZERO
        this % prob_tables(i_E, i_T) % avg_g % xs = ZERO
        this % prob_tables(i_E, i_T) % avg_f % xs = ZERO
        this % prob_tables(i_E, i_T) % avg_x % xs = ZERO
      end do
    end do

  end subroutine flush_histories

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! ACCUM_BATCH adds the single-batch results to the overall accumulators
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine accum_batch(this, i_E)

    class(Isotope), intent(inout) :: this ! cross section object
    integer :: i_E ! tabulated URR energy index
    integer :: i_T ! temperature index
    integer :: i_b ! probability band index

    if (this % E /= this % Etabs(i_E)) return
    do i_T = 1, this % nT
      do i_b = 1, n_bands
        ! accumulate single-batch band means and squared band means
        if (this % prob_tables(i_E, i_T) % t(i_b) % cnt_tmp_sum /= 0) then
          this % prob_tables(i_E, i_T) % t(i_b) % xs_sum &
            & = this % prob_tables(i_E, i_T) % t(i_b) % xs_sum &
            & + this % prob_tables(i_E, i_T) % t(i_b) % xs_tmp_sum
          this % prob_tables(i_E, i_T) % t(i_b) % xs_sum2 &
            & = this % prob_tables(i_E, i_T) % t(i_b) % xs_sum2 &
            & + (this % prob_tables(i_E, i_T) % t(i_b) % xs_tmp_sum &
            & *  this % prob_tables(i_E, i_T) % t(i_b) % xs_tmp_sum)
          this % prob_tables(i_E, i_T) % n(i_b) % xs_sum &
            & = this % prob_tables(i_E, i_T) % n(i_b) % xs_sum &
            & + this % prob_tables(i_E, i_T) % n(i_b) % xs_tmp_sum
          this % prob_tables(i_E, i_T) % n(i_b) % xs_sum2 &
            & = this % prob_tables(i_E, i_T) % n(i_b) % xs_sum2 &
            & + (this % prob_tables(i_E, i_T) % n(i_b) % xs_tmp_sum &
            & *  this % prob_tables(i_E, i_T) % n(i_b) % xs_tmp_sum)
          this % prob_tables(i_E, i_T) % g(i_b) % xs_sum &
            & = this % prob_tables(i_E, i_T) % g(i_b) % xs_sum &
            & + this % prob_tables(i_E, i_T) % g(i_b) % xs_tmp_sum
          this % prob_tables(i_E, i_T) % g(i_b) % xs_sum2 &
            & = this % prob_tables(i_E, i_T) % g(i_b) % xs_sum2 &
            & + (this % prob_tables(i_E, i_T) % g(i_b) % xs_tmp_sum &
            & *  this % prob_tables(i_E, i_T) % g(i_b) % xs_tmp_sum)
          this % prob_tables(i_E, i_T) % f(i_b) % xs_sum &
            & = this % prob_tables(i_E, i_T) % f(i_b) % xs_sum &
            & + this % prob_tables(i_E, i_T) % f(i_b) % xs_tmp_sum
          this % prob_tables(i_E, i_T) % f(i_b) % xs_sum2 &
            & = this % prob_tables(i_E, i_T) % f(i_b) % xs_sum2 &
            & + (this % prob_tables(i_E, i_T) % f(i_b) % xs_tmp_sum &
            & *  this % prob_tables(i_E, i_T) % f(i_b) % xs_tmp_sum)
          this % prob_tables(i_E, i_T) % x(i_b) % xs_sum &
            & = this % prob_tables(i_E, i_T) % x(i_b) % xs_sum &
            & + this % prob_tables(i_E, i_T) % x(i_b) % xs_tmp_sum
          this % prob_tables(i_E, i_T) % x(i_b) % xs_sum2 &
            & = this % prob_tables(i_E, i_T) % x(i_b) % xs_sum2 &
            & + (this % prob_tables(i_E, i_T) % x(i_b) % xs_tmp_sum &
            & *  this % prob_tables(i_E, i_T) % x(i_b) % xs_tmp_sum)

          ! accumulate single-batch band counts
          this % prob_tables(i_E, i_T) % t(i_b) % cnt &
            & = this % prob_tables(i_E, i_T) % t(i_b) % cnt &
            & + this % prob_tables(i_E, i_T) % t(i_b) % cnt_tmp_sum

          ! accumulate squared single-batch band counts
          this % prob_tables(i_E, i_T) % t(i_b) % cnt2 &
            & = this % prob_tables(i_E, i_T) % t(i_b) % cnt2 &
            & + (this % prob_tables(i_E, i_T) % t(i_b) % cnt_tmp_sum &
            & *  this % prob_tables(i_E, i_T) % t(i_b) % cnt_tmp_sum)
        end if
      end do

      ! accumulate single-batch infinite-dilute means
      this % prob_tables(i_E, i_T) % avg_t % xs_sum &
        & = this % prob_tables(i_E, i_T) % avg_t % xs_sum &
        & + this % prob_tables(i_E, i_T) % avg_t % xs_tmp_sum &
        & / dble(histories_avg_urr)
      this % prob_tables(i_E, i_T) % avg_n % xs_sum &
        & = this % prob_tables(i_E, i_T) % avg_n % xs_sum &
        & + this % prob_tables(i_E, i_T) % avg_n % xs_tmp_sum &
        & / dble(histories_avg_urr)
      this % prob_tables(i_E, i_T) % avg_g % xs_sum &
        & = this % prob_tables(i_E, i_T) % avg_g % xs_sum &
        & + this % prob_tables(i_E, i_T) % avg_g % xs_tmp_sum &
        & / dble(histories_avg_urr)
      this % prob_tables(i_E, i_T) % avg_f % xs_sum &
        & = this % prob_tables(i_E, i_T) % avg_f % xs_sum &
        & + this % prob_tables(i_E, i_T) % avg_f % xs_tmp_sum &
        & / dble(histories_avg_urr)
      this % prob_tables(i_E, i_T) % avg_x % xs_sum &
        & = this % prob_tables(i_E, i_T) % avg_x % xs_sum &
        & + this % prob_tables(i_E, i_T) % avg_x % xs_tmp_sum &
        & / dble(histories_avg_urr)

      ! accumulate squared single-batch infinite-dilute means
      this % prob_tables(i_E, i_T) % avg_t % xs_sum2 &
        & = this % prob_tables(i_E, i_T) % avg_t % xs_sum2 &
        & + (this % prob_tables(i_E, i_T) % avg_t % xs_tmp_sum &
        & / dble(histories_avg_urr))&
        & * (this % prob_tables(i_E, i_T) % avg_t % xs_tmp_sum &
        & / dble(histories_avg_urr))
      this % prob_tables(i_E, i_T) % avg_n % xs_sum2 &
        & = this % prob_tables(i_E, i_T) % avg_n % xs_sum2 &
        & + (this % prob_tables(i_E, i_T) % avg_n % xs_tmp_sum &
        & / dble(histories_avg_urr))&
        & * (this % prob_tables(i_E, i_T) % avg_n % xs_tmp_sum &
        & / dble(histories_avg_urr))
      this % prob_tables(i_E, i_T) % avg_g % xs_sum2 &
        & = this % prob_tables(i_E, i_T) % avg_g % xs_sum2 &
        & + (this % prob_tables(i_E, i_T) % avg_g % xs_tmp_sum &
        & / dble(histories_avg_urr))&
        & * (this % prob_tables(i_E, i_T) % avg_g % xs_tmp_sum &
        & / dble(histories_avg_urr))
      this % prob_tables(i_E, i_T) % avg_f % xs_sum2 &
        & = this % prob_tables(i_E, i_T) % avg_f % xs_sum2 &
        & + (this % prob_tables(i_E, i_T) % avg_f % xs_tmp_sum &
        & / dble(histories_avg_urr))&
        & * (this % prob_tables(i_E, i_T) % avg_f % xs_tmp_sum &
        & / dble(histories_avg_urr))
      this % prob_tables(i_E, i_T) % avg_x % xs_sum2 &
        & = this % prob_tables(i_E, i_T) % avg_x % xs_sum2 &
        & + (this % prob_tables(i_E, i_T) % avg_x % xs_tmp_sum &
        & / dble(histories_avg_urr))&
        & * (this % prob_tables(i_E, i_T) % avg_x % xs_tmp_sum &
        & / dble(histories_avg_urr))

    end do

  end subroutine accum_batch

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! CALC_STATS computes batch-based means and standard errors of those means
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine calc_stats(this, i_bat_int)

    class(Isotope), intent(inout), target :: this ! cross section object
    type(ProbabilityTable), pointer :: ptable => null() ! prob table pointer
    integer :: i_E       ! tabulated URR energy index
    integer :: i_T       ! temperature index
    integer :: i_b       ! probability band index
    integer :: i_bat_int ! current batch index
    real(8) :: i_bat ! current batch index for calcs

    i_bat = dble(i_bat_int)

    do i_E = 1, this % ntabs
      if (this % E /= this % Etabs(i_E)) cycle
      do i_T = 1, this % nT
        ptable => this % prob_tables(i_E, i_T)
        do i_b = 1, n_bands
          if (ptable % t(i_b) % cnt > 0) then
            ptable % t(i_b) % xs_mean &
              & = ptable % t(i_b) % xs_sum / ptable % t(i_b) % cnt
            ptable % n(i_b) % xs_mean &
              & = ptable % n(i_b) % xs_sum / ptable % t(i_b) % cnt
            ptable % g(i_b) % xs_mean &
              & = ptable % g(i_b) % xs_sum / ptable % t(i_b) % cnt
            ptable % f(i_b) % xs_mean &
              & = ptable % f(i_b) % xs_sum / ptable % t(i_b) % cnt
            ptable % x(i_b) % xs_mean &
              & = ptable % x(i_b) % xs_sum / ptable % t(i_b) % cnt
          end if
          ptable % t(i_b) % cnt_mean &
            & = ptable % t(i_b) % cnt / (i_bat * histories_avg_urr)
          ptable % n(i_b) % cnt_mean &
            & = ptable % n(i_b) % cnt / (i_bat * histories_avg_urr)
          ptable % g(i_b) % cnt_mean &
            & = ptable % g(i_b) % cnt / (i_bat * histories_avg_urr)
          ptable % f(i_b) % cnt_mean &
            & = ptable % f(i_b) % cnt / (i_bat * histories_avg_urr)
          ptable % x(i_b) % cnt_mean &
            & = ptable % x(i_b) % cnt / (i_bat * histories_avg_urr)
        end do

        ptable % avg_t % xs_mean &
          & = ptable % avg_t % xs_sum / i_bat
        ptable % avg_n % xs_mean &
          & = ptable % avg_n % xs_sum / i_bat
        ptable % avg_g % xs_mean &
          & = ptable % avg_g % xs_sum / i_bat
        ptable % avg_f % xs_mean &
          & = ptable % avg_f % xs_sum / i_bat
        ptable % avg_x % xs_mean &
          & = ptable % avg_x % xs_sum / i_bat

        if (ptable % avg_t % xs_mean /= ZERO &
          & .and. i_bat_int > 1) then
          ! compute standard errors of mean xs values
          ptable % avg_t % xs_sem &
            & = sqrt((ONE / (i_bat - ONE)) &
            & * (ptable % avg_t % xs_sum2 / i_bat &
            & - (ptable % avg_t % xs_mean) &
            & * (ptable % avg_t % xs_mean)))

          ! compute relative uncertainties in xs values
          ptable % avg_t % rel_unc &
            & = ptable % avg_t % xs_sem &
            & / ptable % avg_t % xs_mean
        end if
        if (ptable % avg_n % xs_mean /= ZERO &
          & .and. i_bat_int > 1) then
          ! compute standard errors of mean xs values
          ptable % avg_n % xs_sem &
            & = sqrt((ONE / (i_bat - ONE)) &
            & * (ptable % avg_n % xs_sum2 / i_bat &
            & - (ptable % avg_n % xs_mean) &
            & * (ptable % avg_n % xs_mean)))

          ! compute relative uncertainties in xs values
          ptable % avg_n % rel_unc &
            & = ptable % avg_n % xs_sem &
            & / ptable % avg_n % xs_mean
        end if
        if (ptable % avg_g % xs_mean /= ZERO &
          & .and. i_bat_int > 1) then
          ! compute standard errors of mean xs values
          ptable % avg_g % xs_sem &
            & = sqrt((ONE / (i_bat - ONE)) &
            & * (ptable % avg_g % xs_sum2 / i_bat &
            & - (ptable % avg_g % xs_mean) &
            & * (ptable % avg_g % xs_mean)))

          ! compute relative uncertainties in xs values
          ptable % avg_g % rel_unc &
            & = ptable % avg_g % xs_sem &
            & / ptable % avg_g % xs_mean
        end if
        if (ptable % avg_f % xs_mean > ZERO &
          & .and. i_bat_int > 1 .and. ptable % avg_f % xs_sum2 / i_bat &
          & - ptable % avg_f % xs_mean * ptable % avg_f % xs_mean > ZERO) then
          ! compute standard errors of mean xs values
          ptable % avg_f % xs_sem &
            & = sqrt((ONE / (i_bat - ONE)) &
            & * (ptable % avg_f % xs_sum2 / i_bat &
            & - (ptable % avg_f % xs_mean) &
            & * (ptable % avg_f % xs_mean)))

          ! compute relative uncertainties in xs values
          ptable % avg_f % rel_unc &
            & = ptable % avg_f % xs_sem &
            & / ptable % avg_f % xs_mean
        end if
        if (ptable % avg_x % xs_mean /= ZERO &
          & .and. i_bat_int > 1 .and. competitive) then
          ! compute standard errors of mean xs values
          ptable % avg_x % xs_sem &
            & = sqrt((ONE / (i_bat - ONE)) &
            & * (ptable % avg_x % xs_sum2 / i_bat &
            & - (ptable % avg_x % xs_mean) &
            & * (ptable % avg_x % xs_mean)))

          ! compute relative uncertainties in xs values
          ptable % avg_x % rel_unc &
            & = ptable % avg_x % xs_sem &
            & / ptable % avg_x % xs_mean
        end if
      end do
    end do

  end subroutine calc_stats

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! INTERP_FACTOR computes an interpolation factor
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  function interp_factor(val, val_low, val_up, scheme) result(factor)

    integer :: scheme ! interpolation scheme
    real(8) :: val     ! value we're interpolating to
    real(8) :: val_low ! lower bounding value
    real(8) :: val_up  ! upper bounding value
    real(8) :: factor  ! interpolation factor

    select case(scheme)
    case(LINEAR_LINEAR)

      factor = (val - val_low) / (val_up - val_low)

    case(LOG_LOG)

      if (min(val, val_low, val_up) <= ZERO) &
        & call fatal_error('Argument of log() non-positive')
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

    integer :: scheme ! interpolation scheme
    real(8) :: factor  ! interpolation factor
    real(8) :: val_low ! lower bounding value
    real(8) :: val_up  ! upper bounding value
    real(8) :: val     ! interpolated value

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

  function E_last_rrr(iso, l_val, J_val) result(E_val)

    type(Isotope), pointer :: tope => null() ! nuclide pointer
    integer :: iso   ! isotope index
    integer :: i_res ! RRR resonance index for a given l
    integer :: l_val ! orbital quantum number
    real(8) :: J_val ! total angular momentum quantum number
    real(8) :: E_val ! highest-energy RRR resonance energy for (l,J)
!$omp threadprivate(tope)

    tope => isotopes(iso)

    select case (tope % LRF(tope % i_urr - 1))
    case (SLBW)
      do i_res = size(tope % slbw_resonances(l_val + 1) % E_lam), 1, -1
        if (tope % slbw_resonances(l_val + 1) % AJ(i_res) == J_val &
          .and. tope % slbw_resonances(l_val + 1) % E_lam(i_res) < tope % EL(tope % i_urr)) then
          E_val = tope % slbw_resonances(l_val + 1) % E_lam(i_res)
          exit
        end if
      end do
    case (MLBW)
      do i_res = size(tope % mlbw_resonances(l_val + 1) % E_lam), 1, -1
        if (tope % mlbw_resonances(l_val + 1) % AJ(i_res) == J_val &
          .and. tope % mlbw_resonances(l_val + 1) % E_lam(i_res) < tope % EL(tope % i_urr)) then
          E_val = tope % mlbw_resonances(l_val + 1) % E_lam(i_res)
          exit
        end if
      end do
    case (REICH_MOORE)
      do i_res = size(tope % rm_resonances(l_val + 1) % E_lam), 1, -1
        if (tope % rm_resonances(l_val + 1) % AJ(i_res) == J_val &
          .and. tope % rm_resonances(l_val + 1) % E_lam(i_res) < tope % EL(tope % i_urr)) then
          E_val = tope % rm_resonances(l_val + 1) % E_lam(i_res)
          exit
        end if
      end do
    case default
      call fatal_error('Unrecognized/unsupported RRR formalism')
    end select

  end function E_last_rrr

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! RRR_RES finds the index of the RRR resonance which we need to add the
! contribution of to a URR xs
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  function rrr_res(iso, n_rrr_res, l_val, J_val) result(i_res)

    type(Isotope), pointer :: tope => null() ! nuclide pointer
    integer :: iso       ! isotope index
    integer :: n_rrr_res ! how many RRR resonances to go back
    integer :: cnt_res   ! how many RRR resonances have we gone back
    integer :: i_res     ! index of the RRR resonance cnt_res resonances back
    integer :: l_val     ! orbital quantum number
    real(8) :: J_val ! total angular momentum quantum number
!$omp threadprivate(tope)

    tope => isotopes(iso)

    cnt_res   = 0

    select case (tope % LRF(tope % i_urr - 1))
    case (SLBW)
      do i_res = size(tope % slbw_resonances(l_val + 1) % E_lam), 1, -1
        if (tope % slbw_resonances(l_val + 1) % AJ(i_res) == J_val &
          .and. tope % slbw_resonances(l_val + 1) % E_lam(i_res) < tope % EL(tope % i_urr)) then
          cnt_res = cnt_res + 1
        end if
        if (cnt_res == n_rrr_res) exit
      end do
    case (MLBW)
      do i_res = size(tope % mlbw_resonances(l_val + 1) % E_lam), 1, -1
        if (tope % mlbw_resonances(l_val + 1) % AJ(i_res) == J_val &
          .and. tope % mlbw_resonances(l_val + 1) % E_lam(i_res) < tope % EL(tope % i_urr)) then
          cnt_res = cnt_res + 1
        end if
        if (cnt_res == n_rrr_res) exit
      end do
    case (REICH_MOORE)
      do i_res = size(tope % rm_resonances(l_val + 1) % E_lam), 1, -1
        if (tope % rm_resonances(l_val + 1) % AJ(i_res) == J_val &
          .and. tope % rm_resonances(l_val + 1) % E_lam(i_res) < tope % EL(tope % i_urr)) then
          cnt_res = cnt_res + 1
        end if
        if (cnt_res == n_rrr_res) exit
      end do
    case default
      call fatal_error('Unrecognized/unsupported RRR formalism')
    end select

  end function rrr_res

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! ACCUM_RESONANCES accumulates contribution from an additional resonance
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine accum_resonances(res, t, n, g, f, x)

    type(Resonance) :: res ! resonance object
    type(CrossSection) :: t ! total xs object
    type(CrossSection) :: n ! elastic scattering xs object
    type(CrossSection) :: g ! radiative capture xs object
    type(CrossSection) :: f ! fission xs object
    type(CrossSection) :: x ! competitive inelastic scattering xs object

    call t % accum_resonance(res % dxs_t)
    call n % accum_resonance(res % dxs_n)
    call g % accum_resonance(res % dxs_g)
    call f % accum_resonance(res % dxs_f)
    call x % accum_resonance(res % dxs_x)

  end subroutine accum_resonances

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! RES_CONTRIB add the contribution from an additional resonance
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine res_contrib(this, res, i_E, i_T)

    class(Isotope), intent(inout) :: this ! isotope object
    type(Resonance) :: res ! resonance object
    integer :: i_E ! tabulated URR energy index
    integer :: i_T ! temperature index

    call this % prob_tables(i_E, i_T) % avg_t % accum_resonance(res % dxs_t)
    call this % prob_tables(i_E, i_T) % avg_n % accum_resonance(res % dxs_n)
    call this % prob_tables(i_E, i_T) % avg_g % accum_resonance(res % dxs_g)
    call this % prob_tables(i_E, i_T) % avg_f % accum_resonance(res % dxs_f)
    call this % prob_tables(i_E, i_T) % avg_x % accum_resonance(res % dxs_x)

  end subroutine res_contrib

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! FLUSH_BATCHES zeroes out probability table batch accumulators
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine flush_batches(this)

    class(Isotope), intent(inout) :: this ! isotope object
    integer :: i_E ! tabulated URR energy index
    integer :: i_T ! temperature index
    integer :: i_b ! probability band index

    do i_E = 1, this % ntabs
      do i_T = 1, this % nT
        do i_b = 1, n_bands
          this % prob_tables(i_E, i_T) % t(i_b) % xs_tmp_sum = ZERO
          this % prob_tables(i_E, i_T) % n(i_b) % xs_tmp_sum = ZERO
          this % prob_tables(i_E, i_T) % g(i_b) % xs_tmp_sum = ZERO
          this % prob_tables(i_E, i_T) % f(i_b) % xs_tmp_sum = ZERO
          this % prob_tables(i_E, i_T) % x(i_b) % xs_tmp_sum = ZERO
          this % prob_tables(i_E, i_T) % t(i_b) % cnt_tmp_sum = 0
          this % prob_tables(i_E, i_T) % n(i_b) % cnt_tmp_sum = 0
          this % prob_tables(i_E, i_T) % g(i_b) % cnt_tmp_sum = 0
          this % prob_tables(i_E, i_T) % f(i_b) % cnt_tmp_sum = 0
          this % prob_tables(i_E, i_T) % x(i_b) % cnt_tmp_sum = 0
        end do
        this % prob_tables(i_E, i_T) % avg_t % xs_tmp_sum = ZERO
        this % prob_tables(i_E, i_T) % avg_n % xs_tmp_sum = ZERO
        this % prob_tables(i_E, i_T) % avg_g % xs_tmp_sum = ZERO
        this % prob_tables(i_E, i_T) % avg_f % xs_tmp_sum = ZERO
        this % prob_tables(i_E, i_T) % avg_x % xs_tmp_sum = ZERO
        this % prob_tables(i_E, i_T) % avg_t % cnt_tmp_sum = 0
        this % prob_tables(i_E, i_T) % avg_n % cnt_tmp_sum = 0
        this % prob_tables(i_E, i_T) % avg_g % cnt_tmp_sum = 0
        this % prob_tables(i_E, i_T) % avg_f % cnt_tmp_sum = 0
        this % prob_tables(i_E, i_T) % avg_x % cnt_tmp_sum = 0
      end do
    end do

  end subroutine flush_batches

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! FLUSH_SIGMAS flushes cross section object histories
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine flush_sigmas(t, n, g, f, x)

    type(CrossSection) :: t ! total xs object
    type(CrossSection) :: n ! elastic scattering xs object
    type(CrossSection) :: g ! radiative capture xs object
    type(CrossSection) :: f ! fission xs object
    type(CrossSection) :: x ! competitive inelastic scattering xs object

    t % xs = ZERO
    n % xs = ZERO
    g % xs = ZERO
    f % xs = ZERO
    x % xs = ZERO

  end subroutine flush_sigmas

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! ADD_PARAMETERS adds the URR resonance parameters for a single URR resonance to
! the realization
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine add_parameters(res, iso, i_ens, i_res, i_l, i_J)

    type(Isotope), pointer :: tope => null() ! isotope object pointer
    type(Resonance) :: res ! resonance object
    integer :: iso   ! isotope index
    integer :: i_ens ! resonance ensemble index
    integer :: i_res ! resonance counter
    integer :: i_l   ! orbital quantum number index
    integer :: i_J   ! total angular momentum quantum number
!$omp threadprivate(tope)

    tope => isotopes(iso)

    tope % urr_resonances_tmp(i_ens, i_l, i_J) % E_lam(i_res) = res % E_lam
    tope % urr_resonances_tmp(i_ens, i_l, i_J) % GN(i_res)    = res % Gam_n
    tope % urr_resonances_tmp(i_ens, i_l, i_J) % GG(i_res)    = res % Gam_g
    tope % urr_resonances_tmp(i_ens, i_l, i_J) % GF(i_res)    = res % Gam_f
    tope % urr_resonances_tmp(i_ens, i_l, i_J) % GX(i_res)    = res % Gam_x
    tope % urr_resonances_tmp(i_ens, i_l, i_J) % GT(i_res)    = res % Gam_t

  end subroutine add_parameters

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! SET_PARAMETERS sets the URR resonance parameters for a single resonance
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine set_parameters(res, iso, i_res, i_l, i_J, i_ER)

    type(Isotope), pointer :: tope => null() ! isotope object pointer
    type(Resonance) :: res ! resonance object
    integer :: iso   ! isotope index
    integer :: i_res ! resonance counter
    integer :: i_l   ! orbital quantum number index
    integer :: i_J   ! total angular momentum quantum number
    integer :: i_ER  ! resonance energy region index
!$omp threadprivate(tope)

    tope => isotopes(iso)

    select case (tope % LRU(i_ER))
    ! resolved parameters
    case (1)
      select case (tope % LRF(i_ER))
      case (SLBW)
        res % E_lam = tope % slbw_resonances(i_l) % E_lam(i_res)
        res % Gam_n = tope % slbw_resonances(i_l) % GN(i_res)
        res % Gam_g = tope % slbw_resonances(i_l) % GG(i_res)
        res % Gam_f = tope % slbw_resonances(i_l) % GF(i_res)
        res % Gam_t = tope % slbw_resonances(i_l) % GT(i_res)
        res % Gam_x = res % Gam_t - res % Gam_n - res % Gam_g - res % Gam_f

      case (MLBW)
        res % E_lam = tope % mlbw_resonances(i_l) % E_lam(i_res)
        res % Gam_n = tope % mlbw_resonances(i_l) % GN(i_res)
        res % Gam_g = tope % mlbw_resonances(i_l) % GG(i_res)
        res % Gam_f = tope % mlbw_resonances(i_l) % GF(i_res)
        res % Gam_t = tope % mlbw_resonances(i_l) % GT(i_res)
        res % Gam_x = res % Gam_t - res % Gam_n - res % Gam_g - res % Gam_f

      case (REICH_MOORE)
        res % E_lam = tope % rm_resonances(i_l) % E_lam(i_res)
        res % Gam_n = tope % rm_resonances(i_l) % GN(i_res)
        res % Gam_g = tope % rm_resonances(i_l) % GG(i_res)
        res % Gam_f = tope % rm_resonances(i_l) % GFA(i_res) &
          & + tope % rm_resonances(i_l) % GFB(i_res)
        res % Gam_x = ZERO
        res % Gam_t = res % Gam_n + res % Gam_g + res % Gam_f + res % Gam_x

      case default
        call fatal_error('Unrecognized resolved resonance region formalism')
      end select

    ! unresolved parameters
    case (2)
      res % E_lam = tope % urr_resonances(i_real, i_l, i_J) % E_lam(i_res)
      res % Gam_n = tope % urr_resonances(i_real, i_l, i_J) % GN(i_res)
      res % Gam_g = tope % urr_resonances(i_real, i_l, i_J) % GG(i_res)
      res % Gam_f = tope % urr_resonances(i_real, i_l, i_J) % GF(i_res)
      res % Gam_x = tope % urr_resonances(i_real, i_l, i_J) % GX(i_res)
      res % Gam_t = tope % urr_resonances(i_real, i_l, i_J) % GT(i_res)

    case default
      call fatal_error('Only 1 and 2 are supported ENDF-6 LRU values')
    end select

    tope % local_realization(i_l, i_J) % E_lam(res % i_res) = res % E_lam
    tope % local_realization(i_l, i_J) % Gam_n(res % i_res) = res % Gam_n
    tope % local_realization(i_l, i_J) % Gam_g(res % i_res) = res % Gam_g
    tope % local_realization(i_l, i_J) % Gam_f(res % i_res) = res % Gam_f
    tope % local_realization(i_l, i_J) % Gam_x(res % i_res) = res % Gam_x
    tope % local_realization(i_l, i_J) % Gam_t(res % i_res) = res % Gam_t

  end subroutine set_parameters

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! SET_MEAN_PARAMETERS sets the URR mean resonance parameters at an energy
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine set_mean_parameters(iso, E_res, i_l, i_J)

    type(Isotope), pointer :: tope => null() ! isotope object pointer
    integer :: iso ! isotope index
    integer :: i_E ! tabulated URR parameters energy index
    integer :: i_l ! orbital quantum number index
    integer :: i_J ! total angular momentum quantum number
    real(8) :: E_res ! current resonance (lab) energy (e.g. E_lam)
    real(8) :: fendf ! ENDF6 URR parameters energy interpolation factor
!$omp threadprivate(tope)

    tope => isotopes(iso)

    ! compute interpolation factor
    if (E_res < tope % ES(1)) then
      i_E = 1
    else if (E_res > tope % ES(tope % NE)) then
      i_E = tope % NE - 1
    else
      i_E = binary_search(tope % ES, tope % NE, E_res)
    end if
    fendf = interp_factor(E_res, tope % ES(i_E), tope % ES(i_E + 1),&
      tope % INT)

    ! set current mean unresolved resonance parameters
    tope % D   = interpolator(fendf, &
      & tope % D_mean(i_l) % data(i_J) % data(i_E), &
      & tope % D_mean(i_l) % data(i_J) % data(i_E + 1), tope % INT)
    tope % GN0 = interpolator(fendf, &
      & tope % GN0_mean(i_l) % data(i_J) % data(i_E), &
      & tope % GN0_mean(i_l) % data(i_J) % data(i_E + 1), tope % INT)
    tope % GG  = interpolator(fendf, &
      & tope % GG_mean(i_l) % data(i_J) % data(i_E), &
      & tope % GG_mean(i_l) % data(i_J) % data(i_E + 1), tope % INT)
    if (tope % INT == LINEAR_LINEAR .or. &
      & tope % GF_mean(i_l) % data(i_J) % data(i_E) > ZERO) then
      tope % GF  = interpolator(fendf, &
        & tope % GF_mean(i_l) % data(i_J) % data(i_E), &
        & tope % GF_mean(i_l) % data(i_J) % data(i_E + 1), tope % INT)
    else
      tope % GF = ZERO
    end if
    if (tope % INT == LINEAR_LINEAR .or. &
      & tope % GX_mean(i_l) % data(i_J) % data(i_E) > ZERO) then
      tope % GX  = interpolator(fendf, &
        & tope % GX_mean(i_l) % data(i_J) % data(i_E), &
        & tope % GX_mean(i_l) % data(i_J) % data(i_E + 1), tope % INT)
    else
      tope % GX = ZERO
    end if

  end subroutine set_mean_parameters

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! INTERP_AVG_URR_XS interpolates the averaged, infinite-dilute URR cross
! sections computed via Monte Carlo from mean resonance parameters (i.e. not
! the evaluator-supplied File 3 background cross sections)
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine interp_avg_urr_xs(f, iso, iavg, n_xs, f_xs, g_xs, x_xs)

    type(Isotope), pointer, save :: tope => null() ! isotope object pointer
    integer :: iso
    integer :: iavg ! averaged cross section index
    real(8) :: f
    real(8) :: n_xs
    real(8) :: f_xs
    real(8) :: g_xs
    real(8) :: x_xs
!$omp threadprivate(tope)

    tope => isotopes(iso)

    ! infinite-dilute elastic scattering
    if (tope % avg_urr_n(iavg) > ZERO) then
      n_xs = interpolator(f, &
        & tope % avg_urr_n(iavg), tope % avg_urr_n(iavg + 1), tope % INT)
    else
      call fatal_error('Non-positive (n,n) infinite-dilute cross section')
    end if

    ! infinite-dilute fission
    if (tope % INT == LINEAR_LINEAR .or. tope % avg_urr_f(iavg) > ZERO) then
      f_xs = interpolator(f, &
        & tope % avg_urr_f(iavg), tope % avg_urr_f(iavg + 1), tope % INT)
    else
      f_xs = ZERO
    end if

    ! infinite-dilute capture
    if (tope % INT == LINEAR_LINEAR .or. tope % avg_urr_g(iavg) > ZERO) then
      g_xs = interpolator(f, &
        & tope % avg_urr_g(iavg), tope % avg_urr_g(iavg + 1), tope % INT)
    else
      g_xs = ZERO
    end if

    ! infinite-dilute competitive reaction xs
    if (tope % INT == LINEAR_LINEAR .or. tope % avg_urr_x(iavg) > ZERO) then
      x_xs = interpolator(f, &
        & tope % avg_urr_x(iavg), tope % avg_urr_x(iavg + 1), tope % INT)
    else
      x_xs = ZERO
    end if

  end subroutine interp_avg_urr_xs

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! INTERP_ACE_BACKGROUND interpolates the processed ENDF-6 File 3 background
! cross sections
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine interp_ace_background(iso, i_nuc, n_pts, f, i_grid)

    type(Isotope), pointer :: tope => null() ! isotope object pointer
    type(Nuclide), pointer :: nuc => null() ! nuclide object pointer
    integer :: iso
    integer :: i_nuc
    integer :: n_pts
    integer :: i_grid
    real(8) :: f
!$omp threadprivate(tope, nuc)

    tope => isotopes(iso)
    nuc  => nuclides(i_nuc)

    ! elastic scattering xs
    tope % n_tmp(n_pts) = interpolator(f, &
      & nuc % elastic(i_grid), nuc % elastic(i_grid + 1), tope % INT)

    ! radiative capture xs
    tope % g_tmp(n_pts) = interpolator(f, &
      & nuc % absorption(i_grid) - nuc % fission(i_grid), &
      & nuc % absorption(i_grid + 1) - nuc % fission(i_grid + 1), tope % INT)

    ! fission xs
    if (tope % INT == LINEAR_LINEAR .or. nuc % fission(i_grid) > ZERO) then
      tope % f_tmp(n_pts) = interpolator(f, nuc % fission(i_grid), &
        & nuc % fission(i_grid + 1), tope % INT)
    else
      tope % f_tmp(n_pts) = ZERO
    end if

    ! competitive first level inelastic scattering xs
    if (tope % INT == LINEAR_LINEAR .or. nuc % total(i_grid) &
      & - nuc % absorption(i_grid) - nuc % elastic(i_grid) > ZERO) then
      tope % x_tmp(n_pts) = interpolator(f, &
        &   nuc % total(i_grid) - nuc % absorption(i_grid) &
        & - nuc % elastic(i_grid), &
        &   nuc % total(i_grid + 1) - nuc % absorption(i_grid + 1) &
        & - nuc % elastic(i_grid + 1), tope % INT)
    else
      tope % x_tmp(n_pts) = ZERO
    end if

    ! total xs
    tope % t_tmp(n_pts) = tope % n_tmp(n_pts) &
      & + tope % g_tmp(n_pts) + tope % f_tmp(n_pts)&
      & + tope % x_tmp(n_pts)

  end subroutine interp_ace_background

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! ADD_MF3_BACKGROUND adds the resonance xs component to the evaluator-supplied
! background
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine add_mf3_background(iso, i_nuc, n_pts, t, n, g, f, x)

    type(Isotope), pointer :: tope => null() ! isotope object pointer
    type(CrossSection) :: t ! total xs object
    type(CrossSection) :: n ! elastic scattering xs object
    type(CrossSection) :: g ! radiative capture xs object
    type(CrossSection) :: f ! fission xs object
    type(CrossSection) :: x ! competitive inelastic scattering xs object
    integer :: iso    ! isotope index
    integer :: i_nuc  ! nuclide index
    integer :: i_grid ! background energy grid index
    integer :: n_pts  ! URR pointwise xs energy grid index
    real(8) :: fact         ! File 3 interpolation factor
    real(8) :: capture_xs   ! radiative capture xs
    real(8) :: inelastic_xs ! first level inelastic scattering xs
!$omp threadprivate(tope)

    tope => isotopes(iso)

    if (tope % point_urr_xs) then

      ! elastic scattering xs
      tope % n_tmp(n_pts) = tope % n_tmp(n_pts) &
        & + n % xs

      ! radiative capture xs
      tope % g_tmp(n_pts) = tope % g_tmp(n_pts) &
        & + g % xs

      ! fission xs
      tope % f_tmp(n_pts) = tope % f_tmp(n_pts) &
        & + f % xs

      ! competitive first level inelastic scattering xs
      tope % x_tmp(n_pts) = tope % x_tmp(n_pts) &
        & + x % xs

      ! total xs
      tope % t_tmp(n_pts) = tope % t_tmp(n_pts) &
        & + t % xs

    else

      ! elastic scattering xs
      if (tope % E < tope % MF3_n_e(1)) then
        call fatal_error('Energy is below File 3 elastic energy grid')
      else if (tope % E > tope % MF3_n_e(size(tope % MF3_n_e))) then
        call fatal_error('Energy is above File 3 elastic energy grid')
      else
        i_grid = binary_search(tope % MF3_n_e, size(tope % MF3_n_e),tope % E)
        if (tope % INT == LINEAR_LINEAR .or. &
          & tope % MF3_n_e(i_grid) > ZERO) then
          fact = interp_factor(tope % E, tope % MF3_n_e(i_grid), &
            & tope % MF3_n_e(i_grid + 1), tope % INT)
          micro_xs(i_nuc) % elastic = interpolator(fact, tope % MF3_n(i_grid),&
            & tope % MF3_n(i_grid + 1), tope % INT) + n % xs
        else
          micro_xs(i_nuc) % elastic = ZERO
        end if
      end if
      if (micro_xs(i_nuc) % elastic < ZERO) micro_xs(i_nuc) % elastic = ZERO

      ! competitive reaction xs
      if (tope % E < tope % MF3_x_e(1)) then
        inelastic_xs = ZERO
      else if (tope % E > tope % MF3_x_e(size(tope % MF3_x_e))) then
        call fatal_error('Energy is above File 3 inelastic energy grid')
      else
        i_grid = binary_search(tope % MF3_x_e, size(tope % MF3_x_e), tope % E)
        if (tope % INT == LINEAR_LINEAR .or. &
          & tope % MF3_x_e(i_grid) > ZERO) then
          fact = interp_factor(tope % E, tope % MF3_x_e(i_grid), &
            & tope % MF3_x_e(i_grid + 1), tope % INT)
          inelastic_xs = interpolator(fact, tope % MF3_x(i_grid), &
            & tope % MF3_x(i_grid + 1), tope % INT)
        else
          inelastic_xs = ZERO
        end if
        if (competitive) inelastic_xs = inelastic_xs + x % xs
      end if
      if (inelastic_xs < ZERO) inelastic_xs = ZERO

      ! capture xs
      if (tope % E < tope % MF3_g_e(1)) then
        call fatal_error('Energy is below File 3 capture energy grid')
      else if (tope % E > tope % MF3_g_e(size(tope % MF3_g_e))) then
        call fatal_error('Energy is above File 3 capture energy grid')
      else
        i_grid = binary_search(tope % MF3_g_e, size(tope % MF3_g_e), tope % E)
        if (tope % INT == LINEAR_LINEAR .or. &
          & tope % MF3_g_e(i_grid) > ZERO) then
          fact = interp_factor(tope % E, tope % MF3_g_e(i_grid), &
            & tope % MF3_g_e(i_grid + 1), tope % INT)
          capture_xs = interpolator(fact, tope % MF3_g(i_grid), &
            & tope % MF3_g(i_grid + 1), tope % INT) + g % xs
        else
          capture_xs = g % xs
        end if
      end if

      ! fission xs
      if (.not. (allocated(tope % MF3_f_e))) then
        micro_xs(i_nuc) % fission = ZERO
      else
        if (tope % E < tope % MF3_f_e(1)) then
          micro_xs(i_nuc) % fission = ZERO
        else if (tope % E > tope % MF3_f_e(size(tope % MF3_f_e))) then
          call fatal_error('Energy is above File 3 fission energy grid')
        else
          i_grid = binary_search(tope % MF3_f_e, size(tope % MF3_f_e), tope%E)
          if (tope % INT == LINEAR_LINEAR .or. &
            & tope % MF3_f_e(i_grid) > ZERO) then
            fact = interp_factor(tope % E, tope % MF3_f_e(i_grid), &
              & tope % MF3_f_e(i_grid + 1), tope % INT)
            micro_xs(i_nuc) % fission = interpolator(fact, tope%MF3_f(i_grid),&
              & tope % MF3_f(i_grid + 1), tope % INT) + f % xs
          else
            micro_xs(i_nuc) % fission = f % xs
          end if
        end if
      end if

      ! absorption xs
      micro_xs(i_nuc) % absorption = micro_xs(i_nuc) % fission + capture_xs

      ! total xs
      micro_xs(i_nuc) % total = micro_xs(i_nuc) % elastic &
           + micro_xs(i_nuc) % absorption + inelastic_xs
    end if

  end subroutine add_mf3_background

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! ALLOC_ENERGY_RANGE allocates variables for the resonance energy ranges
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine alloc_energy_range(this)

    class(Isotope), intent(inout) :: this ! isotope object

    allocate(this % NLS(this % NER))
    allocate(this % LRU(this % NER))
    allocate(this % LRF(this % NER))
    allocate(this % NRO(this % NER))
    allocate(this % NAPS(this % NER))
    allocate(this % EL(this % NER))
    allocate(this % EH(this % NER))
    allocate(this % AP(this % NER))
    allocate(this % ac(this % NER))
    allocate(this % SPI(this % NER))

  end subroutine alloc_energy_range

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! DEALLOC_ENERGY_RANGES deallocates variables for the resonance energy ranges
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine dealloc_energy_ranges(this)

    class(Isotope), intent(inout) :: this ! isotope object

    deallocate(this % NLS)
    deallocate(this % LRU)
    deallocate(this % LRF)
    deallocate(this % NRO)
    deallocate(this % NAPS)
    deallocate(this % EL)
    deallocate(this % EH)
    deallocate(this % AP)
    deallocate(this % ac)
    deallocate(this % SPI)

  end subroutine dealloc_energy_ranges

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! DEALLOC_MF3_URR deallocates ENDF-6 File 3 evaluator-supplied background
! cross sections
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine dealloc_MF3(this)

    class(Isotope), intent(inout) :: this ! isotope object

    if (allocated(this % MF3_n_e)) deallocate(this % MF3_n_e)
    if (allocated(this % MF3_f_e)) deallocate(this % MF3_f_e)
    if (allocated(this % MF3_g_e)) deallocate(this % MF3_g_e)
    if (allocated(this % MF3_x_e)) deallocate(this % MF3_x_e)
    if (allocated(this % MF3_n)) deallocate(this % MF3_n)
    if (allocated(this % MF3_f)) deallocate(this % MF3_f)
    if (allocated(this % MF3_g)) deallocate(this % MF3_g)
    if (allocated(this % MF3_x)) deallocate(this % MF3_x)

  end subroutine dealloc_MF3

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! ALLOC_ENSEMBLE_TMP allocates temporary URR resonance ensemble realizations
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine alloc_ensemble_tmp(this, n_reals)

    class(Isotope), intent(inout) :: this ! isotope object
    integer :: n_reals ! number of realizations
    integer :: i_ens   ! realization/ensemble index
    integer :: i_l     ! orbital angular momentum quantum number index
    integer :: i_J     ! total angular momentum quantum number index

    ! allocate temporary URR resonances
    allocate(this % urr_resonances_tmp(n_reals, this % NLS(this % i_urr),&
      & maxval(this % NJS(:))))

    ! loop over realizations
    do i_ens = 1, n_reals

      ! loop over orbital quantum numbers
      do i_l = 1, this % NLS(this % i_urr)

        ! loop over total angular momenta
        do i_J = 1, this % NJS(i_l)

          ! allocate resonance parameters
          call this % urr_resonances_tmp(i_ens, i_l, i_J) &
            & % alloc_slbw_resonances(this % n_lam_tmp)
        end do
      end do
    end do

  end subroutine alloc_ensemble_tmp

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! DEALLOC_ENSEMBLE_TMP deallocates temporary URR resonance ensemble realizations
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine dealloc_ensemble_tmp(this, n_reals)

    class(Isotope), intent(inout) :: this ! isotope object
    integer :: n_reals ! number of realizations
    integer :: i_ens   ! realization/ensemble index
    integer :: i_l     ! orbital angular momentum quantum number index
    integer :: i_J     ! total angular momentum quantum number index

    ! loop over realizations
    do i_ens = 1, n_reals

      ! loop over orbital quantum numbers
      do i_l = 1, this % NLS(this % i_urr)

        ! loop over total angular momenta
        do i_J = 1, this % NJS(i_l)

! TODO: deallocate resonances of the proper formalism once MLBW, RM allowed
          ! allocate resonance parameters
          call this % urr_resonances_tmp(i_ens, i_l, i_J) &
            & % dealloc_slbw_resonances()
        end do
      end do
    end do

    ! deallocate temporary URR resonances
    deallocate(this % urr_resonances_tmp)

  end subroutine dealloc_ensemble_tmp

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! ALLOC_ENSEMBLE allocates URR resonance ensemble realizations
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine alloc_ensemble(this, n_reals)

    class(Isotope), intent(inout) :: this ! isotope object
    integer :: n_reals ! number of realizations
    integer :: i_ens   ! realization/ensemble index
    integer :: i_l     ! orbital angular momentum quantum number index
    integer :: i_J     ! total angular momentum quantum number index

    ! allocate URR resonance realizations
    allocate(this % urr_resonances(n_reals, this % NLS(this % i_urr),&
      & maxval(this % NJS(:))))

    ! loop over realizations
    do i_ens = 1, n_reals

      ! loop over orbital quantum numbers
      do i_l = 1, this % NLS(this % i_urr)

        ! loop over total angular momenta
        do i_J = 1, this % NJS(i_l)

          ! allocate URR resonances
          call this % urr_resonances(i_ens, i_l, i_J) &
            & % alloc_slbw_resonances(this % n_lam(i_ens, i_l, i_J))
        end do
      end do
    end do

  end subroutine alloc_ensemble

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! DEALLOC_ENSEMBLE deallocates a URR resonance ensemble realization
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine dealloc_ensemble(this, n_reals)

    class(Isotope), intent(inout) :: this ! isotope object
    integer :: n_reals ! number of realizations
    integer :: i_ens   ! realization/ensemble index
    integer :: i_l     ! orbital angular momentum quantum number index
    integer :: i_J     ! total angular momentum quantum number index

    ! loop over realizations
    do i_ens = 1, n_reals

      ! loop over orbital quantum numbers
      do i_l = 1, this % NLS(this % i_urr)

        ! loop over total angular momenta
        do i_J = 1, this % NJS(i_l)

! TODO: deallocate resonances of the proper formalism once MLBW, RM allowed
          ! deallocate resonance parameters
          call this % urr_resonances(i_ens, i_l, i_J) &
            & % dealloc_slbw_resonances()
        end do
      end do
    end do

    ! deallocate URR resonance realizations
    deallocate(this % urr_resonances)

  end subroutine dealloc_ensemble

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! ALLOC_SLBW_RESONANCES allocates a vector of SLBW resonances for a given (l,J)
! in the URR case, and for a given number of resonances, NRS, in the RRR case
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine alloc_slbw_resonances(this, NRS)

    class(SLBWResonances), intent(inout) :: this ! resonance vector object
    integer :: NRS

    allocate(this % E_lam(NRS))
    allocate(this % AJ(NRS))
    allocate(this % GN(NRS))
    allocate(this % GG(NRS))
    allocate(this % GF(NRS))
    allocate(this % GX(NRS))
    allocate(this % GT(NRS))

  end subroutine alloc_slbw_resonances

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! DEALLOC_SLBW_RESONANCES deallocates a vector of SLBW resonances for a given
! (l,J) in the URR case, and for a given number of resonances, NRS, in the RRR
! case
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine dealloc_slbw_resonances(this)

    class(SLBWResonances), intent(inout) :: this ! resonance vector object

    deallocate(this % E_lam)
    deallocate(this % AJ)
    deallocate(this % GN)
    deallocate(this % GG)
    deallocate(this % GF)
    deallocate(this % GX)
    deallocate(this % GT)

  end subroutine dealloc_slbw_resonances

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! ALLOC_MLBW_RESONANCES allocates a vector of NRS MLBW resonances
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine alloc_mlbw_resonances(this, NRS)

    class(MLBWResonances), intent(inout) :: this ! resonance vector object
    integer :: NRS

    allocate(this % E_lam(NRS))
    allocate(this % AJ(NRS))
    allocate(this % GN(NRS))
    allocate(this % GG(NRS))
    allocate(this % GF(NRS))
    allocate(this % GX(NRS))
    allocate(this % GT(NRS))

  end subroutine alloc_mlbw_resonances

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! DEALLOC_MLBW_RESONANCES deallocates a vector of NRS MLBW resonances
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine dealloc_mlbw_resonances(this)

    class(MLBWResonances), intent(inout) :: this ! resonance vector object

    deallocate(this % E_lam)
    deallocate(this % AJ)
    deallocate(this % GN)
    deallocate(this % GG)
    deallocate(this % GF)
    deallocate(this % GX)
    deallocate(this % GT)

  end subroutine dealloc_mlbw_resonances

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! ALLOC_RM_RESONANCES allocates a vector of NRS Reich-Moore resonances
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine alloc_rm_resonances(this, NRS)

    class(RMResonances), intent(inout) :: this ! resonance vector object
    integer :: NRS

    allocate(this % E_lam(NRS))
    allocate(this % AJ(NRS))
    allocate(this % GN(NRS))
    allocate(this % GG(NRS))
    allocate(this % GFA(NRS))
    allocate(this % GFB(NRS))

  end subroutine alloc_rm_resonances

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! DEALLOC_RM_RESONANCES deallocates a vector of NRS Reich-Moore resonances
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine dealloc_rm_resonances(this)

    class(RMResonances), intent(inout) :: this ! resonance vector object

    deallocate(this % E_lam)
    deallocate(this % AJ)
    deallocate(this % GN)
    deallocate(this % GG)
    deallocate(this % GFA)
    deallocate(this % GFB)

  end subroutine dealloc_rm_resonances

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! ALLOC_LOCAL_SEQUENCE allocates a local spin sequence of URR resonances about
! E_n
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine alloc_local_sequence(this, n_lam)

    class(LocalSequence), intent(inout) :: this ! resonance vector object
    integer :: n_lam

    allocate(this % E_lam(n_lam))
    allocate(this % Gam_n(n_lam))
    allocate(this % Gam_f(n_lam))
    allocate(this % Gam_g(n_lam))
    allocate(this % Gam_x(n_lam))
    allocate(this % Gam_t(n_lam))

  end subroutine alloc_local_sequence

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! DEALLOC_LOCAL_SEQUENCE deallocates a local spin sequence of URR resonances
! about E_n
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine dealloc_local_sequence(this)

    class(LocalSequence), intent(inout) :: this ! resonance vector object

    deallocate(this % E_lam)
    deallocate(this % Gam_n)
    deallocate(this % Gam_f)
    deallocate(this % Gam_g)
    deallocate(this % Gam_x)
    deallocate(this % Gam_t)

  end subroutine dealloc_local_sequence

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! ALLOC_LOCAL_REALIZATION allocates a local realization of URR resonances
! about E_n
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine alloc_local_realization(this)

    class(Isotope), intent(inout) :: this ! isotope object
    integer :: i_l ! orbital angular momentum quantum number index
    integer :: i_J ! total angular momentum quantum number index

    ! allocate URR spin sequences
    allocate(this % local_realization(this % NLS(this % i_urr),&
      maxval(this % NJS(:))))

    ! loop over orbital quantum numbers
    do i_l = 1, this % NLS(this % i_urr)

      ! loop over total angular momenta
      do i_J = 1, this % NJS(i_l)

        ! allocate spin sequence resonances
        call this % local_realization(i_l, i_J) &
          % alloc_local_sequence(l_waves(i_l))
      end do
    end do

  end subroutine alloc_local_realization

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! DEALLOC_LOCAL_REALIZATION deallocates a local realization of URR resonances
! about E_n
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine dealloc_local_realization(this)

    class(Isotope), intent(inout) :: this ! isotope object
    integer :: i_l ! orbital angular momentum quantum number index
    integer :: i_J ! total angular momentum quantum number index

    ! loop over orbital quantum numbers
    do i_l = 1, this % NLS(this % i_urr)

      ! loop over total angular momenta
      do i_J = 1, this % NJS(i_l)

        ! deallocate spin sequence resonances
        call this % local_realization(i_l, i_J) &
          % dealloc_local_sequence()
      end do
    end do

    ! deallocate URR spin sequences
    deallocate(this % local_realization)

  end subroutine dealloc_local_realization

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! ALLOC_POINTWISE_TMP allocates the temporary pointwise URR energy-cross section
! grids
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine alloc_pointwise_tmp(this)

    class(Isotope), intent(inout) :: this ! isotope object

    allocate(this % E_tmp(this % n_urr_gridpoints))
    allocate(this % n_tmp(this % n_urr_gridpoints))
    allocate(this % g_tmp(this % n_urr_gridpoints))
    allocate(this % f_tmp(this % n_urr_gridpoints))
    allocate(this % x_tmp(this % n_urr_gridpoints))
    allocate(this % t_tmp(this % n_urr_gridpoints))

  end subroutine alloc_pointwise_tmp

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! DEALLOC_POINTWISE_TMP deallocates the temporary pointwise URR energy-cross
! section grids
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine dealloc_pointwise_tmp(this)

    class(Isotope), intent(inout) :: this ! isotope object

    deallocate(this % E_tmp)
    deallocate(this % n_tmp)
    deallocate(this % g_tmp)
    deallocate(this % f_tmp)
    deallocate(this % x_tmp)
    deallocate(this % t_tmp)

  end subroutine dealloc_pointwise_tmp

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! ALLOC_POINTWISE allocates the pointwise URR energy-cross section grids
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine alloc_pointwise(this, n_pts)

    class(Isotope), intent(inout) :: this ! isotope object
    integer :: n_pts ! number of points in grid

    allocate(this % urr_E(n_pts))
    allocate(this % urr_n(n_pts))
    allocate(this % urr_g(n_pts))
    allocate(this % urr_f(n_pts))
    allocate(this % urr_x(n_pts))
    allocate(this % urr_t(n_pts))

  end subroutine alloc_pointwise

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! DEALLOC_POINTWISE deallocates the pointwise URR energy-cross section grids
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine dealloc_pointwise(this)

    class(Isotope), intent(inout) :: this ! isotope object

    deallocate(this % urr_E)
    deallocate(this % urr_n)
    deallocate(this % urr_g)
    deallocate(this % urr_f)
    deallocate(this % urr_x)
    deallocate(this % urr_t)

  end subroutine dealloc_pointwise

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! ALLOC_PROB_TABLES allocates the probability tables for this isotope
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine alloc_prob_tables(this)

    class(Isotope), intent(inout) :: this ! isotope object
    integer :: i_E ! tabulated URR energy index
    integer :: i_T ! temperature index

    if (E_spacing == ENDF6) then
      this % ntabs = this % NE
      allocate(this % Etabs(this % ntabs))
      this % Etabs(:) = this % ES
    else if (E_spacing == LOGARITHMIC) then
      this % ntabs = n_tables
      allocate(this % Etabs(this % ntabs))
      this % Etabs(:) = (/(this % EL(this % i_urr) &
        & * exp(i_E * log(this % EH(this % i_urr) / this % EL(this % i_urr)) &
        & / (this % ntabs - 1)), i_E = 0, this % ntabs - 1)/)
    else if (E_spacing == USER) then
      this % ntabs = n_tables
      allocate(this % Etabs(this % ntabs))
      this % Etabs(:) = Etables
    end if

    allocate(this % prob_tables(this % ntabs, this % nT))

    do i_E = 1, this % ntabs
      do i_T = 1, this % nT
        allocate(this % prob_tables(i_E, i_T) % t(n_bands))
        allocate(this % prob_tables(i_E, i_T) % n(n_bands))
        allocate(this % prob_tables(i_E, i_T) % g(n_bands))
        allocate(this % prob_tables(i_E, i_T) % f(n_bands))
        allocate(this % prob_tables(i_E, i_T) % x(n_bands))
      end do
    end do

  end subroutine alloc_prob_tables

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! DEALLOC_PROB_TABLES deallocates the probability tables for this isotope
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine dealloc_prob_tables(this)

    class(Isotope), intent(inout) :: this ! isotope object
    integer :: i_E ! tabulated URR energy index
    integer :: i_T ! temperature index

    do i_E = 1, this % ntabs
      do i_T = 1, this % nT
        deallocate(this % prob_tables(i_E, i_T) % t)
        deallocate(this % prob_tables(i_E, i_T) % n)
        deallocate(this % prob_tables(i_E, i_T) % g)
        deallocate(this % prob_tables(i_E, i_T) % f)
        deallocate(this % prob_tables(i_E, i_T) % x)
      end do
    end do

    deallocate(this % Etabs)
    deallocate(this % prob_tables)

  end subroutine dealloc_prob_tables

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! FLUSH_PTABLE_STATS zeroes out statistics accumulators for an isotope's
! probability tables
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine flush_ptable_stats(this, i_E)

    class(Isotope), intent(inout) :: this ! isotope object
    integer :: i_E ! tabulated URR energy index
    integer :: i_T ! temperature index
    integer :: i_b ! probability band index

    do i_T = 1, this % nT
      do i_b = 1, n_bands
        call this % prob_tables(i_E, i_T) % t(i_b) % flush_xs_stats()
        call this % prob_tables(i_E, i_T) % n(i_b) % flush_xs_stats()
        call this % prob_tables(i_E, i_T) % g(i_b) % flush_xs_stats()
        call this % prob_tables(i_E, i_T) % f(i_b) % flush_xs_stats()
        call this % prob_tables(i_E, i_T) % x(i_b) % flush_xs_stats()
      end do
      call this % prob_tables(i_E, i_T) % avg_t % flush_xs_stats()
      call this % prob_tables(i_E, i_T) % avg_n % flush_xs_stats()
      call this % prob_tables(i_E, i_T) % avg_g % flush_xs_stats()
      call this % prob_tables(i_E, i_T) % avg_f % flush_xs_stats()
      call this % prob_tables(i_E, i_T) % avg_x % flush_xs_stats()
    end do

  end subroutine flush_ptable_stats

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! FLUSH_XS_STATS clears the batch statistics and accumulators
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine flush_xs_stats(this)

    class(CrossSection), intent(inout) :: this ! cross section object

    ! clear mean, SEM, and overall accumulators
    this % xs          = ZERO
    this % xs_tmp_sum  = ZERO
    this % cnt_tmp_sum = 0
    this % xs_sum      = ZERO
    this % xs_sum2     = ZERO
    this % cnt         = 0
    this % cnt2        = 0
    this % xs_mean     = ZERO
    this % cnt_mean    = ZERO
    this % xs_sem      = ZERO
    this % cnt_sem     = ZERO
    this % rel_unc     = ZERO

  end subroutine flush_xs_stats

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! CHANNEL_RADIUS computes or sets the channel radius depending on ENDF flags
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine channel_radius(this, i_ER)

    class(Isotope), intent(inout) :: this ! isotope object
    integer :: i_ER ! resonance energy range index

    select case (this % NRO(i_ER))

    ! scattering radius is independent of energy
    case (0)

      select case (this % NAPS(i_ER))

      ! use channel radius for penetrabilities and shift factors but
      ! scattering radius for phase shifts
      case (0)
        this % ac(i_ER) = 0.123_8 * this % awr**(ONE/THREE) + 0.08_8

      ! use scattering radius for penetrabilities, shift factors and phase
      ! shifts
      case (1)
        this % ac(i_ER) = this % AP(i_ER)

      ! invalid scattering radius treatment flag
      case default
        call fatal_error('ENDF-6 NAPS flag must be 0 or 1 when NRO is 0')
      end select

    ! scattering radius is energy dependent
    case (1)

      select case (this % NAPS(i_ER))

      ! use channel radius for penetrabilities and shift factors but
      ! scattering radius for phase shifts
      case (0)
        this % ac(i_ER) = 0.123_8 * this % awr**(ONE/THREE) + 0.08_8

      ! use scattering radius for penetrabilities, shift factors and phase
      ! shifts
      case (1)
        this % ac(i_ER) = this % AP(i_ER)

! TODO: understand this and implement it correctly
      ! use energy dependent scattering radius in phase shifts but the energy
      ! independent AP scattering radius value for penetrabilities and shift
      ! factors
      case (2)
        this % ac(i_ER) = this % AP(i_ER)

      ! invalid scattering radius treatment flag
      case default
        call fatal_error('ENDF-6 NAPS flag must be 0, 1, or 2 when NRO is 1')
      end select

    ! invalid energy dependence of scattering radius flag
    case default
      call fatal_error('ENDF-6 NRO flag must be 0 or 1')
    end select

  end subroutine channel_radius

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! DEALLOC_ISOTOPE deallocates a URR isotope object
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine dealloc_isotope(this)

    class(Isotope), intent(inout) :: this ! isotope object
    integer :: i_l  ! orbital angular momentum index
    integer :: i_J  ! total angular momentum index

    ! deallocate mean parameters
    do i_l = 1, this % NLS(this % i_urr)
      do i_J = 1, this % NJS(i_l)
        deallocate(this % D_mean(i_l) % data(i_J) % data)
        deallocate(this % GN0_mean(i_l) % data(i_J) % data)
        deallocate(this % GG_mean(i_l) % data(i_J) % data)
        deallocate(this % GF_mean(i_l) % data(i_J) % data)
        deallocate(this % GX_mean(i_l) % data(i_J) % data)
      end do
      deallocate(this % D_mean(i_l) % data)
      deallocate(this % GN0_mean(i_l) % data)
      deallocate(this % GG_mean(i_l) % data)
      deallocate(this % GF_mean(i_l) % data)
      deallocate(this % GX_mean(i_l) % data)
      deallocate(this % AJ(i_l) % data)
      deallocate(this % DOFX(i_l) % data)
      deallocate(this % DOFN(i_l) % data)
      deallocate(this % DOFG(i_l) % data)
      deallocate(this % DOFF(i_l) % data)
    end do
    deallocate(this % D_mean)
    deallocate(this % GN0_mean)
    deallocate(this % GG_mean)
    deallocate(this % GF_mean)
    deallocate(this % GX_mean)
    deallocate(this % AJ)
    deallocate(this % DOFX)
    deallocate(this % DOFN)
    deallocate(this % DOFG)
    deallocate(this % DOFF)
    deallocate(this % ES)

    ! deallocate probability tables
    if (allocated(this % prob_tables)) call this % dealloc_prob_tables()

    ! deallocate URR resonances
    if (allocated(this % n_lam)) deallocate(this % n_lam)
    if (allocated(this % urr_resonances)) call this % dealloc_ensemble(n_reals)

    ! deallocate RRR resonances
    if (allocated(this % slbw_resonances)) then
      do i_l = 1, this % NLS(this % i_urr - 1)
        call this % slbw_resonances(i_l) % dealloc_slbw_resonances()
      end do
      deallocate(this % slbw_resonances)
    end if
    if (allocated(this % mlbw_resonances)) then
      do i_l = 1, this % NLS(this % i_urr - 1)
        call this % mlbw_resonances(i_l) % dealloc_mlbw_resonances()
      end do
      deallocate(this % mlbw_resonances)
    end if
    if (allocated(this % rm_resonances)) then
      do i_l = 1, this % NLS(this % i_urr - 1)
        call this % rm_resonances(i_l) % dealloc_rm_resonances()
      end do
      deallocate(this % rm_resonances)
    end if

    ! deallocate pointwise data
    if (allocated(this % urr_E)) call this % dealloc_pointwise()

    ! deallocate averaged, infinite-dilute URR cross sections
    deallocate(this % Eavg)
    deallocate(this % avg_urr_t)
    deallocate(this % avg_urr_n)
    deallocate(this % avg_urr_f)
    deallocate(this % avg_urr_g)
    deallocate(this % avg_urr_x)

    ! deallocate ENDF-6 File 3 cross sections
    call this % dealloc_MF3()

    ! deallocate energy range variables and total angular momenta counts
    call this % dealloc_energy_ranges()
    deallocate(this % NJS)

  end subroutine dealloc_isotope

end module unresolved
