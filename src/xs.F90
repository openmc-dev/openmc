module xs

  use error,         only: fatal_error
  use faddeeva,      only: quickw, faddeeva_w
  use fission,       only: nu_total
  use global
  use list_header,   only: ListInt, ListReal
  use random_lcg,    only: prn
  use resonances,    only: BWResonanceListVec,&
                           BWResonanceVec,&
                           BWResonanceVecVec,&
                           RMResonanceVec
  use search,        only: binary_search
  use vectors,       only: VecReal, VecVecReal, VecInt

  implicit none
  private
  public :: background,&
            calculate_prob_band_xs,&
            calculate_urr_xs_otf,&
            calc_urr_xs_otf,&
            competitive,&
            endf_files,&
            E_spacing,&
            E_tabs,&
            formalism,&
            histories_avg_urr,&
            INT_T,&
            i_real_user,&
            Isotope,&
            isotopes,&
            load_prob_tables,&
            load_urr_tables,&
            l_waves,&
            max_batches_avg_urr,&
            min_batches_avg_urr,&
            min_dE_point_urr,&
            n_bands,&
            n_isotopes,&
            n_reals,&
            n_energies,&
            n_temperatures,&
            path_avg_urr_xs,&
            path_endf,&
            path_prob_tables,&
            pointwise_urr,&
            calc_prob_tables,&
            real_freq,&
            represent_params,&
            represent_urr,&
            resonance_ensemble,&
            run_fasturr,&
            T_tabs,&
            tol_avg_urr,&
            tol_point_urr,&
            w_eval,&
            wigner_surmise,&
            write_avg_urr_xs,&
            write_urr_tables

  logical :: competitive      ! use competitve reaction xs resonance structure?
  logical :: load_urr_tables  ! load tables or generate new ones?
  logical :: run_fasturr      ! use special treatment for Fast/URR data?
  logical :: write_avg_urr_xs ! write averaged xs values to a file?
  logical :: write_urr_tables ! write probability tables to a file?
  integer :: background          ! where to get background cross sections
  integer :: E_spacing           ! probability table energy spacing scheme
  integer :: formalism           ! URR resonance formalism
  integer :: histories_avg_urr   ! histories for averaged cross section calc
  integer :: INT_T               ! temperature interpolation scheme
  integer :: i_real              ! index of URR realization used for calc
  integer :: i_real_user         ! user-specified realization index
  integer :: l_waves(4)          ! number of contributing l-wave
  integer :: max_batches_avg_urr ! max batches for averaged cross section calc
  integer :: min_batches_avg_urr ! min batches for averaged cross section calc
  integer :: n_bands             ! number of probability table xs bands
  integer :: n_temperatures      ! number of probability table temperatures
  integer :: n_isotopes          ! number of URR isotopes being processed
  integer :: n_energies          ! number of probability tables (energies)
  integer :: n_reals             ! number of independent URR realizations
  integer :: real_freq           ! frequency of URR realizations
  integer :: represent_params    ! representation of URR parameters
  integer :: represent_urr       ! representation of URR cross sections
  integer :: w_eval              ! W function evaluation method
  real(8) :: min_dE_point_urr    ! min diff between reconstructed energies [eV]
  real(8) :: tol_avg_urr         ! max rel err for inf dil xs calc termination
  real(8) :: tol_point_urr       ! max pointwise xs reconstruction rel err
  real(8), allocatable :: E_tabs(:) ! probability table energies
  real(8), allocatable :: T_tabs(:) ! probability table temperatures
  character(MAX_FILE_LEN) :: path_avg_urr_xs  ! path to averaged URR xs files
  character(MAX_FILE_LEN) :: path_endf        ! path to ENDF-6 files
  character(MAX_FILE_LEN) :: path_prob_tables ! path to probability table files
  character(80), allocatable :: endf_files(:) ! list of ENDF-6 filenames

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
! XSSAMPLE is an object containing a single cross section value for the elastic,
! capture, fission, competitive, and total reactions
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
 
  type xsSample

    real(8) :: t
    real(8) :: n
    real(8) :: g
    real(8) :: f
    real(8) :: x

  end type xsSample
  
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! ISOTOPE is an object containing data for a single isotope with a URR that is
! to be processed
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  type Isotope

    real(8) :: AWR ! weight of nucleus in neutron masses
    real(8) :: T   ! current isotope temperature [K]
    type(ListReal) :: ace_T_list ! list of temperatures isotope has ACE data at
    type(ListInt) :: ace_index_list ! list of indices for different temperatures
    logical :: fissionable ! is isotope fissionable?

    ! Current quantum mechanical variables
    real(8) :: k_n       ! wavenumber for neutron energy
    real(8) :: k_lam     ! wavenumber for resonance energy
    real(8) :: P_l_n     ! penetration for neutron energy
    real(8) :: P_l_lam   ! penetration for resonance energy
    real(8) :: S_l_n     ! shift for neutron energy
    real(8) :: phi_l_n   ! phase shift for neutron energy

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
! TODO: use LRX if ZERO's aren't given for competitive width in ENDF when LRX = 0
    integer :: LFW     ! energy-dependent fission width flag
    integer :: LRX     ! competitive width flag
    integer :: NE      ! number of URR tabulated data energies
    real(8) :: E_ex1 = INF ! first level inelastic scattering excitation energy
    real(8) :: E_ex2 = INF ! second level inelastic scattering excitation energy
    real(8), allocatable :: EL(:)  ! lower energy bound of energy region
    real(8), allocatable :: EH(:)  ! upper energy bound of energy region
    real(8), allocatable :: AP(:)  ! scattering radius
    real(8), allocatable :: ac(:)  ! channel radius
    real(8), allocatable :: SPI(:) ! total spin
    real(8), allocatable :: ES(:)  ! URR tabulated data energies

    ! mean values
    type(VecVecReal), allocatable :: D_mean(:)   ! level spacing
    type(VecVecReal), allocatable :: GN0_mean(:) ! reduced neutron width
    type(VecVecReal), allocatable :: GG_mean(:)  ! radiative capture width
    type(VecVecReal), allocatable :: GF_mean(:)  ! fission width
    type(VecVecReal), allocatable :: GX_mean(:)  ! competitive width
    type(VecReal), allocatable :: AJ(:)   ! total angular momentum
    type(VecReal), allocatable :: DOFN(:) ! # neutron channels
    type(VecReal), allocatable :: DOFG(:) ! # capture channels
    type(VecReal), allocatable :: DOFF(:) ! # fission channels
    type(VecReal), allocatable :: DOFX(:) ! # competitive channels

    ! current values
    real(8) :: E_last = ZERO ! last neutron energy
    real(8) :: E   ! neutron energy
    real(8) :: J   ! total angular momentum
    real(8) :: g_J ! statistical spin factor
    real(8) :: D   ! level spacing
    real(8) :: GN0 ! reduced neutron width
    real(8) :: GG  ! radiative capture width
    real(8) :: GF  ! fission width
    real(8) :: GX  ! competitive width
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
    real(8) :: max_E_urr ! max energy for URR treatment [eV]

    ! probability tables for given (energy, temperature) pairs
    integer :: nE_tabs ! number of probability table energies
    integer :: nT_tabs ! number of probability table temperatures
    integer :: n_bands ! number of probability table bands
    real(8), allocatable :: E_tabs(:) ! probability table energies
    real(8), allocatable :: T_tabs(:) ! probability table temperatures
    type(ProbabilityTable), allocatable :: prob_tables(:,:)
    type(xsSample), allocatable :: xs_samples(:,:)

    ! for each (l, realization), a vector of lists of (l, J) resonances
    type(BWResonanceListVec), allocatable :: urr_resonances_tmp(:,:)

    ! for each (l, realization), a vector of vectors of (l, J) resonances
    type(BWResonanceVecVec), allocatable :: urr_resonances(:,:)

    ! max resonances for a given (l, realization)
    type(VecInt), allocatable :: n_lam(:,:)

    ! vector of Reich-Moore resonances for each l
    type(RMResonanceVec), allocatable :: rm_resonances(:)

    ! vector of BW resonances for each l
    type(BWResonanceVec), allocatable :: bw_resonances(:)

    ! for each (l), a vector of vectors of local (l, J) resonances
    type(BWResonanceVecVec), allocatable :: local_realization(:)

    ! pointwise URR cross section data
    type(ListReal) :: E_tmp ! scratch energy grid values
    type(ListReal) :: n_tmp ! scratch elastic xs
    type(ListReal) :: g_tmp ! scratch capture xs
    type(ListReal) :: f_tmp ! scratch fission xs
    type(ListReal) :: x_tmp ! scratch competitive xs
    type(ListReal) :: t_tmp ! scratch total xs
    real(8), allocatable :: urr_E(:) ! energy grid values
    real(8), allocatable :: urr_n(:) ! elastic xs
    real(8), allocatable :: urr_g(:) ! capture xs
    real(8), allocatable :: urr_f(:) ! fission xs
    real(8), allocatable :: urr_x(:) ! competitive xs
    real(8), allocatable :: urr_t(:) ! total xs

  ! Type-Bound procedures
  contains

    ! allocate resonance energy range variables
    procedure :: alloc_energy_range => alloc_energy_range

    ! deallocate resonance energy range variables
    procedure :: dealloc_energy_ranges => dealloc_energy_ranges

    ! deallocate File 3 average (infinite-dilute) cross sections
    procedure :: dealloc_MF3 => dealloc_MF3

    ! allocate URR resonance ensemble realization
    procedure :: alloc_ensemble => alloc_ensemble

    ! deallocate URR resonance ensemble realization
    procedure :: dealloc_ensemble => dealloc_ensemble

    ! allocate local URR resonance realization
    procedure :: alloc_local_realization => alloc_local_realization

    ! deallocate local URR resonance realization
    procedure :: dealloc_local_realization => dealloc_local_realization

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
  type(xsSample), allocatable :: xs_samples_tmp(:,:)

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! Tabulated Chi-Squared Distribution
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  ! tabulated chi-squared distribution for 1-4 degrees of freedom with expected
  ! bin values taken while preserving equiprobable bins
  real(8), dimension(100,4), parameter :: chi2 = reshape((/&
       5.2361522573878303e-05_8,&
       3.6657014603380927e-04_8,&
       9.9518491948712953e-04_8,&
       1.9386012997931335e-03_8,&
       3.1974134834564794e-03_8,&
       4.7724156249104943e-03_8,&
       6.6646034689871562e-03_8,&
       8.8751764045787822e-03_8,&
       1.1405539948348784e-02_8,&
       1.4257308669262333e-02_8,&
       1.7432309566699937e-02_8,&
       2.0932585917008924e-02_8,&
       2.4760401605548443e-02_8,&
       2.8918245963614714e-02_8,&
       3.3408839132121698e-02_8,&
       3.8235137976558628e-02_8,&
       4.3400342580600923e-02_8,&
       4.8907903348823251e-02_8,&
       5.4761528752270645e-02_8,&
       6.0965193754262480e-02_8,&
       6.7523148957703075e-02_8,&
       7.4439930519469649e-02_8,&
       8.1720370882077772e-02_8,&
       8.9369610378026390e-02_8,&
       9.7393109767770936e-02_8,&
       1.0579666377851940e-01_8,&
       1.1458641571792876e-01_8,&
       1.2376887324425483e-01_8,&
       1.3335092538298998e-01_8,&
       1.4333986088928199e-01_8,&
       1.5374338806581933e-01_8,&
       1.6456965615734245e-01_8,&
       1.7582727845596630e-01_8,&
       1.8752535726568245e-01_8,&
       1.9967351089080942e-01_8,&
       2.1228190283110535e-01_8,&
       2.2536127338660603e-01_8,&
       2.3892297389851302e-01_8,&
       2.5297900387789557e-01_8,&
       2.6754205130376008e-01_8,&
       2.8262553640506494e-01_8,&
       2.9824365927889102e-01_8,&
       3.1441145174029650e-01_8,&
       3.3114483384809063e-01_8,&
       3.4846067560687466e-01_8,&
       3.6637686441014672e-01_8,&
       3.8491237886273233e-01_8,&
       4.0408736970584419e-01_8,&
       4.2392324866638764e-01_8,&
       4.4444278616514943e-01_8,&
       4.6567021895116673e-01_8,&
       4.8763136888124214e-01_8,&
       5.1035377424519790e-01_8,&
       5.3386683524416667e-01_8,&
       5.5820197547795081e-01_8,&
       5.8339282158625361e-01_8,&
       6.0947540353240148e-01_8,&
       6.3648837842569905e-01_8,&
       6.6447328126331806e-01_8,&
       6.9347480655435001e-01_8,&
       7.2354112548443750e-01_8,&
       7.5472424412145056e-01_8,&
       7.8708040918100841e-01_8,&
       8.2067056911047287e-01_8,&
       8.5556089976799887e-01_8,&
       8.9182340583782349e-01_8,&
       9.2953661143024902e-01_8,&
       9.6878635618287412e-01_8,&
       1.0096667167701208e+00_8,&
       1.0522810782469578e+00_8,&
       1.0967433853865443e+00_8,&
       1.1431796114924373e+00_8,&
       1.1917294915911445e+00_8,&
       1.2425485791436743e+00_8,&
       1.2958107014323266e+00_8,&
       1.3517109099512452e+00_8,&
       1.4104690504006043e+00_8,&
       1.4723341150469957e+00_8,&
       1.5375895923398502e+00_8,&
       1.6065601007852524e+00_8,&
       1.6796196952261095e+00_8,&
       1.7572023777114727e+00_8,&
       1.8398155536525733e+00_8,&
       1.9280574813581017e+00_8,&
       2.0226402253972973e+00_8,&
       2.1244203353453672e+00_8,&
       2.2344405927812980e+00_8,&
       2.3539879879111569e+00_8,&
       2.4846761319111303e+00_8,&
       2.6285655985359493e+00_8,&
       2.7883452676527738e+00_8,&
       2.9676159813222007e+00_8,&
       3.1713546653829661e+00_8,&
       3.4067170514840428e+00_8,&
       3.6845270836767279e+00_8,&
       4.0223076143802530e+00_8,&
       4.4512871914763084e+00_8,&
       5.0360207017734400e+00_8,&
       5.9512649086256326e+00_8,&
       8.4491659621504542e+00_8,&
       1.0033501006713585e-02_8,&
       3.0235864759474755e-02_8,&
       5.0644382200351203e-02_8,&
       7.1263304144476777e-02_8,&
       9.2097014254383214e-02_8,&
       1.1315003463415688e-01_8,&
       1.3442703172104936e-01_8,&
       1.5593282249399767e-01_8,&
       1.7767238101947350e-01_8,&
       1.9965084535718647e-01_8,&
       2.2187352484936035e-01_8,&
       2.4434590781962723e-01_8,&
       2.6707366970940932e-01_8,&
       2.9006268168195443e-01_8,&
       3.1331901972664644e-01_8,&
       3.3684897429907912e-01_8,&
       3.6065906053472679e-01_8,&
       3.8475602907845413e-01_8,&
       4.0914687757374879e-01_8,&
       4.3383886286217227e-01_8,&
       4.5883951394450606e-01_8,&
       4.8415664576310247e-01_8,&
       5.0979837386716853e-01_8,&
       5.3577313003122440e-01_8,&
       5.6208967890041917e-01_8,&
       5.8875713574674671e-01_8,&
       6.1578498542415216e-01_8,&
       6.4318310262303713e-01_8,&
       6.7096177353100350e-01_8,&
       6.9913171901966831e-01_8,&
       7.2770411948771163e-01_8,&
       7.5669064150490140e-01_8,&
       7.8610346641513551e-01_8,&
       8.1595532107502389e-01_8,&
       8.4625951092071139e-01_8,&
       8.7702995558135655e-01_8,&
       9.0828122727129657e-01_8,&
       9.4002859223442736e-01_8,&
       9.7228805552877706e-01_8,&
       1.0050764094843230e+00_8,&
       1.0384112861988468e+00_8,&
       1.0723112144860729e+00_8,&
       1.1067956817302635e+00_8,&
       1.1418852011741760e+00_8,&
       1.1776013852112666e+00_8,&
       1.2139670253460237e+00_8,&
       1.2510061795594358e+00_8,&
       1.2887442679197412e+00_8,&
       1.3272081773889526e+00_8,&
       1.3664263769095955e+00_8,&
       1.4064290440029588e+00_8,&
       1.4472482042923114e+00_8,&
       1.4889178855641574e+00_8,&
       1.5314742882274328e+00_8,&
       1.5749559743082053e+00_8,&
       1.6194040774542930e+00_8,&
       1.6648625368154983e+00_8,&
       1.7113783581329287e+00_8,&
       1.7590019059264674e+00_8,&
       1.8077872313378793e+00_8,&
       1.8577924409736915e+00_8,&
       1.9090801130691049e+00_8,&
       1.9617177684434415e+00_8,&
       2.0157784051434957e+00_8,&
       2.0713411073952446e+00_8,&
       2.1284917416161937e+00_8,&
       2.1873237548648219e+00_8,&
       2.2479390943710857e+00_8,&
       2.3104492708727187e+00_8,&
       2.3749765936265970e+00_8,&
       2.4416556114623531e+00_8,&
       2.5106348025719853e+00_8,&
       2.5820785663985397e+00_8,&
       2.6561695848594633e+00_8,&
       2.7331116382691767e+00_8,&
       2.8131329852675315e+00_8,&
       2.8964904480156624e+00_8,&
       2.9834743870013902e+00_8,&
       3.0744148085940641e+00_8,&
       3.1696889297520667e+00_8,&
       3.2697306381410751e+00_8,&
       3.3750424479133900e+00_8,&
       3.4862107856258033e+00_8,&
       3.6039257857375073e+00_8,&
       3.7290072933696399e+00_8,&
       3.8624395681371251e+00_8,&
       4.0054184367488785e+00_8,&
       4.1594166728882609e+00_8,&
       4.3262767786284364e+00_8,&
       4.5083462302929407e+00_8,&
       4.7086809041471618e+00_8,&
       4.9313626468018068e+00_8,&
       5.1820177918730055e+00_8,&
       5.4687119159383473e+00_8,&
       5.8036058655805505e+00_8,&
       6.2063161365943529e+00_8,&
       6.7116592150256986e+00_8,&
       7.3912553622073647e+00_8,&
       8.4377516497364073e+00_8,&
       1.1210340370973686e+01_8,&
       6.8445840470308383e-02_8,&
       1.5100034771915843e-01_8,&
       2.1551987807422432e-01_8,&
       2.7296268032709492e-01_8,&
       3.2622940327935962e-01_8,&
       3.7667790160175518e-01_8,&
       4.2508370358766145e-01_8,&
       4.7194279239893744e-01_8,&
       5.1759698373925178e-01_8,&
       5.6229485133305723e-01_8,&
       6.0622470846289511e-01_8,&
       6.4953392303497304e-01_8,&
       6.9234093087234361e-01_8,&
       7.3474307291383112e-01_8,&
       7.7682191088101582e-01_8,&
       8.1864695147378741e-01_8,&
       8.6027832833757922e-01_8,&
       9.0176877983612724e-01_8,&
       9.4316513811993008e-01_8,&
       9.8450947108920839e-01_8,&
       1.0258399727912586e+00_8,&
       1.0671916682347822e+00_8,&
       1.1085969791461903e+00_8,&
       1.1500861840892049e+00_8,&
       1.1916877973642679e+00_8,&
       1.2334288848024455e+00_8,&
       1.2753353300859189e+00_8,&
       1.3174320619909312e+00_8,&
       1.3597432505811178e+00_8,&
       1.4022924786272399e+00_8,&
       1.4451028932156729e+00_8,&
       1.4881973415202965e+00_8,&
       1.5315984939545633e+00_8,&
       1.5753289573465672e+00_8,&
       1.6194113803319028e+00_8,&
       1.6638685528185568e+00_8,&
       1.7087235011125861e+00_8,&
       1.7539995800922088e+00_8,&
       1.7997205636636548e+00_8,&
       1.8459107346190737e+00_8,&
       1.8925949749385853e+00_8,&
       1.9397988575227763e+00_8,&
       1.9875487403178107e+00_8,&
       2.0358718637868938e+00_8,&
       2.0847964526978844e+00_8,&
       2.1343518232300300e+00_8,&
       2.1845684964574215e+00_8,&
       2.2354783193401633e+00_8,&
       2.2871145944518232e+00_8,&
       2.3395122197887064e+00_8,&
       2.3927078401548001e+00_8,&
       2.4467400117919293e+00_8,&
       2.5016493821335106e+00_8,&
       2.5574788868147769e+00_8,&
       2.6142739663679690e+00_8,&
       2.6720828053882273e+00_8,&
       2.7309565973775838e+00_8,&
       2.7909498389794463e+00_8,&
       2.8521206579161578e+00_8,&
       2.9145311796652655e+00_8,&
       2.9782479387776153e+00_8,&
       3.0433423417829970e+00_8,&
       3.1098911899031334e+00_8,&
       3.1779772713287335e+00_8,&
       3.2476900347077642e+00_8,&
       3.3191263578080776e+00_8,&
       3.3923914281841236e+00_8,&
       3.4675997562325906e+00_8,&
       3.5448763454753069e+00_8,&
       3.6243580505004620e+00_8,&
       3.7061951600825447e+00_8,&
       3.7905532520487486e+00_8,&
       3.8776153780842386e+00_8,&
       3.9675846517429401e+00_8,&
       4.0606873326438402e+00_8,&
       4.1571765258515283e+00_8,&
       4.2573366501392140e+00_8,&
       4.3614888756036256e+00_8,&
       4.4699977948979086e+00_8,&
       4.5832796804675722e+00_8,&
       4.7018128035906628e+00_8,&
       4.8261504664869381e+00_8,&
       4.9569376523570625e+00_8,&
       5.0949325714431399e+00_8,&
       5.2410349415867721e+00_8,&
       5.3963237020682318e+00_8,&
       5.5621082135126638e+00_8,&
       5.7399991881525931e+00_8,&
       5.9320092563219085e+00_8,&
       6.1406994197980982e+00_8,&
       6.3693991105142747e+00_8,&
       6.6225493533220510e+00_8,&
       6.9062624129574450e+00_8,&
       7.2292862961073432e+00_8,&
       7.6047873656073781e+00_8,&
       8.0539619114173657e+00_8,&
       8.6143443538629256e+00_8,&
       9.3629335234829068e+00_8,&
       1.0506089294536801e+01_8,&
       1.3486550432895573e+01_8,&
       1.9559695720196493e-01_8,&
       3.6654661317972215e-01_8,&
       4.8369150691127061e-01_8,&
       5.8196287845694095e-01_8,&
       6.6951354752681236e-01_8,&
       7.4996729294063158e-01_8,&
       8.2531974682625786e-01_8,&
       8.9681085396225890e-01_8,&
       9.6527447379825393e-01_8,&
       1.0313039645081317e+00_8,&
       1.0953399453660284e+00_8,&
       1.1577207743428419e+00_8,&
       1.2187134520889087e+00_8,&
       1.2785335015626327e+00_8,&
       1.3373582809681861e+00_8,&
       1.3953362022539313e+00_8,&
       1.4525932975404519e+00_8,&
       1.5092380114338664e+00_8,&
       1.5653647732450024e+00_8,&
       1.6210567097541417e+00_8,&
       1.6763877397217370e+00_8,&
       1.7314242153557968e+00_8,&
       1.7862262263145512e+00_8,&
       1.8408486486482984e+00_8,&
       1.8953419984428517e+00_8,&
       1.9497531341834025e+00_8,&
       2.0041258407262728e+00_8,&
       2.0585013197703494e+00_8,&
       2.1129186059035727e+00_8,&
       2.1674149230066977e+00_8,&
       2.2220259925977954e+00_8,&
       2.2767863032911895e+00_8,&
       2.3317293487084729e+00_8,&
       2.3868878397761839e+00_8,&
       2.4422938962570533e+00_8,&
       2.4979792215205392e+00_8,&
       2.5539752639006177e+00_8,&
       2.6103133674763810e+00_8,&
       2.6670249147094554e+00_8,&
       2.7241414630623040e+00_8,&
       2.7816948774816499e+00_8,&
       2.8397174604512325e+00_8,&
       2.8982420811839344e+00_8,&
       2.9573023054391534e+00_8,&
       3.0169325273909999e+00_8,&
       3.0771681049573041e+00_8,&
       3.1380455000074954e+00_8,&
       3.1996024249084458e+00_8,&
       3.2618779969343259e+00_8,&
       3.3249129021703441e+00_8,&
       3.3887495706703041e+00_8,&
       3.4534323647957055e+00_8,&
       3.5190077828775581e+00_8,&
       3.5855246805903294e+00_8,&
       3.6530345127414847e+00_8,&
       3.7215915985445833e+00_8,&
       3.7912534138916278e+00_8,&
       3.8620809146701407e+00_8,&
       3.9341388958097676e+00_8,&
       4.0074963915081963e+00_8,&
       4.0822271230123546e+00_8,&
       4.1584100014434586e+00_8,&
       4.2361296945114582e+00_8,&
       4.3154772676087871e+00_8,&
       4.3965509117920645e+00_8,&
       4.4794567736368460e+00_8,&
       4.5643099050101883e+00_8,&
       4.6512353546132843e+00_8,&
       4.7403694278956019e+00_8,&
       4.8318611479259577e+00_8,&
       4.9258739573759858e+00_8,&
       5.0225877114340918e+00_8,&
       5.1222010238797715e+00_8,&
       5.2249340446443524e+00_8,&
       5.3310317682151434e+00_8,&
       5.4407680000031222e+00_8,&
       5.5544501448001302e+00_8,&
       5.6724250313187019e+00_8,&
       5.7950860547996710e+00_8,&
       5.9228820135617131e+00_8,&
       6.0563281468068970e+00_8,&
       6.1960200678016255e+00_8,&
       6.3426515564370485e+00_8,&
       6.4970375722046603e+00_8,&
       6.6601444445014435e+00_8,&
       6.8331301115602772e+00_8,&
       7.0173987177487653e+00_8,&
       7.2146762060187619e+00_8,&
       7.4271174283417203e+00_8,&
       7.6574620270035592e+00_8,&
       7.9092684969487443e+00_8,&
       8.1872789146451925e+00_8,&
       8.4980132789380303e+00_8,&
       8.8507929093144444e+00_8,&
       9.2596300868203301e+00_8,&
       9.7470523027115679e+00_8,&
       1.0352888667619078e+01_8,&
       1.1158690088753749e+01_8,&
       1.2382461545773161e+01_8,&
       1.5538540717378252e+01_8/),(/100,4/))

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
    integer :: i_res       ! l-wave resonance counter
    integer :: n_res       ! number of l-wave resonances to include
    integer :: n_above_urr ! number of resonances abover upper URR energy
    real(8) :: E_res ! current resonance (lab) energy (e.g. E_lam)
    real(8) :: m     ! energy interpolation factor

    tope => isotopes(iso)

    ! allocate vector of linked lists of (l,J) resonances for (l, realization)
    allocate(tope % urr_resonances_tmp(tope % NLS(tope % i_urr), n_reals))
    allocate(tope % n_lam(tope % NLS(tope % i_urr), n_reals))

    res % i_res = 0

    ! loop over independent realizations
    ENSEMBLE_LOOP: do i_ens = 1, n_reals

      ! loop over orbital angular momenta
      ORBITAL_ANG_MOM_LOOP: do i_l = 1, tope % NLS(tope % i_urr)

        ! alloc vector of NJS(l) linked lists of resonances for (l, realization)
        allocate(tope % urr_resonances_tmp(i_l, i_ens) % J(tope % NJS(i_l)))
        allocate(tope % n_lam(i_l, i_ens) % dim1(tope % NJS(i_l)))

        ! set current orbital angular momentum quantum number
        tope % L = i_l - 1

        ! get the number of contributing l-wave resonances for this l
        n_res = n_res_contrib(tope % L)

        ! loop over total angular momenta
        TOTAL_ANG_MOM_LOOP: do i_J = 1, tope % NJS(i_l)

          ! set current total angular momentum quantum number
          tope % J = tope % AJ(i_l) % dim1(i_J)

          ! set current partial width degrees of freedom
          tope % AMUX = int(tope % DOFX(i_l) % dim1(i_J))
          tope % AMUN = int(tope % DOFN(i_l) % dim1(i_J))
          tope % AMUG = int(tope % DOFG(i_l) % dim1(i_J))
          tope % AMUF = int(tope % DOFF(i_l) % dim1(i_J))

          ! set energy of the lowest-lying contributing URR resonance
          tope % D = tope % D_mean(i_l) % dim2(i_J) % dim1(1)
          if (i_l > tope % NLS(tope % i_urr - 1)) then
            ! the URR has more l-states than the RRR; place resonance energy
            ! randomly about lower URR energy bound
            E_res = tope % EL(tope % i_urr)&
                 + (ONE - TWO * prn()) * wigner_surmise(tope % D)
          else
            ! offset first URR resonance energy from the highest-energy RRR
            ! resonance with the same (l,J) spin sequence
            E_res = E_last_rrr(iso, tope % L, tope % J)+wigner_surmise(tope % D)
          end if

          ! point to first resonance for this l-wave, realization, J
          tope % urr_resonances_tmp(i_l, i_ens) % J(i_J) % res&
               => tope % urr_resonances_tmp(i_l, i_ens) % J(i_J) % first

          ! zero resonance counters
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
                 tope % INT)

            ! set current mean unresolved resonance parameters
            tope % D = interpolator(m, &
                 tope % D_mean(i_l) % dim2(i_J) % dim1(i_E), &
                 tope % D_mean(i_l) % dim2(i_J) % dim1(i_E + 1), tope % INT)
            tope % GN0 = interpolator(m, &
                 tope % GN0_mean(i_l) % dim2(i_J) % dim1(i_E), &
                 tope % GN0_mean(i_l) % dim2(i_J) % dim1(i_E + 1), tope % INT)
            tope % GG = interpolator(m, &
                 tope % GG_mean(i_l) % dim2(i_J) % dim1(i_E), &
                 tope % GG_mean(i_l) % dim2(i_J) % dim1(i_E + 1), tope % INT)
            if (tope % INT == LINEAR_LINEAR &
                 .or. tope % GF_mean(i_l) % dim2(i_J) % dim1(i_E) > ZERO) then
              tope % GF = interpolator(m, &
                   tope % GF_mean(i_l) % dim2(i_J) % dim1(i_E), &
                   tope % GF_mean(i_l) % dim2(i_J) % dim1(i_E + 1), tope % INT)
            else
              tope % GF = ZERO
            end if
            if (tope % INT == LINEAR_LINEAR &
              .or. tope % GX_mean(i_l) % dim2(i_J) % dim1(i_E) > ZERO) then
              tope % GX = interpolator(m, &
                   tope % GX_mean(i_l) % dim2(i_J) % dim1(i_E), &
                   tope % GX_mean(i_l) % dim2(i_J) % dim1(i_E + 1), tope % INT)
            else
              tope % GX = ZERO
            end if

            ! sample unresolved resonance parameters for this spin
            ! sequence, at this energy
            res % E_lam = E_res
            tope % E = E_res
            tope % k_n   = wavenumber(tope % AWR, abs(tope % E))
            tope % k_lam = wavenumber(tope % AWR, abs(tope % E))
            tope % P_l_n = penetration(tope % L,&
                 tope % k_n * tope % ac(tope % i_urr))
            tope % P_l_lam = penetration(tope % L,&
                 tope % k_lam * tope % ac(tope % i_urr))
            tope % S_l_n = shift(tope % L,&
                 tope % k_n * tope % ac(tope % i_urr))
            tope % phi_l_n = phase_shift(tope % L,&
                 tope % k_n * tope % AP(tope % i_urr))
            call res % channel_width(iso, i_l, i_J)
            call add_parameters(res, iso, i_ens, i_l, i_J)

            ! add an additional resonance
            E_res = E_res + wigner_surmise(tope % D)

            if (E_res > tope % EH(tope % i_urr)) n_above_urr = n_above_urr + 1

          end do RESONANCE_LOOP
  
          tope % n_lam(i_l, i_ens) % dim1(i_J) = i_res

        end do TOTAL_ANG_MOM_LOOP
      end do ORBITAL_ANG_MOM_LOOP
    end do ENSEMBLE_LOOP

    ! allocate URR resonance ensemble realizations
    call tope % alloc_ensemble()

    ! transfer linked list URR ensembles to permanent array
    do i_ens = 1, n_reals
      do i_l = 1, tope % NLS(tope % i_urr)
        do i_J = 1, tope % NJS(i_l)

          ! point to first resonance for this l-wave, realization, J
          tope % urr_resonances_tmp(i_l, i_ens) % J(i_J) % res&
               => tope % urr_resonances_tmp(i_l, i_ens) % J(i_J) % first
        
          do i_res = 1, tope % n_lam(i_l, i_ens) % dim1(i_J)

            ! resonance energy
            tope % urr_resonances(i_l, i_ens) % J(i_J) % res(i_res) % E_lam&
                 = tope % urr_resonances_tmp(i_l, i_ens) % J(i_J) % res % E_lam
          
            ! total angular momentum, J
            tope % urr_resonances(i_l, i_ens) % J(i_J) % res(i_res) % AJ&
                 = tope % urr_resonances_tmp(i_l, i_ens) % J(i_J) % res % AJ
          
            ! neutron width
            tope % urr_resonances(i_l, i_ens) % J(i_J) % res(i_res) % GN&
                 = tope % urr_resonances_tmp(i_l, i_ens) % J(i_J) % res % GN
          
            ! gamma width
            tope % urr_resonances(i_l, i_ens) % J(i_J) % res(i_res) % GG&
                 = tope % urr_resonances_tmp(i_l, i_ens) % J(i_J) % res % GG
          
            ! fission width
            tope % urr_resonances(i_l, i_ens) % J(i_J) % res(i_res) % GF&
                 = tope % urr_resonances_tmp(i_l, i_ens) % J(i_J) % res % GF
          
            ! competitive width
            tope % urr_resonances(i_l, i_ens) % J(i_J) % res(i_res) % GX&
                 = tope % urr_resonances_tmp(i_l, i_ens) % J(i_J) % res % GX
          
            ! total width
            tope % urr_resonances(i_l, i_ens) % J(i_J) % res(i_res) % GT&
                 = tope % urr_resonances_tmp(i_l, i_ens) % J(i_J) % res % GT
       
            ! point to next resonance
            tope % urr_resonances_tmp(i_l, i_ens) % J(i_J) % res&
                 => tope % urr_resonances_tmp(i_l, i_ens) % J(i_J) % res % next

          end do

          call tope % urr_resonances_tmp(i_l, i_ens) % J(i_J) % clear()

        end do

        deallocate(tope % urr_resonances_tmp(i_l, i_ens) % J)

      end do
    end do

    ! deallocate temporary linked list URR resonance ensemble realizations
    deallocate(tope % urr_resonances_tmp)

  end subroutine resonance_ensemble

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! POINTWISE_URR generates pointwise energy-cross section data in the URR
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine pointwise_urr(iso, T_K)

    type(Isotope), pointer :: tope ! isotope object pointer
    type(Resonance) :: res ! resonance object
    type(CrossSection) :: n ! elastic scattering xs object
    type(CrossSection) :: g ! radiative capture xs object
    type(CrossSection) :: f ! fission xs object
    type(CrossSection) :: x ! competitive xs object
    type(CrossSection) :: t ! total xs object
    integer :: iso       ! isotope index
    integer :: i_l       ! orbital quantum number
    integer :: i_J       ! total angular momentum quantum number
    integer :: n_res     ! number of contributing l-state resonances
    integer :: n_rrr_res ! number of RRR resonances we need to grab
    integer :: i_low     ! index of lowest-lying resonance
    integer :: i_res     ! resonance counter
    integer :: i_rrr_res ! RRR resonance index
    integer :: i_ES      ! index of current URR tabulated energy
    integer :: iE        ! pointwise energy grid index
    real(8) :: T_K          ! isotope temperature [K]
    real(8) :: E_lo         ! lower energy grid point for interpolation
    real(8) :: E_hi         ! upper energy grid point for interpolation
    real(8) :: E_mid        ! trial interpolated energy grid point
    real(8) :: xs_lo        ! total xs at E_lo
    real(8) :: xs_hi        ! total xs at E_hi
    real(8) :: xs_mid       ! total xs value at E_mid
    real(8) :: xs_mid_trial ! interpolated xs at E_mid
    real(8) :: rel_err      ! relative error between interpolated and exact xs
    real(8) :: E_0          ! energy at current list index
    real(8) :: E_1          ! energy at next list index
    real(8) :: E_last       ! energy at last grid index
    integer :: i_list       ! index in list of resonances for all (l,J)
    integer :: i_list_start ! index where to start searching for E_lam

    tope => isotopes(iso)
    tope % T = T_K

    ! only one realization allowed when using a pointwise representation
    i_real = 1

    ! initialize linked list of energies to lower URR bound
    call tope % E_tmp % append(tope % EL(tope % i_urr))

    ! loop over resonances for all spin sequences and add energies
    do i_l = 1, tope % NLS(tope % i_urr)
      do i_J = 1, tope % NJS(i_l)
        if (master)&
             write(*,'(I7,A48,I1,A1,F4.1)') tope % ZAI,&
             ': Generating resonances for (l,J) spin sequence ',&
             i_l - 1, ',', tope % AJ(i_l) % dim1(i_J)
        i_list_start = 1
        do i_res = 1, tope % n_lam(i_l, i_real) % dim1(i_J)
          do i_list = i_list_start, tope % E_tmp % size()
            if (tope%urr_resonances(i_l,i_real)%J(i_J)%res(i_res)%E_lam <=&
                 tope % E_tmp % get_item(i_list)) exit
          end do
          call tope % E_tmp % insert(i_list,&
               tope % urr_resonances(i_l,i_real) % J(i_J) % res(i_res) % E_lam)
          i_list_start = i_list
        end do
      end do
    end do

    ! clean energy grid of duplicates, values outside URR
    i_list = 1
    E_last = min(tope % EH(tope % i_urr), tope % max_E_urr)
    do
      if (i_list == tope % E_tmp % size()) exit
      E_0 = tope % E_tmp % get_item(i_list)
      E_1 = tope % E_tmp % get_item(i_list+1)
      if (E_0 > E_1)&
           call fatal_error('Pointwise URR energy grid not monotonically&
             & increasing')
      if ((E_0 == E_1)&
           .or. (E_0 > E_last)&
           .or. (E_0 < tope % EL(tope % i_urr))) then
        call tope % E_tmp % remove(E_0)
        cycle
      end if
      i_list = i_list + 1
    end do
    E_0 = tope % E_tmp % get_item(i_list)
    if ((E_0 >= E_last) .or. (E_0 < tope % EL(tope % i_urr)))&
         call tope % E_tmp % remove(E_0)
    call tope % E_tmp % append(E_last)

    ! set first two energy points
    iE = 1
    E_lo = tope % E_tmp % get_item(iE)
    E_hi = tope % E_tmp % get_item(iE+1)

    ! calculate xs vals at the first energy
    tope % E = E_lo
    tope % k_n = wavenumber(tope % AWR, abs(tope % E))

    ! reset xs accumulators
    call flush_sigmas(t, n, g, f, x)

    ! Get resonance parameters for a local realization about E_n
    ! loop over orbital quantum numbers
    LOC_ORBITAL_ANG_MOM_LOOP_1: do i_l = 1, tope % NLS(tope % i_urr)

      ! set current orbital angular momentum quantum number
      tope % L = i_l - 1

      ! get the number of contributing l-wave resonances for this l
      n_res = n_res_contrib(tope % L)

      ! loop over total quantum numbers
      LOC_TOTAL_ANG_MOM_LOOP_1: do i_J = 1, tope % NJS(i_l)

        if (tope % E&
             < tope % urr_resonances(i_l, i_real) % J(i_J) % res(1)%E_lam) then
          i_low = 1
        else
          i_low = binary_search(&
               tope % urr_resonances(i_l,i_real) % J(i_J) % res(:) % E_lam,&
               tope % n_lam(i_l, i_real) % dim1(i_J), tope % E)
        end if

        ! set current total angular momentum quantum number
        tope % J = tope % AJ(i_l) % dim1(i_J)

        res % i_res = 0
        ! loop over the addition of resonances to this ladder
        if (i_low - n_res/2 + 1 < 1) then
          ! if we're near the lower end of the URR, need to incorporate
          ! resolved resonance region resonances in order to fix-up
          ! (i.e. smooth out) cross sections at the RRR-URR crossover
          ! energy

          ! if the RRR has resonances with this l-state
          if (i_l <= tope % NLS(tope % i_urr - 1)) then

            ! how many RRR resonances are contributing
            n_rrr_res = abs(i_low - n_res/2)

            ! loop over contributing resolved resonance region resonances
            LOC_RRR_RESONANCES_LOOP_1: do i_res = n_rrr_res, 1, -1
              i_rrr_res = rrr_res(iso, i_res, tope % L, tope % J)

              ! fewer RRR resonances w/ this J value then needed;
              ! just generate URR resonances instead
!TODO: take however many RRR resonances there actually are, even if too few
              if (i_rrr_res == 0) exit

              ! add this resolved resonance
              res % i_res = res % i_res + 1
              call set_parameters(res, iso, i_rrr_res, i_l, i_J,&
                   tope % i_urr - 1)

            end do LOC_RRR_RESONANCES_LOOP_1
          end if

          ! loop over contributing unresolved resonance region resonances
          LOC_URR_RESONANCES_LOOP_1: do i_res = 1, n_res - res % i_res
            res % i_res = res % i_res + 1
            call set_parameters(res, iso, i_res, i_l, i_J, tope % i_urr)
          end do LOC_URR_RESONANCES_LOOP_1

        else
          ! we're firmly in the URR and can ignore anything going on in
          ! the upper resolved resonance region energies
          LOC_URR_LOOP_1: do i_res = i_low - n_res/2 + 1, i_low + n_res/2
            res % i_res = res % i_res + 1
            call set_parameters(res, iso, i_res, i_l, i_J, tope % i_urr)
          end do LOC_URR_LOOP_1
        end if
      end do LOC_TOTAL_ANG_MOM_LOOP_1
    end do LOC_ORBITAL_ANG_MOM_LOOP_1

    ! loop over orbital quantum numbers
    ORBITAL_ANG_MOM_LOOP_1: do i_l = 1, tope % NLS(tope % i_urr)

      ! set current orbital angular momentum quantum number
      tope % L = i_l - 1

      ! penetration
      tope % P_l_n = penetration(tope % L, tope % k_n * tope % ac(tope % i_urr))

      ! resonance energy shift factor
      tope % S_l_n = shift(tope % L, tope % k_n * tope % ac(tope % i_urr))

      ! hard-sphere phase shift
      tope % phi_l_n = phase_shift(tope % L, tope % k_n * tope % AP(tope%i_urr))

      ! get the number of contributing l-wave resonances for this l
      n_res = n_res_contrib(tope % L)

      ! loop over total angular momentum quantum numbers
      TOTAL_ANG_MOM_LOOP_1: do i_J = 1, tope % NJS(i_l)

        ! set current total angular momentum quantum number
        tope % J = tope % AJ(i_l) % dim1(i_J)

        ! compute statistical spin factor
        tope % g_J = (TWO*tope % J + ONE) / (FOUR*tope % SPI(tope%i_urr) + TWO)

        ! loop over resonances localized about E_n
        RESONANCES_LOOP_1: do i_res = 1, n_res

          res % i_res = i_res
          res % E_lam = tope%local_realization(i_l) % J(i_J) % res(i_res)%E_lam
          res % Gam_n = tope%local_realization(i_l) % J(i_J) % res(i_res)%GN
          res % Gam_g = tope%local_realization(i_l) % J(i_J) % res(i_res)%GG
          res % Gam_f = tope%local_realization(i_l) % J(i_J) % res(i_res)%GF
          res % Gam_x = tope%local_realization(i_l) % J(i_J) % res(i_res)%GX
          res % Gam_t = tope%local_realization(i_l) % J(i_J) % res(i_res)%GT

          call res % calc_xs(iso)
          call accum_resonances(res, t, n, g, f, x)

        end do RESONANCES_LOOP_1
      end do TOTAL_ANG_MOM_LOOP_1
    end do ORBITAL_ANG_MOM_LOOP_1

    ! add potential scattering contribution
    call n % potential_xs(iso)
    call t % potential_xs(iso)

    ! combine resonance and File 3 components
    call mf3_self_shielding(iso, n, g, f, x, t)

    ! initialize xs linked lists
    call tope % n_tmp % append(n % xs)
    call tope % g_tmp % append(g % xs)
    call tope % f_tmp % append(f % xs)
    call tope % x_tmp % append(x % xs)
    call tope % t_tmp % append(t % xs)

    xs_lo = tope % t_tmp % get_item(iE)

    ! calculate xs vals at the second energy
    tope % E = E_hi
    tope % k_n = wavenumber(tope % AWR, abs(tope % E))

    ! reset xs accumulators
    call flush_sigmas(t, n, g, f, x)

    ! loop over orbital quantum numbers
    ORBITAL_ANG_MOM_LOOP_2: do i_l = 1, tope % NLS(tope % i_urr)

      ! set current orbital angular momentum quantum number
      tope % L = i_l - 1

      ! penetration
      tope % P_l_n = penetration(tope % L, tope % k_n * tope % ac(tope % i_urr))

      ! resonance energy shift factor
      tope % S_l_n = shift(tope % L, tope % k_n * tope % ac(tope % i_urr))

      ! hard-sphere phase shift
      tope % phi_l_n = phase_shift(tope % L, tope % k_n * tope % AP(tope%i_urr))

      ! get the number of contributing l-wave resonances for this l
      n_res = n_res_contrib(tope % L)

      ! loop over total angular momentum quantum numbers
      TOTAL_ANG_MOM_LOOP_2: do i_J = 1, tope % NJS(i_l)

        ! set current total angular momentum quantum number
        tope % J = tope % AJ(i_l) % dim1(i_J)

        ! compute statistical spin factor
        tope % g_J = (TWO*tope % J + ONE) / (FOUR*tope % SPI(tope%i_urr) + TWO)

        ! loop over resonances localized about E_n
        RESONANCES_LOOP_2: do i_res = 1, n_res

          res % i_res = i_res
          res % E_lam = tope%local_realization(i_l) % J(i_J) % res(i_res)%E_lam
          res % Gam_n = tope%local_realization(i_l) % J(i_J) % res(i_res)%GN
          res % Gam_g = tope%local_realization(i_l) % J(i_J) % res(i_res)%GG
          res % Gam_f = tope%local_realization(i_l) % J(i_J) % res(i_res)%GF
          res % Gam_x = tope%local_realization(i_l) % J(i_J) % res(i_res)%GX
          res % Gam_t = tope%local_realization(i_l) % J(i_J) % res(i_res)%GT

          call res % calc_xs(iso)
          call accum_resonances(res, t, n, g, f, x)

        end do RESONANCES_LOOP_2
      end do TOTAL_ANG_MOM_LOOP_2
    end do ORBITAL_ANG_MOM_LOOP_2

    ! add potential scattering contribution
    call n % potential_xs(iso)
    call t % potential_xs(iso)

    ! combine resonance and File 3 components
    call mf3_self_shielding(iso, n, g, f, x, t)

    xs_hi = t % xs

    ! refine energy grid until total cross section convergence
    i_ES = 2
    do
      if (i_ES <= tope % NE) then
        if (tope % E > tope % ES(i_ES - 1)) then
          if (master) write(*,'(I7,A29,ES12.5,A3)') tope % ZAI,&
               ': Reconstructing URR xs below', tope % ES(i_ES), ' eV'
          i_ES = i_ES + 1
        end if
      end if

      ! compute interpolated energy and xs value
      E_mid = (E_lo + E_hi) / TWO
      xs_mid_trial = (xs_lo + xs_hi) / TWO

      ! calculate xs vals at the interpolated energy
      tope % E = E_mid
      tope % k_n = wavenumber(tope % AWR, abs(tope % E))

      ! reset xs accumulators
      call flush_sigmas(t, n, g, f, x)

      ! loop over orbital quantum numbers
      ORBITAL_ANG_MOM_LOOP_3: do i_l = 1, tope % NLS(tope % i_urr)

        ! set current orbital angular momentum quantum number
        tope % L = i_l - 1

        ! penetration
        tope % P_l_n = penetration(tope % L, tope % k_n * tope % ac(tope % i_urr))

        ! resonance energy shift factor
        tope % S_l_n = shift(tope % L, tope % k_n * tope % ac(tope % i_urr))

        ! hard-sphere phase shift
        tope % phi_l_n = phase_shift(tope % L, tope % k_n * tope % AP(tope%i_urr))

        ! get the number of contributing l-wave resonances for this l
        n_res = n_res_contrib(tope % L)

        ! loop over total angular momentum quantum numbers
        TOTAL_ANG_MOM_LOOP_3: do i_J = 1, tope % NJS(i_l)

          ! set current total angular momentum quantum number
          tope % J = tope % AJ(i_l) % dim1(i_J)

          ! compute statistical spin factor
          tope % g_J = (TWO*tope % J + ONE) / (FOUR*tope % SPI(tope%i_urr) + TWO)

          ! loop over resonances localized about E_n
          RESONANCES_LOOP_3: do i_res = 1, n_res

            res % i_res = i_res
            res % E_lam = tope%local_realization(i_l) % J(i_J) % res(i_res)%E_lam
            res % Gam_n = tope%local_realization(i_l) % J(i_J) % res(i_res)%GN
            res % Gam_g = tope%local_realization(i_l) % J(i_J) % res(i_res)%GG
            res % Gam_f = tope%local_realization(i_l) % J(i_J) % res(i_res)%GF
            res % Gam_x = tope%local_realization(i_l) % J(i_J) % res(i_res)%GX
            res % Gam_t = tope%local_realization(i_l) % J(i_J) % res(i_res)%GT

            call res % calc_xs(iso)
            call accum_resonances(res, t, n, g, f, x)

          end do RESONANCES_LOOP_3
        end do TOTAL_ANG_MOM_LOOP_3
      end do ORBITAL_ANG_MOM_LOOP_3

      ! add potential scattering contribution
      call n % potential_xs(iso)
      call t % potential_xs(iso)

      ! combine resonance and File 3 components
      call mf3_self_shielding(iso, n, g, f, x, t)

      xs_mid = t % xs

      ! compute relative error
      if (xs_mid <= ZERO) then
        rel_err = INF
      else
        rel_err = abs((xs_mid_trial - xs_mid) / xs_mid)
      end if

      ! refine energy mesh or accept the point
      if (rel_err <= tol_point_urr .or. (E_mid - E_lo <= min_dE_point_urr)) then
        iE = iE + 1
        if (.not. tope % E_tmp % contains(E_mid))&
             call tope % E_tmp % insert(iE, E_mid)
        call tope % n_tmp % insert(iE, n % xs)
        call tope % g_tmp % insert(iE, g % xs)
        call tope % f_tmp % insert(iE, f % xs)
        call tope % x_tmp % insert(iE, x % xs)
        call tope % t_tmp % insert(iE, t % xs)

        ! proceed from the midpoint energy just added
        E_lo = E_mid
        E_hi = tope % E_tmp % get_item(iE+1)
        xs_lo = xs_mid

        ! calculate xs vals at the new upper energy
        tope % E = E_hi
        tope % k_n = wavenumber(tope % AWR, abs(tope % E))

        ! reset xs accumulators
        call flush_sigmas(t, n, g, f, x)

        ! Get resonance parameters for a local realization about E_n
        ! loop over orbital quantum numbers
        LOC_ORBITAL_ANG_MOM_LOOP_4: do i_l = 1, tope % NLS(tope % i_urr)

          ! set current orbital angular momentum quantum number
          tope % L = i_l - 1

          ! get the number of contributing l-wave resonances for this l
          n_res = n_res_contrib(tope % L)

          ! loop over total quantum numbers
          LOC_TOTAL_ANG_MOM_LOOP_4: do i_J = 1, tope % NJS(i_l)

            if (tope % E&
                 < tope % urr_resonances(i_l, i_real) % J(i_J) % res(1)%E_lam) then
              i_low = 1
            else
              i_low = binary_search(&
                   tope % urr_resonances(i_l,i_real) % J(i_J) % res(:) % E_lam,&
                   tope % n_lam(i_l, i_real) % dim1(i_J), tope % E)
            end if

            ! set current total angular momentum quantum number
            tope % J = tope % AJ(i_l) % dim1(i_J)

            res % i_res = 0
            ! loop over the addition of resonances to this ladder
            if (i_low - n_res/2 + 1 < 1) then
              ! if we're near the lower end of the URR, need to incorporate
              ! resolved resonance region resonances in order to fix-up
              ! (i.e. smooth out) cross sections at the RRR-URR crossover
              ! energy

              ! if the RRR has resonances with this l-state
              if (i_l <= tope % NLS(tope % i_urr - 1)) then

                ! how many RRR resonances are contributing
                n_rrr_res = abs(i_low - n_res/2)

                ! loop over contributing resolved resonance region resonances
                LOC_RRR_RESONANCES_LOOP_4: do i_res = n_rrr_res, 1, -1
                  i_rrr_res = rrr_res(iso, i_res, tope % L, tope % J)

                  ! fewer RRR resonances w/ this J value then needed;
                  ! just generate URR resonances instead
    !TODO: take however many RRR resonances there actually are, even if too few
                  if (i_rrr_res == 0) exit

                  ! add this resolved resonance
                  res % i_res = res % i_res + 1
                  call set_parameters(res, iso, i_rrr_res, i_l, i_J,&
                       tope % i_urr - 1)

                end do LOC_RRR_RESONANCES_LOOP_4
              end if

              ! loop over contributing unresolved resonance region resonances
              LOC_URR_RESONANCES_LOOP_4: do i_res = 1, n_res - res % i_res
                res % i_res = res % i_res + 1
                call set_parameters(res, iso, i_res, i_l, i_J, tope % i_urr)
              end do LOC_URR_RESONANCES_LOOP_4

            else
              ! we're firmly in the URR and can ignore anything going on in
              ! the upper resolved resonance region energies
              LOC_URR_LOOP_4: do i_res = i_low - n_res/2 + 1, i_low + n_res/2
                res % i_res = res % i_res + 1
                call set_parameters(res, iso, i_res, i_l, i_J, tope % i_urr)
              end do LOC_URR_LOOP_4
            end if
          end do LOC_TOTAL_ANG_MOM_LOOP_4
        end do LOC_ORBITAL_ANG_MOM_LOOP_4

        ! loop over orbital quantum numbers
        ORBITAL_ANG_MOM_LOOP_4: do i_l = 1, tope % NLS(tope % i_urr)

          ! set current orbital angular momentum quantum number
          tope % L = i_l - 1

          ! penetration
          tope % P_l_n = penetration(tope % L, tope % k_n * tope % ac(tope % i_urr))

          ! resonance energy shift factor
          tope % S_l_n = shift(tope % L, tope % k_n * tope % ac(tope % i_urr))

          ! hard-sphere phase shift
          tope % phi_l_n = phase_shift(tope % L, tope % k_n * tope % AP(tope%i_urr))

          ! get the number of contributing l-wave resonances for this l
          n_res = n_res_contrib(tope % L)

          ! loop over total angular momentum quantum numbers
          TOTAL_ANG_MOM_LOOP_4: do i_J = 1, tope % NJS(i_l)

            ! set current total angular momentum quantum number
            tope % J = tope % AJ(i_l) % dim1(i_J)

            ! compute statistical spin factor
            tope % g_J = (TWO*tope % J + ONE) / (FOUR*tope % SPI(tope%i_urr) + TWO)

            ! loop over resonances localized about E_n
            RESONANCES_LOOP_4: do i_res = 1, n_res

              res % i_res = i_res
              res % E_lam = tope%local_realization(i_l) % J(i_J) % res(i_res)%E_lam
              res % Gam_n = tope%local_realization(i_l) % J(i_J) % res(i_res)%GN
              res % Gam_g = tope%local_realization(i_l) % J(i_J) % res(i_res)%GG
              res % Gam_f = tope%local_realization(i_l) % J(i_J) % res(i_res)%GF
              res % Gam_x = tope%local_realization(i_l) % J(i_J) % res(i_res)%GX
              res % Gam_t = tope%local_realization(i_l) % J(i_J) % res(i_res)%GT

              call res % calc_xs(iso)
              call accum_resonances(res, t, n, g, f, x)

            end do RESONANCES_LOOP_4
          end do TOTAL_ANG_MOM_LOOP_4
        end do ORBITAL_ANG_MOM_LOOP_4

        ! add potential scattering contribution
        call n % potential_xs(iso)
        call t % potential_xs(iso)

        ! combine resonance and File 3 components
        call mf3_self_shielding(iso, n, g, f, x, t)

        ! add the point if it's the last, otherwise cycle the loop
        if (E_hi >= E_last) then
          if (.not. tope % E_tmp % contains(E_hi))&
               call tope % E_tmp % insert(iE+1, E_hi)
          call tope % n_tmp % insert(iE+1, n % xs)
          call tope % g_tmp % insert(iE+1, g % xs)
          call tope % f_tmp % insert(iE+1, f % xs)
          call tope % x_tmp % insert(iE+1, x % xs)
          call tope % t_tmp % insert(iE+1, t % xs)
          exit

        else
          xs_hi = t % xs

        end if
 
      else
        E_hi = E_mid
        xs_hi = xs_mid

      end if
    end do

    ! pass temporary, energy-xs linked lists to dynamic vectors
    call tope % alloc_pointwise(tope % E_tmp % size())
    do iE = 1, tope % E_tmp % size()
      tope % urr_E(iE) = tope % E_tmp % get_item(iE)
      tope % urr_n(iE) = tope % n_tmp % get_item(iE)
      tope % urr_g(iE) = tope % g_tmp % get_item(iE)
      tope % urr_f(iE) = tope % f_tmp % get_item(iE)
      tope % urr_x(iE) = tope % x_tmp % get_item(iE)
      tope % urr_t(iE) = tope % t_tmp % get_item(iE)
    end do
    call tope % E_tmp % clear()
    call tope % n_tmp % clear()
    call tope % g_tmp % clear()
    call tope % f_tmp % clear()
    call tope % x_tmp % clear()
    call tope % t_tmp % clear()

    nullify(tope)

  end subroutine pointwise_urr

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! CALC_URR_XS_OTF calculates unresolved resonance region cross sections, at a
! single energy, on-the-fly from resonance parameters for a single realization
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine calc_urr_xs_otf(iso, i_nuc, E, T_K)

    type(Isotope), pointer :: tope ! isotope object pointer
    type(Nuclide), pointer :: nuc  ! nuclide object pointer
    type(Resonance) :: res ! resonance object
    type(CrossSection) :: n ! elastic scattering xs object
    type(CrossSection) :: g ! radiative capture xs object
    type(CrossSection) :: f ! fission xs object
    type(CrossSection) :: x ! competitive xs object
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
    real(8) :: inelastic_xs ! competitive cross section

    tope => isotopes(iso)
    nuc  => nuclides(i_nuc)

    ! set current temperature
    tope % T = T_K

    ! set current energy and wavenumber
    tope % E = E
    tope % k_n = wavenumber(tope % AWR, abs(tope % E))

    ! reset xs accumulators
    call flush_sigmas(t, n, g, f, x)

    ! select which resonance structure realization to use
    if (i_real_user == 0) then
      ! random realization
      i_real = ceiling(prn() * n_reals)
      if (i_real == 0) call fatal_error('i_real is sampled to be 0')
      if (i_real == n_reals + 1)&
           call fatal_error('i_real is sampled to be > n_reals')
    else
      ! user-specified realization
      i_real = i_real_user
    end if
    
    ! Get resonance parameters for a local realization about E_n
    ! loop over orbital quantum numbers
    LOC_ORBITAL_ANG_MOM_LOOP: do i_l = 1, tope % NLS(tope % i_urr)

      ! set current orbital angular momentum quantum number
      tope % L = i_l - 1

      ! get the number of contributing l-wave resonances for this l
      n_res = n_res_contrib(tope % L)

      ! loop over total angular momentum quantum numbers
      LOC_TOTAL_ANG_MOM_LOOP: do i_J = 1, tope % NJS(i_l)

        ! set current total angular momentum quantum number
        tope % J = tope % AJ(i_l) % dim1(i_J)

        ! zero the resonance counter
        res % i_res = 0

        ! find the nearest lower resonance
        if (tope % E&
             < tope % urr_resonances(i_l,i_real) % J(i_J) % res(1) % E_lam) then
          i_low = 1
        else
          i_low = binary_search(&
               tope % urr_resonances(i_l, i_real) % J(i_J) % res(:) % E_lam,&
               tope % n_lam(i_l, i_real) % dim1(i_J), tope % E)
        end if

        ! loop over the addition of resonances to this ladder
        if (i_low - n_res/2 + 1 < 1) then
          ! if we're near the lower end of the URR, need to incorporate
          ! resolved resonance region resonances in order to fix-up
          ! (i.e. smooth out) cross sections at the RRR-URR crossover
          ! energy

          ! if the RRR has resonances with this l-state
          if (i_l <= tope % NLS(tope % i_urr - 1)) then

            ! how many RRR resonances are contributing
            n_rrr_res = abs(i_low - n_res/2)

            ! loop over contributing resolved resonance region resonances
            LOC_RRR_RESONANCES_LOOP: do i_res = n_rrr_res, 1, -1
              i_rrr_res = rrr_res(iso, i_res, tope % L, tope % J)
                  
              ! fewer RRR resonances w/ this J value then needed;
              ! just generate URR resonances instead
!TODO: take however many RRR resonances there actually are, even if too few
              if (i_rrr_res == 0) exit
                  
              ! add this resolved resonance
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
          LOC_URR_LOOP: do i_res = i_low - n_res/2 + 1, i_low + n_res/2
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

      ! penetration
      tope % P_l_n = penetration(tope % L,&
           tope % k_n * tope % ac(tope % i_urr))
      
      ! resonance energy shift factor
      tope % S_l_n = shift(tope % L,&
           tope % k_n * tope % ac(tope % i_urr))
      
      ! hard-sphere phase shift
      tope % phi_l_n = phase_shift(tope % L,&
           tope % k_n * tope % AP(tope % i_urr))

      ! get the number of contributing l-wave resonances for this l
      n_res = n_res_contrib(tope % L)

      ! loop over total angular momentum quantum numbers
      TOTAL_ANG_MOM_LOOP: do i_J = 1, tope % NJS(i_l)

        ! set current total angular momentum quantum number
        tope % J = tope % AJ(i_l) % dim1(i_J)

        ! compute statistical spin factor
        tope % g_J = (TWO * tope % J + ONE) &
             / (FOUR * tope % SPI(tope % i_urr) + TWO)

        ! loop over resonances localized about e_n
        RESONANCES_LOOP: do i_res = 1, n_res

          res % i_res = i_res
          res % E_lam&
               = tope % local_realization(i_l) % J(i_J) % res(i_res) % E_lam
          res % Gam_n&
               = tope % local_realization(i_l) % J(i_J) % res(i_res) % GN
          res % Gam_g&
               = tope % local_realization(i_l) % J(i_J) % res(i_res) % GG
          res % Gam_f&
               = tope % local_realization(i_l) % J(i_J) % res(i_res) % GF
          res % Gam_x&
               = tope % local_realization(i_l) % J(i_J) % res(i_res) % GX
          res % Gam_t&
               = tope % local_realization(i_l) % J(i_J) % res(i_res) % GT

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
      call add_mf3_background(iso, i_nuc, n, g, f, x)

    else if (tope % LSSF == 1) then
      ! multipy the self-shielding factors by the infinite-dilute xs

      ! tabulated unresolved resonance parameters interpolation factor
      if (tope % E < tope % Eavg(1)) then
        iavg = 1
      else if (tope % E > tope % Eavg(tope % nEavg - 1)) then
        iavg = tope % nEavg - 1
      else
        iavg = binary_search(tope % Eavg, tope % nEavg, tope % E)
      end if

      favg = interp_factor(tope % E, &
           tope % Eavg(iavg), tope % Eavg(iavg + 1), tope % INT)

      ! interpolate averaged, infinite-dilute URR cross sections
      call interp_avg_urr_xs(favg, iso, iavg, &
           avg_urr_n_xs, avg_urr_f_xs, avg_urr_g_xs, avg_urr_x_xs)

      ! competitive xs
      ! infinite-dilute treatment of competitive xs
      inelastic_xs = micro_xs(i_nuc) % total &
           - micro_xs(i_nuc) % absorption &
           - micro_xs(i_nuc) % elastic
      if (avg_urr_x_xs > ZERO&
           .and. tope % E <= (ONE + ENDF_PRECISION) * tope % E_ex2&
           .and. competitive)&
           ! self-shielded treatment of competitive inelastic xs
           inelastic_xs = x % xs / avg_urr_x_xs * inelastic_xs

      ! elastic scattering xs
      micro_xs(i_nuc) % elastic = n % xs / avg_urr_n_xs &
           * micro_xs(i_nuc) % elastic

      ! set negative elastic xs and competitive xs to zero
      if (micro_xs(i_nuc) % elastic < ZERO) micro_xs(i_nuc) % elastic = ZERO
      if (inelastic_xs < ZERO) inelastic_xs = ZERO

      ! radiative capture xs
      ! background capture cross section
      capture_xs = micro_xs(i_nuc) % absorption - micro_xs(i_nuc) % fission
      if (avg_urr_g_xs > ZERO)&
           capture_xs = g % xs / avg_urr_g_xs * capture_xs

      ! fission xs
      if (avg_urr_f_xs > ZERO)&
           micro_xs(i_nuc) % fission = f % xs / avg_urr_f_xs&
           * micro_xs(i_nuc) % fission

      ! absorption xs
      micro_xs(i_nuc) % absorption = micro_xs(i_nuc) % fission + capture_xs

      ! total xs
      micro_xs(i_nuc) % total = micro_xs(i_nuc) % elastic&
           + micro_xs(i_nuc) % absorption + inelastic_xs

    else
      call fatal_error('ENDF-6 LSSF not allowed - must be 0 or 1.')

    end if

    ! Determine nu-fission cross section
    if (tope % fissionable)&
      micro_xs(i_nuc) % nu_fission = nu_total(nuc, E / 1.0e6_8)&
           * micro_xs(i_nuc) % fission

    nullify(tope)
    nullify(nuc)

  end subroutine calc_urr_xs_otf

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! PROB_TABLES computes probability tables for the URR
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine calc_prob_tables(iso)

    type(Isotope), pointer :: tope ! isotope object pointer
    type(ProbabilityTable), pointer :: ptable ! prob. table pointer
    type(Resonance) :: res ! resonance object
    character(6)  :: zaid_str ! ZAID number as a string
    integer :: i_b    ! batch index
    integer :: i_band ! probability band index
    integer :: i_E    ! energy grid index
    integer :: i_grid ! File 3 energy grid index
    integer :: i_h    ! history index
    integer :: i_J    ! total angular momentum quantum number index
    integer :: i_l    ! orbital quantum number index
    integer :: i_r    ! resonance index
    integer :: i_T    ! temperature index
    integer :: iso    ! isotope index
    integer :: n_res  ! number of resonances to include for a given l-wave
    integer :: i_mag  ! index in array of xs samples based on total xs magnitude
    integer :: hits_per_band ! number of hits in each equiprobable band
    integer :: avg_unit = 98 ! avg xs output file unit
    integer :: tab_unit = 99 ! tables output file unit
    real(8) :: E        ! neutron lab energy [eV]
    real(8) :: fmf3     ! File 3 energy grid interpolation factor
    real(8) :: xs_t_min ! min realized total xs
    real(8) :: xs_t_max ! max realized total xs

    xs_t_min = 1.0e6_8
    xs_t_max = XS_CUTOFF

    tope => isotopes(iso)

    write(zaid_str, '(I6)') tope % ZAI
    if (master)&
      write(*,*) 'Generating probability tables for ZAID = '//&
      trim(adjustl(zaid_str))

    if (write_urr_tables) then
      open(unit = tab_unit, file = trim(adjustl(zaid_str))//'-urr-tables.dat')
      write(tab_unit, '("ENDF-6 Path:")', advance='no')
      write(tab_unit, *) trim(adjustl(path_endf))//trim(adjustl(endf_files(iso)))
      write(tab_unit, '("Resonance Formalism:",i2)') formalism
      write(tab_unit, '("s, p, d, f Resonances:",1x,4i3)')&
           l_waves(1), l_waves(2), l_waves(3), l_waves(4)
      write(tab_unit, '("Competitive Reaction Structure:",L2)') competitive
      write(tab_unit, '("Parameter Energy Dependence &
           &(1-neutron, 2-resonance):",i2)') represent_params
      write(tab_unit, '("Faddeeva Evaluation:",i2)') w_eval
      write(tab_unit, '("Target Relative Tolerance on All Avg. Partial xs:",es13.6)') tol_avg_urr
      write(tab_unit, '("Energies:",i7)') tope % nE_tabs
      write(tab_unit, '("Temperatures:",i3)') tope % nT_tabs
      write(tab_unit, '("Bands:",i10)') tope % n_bands
    end if

    if (write_avg_urr_xs) then
      open(unit = avg_unit, file = trim(adjustl(zaid_str))//'-avg-urr-xs.dat')
      write(avg_unit, '("ENDF-6 Path:")', advance='no')
      write(avg_unit, *) trim(adjustl(path_endf))//trim(adjustl(endf_files(iso)))
      write(avg_unit, '("Resonance Formalism:",i2)') formalism
      write(avg_unit, '("s, p, d, f Resonances:",1x,4i3)')&
           l_waves(1), l_waves(2), l_waves(3), l_waves(4)
      write(avg_unit, '("Competitive Reaction Structure:",L2)') competitive
      write(avg_unit, '("Parameter Energy Dependence &
           &(1-neutron, 2-resonance):",i2)') represent_params
      write(avg_unit, '("Faddeeva Evaluation:",i2)') w_eval
      write(avg_unit, '("Target Relative Tolerance on All Avg. Partial xs:",ES13.6)') tol_avg_urr
      write(avg_unit, '("Energies:",i3)') tope % nE_tabs
      write(avg_unit, '(A13,A13,A13,A13,A13,A13)') &
           'energy [eV]', 'total', 'elastic', 'capture', 'fission','competitive'
    end if

    ! loop over energy mesh
    ENERGY_LOOP: do i_E = 1, tope % nE_tabs

      tope % E = tope % E_tabs(i_E)
      E = tope % E
      tope % k_n = wavenumber(tope % AWR, abs(tope % E))
      tope % k_lam = tope % k_n 
 
      ! reset accumulator of statistics
      call tope % flush_ptable_stats(i_E)

      i_b = 0

      allocate(tope % xs_samples(histories_avg_urr, tope % nT_tabs))
      tope % xs_samples(:,:) % t = ZERO
      tope % xs_samples(:,:) % n = ZERO
      tope % xs_samples(:,:) % g = ZERO
      tope % xs_samples(:,:) % f = ZERO
      tope % xs_samples(:,:) % x = ZERO

      ! loop over batches until convergence
      BATCH_LOOP: do

        i_b = i_b + 1

        ! reset batch accumulators
        call tope % flush_batches()

        ! loop over realizations
        HISTORY_LOOP: do i_h = 1, histories_avg_urr

          ! reset accumulator of histories
          call tope % flush_histories()

          ! Get resonance parameters for a local realization about E_n
          ! loop over orbital quantum numbers
          LOC_ORBITAL_ANG_MOM_LOOP: do i_l = 1, tope % NLS(tope % i_urr)

            ! set current orbital angular momentum quantum number
            tope % L = i_l - 1

            ! penetration
            tope % P_l_n = penetration(tope % L,&
                 tope % k_n * tope % ac(tope % i_urr))

            ! resonance energy shift factor
            tope % S_l_n = shift(tope % L,&
                 tope % k_n * tope % ac(tope % i_urr))
            
            ! hard-sphere phase shift
            tope % phi_l_n = phase_shift(tope % L,&
                 tope % k_n * tope % AP(tope % i_urr))

           ! get the number of contributing l-wave resonances for this l
            n_res = n_res_contrib(tope % L)

            ! loop over total angular momentum quantum numbers
            LOC_TOTAL_ANG_MOM_LOOP: do i_J = 1, tope % NJS(i_l)

              ! set current total angular momentum quantum number
              tope % J = tope % AJ(i_l) % dim1(i_J)

              ! set current partial width degrees of freedom
              tope % AMUX = int(tope % DOFX(i_l) % dim1(i_J))
              tope % AMUN = int(tope % DOFN(i_l) % dim1(i_J))
              tope % AMUG = int(tope % DOFG(i_l) % dim1(i_J))
              tope % AMUF = int(tope % DOFF(i_l) % dim1(i_J))

              ! zero the resonance counter
              res % i_res = 0

              ! set mean URR parameters to neutron energy
              call set_mean_parameters(iso, E, i_l, i_J)

              ! sample unresolved resonance parameters for this spin
              ! sequence, at this energy
              tope % k_lam = tope % k_n
              tope % P_l_lam = penetration(tope % L,&
                   tope % k_lam * tope % ac(tope % i_urr)) 
              call res % sample_parameters(iso, i_l, i_J)

              ! loop over the addition of resonances to this ladder
              LOC_RESONANCES_LOOP: do i_r = 1, n_res

                if (represent_params == E_RESONANCE) then
                  ! interpolate mean URR parameters to current resonance energy
                  call set_mean_parameters(iso, res % E_lam, i_l, i_J)
                  tope % k_lam = wavenumber(tope % AWR, abs(res % E_lam))
                  tope % P_l_lam = penetration(tope % L,&
                       tope % k_lam * tope % ac(tope % i_urr))
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

            ! penetration
            tope % P_l_n = penetration(tope % L,&
                 tope % k_n * tope % ac(tope % i_urr))
      
            ! resonance energy shift factor
            tope % S_l_n = shift(tope % L,&
                 tope % k_n * tope % ac(tope % i_urr))
      
            ! hard-sphere phase shift
            tope % phi_l_n = phase_shift(tope % L,&
                 tope % k_n * tope % AP(tope % i_urr))

            ! get the number of contributing l-wave resonances for this l
            n_res = n_res_contrib(tope % L)

            ! loop over total angular momentum quantum numbers
            TOTAL_ANG_MOM_LOOP: do i_J = 1, tope % NJS(i_l)

              ! set current total angular momentum quantum number
              tope % J = tope % AJ(i_l) % dim1(i_J)

              ! compute statistical spin factor
              tope % g_J = (TWO * tope % J + ONE) &
                   / (FOUR * tope % SPI(tope % i_urr) + TWO)

              ! loop over resonances localized about e_n
              RESONANCES_LOOP: do i_r = 1, n_res

                res % i_res = i_r
                res % E_lam&
                     = tope % local_realization(i_l) % J(i_J) % res(i_r) % E_lam
                res % Gam_n&
                     = tope % local_realization(i_l) % J(i_J) % res(i_r) % GN
                res % Gam_g&
                     = tope % local_realization(i_l) % J(i_J) % res(i_r) % GG
                res % Gam_f&
                     = tope % local_realization(i_l) % J(i_J) % res(i_r) % GF
                res % Gam_x&
                     = tope % local_realization(i_l) % J(i_J) % res(i_r) % GX
                res % Gam_t&
                     = tope % local_realization(i_l) % J(i_J) % res(i_r) % GT

                TEMPERATURES_LOOP: do i_T = 1, tope % nT_tabs

                  ! set current temperature
                  tope % T = tope % T_tabs(i_T)

                  ! calculate the contribution to the partial cross sections,
                  ! at this energy, from an additional resonance
                  call res % calc_xs(iso)

                  ! add this contribution to the accumulated partial cross
                  ! section values built up from all resonances
! TODO: consider moving t outside of loop
                  call tope % res_contrib(res, i_E, i_T)

                end do TEMPERATURES_LOOP
              end do RESONANCES_LOOP
            end do TOTAL_ANG_MOM_LOOP
          end do ORBITAL_ANG_MOM_LOOP

          TEMPERATURES_LOOPb: do i_T = 1, tope % nT_tabs

            ! add potential scattering contribution
            call tope % prob_tables(i_E, i_T) % avg_n % potential_xs(iso)
            call tope % prob_tables(i_E, i_T) % avg_t % potential_xs(iso)

            ! set negative elastic xs to zero
            if (tope % prob_tables(i_E, i_T) % avg_n % xs < ZERO) then
              tope % prob_tables(i_E, i_T) % avg_t % xs &
                   = tope % prob_tables(i_E, i_T) % avg_t % xs &
                   + abs(tope % prob_tables(i_E, i_T) % avg_n % xs)
              tope % prob_tables(i_E, i_T) % avg_n % xs = ZERO
            end if

            ! use MF3 fission cross section if requested
            if (background == FALSE) then
              continue
            else
              tope % prob_tables(i_E, i_T) % avg_t % xs &
                   = tope % prob_tables(i_E, i_T) % avg_t % xs &
                   - tope % prob_tables(i_E, i_T) % avg_f % xs
              tope % prob_tables(i_E, i_T) % avg_f % xs = ZERO
              if (background == ENDFFILE .and. allocated(tope % MF3_f_e)) then
                if (E >= tope % MF3_f_e(1)) then
                  i_grid = binary_search(tope % MF3_f_e,size(tope % MF3_f_e),E)
                  if (tope % INT == LINEAR_LINEAR &
                       .or. (tope % MF3_f(i_grid) > XS_CUTOFF &
                       .and. tope % MF3_f(i_grid + 1) > XS_CUTOFF)) then
                    fmf3 = interp_factor(E, tope % MF3_f_e(i_grid), &
                         tope % MF3_f_e(i_grid + 1), tope % INT)
                    tope % prob_tables(i_E, i_T) % avg_f % xs &
                         = interpolator(fmf3, tope % MF3_f(i_grid), &
                         tope % MF3_f(i_grid + 1), tope % INT)
                    tope % prob_tables(i_E, i_T) % avg_t % xs &
                         = tope % prob_tables(i_E, i_T) % avg_t % xs &
                         + tope % prob_tables(i_E, i_T) % avg_f % xs
                  end if
                end if
              end if
            end if
            if (tope % prob_tables(i_E, i_T) % avg_f % xs < ZERO)&
                 call fatal_error('Negative fission xs encountered')

            ! use MF3 competitive cross section if requested
            if (background == FALSE) then
              continue
            else
              tope % prob_tables(i_E, i_T) % avg_t % xs &
                   = tope % prob_tables(i_E, i_T) % avg_t % xs &
                   - tope % prob_tables(i_E, i_T) % avg_x % xs
              tope % prob_tables(i_E, i_T) % avg_x % xs = ZERO
              if (background == ENDFFILE .and. allocated(tope % MF3_x_e)) then
                if (E >= tope % MF3_x_e(1)) then
                  i_grid = binary_search(tope % MF3_x_e,size(tope % MF3_x_e),E)
                  if (tope % INT == LINEAR_LINEAR &
                       .or. (tope % MF3_x(i_grid) > XS_CUTOFF &
                       .and. tope % MF3_x(i_grid + 1) > XS_CUTOFF)) then
                    fmf3 = interp_factor(E, tope % MF3_x_e(i_grid), &
                         tope % MF3_x_e(i_grid + 1), tope % INT)
                    tope % prob_tables(i_E, i_T) % avg_x % xs &
                         = interpolator(fmf3, tope % MF3_x(i_grid), &
                         tope % MF3_x(i_grid + 1), tope % INT)
                    tope % prob_tables(i_E, i_T) % avg_t % xs &
                         = tope % prob_tables(i_E, i_T) % avg_t % xs &
                         + tope % prob_tables(i_E, i_T) % avg_x % xs
                  end if
                end if
              end if
            end if
    
            ! set negative competitive xs to zero
            if (tope % prob_tables(i_E, i_T) % avg_x % xs < ZERO) then
              call fatal_error('Negative competitive xs encountered')
              tope % prob_tables(i_E, i_T) % avg_t % xs &
                   = tope % prob_tables(i_E, i_T) % avg_t % xs &
                   + abs(tope % prob_tables(i_E, i_T) % avg_x % xs)
              tope % prob_tables(i_E, i_T) % avg_x % xs = ZERO
            end if

            ! Set min and max xs values encountered
            if (tope % prob_tables(i_E, i_T) % avg_t % xs < xs_t_min)&
                 xs_t_min = tope % prob_tables(i_E, i_T) % avg_t % xs
            if (tope % prob_tables(i_E, i_T) % avg_t % xs > xs_t_max)&
                 xs_t_max = tope % prob_tables(i_E, i_T) % avg_t % xs

            ! accumulate the result of this history
            call tope % accum_history(i_E, i_T)

            ! find index where this sample belongs based on total cross section
            if (i_b == 1 .and. i_h == 1) then
              i_mag = histories_avg_urr
            else
              if (tope % prob_tables(i_E, i_T) % avg_t % xs&
                   > tope % xs_samples(i_b * histories_avg_urr, i_T) % t) then
                i_mag = i_b * histories_avg_urr
              else
               i_mag = binary_search(tope % xs_samples(:, i_T) % t,&
                     i_b * histories_avg_urr,&
                     tope % prob_tables(i_E, i_T) % avg_t % xs)
              end if
            end if

            ! insert total and conditional partial cross sections
            tope%xs_samples(1 : i_mag - 1, i_T) = tope % xs_samples(2:i_mag,i_T)
            tope%xs_samples(i_mag,i_T) % t = tope%prob_tables(i_E,i_T)% avg_t%xs
            tope%xs_samples(i_mag,i_T) % n = tope%prob_tables(i_E,i_T)% avg_n%xs
            tope%xs_samples(i_mag,i_T) % g = tope%prob_tables(i_E,i_T)% avg_g%xs
            tope%xs_samples(i_mag,i_T) % f = tope%prob_tables(i_E,i_T)% avg_f%xs
            tope%xs_samples(i_mag,i_T) % x = tope%prob_tables(i_E,i_T)% avg_x%xs
          
          end do TEMPERATURES_LOOPb
        end do HISTORY_LOOP

        call move_alloc(tope % xs_samples, xs_samples_tmp)
        allocate(tope % xs_samples((i_b + 1) * histories_avg_urr, tope % nT_tabs))
        tope % xs_samples(:,:) % t = ZERO
        tope % xs_samples(:,:) % n = ZERO
        tope % xs_samples(:,:) % g = ZERO
        tope % xs_samples(:,:) % f = ZERO
        tope % xs_samples(:,:) % x = ZERO
        tope % xs_samples(histories_avg_urr+1 : (i_b+1)*histories_avg_urr, :)&
             = xs_samples_tmp
        deallocate(xs_samples_tmp)
        
        ! accumulate the result of this batch
        call tope % accum_batch(i_E)

        ! calculate statistics for this batch
        call tope % calc_stats(i_b)

        if ((i_b > min_batches_avg_urr &
             .and. max(maxval(tope % prob_tables(i_E, :) % avg_t % rel_unc), &
                       maxval(tope % prob_tables(i_E, :) % avg_n % rel_unc), &
                       maxval(tope % prob_tables(i_E, :) % avg_g % rel_unc), &
                       maxval(tope % prob_tables(i_E, :) % avg_f % rel_unc), &
                       maxval(tope % prob_tables(i_E, :) % avg_x % rel_unc)) &
             < tol_avg_urr) .or. i_b == max_batches_avg_urr) exit

      end do BATCH_LOOP

      call move_alloc(tope % xs_samples, xs_samples_tmp)
      allocate(tope % xs_samples((i_b) * histories_avg_urr, tope % nT_tabs))
      tope % xs_samples(:,:)&
           = xs_samples_tmp(histories_avg_urr+1:(i_b+1)*histories_avg_urr,:)
      deallocate(xs_samples_tmp)
      do i_T = 1, tope % nT_tabs
        hits_per_band = nint(i_b * histories_avg_urr / dble(tope % n_bands))
        do i_band = 1, tope % n_bands - 1
          tope % prob_tables(i_E, i_T) % t(i_band) % cnt_mean&
               = dble(hits_per_band) / dble(i_b * histories_avg_urr)
          tope % prob_tables(i_E, i_T) % t(i_band) % xs_mean&
               = sum(tope % xs_samples((i_band - 1) * hits_per_band + 1&
               : i_band * hits_per_band, i_T) % t) / dble(hits_per_band)
          tope % prob_tables(i_E, i_T) % n(i_band) % xs_mean&
               = sum(tope % xs_samples((i_band - 1) * hits_per_band + 1&
               : i_band * hits_per_band, i_T) % n) / dble(hits_per_band)
          tope % prob_tables(i_E, i_T) % g(i_band) % xs_mean&
               = sum(tope % xs_samples((i_band - 1) * hits_per_band + 1&
               : i_band * hits_per_band, i_T) % g) / dble(hits_per_band)
          tope % prob_tables(i_E, i_T) % f(i_band) % xs_mean&
               = sum(tope % xs_samples((i_band - 1) * hits_per_band + 1&
               : i_band * hits_per_band, i_T) % f) / dble(hits_per_band)
          tope % prob_tables(i_E, i_T) % x(i_band) % xs_mean&
               = sum(tope % xs_samples((i_band - 1) * hits_per_band + 1&
               : i_band * hits_per_band, i_T) % x) / dble(hits_per_band)
        end do
        tope % prob_tables(i_E, i_T) % t(tope % n_bands) % cnt_mean&
             = dble(i_b * histories_avg_urr - (tope % n_bands - 1) * hits_per_band)&
             / dble(i_b * histories_avg_urr)
        tope % prob_tables(i_E, i_T) % t(tope % n_bands) % xs_mean&
             = sum(tope % xs_samples((tope % n_bands - 1) * hits_per_band + 1&
             : i_b * histories_avg_urr, i_T) % t)&
             /dble(i_b * histories_avg_urr - (tope % n_bands - 1) * hits_per_band) 
        tope % prob_tables(i_E, i_T) % n(tope % n_bands) % xs_mean&
             = sum(tope % xs_samples((tope % n_bands - 1) * hits_per_band + 1&
             : i_b * histories_avg_urr, i_T) % n)&
             /dble(i_b * histories_avg_urr - (tope % n_bands - 1) * hits_per_band) 
        tope % prob_tables(i_E, i_T) % g(tope % n_bands) % xs_mean&
             = sum(tope % xs_samples((tope % n_bands - 1) * hits_per_band + 1&
             : i_b * histories_avg_urr, i_T) % g)&
             /dble(i_b * histories_avg_urr - (tope % n_bands - 1) * hits_per_band) 
        tope % prob_tables(i_E, i_T) % f(tope % n_bands) % xs_mean&
             = sum(tope % xs_samples((tope % n_bands - 1) * hits_per_band + 1&
             : i_b * histories_avg_urr, i_T) % f)&
             /dble(i_b * histories_avg_urr - (tope % n_bands - 1) * hits_per_band) 
        tope % prob_tables(i_E, i_T) % x(tope % n_bands) % xs_mean&
             = sum(tope % xs_samples((tope % n_bands - 1) * hits_per_band + 1&
             : i_b * histories_avg_urr, i_T) % x)&
             /dble(i_b * histories_avg_urr - (tope % n_bands - 1) * hits_per_band) 
      end do

      deallocate(tope % xs_samples)

      ! write probability tables out to a file
      if (write_urr_tables) then
        if (master)&
          write(*,'(A32,ES10.3,A3)') 'Generated probability tables at', &
               tope % E_tabs(i_E), ' eV'
        do i_T = 1, tope % nT_tabs
          tope % T = tope % T_tabs(i_T)
          write(tab_unit, '(A13,ES13.6)') 'E [eV]', tope % E_tabs(i_E)
          write(tab_unit, '(A13,ES13.6)') 'T [K]', tope % T
          write(tab_unit, '(8A13)') 'lower [b]', 'upper [b]',&
               'prob', 'total', 'elastic', 'capture', 'fission', 'competitive'
          ptable => tope % prob_tables(i_E, i_T)
          do i_band = 1, tope % n_bands
            if (i_band == 1) then
              write(tab_unit, '(ES13.6,A13,6ES13.6)')&
                   xs_t_min, '',&
                   ptable % t(i_band) % cnt_mean,&
                   ptable % t(i_band) % xs_mean,&
                   ptable % n(i_band) % xs_mean,&
                   ptable % g(i_band) % xs_mean,&
                   ptable % f(i_band) % xs_mean,&
                   ptable % x(i_band) % xs_mean
            else if (i_band == tope % n_bands) then
              write(tab_unit, '(A13,7ES13.6)')&
                   '', xs_t_max,&
                   ptable % t(i_band) % cnt_mean,&
                   ptable % t(i_band) % xs_mean,&
                   ptable % n(i_band) % xs_mean,&
                   ptable % g(i_band) % xs_mean,&
                   ptable % f(i_band) % xs_mean,&
                   ptable % x(i_band) % xs_mean
            else
              write(tab_unit, '(2A13,6ES13.6)')&
                   '', '',&
                   ptable % t(i_band) % cnt_mean,&
                   ptable % t(i_band) % xs_mean,&
                   ptable % n(i_band) % xs_mean,&
                   ptable % g(i_band) % xs_mean,&
                   ptable % f(i_band) % xs_mean,&
                   ptable % x(i_band) % xs_mean
            end if
          end do
          write(tab_unit,'(2A13,6ES13.6)')&
               'batch', 'averaged xs', ONE,&
               ptable % avg_t % xs_mean,&
               ptable % avg_n % xs_mean,&
               ptable % avg_g % xs_mean,&
               ptable % avg_f % xs_mean,&
               ptable % avg_x % xs_mean
          write(tab_unit,'(I13,A13,6ES13.6)')&
               i_b, '1sigma', ZERO,&
               ptable % avg_t % xs_sem,&
               ptable % avg_n % xs_sem,&
               ptable % avg_g % xs_sem,&
               ptable % avg_f % xs_sem,&
               ptable % avg_x % xs_sem
        end do
      end if

      ! write averaged URR cross sections out to a file
      if (write_avg_urr_xs) then
        ptable => tope % prob_tables(i_E, 1)
        write(avg_unit, '(6ES13.6)')&
             tope % E_tabs(i_E),&
             ptable % avg_t % xs_mean,&
             ptable % avg_n % xs_mean,&
             ptable % avg_g % xs_mean,&
             ptable % avg_f % xs_mean,&
             ptable % avg_x % xs_mean
      end if

    end do ENERGY_LOOP

    if (write_urr_tables) close(tab_unit)
    if (write_avg_urr_xs) close(avg_unit)

    nullify(ptable)

  end subroutine calc_prob_tables

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! LOAD_PROB_TABLES loads in pre-generated URR probability tables
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine load_prob_tables(i)

    type(Isotope), pointer :: tope ! isotope object pointer
    logical :: file_exists ! does probability table file exist?
    integer :: i             ! isotope index
    integer :: i_E           ! energy index
    integer :: i_T           ! temperature index
    integer :: i_b           ! probability band index
    integer :: tab_unit = 99 ! tables output file unit
    character(7)   :: readable ! is probability table file readable?
    character(6)   :: zaid_str ! ZAID number as a string
    character(120) :: rec      ! file record

    tope => isotopes(i)
    write(zaid_str, '(i6)') tope % ZAI

    ! check that file exists and is readable
    inquire(file = trim(path_prob_tables)//trim(adjustl(zaid_str))//&
         &'-urr-tables.dat', exist = file_exists, read = readable)
    if (.not. file_exists) then
      call fatal_error('Probability table file '//trim(adjustl(zaid_str))//&
           &'-urr-tables.dat does not exist.')
    else if (readable(1:3) == 'NO') then
      call fatal_error('Probability table file '//trim(adjustl(zaid_str))//&
           &'-urr-tables.dat is not readable.  Change file permissions with &
           &chmod command.')
    end if

    ! open probability table file
    if (master)&
         write(*,*) 'Loading probability tables for ZAID ='//zaid_str
    open(unit = tab_unit, file =&
         trim(path_prob_tables)//trim(adjustl(zaid_str))//'-urr-tables.dat')

10  format(A120)
    ! ENDF-6 filepath
    read(tab_unit, 10) rec
    
    ! resonance formalism
    read(tab_unit, 10) rec
     
    ! number of s, p, d, f wave resonances
    read(tab_unit, 10) rec

    ! structured competitive cross section?
    read(tab_unit, 10) rec

    ! parameter energy dependence
    read(tab_unit, 10) rec

    ! Faddeeva function evaluation
    read(tab_unit, 10) rec

    ! averaged partial cross section 1sigma SEM tolerance
    read(tab_unit, 10) rec

    ! number of energies
    read(tab_unit, 10) rec
    read(rec(10:80), '(i71)') tope % nE_tabs

    ! number of temperatures
    read(tab_unit, 10) rec
    read(rec(14:80), '(i67)') tope % nT_tabs

    ! number of probability-xs bands
    read(tab_unit, 10) rec
    read(rec(7:80), '(i74)') tope % n_bands

    ! allocate probability tables
    allocate(tope % E_tabs(tope % nE_tabs))
    allocate(tope % T_tabs(tope % nT_tabs))
    allocate(tope % prob_tables(tope % nE_tabs, tope % nT_tabs))
    do i_E = 1, tope % nE_tabs
      do i_T = 1, tope % nT_tabs
        allocate(tope % prob_tables(i_E, i_T) % t(tope % n_bands))
        allocate(tope % prob_tables(i_E, i_T) % n(tope % n_bands))
        allocate(tope % prob_tables(i_E, i_T) % g(tope % n_bands))
        allocate(tope % prob_tables(i_E, i_T) % f(tope % n_bands))
        allocate(tope % prob_tables(i_E, i_T) % x(tope % n_bands))
      end do
    end do

    ! read probability tables
    do i_E = 1, tope % nE_tabs
      do i_T = 1, tope % nT_tabs
        ! read energy        
        read(tab_unit, 10) rec
        read(rec(14:26), '(es13.6)') tope % E_tabs(i_E)

        ! read temperature
        read(tab_unit, 10) rec
        read(rec(14:26), '(es13.6)') tope % T_tabs(i_T)

        ! read column labels
        read(tab_unit, 10) rec

        ! read cross sections
        do i_b = 1, tope % n_bands
          read(tab_unit, 10) rec
          read(rec(27:104), '(6es13.6)')&
               tope % prob_tables(i_E, i_T) % t(i_b) % cnt_mean,&
               tope % prob_tables(i_E, i_T) % t(i_b) % xs_mean,&
               tope % prob_tables(i_E, i_T) % n(i_b) % xs_mean,&
               tope % prob_tables(i_E, i_T) % g(i_b) % xs_mean,&
               tope % prob_tables(i_E, i_T) % f(i_b) % xs_mean,&
               tope % prob_tables(i_E, i_T) % x(i_b) % xs_mean
        end do

        ! read average cross sections
        read(tab_unit, 10) rec

        ! read 1sigma SEM values
        read(tab_unit, 10) rec

      end do
    end do

  end subroutine load_prob_tables

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

    type(Isotope), pointer :: tope ! isotope object pointer
    type(Nuclide), pointer :: nuc  ! nuclide object pointer
    type(Resonance) :: res ! resonance object
    type(CrossSection) :: n ! elastic scattering xs object
    type(CrossSection) :: g ! radiative capture xs object
    type(CrossSection) :: f ! fission xs object
    type(CrossSection) :: x ! competitive xs object
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
    real(8) :: inelastic_xs ! competitive cross section
    real(8) :: capture_xs   ! radiative capture cross section
    real(8) :: E            ! neutron energy [eV]
    real(8) :: m            ! pointwise xs energy interpolation factor
    real(8) :: favg         ! average cross section interpolation factor
    real(8) :: T_K          ! isotope temperature [K]

    tope => isotopes(iso)
    nuc  => nuclides(i_nuc)

    if (tope % point_urr_xs) then
      i_E = binary_search(tope % urr_E, size(tope % urr_E), E)
      m = interp_factor(E, tope % urr_E(i_E), tope % urr_E(i_E + 1), &
           LINEAR_LINEAR)
      micro_xs(i_nuc) % elastic = interpolator(m, &
           tope % urr_n(i_E), tope % urr_n(i_E + 1), LINEAR_LINEAR)
      micro_xs(i_nuc) % fission = interpolator(m, &
           tope % urr_f(i_E), tope % urr_f(i_E + 1), LINEAR_LINEAR)
      micro_xs(i_nuc) % absorption = interpolator(m, &
           tope % urr_g(i_E) + tope % urr_f(i_E), &
           tope % urr_g(i_E + 1) + tope % urr_f(i_E + 1), LINEAR_LINEAR)
      inelastic_xs = ZERO
      if (E >= tope % E_ex1) inelastic_xs = interpolator(m,&
           tope % urr_x(i_E), tope % urr_x(i_E + 1), LINEAR_LINEAR)
      micro_xs(i_nuc) % total = micro_xs(i_nuc) % elastic &
           + micro_xs(i_nuc) % absorption &
           + inelastic_xs

      ! Determine nu-fission cross section
      if (tope % fissionable) micro_xs(i_nuc) % nu_fission&
           = nu_total(nuc, E / 1.0e6_8) * micro_xs(i_nuc) % fission
      return
    end if

    ! set current temperature and neutron energy
    tope % T = T_K
    tope % E = E
    tope % k_n = wavenumber(tope % AWR, abs(tope % E))
    tope % k_lam = tope % k_n

    ! reset xs accumulators
    call flush_sigmas(t, n, g, f, x)

    ! Get resonance parameters for a local realization about E_n
    ! loop over orbital quantum numbers
    LOC_ORBITAL_ANG_MOM_LOOP: do i_l = 1, tope % NLS(tope % i_urr)

      ! set current orbital angular momentum quantum number
      tope % L = i_l - 1

      ! penetration
      tope % P_l_n = penetration(tope % L,&
           tope % k_n * tope % ac(tope % i_urr))
      
      ! resonance energy shift factor
      tope % S_l_n = shift(tope % L,&
           tope % k_n * tope % ac(tope % i_urr))
          
      ! hard-sphere phase shift
      tope % phi_l_n = phase_shift(tope % L,&
           tope % k_n * tope % AP(tope % i_urr))

      ! get the number of contributing l-wave resonances for this l
      n_res = n_res_contrib(tope % L)

      ! loop over total angular momentum quantum numbers
      LOC_TOTAL_ANG_MOM_LOOP: do i_J = 1, tope % NJS(i_l)

        ! set current total angular momentum quantum number
        tope % J = tope % AJ(i_l) % dim1(i_J)

        ! set current partial width degrees of freedom
        tope % AMUX = int(tope % DOFX(i_l) % dim1(i_J))
        tope % AMUN = int(tope % DOFN(i_l) % dim1(i_J))
        tope % AMUG = int(tope % DOFG(i_l) % dim1(i_J))
        tope % AMUF = int(tope % DOFF(i_l) % dim1(i_J))

        ! zero the resonance counter
        res % i_res = 0

        ! set mean URR parameters to neutron energy
        call set_mean_parameters(iso, E, i_l, i_J)

        ! sample unresolved resonance parameters for this spin
        ! sequence, at this energy
        tope % k_lam = tope % k_n
        tope % P_l_lam = penetration(tope % L,&
             tope % k_lam * tope % ac(tope % i_urr)) 
        call res % sample_parameters(iso, i_l, i_J)

        ! loop over the addition of resonances to this ladder
        LOC_RESONANCES_LOOP: do i_r = 1, n_res

          if (represent_params == E_RESONANCE) then
            ! interpolate mean URR parameters to current resonance energy
            call set_mean_parameters(iso, res % E_lam, i_l, i_J)
            tope % k_lam = wavenumber(tope % AWR, abs(res % E_lam))
            tope % P_l_lam = penetration(tope % L,&
                 tope % k_lam * tope % ac(tope % i_urr))
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

      ! penetration
      tope % P_l_n = penetration(tope % L,&
           tope % k_n * tope % ac(tope % i_urr))
      
      ! resonance energy shift factor
      tope % S_l_n = shift(tope % L,&
           tope % k_n * tope % ac(tope % i_urr))
      
      ! hard-sphere phase shift
      tope % phi_l_n = phase_shift(tope % L,&
           tope % k_n * tope % AP(tope % i_urr))

      ! get the number of contributing l-wave resonances for this l
      n_res = n_res_contrib(tope % L)

      ! loop over total angular momentum quantum numbers
      TOTAL_ANG_MOM_LOOP: do i_J = 1, tope % NJS(i_l)

        ! set current total angular momentum quantum number
        tope % J = tope % AJ(i_l) % dim1(i_J)

        ! compute statistical spin factor
        tope % g_J = (TWO * tope % J + ONE) &
             / (FOUR * tope % SPI(tope % i_urr) + TWO)

        ! loop over resonances localized about e_n
        RESONANCES_LOOP: do i_r = 1, n_res

          res % i_res = i_r
          res % E_lam&
               = tope % local_realization(i_l) % J(i_J) % res(i_r) % E_lam
          res % Gam_n&
               = tope % local_realization(i_l) % J(i_J) % res(i_r) % GN
          res % Gam_g&
               = tope % local_realization(i_l) % J(i_J) % res(i_r) % GG
          res % Gam_f&
               = tope % local_realization(i_l) % J(i_J) % res(i_r) % GF
          res % Gam_x&
               = tope % local_realization(i_l) % J(i_J) % res(i_r) % GX
          res % Gam_t&
               = tope % local_realization(i_l) % J(i_J) % res(i_r) % GT

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
      call add_mf3_background(iso, i_nuc, n, g, f, x)

    ! MF3 contains evaluator-supplied background xs values that we multipy
    ! the self-shielding factors computed from MF2 by
    elseif (tope % LSSF == 1) then

      ! determine energy index
      if (tope % E < tope % Eavg(1)) then
        iavg = 1
      else if (tope % E > tope % Eavg(tope % nEavg - 1)) then
        iavg = tope % nEavg - 1
      else
        iavg = binary_search(tope % Eavg, tope % nEavg, tope % E)
      end if

      ! tabulated unresolved resonance parameters interpolation factor
      favg = interp_factor(E, tope % Eavg(iavg), tope % Eavg(iavg + 1), &
           tope % INT)

      ! interpolate infinite-dilute URR xs
      call interp_avg_urr_xs(favg, iso, iavg, &
           avg_urr_n_xs, avg_urr_f_xs, avg_urr_g_xs, avg_urr_x_xs)

      inelastic_xs = micro_xs(i_nuc) % total &
           - micro_xs(i_nuc) % absorption &
           - micro_xs(i_nuc) % elastic
      if (avg_urr_x_xs > ZERO&
           .and. tope % E <= (ONE + ENDF_PRECISION) * tope % E_ex2&
           .and. competitive)&
           ! self-shielded treatment of competitive inelastic cross section
           inelastic_xs = x % xs / avg_urr_x_xs * inelastic_xs

      micro_xs(i_nuc) % elastic = n % xs / avg_urr_n_xs &
           * micro_xs(i_nuc) % elastic

      ! set negative elastic xs and competitive xs to zero
      if (micro_xs(i_nuc) % elastic < ZERO) micro_xs(i_nuc) % elastic = ZERO
      if (inelastic_xs < ZERO) inelastic_xs = ZERO

      capture_xs = micro_xs(i_nuc) % absorption - micro_xs(i_nuc) % fission
      if (avg_urr_g_xs > ZERO)&
           capture_xs = g % xs / avg_urr_g_xs * capture_xs

      if (avg_urr_f_xs > ZERO) micro_xs(i_nuc) % fission&
           = f % xs / avg_urr_f_xs * micro_xs(i_nuc) % fission

      micro_xs(i_nuc) % absorption = micro_xs(i_nuc) % fission + capture_xs

      micro_xs(i_nuc) % total = micro_xs(i_nuc) % elastic &
           + micro_xs(i_nuc) % absorption + inelastic_xs

    else
      call fatal_error('Self-shielding flag (LSSF) must be 0 or 1')

    end if

    ! Determine nu-fission cross section
    if (tope % fissionable) micro_xs(i_nuc) % nu_fission&
         = nu_total(nuc, E / 1.0e6_8) * micro_xs(i_nuc) % fission

    ! set last neutron energy
    tope % E_last = tope % E

    nullify(tope)
    nullify(nuc)

  end subroutine calculate_urr_xs_otf

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! CALCULATE_PROB_BAND_XS calculates a URR cross section from the
! OpenMC-computed probability tables
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine calculate_prob_band_xs(i_so, i_nuc, E, T)

    type(Isotope),  pointer :: tope
    type(Nuclide),  pointer :: nuc
    integer, intent(in) :: i_so  ! index into isotopes array
    integer, intent(in) :: i_nuc ! index into nuclides array
    real(8), intent(in) :: E     ! energy
    real(8), intent(in) :: T     ! temperature [K]
    integer :: i            ! loop index
    integer :: i_E          ! index for energy
    integer :: i_Tlow       ! index for lower temperature bound
    integer :: i_Tup        ! index for upper temperature bound
    integer :: i_low        ! band index at lower bounding energy
    integer :: i_up         ! band index at upper bounding energy
    integer :: same_nuc_idx ! index of same nuclide
    integer :: iavg         ! average cross section index
    real(8) :: fE           ! energy interpolation factor
    real(8) :: fT           ! temperature interpolation factor
    real(8) :: r            ! pseudo-random number
    real(8) :: xs_n         ! elastic resonance cross section
    real(8) :: xs_g         ! capture resonance cross section
    real(8) :: capture      ! temporary capture cross section
    real(8) :: xs_f         ! fission resonance cross section
    real(8) :: xs_x         ! competitive resonance cross section
    real(8) :: inelast      ! temporary competitive cross section
    real(8) :: avg_urr_n_xs ! infinite-dilute n xs
    real(8) :: avg_urr_f_xs ! infinite-dilute f xs
    real(8) :: avg_urr_g_xs ! infinite-dilute g xs
    real(8) :: avg_urr_x_xs ! infinite-dilute x xs
    real(8) :: favg         ! average cross section interpolation factor
    real(8) :: xsTlow       ! energy-interpolated xs at lower temperature
    real(8) :: xsTup        ! energy-interpolated xs at upper temperature
    logical :: same_nuc ! do we know the xs for this nuclide at this energy?

    tope => isotopes(i_so)
    nuc  => nuclides(i_nuc)

    tope % E = E
    tope % T = T

    ! determine energy table
    i_E = 1
    do
      if (E < tope % E_tabs(i_E + 1)) exit
      i_E = i_E + 1
    end do

    ! tabulated unresolved resonance parameters interpolation factor
    fE = interp_factor(E, tope % E_tabs(i_E), tope % E_tabs(i_E+1), tope % INT)

    if (tope % nT_tabs == 1) then
      i_Tlow = 1
      i_Tup  = 1
      fT = ZERO
    else
      if (T < tope % T_tabs(1)) then
        call fatal_error('Encountered temperature below probability tables.')
        i_Tlow = 1
      else if (T > tope % T_tabs(tope % nT_tabs)) then
        call fatal_error('Encountered temperature above probability tables.')
        i_Tlow = tope % nT_tabs - 1
      else
        i_Tlow = binary_search(tope % T_tabs, tope % nT_tabs, T)
      end if
      i_Tup  = i_Tlow + 1
      fT = interp_factor(T, tope % T_tabs(i_Tlow), tope % T_tabs(i_Tup), INT_T)
    end if
    
    ! if we're dealing with a nuclide that we've previously encountered at
    ! this energy but a different temperature, use the original random number
    ! to preserve correlation of temperature in probability tables
    same_nuc = .false.
    do i = 1, nuc % nuc_list % size()
      if (E /= ZERO .and. E / 1.0E6_8&
           == micro_xs(nuc % nuc_list % get_item(i)) % last_E) then
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

    i_low = ceiling(r * tope % n_bands)
    i_up  = i_low
    ! elastic xs from probability bands
    xsTlow = interpolator(fE, &
         tope % prob_tables(i_E, i_Tlow) % n(i_low) % xs_mean, &
         tope % prob_tables(i_E + 1, i_Tlow) % n(i_up) % xs_mean, tope % INT)
    xsTup = interpolator(fE, &
         tope % prob_tables(i_E, i_Tup) % n(i_low) % xs_mean, &
         tope % prob_tables(i_E + 1, i_Tup) % n(i_up) % xs_mean, tope % INT)
    xs_n = interpolator(fT, xsTlow, xsTup, INT_T)

    ! fission xs from probability bands
    if ((tope % INT == LINEAR_LINEAR .and. INT_T == LINEAR_LINEAR) .or. &
         (tope % INT == LINEAR_LINEAR .and. INT_T == STATISTICAL) .or. &
         (tope % prob_tables(i_E, i_Tlow) % f(i_low) % xs_mean > ZERO .and. &
         tope % prob_tables(i_E+1, i_Tlow) % f(i_up) % xs_mean > ZERO)) then
      xsTlow = interpolator(fE, &
           tope % prob_tables(i_E, i_Tlow) % f(i_low) % xs_mean, &
           tope % prob_tables(i_E + 1, i_Tlow) % f(i_up) % xs_mean, tope % INT)
      xsTup = interpolator(fE, &
           tope % prob_tables(i_E, i_Tup) % f(i_low) % xs_mean, &
           tope % prob_tables(i_E + 1, i_Tup) % f(i_up) % xs_mean, tope % INT)
      xs_f = interpolator(fT, xsTlow, xsTup, INT_T)
    else
      xs_f = ZERO
    end if

    ! capture xs from probability bands
    if ((tope % INT == LINEAR_LINEAR .and. INT_T == LINEAR_LINEAR) .or. &
         (tope % INT == LINEAR_LINEAR .and. INT_T == STATISTICAL) .or. &
         (tope % prob_tables(i_E, i_Tlow) % g(i_low) % xs_mean > ZERO .and. &
         tope % prob_tables(i_E+1, i_Tlow) % g(i_up) % xs_mean > ZERO)) then
      xsTlow = interpolator(fE, &
           tope % prob_tables(i_E, i_Tlow) % g(i_low) % xs_mean, &
           tope % prob_tables(i_E + 1, i_Tlow) % g(i_up) % xs_mean, tope % INT)
      xsTup = interpolator(fE, &
           tope % prob_tables(i_E, i_Tup) % g(i_low) % xs_mean, &
           tope % prob_tables(i_E + 1, i_Tup) % g(i_up) % xs_mean, tope % INT)
      xs_g = interpolator(fT, xsTlow, xsTup, INT_T)
    else
      xs_g = ZERO
    end if

    ! competitive xs from probability bands
    if ((tope % INT == LINEAR_LINEAR .and. INT_T == LINEAR_LINEAR) .or. &
         (tope % INT == LINEAR_LINEAR .and. INT_T == STATISTICAL) .or. &
         (tope % prob_tables(i_E, i_Tlow) % x(i_low) % xs_mean > ZERO .and. &
         tope % prob_tables(i_E+1, i_Tlow) % x(i_up) % xs_mean > ZERO)) then
      xsTlow = interpolator(fE, &
           tope % prob_tables(i_E, i_Tlow) % x(i_low) % xs_mean, &
           tope % prob_tables(i_E + 1, i_Tlow) % x(i_up) % xs_mean, tope % INT)
      xsTup = interpolator(fE, &
           tope % prob_tables(i_E, i_Tup) % x(i_low) % xs_mean, &
           tope % prob_tables(i_E + 1, i_Tup) % x(i_up) % xs_mean, tope % INT)
      xs_x = interpolator(fT, xsTlow, xsTup, INT_T)
    else
      xs_x = ZERO
    end if

    inelast = micro_xs(i_nuc) % total - micro_xs(i_nuc) % elastic&
         - micro_xs(i_nuc) % absorption
    if (tope % LSSF == 0) then
! TODO: consider moving addition of File 3 background to prob band calculation
      if (competitive) inelast = inelast + xs_x
      capture = micro_xs(i_nuc) % absorption - micro_xs(i_nuc) % fission + xs_g
      micro_xs(i_nuc) % elastic = micro_xs(i_nuc) % elastic + xs_n
      micro_xs(i_nuc) % fission = micro_xs(i_nuc) % fission + xs_f
! TODO: consider moving multiplication by File 3 background to prob band calculation
    else if (tope % LSSF == 1) then
      ! multipy the self-shielding factors by the infinite-dilute xs
      ! averaged cross sections interpolation factor
      iavg = binary_search(tope % Eavg, tope % nEavg, tope % E)
      favg = interp_factor(tope % E, &
           tope % Eavg(iavg), tope % Eavg(iavg + 1), tope % INT)

      ! interpolate averaged, infinite-dilute URR cross sections
      call interp_avg_urr_xs(favg, i_so, iavg, &
           avg_urr_n_xs, avg_urr_f_xs, avg_urr_g_xs, avg_urr_x_xs)

      ! competitive xs
      ! infinite-dilute treatment of competitive xs
      if (avg_urr_x_xs > ZERO&
           .and. tope % E <= (ONE + ENDF_PRECISION) * tope % E_ex2&
           .and. competitive)&
           ! self-shielded treatment of competitive inelastic xs
           inelast = xs_x / avg_urr_x_xs * inelast

      ! elastic scattering xs
      micro_xs(i_nuc) % elastic = xs_n / avg_urr_n_xs &
           * micro_xs(i_nuc) % elastic

      ! set negative elastic xs and competitive xs to zero
      if (micro_xs(i_nuc) % elastic < ZERO) micro_xs(i_nuc) % elastic = ZERO
      if (inelast < ZERO) inelast = ZERO

      ! radiative capture xs
      ! background capture cross section
      capture = micro_xs(i_nuc) % absorption - micro_xs(i_nuc) % fission
      if (avg_urr_g_xs > ZERO) capture = xs_g / avg_urr_g_xs * capture

      ! fission xs
      if (avg_urr_f_xs > ZERO) micro_xs(i_nuc) % fission = xs_f / avg_urr_f_xs&
           * micro_xs(i_nuc) % fission

    end if

    micro_xs(i_nuc) % absorption = capture + micro_xs(i_nuc) % fission
    micro_xs(i_nuc) % total = micro_xs(i_nuc) % elastic &
         + micro_xs(i_nuc) % absorption + inelast

    ! Determine nu-fission cross section
    if (nuc % fissionable) micro_xs(i_nuc) % nu_fission&
         = nu_total(nuc, E / 1.0E6_8) * micro_xs(i_nuc) % fission

    tope % E_last = E

    nullify(tope)
    nullify(nuc)

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
    type(Isotope), pointer :: tope ! isotope pointer
    integer :: iso   ! isotope index
    integer :: n_res ! number of resonances to include for a given l-wave
    integer :: i_l   ! orbital angular momentum quantum number index
    integer :: i_J   ! total angular momentum quantum number index

    tope => isotopes(iso)

    if (tope % E == tope % E_last) then
      ! Energy hasn't changed since last realization, so use the same one
      if (this % i_res == 0) then
        this % E_lam = tope % local_realization(i_l) % J(i_J) % res(1) % E_lam
      else
        this % E_lam&
             = tope % local_realization(i_l) % J(i_J) % res(this%i_res) % E_lam
      end if

    else
      ! sample a level spacing from the Wigner distribution
      this % D_lJ = wigner_surmise(tope % D)

      if (this % i_res == 0) then
        ! set lowest-energy resonance for this ladder well below the energy grid
        ! point such that the ladder spans a sufficient energy range
        n_res = n_res_contrib(tope % L)
        this % E_lam = (tope % E - n_res/2 * tope % D) &
             + (ONE - TWO * prn()) * this % D_lJ
      
      else
        ! add subsequent resonance energies at the sampled spacing above the
        ! last resonance
        this % E_lam = this % E_lam + this % D_lJ
        tope % local_realization(i_l) % J(i_J) % res(this % i_res) % E_lam&
             = this % E_lam

      end if
    end if

    nullify(tope)

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

  function wigner_surmise(D_avg) result(D_samp)

    real(8) :: D_avg  ! mean level spacing
    real(8) :: D_samp ! sampled level spacing

    ! sample a level spacing by directly inverting the Wigner distribution CDF
    D_samp = D_avg * sqrt(-FOUR * log(prn()) / PI)

  end function wigner_surmise

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! CHANNEL_WIDTH samples the channel partial widths at E_lambda when generating
! a full resonance ensemble or at E_n when generating localized parameters for
! an on-the-fly cross section calculation
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine channel_width(this, iso, i_l, i_J)

    class(Resonance), intent(inout) :: this ! pseudo-resonance object
    type(Isotope), pointer :: tope ! isotope object pointer
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
      if (this % i_res == 0) return
      ! Energy hasn't changed since last realization, so use the same one
      this % Gam_n&
           = tope % local_realization(i_l) % J(i_J) % res(this % i_res) % GN
      this % Gam_f&
           = tope % local_realization(i_l) % J(i_J) % res(this % i_res) % GF
      this % Gam_g&
           = tope % local_realization(i_l) % J(i_J) % res(this % i_res) % GG
      this % Gam_x&
           = tope % local_realization(i_l) % J(i_J) % res(this % i_res) % GX
      this % Gam_t&
           = tope % local_realization(i_l) % J(i_J) % res(this % i_res) % GT

    else
      ! sample indices into table of equiprobable chi-squared function values
      i_tabn = ceiling(prn() * 100.0_8)
      i_tabf = ceiling(prn() * 100.0_8)
      i_tabg = ceiling(prn() * 100.0_8)
      i_tabx = ceiling(prn() * 100.0_8)

      ! use the sampled tabulated chi-squared values to calculate sample widths
      ! neutron width
      if (tope % AMUN > 0) then
        ! compute factors needed to go from the mean reduced width that is
        ! provided by ENDF for elastic scattering to a partial width
        ! (use absolute value energies when handling  bound levels which have
        ! negative resonance energies - this is ENDF-6 convention, not theory)
        if (represent_params == E_NEUTRON) then
          rho = tope % k_n * tope % ac(tope % i_urr)
          nu  = tope % P_l_n / rho
          this % Gam_n = tope % GN0 * sqrt(abs(tope % E)) * nu &
               * chi2(i_tabn, tope % AMUN)
        else if (represent_params == E_RESONANCE) then
          rho = tope % k_lam * tope % ac(tope % i_urr)
          nu  = tope % P_l_lam / rho
          this % Gam_n = tope % GN0 * sqrt(abs(this % E_lam)) * nu &
               * chi2(i_tabn, tope % AMUN)
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
        tope % local_realization(i_l) % J(i_J) % res(this % i_res) % GN&
             = this % Gam_n
        tope % local_realization(i_l) % J(i_J) % res(this % i_res) % GF&
             = this % Gam_f
        tope % local_realization(i_l) % J(i_J) % res(this % i_res) % GG&
             = this % Gam_g
        tope % local_realization(i_l) % J(i_J) % res(this % i_res) % GX&
             = this % Gam_x
        tope % local_realization(i_l) % J(i_J) % res(this % i_res) % GT&
             = this % Gam_t

      end if
    end if

    nullify(tope)

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
    type(Isotope), pointer :: tope ! nuclide pointer
    integer :: iso ! isotope index
    real(8) :: k_n_x   ! center-of-mass neutron wavenumber at E_n - QI
    real(8) :: k_lam_x ! center-of-mass neutron wavenumber at |E_lam - QI|
    real(8) :: E_shift ! shifted resonance energy in the lab system
    real(8) :: theta   ! total width / Doppler width
    real(8) :: x       ! derived variable
    real(8) :: Gam_t_n ! sampled energy-dependent total width at E_n
    real(8) :: Gam_n_n ! sampled energy-dependent neutron width at E_n
    real(8) :: Gam_x_n ! sampled energy-dependent competitive width at E_n
    real(8) :: sig_lam ! peak resonance cross section
    real(8) :: sig_lam_Gam_t_n_psi ! derived variable
    real(8) :: S_l_lam ! resonance energy shift factor at E_lam

    tope => isotopes(iso)

    ! if URR parameters have resonance energy dependence
    if (represent_params == E_RESONANCE) then

      tope % k_lam = wavenumber(tope % AWR, abs(this % E_lam))
      tope % P_l_lam = penetration(tope % L,&
           tope % k_lam * tope % ac(tope % i_urr))
      S_l_lam = shift(tope % L, tope % k_lam * tope % ac(tope % i_urr))
      E_shift = this % E_lam &
           + this % Gam_n * (S_l_lam - tope % S_l_n) &
           / (TWO * tope % P_l_lam)

      Gam_n_n = this % Gam_n * tope % P_l_n / tope % P_l_lam

      if (tope % E > (ONE + ENDF_PRECISION) * tope % E_ex2) then
        ! two competitive rxns possible, can't calculate an energy-dependent
        ! width because it depends on the two (unprovided) rxn partial widths
        Gam_x_n = this % Gam_x
      else if (tope % E >= tope % E_ex1) then
        ! can compute an energy-dependent width for the one competitive reaction
        k_n_x = wavenumber(tope % AWR, tope % E - tope % E_ex1)
        k_lam_x = wavenumber(tope % AWR, abs(this % E_lam - tope % E_ex1))
        Gam_x_n = this % Gam_x &
             * penetration(tope % L, k_n_x   * tope % ac(tope % i_urr)) &
             / penetration(tope % L, k_lam_x * tope % ac(tope % i_urr))
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

    sig_lam = FOUR * PI / (tope % k_n * tope % k_n) * tope % g_J &
         * Gam_n_n / Gam_t_n

    ! this particular form comes from the NJOY2012 manual
    this % dxs_n = sig_lam * &
         ((cos(TWO * tope % phi_l_n) &
         - (ONE - Gam_n_n / Gam_t_n)) * psi(tope % T, theta, x) &
         + sin(TWO * tope % phi_l_n) * chi(tope % T, theta, x))

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

    if (Gam_x_n > ZERO .and. tope % E >= tope % E_ex1) then
      this % dxs_x = sig_lam_Gam_t_n_psi * Gam_x_n
    else
      this % dxs_x = ZERO
    end if

    this % dxs_t = this % dxs_n + this % dxs_g + this % dxs_f + this % dxs_x

    nullify(tope)

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
    type(Isotope), pointer :: tope ! nuclide pointer
    integer :: iso ! isotope index
    real(8) :: k_n_x   ! center-of-mass neutron wavenumber at E_n - QI
    real(8) :: k_lam_x ! center-of-mass neutron wavenumber at |E_lam - QI|
    real(8) :: E_shift ! shifted resonance energy in the lab system
    real(8) :: theta   ! total width / Doppler width
    real(8) :: x       ! derived variable
    real(8) :: Gam_t_n ! sampled energy-dependent total width at E_n
    real(8) :: Gam_n_n ! sampled energy-dependent neutron width at E_n
    real(8) :: Gam_x_n ! sampled energy-dependent competitive width at E_n
    real(8) :: sig_lam ! peak resonance cross section
    real(8) :: sig_lam_Gam_t_n_psi ! derived variable
    real(8) :: S_l_lam ! resonance energy shift factor at E_lam

    tope => isotopes(iso)

    if (tope % J == tope % SPI(tope % i_urr)) then
      call fatal_error('Computing MLBW elastic scattering cross section&
           & for resonance with total orbital angular momentum quantum&
           & number, J, equal to the nuclear spin, I')
    end if

    ! if URR parameters have resonance energy dependence
    if (represent_params == E_RESONANCE) then

      tope % k_lam = wavenumber(tope % AWR, abs(this % E_lam))
      tope % P_l_lam = penetration(tope % L,&
           tope % k_lam * tope % ac(tope % i_urr))
      S_l_lam = shift(tope % L, tope % k_lam * tope % ac(tope % i_urr))
      E_shift = this % E_lam &
           + this % Gam_n * (S_l_lam - tope % S_l_n) &
           / (TWO * tope % P_l_lam)

      Gam_n_n = this % Gam_n * tope % P_l_n / tope % P_l_lam

      if (tope % E > (ONE + ENDF_PRECISION) * tope % E_ex2) then
        ! two competitive rxns possible, can't calculate an energy-dependent
        ! width because it depends on the two (unprovided) rxn partial widths
        Gam_x_n = this % Gam_x
      else if (tope % E >= tope % E_ex1) then
        ! can compute an energy-dependent width for the one competitive reaction
        k_n_x = wavenumber(tope % AWR, tope % E - tope % E_ex1)
        k_lam_x = wavenumber(tope % AWR, abs(this % E_lam - tope % E_ex1))
        Gam_x_n = this % Gam_x &
             * penetration(tope % L, k_n_x   * tope % ac(tope % i_urr)) &
             / penetration(tope % L, k_lam_x * tope % ac(tope % i_urr))
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

    sig_lam = FOUR * PI / (tope % k_n * tope % k_n) * tope % g_J &
         * Gam_n_n / Gam_t_n

    ! this particular form comes from the NJOY2012 manual
    this % dxs_n = sig_lam *&
         ((cos(TWO * tope % phi_l_n)&
         - (ONE - Gam_n_n / Gam_t_n)&
         + HALF * G_func(iso, E_shift, Gam_n_n, Gam_t_n, this % i_res)&
         / Gam_n_n) * psi(tope % T, theta, x)&
         + (sin(TWO * tope % phi_l_n)&
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

    if (Gam_x_n > ZERO .and. tope % E >= tope % E_ex1) then
      this % dxs_x = sig_lam_Gam_t_n_psi * Gam_x_n
    else
      this % dxs_x = ZERO
    end if

    this % dxs_t = this % dxs_n + this % dxs_g + this % dxs_f + this % dxs_x

    nullify(tope)

  end subroutine mlbw_xs

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! G_FUNC calculates a value for the G-function appearing in the NJOY-2012 form
! of the MLBW resonance formalae
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  function G_func(iso, E_res, G_n, G_t, i_res) result(G_val)

    type(Isotope), pointer :: tope ! isotope pointer
    integer :: iso   ! isotope index
    integer :: i_l   ! orbital angular momentum quantum number index
    integer :: i_J   ! total angular momentum quantum number index
    integer :: i_r   ! spin sequence resonance index
    integer :: i_res ! index of current nearest resonance
    real(8) :: G_val ! G-function value
    real(8) :: E_res ! energy of the resonance contributing to the xs at E_n
    real(8) :: G_n   ! neutron width of the resonance at E_res
    real(8) :: G_t   ! total width of the resonance at E_res
    real(8) :: k_n_x   ! center-of-mass neutron wavenumber at E_n - QI
    real(8) :: k_lam   ! center-of-mass neutron wavenumber at E_lam
    real(8) :: k_lam_x ! center-of-mass neutron wavenumber at |E_lam - QI|
    real(8) :: E_shift ! shifted resonance energy in the lab system
    real(8) :: Gam_t_n ! sampled energy-dependent total width at E_n
    real(8) :: Gam_n_n ! sampled energy-dependent neutron width at E_n
    real(8) :: Gam_x_n ! sampled energy-dependent competitive width at E_n
    real(8) :: P_l_lam ! penetration at resonance energy
    real(8) :: S_l_lam ! resonance energy shift factor

    tope => isotopes(iso)

    G_val = ZERO

    i_l = tope % L + 1

    TOTAL_ANG_MOM_LOOP: do i_J = 1, tope % NJS(i_l)
      RESONANCE_LOOP: do i_r = 1, n_res_contrib(tope % L)
        if (i_r == i_res .and. tope % J == tope % AJ(i_l) % dim1(i_J)) cycle

        ! if URR parameters have resonance energy dependence
        if (represent_params == E_RESONANCE) then
          ! absolute value of energy in order to handle bound levels which have
          ! negative resonance energies
          k_lam = wavenumber(tope % AWR,&
               abs(tope % local_realization(i_l) % J(i_J) % res(i_r) % E_lam))
          P_l_lam = penetration(tope % L, k_lam * tope % ac(tope % i_urr))
          S_l_lam = shift(tope % L, k_lam * tope % ac(tope % i_urr))

          E_shift = tope % local_realization(i_l) % J(i_J) % res(i_r) % E_lam&
               + tope % local_realization(i_l) % J(i_J) % res(i_r) % GN&
               * (S_l_lam - tope % S_l_n) / (TWO * P_l_lam)

          Gam_n_n = tope % local_realization(i_l) % J(i_J) % res(i_r) % GN&
               * tope % P_l_n / P_l_lam

          if (tope % E > (ONE + ENDF_PRECISION) * tope % E_ex2) then
            ! two competitive reactions possible;
            ! can't calculate an energy-dependent width because it depends on
            ! the two (unprovided) reaction partial widths
            Gam_x_n = tope % local_realization(i_l) % J(i_J) % res(i_r) % GX
          else if (tope % E >= tope % E_ex1) then
            ! compute an energy-dependent width for the one competitive reaction
            k_n_x = wavenumber(tope % AWR, tope % E - tope % E_ex1)
            k_lam_x = wavenumber(tope % AWR,&
                 abs(tope % local_realization(i_l) % J(i_J) % res(i_r) % E_lam&
                 - tope % E_ex1))
            Gam_x_n = tope % local_realization(i_l) % J(i_J) % res(i_r) % GX&
                 * penetration(tope % L, k_n_x   * tope % ac(tope % i_urr))&
                 / penetration(tope % L, k_lam_x * tope % ac(tope % i_urr))
          else
            Gam_x_n = ZERO
          end if

        else

          ! assume all URR parameters already have neutron energy dependence
          E_shift = tope % local_realization(i_l) % J(i_J) % res(i_r) % E_lam
          Gam_n_n = tope % local_realization(i_l) % J(i_J) % res(i_r) % GN
          Gam_x_n = tope % local_realization(i_l) % J(i_J) % res(i_r) % GX

        end if

        Gam_t_n = tope % local_realization(i_l) % J(i_J) % res(i_r) % GT&
             - tope % local_realization(i_l) % J(i_J) % res(i_r) % GN&
             - tope % local_realization(i_l) % J(i_J) % res(i_r) % GX&
             + Gam_n_n&
             + Gam_x_n

        G_val = G_val + G_n * Gam_n_n * (G_t + Gam_t_n)&
             / ((E_res - E_shift) * (E_res - E_shift)&
             + (G_t + Gam_t_n) * (G_t + Gam_t_n) / 4.0_8)

      end do RESONANCE_LOOP
    end do TOTAL_ANG_MOM_LOOP

    nullify(tope)

  end function G_func

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! H_FUNC calculates a value for the H-function appearing in the NJOY-2012 form
! of the MLBW resonance formalae
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  function H_func(iso, E_res, G_n, G_t, i_res) result(H_val)

    type(Isotope), pointer :: tope ! isotope pointer
    integer :: iso   ! isotope index
    integer :: i_l   ! orbital angular momentum quantum number index
    integer :: i_J   ! total angular momentum quantum number index
    integer :: i_r   ! spin sequence resonance index
    integer :: i_res ! index of current nearest resonance
    real(8) :: H_val ! G-function value
    real(8) :: E_res ! energy of the resonance contributing to the xs at E_n
    real(8) :: G_n   ! neutron width of the resonance at E_res
    real(8) :: G_t   ! total width of the resonance at E_res
    real(8) :: k_n_x   ! center-of-mass neutron wavenumber at E_n - QI
    real(8) :: k_lam   ! center-of-mass neutron wavenumber at E_lam
    real(8) :: k_lam_x ! center-of-mass neutron wavenumber at |E_lam - QI|
    real(8) :: E_shift ! shifted resonance energy in the lab system
    real(8) :: Gam_t_n ! sampled energy-dependent total width at E_n
    real(8) :: Gam_n_n ! sampled energy-dependent neutron width at E_n
    real(8) :: Gam_x_n ! sampled energy-dependent competitive width at E_n
    real(8) :: P_l_lam ! penetration at resonance energy
    real(8) :: S_l_lam ! resonance energy shift factor

    tope => isotopes(iso)

    H_val = ZERO

    i_l = tope % L + 1

    TOTAL_ANG_MOM_LOOP: do i_J = 1, tope % NJS(i_l)
      RESONANCE_LOOP: do i_r = 1, n_res_contrib(tope % L)
        if (i_r == i_res .and. tope % J == tope % AJ(i_l) % dim1(i_J)) cycle

        ! if URR parameters have resonance energy dependence
        if (represent_params == E_RESONANCE) then
          ! absolute value of energy in order to handle bound levels which have
          ! negative resonance energies
          k_lam = wavenumber(tope % AWR,&
               abs(tope % local_realization(i_l) % J(i_J) % res(i_r) % E_lam))
          P_l_lam = penetration(tope % L, k_lam * tope % ac(tope % i_urr))
          S_l_lam = shift(tope % L, k_lam * tope % ac(tope % i_urr))

          E_shift = tope % local_realization(i_l) % J(i_J) % res(i_r) % E_lam&
               + tope % local_realization(i_l) % J(i_J) % res(i_r) % GN&
               * (S_l_lam - tope % S_l_n) / (TWO * P_l_lam)

          Gam_n_n = tope % local_realization(i_l) % J(i_J) % res(i_r) % GN&
               * tope % P_l_n / P_l_lam

          if (tope % E > (ONE + ENDF_PRECISION) * tope % E_ex2) then
            ! two competitive reactions possible;
            ! can't calculate an energy-dependent width because it depends on
            ! the two (unprovided) reaction partial widths
            Gam_x_n = tope % local_realization(i_l) % J(i_J) % res(i_r) % GX
          else if (tope % E >= tope % E_ex1) then
            ! compute an energy-dependent width for the one competitive reaction
            k_n_x = wavenumber(tope % AWR, tope % E - tope % E_ex1)
            k_lam_x = wavenumber(tope % AWR,&
                 abs(tope % local_realization(i_l) % J(i_J) % res(i_r) % E_lam&
                 - tope % E_ex1))
            Gam_x_n = tope % local_realization(i_l) % J(i_J) % res(i_r) % GX&
                 * penetration(tope % L, k_n_x   * tope % ac(tope % i_urr))&
                 / penetration(tope % L, k_lam_x * tope % ac(tope % i_urr))
          else
            Gam_x_n = ZERO
          end if

        else

          ! assume all URR parameters already have neutron energy dependence
          E_shift = tope % local_realization(i_l) % J(i_J) % res(i_r) % E_lam
          Gam_n_n = tope % local_realization(i_l) % J(i_J) % res(i_r) % GN
          Gam_x_n = tope % local_realization(i_l) % J(i_J) % res(i_r) % GX

        end if

        Gam_t_n = tope % local_realization(i_l) % J(i_J) % res(i_r) % GT&
             - tope % local_realization(i_l) % J(i_J) % res(i_r) % GN&
             - tope % local_realization(i_l) % J(i_J) % res(i_r) % GX&
             + Gam_n_n&
             + Gam_x_n

        H_val = H_val + G_n * Gam_n_n * (E_res - E_shift)&
             / ((E_res - E_shift) * (E_res - E_shift)&
             + (G_t + Gam_t_n) * (G_t + Gam_t_n) / 4.0_8)

      end do RESONANCE_LOOP
    end do TOTAL_ANG_MOM_LOOP

    nullify(tope)

  end function H_func

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! PENETRATION calculates hard sphere penetrability factors
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  function penetration(L, rho) result(P)

    integer, intent(in) :: L ! current orbital quantum number
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
           / (225.0_8 + rho2 * (45.0_8 + rho2 * (6.0_8 + rho2)))

    case(4)
      P = rho * rho2 * rho2 * rho2 * rho2 &
           / (11025.0_8 + rho2 * (1575.0_8 + rho2 &
           * (135.0_8 + rho2 * (10.0_8 + rho2))))

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
           / (105.0_8 - 45.0_8 * rho2 + rho4))

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
           / (225.0_8 + rho2 * (45.0_8 + rho2 * (6.0_8 + rho2)))

    case(4)
      S = -(44100.0_8 + rho2 * (4725.0_8 + rho2 * (270.0_8 + 10.0_8 * rho2)))&
           / (11025.0_8 + rho2 * (1575.0_8 + rho2 * (135.0_8 &
           + rho2 * (10.0_8 + rho2))))

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

      case (MIT_W)
        ! call S.G. Johnson's Faddeeva evaluation
        relerr = 1.0e-6
        w_val = faddeeva_w(cmplx(theta * x * HALF, theta * HALF, 8), relerr)
        psi_val = SQRT_PI * HALF * theta&
             * real(real(w_val, 8), 8)

      case (QUICK_W)
        ! QUICKW Faddeeva evaluation from Argonne (also used in NJOY)
        psi_val = SQRT_PI * HALF * theta&
             * real(real(quickw(cmplx(theta * x * HALF, theta * HALF, 8)),8),8)

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

      case (MIT_W)
        ! S.G. Johnson's Faddeeva evaluation
        relerr = 1.0e-6
        w_val = faddeeva_w(cmplx(theta * x * HALF, theta * HALF, 8), relerr)
        chi_val = SQRT_PI * HALF * theta&
             * real(aimag(w_val), 8)

      case (QUICK_W)
        ! QUICKW Faddeeva evaluation from Argonne (also used in NJOY)
        chi_val = SQRT_PI * HALF * theta&
             * real(aimag(quickw(cmplx(theta * x * HALF, theta * HALF, 8))), 8)

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
    type(Isotope), pointer :: tope ! nuclide pointer
    integer :: iso ! isotope index
    integer :: i_l ! orbital quantum number index
    real(8) :: sig_pot ! potential scattering cross section
    real(8) :: phi_l_n ! hard-sphere phase shift

    ! set nuclide variables
    tope => isotopes(iso)

    ! compute potential scattering xs by adding contribution from each l-wave
    sig_pot = ZERO
    do i_l = 0, tope % NLS(tope % i_urr) - 1
      phi_l_n = phase_shift(i_l, tope % k_n * tope % AP(tope % i_urr)) 
      sig_pot = sig_pot&
           + FOUR * PI / (tope % k_n * tope % k_n) * (TWO * dble(i_l) + ONE)&
           * sin(phi_l_n) * sin(phi_l_n)
    end do

    ! add the potential scattering xs to this xs
    this % xs = this % xs + sig_pot

    nullify(tope)

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

    ! accumulate infinite-dilute xs values for this realization
    this % prob_tables(i_E, i_T) % avg_t % xs_tmp_sum &
         = this % prob_tables(i_E, i_T) % avg_t % xs_tmp_sum &
         + this % prob_tables(i_E, i_T) % avg_t % xs
    this % prob_tables(i_E, i_T) % avg_n % xs_tmp_sum &
         = this % prob_tables(i_E, i_T) % avg_n % xs_tmp_sum &
         + this % prob_tables(i_E, i_T) % avg_n % xs
    this % prob_tables(i_E, i_T) % avg_g % xs_tmp_sum &
         = this % prob_tables(i_E, i_T) % avg_g % xs_tmp_sum &
         + this % prob_tables(i_E, i_T) % avg_g % xs
    this % prob_tables(i_E, i_T) % avg_f % xs_tmp_sum &
         = this % prob_tables(i_E, i_T) % avg_f % xs_tmp_sum &
         + this % prob_tables(i_E, i_T) % avg_f % xs
    this % prob_tables(i_E, i_T) % avg_x % xs_tmp_sum &
         = this % prob_tables(i_E, i_T) % avg_x % xs_tmp_sum &
         + this % prob_tables(i_E, i_T) % avg_x % xs

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

    do i_E = 1, this % nE_tabs
      do i_T = 1, this % nT_tabs
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

    if (this % E /= this % E_tabs(i_E)) return
    do i_T = 1, this % nT_tabs
      ! accumulate single-batch infinite-dilute means
      this % prob_tables(i_E, i_T) % avg_t % xs_sum &
           = this % prob_tables(i_E, i_T) % avg_t % xs_sum &
           + this % prob_tables(i_E, i_T) % avg_t % xs_tmp_sum &
           / dble(histories_avg_urr)
      this % prob_tables(i_E, i_T) % avg_n % xs_sum &
           = this % prob_tables(i_E, i_T) % avg_n % xs_sum &
           + this % prob_tables(i_E, i_T) % avg_n % xs_tmp_sum &
           / dble(histories_avg_urr)
      this % prob_tables(i_E, i_T) % avg_g % xs_sum &
           = this % prob_tables(i_E, i_T) % avg_g % xs_sum &
           + this % prob_tables(i_E, i_T) % avg_g % xs_tmp_sum &
           / dble(histories_avg_urr)
      this % prob_tables(i_E, i_T) % avg_f % xs_sum &
           = this % prob_tables(i_E, i_T) % avg_f % xs_sum &
           + this % prob_tables(i_E, i_T) % avg_f % xs_tmp_sum &
           / dble(histories_avg_urr)
      this % prob_tables(i_E, i_T) % avg_x % xs_sum &
           = this % prob_tables(i_E, i_T) % avg_x % xs_sum &
           + this % prob_tables(i_E, i_T) % avg_x % xs_tmp_sum &
           / dble(histories_avg_urr)

      ! accumulate squared single-batch infinite-dilute means
      this % prob_tables(i_E, i_T) % avg_t % xs_sum2 &
           = this % prob_tables(i_E, i_T) % avg_t % xs_sum2 &
           + (this % prob_tables(i_E, i_T) % avg_t % xs_tmp_sum &
           / dble(histories_avg_urr))&
           * (this % prob_tables(i_E, i_T) % avg_t % xs_tmp_sum &
           / dble(histories_avg_urr))
      this % prob_tables(i_E, i_T) % avg_n % xs_sum2 &
           = this % prob_tables(i_E, i_T) % avg_n % xs_sum2 &
           + (this % prob_tables(i_E, i_T) % avg_n % xs_tmp_sum &
           / dble(histories_avg_urr))&
           * (this % prob_tables(i_E, i_T) % avg_n % xs_tmp_sum &
           / dble(histories_avg_urr))
      this % prob_tables(i_E, i_T) % avg_g % xs_sum2 &
           = this % prob_tables(i_E, i_T) % avg_g % xs_sum2 &
           + (this % prob_tables(i_E, i_T) % avg_g % xs_tmp_sum &
           / dble(histories_avg_urr))&
           * (this % prob_tables(i_E, i_T) % avg_g % xs_tmp_sum &
           / dble(histories_avg_urr))
      this % prob_tables(i_E, i_T) % avg_f % xs_sum2 &
           = this % prob_tables(i_E, i_T) % avg_f % xs_sum2 &
           + (this % prob_tables(i_E, i_T) % avg_f % xs_tmp_sum &
           / dble(histories_avg_urr))&
           * (this % prob_tables(i_E, i_T) % avg_f % xs_tmp_sum &
           / dble(histories_avg_urr))
      this % prob_tables(i_E, i_T) % avg_x % xs_sum2 &
           = this % prob_tables(i_E, i_T) % avg_x % xs_sum2 &
           + (this % prob_tables(i_E, i_T) % avg_x % xs_tmp_sum &
           / dble(histories_avg_urr))&
           * (this % prob_tables(i_E, i_T) % avg_x % xs_tmp_sum &
           / dble(histories_avg_urr))

    end do

  end subroutine accum_batch

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! CALC_STATS computes batch-based means and standard errors of those means
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine calc_stats(this, i_bat_int)

    class(Isotope), intent(inout), target :: this ! cross section object
    type(ProbabilityTable), pointer :: ptable ! prob table pointer
    integer :: i_E       ! tabulated URR energy index
    integer :: i_T       ! temperature index
    integer :: i_bat_int ! current batch index
    real(8) :: i_bat ! current batch index for calcs

    i_bat = dble(i_bat_int)

    do i_E = 1, this % nE_tabs
      if (this % E /= this % E_tabs(i_E)) cycle
      do i_T = 1, this % nT_tabs
        ptable => this % prob_tables(i_E, i_T)
        ptable % avg_t % xs_mean = ptable % avg_t % xs_sum / i_bat
        ptable % avg_n % xs_mean = ptable % avg_n % xs_sum / i_bat
        ptable % avg_g % xs_mean = ptable % avg_g % xs_sum / i_bat
        ptable % avg_f % xs_mean = ptable % avg_f % xs_sum / i_bat
        ptable % avg_x % xs_mean = ptable % avg_x % xs_sum / i_bat

        if (ptable % avg_t % xs_mean > ZERO .and. i_bat_int > 1) then
          ! compute standard errors of mean xs values
          ptable % avg_t % xs_sem&
               = sqrt((ONE / (i_bat - ONE)) * (ptable % avg_t % xs_sum2 / i_bat&
               - (ptable % avg_t % xs_mean) * (ptable % avg_t % xs_mean)))

          ! compute relative uncertainties in xs values
          ptable % avg_t % rel_unc &
               = ptable % avg_t % xs_sem / ptable % avg_t % xs_mean
        end if
        if (ptable % avg_n % xs_mean > ZERO .and. i_bat_int > 1) then
          ! compute standard errors of mean xs values
          ptable % avg_n % xs_sem &
               = sqrt((ONE / (i_bat - ONE)) * (ptable % avg_n % xs_sum2 / i_bat&
               - (ptable % avg_n % xs_mean) * (ptable % avg_n % xs_mean)))

          ! compute relative uncertainties in xs values
          ptable % avg_n % rel_unc &
               = ptable % avg_n % xs_sem / ptable % avg_n % xs_mean
        end if
        if (ptable % avg_g % xs_mean > ZERO .and. i_bat_int > 1) then
          ! compute standard errors of mean xs values
          ptable % avg_g % xs_sem &
               = sqrt((ONE / (i_bat - ONE)) * (ptable % avg_g % xs_sum2 / i_bat&
               - (ptable % avg_g % xs_mean) * (ptable % avg_g % xs_mean)))

          ! compute relative uncertainties in xs values
          ptable % avg_g % rel_unc &
               = ptable % avg_g % xs_sem / ptable % avg_g % xs_mean
        end if
        if (ptable % avg_f % xs_mean > ZERO .and. i_bat_int > 1&
             .and. (ptable % avg_f % xs_sum2 / i_bat &
             - ptable % avg_f % xs_mean * ptable % avg_f % xs_mean) > ZERO) then
          ! compute standard errors of mean xs values
          ptable % avg_f % xs_sem &
               = sqrt((ONE / (i_bat - ONE)) * (ptable % avg_f % xs_sum2 / i_bat&
               - (ptable % avg_f % xs_mean) * (ptable % avg_f % xs_mean)))

          ! compute relative uncertainties in xs values
          ptable % avg_f % rel_unc &
               = ptable % avg_f % xs_sem / ptable % avg_f % xs_mean
        end if
        if (ptable % avg_x % xs_mean > ZERO &
             .and. i_bat_int > 1 .and. competitive) then
          ! compute standard errors of mean xs values
          ptable % avg_x % xs_sem &
               = sqrt((ONE / (i_bat - ONE)) * (ptable % avg_x % xs_sum2 / i_bat&
               - (ptable % avg_x % xs_mean) * (ptable % avg_x % xs_mean)))

          ! compute relative uncertainties in xs values
          ptable % avg_x % rel_unc &
               = ptable % avg_x % xs_sem / ptable % avg_x % xs_mean
        end if
      end do
    end do

    nullify(ptable)

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
    real(8) :: xi      ! prn for statistical interpolation
    real(8) :: prob    ! probability for statistical interpolation

    select case(scheme)
    case(HISTOGRAM)
      factor = ONE

    case(LINEAR_LINEAR)
      factor = (val - val_low) / (val_up - val_low)

    case(LINEAR_LOG)
      factor = (val - val_low) / (val_up - val_low)

    case(LOG_LINEAR)
      factor = log(val / val_low) / log(val_up / val_low)

    case(LOG_LOG)
      factor = log(val / val_low) / log(val_up / val_low)

    case(SQRT_LINEAR)
      factor = (sqrt(val) - sqrt(val_low)) / (sqrt(val_up) - sqrt(val_low))

    case(SQRT_LOG)
      factor = (sqrt(val) - sqrt(val_low)) / (sqrt(val_up) - sqrt(val_low))

    case(STATISTICAL)
      xi = prn()
      prob = (val - val_low) / (val_up - val_low)
      if (xi > prob) then
        factor = ZERO
      else
        factor = ONE
      end if

    case default
      call fatal_error('Interpolation scheme not recognized')

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
    case(HISTOGRAM)
      val = val_low

    case(LINEAR_LINEAR)
      val = val_low + factor * (val_up - val_low)

    case(LINEAR_LOG)
      val = val_low * exp(factor * log(val_up / val_low))

    case(LOG_LINEAR)
      val = val_low + factor * (val_up - val_low)

    case(LOG_LOG)
      val = val_low * exp(factor * log(val_up / val_low))

    case(SQRT_LINEAR)
      val = val_low + factor * (val_up - val_low)

    case(SQRT_LOG)
      val = val_low * exp(factor * log(val_up / val_low))

    case(STATISTICAL)
      val = val_low + factor * (val_up - val_low)

    case default
      call fatal_error('Interpolation scheme not recognized')

    end select

  end function interpolator

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! E_LAST_RRR determines the energy of the highest-energy resolved resonance
! region resonance for a given (l,J) spin sequence
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  function E_last_rrr(iso, l_val, J_val) result(E_val)

    type(Isotope), pointer :: tope ! nuclide pointer
    integer :: iso   ! isotope index
    integer :: i_res ! RRR resonance index for a given l
    integer :: l_val ! orbital quantum number
    real(8) :: J_val ! total angular momentum quantum number
    real(8) :: E_val ! highest-energy RRR resonance energy for (l,J)

    tope => isotopes(iso)

    select case (tope % LRF(tope % i_urr - 1))
    
    case (SLBW)
      do i_res = size(tope % bw_resonances(l_val + 1) % res(:)), 1, -1
        if (tope % bw_resonances(l_val + 1) % res(i_res) % AJ== J_val&
             .and. tope % bw_resonances(l_val + 1) % res(i_res) % E_lam&
             < tope % EL(tope % i_urr)) then
          E_val = tope % bw_resonances(l_val + 1) % res(i_res) % E_lam
          exit
        end if
      end do
    
    case (MLBW)
      do i_res = size(tope % bw_resonances(l_val + 1) % res(:)), 1, -1
        if (tope % bw_resonances(l_val + 1) % res(i_res) % AJ == J_val&
             .and. tope % bw_resonances(l_val + 1) % res(i_res) % E_lam&
             < tope % EL(tope % i_urr)) then
          E_val = tope % bw_resonances(l_val + 1) % res(i_res) % E_lam
          exit
        end if
      end do
    
    case (REICH_MOORE)
      do i_res = size(tope % rm_resonances(l_val + 1) % res(:)), 1, -1
        if (tope % rm_resonances(l_val + 1) % res(i_res) % AJ == J_val&
             .and. tope % rm_resonances(l_val + 1) % res(i_res) % E_lam&
             < tope % EL(tope % i_urr)) then
          E_val = tope % rm_resonances(l_val + 1) % res(i_res) % E_lam
          exit
        end if
      end do
    
    case default
      call fatal_error('Unrecognized/unsupported RRR formalism')
    
    end select

    nullify(tope)

  end function E_last_rrr

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! RRR_RES finds the index of the RRR resonance which we need to add the
! contribution of to a URR xs
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  function rrr_res(iso, n_rrr_res, l_val, J_val) result(i_res)

    type(Isotope), pointer :: tope ! nuclide pointer
    integer :: iso       ! isotope index
    integer :: n_rrr_res ! how many RRR resonances to go back
    integer :: cnt_res   ! how many RRR resonances have we gone back
    integer :: i_res     ! index of the RRR resonance cnt_res resonances back
    integer :: l_val     ! orbital quantum number
    real(8) :: J_val ! total angular momentum quantum number

    tope => isotopes(iso)

    cnt_res   = 0

    select case (tope % LRF(tope % i_urr - 1))

    case (SLBW)
      do i_res = size(tope % bw_resonances(l_val + 1) % res(:)), 1, -1
        if (tope % bw_resonances(l_val + 1) % res(i_res) % AJ == J_val &
             .and. tope % bw_resonances(l_val + 1) % res(i_res) % E_lam&
             < tope % EL(tope % i_urr)) cnt_res = cnt_res + 1
        if (cnt_res == n_rrr_res) exit
      end do

    case (MLBW)
      do i_res = size(tope % bw_resonances(l_val + 1) % res(:)), 1, -1
        if (tope % bw_resonances(l_val + 1) % res(i_res) % AJ == J_val &
             .and. tope % bw_resonances(l_val + 1) % res(i_res) % E_lam&
             < tope % EL(tope % i_urr)) cnt_res = cnt_res + 1
        if (cnt_res == n_rrr_res) exit
      end do
    
    case (REICH_MOORE)
      do i_res = size(tope % rm_resonances(l_val + 1) % res(:)), 1, -1
        if (tope % rm_resonances(l_val + 1) % res(i_res) % AJ == J_val &
             .and. tope % rm_resonances(l_val + 1) % res(i_res) % E_lam&
             < tope % EL(tope % i_urr)) cnt_res = cnt_res + 1
        if (cnt_res == n_rrr_res) exit
      end do
    
    case default
      call fatal_error('Unrecognized/unsupported RRR formalism')
    
    end select

    ! if there aren't enough contributing RRR resonances
    if (cnt_res < n_rrr_res) i_res = 0

    nullify(tope)

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
    type(CrossSection) :: x ! competitive xs object

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

    do i_E = 1, this % nE_tabs
      do i_T = 1, this % nT_tabs
        do i_b = 1, this % n_bands
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
    type(CrossSection) :: x ! competitive xs object

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

  subroutine add_parameters(res, iso, i_ens, i_l, i_J)

    type(Isotope), pointer :: tope ! isotope object pointer
    type(Resonance) :: res ! resonance object
    integer :: iso   ! isotope index
    integer :: i_ens ! resonance ensemble index
    integer :: i_l   ! orbital angular momentum index
    integer :: i_J   ! total angular momentum quantum number index

    tope => isotopes(iso)

    tope % urr_resonances_tmp(i_l, i_ens) % J(i_J) % res % E_lam = res % E_lam
    tope % urr_resonances_tmp(i_l, i_ens) % J(i_J) % res % AJ = tope % J
    tope % urr_resonances_tmp(i_l, i_ens) % J(i_J) % res % GN = res % Gam_n
    tope % urr_resonances_tmp(i_l, i_ens) % J(i_J) % res % GG = res % Gam_g
    tope % urr_resonances_tmp(i_l, i_ens) % J(i_J) % res % GF = res % Gam_f
    tope % urr_resonances_tmp(i_l, i_ens) % J(i_J) % res % GX = res % Gam_x
    tope % urr_resonances_tmp(i_l, i_ens) % J(i_J) % res % GT = res % Gam_t

    ! point to next resonance
    allocate(tope % urr_resonances_tmp(i_l, i_ens) % J(i_J) % res % next)
    tope % urr_resonances_tmp(i_l, i_ens) % J(i_J) % res&
         => tope % urr_resonances_tmp(i_l, i_ens) % J(i_J) % res % next

    nullify(tope)

  end subroutine add_parameters

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! SET_PARAMETERS sets the URR resonance parameters for a single resonance
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine set_parameters(res, iso, i_res, i_l, i_J, i_ER)

    type(Isotope), pointer :: tope ! isotope object pointer
    type(Resonance) :: res ! resonance object
    integer :: iso   ! isotope index
    integer :: i_res ! resonance counter
    integer :: i_l   ! orbital quantum number index
    integer :: i_J   ! total angular momentum quantum number
    integer :: i_ER  ! resonance energy region index

    tope => isotopes(iso)

    select case (tope % LRU(i_ER))

    ! resolved parameters
    case (1)
      select case (tope % LRF(i_ER))
      case (SLBW)
        if (i_res > size(tope % bw_resonances(i_l) % res(:))) then
          i_res = size(tope % bw_resonances(i_l) % res(:))
          res % E_lam = res % E_lam&
               + (tope % bw_resonances(i_l) % res(i_res) % E_lam&
               -  tope % bw_resonances(i_l) % res(i_res - 1) % E_lam)
        else
          res % E_lam = tope % bw_resonances(i_l) % res(i_res) % E_lam
        end if

        res % Gam_n = tope % bw_resonances(i_l) % res(i_res) % GN
        res % Gam_g = tope % bw_resonances(i_l) % res(i_res) % GG
        res % Gam_f = tope % bw_resonances(i_l) % res(i_res) % GF
        res % Gam_t = tope % bw_resonances(i_l) % res(i_res) % GT
        res % Gam_x = res % Gam_t - res % Gam_n - res % Gam_g - res % Gam_f

      case (MLBW)
        if (i_res > size(tope % bw_resonances(i_l) % res(:))) then
          i_res = size(tope % bw_resonances(i_l) % res(:))
          res % E_lam = res % E_lam&
               + (tope % bw_resonances(i_l) % res(i_res) % E_lam&
               -  tope % bw_resonances(i_l) % res(i_res - 1) % E_lam)
        else
          res % E_lam = tope % bw_resonances(i_l) % res(i_res) % E_lam
        end if
        
        res % Gam_n = tope % bw_resonances(i_l) % res(i_res) % GN
        res % Gam_g = tope % bw_resonances(i_l) % res(i_res) % GG
        res % Gam_f = tope % bw_resonances(i_l) % res(i_res) % GF
        res % Gam_t = tope % bw_resonances(i_l) % res(i_res) % GT
        res % Gam_x = res % Gam_t - res % Gam_n - res % Gam_g - res % Gam_f

      case (REICH_MOORE)
        if (i_res > size(tope % rm_resonances(i_l) % res(:))) then
          i_res = size(tope % rm_resonances(i_l) % res(:))
          res % E_lam = res % E_lam&
               + (tope % rm_resonances(i_l) % res(i_res) % E_lam&
               -  tope % rm_resonances(i_l) % res(i_res - 1) % E_lam)
        else
          res % E_lam = tope % rm_resonances(i_l) % res(i_res) % E_lam
        end if

        res % Gam_n = tope % rm_resonances(i_l) % res(i_res) % GN
        res % Gam_g = tope % rm_resonances(i_l) % res(i_res) % GG
        res % Gam_f = tope % rm_resonances(i_l) % res(i_res) % GFA&
             + tope % rm_resonances(i_l) % res(i_res) % GFB
        res % Gam_x = ZERO
        res % Gam_t = res % Gam_n + res % Gam_g + res % Gam_f + res % Gam_x

      case default
        call fatal_error('Unrecognized resolved resonance region formalism')

      end select

    ! unresolved parameters
    case (2)
      res % E_lam&
           = tope % urr_resonances(i_l, i_real) % J(i_J) % res(i_res) % E_lam
      res % Gam_n&
           = tope % urr_resonances(i_l, i_real) % J(i_J) % res(i_res) % GN
      res % Gam_g&
           = tope % urr_resonances(i_l, i_real) % J(i_J) % res(i_res) % GG
      res % Gam_f&
           = tope % urr_resonances(i_l, i_real) % J(i_J) % res(i_res) % GF
      res % Gam_x&
           = tope % urr_resonances(i_l, i_real) % J(i_J) % res(i_res) % GX
      res % Gam_t&
           = tope % urr_resonances(i_l, i_real) % J(i_J) % res(i_res) % GT

    case default
      call fatal_error('Only 1 and 2 are supported ENDF-6 LRU values')

    end select

    tope % local_realization(i_l) % J(i_J) % res(res % i_res) % E_lam&
         = res % E_lam
    tope % local_realization(i_l) % J(i_J) % res(res % i_res) % GN = res % Gam_n
    tope % local_realization(i_l) % J(i_J) % res(res % i_res) % GG = res % Gam_g
    tope % local_realization(i_l) % J(i_J) % res(res % i_res) % GF = res % Gam_f
    tope % local_realization(i_l) % J(i_J) % res(res % i_res) % GX = res % Gam_x
    tope % local_realization(i_l) % J(i_J) % res(res % i_res) % GT = res % Gam_t

    nullify(tope)

  end subroutine set_parameters

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! SET_MEAN_PARAMETERS sets the URR mean resonance parameters at an energy
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine set_mean_parameters(iso, E_res, i_l, i_J)

    type(Isotope), pointer :: tope ! isotope object pointer
    integer :: iso ! isotope index
    integer :: i_E ! tabulated URR parameters energy index
    integer :: i_l ! orbital quantum number index
    integer :: i_J ! total angular momentum quantum number
    real(8) :: E_res ! current resonance (lab) energy (e.g. E_lam)
    real(8) :: fendf ! ENDF6 URR parameters energy interpolation factor

    tope => isotopes(iso)

    ! compute interpolation factor
    if (E_res < tope % ES(1)) then
      i_E = 1
    else if (E_res > tope % ES(tope % NE - 1)) then
      i_E = tope % NE - 1
    else
      i_E = binary_search(tope % ES, tope % NE, E_res)
    end if
    fendf = interp_factor(E_res, tope % ES(i_E), tope % ES(i_E + 1),&
      tope % INT)

    ! set current mean unresolved resonance parameters
    tope % D   = interpolator(fendf, &
         tope % D_mean(i_l) % dim2(i_J) % dim1(i_E), &
         tope % D_mean(i_l) % dim2(i_J) % dim1(i_E + 1), tope % INT)
    tope % GN0 = interpolator(fendf, &
         tope % GN0_mean(i_l) % dim2(i_J) % dim1(i_E), &
         tope % GN0_mean(i_l) % dim2(i_J) % dim1(i_E + 1), tope % INT)
    tope % GG  = interpolator(fendf, &
         tope % GG_mean(i_l) % dim2(i_J) % dim1(i_E), &
         tope % GG_mean(i_l) % dim2(i_J) % dim1(i_E + 1), tope % INT)
    if (tope % INT == LINEAR_LINEAR .or. &
         tope % GF_mean(i_l) % dim2(i_J) % dim1(i_E) > ZERO) then
      tope % GF  = interpolator(fendf, &
           tope % GF_mean(i_l) % dim2(i_J) % dim1(i_E), &
           tope % GF_mean(i_l) % dim2(i_J) % dim1(i_E + 1), tope % INT)
    else
      tope % GF = ZERO
    end if
    if (tope % INT == LINEAR_LINEAR .or. &
         tope % GX_mean(i_l) % dim2(i_J) % dim1(i_E) > ZERO) then
      tope % GX  = interpolator(fendf, &
           tope % GX_mean(i_l) % dim2(i_J) % dim1(i_E), &
           tope % GX_mean(i_l) % dim2(i_J) % dim1(i_E + 1), tope % INT)
    else
      tope % GX = ZERO
    end if

    nullify(tope)

  end subroutine set_mean_parameters

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! INTERP_AVG_URR_XS interpolates the averaged, infinite-dilute URR cross
! sections computed via Monte Carlo from mean resonance parameters (i.e. not
! the evaluator-supplied File 3 background cross sections)
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine interp_avg_urr_xs(f, iso, iavg, n_xs, f_xs, g_xs, x_xs)

    type(Isotope), pointer :: tope ! isotope object pointer
    integer :: iso
    integer :: iavg ! averaged cross section index
    real(8) :: f
    real(8) :: n_xs
    real(8) :: f_xs
    real(8) :: g_xs
    real(8) :: x_xs

    tope => isotopes(iso)

    ! infinite-dilute elastic scattering
    if (tope % avg_urr_n(iavg) > ZERO&
         .and. tope % avg_urr_n(iavg + 1) > ZERO) then
      n_xs = interpolator(f, &
           tope % avg_urr_n(iavg), tope % avg_urr_n(iavg + 1), tope % INT)
    else
      call fatal_error('Non-positive (n,n) infinite-dilute cross section')
    end if

    ! infinite-dilute fission
    if (tope % INT == LINEAR_LINEAR&
         .or. (tope % avg_urr_f(iavg) > ZERO&
         .and. tope % avg_urr_f(iavg + 1) > ZERO)) then
      f_xs = interpolator(f, &
           tope % avg_urr_f(iavg), tope % avg_urr_f(iavg + 1), tope % INT)
    else
      f_xs = ZERO
    end if

    ! infinite-dilute capture
    if (tope % INT == LINEAR_LINEAR&
         .or. (tope % avg_urr_g(iavg) > ZERO&
         .and. tope % avg_urr_g(iavg + 1) > ZERO)) then
      g_xs = interpolator(f, &
           tope % avg_urr_g(iavg), tope % avg_urr_g(iavg + 1), tope % INT)
    else
      g_xs = ZERO
    end if

    ! infinite-dilute competitive reaction xs
    if (tope % INT == LINEAR_LINEAR&
         .or. (tope % avg_urr_x(iavg) > ZERO&
         .and. tope % avg_urr_x(iavg + 1) > ZERO)) then
      x_xs = interpolator(f, &
           tope % avg_urr_x(iavg), tope % avg_urr_x(iavg + 1), tope % INT)
    else
      x_xs = ZERO
    end if

    nullify(tope)

  end subroutine interp_avg_urr_xs

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! MF3_SELF_SHIELDING interpolates the evaluator-supplied ENDF-6 File 3
! background and either adds a resonance component (LSSF = 0) or interpolates
! the computed average URR xs value and then multiplies the File 3 xs by
! a resonance component divided by the average, which is a self-shielding factor
! (LSSF = 1)
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine mf3_self_shielding(iso, n, g, f, x, t)
    type(Isotope), pointer :: tope ! isotope object pointer
    integer, intent(in) :: iso   ! isotope index
    type(CrossSection), intent(inout) :: n ! combined elastic cross section
    type(CrossSection), intent(inout) :: g ! combined capture cross section
    type(CrossSection), intent(inout) :: f ! combined fission cross section
    type(CrossSection), intent(inout) :: x ! combined competitive cross section
    type(CrossSection), intent(inout) :: t ! combined total cross section
    integer :: i_mf3 ! index in a File 3 grid
    integer :: i_avg ! index in the averaged xs grid
    real(8) :: avg_n ! averaged elastic xs
    real(8) :: avg_g ! averaged capture xs
    real(8) :: avg_f ! averaged fission xs
    real(8) :: avg_x ! averaged competitive xs
    real(8) :: mf3_n ! MF3 elastic xs
    real(8) :: mf3_g ! MF3 capture xs
    real(8) :: mf3_f ! MF3 fission xs
    real(8) :: mf3_x ! MF3 competitive xs
    real(8) :: f_mf3 ! interpolation factor in a File 3 grid
    real(8) :: f_avg ! interpolation factor in the averaged xs grid

    tope => isotopes(iso)

    ! elastic scattering xs
    if (tope % E < tope % MF3_n_e(1)) then
      call fatal_error('Energy is below File 3 elastic energy grid')
    else if (tope % E > tope % MF3_n_e(size(tope % MF3_n_e))) then
      call fatal_error('Energy is above File 3 elastic energy grid')
    else
      i_mf3 = binary_search(tope % MF3_n_e, size(tope % MF3_n_e), tope % E)
      if (tope % INT == LINEAR_LINEAR&
           .or. (tope % MF3_n(i_mf3) > ZERO&
           .and. tope % MF3_n(i_mf3 + 1) > ZERO)) then
        f_mf3 = interp_factor(tope % E, tope % MF3_n_e(i_mf3),&
             tope % MF3_n_e(i_mf3 + 1), tope % INT)
        mf3_n = interpolator(f_mf3, tope % MF3_n(i_mf3),&
             tope % MF3_n(i_mf3 + 1), tope % INT)
      else
        mf3_n = ZERO
      end if
    end if
    if (mf3_n < ZERO) mf3_n = ZERO

    ! competitive reaction xs
    if (tope % E < tope % MF3_x_e(1)) then
      mf3_x = ZERO
    else if (tope % E > tope % MF3_x_e(size(tope % MF3_x_e))) then
      call fatal_error('Energy is above File 3 competitive energy grid')
    else
      i_mf3 = binary_search(tope % MF3_x_e, size(tope % MF3_x_e), tope % E)
      if (tope % INT == LINEAR_LINEAR&
           .or. (tope % MF3_x(i_mf3) > ZERO&
           .and. tope % MF3_x(i_mf3 + 1) > ZERO)) then
        f_mf3 = interp_factor(tope % E, tope % MF3_x_e(i_mf3),&
             tope % MF3_x_e(i_mf3 + 1), tope % INT)
        mf3_x = interpolator(f_mf3, tope % MF3_x(i_mf3),&
             tope % MF3_x(i_mf3 + 1), tope % INT)
      else
        mf3_x = ZERO
      end if
    end if
    if (mf3_x < ZERO) mf3_x = ZERO

    ! capture xs
    if (tope % E < tope % MF3_g_e(1)) then
      call fatal_error('Energy is below File 3 capture energy grid')
    else if (tope % E > tope % MF3_g_e(size(tope % MF3_g_e))) then
      call fatal_error('Energy is above File 3 capture energy grid')
    else
      i_mf3 = binary_search(tope % MF3_g_e, size(tope % MF3_g_e), tope % E)
      if (tope % INT == LINEAR_LINEAR&
           .or. (tope % MF3_g(i_mf3) > ZERO&
           .and. tope % MF3_g(i_mf3 + 1) > ZERO)) then
        f_mf3 = interp_factor(tope % E, tope % MF3_g_e(i_mf3),&
             tope % MF3_g_e(i_mf3 + 1), tope % INT)
        mf3_g = interpolator(f_mf3, tope % MF3_g(i_mf3),&
             tope % MF3_g(i_mf3 + 1), tope % INT)
      else
        mf3_g = ZERO
      end if
    end if
    if (mf3_g < ZERO) mf3_g = ZERO

    ! fission xs
    if (.not. (allocated(tope % MF3_f_e))) then
      mf3_f = ZERO
    else
      if (tope % E < tope % MF3_f_e(1)) then
        mf3_f = ZERO
      else if (tope % E > tope % MF3_f_e(size(tope % MF3_f_e))) then
        call fatal_error('Energy is above File 3 fission energy grid')
      else
        i_mf3 = binary_search(tope % MF3_f_e, size(tope % MF3_f_e), tope % E)
        if (tope % INT == LINEAR_LINEAR&
             .or. (tope % MF3_f(i_mf3) > ZERO&
             .and. tope % MF3_f(i_mf3 + 1) > ZERO)) then
          f_mf3 = interp_factor(tope % E, tope % MF3_f_e(i_mf3),&
               tope % MF3_f_e(i_mf3 + 1), tope % INT)
          mf3_f = interpolator(f_mf3, tope % MF3_f(i_mf3),&
               tope % MF3_f(i_mf3 + 1), tope % INT)
        else
          mf3_f = ZERO
        end if
      end if
    end if
    if (mf3_f < ZERO) mf3_f = ZERO

    ! add resonance component to File 3 xs or use it as a self-shielding factor
    if (tope % LSSF == 0) then
      n % xs = mf3_n + n % xs
      g % xs = mf3_g + g % xs
      f % xs = mf3_f + f % xs
      if (competitive) x % xs = mf3_x + x % xs

    else if (tope % LSSF == 1) then

      if (tope % E < tope % Eavg(1)) then
        i_avg = 1
      else if (tope % E > tope % Eavg(tope % nEavg - 1)) then
        i_avg = tope % nEavg - 1
      else
        i_avg = binary_search(tope % Eavg, tope % nEavg, tope % E)
      end if

      f_avg = interp_factor(tope % E,&
           tope % Eavg(i_avg), tope % Eavg(i_avg + 1), tope % INT)

      call interp_avg_urr_xs(f_avg, iso, i_avg, avg_n, avg_f, avg_g, avg_x)

      ! competitive xs
      if (avg_x > ZERO&
           .and. tope % E <= (ONE + ENDF_PRECISION) * tope % E_ex2&
           .and. competitive)&
           x % xs = mf3_x * x % xs / avg_x

      ! elastic scattering xs
      n % xs = mf3_n * n % xs / avg_n

      ! set negative elastic xs and competitive xs to zero
      if (n % xs < ZERO) n % xs = ZERO
      if (x % xs < ZERO) x % xs = ZERO

      ! radiative capture xs
      if (avg_g > ZERO) g % xs = mf3_g * g % xs / avg_g

      ! fission xs
      if (avg_f > ZERO) f % xs = mf3_f * f % xs / avg_f

    else
      call fatal_error('ENDF-6 LSSF not allowed - must be 0 or 1.')

    end if

    t % xs = n % xs + g % xs + f % xs + x % xs

    nullify(tope)

  end subroutine mf3_self_shielding

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! ADD_MF3_BACKGROUND adds the resonance xs component to the evaluator-supplied
! background
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine add_mf3_background(iso, i_nuc, n, g, f, x)

    type(Isotope), pointer :: tope ! isotope object pointer
    type(CrossSection) :: n ! elastic scattering xs object
    type(CrossSection) :: g ! radiative capture xs object
    type(CrossSection) :: f ! fission xs object
    type(CrossSection) :: x ! competitive xs object
    integer :: iso    ! isotope index
    integer :: i_nuc  ! nuclide index
    integer :: i_grid ! background energy grid index
    real(8) :: fmf3         ! File 3 interpolation factor
    real(8) :: capture_xs   ! radiative capture xs
    real(8) :: inelastic_xs ! competitive xs

    tope => isotopes(iso)

    ! elastic scattering xs
    if (tope % E < tope % MF3_n_e(1)) then
      call fatal_error('Energy is below File 3 elastic energy grid')
    else if (tope % E > tope % MF3_n_e(size(tope % MF3_n_e))) then
      call fatal_error('Energy is above File 3 elastic energy grid')
    else
      i_grid = binary_search(tope % MF3_n_e, size(tope % MF3_n_e),tope % E)
      if (tope % INT == LINEAR_LINEAR &
           .or. (tope % MF3_n(i_grid) > ZERO &
           .and. tope % MF3_n(i_grid + 1) > ZERO)) then
        fmf3 = interp_factor(tope % E, tope % MF3_n_e(i_grid), &
             tope % MF3_n_e(i_grid + 1), tope % INT)
        micro_xs(i_nuc) % elastic = interpolator(fmf3, tope % MF3_n(i_grid),&
             tope % MF3_n(i_grid + 1), tope % INT) + n % xs
      else
        micro_xs(i_nuc) % elastic = n % xs
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
      if (tope % INT == LINEAR_LINEAR &
           .or. (tope % MF3_x(i_grid) > ZERO &
           .and. tope % MF3_x(i_grid + 1) > ZERO)) then
        fmf3 = interp_factor(tope % E, tope % MF3_x_e(i_grid), &
             tope % MF3_x_e(i_grid + 1), tope % INT)
        inelastic_xs = interpolator(fmf3, tope % MF3_x(i_grid), &
             tope % MF3_x(i_grid + 1), tope % INT)
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
      if (tope % INT == LINEAR_LINEAR &
           .or. (tope % MF3_g(i_grid) > ZERO &
           .and. tope % MF3_g(i_grid + 1) > ZERO)) then
        fmf3 = interp_factor(tope % E, tope % MF3_g_e(i_grid), &
             tope % MF3_g_e(i_grid + 1), tope % INT)
        capture_xs = interpolator(fmf3, tope % MF3_g(i_grid), &
             tope % MF3_g(i_grid + 1), tope % INT) + g % xs
      else
        capture_xs = g % xs
      end if
    end if

    ! fission xs
    if (.not. (allocated(tope % MF3_f_e))) then
      micro_xs(i_nuc) % fission = f % xs
    else
      if (tope % E < tope % MF3_f_e(1)) then
        micro_xs(i_nuc) % fission = ZERO
      else if (tope % E > tope % MF3_f_e(size(tope % MF3_f_e))) then
        call fatal_error('Energy is above File 3 fission energy grid')
      else
        i_grid = binary_search(tope % MF3_f_e, size(tope % MF3_f_e), tope%E)
        if (tope % INT == LINEAR_LINEAR &
             .or. (tope % MF3_f(i_grid) > ZERO &
             .and. tope % MF3_f(i_grid + 1) > ZERO)) then
          fmf3 = interp_factor(tope % E, tope % MF3_f_e(i_grid), &
               tope % MF3_f_e(i_grid + 1), tope % INT)
          micro_xs(i_nuc) % fission = interpolator(fmf3, tope%MF3_f(i_grid),&
               tope % MF3_f(i_grid + 1), tope % INT) + f % xs
        else
          micro_xs(i_nuc) % fission = f % xs
        end if
      end if
    end if

    micro_xs(i_nuc) % absorption = micro_xs(i_nuc) % fission + capture_xs
    micro_xs(i_nuc) % total = micro_xs(i_nuc) % elastic &
         + micro_xs(i_nuc) % absorption + inelastic_xs

    nullify(tope)

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
! ALLOC_ENSEMBLE allocates URR resonance ensemble realizations
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine alloc_ensemble(this)

    class(Isotope), intent(inout) :: this ! isotope object
    integer :: i_ens ! realization/ensemble index
    integer :: i_l   ! orbital angular momentum index
    integer :: i_J   ! total angular momentum index

    ! allocate URR resonance realizations
    allocate(this % urr_resonances(this % NLS(this % i_urr), n_reals))

    ! loop over realizations
    do i_ens = 1, n_reals

      ! loop over orbital quantum numbers
      do i_l = 1, this % NLS(this % i_urr)

        allocate(this % urr_resonances(i_l, i_ens) % J(this % NJS(i_l)))

        ! loop over total angular momenta
        do i_J = 1, this % NJS(i_l)

          ! allocate URR resonances
          allocate(this % urr_resonances(i_l, i_ens) % J(i_J) %&
               res(this % n_lam(i_l, i_ens) % dim1(i_J)))

        end do
      end do
    end do

  end subroutine alloc_ensemble

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! DEALLOC_ENSEMBLE deallocates a URR resonance ensemble realization
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine dealloc_ensemble(this)

    class(Isotope), intent(inout) :: this ! isotope object
    integer :: i_ens ! realization/ensemble index
    integer :: i_l   ! orbital angular momentum quantum number index
    integer :: i_J   ! total angular momentum index

    ! loop over realizations
    do i_ens = 1, n_reals

      ! loop over orbital quantum numbers
      do i_l = 1, this % NLS(this % i_urr)

        ! loop over total angular momenta
        do i_J = 1, this % NJS(i_l)
! TODO: deallocate resonances of the proper formalism once MLBW, RM allowed
          ! deallocate resonance parameters
          deallocate(this % urr_resonances(i_l, i_ens) % J(i_J) % res)
      
        end do

        deallocate(this % urr_resonances(i_l, i_ens) % J)

      end do
    end do

    ! deallocate URR resonance realizations
    deallocate(this % urr_resonances)

  end subroutine dealloc_ensemble

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
      this % nE_tabs = this % NE
      allocate(this % E_tabs(this % nE_tabs))
      this % E_tabs(:) = this % ES
    else if (E_spacing == LOGARITHMIC) then
      this % nE_tabs = n_energies
      allocate(this % E_tabs(this % nE_tabs))
      this % E_tabs(:) = (/(this % EL(this % i_urr) &
           * exp(i_E * log(this % EH(this % i_urr) / this % EL(this % i_urr)) &
           / (this % nE_tabs - 1)), i_E = 0, this % nE_tabs - 1)/)
    else if (E_spacing == USER) then
      this % nE_tabs = n_energies
      allocate(this % E_tabs(this % nE_tabs))
      this % E_tabs(:) = E_tabs
    end if

    this % nT_tabs = n_temperatures
    allocate(this % T_tabs(this % nT_tabs))
    this % T_tabs(:) = T_tabs
    this % n_bands = n_bands

    allocate(this % prob_tables(this % nE_tabs, this % nT_tabs))

    do i_E = 1, this % nE_tabs
      do i_T = 1, this % nT_tabs
        allocate(this % prob_tables(i_E, i_T) % t(this % n_bands))
        allocate(this % prob_tables(i_E, i_T) % n(this % n_bands))
        allocate(this % prob_tables(i_E, i_T) % g(this % n_bands))
        allocate(this % prob_tables(i_E, i_T) % f(this % n_bands))
        allocate(this % prob_tables(i_E, i_T) % x(this % n_bands))
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

    do i_E = 1, this % nE_tabs
      do i_T = 1, this % nT_tabs
        deallocate(this % prob_tables(i_E, i_T) % t)
        deallocate(this % prob_tables(i_E, i_T) % n)
        deallocate(this % prob_tables(i_E, i_T) % g)
        deallocate(this % prob_tables(i_E, i_T) % f)
        deallocate(this % prob_tables(i_E, i_T) % x)
      end do
    end do

    deallocate(this % E_tabs)
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

    do i_T = 1, this % nT_tabs
      do i_b = 1, this % n_bands
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
    integer :: i_l   ! orbital angular momentum index
    integer :: i_J   ! total angular momentum index
    integer :: i_ens ! realization index

    ! deallocate mean parameters
    do i_l = 1, this % NLS(this % i_urr)
      do i_J = 1, this % NJS(i_l)
        deallocate(this % D_mean(i_l) % dim2(i_J) % dim1)
        deallocate(this % GN0_mean(i_l) % dim2(i_J) % dim1)
        deallocate(this % GG_mean(i_l) % dim2(i_J) % dim1)
        deallocate(this % GF_mean(i_l) % dim2(i_J) % dim1)
        deallocate(this % GX_mean(i_l) % dim2(i_J) % dim1)
      end do
      if (allocated(this % n_lam)) then
        do i_ens = 1, n_reals
          deallocate(this % n_lam(i_l, i_ens) % dim1)
        end do
      end if
      deallocate(this % D_mean(i_l) % dim2)
      deallocate(this % GN0_mean(i_l) % dim2)
      deallocate(this % GG_mean(i_l) % dim2)
      deallocate(this % GF_mean(i_l) % dim2)
      deallocate(this % GX_mean(i_l) % dim2)
      deallocate(this % AJ(i_l) % dim1)
      deallocate(this % DOFX(i_l) % dim1)
      deallocate(this % DOFN(i_l) % dim1)
      deallocate(this % DOFG(i_l) % dim1)
      deallocate(this % DOFF(i_l) % dim1)
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
    if (allocated(this % urr_resonances)) call this % dealloc_ensemble()

    ! deallocate RRR resonances
    if (allocated(this % bw_resonances)) then
      do i_l = 1, this % NLS(this % i_urr - 1)
        call this % bw_resonances(i_l) % clear()
      end do
      deallocate(this % bw_resonances)
    end if
    if (allocated(this % rm_resonances)) then
      do i_l = 1, this % NLS(this % i_urr - 1)
        call this % rm_resonances(i_l) % clear()
      end do
      deallocate(this % rm_resonances)
    end if

    ! deallocate pointwise data
    if (allocated(this % urr_E)) call this % dealloc_pointwise()

    ! deallocate averaged, infinite-dilute URR cross sections
    if (allocated(this % Eavg)) deallocate(this % Eavg)
    if (allocated(this % avg_urr_t)) deallocate(this % avg_urr_t)
    if (allocated(this % avg_urr_n)) deallocate(this % avg_urr_n)
    if (allocated(this % avg_urr_f)) deallocate(this % avg_urr_f)
    if (allocated(this % avg_urr_g)) deallocate(this % avg_urr_g)
    if (allocated(this % avg_urr_x)) deallocate(this % avg_urr_x)

    ! deallocate ENDF-6 File 3 cross sections
    call this % dealloc_MF3()

    ! deallocate energy range variables and total angular momenta counts
    call this % dealloc_energy_ranges()
    deallocate(this % NJS)

  end subroutine dealloc_isotope

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! ALLOC_LOCAL_REALIZATION allocates an NLS-length vector of NJS(l)-length
! vectors of NRS(l,J)-length vectors of spin group resonances
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine alloc_local_realization(this)
  
    class(Isotope), intent(inout) :: this ! isotope object
    integer :: i_l ! orbital angular momentum index 
    integer :: i_J ! total angular momentum index

    ! allocate NLS-length vec of NJS(l)-length vec of NRS(l,J)-length vec
    allocate(this % local_realization(this % NLS(this % i_urr)))
    
    ! loop over orbital quantum numbers
    do i_l = 1, this % NLS(this % i_urr)

      ! allocate each vector of (l,J) vectors
      allocate(this % local_realization(i_l) % J(this % NJS(i_l)))

      ! loop of total quantum numbers
      do i_J = 1, this % NJS(i_l)

        ! allocate each vector of (l,J) resonances
        allocate(this%local_realization(i_l)%J(i_J)%res(n_res_contrib(i_l-1)))

      end do

    end do

  end subroutine alloc_local_realization

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! DEALLOC_LOCAL_REALIZATION deallocates an NLS-length vector of NJS(l)-length
! vectors of NRS(l,J)-length vectors of spin group resonances
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine dealloc_local_realization(this)
  
    class(Isotope), intent(inout) :: this ! isotope object
    integer :: i_l ! orbital angular momentum index 
    integer :: i_J ! total angular momentum index

    ! loop over orbital quantum numbers
    do i_l = 1, this % NLS(this % i_urr)

      ! loop of total quantum numbers
      do i_J = 1, this % NJS(i_l)

        ! deallocate each vector of (l,J) resonances
        call this % local_realization(i_l) % J(i_J) % clear()

      end do
      
      ! deallocate each vector of (l,J) vectors
      deallocate(this % local_realization(i_l) % J)

    end do

    ! deallocate NLS-length vec of NJS(l)-length vec of NRS(l,J)-length vec
    deallocate(this % local_realization)

  end subroutine dealloc_local_realization

end module xs
