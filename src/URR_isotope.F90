module URR_isotope

  implicit none
  private
  public :: Isotope,&
            isotopes

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! ENDF6FLAGS
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  type ENDF6Flags

     
     
  end type ENDF6Flags

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

end module URR_isotope
