module ace_header

  use constants,     only: MAX_FILE_LEN, ONE, THREE
  use endf_header,   only: Tab1
  use list_header,   only: ListInt
  use vector_header, only: Vector, JaggedArray

  implicit none

!===============================================================================
! DISTANGLE contains data for a tabular secondary angle distribution whether it
! be tabular or 32 equiprobable cosine bins
!===============================================================================

  type DistAngle
    integer              :: n_energy    ! # of incoming energies
    real(8), allocatable :: energy(:)   ! incoming energy grid
    integer, allocatable :: type(:)     ! type of distribution
    integer, allocatable :: location(:) ! location of each table
    real(8), allocatable :: data(:)     ! angular distribution data

    ! Type-Bound procedures
    contains
      procedure :: clear => distangle_clear ! Deallocates DistAngle
  end type DistAngle

!===============================================================================
! DISTENERGY contains data for a secondary energy distribution for all
! scattering laws
!===============================================================================

  type DistEnergy
    integer    :: law                 ! secondary distribution law
    type(Tab1) :: p_valid             ! probability of law validity
    real(8), allocatable :: data(:)   ! energy distribution data

    ! For reactions that may have multiple energy distributions such as (n,2n),
    ! this pointer allows multiple laws to be stored
    type(DistEnergy), pointer :: next => null()

    ! Type-Bound procedures
    contains
      procedure :: clear => distenergy_clear ! Deallocates DistEnergy
  end type DistEnergy

!===============================================================================
! REACTION contains the cross-section and secondary energy and angle
! distributions for a single reaction in a continuous-energy ACE-format table
!===============================================================================

  type Reaction
    integer :: MT                      ! ENDF MT value
    real(8) :: Q_value                 ! Reaction Q value
    integer :: multiplicity            ! Number of secondary particles released
    integer :: threshold               ! Energy grid index of threshold
    logical :: scatter_in_cm           ! scattering system in center-of-mass?
    real(8), allocatable :: sigma(:)   ! Cross section values
    logical :: has_angle_dist          ! Angle distribution present?
    logical :: has_energy_dist         ! Energy distribution present?
    type(DistAngle)           :: adist ! Secondary angular distribution
    type(DistEnergy), pointer :: edist => null() ! Secondary energy distribution

    ! Type-Bound procedures
    contains
      procedure :: clear => reaction_clear ! Deallocates Reaction
  end type Reaction

!===============================================================================
! URRDATA contains probability tables for the unresolved resonance range.
!===============================================================================

  type UrrData
    integer :: n_energy        ! # of incident neutron energies
    integer :: n_prob          ! # of probabilities
    integer :: interp          ! inteprolation (2=lin-lin, 5=log-log)
    integer :: inelastic_flag  ! inelastic competition flag
    integer :: absorption_flag ! other absorption flag
    logical :: multiply_smooth ! multiply by smooth cross section?
    real(8), allocatable :: energy(:)   ! incident energies
    real(8), allocatable :: prob(:,:,:) ! actual probabibility tables

    ! Type-Bound procedures
    contains
      procedure :: clear => urrdata_clear ! Deallocates UrrData
  end type UrrData

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!                                                                               
! URR_RESONANCES is an object containing a vector of URR resonances' information
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  type URRResonances

     real(8), allocatable :: E_lam(:)
     real(8), allocatable :: GN(:)
     real(8), allocatable :: GG(:)
     real(8), allocatable :: GF(:)
     real(8), allocatable :: GX(:)
     real(8), allocatable :: GT(:)

     ! type-bound procedures
     contains

       ! allocate vector of URR resonances (for a J, for a given (i_lam,i_l))
       procedure :: alloc_resonances => resonances_alloc

       ! deallocate vector of URR resonances (for a J, for a given (i_lam,i_l))
       procedure :: dealloc_resonances => resonances_dealloc

  end type URRResonances

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!                                                                               
! REICHMOORERESONANCES is an object containing a vector of Reich-Moore resonance
! data
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  type ReichMooreResonances

     real(8), allocatable :: E_lam(:)
     real(8), allocatable :: AJ(:)
     real(8), allocatable :: GN(:)
     real(8), allocatable :: GG(:)
     real(8), allocatable :: GFA(:)
     real(8), allocatable :: GFB(:)

  end type ReichMooreResonances

!===============================================================================
! NUCLIDE contains all the data for an ACE-format continuous-energy cross
! section. The ACE format (A Compact ENDF format) is used in MCNP and several
! other Monte Carlo codes.
!===============================================================================

  type Nuclide
    character(10) :: name    ! name of nuclide, e.g. 92235.03c
    integer       :: zaid    ! Z and A identifier, e.g. 92235
    integer       :: listing ! index in xs_listings
    real(8)       :: awr     ! weight of nucleus in neutron masses
    real(8)       :: T       ! termperature in K
    real(8)       :: kT      ! temperature in MeV (k*T)

    ! Linked list of indices in nuclides array of instances of this same nuclide
    type(ListInt) :: nuc_list

    ! Energy grid information
    integer :: n_grid                     ! # of nuclide grid points
    integer, allocatable :: grid_index(:) ! pointers to union grid
    real(8), allocatable :: energy(:)     ! energy values corresponding to xs

    ! Microscopic cross sections
    real(8), allocatable :: total(:)      ! total cross section
    real(8), allocatable :: elastic(:)    ! elastic scattering
    real(8), allocatable :: fission(:)    ! fission
    real(8), allocatable :: nu_fission(:) ! neutron production
    real(8), allocatable :: absorption(:) ! absorption (MT > 100)
    real(8), allocatable :: heating(:)    ! heating

    ! Resonance scattering info
    logical              :: resonant = .false. ! resonant scatterer?
    character(10)        :: name_0K = '' ! name of 0K nuclide, e.g. 92235.00c
    character(16)        :: scheme ! target velocity sampling scheme
    integer              :: n_grid_0K ! number of 0K energy grid points
    real(8), allocatable :: energy_0K(:)  ! energy grid for 0K xs
    real(8), allocatable :: elastic_0K(:) ! Microscopic elastic cross section
    real(8), allocatable :: xs_cdf(:) ! CDF of v_rel times cross section
    real(8)              :: E_min ! lower cutoff energy for res scattering
    real(8)              :: E_max ! upper cutoff energy for res scattering

    ! Fission information
    logical :: fissionable         ! nuclide is fissionable?
    logical :: has_partial_fission ! nuclide has partial fission reactions?
    integer :: n_fission           ! # of fission reactions
    integer, allocatable :: index_fission(:) ! indices in reactions

    ! Total fission neutron emission
    integer :: nu_t_type
    real(8), allocatable :: nu_t_data(:)

    ! Prompt fission neutron emission
    integer :: nu_p_type
    real(8), allocatable :: nu_p_data(:)

    ! Delayed fission neutron emission
    integer :: nu_d_type
    integer :: n_precursor ! # of delayed neutron precursors
    real(8), allocatable :: nu_d_data(:)
    real(8), allocatable :: nu_d_precursor_data(:)
    type(DistEnergy), pointer :: nu_d_edist(:) => null()

    ! URR treatment parameters and indices
    logical :: urr_present
    logical :: otf_urr_xs   ! do an on-the-fly URR cross section calculation?
    logical :: avg_urr_xs   ! do an average URR cross section calculation?
    logical :: point_urr_xs ! calculate pointwise URR cross sections?
    integer :: urr_inelastic_index
    integer :: i_urr ! energy range index of unresolved resonance region
    type(UrrData), pointer :: urr_data => null()

    ! ENDF-6 nuclear data
    integer, allocatable :: NLS(:)  ! number of orbital quantum numbers
    integer, allocatable :: NJS(:)  ! number of J values for each l
    integer, allocatable :: LRU(:)  ! resolved (1) or unresolved (2)?
    integer, allocatable :: LRF(:)  ! resonance formalism #
    integer, allocatable :: NRO(:)  ! AP energy-dependence flag
    integer, allocatable :: NAPS(:) ! channel radius handling flag
    real(8), allocatable :: EL(:)  ! lower energy bound of energy region
    real(8), allocatable :: EH(:)  ! upper energy bound of energy region
    real(8), allocatable :: AP(:)  ! scattering radius
    real(8), allocatable :: ac(:)  ! channel radius
    real(8), allocatable :: SPI(:) ! total spin
    real(8), allocatable :: ES(:)  ! URR tabulated data energies
    integer :: MAT  ! nuclide MAT number
    integer :: LRP  ! resonance parameter flag
    integer :: NER  ! number of resonance energy ranges
    integer :: LSSF ! self-shielding factor flag
    integer :: INT  ! interpolation scheme #
! TODO: use LRX if ZERO's aren't given for inelastic width in ENDF when LRX = 0
    integer :: LRX  ! competitive inelastic width flag
    integer :: NE   ! number of URR tabulated data energies

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
    integer :: AMUF ! number of fission channels(degrees of freedom)
    integer :: AMUX ! number of competitive channels (degrees of freedom)

    ! average (infinite-dilute) cross sections values
    real(8), allocatable :: avg_urr_n(:)
    real(8), allocatable :: avg_urr_f(:)
    real(8), allocatable :: avg_urr_g(:)
    real(8), allocatable :: avg_urr_x(:)

    ! URR resonance realization (vector of resonances for a value
    ! of J for a given (i_lam, L))
    type(URRResonances), allocatable :: urr_resonances(:,:)

    ! set of Reich-Moore resonances (vector of resonances for each l
    type(ReichMooreResonances), allocatable :: rm_resonances(:)

    ! pointwise URR cross section data
    real(8), allocatable :: urr_energy_tmp(:)    ! energy grid values
    real(8), allocatable :: urr_elastic_tmp(:)   ! elastic scattering
    real(8), allocatable :: urr_capture_tmp(:)   ! capture
    real(8), allocatable :: urr_fission_tmp(:)   ! fission
    real(8), allocatable :: urr_inelastic_tmp(:) ! first level inelastic scattering
    real(8), allocatable :: urr_total_tmp(:)     ! total
    real(8), allocatable :: urr_energy(:)        ! energy grid values
    real(8), allocatable :: urr_elastic(:)       ! elastic scattering
    real(8), allocatable :: urr_capture(:)       ! capture
    real(8), allocatable :: urr_fission(:)       ! fission
    real(8), allocatable :: urr_inelastic(:)     ! first level inelastic scattering
    real(8), allocatable :: urr_total(:)         ! total

    ! pointwise URR cross section parameters
    integer :: n_urr_resonances = 1000000 ! max URR resonances for a given (l,J)
    integer :: n_urr_gridpoints = 100000000 ! max URR energy-cross section gridpoints
    real(8) :: urr_dE = 0.1_8 ! difference between URR energy grid points [eV]

    ! Reactions
    integer :: n_reaction ! # of reactions
    type(Reaction), pointer :: reactions(:) => null()

    ! Type-Bound procedures
    contains

      ! allocate resonance energy range variables
      procedure :: alloc_energy_range => energy_range_alloc

      ! deallocate resonance energy range variables
      procedure :: dealloc_energy_range => energy_range_dealloc

      ! allocate average (infinite-dilute) cross sections
      procedure :: alloc_avg_urr => avg_urr_alloc

      ! deallocate average (infinite-dilute) cross sections
      procedure :: dealloc_avg_urr => avg_urr_dealloc

      ! allocate URR resonance ensemble realization
      procedure :: alloc_ensemble => ensemble_alloc

      ! deallocate URR resonance ensemble realization
      procedure :: dealloc_ensemble => ensemble_dealloc

      ! allocate temporary pointwise URR cross sections
      procedure :: alloc_pointwise_tmp => pointwise_tmp_alloc

      ! deallocate temporary pointwise URR cross sections
      procedure :: dealloc_pointwise_tmp => pointwise_tmp_dealloc

      ! allocate pointwise URR cross sections
      procedure :: alloc_pointwise => pointwise_alloc

      ! deallocate pointwise URR cross sections
      procedure :: dealloc_pointwise => pointwise_dealloc

      ! Deallocates Nuclide
      procedure :: clear => nuclide_clear

      ! pre-process data that needs it
      procedure :: pre_process => process_pre

      ! set channel radius
      procedure :: channel_radius => radius_channel

  end type Nuclide

!===============================================================================
! NUCLIDE0K temporarily contains all 0K cross section data and other parameters
! needed to treat resonance scattering before transferring them to NUCLIDE
!===============================================================================

  type Nuclide0K

    character(10) :: nuclide            ! name of nuclide, e.g. U-238
    character(16) :: scheme = 'ares'    ! target velocity sampling scheme
    character(10) :: name               ! name of nuclide, e.g. 92235.03c
    character(10) :: name_0K            ! name of 0K nuclide, e.g. 92235.00c
    real(8)       :: E_min = 0.01e-6    ! lower cutoff energy for res scattering
    real(8)       :: E_max = 1000.0e-6  ! upper cutoff energy for res scattering

  end type Nuclide0K

!===============================================================================
! DISTENERGYSAB contains the secondary energy/angle distributions for inelastic
! thermal scattering collisions which utilize a continuous secondary energy
! representation.
!===============================================================================

  type DistEnergySab
    integer              :: n_e_out
    real(8), allocatable :: e_out(:)
    real(8), allocatable :: e_out_pdf(:)
    real(8), allocatable :: e_out_cdf(:)
    real(8), allocatable :: mu(:,:)
  end type DistEnergySab

!===============================================================================
! SALPHABETA contains S(a,b) data for thermal neutron scattering, typically off
! of light isotopes such as water, graphite, Be, etc
!===============================================================================

  type SAlphaBeta
    character(10) :: name     ! name of table, e.g. lwtr.10t
    real(8)       :: awr      ! weight of nucleus in neutron masses
    real(8)       :: kT       ! temperature in MeV (k*T)
    integer       :: n_zaid   ! Number of valid zaids
    integer, allocatable :: zaid(:) ! List of valid Z and A identifiers, e.g. 6012

    ! threshold for S(a,b) treatment (usually ~4 eV)
    real(8) :: threshold_inelastic
    real(8) :: threshold_elastic = 0.0

    ! Inelastic scattering data
    integer :: n_inelastic_e_in  ! # of incoming E for inelastic
    integer :: n_inelastic_e_out ! # of outgoing E for inelastic
    integer :: n_inelastic_mu    ! # of outgoing angles for inelastic
    integer :: secondary_mode    ! secondary mode (equal/skewed/continuous)
    real(8), allocatable :: inelastic_e_in(:)
    real(8), allocatable :: inelastic_sigma(:)
    ! The following are used only if secondary_mode is 0 or 1
    real(8), allocatable :: inelastic_e_out(:,:)
    real(8), allocatable :: inelastic_mu(:,:,:)
    ! The following is used only if secondary_mode is 3
    ! The different implementation is necessary because the continuous
    ! representation has a variable number of outgoing energy points for each
    ! incoming energy
    type(DistEnergySab), allocatable :: inelastic_data(:) ! One for each Ein

    ! Elastic scattering data
    integer :: elastic_mode   ! elastic mode (discrete/exact)
    integer :: n_elastic_e_in ! # of incoming E for elastic
    integer :: n_elastic_mu   ! # of outgoing angles for elastic
    real(8), allocatable :: elastic_e_in(:)
    real(8), allocatable :: elastic_P(:)
    real(8), allocatable :: elastic_mu(:,:)
  end type SAlphaBeta

!===============================================================================
! XSLISTING contains data read from a cross_sections.xml file
!===============================================================================

  type XsListing
    character(12) :: name       ! table name, e.g. 92235.70c
    character(12) :: alias      ! table alias, e.g. U-235.70c
    integer       :: type       ! type of table (cont-E neutron, S(A,b), etc)
    integer       :: zaid       ! ZAID identifier = 1000*Z + A
    integer       :: filetype   ! ASCII or BINARY
    integer       :: location   ! location of table within library
    integer       :: recl       ! record length for library
    integer       :: entries    ! number of entries per record
    real(8)       :: awr        ! atomic weight ratio (# of neutron masses)
    real(8)       :: kT         ! Boltzmann constant * temperature (MeV)
    logical       :: metastable ! is this nuclide metastable?
    character(MAX_FILE_LEN) :: path ! path to library containing table
  end type XsListing

!===============================================================================
! NUCLIDEMICROXS contains cached microscopic cross sections for a
! particular nuclide at the current energy
!===============================================================================

  type NuclideMicroXS
    integer :: index_grid      ! index on nuclide energy grid
    integer :: index_temp      ! temperature index for nuclide
    real(8) :: last_E = 0.0    ! last evaluated energy
    real(8) :: interp_factor   ! interpolation factor on nuc. energy grid
    real(8) :: total           ! microscropic total xs
    real(8) :: elastic         ! microscopic elastic scattering xs
    real(8) :: absorption      ! microscopic absorption xs
    real(8) :: fission         ! microscopic fission xs
    real(8) :: nu_fission      ! microscopic production xs
    real(8) :: kappa_fission   ! microscopic energy-released from fission

    ! Information for S(a,b) use
    integer :: index_sab          ! index in sab_tables (zero means no table)
    integer :: last_index_sab = 0 ! index in sab_tables last used by this nuclide
    real(8) :: elastic_sab        ! microscopic elastic scattering on S(a,b) table

    ! Information for URR probability table use
    logical :: use_ptable  ! in URR range with probability tables?
    real(8) :: last_prn
  end type NuclideMicroXS

!===============================================================================
! MATERIALMACROXS contains cached macroscopic cross sections for the material a
! particle is traveling through
!===============================================================================

  type MaterialMacroXS
    real(8) :: total         ! macroscopic total xs
    real(8) :: elastic       ! macroscopic elastic scattering xs
    real(8) :: absorption    ! macroscopic absorption xs
    real(8) :: fission       ! macroscopic fission xs
    real(8) :: nu_fission    ! macroscopic production xs
    real(8) :: kappa_fission ! macroscopic energy-released from fission
  end type MaterialMacroXS

contains

!===============================================================================
! DISTANGLE_CLEAR resets and deallocates data in Reaction.
!===============================================================================

  subroutine distangle_clear(this)

    class(DistAngle), intent(inout) :: this ! The DistAngle object to clear

    if (allocated(this % energy)) &
      deallocate(this % energy, this % type, this % location, this % data)

  end subroutine distangle_clear

!===============================================================================
! DISTENERGY_CLEAR resets and deallocates data in DistEnergy.
!===============================================================================

  recursive subroutine distenergy_clear(this)

    class(DistEnergy), intent(inout) :: this ! The DistEnergy object to clear

    ! Clear p_valid
    call this % p_valid % clear()

    if (allocated(this % data)) &
      deallocate(this % data)

    if (associated(this % next)) then
      ! recursively clear this item
      call this % next % clear()
      deallocate(this % next)
    end if

  end subroutine distenergy_clear

!===============================================================================
! REACTION_CLEAR resets and deallocates data in Reaction.
!===============================================================================

  subroutine reaction_clear(this)

    class(Reaction), intent(inout) :: this ! The Reaction object to clear

    if (allocated(this % sigma)) &
      deallocate(this % sigma)

    if (associated(this % edist)) then
      call this % edist % clear()
      deallocate(this % edist)
    end if

    call this % adist % clear()

  end subroutine reaction_clear

!===============================================================================
! URRDATA_CLEAR resets and deallocates data in Reaction.
!===============================================================================

  subroutine urrdata_clear(this)

    class(UrrData), intent(inout) :: this ! The UrrData object to clear

    if (allocated(this % energy)) &
      deallocate(this % energy, this % prob)

  end subroutine urrdata_clear

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! ENERGY_RANGE_ALLOC allocated variables for the resonance energy ranges
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine energy_range_alloc(this)

    class(Nuclide), intent(inout) :: this ! nuclide object

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

  end subroutine energy_range_alloc

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! !TODO: ENERGY_RANGE_DEALLOC deallocates variables for the resonance energy ranges
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine energy_range_dealloc(this)

    class(Nuclide), intent(inout) :: this ! nuclide object

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

  end subroutine energy_range_dealloc

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! AVG_URR_ALLOC allocates average (infinite-dilute) cross sections
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine avg_urr_alloc(this)

    class(Nuclide), intent(inout) :: this ! nuclide object

    allocate(this % avg_urr_n(this % NE))
    allocate(this % avg_urr_f(this % NE))
    allocate(this % avg_urr_g(this % NE))
    allocate(this % avg_urr_x(this % NE))

  end subroutine avg_urr_alloc

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! !TODO: AVG_URR_DEALLOC deallocates average (infinite-dilute) cross sections
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine avg_urr_dealloc(this)

    class(Nuclide), intent(inout) :: this ! nuclide object

    deallocate(this % avg_urr_n)
    deallocate(this % avg_urr_f)
    deallocate(this % avg_urr_g)
    deallocate(this % avg_urr_x)

  end subroutine avg_urr_dealloc

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! ENSEMBLE_ALLOC allocates a URR resonance ensemble realization
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine ensemble_alloc(this)

    class(Nuclide), intent(inout) :: this ! nuclide object
    integer :: i_l   ! orbital angular momentum quantum number index
    integer :: i_lam ! resonance index

    ! allocate energies and orbital quantum numbers for resonances
    allocate(this % urr_resonances(this % n_urr_resonances, &
      & this % NLS(this % i_urr)))

    ! loop over orbital quantum numbers
    do i_l = 1, this % NLS(this % i_urr)

      ! loop over resonances
      do i_lam = 1, this % n_urr_resonances

        ! allocate resonance parameters
        call this % urr_resonances(i_lam, i_l) % alloc_resonances(this % NJS(i_l))

      end do
    end do

  end subroutine ensemble_alloc

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! TODO: ! ENSEMBLE_DEALLOC deallocates a URR resonance ensemble realization
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine ensemble_dealloc(this)

    class(Nuclide), intent(inout) :: this ! nuclide object

    

  end subroutine ensemble_dealloc

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! RESONANCES_ALLOC allocates a vector of URR resonances for a given J, for a
! given (i_lam, i_l)
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine resonances_alloc(this, N_J)

    class(URRResonances), intent(inout) :: this ! resonance vector object
    integer :: N_J

    allocate(this % E_lam(N_J))
    allocate(this % GN(N_J))
    allocate(this % GG(N_J))
    allocate(this % GF(N_J))
    allocate(this % GX(N_J))
    allocate(this % GT(N_J))

  end subroutine resonances_alloc

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! TODO: ! RESONANCES_DEALLOC deallocates a vector of URR resonances for a given J, for a
! given (i_lam, i_l)
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine resonances_dealloc(this)

    class(URRResonances), intent(inout) :: this ! resonance vector object

    

  end subroutine resonances_dealloc

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! POINTWISE_TMP_ALLOC allocates the temporary pointwise URR energy-cross section
! grids
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine pointwise_tmp_alloc(this)

    class(Nuclide), intent(inout) :: this ! nuclide object

    allocate(this % urr_energy_tmp(this % n_urr_gridpoints))
    allocate(this % urr_elastic_tmp(this % n_urr_gridpoints))
    allocate(this % urr_capture_tmp(this % n_urr_gridpoints))
    allocate(this % urr_fission_tmp(this % n_urr_gridpoints))
    allocate(this % urr_inelastic_tmp(this % n_urr_gridpoints))
    allocate(this % urr_total_tmp(this % n_urr_gridpoints))

  end subroutine pointwise_tmp_alloc

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! POINTWISE_TMP_DEALLOC deallocates the temporary pointwise URR energy-cross
! section grids
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine pointwise_tmp_dealloc(this)

    class(Nuclide), intent(inout) :: this ! nuclide object

    deallocate(this % urr_energy_tmp)
    deallocate(this % urr_elastic_tmp)
    deallocate(this % urr_capture_tmp)
    deallocate(this % urr_fission_tmp)
    deallocate(this % urr_inelastic_tmp)
    deallocate(this % urr_total_tmp)

  end subroutine pointwise_tmp_dealloc

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! POINTWISE_ALLOC allocates the pointwise URR energy-cross section grids
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine pointwise_alloc(this, n_pts)

    class(Nuclide), intent(inout) :: this ! nuclide object
    integer :: n_pts ! number of points in grid

    allocate(this % urr_energy(n_pts))
    allocate(this % urr_elastic(n_pts))
    allocate(this % urr_capture(n_pts))
    allocate(this % urr_fission(n_pts))
    allocate(this % urr_inelastic(n_pts))
    allocate(this % urr_total(n_pts))

  end subroutine pointwise_alloc

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! TODO: use this : POINTWISE_DEALLOC deallocates the pointwise URR energy-cross section grids
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine pointwise_dealloc(this)

    class(Nuclide), intent(inout) :: this ! nuclide object

    deallocate(this % urr_energy)
    deallocate(this % urr_elastic)
    deallocate(this % urr_capture)
    deallocate(this % urr_fission)
    deallocate(this % urr_inelastic)
    deallocate(this % urr_total)

  end subroutine pointwise_dealloc

!===============================================================================
! !TODO: deallocate URR structures NUCLIDE_CLEAR resets and deallocates data in Nuclide.
!===============================================================================

    subroutine nuclide_clear(this)

      class(Nuclide), intent(inout) :: this ! The Nuclide object to clear

      integer :: i ! Loop counter

      if (allocated(this % grid_index)) &
           deallocate(this % grid_index)

      if (allocated(this % energy)) &
           deallocate(this % energy, this % total, this % elastic, &
           & this % fission, this % nu_fission, this % absorption)

      if (allocated(this % energy_0K)) &
           deallocate(this % energy_0K)

      if (allocated(this % elastic_0K)) &
           deallocate(this % elastic_0K)

      if (allocated(this % xs_cdf)) &
           deallocate(this % xs_cdf)

      if (allocated(this % heating)) &
           deallocate(this % heating)

      if (allocated(this % index_fission)) deallocate(this % index_fission)

      if (allocated(this % nu_t_data)) deallocate(this % nu_t_data)
      if (allocated(this % nu_p_data)) deallocate(this % nu_p_data)
      if (allocated(this % nu_d_data)) deallocate(this % nu_d_data)

      if (allocated(this % nu_d_precursor_data)) &
           deallocate(this % nu_d_precursor_data)

      if (associated(this % nu_d_edist)) then
        do i = 1, size(this % nu_d_edist)
          call this % nu_d_edist(i) % clear()
        end do
        deallocate(this % nu_d_edist)
      end if

      if (associated(this % urr_data)) then
        call this % urr_data % clear()
        deallocate(this % urr_data)
      end if

      if (associated(this % reactions)) then
        do i = 1, size(this % reactions)
          call this % reactions(i) % clear()
        end do
        deallocate(this % reactions)
      end if

      call this % nuc_list % clear()

    end subroutine nuclide_clear

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! PROCESS_PRE pre-processes any nuclear data that needs it
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine process_pre(this, i_ER)

    class(Nuclide), intent(inout) :: this ! nuclide object
    integer :: i_ER ! resonance energy range index

    ! set the channel radius
    call this % channel_radius(i_ER)

  end subroutine process_pre

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! RADIUS_CHANNEL computes or sets the channel radius depending on ENDF flags
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine radius_channel(this, i_ER)

    class(Nuclide), intent(inout) :: this ! nuclide object
    integer :: i_ER ! resonance energy range index

    select case(this % NRO(i_ER))

    ! scattering radius is independent of energy
    case(0)

      select case(this % NAPS(i_ER))

        ! use channel radius for penetrabilities and shift factors but
        ! scattering radius for phase shifts
        case(0)
          this % ac(i_ER) = 0.123_8 * this % awr**(ONE/THREE) + 0.08_8

        ! use scattering radius for penetrabilities, shift factors and phase
        ! shifts
        case(1)
          this % ac(i_ER) = this % AP(i_ER)

        ! invalid scattering radius treatment flag
        case default
          print*, 'ENDF NAPS flag must be 0 or 1 when NRO is 0'
          stop
      end select

    ! scattering radius is energy dependent
    case(1)

      select case(this % NAPS(i_ER))

        ! use channel radius for penetrabilities and shift factors but
        ! scattering radius for phase shifts
        case(0)
          this % ac(i_ER) = 0.123_8 * this % awr**(ONE/THREE) + 0.08_8

        ! use scattering radius for penetrabilities, shift factors and phase
        ! shifts
        case(1)
          this % ac(i_ER) = this % AP(i_ER)

! TODO: understand this and implement it correctly
        ! use energy dependent scattering radius in phase shifts but the energy
        ! independent AP scattering radius value for penetrabilities and shift
        ! factors
        case(2)
          this % ac(i_ER) = this % AP(i_ER)

        ! invalid scattering radius treatment flag
        case default
          print*, 'ENDF NAPS flag must be 0, 1, or 2 when NRO is 1'
          stop
      end select

    ! invalid energy dependence of scattering radius flag
    case default
      print*, 'ENDF NRO flag must be 0 or 1'
      stop
  
    end select

  end subroutine radius_channel

end module ace_header
