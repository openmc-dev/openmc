module ace_header

  use constants,     only: MAX_FILE_LEN, ONE, THREE
  use endf_header,   only: Tab1
  use vector_header, only: URRVector

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

    ! Unresolved resonance information
    logical                :: otf_urr = .false.
    integer                :: n_resonances = 10
    logical                :: urr_present
    integer                :: urr_inelastic
    type(UrrData), pointer :: urr_data => null()

    ! ENDF URR data
    real(8) :: EL           ! lower energy bound of URR
    real(8) :: EH           ! upper energy bound of URR
    integer :: NRO          ! AP energy-dependence flag
    integer :: NAPS         ! channel radius handling flag
    real(8) :: SPI          ! total spin
    real(8) :: AP           ! scattering radius
    real(8) :: ac           ! channel radius
    integer :: LSSF         ! self-shielding factor flag
    integer :: NLS          ! # of orbital quantum #'s
    integer :: L            ! current orbital quantum #
    integer :: LRX          ! competitive inelastic width flag
    integer :: LRP          ! resonance parameter flag
    integer :: LRU
    integer :: LRF
    integer :: INT

! TODO: only accept LRP of 1
! TODO: use LRX if ZERO's aren't given for inelastic width in ENDF when LRX = 0

    ! # of total angular momenta for each orbital quantum #
    integer,         allocatable :: NJS(:)
    type(URRVector), allocatable :: J_grid(:) ! values at each orbital quantum #
    real(8)                      :: J         ! current total angular momentum

    ! degrees of freedom at each orbital quantum # (one value for each J value)
    type(URRVector), allocatable :: AMUX_grid(:) ! competitive reaction values
    integer                      :: AMUX         ! current value
    type(URRVector), allocatable :: AMUN_grid(:) ! elastic neutron scatter values
    integer                      :: AMUN         ! current value
    type(URRVector), allocatable :: AMUG_grid(:) ! capture values
    integer                      :: AMUG         ! current value
    type(URRVector), allocatable :: AMUF_grid(:) ! fission values
    integer                      :: AMUF         ! current value

    ! energy grid of unresolved resonance parameters
    integer                      :: NE
    real(8), allocatable         :: ES(:)
    real(8)                      :: E          ! neutron energy

    ! mean level spacing for a value of J for a given (E,l)
    type(URRVector), allocatable :: D_means(:,:)
    real(8)                      :: D          ! current value

    ! mean neutron width for a value of J for a given (E,l)
    type(URRVector), allocatable :: Gam_n_means(:,:)
    real(8)                      :: GN0        ! current value

    ! mean radiative width for a value of J for a given (E,l)
    type(URRVector), allocatable :: Gam_gam_means(:,:)
    real(8)                      :: GG         ! current mean value

    ! mean fission width for a value of J for a given (E,l)
    type(URRVector), allocatable :: Gam_f_means(:,:)
    real(8)                      :: GF         ! current value

    ! mean competitive width for a value of J for a given (E,l)
    type(URRVector), allocatable :: Gam_x_means(:,:)
    real(8)                      :: GX         ! current value

    ! Reactions
    integer :: n_reaction ! # of reactions
    type(Reaction), pointer :: reactions(:) => null()

    ! Type-Bound procedures
    contains

      ! allocate a spin sequence for the nuclide
      procedure :: alloc_lJ => alloc_spin_seq

      ! Deallocates Nuclide
      procedure :: clear => nuclide_clear

      ! pre-process data that needs it
      procedure :: pprocess => pre_process

      ! set channel radius
      procedure :: rad_channel => channel_radius

  end type Nuclide

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
! ALLOC_SPIN_SEQ allocates memory for a single spin sequence
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine alloc_spin_seq(this)
    
    integer :: i_L
    integer :: i_E

    class(Nuclide), intent(inout) :: this ! nuclide object

    ! allocate energy grid
    allocate(this % ES(this % NE))

    ! allocate total angular momentum quantum #'s
    allocate(this % NJS(this % NLS))
    allocate(this % J_grid(this % NLS))

    ! allocate degress of freedom for partial widths
    allocate(this % AMUX_grid(this % NLS))
    allocate(this % AMUN_grid(this % NLS))
    allocate(this % AMUG_grid(this % NLS))
    allocate(this % AMUF_grid(this % NLS))

    ! allocate mean widths and spacings
    allocate(this % D_means(this % NE, this % NLS))
    allocate(this % Gam_n_means(this % NE, this % NLS))
    allocate(this % Gam_gam_means(this % NE, this % NLS))
    allocate(this % Gam_f_means(this % NE, this % NLS))
    allocate(this % Gam_x_means(this % NE, this % NLS))

    ! allocate space for the different spin sequences (i.e. (l,J) pairs)
    do i_L = 1, this % NLS

      allocate(this % J_grid(i_L)    % vals(this % NJS(i_L)))
      allocate(this % AMUX_grid(i_L) % vals(this % NJS(i_L)))
      allocate(this % AMUN_grid(i_L) % vals(this % NJS(i_L)))
      allocate(this % AMUG_grid(i_L) % vals(this % NJS(i_L)))
      allocate(this % AMUF_grid(i_L) % vals(this % NJS(i_L)))

      do i_E = 1, this % NE
        allocate(this % D_means(i_E, i_L)       % vals(this % NJS(i_L)))
        allocate(this % Gam_n_means(i_E, i_L)   % vals(this % NJS(i_L)))
        allocate(this % Gam_gam_means(i_E, i_L) % vals(this % NJS(i_L)))
        allocate(this % Gam_f_means(i_E, i_L)   % vals(this % NJS(i_L)))
        allocate(this % Gam_x_means(i_E, i_L)   % vals(this % NJS(i_L)))
      end do
    end do

  end subroutine alloc_spin_seq

!===============================================================================
! NUCLIDE_CLEAR resets and deallocates data in Nuclide.
!===============================================================================

    subroutine nuclide_clear(this)

      class(Nuclide), intent(inout) :: this ! The Nuclide object to clear

      integer :: i ! Loop counter
      integer :: i_L ! orbital quantum # index
      integer :: i_E ! energy grid index

      if (allocated(this % grid_index)) deallocate(this % grid_index)

      if (allocated(this % energy)) &
           deallocate(this % energy, this % total, this % elastic, &
           this % fission, this % nu_fission, this % absorption)
      if (allocated(this % heating)) deallocate(this % heating)

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

      ! deallocate energy grid
      if (allocated(this % ES)) deallocate(this % ES)
      
      ! deallocate space for the different spin sequences (i.e. (l,J) pairs)
      do i_L = 1, this % NLS
        
        if (allocated(this % J_grid(i_L) % vals)) &
          & deallocate(this % J_grid(i_L) % vals)
        
        if (allocated(this % AMUX_grid(i_L) % vals)) &
          & deallocate(this % AMUX_grid(i_L) % vals)
        
        if (allocated(this % AMUN_grid(i_L) % vals)) &
          & deallocate(this % AMUN_grid(i_L) % vals)
        
        if (allocated(this % AMUG_grid(i_L) % vals)) &
          & deallocate(this % AMUG_grid(i_L) % vals)
        
        if (allocated(this % AMUF_grid(i_L) % vals)) &
          & deallocate(this % AMUF_grid(i_L) % vals)
        
        do i_E = 1, this % NE
          if (allocated(this % D_means(i_E, i_L) % vals)) &
            & deallocate(this % D_means(i_E, i_L) % vals)
          if (allocated(this % Gam_n_means(i_E, i_L) % vals)) &
            & deallocate(this % Gam_n_means(i_E, i_L) % vals)
          if (allocated(this % Gam_gam_means(i_E, i_L) % vals)) &
            & deallocate(this % Gam_gam_means(i_E, i_L) % vals)
          if (allocated(this % Gam_f_means(i_E, i_L) % vals)) &
            & deallocate(this % Gam_f_means(i_E, i_L) % vals)
          if (allocated(this % Gam_x_means(i_E, i_L) % vals)) &
            & deallocate(this % Gam_x_means(i_E, i_L) % vals)
        end do
      end do
      
      ! deallocate total angular momentum quantum #'s
      if (allocated(this % NJS))     deallocate(this % NJS)
      if (allocated(this % J_grid))  deallocate(this % J_grid)
      
      ! deallocate mean unresolved resonance parameters
      if (allocated(this % AMUX_grid))     deallocate(this % AMUX_grid)
      if (allocated(this % AMUN_grid))     deallocate(this % AMUN_grid)
      if (allocated(this % AMUG_grid))     deallocate(this % AMUG_grid)
      if (allocated(this % AMUF_grid))     deallocate(this % AMUF_grid)
      if (allocated(this % D_means))       deallocate(this % D_means)
      if (allocated(this % Gam_n_means))   deallocate(this % Gam_n_means)
      if (allocated(this % Gam_gam_means)) deallocate(this % Gam_gam_means)
      if (allocated(this % Gam_f_means))   deallocate(this % Gam_f_means)
      if (allocated(this % Gam_x_means))   deallocate(this % Gam_x_means)
      
    end subroutine nuclide_clear

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! PRE_PROCESS pre-processes any nuclear data that needs it
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine pre_process(this)

    class(Nuclide), intent(inout) :: this ! nuclide object

    ! set the channel radius
    call this % rad_channel

  end subroutine pre_process

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!
! CHANNEL_RADIUS computes or sets the channel radius depending on ENDF flags
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  subroutine channel_radius(this)

    class(Nuclide), intent(inout) :: this ! nuclide object

    select case(this % NRO)

    ! scattering radius is independent of energy
    case(0)

      select case(this % NAPS)

        ! use channel radius for penetrabilities and shift factors but
        ! scattering radius for phase shifts
        case(0)
          this % ac = 0.123_8 * this % awr**(ONE/THREE) + 0.08_8

        ! use scattering radius for penetrabilities, shift factors and phase
        ! shifts
        case(1)
          this % ac = this % AP

        ! invalid scattering radius treatment flag
        case default
          print*, 'ENDF NAPS flag must be 0 or 1 when NRO is 0'
          stop
      end select

    ! scattering radius is energy dependent
    case(1)

      select case(this % NAPS)

        ! use channel radius for penetrabilities and shift factors but
        ! scattering radius for phase shifts
        case(0)
          this % ac = 0.123_8 * this % awr**(ONE/THREE) + 0.08_8

        ! use scattering radius for penetrabilities, shift factors and phase
        ! shifts
        case(1)
          this % ac = this % AP

! TODO: understand this and implement it correctly
        ! use energy dependent scattering radius in phase shifts but the energy
        ! independent AP scattering radius value for penetrabilities and shift
        ! factors
        case(2)
          this % ac = this % AP

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

  end subroutine channel_radius

end module ace_header
