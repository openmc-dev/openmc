module nuclide_header

  use ace_header
  use constants
  use list_header,     only: ListInt
  ! use math,            only: calc_pn, calc_rn!, expand_harmonic
  !use scattdata_header
  ! use xml_interface

  implicit none

!===============================================================================
! NUCLIDE_BASE contains the base nuclidic data for a nuclide, which does not depend
! upon how the nuclear data is represented (i.e., CE, or any variant of MG).
! The extended types, Nuclide_CE and Nuclide_MG deal with the rest
!===============================================================================

  type, abstract :: Nuclide_Base
    character(12) :: name    ! name of nuclide, e.g. 92235.03c
    integer       :: zaid    ! Z and A identifier, e.g. 92235
    real(8)       :: awr     ! Atomic Weight Ratio
    integer       :: listing ! index in xs_listings
    real(8)       :: kT      ! temperature in MeV (k*T)

    ! Linked list of indices in nuclides array of instances of this same nuclide
    type(ListInt) :: nuc_list

    ! Fission information
    logical :: fissionable         ! nuclide is fissionable?

    contains
      procedure(nuclide_base_clear_), deferred, pass :: clear ! Deallocates Nuclide
  end type Nuclide_Base

  type, extends(Nuclide_Base) :: Nuclide_CE
    ! Energy grid information
    integer :: n_grid                     ! # of nuclide grid points
    integer, allocatable :: grid_index(:) ! log grid mapping indices
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

    ! Unresolved resonance data
    logical                :: urr_present
    integer                :: urr_inelastic
    type(UrrData), pointer :: urr_data => null()

    ! Reactions
    integer :: n_reaction ! # of reactions
    type(Reaction), pointer :: reactions(:) => null()

    ! Type-Bound procedures
    contains
      procedure, pass :: clear => nuclide_ce_clear
  end type Nuclide_CE

!===============================================================================
! NUCLIDECONTAINER pointer array for storing Nuclides
!===============================================================================

  type NuclideContainer
    class(Nuclide_Base), pointer :: obj
  end type NuclideContainer

!===============================================================================
! NUCLIDE0K temporarily contains all 0K cross section data and other parameters
! needed to treat resonance scattering before transferring them to NUCLIDE_CE
!===============================================================================

  type Nuclide0K

    character(10) :: nuclide             ! name of nuclide, e.g. U-238
    character(16) :: scheme = 'ares'     ! target velocity sampling scheme
    character(10) :: name                ! name of nuclide, e.g. 92235.03c
    character(10) :: name_0K             ! name of 0K nuclide, e.g. 92235.00c
    real(8)       :: E_min = 0.01e-6_8   ! lower cutoff energy for res scattering
    real(8)       :: E_max = 1000.0e-6_8 ! upper cutoff energy for res scattering

  end type Nuclide0K

!===============================================================================
! NUCLIDEMICROXS contains cached microscopic cross sections for a
! particular nuclide at the current energy
!===============================================================================

  type NuclideMicroXS
    integer :: index_grid      ! index on nuclide energy grid
    integer :: index_temp      ! temperature index for nuclide
    real(8) :: last_E = ZERO   ! last evaluated energy
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

!===============================================================================
! XSLISTING contains data read from a CE or MG cross_sections.xml file
! (or equivalent)
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

  contains


!===============================================================================
! NUCLIDE_*_CLEAR resets and deallocates data in Nuclide_Base, Nuclide_Iso
! or Nuclide_Angle
!===============================================================================

    subroutine nuclide_base_clear_(this)

      class(Nuclide_Base), intent(inout) :: this

      call this % nuc_list % clear()
    end subroutine nuclide_base_clear_

    subroutine nuclide_ce_clear(this)

      class(Nuclide_CE), intent(inout) :: this ! The Nuclide object to clear

      integer :: i ! Loop counter

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

      call nuclide_base_clear_(this)

    end subroutine nuclide_ce_clear

  end module nuclide_header