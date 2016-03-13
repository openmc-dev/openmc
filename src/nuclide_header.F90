module nuclide_header

  use, intrinsic :: ISO_FORTRAN_ENV

  use ace_header
  use constants
  use endf,        only: reaction_name
  use error,       only: fatal_error
  use list_header, only: ListInt
  use scattdata_header
  use string
  use xml_interface

  implicit none

!===============================================================================
! Nuclide contains the base nuclidic data for a nuclide described as needed
! for continuous-energy neutron transport.
!===============================================================================

  type :: Nuclide
    ! Nuclide meta-data
    character(12) :: name    ! name of nuclide, e.g. 92235.03c
    integer       :: zaid    ! Z and A identifier, e.g. 92235
    real(8)       :: awr     ! Atomic Weight Ratio
    integer       :: listing ! index in xs_listings
    real(8)       :: kT      ! temperature in MeV (k*T)

    ! Fission information
    logical :: fissionable   ! nuclide is fissionable?

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
    type(AngleEnergyContainer), allocatable :: nu_d_edist(:)

    ! Unresolved resonance data
    logical                :: urr_present
    integer                :: urr_inelastic
    type(UrrData), pointer :: urr_data => null()

    ! Reactions
    integer :: n_reaction ! # of reactions
    type(Reaction), allocatable :: reactions(:)
    type(DictIntInt) :: reaction_index ! map MT values to index in reactions
                                       ! array; used at tally-time

  contains
    procedure :: clear => nuclide_clear
    procedure :: print => nuclide_print
  end type Nuclide

!===============================================================================
! NUCLIDE0K temporarily contains all 0K cross section data and other parameters
! needed to treat resonance scattering before transferring them to Nuclide
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
! NUCLIDE_CLEAR resets and deallocates data in Nuclide, NuclideIso
! or NuclideAngle
!===============================================================================

    subroutine nuclide_clear(this)

      class(Nuclide), intent(inout) :: this ! The Nuclide object to clear

      integer :: i ! Loop counter

      if (associated(this % urr_data)) deallocate(this % urr_data)

      if (allocated(this % reactions)) then
        do i = 1, size(this % reactions)
          call this % reactions(i) % clear()
        end do
      end if

      call this % reaction_index % clear()

    end subroutine nuclide_clear

!===============================================================================
! NUCLIDE_PRINT displays information about a continuous-energy neutron
! cross_section table and its reactions and secondary angle/energy distributions
!===============================================================================

    subroutine nuclide_print(this, unit)
      class(Nuclide), intent(in) :: this
      integer, intent(in), optional :: unit

      integer :: i                 ! loop index over nuclides
      integer :: unit_             ! unit to write to
      integer :: size_xs           ! memory used for cross-sections (bytes)
      integer :: size_urr          ! memory used for probability tables (bytes)
      type(UrrData),  pointer :: urr

      ! set default unit for writing information
      if (present(unit)) then
        unit_ = unit
      else
        unit_ = OUTPUT_UNIT
      end if

      ! Initialize totals
      size_urr = 0
      size_xs = 0

      ! Basic nuclide information
      write(unit_,*) 'Nuclide ' // trim(this % name)
      write(unit_,*) '  zaid = ' // trim(to_str(this % zaid))
      write(unit_,*) '  awr = ' // trim(to_str(this % awr))
      write(unit_,*) '  kT = ' // trim(to_str(this % kT))
      write(unit_,*) '  # of grid points = ' // trim(to_str(this % n_grid))
      write(unit_,*) '  Fissionable = ', this % fissionable
      write(unit_,*) '  # of fission reactions = ' // trim(to_str(this % n_fission))
      write(unit_,*) '  # of reactions = ' // trim(to_str(this % n_reaction))

      ! Information on each reaction
      write(unit_,*) '  Reaction     Q-value  COM    IE'
      do i = 1, this % n_reaction
        associate (rxn => this % reactions(i))
          write(unit_,'(3X,A11,1X,F8.3,3X,L1,3X,I6)') &
               reaction_name(rxn % MT), rxn % Q_value, rxn % scatter_in_cm, &
               rxn % threshold

          ! Accumulate data size
          size_xs = size_xs + (this % n_grid - rxn%threshold + 1) * 8
        end associate
      end do

      ! Add memory required for summary reactions (total, absorption, fission,
      ! nu-fission)
      size_xs = 8 * this % n_grid * 4

      ! Write information about URR probability tables
      size_urr = 0
      if (this % urr_present) then
        urr => this % urr_data
        write(unit_,*) '  Unresolved resonance probability table:'
        write(unit_,*) '    # of energies = ' // trim(to_str(urr % n_energy))
        write(unit_,*) '    # of probabilities = ' // trim(to_str(urr % n_prob))
        write(unit_,*) '    Interpolation =  ' // trim(to_str(urr % interp))
        write(unit_,*) '    Inelastic flag = ' // trim(to_str(urr % inelastic_flag))
        write(unit_,*) '    Absorption flag = ' // trim(to_str(urr % absorption_flag))
        write(unit_,*) '    Multiply by smooth? ', urr % multiply_smooth
        write(unit_,*) '    Min energy = ', trim(to_str(urr % energy(1)))
        write(unit_,*) '    Max energy = ', trim(to_str(urr % energy(urr % n_energy)))

        ! Calculate memory used by probability tables and add to total
        size_urr = urr % n_energy * (urr % n_prob * 6 + 1) * 8
      end if

      ! Write memory used
      write(unit_,*) '  Memory Requirements'
      write(unit_,*) '    Cross sections = ' // trim(to_str(size_xs)) // ' bytes'
      write(unit_,*) '    Probability Tables = ' // &
           trim(to_str(size_urr)) // ' bytes'

      ! Blank line at end of nuclide
      write(unit_,*)
    end subroutine nuclide_print

end module nuclide_header
