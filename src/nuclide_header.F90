module nuclide_header

  use, intrinsic :: ISO_FORTRAN_ENV
  use, intrinsic :: ISO_C_BINDING

  use hdf5,                   only: HID_T, HSIZE_T, SIZE_T

  use algorithm,              only: sort, find, binary_search
  use constants
  use dict_header,            only: DictIntInt, DictCharInt
  use endf,                   only: reaction_name, is_fission, is_disappearance
  use endf_header,            only: Function1D, Polynomial, Tabulated1D
  use error
  use hdf5_interface
  use list_header,            only: ListInt
  use math,                   only: faddeeva, w_derivative, &
                                    broaden_wmp_polynomials
  use multipole_header,       only: FORM_RM, FORM_MLBW, MP_EA, RM_RT, RM_RA, &
                                    RM_RF, MLBW_RT, MLBW_RX, MLBW_RA, MLBW_RF, &
                                    FIT_T, FIT_A, FIT_F, MultipoleArray
  use message_passing
  use multipole_header,       only: MultipoleArray
  use product_header,         only: AngleEnergyContainer
  use random_lcg,             only: prn, future_prn, prn_set_stream
  use reaction_header,        only: Reaction
  use sab_header,             only: SAlphaBeta, sab_tables
  use secondary_uncorrelated, only: UncorrelatedAngleEnergy
  use settings
  use stl_vector,             only: VectorInt, VectorReal
  use string
  use simulation_header,      only: need_depletion_rx
  use urr_header,             only: UrrData

  implicit none

!===============================================================================
! Nuclide contains the base nuclidic data for a nuclide described as needed
! for continuous-energy neutron transport.
!===============================================================================

  type EnergyGrid
    integer, allocatable :: grid_index(:) ! log grid mapping indices
    real(8), allocatable :: energy(:)     ! energy values corresponding to xs
  end type EnergyGrid

  ! Positions for first dimension of Nuclide % xs
  integer, parameter :: &
       XS_TOTAL      = 1, &
       XS_ABSORPTION = 2, &
       XS_FISSION    = 3, &
       XS_NU_FISSION = 4

  ! The array within SumXS is of shape (4, n_energy) where the first dimension
  ! corresponds to the following values: 1) total, 2) absorption (MT > 100), 3)
  ! fission, 4) neutron production
  type SumXS
    real(8), allocatable :: value(:,:)
  end type SumXS

  type :: Nuclide
    ! Nuclide meta-data
    character(20) :: name    ! name of nuclide, e.g. U235.71c
    integer       :: Z       ! atomic number
    integer       :: A       ! mass number
    integer       :: metastable ! metastable state
    real(8)       :: awr     ! Atomic Weight Ratio
    integer       :: i_nuclide ! The nuclides index in the nuclides array
    real(8), allocatable :: kTs(:) ! temperature in eV (k*T)

    ! Fission information
    logical :: fissionable = .false.  ! nuclide is fissionable?

    ! Energy grid for each temperature
    type(EnergyGrid), allocatable :: grid(:)

    ! Microscopic cross sections
    type(SumXS), allocatable :: xs(:)

    ! Resonance scattering info
    logical              :: resonant = .false. ! resonant scatterer?
    real(8), allocatable :: energy_0K(:)  ! energy grid for 0K xs
    real(8), allocatable :: elastic_0K(:) ! Microscopic elastic cross section
    real(8), allocatable :: xs_cdf(:) ! CDF of v_rel times cross section

    ! Fission information
    logical :: has_partial_fission = .false. ! nuclide has partial fission reactions?
    integer :: n_fission = 0                 ! # of fission reactions
    integer :: n_precursor = 0               ! # of delayed neutron precursors
    integer, allocatable :: index_fission(:) ! indices in reactions
    class(Function1D), allocatable :: total_nu

    ! Unresolved resonance data
    logical                :: urr_present = .false.
    integer                :: urr_inelastic
    type(UrrData), allocatable :: urr_data(:)

    ! Multipole data
    logical                       :: mp_present = .false.
    type(MultipoleArray), pointer :: multipole => null()

    ! Reactions
    type(Reaction), allocatable :: reactions(:)
    integer, allocatable :: index_inelastic_scatter(:)

    ! Array that maps MT values to index in reactions; used at tally-time. Note
    ! that ENDF-102 does not assign any MT values above 891.
    integer :: reaction_index(891)

    ! Fission energy release
    class(Function1D), allocatable :: fission_q_prompt ! prompt neutrons, gammas
    class(Function1D), allocatable :: fission_q_recov  ! neutrons, gammas, betas

  contains
    procedure :: assign_0K_elastic_scattering
    procedure :: clear => nuclide_clear
    procedure :: from_hdf5 => nuclide_from_hdf5
    procedure :: init_grid => nuclide_init_grid
    procedure :: nu    => nuclide_nu
    procedure, private :: create_derived => nuclide_create_derived
    procedure :: calculate_xs => nuclide_calculate_xs
    procedure :: calculate_elastic_xs => nuclide_calculate_elastic_xs
  end type Nuclide

!===============================================================================
! NUCLIDEMICROXS contains cached microscopic cross sections for a particular
! nuclide at the current energy
!===============================================================================

  ! Arbitrary value to indicate invalid cache state for elastic scattering
  ! (NuclideMicroXS % elastic)
  real(8), parameter :: CACHE_INVALID = dble(Z"FFE0000000000000")

  type NuclideMicroXS
    ! Microscopic cross sections in barns
    real(8) :: total
    real(8) :: absorption       ! absorption (disappearance)
    real(8) :: fission          ! fission
    real(8) :: nu_fission       ! neutron production from fission

    real(8) :: elastic          ! If sab_frac is not 1 or 0, then this value is
                                !   averaged over bound and non-bound nuclei
    real(8) :: thermal          ! Bound thermal elastic & inelastic scattering
    real(8) :: thermal_elastic  ! Bound thermal elastic scattering

    ! Cross sections for depletion reactions (note that these are not stored in
    ! macroscopic cache)
    real(8) :: reaction(size(DEPLETION_RX))

    ! Indicies and factors needed to compute cross sections from the data tables
    integer :: index_grid        ! Index on nuclide energy grid
    integer :: index_temp        ! Temperature index for nuclide
    real(8) :: interp_factor     ! Interpolation factor on nuc. energy grid
    integer :: index_sab = NONE  ! Index in sab_tables
    integer :: index_temp_sab    ! Temperature index for sab_tables
    real(8) :: sab_frac          ! Fraction of atoms affected by S(a,b)
    logical :: use_ptable        ! In URR range with probability tables?

    ! Energy and temperature last used to evaluate these cross sections.  If
    ! these values have changed, then the cross sections must be re-evaluated.
    real(8) :: last_E = ZERO       ! Last evaluated energy
    real(8) :: last_sqrtkT = ZERO  ! Last temperature in sqrt(Boltzmann
                                   !   constant * temperature (eV))
  end type NuclideMicroXS

!===============================================================================
! MATERIALMACROXS contains cached macroscopic cross sections for the material a
! particle is traveling through
!===============================================================================

  type MaterialMacroXS
    real(8) :: total         ! macroscopic total xs
    real(8) :: absorption    ! macroscopic absorption xs
    real(8) :: fission       ! macroscopic fission xs
    real(8) :: nu_fission    ! macroscopic production xs
  end type MaterialMacroXS

!===============================================================================
! LIBRARY contains data read from a cross_sections.xml file
!===============================================================================

  type Library
    integer :: type
    character(MAX_WORD_LEN), allocatable :: materials(:)
    character(MAX_FILE_LEN) :: path
  end type Library

  ! Cross section libraries
  type(Library), allocatable :: libraries(:)
  type(DictCharInt) :: library_dict

  ! Nuclear data for each nuclide
  type(Nuclide), allocatable, target :: nuclides(:)
  integer(C_INT), bind(C) :: n_nuclides
  type(DictCharInt) :: nuclide_dict

  ! Cross section caches
  type(NuclideMicroXS), allocatable :: micro_xs(:)  ! Cache for each nuclide
  type(MaterialMacroXS)             :: material_xs  ! Cache for current material
!$omp threadprivate(micro_xs, material_xs)

  ! Minimum/maximum energies
  real(8) :: energy_min_neutron = ZERO
  real(8) :: energy_max_neutron = INFINITY

contains

!===============================================================================
! ASSIGN_0K_ELASTIC_SCATTERING
!===============================================================================

  subroutine assign_0K_elastic_scattering(this)
    class(Nuclide), intent(inout) :: this

    integer :: i
    real(8) :: xs_cdf_sum

    this % resonant = .false.
    if (allocated(res_scat_nuclides)) then
      ! If resonant nuclides were specified, check the list explicitly
      do i = 1, size(res_scat_nuclides)
        if (this % name == res_scat_nuclides(i)) then
          this % resonant = .true.

          ! Make sure nuclide has 0K data
          if (.not. allocated(this % energy_0K)) then
            call fatal_error("Cannot treat " // trim(this % name) // " as a &
                 &resonant scatterer because 0 K elastic scattering data is &
                 &not present.")
          end if

          exit
        end if
      end do
    else
      ! Otherwise, assume that any that have 0 K elastic scattering data are
      ! resonant
      this % resonant = allocated(this % energy_0K)
    end if

    if (this % resonant) then
      ! Build CDF for 0K elastic scattering
      xs_cdf_sum = ZERO
      allocate(this % xs_cdf(0:size(this % energy_0K)))
      this % xs_cdf(0) = ZERO

      associate (E => this % energy_0K, xs => this % elastic_0K)
        do i = 1, size(E) - 1
          ! Negative cross sections result in a CDF that is not monotonically
          ! increasing. Set all negative xs values to zero.
          if (xs(i) < ZERO) xs(i) = ZERO

          ! build xs cdf
          xs_cdf_sum = xs_cdf_sum + (sqrt(E(i))*xs(i) + sqrt(E(i+1))*xs(i+1))&
               / TWO * (E(i+1) - E(i))
          this % xs_cdf(i) = xs_cdf_sum
        end do
      end associate
    end if
  end subroutine assign_0K_elastic_scattering

!===============================================================================
! NUCLIDE_CLEAR resets and deallocates data in Nuclide
!===============================================================================

  subroutine nuclide_clear(this)
    class(Nuclide), intent(inout) :: this ! The Nuclide object to clear

    if (associated(this % multipole)) deallocate(this % multipole)

  end subroutine nuclide_clear

  subroutine nuclide_from_hdf5(this, group_id, temperature, method, tolerance, &
                               minmax, master, i_nuclide)
    class(Nuclide),   intent(inout) :: this
    integer(HID_T),   intent(in)    :: group_id
    type(VectorReal), intent(in)    :: temperature ! list of desired temperatures
    integer,          intent(inout) :: method
    real(8),          intent(in)    :: tolerance
    real(8),          intent(in)    :: minmax(2)  ! range of temperatures
    logical,          intent(in)    :: master     ! if this is the master proc
    integer,          intent(in)    :: i_nuclide  ! Nuclide index in nuclides

    integer :: i
    integer :: i_closest
    integer :: n_temperature
    integer(HID_T) :: urr_group, nu_group
    integer(HID_T) :: energy_group, energy_dset
    integer(HID_T) :: kT_group
    integer(HID_T) :: rxs_group
    integer(HID_T) :: rx_group
    integer(HID_T) :: xs, temp_group
    integer(HID_T) :: total_nu
    integer(HID_T) :: fer_group                 ! fission_energy_release group
    integer(HID_T) :: fer_dset
    integer(SIZE_T) :: name_len
    integer(HSIZE_T) :: j
    integer(HSIZE_T) :: dims(1)
    character(MAX_WORD_LEN) :: temp_str
    character(MAX_WORD_LEN), allocatable :: dset_names(:)
    character(MAX_WORD_LEN), allocatable :: grp_names(:)
    real(8), allocatable :: temps_available(:) ! temperatures available
    real(8) :: temp_desired
    real(8) :: temp_actual
    type(VectorInt) :: MTs
    type(VectorInt) :: temps_to_read
    type(VectorInt) :: index_inelastic_scatter_vector

    ! Get name of nuclide from group
    name_len = len(this % name)
    this % name = get_name(group_id, name_len)

    ! Get rid of leading '/'
    this % name = trim(this % name(2:))

    call read_attribute(this % Z, group_id, 'Z')
    call read_attribute(this % A, group_id, 'A')
    call read_attribute(this % metastable, group_id, 'metastable')
    call read_attribute(this % awr, group_id, 'atomic_weight_ratio')
    kT_group = open_group(group_id, 'kTs')

    ! Determine temperatures available
    call get_datasets(kT_group, dset_names)
    allocate(temps_available(size(dset_names)))
    do i = 1, size(dset_names)
      ! Read temperature value
      call read_dataset(temps_available(i), kT_group, trim(dset_names(i)))
      temps_available(i) = temps_available(i) / K_BOLTZMANN
    end do
    call sort(temps_available)

    ! If only one temperature is available, revert to nearest temperature
    if (size(temps_available) == 1 .and. method == TEMPERATURE_INTERPOLATION) then
      if (master) then
        call warning("Cross sections for " // trim(this % name) // " are only &
             &available at one temperature. Reverting to nearest temperature &
             &method.")
      end if
      method = TEMPERATURE_NEAREST
    end if

    ! Determine actual temperatures to read -- start by checking whether a
    ! temperature range was given, in which case all temperatures in the range
    ! are loaded irrespective of what temperatures actually appear in the model
    if (minmax(2) > ZERO) then
      do i = 1, size(temps_available)
        temp_actual = temps_available(i)
        if (minmax(1) <= temp_actual .and. temp_actual <= minmax(2)) then
          call temps_to_read % push_back(nint(temp_actual))
        end if
      end do
    end if

    select case (method)
    case (TEMPERATURE_NEAREST)
      ! Find nearest temperatures
      do i = 1, temperature % size()
        temp_desired = temperature % data(i)
        i_closest = minloc(abs(temps_available - temp_desired), dim=1)
        temp_actual = temps_available(i_closest)
        if (abs(temp_actual - temp_desired) < tolerance) then
          if (find(temps_to_read, nint(temp_actual)) == -1) then
            call temps_to_read % push_back(nint(temp_actual))

            ! Write warning for resonance scattering data if 0K is not available
            if (abs(temp_actual - temp_desired) > 0 .and. temp_desired == 0 &
                 .and. master) then
              call warning(trim(this % name) // " does not contain 0K data &
                   &needed for resonance scattering options selected. Using &
                   &data at " // trim(to_str(temp_actual)) &
                   // " K instead.")
            end if
          end if
        else
          call fatal_error("Nuclear data library does not contain cross &
               &sections for " // trim(this % name) // " at or near " // &
               trim(to_str(nint(temp_desired))) // " K.")
        end if
      end do

    case (TEMPERATURE_INTERPOLATION)
      ! If temperature interpolation or multipole is selected, get a list of
      ! bounding temperatures for each actual temperature present in the model
      TEMP_LOOP: do i = 1, temperature % size()
        temp_desired = temperature % data(i)

        do j = 1, size(temps_available) - 1
          if (temps_available(j) <= temp_desired .and. &
               temp_desired < temps_available(j + 1)) then
            if (find(temps_to_read, nint(temps_available(j))) == -1) then
              call temps_to_read % push_back(nint(temps_available(j)))
            end if
            if (find(temps_to_read, nint(temps_available(j + 1))) == -1) then
              call temps_to_read % push_back(nint(temps_available(j + 1)))
            end if
            cycle TEMP_LOOP
          end if
        end do

        call fatal_error("Nuclear data library does not contain cross sections &
             &for " // trim(this % name) // " at temperatures that bound " // &
             trim(to_str(nint(temp_desired))) // " K.")
      end do TEMP_LOOP

    end select

    ! Sort temperatures to read
    call sort(temps_to_read)

    n_temperature = temps_to_read % size()
    allocate(this % kTs(n_temperature))
    allocate(this % grid(n_temperature))

    ! Get kT values
    do i = 1, n_temperature
      ! Get temperature as a string
      temp_str = trim(to_str(temps_to_read % data(i))) // "K"

      ! Read exact temperature value
      call read_dataset(this % kTs(i), kT_group, trim(temp_str))
    end do
    call close_group(kT_group)

    ! Read energy grid
    energy_group = open_group(group_id, 'energy')
    do i = 1, n_temperature
      temp_str = trim(to_str(temps_to_read % data(i))) // "K"
      energy_dset = open_dataset(energy_group, temp_str)
      call get_shape(energy_dset, dims)
      allocate(this % grid(i) % energy(int(dims(1), 4)))
      call read_dataset(this % grid(i) % energy, energy_dset)
      call close_dataset(energy_dset)
    end do

    ! Check for 0K energy grid
    if (object_exists(energy_group, '0K')) then
      energy_dset = open_dataset(energy_group, '0K')
      call get_shape(energy_dset, dims)
      allocate(this % energy_0K(int(dims(1), 4)))
      call read_dataset(this % energy_0K, energy_dset)
      call close_dataset(energy_dset)
    end if
    call close_group(energy_group)

    ! Get MT values based on group names
    rxs_group = open_group(group_id, 'reactions')
    call get_groups(rxs_group, grp_names)
    do j = 1, size(grp_names)
      if (starts_with(grp_names(j), "reaction_")) then
        call MTs % push_back(int(str_to_int(grp_names(j)(10:12))))
      end if
    end do

    ! Read reactions
    allocate(this % reactions(MTs % size()))
    do i = 1, size(this % reactions)
      rx_group = open_group(rxs_group, 'reaction_' // trim(&
           zero_padded(MTs % data(i), 3)))

      call this % reactions(i) % from_hdf5(rx_group, temps_to_read)

      ! Check for 0K elastic scattering
      if (this % reactions(i) % MT == 2) then
        if (object_exists(rx_group, '0K')) then
          temp_group = open_group(rx_group, '0K')
          xs = open_dataset(temp_group, 'xs')
          call get_shape(xs, dims)
          allocate(this % elastic_0K(int(dims(1), 4)))
          call read_dataset(this % elastic_0K, xs)
          call close_dataset(xs)
          call close_group(temp_group)
        end if
      end if

      ! Add the reaction index to the scattering array if this is an inelastic
      ! scatter reaction
      if (MTs % data(i) /= N_FISSION .and. MTs % data(i) /= N_F .and. &
           MTs % data(i) /= N_NF .and. MTs % data(i) /= N_2NF .and. &
           MTs % data(i) /= N_3NF .and. MTs % data(i) < 200 .and. &
           MTs % data(i) /= N_LEVEL .and. MTs % data(i) /= ELASTIC) then

        call index_inelastic_scatter_vector % push_back(i)
      end if

      call close_group(rx_group)
    end do
    call close_group(rxs_group)

    ! Recast to a regular array to save space
    allocate(this % index_inelastic_scatter( &
         index_inelastic_scatter_vector % size()))
    this % index_inelastic_scatter = index_inelastic_scatter_vector % data(:)
    call index_inelastic_scatter_vector % clear()

    ! Read unresolved resonance probability tables if present
    if (object_exists(group_id, 'urr')) then
      this % urr_present = .true.
      allocate(this % urr_data(n_temperature))

      do i = 1, n_temperature
        ! Get temperature as a string
        temp_str = trim(to_str(temps_to_read % data(i))) // "K"

        ! Read probability tables for i-th temperature
        urr_group = open_group(group_id, 'urr/' // trim(temp_str))
        call this % urr_data(i) % from_hdf5(urr_group)
        call close_group(urr_group)

        ! Check for negative values
        if (any(this % urr_data(i) % prob < ZERO) .and. master) then
          call warning("Negative value(s) found on probability table &
               &for nuclide " // this % name // " at " // trim(temp_str))
        end if
      end do

      ! if the inelastic competition flag indicates that the inelastic cross
      ! section should be determined from a normal reaction cross section, we
      ! need to get the index of the reaction
      if (n_temperature > 0) then
        if (this % urr_data(1) % inelastic_flag > 0) then
          do i = 1, size(this % reactions)
            if (this % reactions(i) % MT == this % urr_data(1) % inelastic_flag) then
              this % urr_inelastic = i
            end if
          end do

          ! Abort if no corresponding inelastic reaction was found
          if (this % urr_inelastic == NONE) then
            call fatal_error("Could not find inelastic reaction specified on &
                 &unresolved resonance probability table.")
          end if
        end if
      end if
    end if

    ! Check for nu-total
    if (object_exists(group_id, 'total_nu')) then
      nu_group = open_group(group_id, 'total_nu')

      ! Read total nu data
      total_nu = open_dataset(nu_group, 'yield')
      call read_attribute(temp_str, total_nu, 'type')
      select case (temp_str)
      case ('Tabulated1D')
        allocate(Tabulated1D :: this % total_nu)
      case ('Polynomial')
        allocate(Polynomial :: this % total_nu)
      end select
      call this % total_nu % from_hdf5(total_nu)
      call close_dataset(total_nu)

      call close_group(nu_group)
    end if

    ! Read fission energy release data if present
    if (object_exists(group_id, 'fission_energy_release')) then
      fer_group = open_group(group_id, 'fission_energy_release')

      ! Check to see if this is polynomial or tabulated data
      fer_dset = open_dataset(fer_group, 'q_prompt')
      call read_attribute(temp_str, fer_dset, 'type')
      if (temp_str == 'Polynomial') then
        ! Read the prompt Q-value
        allocate(Polynomial :: this % fission_q_prompt)
        call this % fission_q_prompt % from_hdf5(fer_dset)
        call close_dataset(fer_dset)

        ! Read the recoverable energy Q-value
        allocate(Polynomial :: this % fission_q_recov)
        fer_dset = open_dataset(fer_group, 'q_recoverable')
        call this % fission_q_recov % from_hdf5(fer_dset)
        call close_dataset(fer_dset)
      else if (temp_str == 'Tabulated1D') then
        ! Read the prompt Q-value
        allocate(Tabulated1D :: this % fission_q_prompt)
        call this % fission_q_prompt % from_hdf5(fer_dset)
        call close_dataset(fer_dset)

        ! Read the recoverable energy Q-value
        allocate(Tabulated1D :: this % fission_q_recov)
        fer_dset = open_dataset(fer_group, 'q_recoverable')
        call this % fission_q_recov % from_hdf5(fer_dset)
        call close_dataset(fer_dset)
      else
        call fatal_error('Unrecognized fission energy release format.')
      end if
      call close_group(fer_group)
    end if

    ! Create derived cross section data
    call this % create_derived()

    ! Finalize with the nuclide index
    this % i_nuclide = i_nuclide

  end subroutine nuclide_from_hdf5

  subroutine nuclide_create_derived(this)
    class(Nuclide), intent(inout) :: this

    integer :: i, j, k
    integer :: t
    integer :: m
    integer :: n
    integer :: n_grid
    integer :: i_fission
    integer :: n_temperature
    type(VectorInt) :: MTs

    n_temperature = size(this % kTs)
    allocate(this % xs(n_temperature))
    this % reaction_index(:) = 0
    do i = 1, n_temperature
      ! Allocate and initialize derived cross sections
      n_grid = size(this % grid(i) % energy)
      allocate(this % xs(i) % value(4,n_grid))
      this % xs(i) % value(:,:) = ZERO
    end do

    i_fission = 0

    do i = 1, size(this % reactions)
      call MTs % push_back(this % reactions(i) % MT)
      this % reaction_index(this % reactions(i) % MT) = i

      associate (rx => this % reactions(i))
        ! Skip total inelastic level scattering, gas production cross sections
        ! (MT=200+), etc.
        if (rx % MT == N_LEVEL .or. rx % MT == N_NONELASTIC) cycle
        if (rx % MT > N_5N2P .and. rx % MT < N_P0) cycle

        ! Skip level cross sections if total is available
        if (rx % MT >= N_P0 .and. rx % MT <= N_PC .and. find(MTs, N_P) /= -1) cycle
        if (rx % MT >= N_D0 .and. rx % MT <= N_DC .and. find(MTs, N_D) /= -1) cycle
        if (rx % MT >= N_T0 .and. rx % MT <= N_TC .and. find(MTs, N_T) /= -1) cycle
        if (rx % MT >= N_3HE0 .and. rx % MT <= N_3HEC .and. find(MTs, N_3HE) /= -1) cycle
        if (rx % MT >= N_A0 .and. rx % MT <= N_AC .and. find(MTs, N_A) /= -1) cycle
        if (rx % MT >= N_2N0 .and. rx % MT <= N_2NC .and. find(MTs, N_2N) /= -1) cycle

        do t = 1, n_temperature
          j = rx % xs(t) % threshold
          n = size(rx % xs(t) % value)

          ! Add contribution to total cross section
          this % xs(t) % value(XS_TOTAL,j:j+n-1) = this % xs(t) % &
               value(XS_TOTAL,j:j+n-1) + rx % xs(t) % value

          ! Add contribution to absorption cross section
          if (is_disappearance(rx % MT)) then
            this % xs(t) % value(XS_ABSORPTION,j:j+n-1) = this % xs(t) % &
                 value(XS_ABSORPTION,j:j+n-1) + rx % xs(t) % value
          end if

          ! Information about fission reactions
          if (t == 1) then
            if (rx % MT == N_FISSION) then
              allocate(this % index_fission(1))
            elseif (rx % MT == N_F) then
              allocate(this % index_fission(PARTIAL_FISSION_MAX))
              this % has_partial_fission = .true.
            end if
          end if

          ! Add contribution to fission cross section
          if (is_fission(rx % MT)) then
            this % fissionable = .true.
            this % xs(t) % value(XS_FISSION,j:j+n-1) = this % xs(t) % &
                 value(XS_FISSION,j:j+n-1) + rx % xs(t) % value

            ! Also need to add fission cross sections to absorption
            this % xs(t) % value(XS_ABSORPTION,j:j+n-1) = this % xs(t) % &
                 value(XS_ABSORPTION,j:j+n-1) + rx % xs(t) % value

            ! If total fission reaction is present, there's no need to store the
            ! reaction cross-section since it was copied to this % fission
            if (rx % MT == N_FISSION) deallocate(rx % xs(t) % value)

            ! Keep track of this reaction for easy searching later
            if (t == 1) then
              i_fission = i_fission + 1
              this % index_fission(i_fission) = i
              this % n_fission = this % n_fission + 1

              ! <<<<<<<<<<<<<<<<<<<<<<<<<<<< REMOVE THIS <<<<<<<<<<<<<<<<<<<<<<<<<
              ! Before the secondary distribution refactor, when the angle/energy
              ! distribution was uncorrelated, no angle was actually sampled. With
              ! the refactor, an angle is always sampled for an uncorrelated
              ! distribution even when no angle distribution exists in the ACE file
              ! (isotropic is assumed). To preserve the RNG stream, we explicitly
              ! mark fission reactions so that we avoid the angle sampling.
              do k = 1, size(rx % products)
                if (rx % products(k) % particle == NEUTRON) then
                  do m = 1, size(rx % products(k) % distribution)
                    associate (aedist => rx % products(k) % distribution(m) % obj)
                      select type (aedist)
                      type is (UncorrelatedAngleEnergy)
                        aedist % fission = .true.
                      end select
                    end associate
                  end do
                end if
              end do
              ! <<<<<<<<<<<<<<<<<<<<<<<<<<<< REMOVE THIS <<<<<<<<<<<<<<<<<<<<<<<<<
            end if
          end if  ! fission
        end do  ! temperature
      end associate  ! rx
    end do ! reactions

    ! Determine number of delayed neutron precursors
    if (this % fissionable) then
      do i = 1, size(this % reactions(this % index_fission(1)) % products)
        if (this % reactions(this % index_fission(1)) % products(i) % &
             emission_mode == EMISSION_DELAYED) then
          this % n_precursor = this % n_precursor + 1
        end if
      end do
    end if

    ! Calculate nu-fission cross section
    do t = 1, n_temperature
      if (this % fissionable) then
        do i = 1, n_grid
          this % xs(t) % value(XS_NU_FISSION,i) = &
               this % nu(this % grid(t) % energy(i), EMISSION_TOTAL) * &
               this % xs(t) % value(XS_FISSION,i)
        end do
      end if
    end do
  end subroutine nuclide_create_derived

!===============================================================================
! NUCLIDE_NU is an interface to the number of fission neutrons produced
!===============================================================================

  pure function nuclide_nu(this, E, emission_mode, group) result(nu)
    class(Nuclide),    intent(in) :: this
    real(8),           intent(in) :: E
    integer,           intent(in) :: emission_mode
    integer, optional, intent(in) :: group
    real(8)                       :: nu

    integer :: i

    if (.not. this % fissionable) then
      nu = ZERO
      return
    end if

    select case (emission_mode)
    case (EMISSION_PROMPT)
      associate (product => this % reactions(this % index_fission(1)) % products(1))
        nu = product % yield % evaluate(E)
      end associate

    case (EMISSION_DELAYED)
      if (this % n_precursor > 0) then
        if (present(group) .and. group < &
             size(this % reactions(this % index_fission(1)) % products)) then
          ! If delayed group specified, determine yield immediately
          associate(p => this % reactions(this % index_fission(1)) % products(1 + group))
            nu = p % yield % evaluate(E)
          end associate

        else
          nu = ZERO

          associate (rx => this % reactions(this % index_fission(1)))
            do i = 2, size(rx % products)
              associate (product => rx % products(i))
                ! Skip any non-neutron products
                if (product % particle /= NEUTRON) exit

                ! Evaluate yield
                if (product % emission_mode == EMISSION_DELAYED) then
                  nu = nu + product % yield % evaluate(E)
                end if
              end associate
            end do
          end associate
        end if
      else
        nu = ZERO
      end if

    case (EMISSION_TOTAL)
      if (allocated(this % total_nu)) then
        nu = this % total_nu % evaluate(E)
      else
        associate (product => this % reactions(this % index_fission(1)) % products(1))
          nu = product % yield % evaluate(E)
        end associate
      end if
    end select

  end function nuclide_nu

  subroutine nuclide_init_grid(this, E_min, E_max, M)
    class(Nuclide), intent(inout) :: this
    real(8), intent(in) :: E_min           ! Minimum energy in MeV
    real(8), intent(in) :: E_max           ! Maximum energy in MeV
    integer, intent(in) :: M               ! Number of equally log-spaced bins

    integer :: i, j, k               ! Loop indices
    integer :: t                     ! temperature index
    real(8) :: spacing
    real(8), allocatable :: umesh(:) ! Equally log-spaced energy grid

    ! Determine equal-logarithmic energy spacing
    spacing = log(E_max/E_min)/M

    ! Create equally log-spaced energy grid
    allocate(umesh(0:M))
    umesh(:) = [(i*spacing, i=0, M)]

    do t = 1, size(this % grid)
      ! Allocate logarithmic mapping for nuclide
      allocate(this % grid(t) % grid_index(0:M))

      ! Determine corresponding indices in nuclide grid to energies on
      ! equal-logarithmic grid
      j = 1
      do k = 0, M
        do while (log(this % grid(t) % energy(j + 1)/E_min) <= umesh(k))
          ! Ensure that for isotopes where maxval(this % energy) << E_max
          ! that there are no out-of-bounds issues.
          if (j + 1 == size(this % grid(t) % energy)) exit
          j = j + 1
        end do
        this % grid(t) % grid_index(k) = j
      end do
    end do

  end subroutine nuclide_init_grid

!===============================================================================
! NUCLIDE_CALCULATE_XS determines microscopic cross sections for the nuclide
! at the energy of the given particle
!===============================================================================

  subroutine nuclide_calculate_xs(this, i_sab, E, i_log_union, sqrtkT, &
                                  sab_frac, micro_xs)
    class(Nuclide), intent(in) :: this ! Nuclide object
    integer, intent(in) :: i_sab       ! index into sab_tables array
    real(8), intent(in) :: E           ! energy
    integer, intent(in) :: i_log_union ! index into logarithmic mapping array or
                                       ! material union energy grid
    real(8), intent(in) :: sqrtkT      ! square root of kT, material dependent
    real(8), intent(in) :: sab_frac    ! fraction of atoms affected by S(a,b)
    type(NuclideMicroXS), intent(inout) :: micro_xs ! Cross section cache

    logical :: use_mp ! true if XS can be calculated with windowed multipole
    integer :: i_temp ! index for temperature
    integer :: i_grid ! index on nuclide energy grid
    integer :: i_low  ! lower logarithmic mapping index
    integer :: i_high ! upper logarithmic mapping index
    integer :: i_rxn  ! reaction index
    integer :: j      ! index in DEPLETION_RX
    real(8) :: f      ! interp factor on nuclide energy grid
    real(8) :: kT     ! temperature in eV
    real(8) :: sig_t, sig_a, sig_f ! Intermediate multipole variables

    ! Initialize cached cross sections to zero
    micro_xs % elastic         = CACHE_INVALID
    micro_xs % thermal         = ZERO
    micro_xs % thermal_elastic = ZERO

    ! Check to see if there is multipole data present at this energy
    use_mp = .false.
    if (this % mp_present) then
      if (E >= this % multipole % start_E .and. &
           E <= this % multipole % end_E) then
        use_mp = .true.
      end if
    end if

    ! Evaluate multipole or interpolate
    if (use_mp) then
      ! Call multipole kernel
      call multipole_eval(this % multipole, E, sqrtkT, sig_t, sig_a, sig_f)

      micro_xs % total = sig_t
      micro_xs % absorption = sig_a
      micro_xs % fission = sig_f

      if (this % fissionable) then
        micro_xs % nu_fission = sig_f * this % nu(E, EMISSION_TOTAL)
      else
        micro_xs % nu_fission = ZERO
      end if

      if (need_depletion_rx) then
        ! Initialize all reaction cross sections to zero
        micro_xs % reaction(:) = ZERO

        ! Only non-zero reaction is (n,gamma)
        micro_xs % reaction(4) = sig_a - sig_f
      end if

      ! Ensure these values are set
      ! Note, the only time either is used is in one of 4 places:
      ! 1. physics.F90 - scatter - For inelastic scatter.
      ! 2. physics.F90 - sample_fission - For partial fissions.
      ! 3. tally.F90 - score_general - For tallying on MTxxx reactions.
      ! 4. nuclide_header.F90 - calculate_urr_xs - For unresolved purposes.
      ! It is worth noting that none of these occur in the resolved
      ! resonance range, so the value here does not matter.  index_temp is
      ! set to -1 to force a segfault in case a developer messes up and tries
      ! to use it with multipole.
      micro_xs % index_temp    = -1
      micro_xs % index_grid    = 0
      micro_xs % interp_factor = ZERO

    else
      ! Find the appropriate temperature index.
      kT = sqrtkT**2
      select case (temperature_method)
      case (TEMPERATURE_NEAREST)
        i_temp = minloc(abs(this % kTs - kT), dim=1)

      case (TEMPERATURE_INTERPOLATION)
        ! Find temperatures that bound the actual temperature
        do i_temp = 1, size(this % kTs) - 1
          if (this % kTs(i_temp) <= kT .and. kT < this % kTs(i_temp + 1)) exit
        end do

        ! Randomly sample between temperature i and i+1
        f = (kT - this % kTs(i_temp)) / &
             (this % kTs(i_temp + 1) - this % kTs(i_temp))
        if (f > prn()) i_temp = i_temp + 1
      end select

      associate (grid => this % grid(i_temp), xs => this % xs(i_temp))
        ! Determine the energy grid index using a logarithmic mapping to
        ! reduce the energy range over which a binary search needs to be
        ! performed

        if (E < grid % energy(1)) then
          i_grid = 1
        elseif (E > grid % energy(size(grid % energy))) then
          i_grid = size(grid % energy) - 1
        else
          ! Determine bounding indices based on which equal log-spaced
          ! interval the energy is in
          i_low  = grid % grid_index(i_log_union)
          i_high = grid % grid_index(i_log_union + 1) + 1

          ! Perform binary search over reduced range
          i_grid = binary_search(grid % energy(i_low:i_high), &
               i_high - i_low + 1, E) + i_low - 1
        end if

        ! check for rare case where two energy points are the same
        if (grid % energy(i_grid) == grid % energy(i_grid + 1)) &
             i_grid = i_grid + 1

        ! calculate interpolation factor
        f = (E - grid % energy(i_grid)) / &
             (grid % energy(i_grid + 1) - grid % energy(i_grid))

        micro_xs % index_temp    = i_temp
        micro_xs % index_grid    = i_grid
        micro_xs % interp_factor = f

        ! Calculate microscopic nuclide total cross section
        micro_xs % total = (ONE - f) * xs % value(XS_TOTAL,i_grid) &
             + f * xs % value(XS_TOTAL,i_grid + 1)

        ! Calculate microscopic nuclide absorption cross section
        micro_xs % absorption = (ONE - f) * xs % value(XS_ABSORPTION, &
             i_grid) + f * xs % value(XS_ABSORPTION,i_grid + 1)

        if (this % fissionable) then
          ! Calculate microscopic nuclide total cross section
          micro_xs % fission = (ONE - f) * xs % value(XS_FISSION,i_grid) &
               + f * xs % value(XS_FISSION,i_grid + 1)

          ! Calculate microscopic nuclide nu-fission cross section
          micro_xs % nu_fission = (ONE - f) * xs % value(XS_NU_FISSION, &
               i_grid) + f * xs % value(XS_NU_FISSION,i_grid + 1)
        else
          micro_xs % fission         = ZERO
          micro_xs % nu_fission      = ZERO
        end if
      end associate

      ! Depletion-related reactions
      if (need_depletion_rx) then
        do j = 1, 6
          ! Initialize reaction xs to zero
          micro_xs % reaction(j) = ZERO

          ! If reaction is present and energy is greater than threshold, set
          ! the reaction xs appropriately
          i_rxn = this % reaction_index(DEPLETION_RX(j))
          if (i_rxn > 0) then
            associate (xs => this % reactions(i_rxn) % xs(i_temp))
              if (i_grid >= xs % threshold) then
                micro_xs % reaction(j) = (ONE - f) * &
                     xs % value(i_grid - xs % threshold + 1) + &
                     f * xs % value(i_grid - xs % threshold + 2)
              end if
            end associate
          end if
        end do
      end if

    end if

    ! Initialize sab treatment to false
    micro_xs % index_sab = NONE
    micro_xs % sab_frac = ZERO

    ! Initialize URR probability table treatment to false
    micro_xs % use_ptable = .false.

    ! If there is S(a,b) data for this nuclide, we need to set the sab_scatter
    ! and sab_elastic cross sections and correct the total and elastic cross
    ! sections.

    if (i_sab > 0) then
      call calculate_sab_xs(this, i_sab, E, sqrtkT, sab_frac, micro_xs)
    end if

    ! If the particle is in the unresolved resonance range and there are
    ! probability tables, we need to determine cross sections from the table

    if (urr_ptables_on .and. this % urr_present .and. .not. use_mp) then
      if (E > this % urr_data(i_temp) % energy(1) .and. E < this % &
           urr_data(i_temp) % energy(this % urr_data(i_temp) % n_energy)) then
        call calculate_urr_xs(this, i_temp, E, micro_xs)
      end if
    end if

    micro_xs % last_E = E
    micro_xs % last_sqrtkT = sqrtkT

  end subroutine nuclide_calculate_xs

!===============================================================================
! NUCLIDE_CALCULATE_ELASTIC_XS precalculates the free atom elastic scattering
! cross section. Normally it is not needed until a collision actually occurs in
! a material. However, in the thermal and unresolved resonance regions, we have
! to calculate it early to adjust the total cross section correctly.
!===============================================================================

  subroutine nuclide_calculate_elastic_xs(this, micro_xs)
    class(Nuclide),       intent(in)    :: this
    type(NuclideMicroXS), intent(inout) :: micro_xs ! Cross section cache

    integer :: i_temp
    integer :: i_grid
    real(8) :: f

    ! Get temperature index, grid index, and interpolation factor
    i_temp =  micro_xs % index_temp
    i_grid =  micro_xs % index_grid
    f      =  micro_xs % interp_factor

    if (i_temp > 0) then
      associate (xs => this % reactions(1) % xs(i_temp) % value)
        micro_xs % elastic = (ONE - f) * xs(i_grid) + f * xs(i_grid + 1)
      end associate
    else
      ! For multipole, elastic is total - absorption
      micro_xs % elastic = micro_xs % total - micro_xs % absorption
    end if
  end subroutine nuclide_calculate_elastic_xs

!===============================================================================
! CALCULATE_SAB_XS determines the elastic and inelastic scattering
! cross-sections in the thermal energy range. These cross sections replace a
! fraction of whatever data were taken from the normal Nuclide table.
!===============================================================================

  subroutine calculate_sab_xs(this, i_sab, E, sqrtkT, sab_frac, micro_xs)
    class(Nuclide), intent(in) :: this ! Nuclide object
    integer, intent(in) :: i_sab     ! index into sab_tables array
    real(8), intent(in) :: E         ! energy
    real(8), intent(in) :: sqrtkT    ! temperature
    real(8), intent(in) :: sab_frac  ! fraction of atoms affected by S(a,b)
    type(NuclideMicroXS), intent(inout) :: micro_xs ! Cross section cache

    integer :: i_temp    ! temperature index
    real(8) :: inelastic ! S(a,b) inelastic cross section
    real(8) :: elastic   ! S(a,b) elastic cross section

    ! Set flag that S(a,b) treatment should be used for scattering
    micro_xs % index_sab = i_sab

    ! Calculate the S(a,b) cross section
    call sab_tables(i_sab) % calculate_xs(E, sqrtkT, i_temp, elastic, inelastic)

    ! Store the S(a,b) cross sections.
    micro_xs % thermal = sab_frac * (elastic + inelastic)
    micro_xs % thermal_elastic = sab_frac * elastic

    ! Calculate free atom elastic cross section
    call this % calculate_elastic_xs(micro_xs)

    ! Correct total and elastic cross sections
    micro_xs % total = micro_xs % total + micro_xs % thermal - &
         sab_frac *  micro_xs % elastic
    micro_xs % elastic = micro_xs % thermal + (ONE - sab_frac) * &
         micro_xs % elastic

    ! Save temperature index and thermal fraction
    micro_xs % index_temp_sab = i_temp
    micro_xs % sab_frac = sab_frac

  end subroutine calculate_sab_xs

!===============================================================================
! MULTIPOLE_EVAL evaluates the windowed multipole equations for cross
! sections in the resolved resonance regions
!===============================================================================

  subroutine multipole_eval(multipole, E, sqrtkT, sig_t, sig_a, sig_f)
    type(MultipoleArray), intent(in) :: multipole ! The windowed multipole
                                                  !  object to process.
    real(8), intent(in)              :: E         ! The energy at which to
                                                  !  evaluate the cross section
    real(8), intent(in)              :: sqrtkT    ! The temperature in the form
                                                  !  sqrt(kT), at which
                                                  !  to evaluate the XS.
    real(8), intent(out)             :: sig_t     ! Total cross section
    real(8), intent(out)             :: sig_a     ! Absorption cross section
    real(8), intent(out)             :: sig_f     ! Fission cross section
    complex(8) :: psi_chi  ! The value of the psi-chi function for the
                           !  asymptotic form
    complex(8) :: c_temp   ! complex temporary variable
    complex(8) :: w_val    ! The faddeeva function evaluated at Z
    complex(8) :: Z        ! sqrt(atomic weight ratio / kT) * (sqrt(E) - pole)
    complex(8) :: sig_t_factor(multipole % num_l)
    real(8) :: broadened_polynomials(multipole % fit_order + 1)
    real(8) :: sqrtE       ! sqrt(E), eV
    real(8) :: invE        ! 1/E, eV
    real(8) :: dopp        ! sqrt(atomic weight ratio / kT) = 1 / (2 sqrt(xi))
    real(8) :: temp        ! real temporary value
    integer :: i_pole      ! index of pole
    integer :: i_poly      ! index of curvefit
    integer :: i_window    ! index of window
    integer :: startw      ! window start pointer (for poles)
    integer :: endw        ! window end pointer

    ! ==========================================================================
    ! Bookkeeping

    ! Define some frequently used variables.
    sqrtE = sqrt(E)
    invE = ONE / E

    ! Locate us.
    i_window = floor((sqrtE - sqrt(multipole % start_E)) / multipole % spacing &
         + ONE)
    startw = multipole % w_start(i_window)
    endw = multipole % w_end(i_window)

    ! Fill in factors.
    if (startw <= endw) then
      call compute_sig_t_factor(multipole, sqrtE, sig_t_factor)
    end if

    ! Initialize the ouptut cross sections.
    sig_t = ZERO
    sig_a = ZERO
    sig_f = ZERO

    ! ==========================================================================
    ! Add the contribution from the curvefit polynomial.

    if (sqrtkT /= ZERO .and. multipole % broaden_poly(i_window) == 1) then
      ! Broaden the curvefit.
      dopp = multipole % sqrtAWR / sqrtkT
      call broaden_wmp_polynomials(E, dopp, multipole % fit_order + 1, &
           broadened_polynomials)
      do i_poly = 1, multipole % fit_order+1
        sig_t = sig_t + multipole % curvefit(FIT_T, i_poly, i_window) &
             * broadened_polynomials(i_poly)
        sig_a = sig_a + multipole % curvefit(FIT_A, i_poly, i_window) &
             * broadened_polynomials(i_poly)
        if (multipole % fissionable) then
          sig_f = sig_f + multipole % curvefit(FIT_F, i_poly, i_window) &
               * broadened_polynomials(i_poly)
        end if
      end do
    else ! Evaluate as if it were a polynomial
      temp = invE
      do i_poly = 1, multipole % fit_order+1
        sig_t = sig_t + multipole % curvefit(FIT_T, i_poly, i_window) * temp
        sig_a = sig_a + multipole % curvefit(FIT_A, i_poly, i_window) * temp
        if (multipole % fissionable) then
          sig_f = sig_f + multipole % curvefit(FIT_F, i_poly, i_window) * temp
        end if
        temp = temp * sqrtE
      end do
    end if

    ! ==========================================================================
    ! Add the contribution from the poles in this window.

    if (sqrtkT == ZERO) then
      ! If at 0K, use asymptotic form.
      do i_pole = startw, endw
        psi_chi = -ONEI / (multipole % data(MP_EA, i_pole) - sqrtE)
        c_temp = psi_chi / E
        if (multipole % formalism == FORM_MLBW) then
          sig_t = sig_t + real(multipole % data(MLBW_RT, i_pole) * c_temp * &
                               sig_t_factor(multipole % l_value(i_pole))) &
                        + real(multipole % data(MLBW_RX, i_pole) * c_temp)
          sig_a = sig_a + real(multipole % data(MLBW_RA, i_pole) * c_temp)
          if (multipole % fissionable) then
            sig_f = sig_f + real(multipole % data(MLBW_RF, i_pole) * c_temp)
          end if
        else if (multipole % formalism == FORM_RM) then
          sig_t = sig_t + real(multipole % data(RM_RT, i_pole) * c_temp * &
                               sig_t_factor(multipole % l_value(i_pole)))
          sig_a = sig_a + real(multipole % data(RM_RA, i_pole) * c_temp)
          if (multipole % fissionable) then
            sig_f = sig_f + real(multipole % data(RM_RF, i_pole) * c_temp)
          end if
        end if
      end do
    else
      ! At temperature, use Faddeeva function-based form.
      dopp = multipole % sqrtAWR / sqrtkT
      if (endw >= startw) then
        do i_pole = startw, endw
          Z = (sqrtE - multipole % data(MP_EA, i_pole)) * dopp
          w_val = faddeeva(Z) * dopp * invE * SQRT_PI
          if (multipole % formalism == FORM_MLBW) then
            sig_t = sig_t + real((multipole % data(MLBW_RT, i_pole) * &
                            sig_t_factor(multipole % l_value(i_pole)) + &
                            multipole % data(MLBW_RX, i_pole)) * w_val)
            sig_a = sig_a + real(multipole % data(MLBW_RA, i_pole) * w_val)
            if (multipole % fissionable) then
              sig_f = sig_f + real(multipole % data(MLBW_RF, i_pole) * w_val)
            end if
          else if (multipole % formalism == FORM_RM) then
            sig_t = sig_t + real(multipole % data(RM_RT, i_pole) * w_val * &
                                 sig_t_factor(multipole % l_value(i_pole)))
            sig_a = sig_a + real(multipole % data(RM_RA, i_pole) * w_val)
            if (multipole % fissionable) then
              sig_f = sig_f + real(multipole % data(RM_RF, i_pole) * w_val)
            end if
          end if
        end do
      end if
    end if
  end subroutine multipole_eval

!===============================================================================
! MULTIPOLE_DERIV_EVAL evaluates the windowed multipole equations for the
! derivative of cross sections in the resolved resonance regions with respect to
! temperature.
!===============================================================================

  subroutine multipole_deriv_eval(multipole, E, sqrtkT, sig_t, sig_a, sig_f)
    type(MultipoleArray), intent(in) :: multipole ! The windowed multipole
                                                  !  object to process.
    real(8), intent(in)              :: E         ! The energy at which to
                                                  !  evaluate the cross section
    real(8), intent(in)              :: sqrtkT    ! The temperature in the form
                                                  !  sqrt(kT), at which to
                                                  !  evaluate the XS.
    real(8), intent(out)             :: sig_t     ! Total cross section
    real(8), intent(out)             :: sig_a     ! Absorption cross section
    real(8), intent(out)             :: sig_f     ! Fission cross section
    complex(8) :: w_val    ! The faddeeva function evaluated at Z
    complex(8) :: Z        ! sqrt(atomic weight ratio / kT) * (sqrt(E) - pole)
    complex(8) :: sig_t_factor(multipole % num_l)
    real(8) :: sqrtE       ! sqrt(E), eV
    real(8) :: invE        ! 1/E, eV
    real(8) :: dopp        ! sqrt(atomic weight ratio / kT)
    integer :: i_pole      ! index of pole
    integer :: i_window    ! index of window
    integer :: startw      ! window start pointer (for poles)
    integer :: endw        ! window end pointer
    real(8) :: T

    ! ==========================================================================
    ! Bookkeeping

    ! Define some frequently used variables.
    sqrtE = sqrt(E)
    invE = ONE / E
    T = sqrtkT**2 / K_BOLTZMANN

    if (sqrtkT == ZERO) call fatal_error("Windowed multipole temperature &
         &derivatives are not implemented for 0 Kelvin cross sections.")

    ! Locate us
    i_window = floor((sqrtE - sqrt(multipole % start_E)) / multipole % spacing &
         + ONE)
    startw = multipole % w_start(i_window)
    endw = multipole % w_end(i_window)

    ! Fill in factors.
    if (startw <= endw) then
      call compute_sig_t_factor(multipole, sqrtE, sig_t_factor)
    end if

    ! Initialize the ouptut cross sections.
    sig_t = ZERO
    sig_a = ZERO
    sig_f = ZERO

    ! TODO Polynomials: Some of the curvefit polynomials Doppler broaden so
    ! rigorously we should be computing the derivative of those.  But in
    ! practice, those derivatives are only large at very low energy and they
    ! have no effect on reactor calculations.

    ! ==========================================================================
    ! Add the contribution from the poles in this window.

    dopp = multipole % sqrtAWR / sqrtkT
    if (endw >= startw) then
      do i_pole = startw, endw
        Z = (sqrtE - multipole % data(MP_EA, i_pole)) * dopp
        w_val = -invE * SQRT_PI * HALF * w_derivative(Z, 2)
        if (multipole % formalism == FORM_MLBW) then
          sig_t = sig_t + real((multipole % data(MLBW_RT, i_pole) * &
                          sig_t_factor(multipole%l_value(i_pole)) + &
                          multipole % data(MLBW_RX, i_pole)) * w_val)
          sig_a = sig_a + real(multipole % data(MLBW_RA, i_pole) * w_val)
          if (multipole % fissionable) then
            sig_f = sig_f + real(multipole % data(MLBW_RF, i_pole) * w_val)
          end if
        else if (multipole % formalism == FORM_RM) then
          sig_t = sig_t + real(multipole % data(RM_RT, i_pole) * w_val * &
                               sig_t_factor(multipole % l_value(i_pole)))
          sig_a = sig_a + real(multipole % data(RM_RA, i_pole) * w_val)
          if (multipole % fissionable) then
            sig_f = sig_f + real(multipole % data(RM_RF, i_pole) * w_val)
          end if
        end if
      end do
      sig_t = -HALF*multipole % sqrtAWR / sqrt(K_BOLTZMANN) * T**(-1.5) * sig_t
      sig_a = -HALF*multipole % sqrtAWR / sqrt(K_BOLTZMANN) * T**(-1.5) * sig_a
      sig_f = -HALF*multipole % sqrtAWR / sqrt(K_BOLTZMANN) * T**(-1.5) * sig_f
    end if
  end subroutine multipole_deriv_eval

!===============================================================================
! COMPUTE_SIG_T_FACTOR calculates the sig_t_factor, a factor inside of the sig_t
! equation not present in the sig_a and sig_f equations.
!===============================================================================

  subroutine compute_sig_t_factor(multipole, sqrtE, sig_t_factor)
    type(MultipoleArray), intent(in)  :: multipole
    real(8),              intent(in)  :: sqrtE
    complex(8),           intent(out) :: sig_t_factor(multipole % num_l)

    integer :: iL
    real(8) :: twophi(multipole % num_l)
    real(8) :: arg

    do iL = 1, multipole % num_l
      twophi(iL) = multipole % pseudo_k0RS(iL) * sqrtE
      if (iL == 2) then
        twophi(iL) = twophi(iL) - atan(twophi(iL))
      else if (iL == 3) then
        arg = 3.0_8 * twophi(iL) / (3.0_8 - twophi(iL)**2)
        twophi(iL) = twophi(iL) - atan(arg)
      else if (iL == 4) then
        arg = twophi(iL) * (15.0_8 - twophi(iL)**2) &
             / (15.0_8 - 6.0_8 * twophi(iL)**2)
        twophi(iL) = twophi(iL) - atan(arg)
      end if
    end do

    twophi = 2.0_8 * twophi
    sig_t_factor = cmplx(cos(twophi), -sin(twophi), KIND=8)
  end subroutine compute_sig_t_factor

!===============================================================================
! 0K_ELASTIC_XS determines the microscopic 0K elastic cross section
! for a given nuclide at the trial relative energy used in resonance scattering
!===============================================================================

  pure function elastic_xs_0K(E, nuc) result(xs_out)
    real(8),       intent(in) :: E      ! trial energy
    type(Nuclide), intent(in) :: nuc    ! target nuclide at temperature
    real(8)                   :: xs_out ! 0K xs at trial energy

    integer :: i_grid ! index on nuclide energy grid
    integer :: n_grid
    real(8) :: f      ! interp factor on nuclide energy grid

    ! Determine index on nuclide energy grid
    n_grid = size(nuc % energy_0K)
    if (E < nuc % energy_0K(1)) then
      i_grid = 1
    elseif (E > nuc % energy_0K(n_grid)) then
      i_grid = n_grid - 1
    else
      i_grid = binary_search(nuc % energy_0K, n_grid, E)
    end if

    ! check for rare case where two energy points are the same
    if (nuc % energy_0K(i_grid) == nuc % energy_0K(i_grid+1)) then
      i_grid = i_grid + 1
    end if

    ! calculate interpolation factor
    f = (E - nuc % energy_0K(i_grid)) &
         & / (nuc % energy_0K(i_grid + 1) - nuc % energy_0K(i_grid))

    ! Calculate microscopic nuclide elastic cross section
    xs_out = (ONE - f) * nuc % elastic_0K(i_grid) &
         & + f * nuc % elastic_0K(i_grid + 1)

  end function elastic_xs_0K

!===============================================================================
! CALCULATE_URR_XS determines cross sections in the unresolved resonance range
! from probability tables
!===============================================================================

  subroutine calculate_urr_xs(this, i_temp, E, micro_xs)
    class(Nuclide), intent(in)        :: this     ! Nuclide object
    integer, intent(in)               :: i_temp   ! temperature index
    real(8), intent(in)               :: E        ! energy
    type(NuclideMicroXS), intent(inout) :: micro_xs ! Cross section cache

    integer :: i_energy     ! index for energy
    integer :: i_low        ! band index at lower bounding energy
    integer :: i_up         ! band index at upper bounding energy
    real(8) :: f            ! interpolation factor
    real(8) :: r            ! pseudo-random number
    real(8) :: elastic      ! elastic cross section
    real(8) :: capture      ! (n,gamma) cross section
    real(8) :: fission      ! fission cross section
    real(8) :: inelastic    ! inelastic cross section

    micro_xs % use_ptable = .true.

    associate (urr => this % urr_data(i_temp))
      ! determine energy table
      i_energy = 1
      do
        if (E < urr % energy(i_energy + 1)) exit
        i_energy = i_energy + 1
      end do

      ! determine interpolation factor on table
      f = (E - urr % energy(i_energy)) / &
           (urr % energy(i_energy + 1) - urr % energy(i_energy))

      ! sample probability table using the cumulative distribution

      ! Random numbers for xs calculation are sampled from a separated stream.
      ! This guarantees the randomness and, at the same time, makes sure we reuse
      ! random number for the same nuclide at different temperatures, therefore
      ! preserving correlation of temperature in probability tables.
      call prn_set_stream(STREAM_URR_PTABLE)
      r = future_prn(int(this % i_nuclide, 8))
      call prn_set_stream(STREAM_TRACKING)

      i_low = 1
      do
        if (urr % prob(i_energy, URR_CUM_PROB, i_low) > r) exit
        i_low = i_low + 1
      end do
      i_up = 1
      do
        if (urr % prob(i_energy + 1, URR_CUM_PROB, i_up) > r) exit
        i_up = i_up + 1
      end do

      ! determine elastic, fission, and capture cross sections from probability
      ! table
      if (urr % interp == LINEAR_LINEAR) then
        elastic = (ONE - f) * urr % prob(i_energy, URR_ELASTIC, i_low) + &
             f * urr % prob(i_energy + 1, URR_ELASTIC, i_up)
        fission = (ONE - f) * urr % prob(i_energy, URR_FISSION, i_low) + &
             f * urr % prob(i_energy + 1, URR_FISSION, i_up)
        capture = (ONE - f) * urr % prob(i_energy, URR_N_GAMMA, i_low) + &
             f * urr % prob(i_energy + 1, URR_N_GAMMA, i_up)
      elseif (urr % interp == LOG_LOG) then
        ! Get logarithmic interpolation factor
        f = log(E / urr % energy(i_energy)) / &
             log(urr % energy(i_energy + 1) / urr % energy(i_energy))

        ! Calculate elastic cross section/factor
        elastic = ZERO
        if (urr % prob(i_energy, URR_ELASTIC, i_low) > ZERO .and. &
             urr % prob(i_energy + 1, URR_ELASTIC, i_up) > ZERO) then
          elastic = exp((ONE - f) * log(urr % prob(i_energy, URR_ELASTIC, &
               i_low)) + f * log(urr % prob(i_energy + 1, URR_ELASTIC, &
               i_up)))
        end if

        ! Calculate fission cross section/factor
        fission = ZERO
        if (urr % prob(i_energy, URR_FISSION, i_low) > ZERO .and. &
             urr % prob(i_energy + 1, URR_FISSION, i_up) > ZERO) then
          fission = exp((ONE - f) * log(urr % prob(i_energy, URR_FISSION, &
               i_low)) + f * log(urr % prob(i_energy + 1, URR_FISSION, &
               i_up)))
        end if

        ! Calculate capture cross section/factor
        capture = ZERO
        if (urr % prob(i_energy, URR_N_GAMMA, i_low) > ZERO .and. &
             urr % prob(i_energy + 1, URR_N_GAMMA, i_up) > ZERO) then
          capture = exp((ONE - f) * log(urr % prob(i_energy, URR_N_GAMMA, &
               i_low)) + f * log(urr % prob(i_energy + 1, URR_N_GAMMA, &
               i_up)))
        end if
      end if

      ! Determine treatment of inelastic scattering
      inelastic = ZERO
      if (urr % inelastic_flag > 0) then
        ! Get index on energy grid and interpolation factor
        i_energy = micro_xs % index_grid
        f = micro_xs % interp_factor

        ! Determine inelastic scattering cross section
        associate (xs => this % reactions(this % urr_inelastic) % xs(i_temp))
          if (i_energy >= xs % threshold) then
            inelastic = (ONE - f) * xs % value(i_energy - xs % threshold + 1) + &
                 f * xs % value(i_energy - xs % threshold + 2)
          end if
        end associate
      end if

      ! Multiply by smooth cross-section if needed
      if (urr % multiply_smooth) then
        call this % calculate_elastic_xs(micro_xs)
        elastic = elastic * micro_xs % elastic
        capture = capture * (micro_xs % absorption - micro_xs % fission)
        fission = fission * micro_xs % fission
      end if

      ! Check for negative values
      if (elastic < ZERO) elastic = ZERO
      if (fission < ZERO) fission = ZERO
      if (capture < ZERO) capture = ZERO

      ! Set elastic, absorption, fission, and total cross sections. Note that the
      ! total cross section is calculated as sum of partials rather than using the
      ! table-provided value
      micro_xs % elastic = elastic
      micro_xs % absorption = capture + fission
      micro_xs % fission = fission
      micro_xs % total = elastic + inelastic + capture + fission

      ! Determine nu-fission cross section
      if (this % fissionable) then
        micro_xs % nu_fission = this % nu(E, EMISSION_TOTAL) * &
             micro_xs % fission
      end if
    end associate

  end subroutine calculate_urr_xs

!===============================================================================
! CHECK_DATA_VERSION checks for the right version of nuclear data within HDF5
! files
!===============================================================================

  subroutine check_data_version(file_id)
    integer(HID_T), intent(in) :: file_id

    integer, allocatable :: version(:)

    if (attribute_exists(file_id, 'version')) then
      call read_attribute(version, file_id, 'version')
      if (version(1) /= HDF5_VERSION(1)) then
        call fatal_error("HDF5 data format uses version " // trim(to_str(&
             version(1))) // "." // trim(to_str(version(2))) // " whereas &
             &your installation of OpenMC expects version " // trim(to_str(&
             HDF5_VERSION(1))) // ".x data.")
      end if
    else
      call fatal_error("HDF5 data does not indicate a version. Your &
           &installation of OpenMC expects version " // trim(to_str(&
           HDF5_VERSION(1))) // ".x data.")
    end if
  end subroutine check_data_version

!===============================================================================
! FREE_MEMORY_NUCLIDE deallocates global arrays defined in this module
!===============================================================================

  subroutine free_memory_nuclide()
    integer :: i

    ! Deallocate cross section data, listings, and cache
    if (allocated(nuclides)) then
      ! First call the clear routines
      do i = 1, size(nuclides)
        call nuclides(i) % clear()
      end do
      deallocate(nuclides)
    end if
    n_nuclides = 0

    if (allocated(libraries)) deallocate(libraries)

    call nuclide_dict % clear()
    call library_dict % clear()

  end subroutine free_memory_nuclide

!===============================================================================
!                               C API FUNCTIONS
!===============================================================================

  function openmc_get_nuclide_index(name, index) result(err) bind(C)
    ! Return the index in the nuclides array of a nuclide with a given name
    character(kind=C_CHAR), intent(in) :: name(*)
    integer(C_INT), intent(out) :: index
    integer(C_INT) :: err

    character(:), allocatable :: name_

    ! Copy array of C_CHARs to normal Fortran string
    name_ = to_f_string(name)

    if (allocated(nuclides)) then
      if (nuclide_dict % has(to_lower(name_))) then
        index = nuclide_dict % get(to_lower(name_))
        err = 0
      else
        err = E_DATA
        call set_errmsg("No nuclide named '" // trim(name_) // &
             "' has been loaded.")
      end if
    else
      err = E_ALLOCATE
      call set_errmsg("Memory for nuclides has not been allocated.")
    end if
  end function openmc_get_nuclide_index


  function openmc_load_nuclide(name) result(err) bind(C)
    ! Load a nuclide from the cross section library
    character(kind=C_CHAR), intent(in) :: name(*)
    integer(C_INT) :: err

    integer :: i_library
    integer :: n
    integer(HID_T) :: file_id
    integer(HID_T) :: group_id
    character(:), allocatable :: name_
    real(8) :: minmax(2) = [ZERO, INFINITY]
    type(VectorReal) :: temperature
    type(Nuclide), allocatable :: new_nuclides(:)

    ! Copy array of C_CHARs to normal Fortran string
    name_ = to_f_string(name)

    err = 0
    if (.not. nuclide_dict % has(to_lower(name_))) then
      if (library_dict % has(to_lower(name_))) then
        ! allocate extra space in nuclides array
        n = n_nuclides
        allocate(new_nuclides(n + 1))
        new_nuclides(1:n) = nuclides(:)
        call move_alloc(FROM=new_nuclides, TO=nuclides)
        n = n + 1

        i_library = library_dict % get(to_lower(name_))

        ! Open file and make sure version is sufficient
        file_id = file_open(libraries(i_library) % path, 'r')
        call check_data_version(file_id)

        ! Read nuclide data from HDF5
        group_id = open_group(file_id, name_)
        call nuclides(n) % from_hdf5(group_id, temperature, &
             temperature_method, temperature_tolerance, minmax, &
             master, n)
        call close_group(group_id)
        call file_close(file_id)

        ! Add entry to nuclide dictionary
        call nuclide_dict % set(to_lower(name_), n)
        n_nuclides = n

        ! Assign resonant scattering data
        if (res_scat_on) call nuclides(n) % assign_0K_elastic_scattering()

        ! Initialize nuclide grid
        call nuclides(n) % init_grid(energy_min_neutron, &
             energy_max_neutron, n_log_bins)
      else
        err = E_DATA
        call set_errmsg("Nuclide '" // trim(name_) // "' is not present &
             &in library.")
      end if
    end if

  end function openmc_load_nuclide


  function openmc_nuclide_name(index, name) result(err) bind(C)
    ! Return the name of a nuclide with a given index
    integer(C_INT), value, intent(in) :: index
    type(c_ptr), intent(out) :: name
    integer(C_INT) :: err

    character(C_CHAR), pointer :: name_

    err = E_UNASSIGNED
    if (allocated(nuclides)) then
      if (index >= 1 .and. index <= size(nuclides)) then
        name_ => nuclides(index) % name(1:1)
        name = C_LOC(name_)
        err = 0
      else
        err = E_OUT_OF_BOUNDS
        call set_errmsg("Index in nuclides array is out of bounds.")
      end if
    else
      err = E_ALLOCATE
      call set_errmsg("Memory for nuclides has not been allocated yet.")
    end if
  end function openmc_nuclide_name

end module nuclide_header
