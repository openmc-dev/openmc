module nuclide_header

  use, intrinsic :: ISO_FORTRAN_ENV
  use, intrinsic :: ISO_C_BINDING

  use hdf5, only: HID_T, HSIZE_T, SIZE_T

  use algorithm, only: sort, find
  use constants
  use dict_header, only: DictIntInt
  use endf,        only: reaction_name, is_fission, is_disappearance
  use endf_header, only: Function1D, Polynomial, Tabulated1D
  use error,       only: fatal_error, warning
  use hdf5_interface, only: read_attribute, open_group, close_group, &
       open_dataset, read_dataset, close_dataset, get_shape, get_datasets, &
       object_exists, get_name, get_groups
  use list_header, only: ListInt
  use math,        only: evaluate_legendre
  use multipole_header, only: MultipoleArray
  use product_header, only: AngleEnergyContainer
  use reaction_header, only: Reaction
  use secondary_uncorrelated, only: UncorrelatedAngleEnergy
  use stl_vector,  only: VectorInt, VectorReal
  use string
  use urr_header, only: UrrData
  use xml_interface

  implicit none

!===============================================================================
! Nuclide contains the base nuclidic data for a nuclide described as needed
! for continuous-energy neutron transport.
!===============================================================================

  type EnergyGrid
    integer, allocatable :: grid_index(:) ! log grid mapping indices
    real(8), allocatable :: energy(:)     ! energy values corresponding to xs
  end type EnergyGrid

  type SumXS
    real(8), allocatable :: total(:)      ! total cross section
    real(8), allocatable :: elastic(:)    ! elastic scattering
    real(8), allocatable :: fission(:)    ! fission
    real(8), allocatable :: nu_fission(:) ! neutron production
    real(8), allocatable :: absorption(:) ! absorption (MT > 100)
    real(8), allocatable :: heating(:)    ! heating
  end type SumXS

  type :: Nuclide
    ! Nuclide meta-data
    character(20) :: name    ! name of nuclide, e.g. U235.71c
    integer       :: Z       ! atomic number
    integer       :: A       ! mass number
    integer       :: metastable ! metastable state
    real(8)       :: awr     ! Atomic Weight Ratio
    real(8), allocatable :: kTs(:) ! temperature in eV (k*T)

    ! Fission information
    logical :: fissionable = .false.  ! nuclide is fissionable?

    ! Energy grid for each temperature
    type(EnergyGrid), allocatable :: grid(:)

    ! Microscopic cross sections
    type(SumXS), allocatable :: sum_xs(:)

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
    type(DictIntInt) :: reaction_index ! map MT values to index in reactions
                                       ! array; used at tally-time

    ! Fission energy release
    class(Function1D), allocatable :: fission_q_prompt ! prompt neutrons, gammas
    class(Function1D), allocatable :: fission_q_recov  ! neutrons, gammas, betas

  contains
    procedure :: clear => nuclide_clear
    procedure :: from_hdf5 => nuclide_from_hdf5
    procedure :: nu    => nuclide_nu
    procedure, private :: create_derived => nuclide_create_derived
  end type Nuclide

!===============================================================================
! NUCLIDE0K temporarily contains all 0K cross section data and other parameters
! needed to treat resonance scattering before transferring them to Nuclide
!===============================================================================

  type Nuclide0K
    character(10) :: nuclide             ! name of nuclide, e.g. U238
    character(16) :: scheme = 'ares'     ! target velocity sampling scheme
    real(8)       :: E_min = 0.01_8   ! lower cutoff energy for res scattering
    real(8)       :: E_max = 1000.0_8 ! upper cutoff energy for res scattering
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
    integer :: index_temp_sab     ! temperature index for sab_tables
    real(8) :: elastic_sab        ! microscopic elastic scattering on S(a,b) table

    ! Information for URR probability table use
    logical :: use_ptable  ! in URR range with probability tables?

    ! Information for Doppler broadening
    real(8) :: last_sqrtkT = ZERO  ! Last temperature in sqrt(Boltzmann
                                   ! constant * temperature (eV))
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
! LIBRARY contains data read from a cross_sections.xml file
!===============================================================================

  type Library
    integer :: type
    character(MAX_WORD_LEN), allocatable :: materials(:)
    character(MAX_FILE_LEN) :: path
  end type Library

  contains

!===============================================================================
! NUCLIDE_CLEAR resets and deallocates data in Nuclide
!===============================================================================

  subroutine nuclide_clear(this)
    class(Nuclide), intent(inout) :: this ! The Nuclide object to clear

    if (associated(this % multipole)) deallocate(this % multipole)

  end subroutine nuclide_clear

  subroutine nuclide_from_hdf5(this, group_id, temperature, method, tolerance, &
                               master)
    class(Nuclide),   intent(inout) :: this
    integer(HID_T),   intent(in)    :: group_id
    type(VectorReal), intent(in)   :: temperature ! list of desired temperatures
    integer,          intent(inout) :: method
    real(8),          intent(in)    :: tolerance
    logical,          intent(in)    :: master     ! if this is the master proc

    integer :: i
    integer :: i_closest
    integer :: n_temperature
    integer(HID_T) :: urr_group, nu_group
    integer(HID_T) :: energy_group, energy_dset
    integer(HID_T) :: kT_group
    integer(HID_T) :: rxs_group
    integer(HID_T) :: rx_group
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
    if (size(temps_available) == 1 .and. &
         method == TEMPERATURE_INTERPOLATION) then
      call warning("Cross sections for " // trim(this % name) // " are only &
           &available at one temperature. Reverting to nearest temperature &
           &method.")
      method = TEMPERATURE_NEAREST
    end if

    ! Determine actual temperatures to read
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

    do i = 1, n_temperature
      ! Get temperature as a string
      temp_str = trim(to_str(temps_to_read % data(i))) // "K"

      ! Read exact temperature value
      call read_dataset(this % kTs(i), kT_group, trim(temp_str))

      ! Read energy grid
      energy_group = open_group(group_id, 'energy')
      energy_dset = open_dataset(energy_group, temp_str)
      call get_shape(energy_dset, dims)
      allocate(this % grid(i) % energy(int(dims(1), 4)))
      call read_dataset(this % grid(i) % energy, energy_dset)
      call close_dataset(energy_dset)
      call close_group(energy_group)
    end do

    call close_group(kT_group)

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
      call close_group(rx_group)
    end do
    call close_group(rxs_group)

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
        if (any(this % urr_data(i) % prob < ZERO)) then
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
    allocate(this % sum_xs(n_temperature))

    do i = 1, n_temperature
      ! Allocate and initialize derived cross sections
      n_grid = size(this % grid(i) % energy)
      allocate(this % sum_xs(i) % total(n_grid))
      allocate(this % sum_xs(i) % elastic(n_grid))
      allocate(this % sum_xs(i) % fission(n_grid))
      allocate(this % sum_xs(i) % nu_fission(n_grid))
      allocate(this % sum_xs(i) % absorption(n_grid))
      this % sum_xs(i) % total(:) = ZERO
      this % sum_xs(i) % elastic(:) = ZERO
      this % sum_xs(i) % fission(:) = ZERO
      this % sum_xs(i) % nu_fission(:) = ZERO
      this % sum_xs(i) % absorption(:) = ZERO
    end do

    i_fission = 0

    do i = 1, size(this % reactions)
      call MTs % push_back(this % reactions(i) % MT)
      call this % reaction_index % add_key(this % reactions(i) % MT, i)

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

          ! Copy elastic
          if (rx % MT == ELASTIC) this % sum_xs(t) % elastic(:) = rx % xs(t) % value

          ! Add contribution to total cross section
          this % sum_xs(t) % total(j:j+n-1) = this % sum_xs(t) % total(j:j+n-1) + &
               rx % xs(t) % value

          ! Add contribution to absorption cross section
          if (is_disappearance(rx % MT)) then
            this % sum_xs(t) % absorption(j:j+n-1) = this % sum_xs(t) % &
                 absorption(j:j+n-1) + rx % xs(t) % value
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
            this % sum_xs(t) % fission(j:j+n-1) = this % sum_xs(t) % &
                 fission(j:j+n-1) + rx % xs(t) % value

            ! Also need to add fission cross sections to absorption
            this % sum_xs(t) % absorption(j:j+n-1) = this % sum_xs(t) % &
                 absorption(j:j+n-1) + rx % xs(t) % value

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
        do i = 1, size(this % sum_xs(t) % fission)
          this % sum_xs(t) % nu_fission(i) = this % nu(this % grid(t) % energy(i), &
               EMISSION_TOTAL) * this % sum_xs(t) % fission(i)
        end do
      else
        this % sum_xs(t) % nu_fission(:) = ZERO
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
        if (present(group)) then
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

end module nuclide_header
