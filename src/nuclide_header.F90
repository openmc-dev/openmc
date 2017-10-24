module nuclide_header

  use, intrinsic :: ISO_FORTRAN_ENV
  use, intrinsic :: ISO_C_BINDING

  use hdf5, only: HID_T, HSIZE_T, SIZE_T

  use algorithm, only: sort, find
  use constants
  use dict_header, only: DictIntInt, DictCharInt
  use endf,        only: reaction_name, is_fission, is_disappearance
  use endf_header, only: Function1D, Polynomial, Tabulated1D
  use error
  use hdf5_interface
  use list_header, only: ListInt
  use math,        only: evaluate_legendre
  use message_passing
  use multipole_header, only: MultipoleArray
  use product_header, only: AngleEnergyContainer
  use reaction_header, only: Reaction
  use secondary_uncorrelated, only: UncorrelatedAngleEnergy
  use settings
  use stl_vector,  only: VectorInt, VectorReal
  use string
  use urr_header, only: UrrData

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
    type(DictIntInt) :: reaction_index ! map MT values to index in reactions
                                       ! array; used at tally-time

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
  end type Nuclide

!===============================================================================
! NUCLIDEMICROXS contains cached microscopic cross sections for a particular
! nuclide at the current energy
!===============================================================================

  type NuclideMicroXS
    ! Microscopic cross sections in barns
    real(8) :: total
    real(8) :: elastic          ! If sab_frac is not 1 or 0, then this value is
                                !   averaged over bound and non-bound nuclei
    real(8) :: absorption
    real(8) :: fission
    real(8) :: nu_fission
    real(8) :: thermal          ! Bound thermal elastic & inelastic scattering
    real(8) :: thermal_elastic  ! Bound thermal elastic scattering

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
                               minmax, master)
    class(Nuclide),   intent(inout) :: this
    integer(HID_T),   intent(in)    :: group_id
    type(VectorReal), intent(in)    :: temperature ! list of desired temperatures
    integer,          intent(inout) :: method
    real(8),          intent(in)    :: tolerance
    real(8),          intent(in)    :: minmax(2)  ! range of temperatures
    logical,          intent(in)    :: master     ! if this is the master proc

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
      call this % reaction_index % set(this % reactions(i) % MT, i)

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
             master)
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
