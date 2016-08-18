module nuclide_header

  use, intrinsic :: ISO_FORTRAN_ENV
  use, intrinsic :: ISO_C_BINDING

  use hdf5, only: HID_T, HSIZE_T, SIZE_T, h5iget_name_f, h5gget_info_f, &
                  h5lget_name_by_idx_f, H5_INDEX_NAME_F, H5_ITER_INC_F
  use h5lt, only: h5ltpath_valid_f

  use constants
  use dict_header, only: DictIntInt
  use endf,        only: reaction_name, is_fission, is_disappearance
  use endf_header, only: Function1D, Polynomial, Tabulated1D
  use error,       only: fatal_error, warning
  use hdf5_interface, only: read_attribute, open_group, close_group, &
       open_dataset, read_dataset, close_dataset, get_shape
  use list_header, only: ListInt
  use math,        only: evaluate_legendre
  use multipole_header, only: MultipoleArray
  use product_header, only: AngleEnergyContainer
  use reaction_header, only: Reaction
  use secondary_uncorrelated, only: UncorrelatedAngleEnergy
  use stl_vector,  only: VectorInt
  use string
  use urr_header, only: UrrData
  use xml_interface

  implicit none

!===============================================================================
! Nuclide contains the base nuclidic data for a nuclide described as needed
! for continuous-energy neutron transport.
!===============================================================================

  type :: Nuclide
    ! Nuclide meta-data
    character(20) :: name    ! name of nuclide, e.g. U235.71c
    integer       :: zaid    ! Z and A identifier, e.g. 92235
    integer       :: metastable ! metastable state
    real(8)       :: awr     ! Atomic Weight Ratio
    real(8)       :: kT      ! temperature in MeV (k*T)

    ! Fission information
    logical :: fissionable = .false.  ! nuclide is fissionable?

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
    logical :: has_partial_fission = .false. ! nuclide has partial fission reactions?
    integer :: n_fission = 0                 ! # of fission reactions
    integer :: n_precursor = 0               ! # of delayed neutron precursors
    integer, allocatable :: index_fission(:) ! indices in reactions
    class(Function1D), allocatable :: total_nu

    ! Unresolved resonance data
    logical                :: urr_present = .false.
    integer                :: urr_inelastic
    type(UrrData), pointer :: urr_data => null()

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
    procedure :: print => nuclide_print
    procedure :: from_hdf5 => nuclide_from_hdf5
    procedure :: nu    => nuclide_nu
    procedure, private :: create_derived => nuclide_create_derived
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

    ! Information for Doppler broadening
    real(8) :: last_sqrtkT = ZERO  ! Last temperature in sqrt(Boltzmann
                                   ! constant * temperature (MeV))
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

    if (associated(this % urr_data)) deallocate(this % urr_data)
    if (associated(this % multipole)) deallocate(this % multipole)

  end subroutine nuclide_clear

  subroutine nuclide_from_hdf5(this, group_id)
    class(Nuclide), intent(inout) :: this
    integer(HID_T),   intent(in)    :: group_id

    integer :: i
    integer :: Z
    integer :: A
    integer :: storage_type
    integer :: max_corder
    integer :: n_links
    integer :: hdf5_err
    integer(HID_T) :: urr_group, nu_group
    integer(HID_T) :: energy_dset
    integer(HID_T) :: rxs_group
    integer(HID_T) :: rx_group
    integer(HID_T) :: total_nu
    integer(HID_T) :: fer_group                 ! fission_energy_release group
    integer(HID_T) :: fer_dset
    integer(SIZE_T) :: name_len, name_file_len
    integer(HSIZE_T) :: j
    integer(HSIZE_T) :: dims(1)
    character(MAX_WORD_LEN) :: temp
    type(VectorInt) :: MTs
    logical :: exists

    ! Get name of nuclide from group
    name_len = len(this % name)
    call h5iget_name_f(group_id, this % name, name_len, name_file_len, hdf5_err)

    ! Get rid of leading '/'
    this % name = trim(this % name(2:))

    call read_attribute(Z, group_id, 'Z')
    call read_attribute(A, group_id, 'A')
    call read_attribute(this % metastable, group_id, 'metastable')
    this % zaid = 1000*Z + A + 400*this % metastable
    call read_attribute(this % awr, group_id, 'atomic_weight_ratio')
    call read_attribute(this % kT, group_id, 'temperature')

    ! Read energy grid
    energy_dset = open_dataset(group_id, 'energy')
    call get_shape(energy_dset, dims)
    this % n_grid = int(dims(1), 4)
    allocate(this % energy(this % n_grid))
    call read_dataset(this % energy, energy_dset)
    call close_dataset(energy_dset)

    ! Get MT values based on group names
    rxs_group = open_group(group_id, 'reactions')
    call h5gget_info_f(rxs_group, storage_type, n_links, max_corder, hdf5_err)
    do j = 0, n_links - 1
      call h5lget_name_by_idx_f(rxs_group, ".", H5_INDEX_NAME_F, H5_ITER_INC_F, &
           j, temp, hdf5_err, name_len)
      if (starts_with(temp, "reaction_")) then
        call MTs % push_back(int(str_to_int(temp(10:12))))
      end if
    end do

    ! Read reactions
    allocate(this % reactions(MTs % size()))
    do i = 1, size(this % reactions)
      rx_group = open_group(rxs_group, 'reaction_' // trim(&
           zero_padded(MTs % data(i), 3)))
      call this % reactions(i) % from_hdf5(rx_group)
      call close_group(rx_group)
    end do
    call close_group(rxs_group)

    ! Read unresolved resonance probability tables if present
    call h5ltpath_valid_f(group_id, 'urr', .true., exists, hdf5_err)
    if (exists) then
      this % urr_present = .true.
      allocate(this % urr_data)
      urr_group = open_group(group_id, 'urr')
      call this % urr_data % from_hdf5(urr_group)

      ! if the inelastic competition flag indicates that the inelastic cross
      ! section should be determined from a normal reaction cross section, we
      ! need to get the index of the reaction
      if (this % urr_data % inelastic_flag > 0) then
        do i = 1, size(this % reactions)
          if (this % reactions(i) % MT == this % urr_data % inelastic_flag) then
            this % urr_inelastic = i
          end if
        end do

        ! Abort if no corresponding inelastic reaction was found
        if (this % urr_inelastic == NONE) then
          call fatal_error("Could not find inelastic reaction specified on &
               &unresolved resonance probability table.")
        end if
      end if

      ! Check for negative values
      if (any(this % urr_data % prob < ZERO)) then
        call warning("Negative value(s) found on probability table &
             &for nuclide " // this % name)
      end if
    end if

    ! Check for nu-total
    call h5ltpath_valid_f(group_id, 'total_nu', .true., exists, hdf5_err)
    if (exists) then
      nu_group = open_group(group_id, 'total_nu')

      ! Read total nu data
      total_nu = open_dataset(nu_group, 'yield')
      call read_attribute(temp, total_nu, 'type')
      select case (temp)
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
    call h5ltpath_valid_f(group_id, 'fission_energy_release', .true., exists, &
                          hdf5_err)
    if (exists) then
      fer_group = open_group(group_id, 'fission_energy_release')

      ! Check to see if this is polynomial or tabulated data
      fer_dset = open_dataset(fer_group, 'q_prompt')
      call read_attribute(temp, fer_dset, 'type')
      if (temp == 'Polynomial') then
        ! Read the prompt Q-value
        allocate(Polynomial :: this % fission_q_prompt)
        call this % fission_q_prompt % from_hdf5(fer_dset)
        call close_dataset(fer_dset)

        ! Read the recoverable energy Q-value
        allocate(Polynomial :: this % fission_q_recov)
        fer_dset = open_dataset(fer_group, 'q_recoverable')
        call this % fission_q_recov % from_hdf5(fer_dset)
        call close_dataset(fer_dset)
      else if (temp == 'Tabulated1D') then
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

    integer :: i
    integer :: j
    integer :: k
    integer :: m
    integer :: n
    integer :: i_fission
    type(ListInt) :: MTs

    ! Allocate and initialize derived cross sections
    allocate(this % total(this % n_grid))
    allocate(this % elastic(this % n_grid))
    allocate(this % fission(this % n_grid))
    allocate(this % nu_fission(this % n_grid))
    allocate(this % absorption(this % n_grid))
    this % total(:) = ZERO
    this % elastic(:) = ZERO
    this % fission(:) = ZERO
    this % nu_fission(:) = ZERO
    this % absorption(:) = ZERO

    i_fission = 0

    do i = 1, size(this % reactions)
      call MTs % append(this % reactions(i) % MT)
      call this % reaction_index % add_key(this % reactions(i) % MT, i)

      associate (rx => this % reactions(i))
        j = rx % threshold
        n = size(rx % sigma)

        ! Skip total inelastic level scattering, gas production cross sections
        ! (MT=200+), etc.
        if (rx % MT == N_LEVEL .or. rx % MT == N_NONELASTIC) cycle
        if (rx % MT > N_5N2P .and. rx % MT < N_P0) cycle

        ! Skip level cross sections if total is available
        if (rx % MT >= N_P0 .and. rx % MT <= N_PC .and. MTs % contains(N_P)) cycle
        if (rx % MT >= N_D0 .and. rx % MT <= N_DC .and. MTs % contains(N_D)) cycle
        if (rx % MT >= N_T0 .and. rx % MT <= N_TC .and. MTs % contains(N_T)) cycle
        if (rx % MT >= N_3HE0 .and. rx % MT <= N_3HEC .and. MTs % contains(N_3HE)) cycle
        if (rx % MT >= N_A0 .and. rx % MT <= N_AC .and. MTs % contains(N_A)) cycle
        if (rx % MT >= N_2N0 .and. rx % MT <= N_2NC .and. MTs % contains(N_2N)) cycle

        ! Copy elastic
        if (rx % MT == ELASTIC) this % elastic(:) = rx % sigma

        ! Add contribution to total cross section
        this % total(j:j+n-1) = this % total(j:j+n-1) + rx % sigma

        ! Add contribution to absorption cross section
        if (is_disappearance(rx % MT)) then
          this % absorption(j:j+n-1) = this % absorption(j:j+n-1) + rx % sigma
        end if

        ! Information about fission reactions
        if (rx % MT == N_FISSION) then
          allocate(this % index_fission(1))
        elseif (rx % MT == N_F) then
          allocate(this % index_fission(PARTIAL_FISSION_MAX))
          this % has_partial_fission = .true.
        end if

        ! Add contribution to fission cross section
        if (is_fission(rx % MT)) then
          this % fissionable = .true.
          this % fission(j:j+n-1) = this % fission(j:j+n-1) + rx % sigma

          ! Also need to add fission cross sections to absorption
          this % absorption(j:j+n-1) = this % absorption(j:j+n-1) + rx % sigma

          ! If total fission reaction is present, there's no need to store the
          ! reaction cross-section since it was copied to this % fission
          if (rx % MT == N_FISSION) deallocate(rx % sigma)

          ! Keep track of this reaction for easy searching later
          i_fission = i_fission + 1
          this % index_fission(i_fission) = i
          this % n_fission = this % n_fission + 1

          ! <<<<<<<<<<<<<<<<<<<<<<<<<<<< REMOVE THIS <<<<<<<<<<<<<<<<<<<<<<<<<<<
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
          ! <<<<<<<<<<<<<<<<<<<<<<<<<<<< REMOVE THIS <<<<<<<<<<<<<<<<<<<<<<<<<<<
        end if
      end associate
    end do

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
    if (this % fissionable) then
      do i = 1, size(this % energy)
        this % nu_fission(i) = this % nu(this % energy(i), EMISSION_TOTAL) * &
             this % fission(i)
      end do
    else
      this % nu_fission(:) = ZERO
    end if

    ! Clear MTs set
    call MTs % clear()
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


!===============================================================================
! NUCLIDE*_PRINT displays information about a continuous-energy neutron
! cross_section table and its reactions and secondary angle/energy distributions
!===============================================================================

  subroutine nuclide_print(this, unit)
    class(Nuclide), intent(in) :: this
    integer, intent(in), optional :: unit

    integer :: i                 ! loop index over nuclides
    integer :: unit_             ! unit to write to
    integer :: size_xs           ! memory used for cross-sections (bytes)
    integer :: size_urr          ! memory used for probability tables (bytes)

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
    write(unit_,*) '  # of reactions = ' // trim(to_str(size(this % reactions)))

    ! Information on each reaction
    write(unit_,*) '  Reaction     Q-value  COM    IE'
    do i = 1, size(this % reactions)
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
      associate(urr => this % urr_data)
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
      end associate
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
