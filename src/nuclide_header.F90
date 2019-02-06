module nuclide_header

  use, intrinsic :: ISO_FORTRAN_ENV
  use, intrinsic :: ISO_C_BINDING

  use algorithm,              only: sort, find, binary_search
  use constants
  use endf,                   only: is_fission, is_disappearance
  use endf_header,            only: Function1D, Polynomial, Tabulated1D
  use error
  use hdf5_interface
  use message_passing
  use reaction_header,        only: Reaction
  use settings
  use stl_vector,             only: VectorInt, VectorReal
  use string

  implicit none

!===============================================================================
! Nuclide contains the base nuclidic data for a nuclide described as needed
! for continuous-energy neutron transport.
!===============================================================================

  type EnergyGrid
    real(8), allocatable :: energy(:)     ! energy values corresponding to xs
  end type EnergyGrid

  ! Positions for first dimension of Nuclide % xs
  integer, parameter :: &
       XS_TOTAL       = 1, &
       XS_ABSORPTION  = 2, &
       XS_FISSION     = 3, &
       XS_NU_FISSION  = 4, &
       XS_PHOTON_PROD = 5

  ! The array within SumXS is of shape (5, n_energy) where the first dimension
  ! corresponds to the following values: 1) total, 2) absorption (MT > 100), 3)
  ! fission, 4) neutron production, 5) photon production
  type SumXS
    real(8), allocatable :: value(:,:)
  end type SumXS

  type :: Nuclide
    ! Nuclide meta-data
    character(20) :: name    ! name of nuclide, e.g. U235
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

    ! Fission information
    integer, allocatable :: index_fission(:) ! indices in reactions

    ! Reactions
    type(Reaction), allocatable :: reactions(:)

    ! Array that maps MT values to index in reactions; used at tally-time. Note
    ! that ENDF-102 does not assign any MT values above 891.
    integer :: reaction_index(891)

    type(C_PTR) :: ptr

  contains
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
  real(8), parameter :: CACHE_INVALID = -1

  type, bind(C) :: NuclideMicroXS
    ! Microscopic cross sections in barns
    real(C_DOUBLE) :: total
    real(C_DOUBLE) :: absorption       ! absorption (disappearance)
    real(C_DOUBLE) :: fission          ! fission
    real(C_DOUBLE) :: nu_fission       ! neutron production from fission

    real(C_DOUBLE) :: elastic          ! If sab_frac is not 1 or 0, then this value is
                                !   averaged over bound and non-bound nuclei
    real(C_DOUBLE) :: thermal          ! Bound thermal elastic & inelastic scattering
    real(C_DOUBLE) :: thermal_elastic  ! Bound thermal elastic scattering
    real(C_DOUBLE) :: photon_prod      ! microscopic photon production xs

    ! Cross sections for depletion reactions (note that these are not stored in
    ! macroscopic cache)
    real(C_DOUBLE) :: reaction(size(DEPLETION_RX))

    ! Indicies and factors needed to compute cross sections from the data tables
    integer(C_INT) :: index_grid        ! Index on nuclide energy grid
    integer(C_INT) :: index_temp        ! Temperature index for nuclide
    real(C_DOUBLE) :: interp_factor     ! Interpolation factor on nuc. energy grid
    integer(C_INT) :: index_sab = NONE  ! Index in sab_tables
    integer(C_INT) :: index_temp_sab    ! Temperature index for sab_tables
    real(C_DOUBLE) :: sab_frac          ! Fraction of atoms affected by S(a,b)
    logical(C_BOOL) :: use_ptable        ! In URR range with probability tables?

    ! Energy and temperature last used to evaluate these cross sections.  If
    ! these values have changed, then the cross sections must be re-evaluated.
    real(C_DOUBLE) :: last_E = ZERO       ! Last evaluated energy
    real(C_DOUBLE) :: last_sqrtkT = ZERO  ! Last temperature in sqrt(Boltzmann
                                   !   constant * temperature (eV))
  end type NuclideMicroXS

!===============================================================================
! MATERIALMACROXS contains cached macroscopic cross sections for the material a
! particle is traveling through
!===============================================================================

  type, bind(C) :: MaterialMacroXS
    real(C_DOUBLE) :: total         ! macroscopic total xs
    real(C_DOUBLE) :: absorption    ! macroscopic absorption xs
    real(C_DOUBLE) :: fission       ! macroscopic fission xs
    real(C_DOUBLE) :: nu_fission    ! macroscopic production xs
    real(C_DOUBLE) :: photon_prod   ! macroscopic photon production xs

    ! Photon cross sections
    real(C_DOUBLE) :: coherent        ! macroscopic coherent xs
    real(C_DOUBLE) :: incoherent      ! macroscopic incoherent xs
    real(C_DOUBLE) :: photoelectric   ! macroscopic photoelectric xs
    real(C_DOUBLE) :: pair_production ! macroscopic pair production xs
  end type MaterialMacroXS

  ! Nuclear data for each nuclide
  type(Nuclide), allocatable, target :: nuclides(:)
  integer(C_INT), bind(C) :: n_nuclides

  ! Cross section caches
  type(NuclideMicroXS), allocatable, target :: micro_xs(:)  ! Cache for each nuclide
  type(MaterialMacroXS), bind(C)            :: material_xs  ! Cache for current material
!$omp threadprivate(micro_xs, material_xs)

  interface
    function library_present_c(type, name) result(b) bind(C, name='library_present')
      import C_INT, C_CHAR, C_BOOL
      integer(C_INT), value :: type
      character(kind=C_CHAR), intent(in) :: name(*)
      logical(C_BOOL) :: b
    end function

    function library_path_c(type, name) result(path) bind(C, name='library_path')
      import C_INT, C_CHAR, C_PTR
      integer(C_INT), value :: type
      character(kind=C_CHAR), intent(in) :: name(*)
      type(C_PTR) :: path
    end function

    subroutine multipole_deriv_eval(ptr, E, sqrtkT, sig_s, sig_a, sig_f) bind(C)
      import C_PTR, C_DOUBLE
      type(C_PTR), value :: ptr
      real(C_DOUBLE), value :: E
      real(C_DOUBLE), value :: sqrtkT
      real(C_DOUBLE), intent(out) :: sig_s
      real(C_DOUBLE), intent(out) :: sig_a
      real(C_DOUBLE), intent(out) :: sig_f
    end subroutine

    function multipole_in_range(ptr, E) result(b) bind(C)
      import C_PTR, C_DOUBLE, C_BOOL
      type(C_PTR), value :: ptr
      real(C_DOUBLE), value :: E
      logical(C_BOOL) :: b
    end function

    function nuclide_awr(i_nuc) result(awr) bind(C)
      import C_INT, C_DOUBLE
      integer(C_INT), value :: i_nuc
      real(C_DOUBLE) :: awr
    end function

    function nuclide_fission_q_prompt(ptr, E) result(q) bind(C)
      import C_PTR, C_DOUBLE
      type(C_PTR), value :: ptr
      real(C_DOUBLE), value :: E
      real(C_DOUBLE) :: q
    end function

    function nuclide_fission_q_recov(ptr, E) result(q) bind(C)
      import C_PTR, C_DOUBLE
      type(C_PTR), value :: ptr
      real(C_DOUBLE), value :: E
      real(C_DOUBLE) :: q
    end function

    function nuclide_map_get(name) result(idx) bind(C)
      import C_CHAR, C_INT
      character(kind=C_CHAR), intent(in) :: name(*)
      integer(C_INT) :: idx
    end function
  end interface

contains

  function library_path(type, name) result(path)
    integer, intent(in) :: type
    character(len=*), intent(in) :: name
    character(MAX_FILE_LEN) :: path

    type(C_PTR) :: ptr
    character(kind=C_CHAR), pointer :: string(:)

    ptr = library_path_c(type, to_c_string(name))
    call c_f_pointer(ptr, string, [255])
    path = to_f_string(string)
  end function

  function library_present(type, name) result(b)
    integer, intent(in) :: type
    character(len=*), intent(in) :: name
    logical :: b

    b = library_present_c(type, to_c_string(name))
  end function

  function micro_xs_ptr() result(ptr) bind(C)
    type(C_PTR) :: ptr
    ptr = C_LOC(micro_xs(1))
  end function

  subroutine nuclide_from_hdf5(group_id, ptr, temps, n, i_nuclide) bind(C)
    integer(HID_T), value :: group_id
    type(C_PTR), value :: ptr
    type(C_PTR), value :: temps
    integer(C_INT), value :: n
    integer(C_INT), value :: i_nuclide

    real(C_DOUBLE), pointer :: temperature(:) ! list of desired temperatures

    integer :: i
    integer :: i_closest
    integer :: n_temperature
    integer(HID_T) :: energy_group, energy_dset
    integer(HID_T) :: kT_group
    integer(HID_T) :: rxs_group
    integer(HID_T) :: rx_group
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

    ! Get array passed
    call c_f_pointer(temps, temperature, [n])

    associate (this => nuclides(i_nuclide))
    this % ptr = ptr

    ! Get name of nuclide from group
    this % name = get_name(group_id)

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
    if (size(temps_available) == 1 .and. temperature_method == TEMPERATURE_INTERPOLATION) then
      if (master) then
        call warning("Cross sections for " // trim(this % name) // " are only &
             &available at one temperature. Reverting to nearest temperature &
             &method.")
      end if
      temperature_method = TEMPERATURE_NEAREST
    end if

    ! Determine actual temperatures to read -- start by checking whether a
    ! temperature range was given, in which case all temperatures in the range
    ! are loaded irrespective of what temperatures actually appear in the model
    if (temperature_range(2) > ZERO) then
      do i = 1, size(temps_available)
        temp_actual = temps_available(i)
        if (temperature_range(1) <= temp_actual .and. temp_actual <= temperature_range(2)) then
          call temps_to_read % push_back(nint(temp_actual))
        end if
      end do
    end if

    select case (temperature_method)
    case (TEMPERATURE_NEAREST)
      ! Find nearest temperatures
      do i = 1, n
        temp_desired = temperature(i)
        i_closest = minloc(abs(temps_available - temp_desired), dim=1)
        temp_actual = temps_available(i_closest)
        if (abs(temp_actual - temp_desired) < temperature_tolerance) then
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
      TEMP_LOOP: do i = 1, n
        temp_desired = temperature(i)

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

      ! Set pointer for each reaction
      call this % reactions(i) % init(this % ptr, i)

      call close_group(rx_group)
    end do
    call close_group(rxs_group)

    ! Create derived cross section data
    call this % create_derived()

    ! Finalize with the nuclide index
    this % i_nuclide = i_nuclide
    end associate

  end subroutine nuclide_from_hdf5

  subroutine nuclide_create_derived(this)
    class(Nuclide), intent(inout) :: this

    integer :: i, j, k, l
    integer :: t
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
      allocate(this % xs(i) % value(5,n_grid))
      this % xs(i) % value(:,:) = ZERO
    end do

    i_fission = 0

    do i = 1, size(this % reactions)
      call MTs % push_back(this % reactions(i) % MT)
      this % reaction_index(this % reactions(i) % MT) = i

      associate (rx => this % reactions(i))
        do t = 1, n_temperature
          j = rx % xs_threshold(t)
          n = rx % xs_size(t)

          ! Calculate photon production cross section
          do k = 1, rx % products_size()
            if (rx % product_particle(k) == PHOTON) then
              do l = 1, n
                this % xs(t) % value(XS_PHOTON_PROD,l+j-1) = &
                     this % xs(t) % value(XS_PHOTON_PROD,l+j-1) + &
                     rx % xs(t, l) * rx % product_yield(k, &
                     this % grid(t) % energy(l+j-1))
              end do
            end if
          end do

          ! Skip gas production cross sections (MT=200+), etc.
          if (rx % MT > N_5N2P .and. rx % MT < N_P0) cycle

          ! Skip any reaction that has been marked as redundant
          if (rx % redundant) cycle

          ! Add contribution to total cross section
          do k = j, j + n - 1
            this % xs(t) % value(XS_TOTAL,k) = this % xs(t) % &
                 value(XS_TOTAL,k) + rx % xs(t, k - j + 1)
          end do

          ! Add contribution to absorption cross section
          if (is_disappearance(rx % MT)) then
            do k = j, j + n - 1
              this % xs(t) % value(XS_ABSORPTION,k) = this % xs(t) % &
                   value(XS_ABSORPTION,k) + rx % xs(t, k - j + 1)
            end do
          end if

          ! Information about fission reactions
          if (t == 1) then
            if (rx % MT == N_FISSION) then
              allocate(this % index_fission(1))
            elseif (rx % MT == N_F) then
              allocate(this % index_fission(PARTIAL_FISSION_MAX))
            end if
          end if

          ! Add contribution to fission cross section
          if (is_fission(rx % MT)) then
            this % fissionable = .true.
            do k = j, j + n - 1
              this % xs(t) % value(XS_FISSION,k) = this % xs(t) % &
                   value(XS_FISSION,k) + rx % xs(t, k - j + 1)

              ! Also need to add fission cross sections to absorption
              this % xs(t) % value(XS_ABSORPTION,k) = this % xs(t) % &
                   value(XS_ABSORPTION,k) + rx % xs(t, k - j + 1)
            end do

            ! Keep track of this reaction for easy searching later
            if (t == 1) then
              i_fission = i_fission + 1
              this % index_fission(i_fission) = i
            end if
          end if  ! fission
        end do  ! temperature
      end associate  ! rx
    end do ! reactions

    ! Calculate nu-fission cross section
    do t = 1, n_temperature
      if (this % fissionable) then
        n_grid = size(this % grid(t) % energy)
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

    interface
      pure function nuclide_nu_c(ptr, E, emission_mode, group) result(nu) bind(C)
        import C_PTR, C_DOUBLE, C_INT
        type(C_PTR), value, intent(in) :: ptr
        real(C_DOUBLE), value, intent(in) :: E
        integer(C_INT), value, intent(in) :: emission_mode
        integer(C_INT), value, intent(in) :: group
        real(C_DOUBLE) :: nu
      end function
    end interface

    if (present(group)) then
      nu = nuclide_nu_c(this % ptr, E, emission_mode, group)
    else
      nu = nuclide_nu_c(this % ptr, E, emission_mode, 0)
    end if
  end function nuclide_nu

  subroutine nuclide_init_grid(this)
    class(nuclide), intent(inout) :: this
    interface
      subroutine nuclide_init_grid_c(ptr) bind(C)
        import C_PTR
        type(C_PTR), value :: ptr
      end subroutine
    end interface
    call nuclide_init_grid_c(this % ptr)
  end subroutine

!===============================================================================
! NUCLIDE_CALCULATE_XS determines microscopic cross sections for the nuclide
! at the energy of the given particle
!===============================================================================

  subroutine nuclide_calculate_xs(this, i_sab, E, i_log_union, sqrtkT, &
                                  sab_frac)
    class(Nuclide), intent(in) :: this ! Nuclide object
    integer, intent(in) :: i_sab       ! index into sab_tables array
    real(8), intent(in) :: E           ! energy
    integer, intent(in) :: i_log_union ! index into logarithmic mapping array or
                                       ! material union energy grid
    real(8), intent(in) :: sqrtkT      ! square root of kT, material dependent
    real(8), intent(in) :: sab_frac    ! fraction of atoms affected by S(a,b)

    interface
      subroutine nuclide_calculate_xs_c(ptr, i_sab, E, i_log_union, sqrtkT, sab_frac) bind(C)
        import C_PTR, C_INT, C_DOUBLE
        type(C_PTR), value :: ptr
        integer(C_INT), value :: i_sab
        real(C_DOUBLE), value :: E
        integer(C_INT), value :: i_log_union
        real(C_DOUBLE), value :: sqrtkT
        real(C_DOUBLE), value :: sab_frac
      end subroutine
    end interface

    call nuclide_calculate_xs_c(this % ptr, i_sab - 1, E, i_log_union, sqrtkT, sab_frac)
  end subroutine nuclide_calculate_xs

!===============================================================================
! NUCLIDE_CALCULATE_ELASTIC_XS precalculates the free atom elastic scattering
! cross section. Normally it is not needed until a collision actually occurs in
! a material. However, in the thermal and unresolved resonance regions, we have
! to calculate it early to adjust the total cross section correctly.
!===============================================================================

  subroutine nuclide_calculate_elastic_xs(this)
    class(Nuclide),       intent(in)    :: this
    interface
      subroutine nuclide_calculate_elastic_xs_c(ptr) bind(C)
        import C_PTR
        type(C_PTR), value :: ptr
      end subroutine
    end interface
    call nuclide_calculate_elastic_xs_c(this % ptr)
  end subroutine nuclide_calculate_elastic_xs

!===============================================================================
! FREE_MEMORY_NUCLIDE deallocates global arrays defined in this module
!===============================================================================

  subroutine free_memory_nuclide()
    interface
      subroutine library_clear() bind(C)
      end subroutine

      subroutine nuclides_clear() bind(C)
      end subroutine
    end interface

    ! Deallocate cross section data, listings, and cache
    if (allocated(nuclides)) then
      ! First call the clear routines
      deallocate(nuclides)
      call nuclides_clear()
    end if
    n_nuclides = 0

    call library_clear()

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

    err = 0
    if (allocated(nuclides)) then
      index = nuclide_map_get(name)
      if (index == -1) then
        err = E_DATA
        call set_errmsg("No nuclide named '" // trim(name_) // &
             "' has been loaded.")
      end if
    else
      err = E_ALLOCATE
      call set_errmsg("Memory for nuclides has not been allocated.")
    end if
  end function openmc_get_nuclide_index


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

  subroutine extend_nuclides() bind(C)
    type(Nuclide), allocatable :: new_nuclides(:)

    ! allocate extra space in nuclides array
    allocate(new_nuclides(n_nuclides + 1))
    new_nuclides(1:n_nuclides) = nuclides(:)
    call move_alloc(FROM=new_nuclides, TO=nuclides)
    n_nuclides = n_nuclides + 1
  end subroutine

end module nuclide_header
