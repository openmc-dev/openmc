module sab_header

  use, intrinsic :: ISO_C_BINDING
  use, intrinsic :: ISO_FORTRAN_ENV

  use algorithm, only: find, sort, binary_search
  use constants
  use dict_header, only: DictIntInt, DictCharInt
  use distribution_univariate, only: Tabular
  use error,       only: warning, fatal_error
  use hdf5, only: HID_T, HSIZE_T, SIZE_T
  use hdf5_interface, only: read_attribute, get_shape, open_group, close_group, &
       open_dataset, read_dataset, close_dataset, get_datasets, object_exists, &
       get_name
  use random_lcg,  only: prn
  use secondary_correlated, only: CorrelatedAngleEnergy
  use settings
  use stl_vector, only: VectorInt, VectorReal
  use string, only: to_str, str_to_int

  implicit none

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

  type SabData
    ! threshold for S(a,b) treatment (usually ~4 eV)
    real(8) :: threshold_inelastic
    real(8) :: threshold_elastic = ZERO

    ! Inelastic scattering data
    integer :: n_inelastic_e_in  ! # of incoming E for inelastic
    integer :: n_inelastic_e_out ! # of outgoing E for inelastic
    integer :: n_inelastic_mu    ! # of outgoing angles for inelastic
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
  end type SabData

  type SAlphaBeta
    character(150) :: name     ! name of table, e.g. lwtr.10t
    real(8)        :: awr      ! weight of nucleus in neutron masses
    real(8), allocatable :: kTs(:)  ! temperatures in eV (k*T)
    character(10), allocatable :: nuclides(:) ! List of valid nuclides
    integer :: secondary_mode    ! secondary mode (equal/skewed/continuous)

    ! cross sections and distributions at each temperature
    type(SabData), allocatable :: data(:)
  contains
    procedure :: from_hdf5 => salphabeta_from_hdf5
    procedure :: calculate_xs => sab_calculate_xs
  end type SAlphaBeta

  ! S(a,b) tables
  type(SAlphaBeta), allocatable, target :: sab_tables(:)
  integer(C_INT), bind(C) :: n_sab_tables
  type(DictCharInt) :: sab_dict

contains

  subroutine salphabeta_from_hdf5(this, group_id, temperature, method, &
       tolerance, minmax)
    class(SAlphaBeta), intent(inout) :: this
    integer(HID_T),    intent(in)    :: group_id
    type(VectorReal),  intent(in)    :: temperature ! list of temperatures
    integer,           intent(in)    :: method
    real(8),           intent(in)    :: tolerance
    real(8),           intent(in)    :: minmax(2)

    integer :: i, j
    integer :: t
    integer :: n_energy, n_energy_out, n_mu
    integer :: i_closest
    integer :: n_temperature
    integer(SIZE_T) :: name_len
    integer(HID_T) :: T_group
    integer(HID_T) :: elastic_group
    integer(HID_T) :: inelastic_group
    integer(HID_T) :: dset_id
    integer(HID_T) :: kT_group
    integer(HSIZE_T) :: dims2(2)
    integer(HSIZE_T) :: dims3(3)
    real(8), allocatable :: temp(:,:)
    character(20) :: type
    type(CorrelatedAngleEnergy) :: correlated_dist

    character(MAX_WORD_LEN) :: temp_str
    character(MAX_WORD_LEN), allocatable :: dset_names(:)
    real(8), allocatable :: temps_available(:) ! temperatures available
    real(8) :: temp_desired
    real(8) :: temp_actual
    type(VectorInt) :: temps_to_read

    ! Get name of table from group
    name_len = len(this % name)
    this % name = get_name(group_id, name_len)

    ! Get rid of leading '/'
    this % name = trim(this % name(2:))

    call read_attribute(this % awr, group_id, 'atomic_weight_ratio')
    call read_attribute(this % nuclides, group_id, 'nuclides')
    call read_attribute(type, group_id, 'secondary_mode')
    select case (type)
    case ('equal')
      this % secondary_mode = SAB_SECONDARY_EQUAL
    case ('skewed')
      this % secondary_mode = SAB_SECONDARY_SKEWED
    case ('continuous')
      this % secondary_mode = SAB_SECONDARY_CONT
    end select

    ! Read temperatures
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
      ! Determine actual temperatures to read
      do i = 1, temperature % size()
        temp_desired = temperature % data(i)
        i_closest = minloc(abs(temps_available - temp_desired), dim=1)
        temp_actual = temps_available(i_closest)
        if (abs(temp_actual - temp_desired) < tolerance) then
          if (find(temps_to_read, nint(temp_actual)) == -1) then
            call temps_to_read % push_back(nint(temp_actual))
          end if
        else
          call fatal_error("Nuclear data library does not contain cross sections &
               &for " // trim(this % name) // " at or near " // &
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
    allocate(this % data(n_temperature))

    do t = 1, n_temperature
      ! Get temperature as a string
      temp_str = trim(to_str(temps_to_read % data(t))) // "K"

      ! Read exact temperature value
      call read_dataset(this % kTs(t), kT_group, temp_str)

      ! Open group for temperature i
      T_group = open_group(group_id, temp_str)

      ! Coherent elastic data
      if (object_exists(T_group, 'elastic')) then
        ! Read cross section data
        elastic_group = open_group(T_group, 'elastic')
        dset_id = open_dataset(elastic_group, 'xs')
        call read_attribute(type, dset_id, 'type')
        call get_shape(dset_id, dims2)
        allocate(temp(dims2(1), dims2(2)))
        call read_dataset(temp, dset_id)
        call close_dataset(dset_id)

        ! Set cross section data and type
        this % data(t) % n_elastic_e_in = int(dims2(1), 4)
        allocate(this % data(t) % elastic_e_in(this % data(t) % n_elastic_e_in))
        allocate(this % data(t) % elastic_P(this % data(t) % n_elastic_e_in))
        this % data(t) % elastic_e_in(:) = temp(:, 1)
        this % data(t) % elastic_P(:) = temp(:, 2)
        select case (type)
        case ('tab1')
          this % data(t) % elastic_mode = SAB_ELASTIC_DISCRETE
        case ('bragg')
          this % data(t) % elastic_mode = SAB_ELASTIC_EXACT
        end select
        deallocate(temp)

        ! Set elastic threshold
        this % data(t) % threshold_elastic = this % data(t) % elastic_e_in(&
             this % data(t) % n_elastic_e_in)

        ! Read angle distribution
        if (this % data(t) % elastic_mode /= SAB_ELASTIC_EXACT) then
          dset_id = open_dataset(elastic_group, 'mu_out')
          call get_shape(dset_id, dims2)
          this % data(t) % n_elastic_mu = int(dims2(1), 4)
          allocate(this % data(t) % elastic_mu(dims2(1), dims2(2)))
          call read_dataset(this % data(t) % elastic_mu, dset_id)
          call close_dataset(dset_id)
        end if

        call close_group(elastic_group)
      end if

      ! Inelastic data
      if (object_exists(T_group, 'inelastic')) then
        ! Read type of inelastic data
        inelastic_group = open_group(T_group, 'inelastic')

        ! Read cross section data
        dset_id = open_dataset(inelastic_group, 'xs')
        call get_shape(dset_id, dims2)
        allocate(temp(dims2(1), dims2(2)))
        call read_dataset(temp, dset_id)
        call close_dataset(dset_id)

        ! Set cross section data
        this % data(t) % n_inelastic_e_in = int(dims2(1), 4)
        allocate(this % data(t) % inelastic_e_in(this % data(t) % n_inelastic_e_in))
        allocate(this % data(t) % inelastic_sigma(this % data(t) % n_inelastic_e_in))
        this % data(t) % inelastic_e_in(:) = temp(:, 1)
        this % data(t) % inelastic_sigma(:) = temp(:, 2)
        deallocate(temp)

        ! Set inelastic threshold
        this % data(t) % threshold_inelastic = this % data(t) % inelastic_e_in(&
             this % data(t) % n_inelastic_e_in)

        if (this % secondary_mode /= SAB_SECONDARY_CONT) then
          ! Read energy distribution
          dset_id = open_dataset(inelastic_group, 'energy_out')
          call get_shape(dset_id, dims2)
          this % data(t) % n_inelastic_e_out = int(dims2(1), 4)
          allocate(this % data(t) % inelastic_e_out(dims2(1), dims2(2)))
          call read_dataset(this % data(t) % inelastic_e_out, dset_id)
          call close_dataset(dset_id)

          ! Read angle distribution
          dset_id = open_dataset(inelastic_group, 'mu_out')
          call get_shape(dset_id, dims3)
          this % data(t) % n_inelastic_mu = int(dims3(1), 4)
          allocate(this % data(t) % inelastic_mu(dims3(1), dims3(2), dims3(3)))
          call read_dataset(this % data(t) % inelastic_mu, dset_id)
          call close_dataset(dset_id)
        else
          ! Read correlated angle-energy distribution
          call correlated_dist % from_hdf5(inelastic_group)

          ! Convert to S(a,b) native format
          n_energy = size(correlated_dist % energy)
          allocate(this % data(t) % inelastic_data(n_energy))
          do i = 1, n_energy
            associate (edist => correlated_dist % distribution(i))
              ! Get number of outgoing energies for incoming energy i
              n_energy_out = size(edist % e_out)
              this % data(t) % inelastic_data(i) % n_e_out = n_energy_out
              allocate(this % data(t) % inelastic_data(i) % e_out(n_energy_out))
              allocate(this % data(t) % inelastic_data(i) % e_out_pdf(n_energy_out))
              allocate(this % data(t) % inelastic_data(i) % e_out_cdf(n_energy_out))

              ! Copy outgoing energy distribution
              this % data(t) % inelastic_data(i) % e_out(:) = edist % e_out
              this % data(t) % inelastic_data(i) % e_out_pdf(:) = edist % p
              this % data(t) % inelastic_data(i) % e_out_cdf(:) = edist % c

              do j = 1, n_energy_out
                select type (adist => edist % angle(j) % obj)
                type is (Tabular)
                  ! On first pass, allocate space for angles
                  if (j == 1) then
                    n_mu = size(adist % x)
                    this % data(t) % n_inelastic_mu = n_mu
                    allocate(this % data(t) % inelastic_data(i) % mu(&
                         n_mu, n_energy_out))
                  end if

                  ! Copy outgoing angles
                  this % data(t) % inelastic_data(i) % mu(:, j) = adist % x
                end select
              end do
            end associate
          end do

          ! Clear data on correlated angle-energy object
          deallocate(correlated_dist % breakpoints)
          deallocate(correlated_dist % interpolation)
          deallocate(correlated_dist % energy)
          deallocate(correlated_dist % distribution)
        end if

        call close_group(inelastic_group)
      end if
      call close_group(T_group)
    end do

    call close_group(kT_group)
  end subroutine salphabeta_from_hdf5

!===============================================================================
! SAB_CALCULATE_XS determines the elastic and inelastic scattering
! cross-sections in the thermal energy range.
!===============================================================================

  subroutine sab_calculate_xs(this, E, sqrtkT, i_temp, elastic, inelastic)
    class(SAlphaBeta), intent(in) :: this ! S(a,b) object
    real(8), intent(in) :: E          ! energy
    real(8), intent(in) :: sqrtkT     ! temperature
    integer, intent(out) :: i_temp    ! index in the S(a,b)'s temperature
    real(8), intent(out) :: elastic   ! thermal elastic cross section
    real(8), intent(out) :: inelastic ! thermal inelastic cross section

    integer :: i_grid    ! index on S(a,b) energy grid
    real(8) :: f         ! interp factor on S(a,b) energy grid
    real(8) :: kT

    ! Determine temperature for S(a,b) table
    kT = sqrtkT**2
    if (temperature_method == TEMPERATURE_NEAREST) then
      ! If using nearest temperature, do linear search on temperature
      do i_temp = 1, size(this % kTs)
        if (abs(this % kTs(i_temp) - kT) < &
             K_BOLTZMANN*temperature_tolerance) exit
      end do
    else
      ! Find temperatures that bound the actual temperature
      do i_temp = 1, size(this % kTs) - 1
        if (this % kTs(i_temp) <= kT .and. &
             kT < this % kTs(i_temp + 1)) exit
      end do

      ! Randomly sample between temperature i and i+1
      f = (kT - this % kTs(i_temp)) / &
           (this % kTs(i_temp + 1) - this % kTs(i_temp))
      if (f > prn()) i_temp = i_temp + 1
    end if


    ! Get pointer to S(a,b) table
    associate (sab => this % data(i_temp))

      ! Get index and interpolation factor for inelastic grid
      if (E < sab % inelastic_e_in(1)) then
        i_grid = 1
        f = ZERO
      else
        i_grid = binary_search(sab % inelastic_e_in, sab % n_inelastic_e_in, E)
        f = (E - sab%inelastic_e_in(i_grid)) / &
             (sab%inelastic_e_in(i_grid+1) - sab%inelastic_e_in(i_grid))
      end if

      ! Calculate S(a,b) inelastic scattering cross section
      inelastic = (ONE - f) * sab % inelastic_sigma(i_grid) + &
           f * sab % inelastic_sigma(i_grid + 1)

      ! Check for elastic data
      if (E < sab % threshold_elastic) then
        ! Determine whether elastic scattering is given in the coherent or
        ! incoherent approximation. For coherent, the cross section is
        ! represented as P/E whereas for incoherent, it is simply P

        if (sab % elastic_mode == SAB_ELASTIC_EXACT) then
          if (E < sab % elastic_e_in(1)) then
            ! If energy is below that of the lowest Bragg peak, the elastic
            ! cross section will be zero
            elastic = ZERO
          else
            i_grid = binary_search(sab % elastic_e_in, &
                 sab % n_elastic_e_in, E)
            elastic = sab % elastic_P(i_grid) / E
          end if
        else
          ! Determine index on elastic energy grid
          if (E < sab % elastic_e_in(1)) then
            i_grid = 1
          else
            i_grid = binary_search(sab % elastic_e_in, &
                 sab % n_elastic_e_in, E)
          end if

          ! Get interpolation factor for elastic grid
          f = (E - sab%elastic_e_in(i_grid))/(sab%elastic_e_in(i_grid+1) - &
               sab%elastic_e_in(i_grid))

          ! Calculate S(a,b) elastic scattering cross section
          elastic = (ONE - f) * sab % elastic_P(i_grid) + &
               f * sab % elastic_P(i_grid + 1)
        end if
      else
        ! No elastic data
        elastic = ZERO
      end if
    end associate

  end subroutine sab_calculate_xs


!===============================================================================
! FREE_MEMORY_SAB deallocates global arrays defined in this module
!===============================================================================

  subroutine free_memory_sab()
    n_sab_tables = 0
    if (allocated(sab_tables)) deallocate(sab_tables)
    call sab_dict % clear()
  end subroutine free_memory_sab

end module sab_header
