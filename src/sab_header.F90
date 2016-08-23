module sab_header

  use, intrinsic :: ISO_FORTRAN_ENV

  use constants
  use distribution_univariate, only: Tabular
  use error,       only: warning
  use hdf5, only: HID_T, HSIZE_T, SIZE_T
  use h5lt, only: h5ltpath_valid_f, h5iget_name_f
  use hdf5_interface, only: read_attribute, get_shape, open_group, close_group, &
       open_dataset, read_dataset, close_dataset, get_datasets
  use secondary_correlated, only: CorrelatedAngleEnergy
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

  type SAlphaBeta
    character(100) :: name     ! name of table, e.g. lwtr.10t
    real(8)        :: awr      ! weight of nucleus in neutron masses
    real(8)        :: kT       ! temperature in MeV (k*T)
    integer        :: n_zaid   ! Number of valid zaids
    integer, allocatable :: zaid(:) ! List of valid Z and A identifiers, e.g. 6012

    ! threshold for S(a,b) treatment (usually ~4 eV)
    real(8) :: threshold_inelastic
    real(8) :: threshold_elastic = ZERO

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
  contains
    procedure :: print => salphabeta_print
    procedure :: from_hdf5 => salphabeta_from_hdf5
  end type SAlphaBeta

contains

!===============================================================================
! PRINT_SAB_TABLE displays information about a S(a,b) table containing data
! describing thermal scattering from bound materials such as hydrogen in water.
!===============================================================================

  subroutine salphabeta_print(this, unit)
    class(SAlphaBeta), intent(in)  :: this
    integer, intent(in), optional :: unit

    integer :: size_sab   ! memory used by S(a,b) table
    integer :: unit_      ! unit to write to
    integer :: i          ! Loop counter for parsing through this % zaid
    integer :: char_count ! Counter for the number of characters on a line

    ! set default unit for writing information
    if (present(unit)) then
      unit_ = unit
    else
      unit_ = OUTPUT_UNIT
    end if

    ! Basic S(a,b) table information
    write(unit_,*) 'S(a,b) Table ' // trim(this % name)
    write(unit_,'(A)',advance="no") '   zaids = '
    ! Initialize the counter based on the above string
    char_count = 11
    do i = 1, this % n_zaid
      ! Deal with a line thats too long
      if (char_count >= 73) then  ! 73 = 80 - (5 ZAID chars + 1 space + 1 comma)
        ! End the line
        write(unit_,*) ""
        ! Add 11 leading blanks
        write(unit_,'(A)', advance="no") "           "
        ! reset the counter to 11
        char_count = 11
      end if
      if (i < this % n_zaid) then
        ! Include a comma
        write(unit_,'(A)',advance="no") trim(to_str(this % zaid(i))) // ", "
        char_count = char_count + len(trim(to_str(this % zaid(i)))) + 2
      else
        ! Don't include a comma, since we are all done
        write(unit_,'(A)',advance="no") trim(to_str(this % zaid(i)))
      end if

    end do
    write(unit_,*) "" ! Move to next line
    write(unit_,*) '  awr = ' // trim(to_str(this % awr))
    write(unit_,*) '  kT = ' // trim(to_str(this % kT))

    ! Inelastic data
    write(unit_,*) '  # of Incoming Energies (Inelastic) = ' // &
         trim(to_str(this % n_inelastic_e_in))
    write(unit_,*) '  # of Outgoing Energies (Inelastic) = ' // &
         trim(to_str(this % n_inelastic_e_out))
    write(unit_,*) '  # of Outgoing Angles (Inelastic) = ' // &
         trim(to_str(this % n_inelastic_mu))
    write(unit_,*) '  Threshold for Inelastic = ' // &
         trim(to_str(this % threshold_inelastic))

    ! Elastic data
    if (this % n_elastic_e_in > 0) then
      write(unit_,*) '  # of Incoming Energies (Elastic) = ' // &
           trim(to_str(this % n_elastic_e_in))
      write(unit_,*) '  # of Outgoing Angles (Elastic) = ' // &
           trim(to_str(this % n_elastic_mu))
      write(unit_,*) '  Threshold for Elastic = ' // &
           trim(to_str(this % threshold_elastic))
    end if

    ! Determine memory used by S(a,b) table and write out
    size_sab = 8 * (this % n_inelastic_e_in * (2 + this % n_inelastic_e_out * &
         (1 + this % n_inelastic_mu)) + this % n_elastic_e_in * &
         (2 + this % n_elastic_mu))
    write(unit_,*) '  Memory Used = ' // trim(to_str(size_sab)) // ' bytes'

    ! Blank line at end
    write(unit_,*)

  end subroutine salphabeta_print

  subroutine salphabeta_from_hdf5(this, group_id, temperature)
    class(SAlphaBeta), intent(inout) :: this
    integer(HID_T),    intent(in)    :: group_id
    character(6),      intent(in)    :: temperature

    integer :: i, j
    integer :: n_energy, n_energy_out, n_mu
    integer :: hdf5_err
    integer(SIZE_T) :: name_len, name_file_len
    integer(HID_T) :: T_group
    integer(HID_T) :: elastic_group
    integer(HID_T) :: inelastic_group
    integer(HID_T) :: dset_id
    integer(HID_T) :: kT_group, kT_dset
    integer(HSIZE_T) :: dims2(2)
    integer(HSIZE_T) :: dims3(3)
    real(8), allocatable :: temp(:,:)
    character(20) :: type
    logical :: exists
    type(CorrelatedAngleEnergy) :: correlated_dist
    character(MAX_FILE_LEN), allocatable :: temperatures(:)
    character(MAX_FILE_LEN) :: temp_str
    integer, allocatable :: temperatures_integer(:)
    integer :: temperature_delta
    character(6) :: my_temperature
    integer :: temperature_integer

    ! Get name of table from group
    name_len = len(this % name)
    call h5iget_name_f(group_id, this % name, name_len, name_file_len, hdf5_err)

    ! Get rid of leading '/'
    this % name = trim(this % name(2:))

    call read_attribute(this % awr, group_id, 'atomic_weight_ratio')
    call read_attribute(this % zaid, group_id, 'zaids')
    call read_attribute(type, group_id, 'secondary_mode')
    select case (type)
    case ('equal')
      this % secondary_mode = SAB_SECONDARY_EQUAL
    case ('skewed')
      this % secondary_mode = SAB_SECONDARY_SKEWED
    case ('continuous')
      this % secondary_mode = SAB_SECONDARY_CONT
    end select
    this % n_zaid = size(this % zaid)
    kT_group = open_group(group_id, 'kTs')
    ! Before accessing the temperature data, see if the user-provied temperature
    ! exists.  We can find this out by looking at the datasets within kT_group
    temp_str = adjustr(trim(temperature))
    temperature_integer = str_to_int(temp_str(1:len(temp_str) - 1))
    call get_datasets(kT_group, temperatures)
    allocate(temperatures_integer(size(temperatures)))
    do i = 1, size(temperatures)
      temp_str = adjustr(trim(temperatures(i)))
      temperatures_integer(i) = str_to_int(temp_str(1:len(temp_str) - 1))
    end do
    j = 1
    temperature_delta = temperature_integer - temperatures_integer(j)
    do i = 2, size(temperatures)
      if (abs(temperature_integer - temperatures_integer(i)) < temperature_delta) &
           j = i
    end do
    ! Now print a warning if there is no matching temperature and then use the
    ! closest temperature
    my_temperature = temperatures(j)
    if (temperature /= my_temperature) then
      call warning(trim(this % name) // " does not contain data at a &
                   &temperature of " // trim(temperature) // "; using the &
                   &nearest available temperature of " // trim(my_temperature))
    end if

    kT_dset = open_dataset(kT_group, my_temperature)
    call read_dataset(this % kT, kT_dset)
    call close_dataset(kT_dset)
    call close_group(kT_group)

    ! Open my_temperature group
    T_group = open_group(group_id, my_temperature)

    ! Coherent elastic data
    call h5ltpath_valid_f(T_group, 'elastic', .true., exists, hdf5_err)
    if (exists) then
      ! Read cross section data
      elastic_group = open_group(T_group, 'elastic')
      dset_id = open_dataset(elastic_group, 'xs')
      call read_attribute(type, dset_id, 'type')
      call get_shape(dset_id, dims2)
      allocate(temp(dims2(1), dims2(2)))
      call read_dataset(temp, dset_id)
      call close_dataset(dset_id)

      ! Set cross section data and type
      this % n_elastic_e_in = int(dims2(1), 4)
      allocate(this % elastic_e_in(this % n_elastic_e_in))
      allocate(this % elastic_P(this % n_elastic_e_in))
      this % elastic_e_in(:) = temp(:, 1)
      this % elastic_P(:) = temp(:, 2)
      select case (type)
      case ('tab1')
        this % elastic_mode = SAB_ELASTIC_DISCRETE
      case ('bragg')
        this % elastic_mode = SAB_ELASTIC_EXACT
      end select
      deallocate(temp)

      ! Set elastic threshold
      this % threshold_elastic = this % elastic_e_in(this % n_elastic_e_in)

      ! Read angle distribution
      if (this % elastic_mode /= SAB_ELASTIC_EXACT) then
        dset_id = open_dataset(elastic_group, 'mu_out')
        call get_shape(dset_id, dims2)
        this % n_elastic_mu = int(dims2(1), 4)
        allocate(this % elastic_mu(dims2(1), dims2(2)))
        call read_dataset(this % elastic_mu, dset_id)
        call close_dataset(dset_id)
      end if

      call close_group(elastic_group)
    end if

    ! Inelastic data
    call h5ltpath_valid_f(T_group, 'inelastic', .true., exists, hdf5_err)
    if (exists) then
      ! Read type of inelastic data
      inelastic_group = open_group(T_group, 'inelastic')

      ! Read cross section data
      dset_id = open_dataset(inelastic_group, 'xs')
      call get_shape(dset_id, dims2)
      allocate(temp(dims2(1), dims2(2)))
      call read_dataset(temp, dset_id)
      call close_dataset(dset_id)

      ! Set cross section data
      this % n_inelastic_e_in = int(dims2(1), 4)
      allocate(this % inelastic_e_in(this % n_inelastic_e_in))
      allocate(this % inelastic_sigma(this % n_inelastic_e_in))
      this % inelastic_e_in(:) = temp(:, 1)
      this % inelastic_sigma(:) = temp(:, 2)
      deallocate(temp)

      ! Set inelastic threshold
      this % threshold_inelastic = this % inelastic_e_in(this % n_inelastic_e_in)

      if (this % secondary_mode /= SAB_SECONDARY_CONT) then
        ! Read energy distribution
        dset_id = open_dataset(inelastic_group, 'energy_out')
        call get_shape(dset_id, dims2)
        this % n_inelastic_e_out = int(dims2(1), 4)
        allocate(this % inelastic_e_out(dims2(1), dims2(2)))
        call read_dataset(this % inelastic_e_out, dset_id)
        call close_dataset(dset_id)

        ! Read angle distribution
        dset_id = open_dataset(inelastic_group, 'mu_out')
        call get_shape(dset_id, dims3)
        this % n_inelastic_mu = int(dims3(1), 4)
        allocate(this % inelastic_mu(dims3(1), dims3(2), dims3(3)))
        call read_dataset(this % inelastic_mu, dset_id)
        call close_dataset(dset_id)
      else
        ! Read correlated angle-energy distribution
        call correlated_dist % from_hdf5(inelastic_group)

        ! Convert to S(a,b) native format
        n_energy = size(correlated_dist % energy)
        allocate(this % inelastic_data(n_energy))
        do i = 1, n_energy
          associate (edist => correlated_dist % distribution(i))
            ! Get number of outgoing energies for incoming energy i
            n_energy_out = size(edist % e_out)
            this % inelastic_data(i) % n_e_out = n_energy_out
            allocate(this % inelastic_data(i) % e_out(n_energy_out))
            allocate(this % inelastic_data(i) % e_out_pdf(n_energy_out))
            allocate(this % inelastic_data(i) % e_out_cdf(n_energy_out))

            ! Copy outgoing energy distribution
            this % inelastic_data(i) % e_out(:) = edist % e_out
            this % inelastic_data(i) % e_out_pdf(:) = edist % p
            this % inelastic_data(i) % e_out_cdf(:) = edist % c

            do j = 1, n_energy_out
              select type (adist => edist % angle(j) % obj)
              type is (Tabular)
                ! On first pass, allocate space for angles
                if (j == 1) then
                  n_mu = size(adist % x)
                  this % n_inelastic_mu = n_mu
                  allocate(this % inelastic_data(i) % mu(n_mu, n_energy_out))
                end if

                ! Copy outgoing angles
                this % inelastic_data(i) % mu(:, j) = adist % x
              end select
            end do
          end associate
        end do
      end if

      call close_group(inelastic_group)
    end if
    call close_group(T_group)
  end subroutine salphabeta_from_hdf5

end module sab_header
