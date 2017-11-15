module photon_header

  use hdf5, only: HID_T, HSIZE_T, SIZE_T

  use algorithm,   only: binary_search
  use constants,   only: ZERO, HALF, SUBSHELLS
  use dict_header, only: DictIntInt
  use endf_header, only: Tabulated1D
  use hdf5_interface

  real(8), allocatable :: compton_profile_pz(:)
  real(8), allocatable :: ttb_energy_electron(:) ! incident electron energy grid
  real(8), allocatable :: ttb_energy_photon(:)   ! reduced photon energy grid
  real(8)              :: ttb_energy_cutoff

  type ElectronSubshell
    integer :: index_subshell  ! index in SUBSHELLS
    integer :: threshold
    real(8) :: n_electrons
    real(8) :: binding_energy
    real(8), allocatable :: cross_section(:)

    ! Transition data
    integer :: n_transitions
    integer, allocatable :: transition_subshells(:,:)
    real(8), allocatable :: transition_energy(:)
    real(8), allocatable :: transition_probability(:)
  end type ElectronSubshell

  type PhotonInteraction
    character(3) :: name  ! atomic symbol, e.g. 'Zr'
    integer      :: Z     ! atomic number

    ! Microscopic cross sections
    real(8), allocatable :: energy(:)
    real(8), allocatable :: coherent(:)
    real(8), allocatable :: incoherent(:)
    real(8), allocatable :: photoelectric_total(:)
    real(8), allocatable :: pair_production_total(:)
    real(8), allocatable :: pair_production_electron(:)
    real(8), allocatable :: pair_production_nuclear(:)

    ! Form factors
    type(Tabulated1D) :: incoherent_form_factor
    type(Tabulated1D) :: coherent_int_form_factor
    type(Tabulated1D) :: coherent_anomalous_real
    type(Tabulated1D) :: coherent_anomalous_imag

    ! Photoionization and atomic relaxation data
    type(DictIntInt) :: shell_dict ! Given a shell designator, e.g. 3, this
                                   ! dictionary gives an index in shells(:)
    type(ElectronSubshell), allocatable :: shells(:)

    ! Compton profile data
    real(8), allocatable :: profile_pdf(:,:)
    real(8), allocatable :: profile_cdf(:,:)
    real(8), allocatable :: binding_energy(:)
    real(8), allocatable :: electron_pdf(:)

    ! Stopping power data
    real(8) :: density
    real(8), allocatable :: stopping_power_collision(:)
    real(8), allocatable :: stopping_power_radiative(:)

    ! Bremsstrahlung scaled DCS
    real(8), allocatable :: dcs(:,:)

  contains
    procedure :: from_hdf5 => photon_from_hdf5
  end type PhotonInteraction

  type Bremsstrahlung
    integer :: i_material ! Index in materials array

    real(8), allocatable :: yield(:,:) ! Photon number yield
    real(8), allocatable :: dcs(:,:)   ! Scaled bremsstrahlung DCS
    real(8), allocatable :: cdf(:,:)   ! Bremsstrahlung energy CDF

  contains
    procedure :: init => bremsstrahlung_init
  end type Bremsstrahlung

!===============================================================================
! ELEMENTMICROXS contains cached microscopic photon cross sections for a
! particular element at the current energy
!===============================================================================

  type ElementMicroXS
    integer :: index_grid      ! index on element energy grid
    real(8) :: last_E = ZERO   ! last evaluated energy
    real(8) :: interp_factor   ! interpolation factor on energy grid
    real(8) :: total           ! microscropic total photon xs
    real(8) :: coherent        ! microscopic coherent xs
    real(8) :: incoherent      ! microscopic incoherent xs
    real(8) :: photoelectric   ! microscopic photoelectric xs
    real(8) :: pair_production ! microscopic pair production xs
  end type ElementMicroXS

contains

  subroutine photon_from_hdf5(this, group_id)
    class(PhotonInteraction), intent(inout) :: this
    integer(HID_T), intent(in) :: group_id

    integer          :: i, j
    integer(HID_T)   :: rgroup, tgroup
    integer(HID_T)   :: dset_id
    integer(HSIZE_T) :: dims(1), dims2(2)
    integer(SIZE_T)  :: name_len
    integer          :: n_energy
    integer          :: n_shell
    integer          :: n_profile
    integer          :: n_transition
    integer          :: n_k
    integer          :: n_e
    character(3), allocatable :: designators(:)
    real(8)          :: c
    real(8), allocatable :: matrix(:,:)

    ! Get name of nuclide from group
    name_len = len(this % name)
    this % name = get_name(group_id, name_len)

    ! Get rid of leading '/'
    this % name = trim(this % name(2:))

    ! Get atomic number
    call read_attribute(this % Z, group_id, 'Z')

    ! Determine number of energies and read energy grid
    dset_id = open_dataset(group_id, 'energy')
    call get_shape(dset_id, dims)
    n_energy = int(dims(1), 4)
    allocate(this % energy(dims(1)))
    call read_dataset(this % energy, dset_id)
    call close_dataset(dset_id)

    ! Allocate arrays
    allocate(this % coherent(n_energy))
    allocate(this % incoherent(n_energy))
    allocate(this % pair_production_total(n_energy))
    allocate(this % pair_production_nuclear(n_energy))
    allocate(this % pair_production_electron(n_energy))
    allocate(this % photoelectric_total(n_energy))

    ! Read coherent scattering
    rgroup = open_group(group_id, 'coherent')
    call read_dataset(this % coherent, rgroup, 'xs')

    dset_id = open_dataset(rgroup, 'integrated_scattering_factor')
    call this % coherent_int_form_factor % from_hdf5(dset_id)
    call close_dataset(dset_id)

    dset_id = open_dataset(rgroup, 'anomalous_real')
    call this % coherent_anomalous_real % from_hdf5(dset_id)
    call close_dataset(dset_id)

    dset_id = open_dataset(rgroup, 'anomalous_imag')
    call this % coherent_anomalous_imag % from_hdf5(dset_id)
    call close_dataset(dset_id)
    call close_group(rgroup)

    ! Read incoherent scattering
    rgroup = open_group(group_id, 'incoherent')
    call read_dataset(this % incoherent, rgroup, 'xs')
    dset_id = open_dataset(rgroup, 'scattering_factor')
    call this % incoherent_form_factor % from_hdf5(dset_id)
    call close_dataset(dset_id)
    call close_group(rgroup)

    ! Read pair production
    rgroup = open_group(group_id, 'pair_production')
    call read_dataset(this % pair_production_nuclear, rgroup, 'xs')
    call close_group(rgroup)

    ! Read pair production
    rgroup = open_group(group_id, 'triplet_production')
    call read_dataset(this % pair_production_electron, rgroup, 'xs')
    call close_group(rgroup)

    ! Read photoelectric
    rgroup = open_group(group_id, 'photoelectric')
    call read_dataset(this % photoelectric_total, rgroup, 'xs')
    call close_group(rgroup)

    ! Read subshell photoionization cross section and atomic relaxation data
    rgroup = open_group(group_id, 'subshells')
    call read_attribute(designators, rgroup, 'designators')
    n_shell = size(designators)
    allocate(this % shells(n_shell))
    do i = 1, n_shell
      ! Create mapping from designator to index
      do j = 1, size(SUBSHELLS)
        if (designators(i) == SUBSHELLS(j)) then
          call this % shell_dict % add_key(j, i)
          this % shells(i) % index_subshell = j
          exit
        end if
      end do

      ! Read binding energy and number of electrons
      tgroup = open_group(rgroup, trim(designators(i)))
      call read_attribute(this % shells(i) % binding_energy, tgroup, &
           'binding_energy')
      call read_attribute(this % shells(i) % n_electrons, tgroup, &
           'num_electrons')

      ! Read subshell cross section
      dset_id = open_dataset(tgroup, 'xs')
      call read_attribute(j, dset_id, 'threshold_idx')
      this % shells(i) % threshold = j
      allocate(this % shells(i) % cross_section(n_energy - j))
      call read_dataset(this % shells(i) % cross_section, dset_id)
      call close_dataset(dset_id)
      where (this % shells(i) % cross_section > ZERO)
        this % shells(i) % cross_section = log(this % shells(i) % cross_section)
      elsewhere
        this % shells(i) % cross_section = -500.0_8
      end where

      if (object_exists(tgroup, 'transitions')) then
        dset_id = open_dataset(tgroup, 'transitions')
        call get_shape(dset_id, dims2)
        n_transition = int(dims2(2), 4)
        this % shells(i) % n_transitions = n_transition
        if (n_transition > 0) then
          allocate(this % shells(i) % transition_subshells(2, n_transition))
          allocate(this % shells(i) % transition_energy(n_transition))
          allocate(this % shells(i) % transition_probability(n_transition))

          allocate(matrix(dims2(1), dims2(2)))
          call read_dataset(matrix, dset_id)

          this % shells(i) % transition_subshells(:,:) = int(matrix(1:2, :), 4)
          this % shells(i) % transition_energy(:) = matrix(3, :)
          this % shells(i) % transition_probability(:) = matrix(4, :)
          deallocate(matrix)
        end if
        call close_dataset(dset_id)
      else
        this % shells(i) % n_transitions = 0
      end if
      call close_group(tgroup)
    end do
    call close_group(rgroup)
    deallocate(designators)

    ! Determine number of electron shells
    rgroup = open_group(group_id, 'compton_profiles')

    ! Determine number of shells
    dset_id = open_dataset(rgroup, 'num_electrons')
    call get_shape(dset_id, dims)
    n_shell = int(dims(1), 4)

    ! Read electron shell PDF and binding energies
    allocate(this % electron_pdf(n_shell), this % binding_energy(n_shell))
    call read_dataset(this % electron_pdf, dset_id)
    call close_dataset(dset_id)
    call read_dataset(this % binding_energy, rgroup, 'binding_energy')
    this % electron_pdf(:) = this % electron_pdf / sum(this % electron_pdf)

    ! Read Compton profiles
    dset_id = open_dataset(rgroup, 'J')
    call get_shape(dset_id, dims2)
    n_profile = int(dims2(1), 4)
    allocate(this % profile_pdf(n_profile, n_shell))
    call read_dataset(this % profile_pdf, dset_id)
    call close_dataset(dset_id)

    ! Get Compton profile momentum grid
    if (.not. allocated(compton_profile_pz)) then
      allocate(compton_profile_pz(n_profile))
      call read_dataset(compton_profile_pz, rgroup, 'pz')
    end if
    call close_group(rgroup)

    ! Create Compton profile CDF
    allocate(this % profile_cdf(n_profile, n_shell))
    do i = 1, n_shell
      c = ZERO
      this % profile_cdf(1,i) = ZERO
      do j = 1, n_profile - 1
        c = c + HALF*(compton_profile_pz(j+1) - compton_profile_pz(j)) * &
             (this%profile_pdf(j,i) + this%profile_pdf(j+1,i))
        this % profile_cdf(j+1,i) = c
      end do
    end do

    ! Calculate total pair production
    this % pair_production_total(:) = this % pair_production_nuclear + &
         this % pair_production_electron

    if (electron_treatment == ELECTRON_TTB) then
      ! TODO: read in
      ttb_energy_cutoff = 1.e-3

      ! Read bremsstrahlung scaled DCS
      rgroup = open_group(group_id, 'bremsstrahlung')
      dset_id = open_dataset(rgroup, 'dcs')
      call get_shape(dset_id, dims2)
      n_k = int(dims2(1), 4)
      n_e = int(dims2(2), 4)
      allocate(this % dcs(n_k, n_e))
      call read_dataset(this % dcs, dset_id)
      call close_dataset(dset_id)

      ! Get energy grids used for bremsstrahlung DCS and for stopping powers
      if (.not. allocated(ttb_energy_electron)) then
        allocate(ttb_energy_electron(n_e))
        call read_dataset(ttb_energy_electron, rgroup, 'electron_energy')
      end if
      if (.not. allocated(ttb_energy_photon)) then
        allocate(ttb_energy_photon(n_k))
        call read_dataset(ttb_energy_photon, rgroup, 'photon_energy')
      end if
      call close_group(rgroup)

      ! Read stopping power data
      if (this % Z < 99) then
        rgroup = open_group(group_id, 'stopping_powers')
        allocate(this % stopping_power_collision(n_e))
        allocate(this % stopping_power_radiative(n_e))
        call read_dataset(this % stopping_power_collision, rgroup, 's_collision')
        call read_dataset(this % stopping_power_radiative, rgroup, 's_radiative')
        call read_attribute(this % density, rgroup, 'density')
        call close_group(rgroup)
      end if
    end if

    ! Take logarithm of energies and cross sections since they are log-log
    ! interpolated
    this % energy = log(this % energy)

    where (this % coherent > ZERO)
      this % coherent = log(this % coherent)
    elsewhere
      this % coherent = -500.0_8
    end where

    where (this % incoherent > ZERO)
      this % incoherent = log(this % incoherent)
    elsewhere
      this % incoherent = -500.0_8
    end where

    where (this % photoelectric_total > ZERO)
      this % photoelectric_total = log(this % photoelectric_total)
    elsewhere
      this % photoelectric_total = -500.0_8
    end where

    where (this % pair_production_total > ZERO)
      this % pair_production_total = log(this % pair_production_total)
    elsewhere
      this % pair_production_total = -500.0_8
    end where

    if (electron_treatment == ELECTRON_TTB) then
      ttb_energy_electron = log(ttb_energy_electron)
    end if

  end subroutine photon_from_hdf5

  subroutine bremsstrahlung_init(this, i_material)
    class(Bremsstrahlung), intent(inout) :: this
    integer, intent(in) :: i_material

    integer                 :: i, j
    integer                 :: i_k
    integer                 :: n_e, n_k
    real(8)                 :: e
    real(8)                 :: c
    real(8)                 :: k, k_l, k_r, k_c
    real(8)                 :: x_l, x_r, x_c
    real(8)                 :: Z_eq_sq
    real(8)                 :: beta
    real(8), allocatable    :: mass_fraction(:)
    real(8), allocatable    :: stopping_power(:)
    real(8), allocatable    :: mfp_inv(:)
    type(Material), pointer :: mat
    type(PhotonInteraction), pointer :: elm

    ! Get pointer to this material
    mat => materials(i_material)
    this % i_material = i_material

    ! Allocate and initialize arrays
    n_k = size(ttb_energy_photon)
    n_e = size(ttb_energy_electron)
    allocate(mass_fraction(size(mat % element)))
    allocate(stopping_power(n_e))
    allocate(mfp_inv(n_e))
    allocate(this % yield(n_e))
    allocate(this % dcs(n_k, n_e))
    allocate(this % cdf(n_k, n_e))
    stopping_power(:) = ZERO
    mfp_inv(:) = ZERO
    this % dcs(:,:) = ZERO
    this % cdf(:,:) = ZERO

    ! TODO
    ! Calculate the "equivalent" atomic number Zeq and get the mass fraction of
    ! each element
    Z_eq_sq = 0
    do i = 1, mat % n_nuclides
      ! Zeq**2 = (atomic fraction of x)*Zx**2 + (atomic fraction of y)*Zy**2 + ...
    end do

    ! Calculate the molecular DCS and the molecular total stopping power using
    ! Bragg's additivity rule. Note: the collision stopping power cannot be
    ! accurately calculated using Bragg's additivity rule since the mean
    ! excitation energies and the density effect corrections cannot simply be
    ! summed together. Bragg's additivity rule fails especially when a
    ! higher-density compound is composed of elements that are in lower-density
    ! form at normal temperature and pressure (at which the NIST stopping
    ! powers are given). It will be used to approximate the collision stopping
    ! powers for now, but should be fixed in the future.
    do i = 1, size(mat % element)
      ! Get pointer to current element
      elm => mat % element(i)

      ! TODO: n_atoms
      ! Accumulate material DCS
      this % dcs = this % dcs + n_atoms(i) * elm % Z**2 / Z_eq_sq * elm % dcs

      ! Accumulate material total stopping power
      stopping_power = stopping_power + mass_fraction(i) * elm % density * &
           (elm % stopping_power_collision + elm % stopping_power_radiative)
    end do

    ! Calculate inverse bremsstrahlung mean free path
    do i = 1, n_e
      e = exp(ttb_energy_electron(i))
      if (e <= ttb_energy_cutoff) cycle

      ! Ratio of the velocity of the charged particle to the speed of light
      beta = sqrt(e*(e + TWO*MASS_ELECTRON/1.e6)) / (e + MASS_ELECTRON/1.e6)

      ! Integration lower bound
      k_c = ttb_energy_cutoff / e

      ! Find the upper bounding index of the reduced photon cutoff energy
      i_k = binary_search(ttb_energy_photon, n_k, k_c) + 1

      ! Get the interpolation bounds
      k_l = ttb_energy_photon(i_k-1)
      k_r = ttb_energy_photon(i_k)
      x_l = this % dcs(i_k-1, i)
      x_r = this % dcs(i_k, i)

      ! Use linear interpolation in reduced photon energy k to find value of
      ! the DCS at the cutoff energy
      x_c = (x_l * (k_r - k_c) + x_r * (k_c - k_l)) / (k_r - k_l)

      ! Integrate using the trapezoidal rule in log-log space
      c = HALF * (log(k_r) - log(k_c)) * (x_c + x_r)
      this % cdf(i_k,i) = c
      do j = i_k, n_k - 1
        c = c + HALF * (log(ttb_energy_photon(j+1)) - &
             log(ttb_energy_photon(j))) * (this % dcs(j,i) + this % dcs(j+1,i))
        this % cdf(j+1,i) = c

      ! TODO: density in atom/cm^3
      mfp_inv(i) = c * mat % density * Z_eq_sq / beta**2 * 1.0e-27_8
    end do

    ! TODO:
    ! Calculate photon number yield
    ! cs = cubic_spline(exp(ttb_energy_electron), mfp_inv / stopping_power)
    ! do i = 1, n_e
    !   yield(i) = cs % integrate(ttb_energy_cutoff, exp(ttb_energy_electron(i)))
    ! Use logarithm of number yield since it is log-log interpolated
    ! yield = log(yield)

  end subroutine bremsstrahlung_init

end module photon_header
