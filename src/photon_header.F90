module photon_header

  use, intrinsic :: ISO_C_BINDING

  use algorithm,        only: binary_search
  use constants
  use dict_header,      only: DictIntInt, DictCharInt
  use endf_header,      only: Tabulated1D
  use hdf5_interface
  use nuclide_header,   only: nuclides
  use settings

  real(8), allocatable :: compton_profile_pz(:)
  real(8), allocatable :: ttb_e_grid(:) ! energy T of incident electron
  real(8), allocatable :: ttb_k_grid(:) ! reduced energy W/T of emitted photon

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
    real(8) :: I ! mean excitation energy
    real(8), allocatable :: stopping_power_collision(:)
    real(8), allocatable :: stopping_power_radiative(:)

    ! Bremsstrahlung scaled DCS
    real(8), allocatable :: dcs(:,:)

  contains
    procedure :: from_hdf5 => photon_from_hdf5
    procedure :: calculate_xs => photon_calculate_xs
  end type PhotonInteraction

  type BremsstrahlungData
    real(8), allocatable :: pdf(:,:) ! Bremsstrahlung energy PDF
    real(8), allocatable :: cdf(:,:) ! Bremsstrahlung energy CDF
    real(8), allocatable :: yield(:) ! Photon number yield
  end type BremsstrahlungData

  type Bremsstrahlung
    type(BremsstrahlungData) :: electron
    type(BremsstrahlungData) :: positron
  end type Bremsstrahlung

  type(PhotonInteraction), allocatable, target :: elements(:) ! Photon cross sections
  integer :: n_elements       ! Number of photon cross section tables

  type(DictCharInt) :: element_dict

  type(Bremsstrahlung), allocatable, target :: ttb(:) ! Bremsstrahlung data

!===============================================================================
! ELEMENTMICROXS contains cached microscopic photon cross sections for a
! particular element at the current energy
!===============================================================================

  type, bind(C) :: ElementMicroXS
    integer(C_INT) :: index_grid      ! index on element energy grid
    real(C_DOUBLE) :: last_E = ZERO   ! last evaluated energy
    real(C_DOUBLE) :: interp_factor   ! interpolation factor on energy grid
    real(C_DOUBLE) :: total           ! microscropic total photon xs
    real(C_DOUBLE) :: coherent        ! microscopic coherent xs
    real(C_DOUBLE) :: incoherent      ! microscopic incoherent xs
    real(C_DOUBLE) :: photoelectric   ! microscopic photoelectric xs
    real(C_DOUBLE) :: pair_production ! microscopic pair production xs
  end type ElementMicroXS

  type(ElementMicroXS), allocatable, target :: micro_photon_xs(:) ! Cache for each element
!$omp threadprivate(micro_photon_xs)

contains

  subroutine photon_from_hdf5(this, group_id)
    class(PhotonInteraction), intent(inout) :: this
    integer(HID_T), intent(in) :: group_id

    integer          :: i, j
    integer(HID_T)   :: rgroup, tgroup
    integer(HID_T)   :: dset_id
    integer(HSIZE_T) :: dims(1), dims2(2)
    integer          :: n_energy
    integer          :: n_shell
    integer          :: n_profile
    integer          :: n_transition
    integer          :: n_k
    integer          :: n_e
    character(3), allocatable :: designators(:)
    real(8)          :: c
    real(8)          :: f
    real(8)          :: y
    real(8), allocatable :: electron_energy(:)
    real(8), allocatable :: matrix(:,:)
    real(8), allocatable :: dcs(:,:)

    interface
      subroutine photon_from_hdf5_c(group) bind(C)
        import HID_T
        integer(HID_T), value :: group
      end subroutine
    end interface

    ! Read element data on C++ side
    call photon_from_hdf5_c(group_id)

    ! Get name of nuclide from group
    this % name = get_name(group_id)

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

    if (object_exists(group_id, 'anomalous_real')) then
      dset_id = open_dataset(rgroup, 'anomalous_real')
      call this % coherent_anomalous_real % from_hdf5(dset_id)
      call close_dataset(dset_id)
    end if

    if (object_exists(group_id, 'anomalous_imag')) then
      dset_id = open_dataset(rgroup, 'anomalous_imag')
      call this % coherent_anomalous_imag % from_hdf5(dset_id)
      call close_dataset(dset_id)
      call close_group(rgroup)
    end if

    ! Read incoherent scattering
    rgroup = open_group(group_id, 'incoherent')
    call read_dataset(this % incoherent, rgroup, 'xs')
    dset_id = open_dataset(rgroup, 'scattering_factor')
    call this % incoherent_form_factor % from_hdf5(dset_id)
    call close_dataset(dset_id)
    call close_group(rgroup)

    ! Read pair production
    rgroup = open_group(group_id, 'pair_production_electron')
    call read_dataset(this % pair_production_electron, rgroup, 'xs')
    call close_group(rgroup)

    ! Read pair production
    if (object_exists(group_id, 'pair_production_nuclear')) then
      rgroup = open_group(group_id, 'pair_production_nuclear')
      call read_dataset(this % pair_production_nuclear, rgroup, 'xs')
      call close_group(rgroup)
    else
      this % pair_production_nuclear(:) = ZERO
    end if

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
          call this % shell_dict % set(j, i)
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
          this % shells(i) % transition_probability(:) = matrix(4, :) &
               / sum(matrix(4, :))
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
      allocate(electron_energy(n_e))
      call read_dataset(electron_energy, rgroup, 'electron_energy')
      if (.not. allocated(ttb_k_grid)) then
        allocate(ttb_k_grid(n_k))
        call read_dataset(ttb_k_grid, rgroup, 'photon_energy')
      end if
      call close_group(rgroup)

      ! Read stopping power data
      if (this % Z < 99) then
        rgroup = open_group(group_id, 'stopping_powers')
        allocate(this % stopping_power_collision(n_e))
        allocate(this % stopping_power_radiative(n_e))
        call read_dataset(this % stopping_power_collision, rgroup, 's_collision')
        call read_dataset(this % stopping_power_radiative, rgroup, 's_radiative')
        call read_attribute(this % I, rgroup, 'I')
        call close_group(rgroup)
      end if

      ! Truncate the bremsstrahlung data at the cutoff energy
      if (energy_cutoff(PHOTON) > electron_energy(1)) then
        i_grid = binary_search(electron_energy, n_e, energy_cutoff(PHOTON))

        ! calculate interpolation factor
        f = (log(energy_cutoff(PHOTON)) - log(electron_energy(i_grid))) / &
             (log(electron_energy(i_grid+1)) - log(electron_energy(i_grid)))

        ! Interpolate collision stopping power at the cutoff energy and
        ! truncate
        y = exp(log(this % stopping_power_collision(i_grid)) + &
             f*(log(this % stopping_power_collision(i_grid+1)) - &
             log(this % stopping_power_collision(i_grid))))
        this % stopping_power_collision = &
             [y, this % stopping_power_collision(i_grid+1:n_e)]

        ! Interpolate radiative stopping power at the cutoff energy and
        ! truncate
        y = exp(log(this % stopping_power_radiative(i_grid)) + &
             f*(log(this % stopping_power_radiative(i_grid+1)) - &
             log(this % stopping_power_radiative(i_grid))))
        this % stopping_power_radiative = &
             [y, this % stopping_power_radiative(i_grid+1:n_e)]

        ! Interpolate bremsstrahlung DCS at the cutoff energy and truncate
        allocate(dcs(n_k, n_e-i_grid+1))
        do i = 1, n_k
          y = exp(log(this % dcs(i,i_grid)) + &
               f*(log(this % dcs(i,i_grid+1)) - log(this % dcs(i,i_grid))))
          dcs(i,:) = [y, this % dcs(i,i_grid+1:n_e)]
        end do
        call move_alloc(dcs, this % dcs)

        electron_energy = [energy_cutoff(PHOTON), electron_energy(i_grid+1:n_e)]
      end if

      ! Set incident particle energy grid
      if (.not. allocated(ttb_e_grid)) then
        call move_alloc(electron_energy, ttb_e_grid)
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

  end subroutine photon_from_hdf5

!===============================================================================
! CALCULATE_ELEMENT_XS determines microscopic photon cross sections for an
! element of a given index in the elements array at the energy of the given
! particle
!===============================================================================

  subroutine photon_calculate_xs(this, E, xs)
    class(PhotonInteraction), intent(in)    :: this ! index into elements array
    real(8),                  intent(in)    :: E ! energy
    type(ElementMicroXS),     intent(inout) :: xs

    integer :: i_grid  ! index on element energy grid
    integer :: i_shell ! index in subshells
    integer :: i_start ! threshold index
    integer :: n_grid  ! number of grid points
    real(8) :: f       ! interp factor on element energy grid
    real(8) :: log_E   ! logarithm of the energy

    ! Perform binary search on the element energy grid in order to determine
    ! which points to interpolate between
    n_grid = size(this % energy)
    log_E = log(E)
    if (log_E <= this % energy(1)) then
      i_grid = 1
    elseif (log_E > this % energy(n_grid)) then
      i_grid = n_grid - 1
    else
      i_grid = binary_search(this % energy, n_grid, log_E)
    end if

    ! check for case where two energy points are the same
    if (this % energy(i_grid) == this % energy(i_grid+1)) i_grid = i_grid + 1

    ! calculate interpolation factor
    f = (log_E - this % energy(i_grid)) / &
         (this % energy(i_grid+1) - this % energy(i_grid))

    xs % index_grid    = i_grid
    xs % interp_factor = f

    ! Calculate microscopic coherent cross section
    xs % coherent = exp(this % coherent(i_grid) + f * &
         (this % coherent(i_grid+1) - this % coherent(i_grid)))

    ! Calculate microscopic incoherent cross section
    xs % incoherent = exp(this % incoherent(i_grid) + &
         f*(this % incoherent(i_grid+1) - this % incoherent(i_grid)))

    ! Calculate microscopic photoelectric cross section
    xs % photoelectric = ZERO
    do i_shell = 1, size(this % shells)
      ! Check threshold of reaction
      i_start = this % shells(i_shell) % threshold
      if (i_grid <= i_start) cycle

      ! Evaluation subshell photoionization cross section
      xs % photoelectric = xs % photoelectric + &
           exp(this % shells(i_shell) % cross_section(i_grid-i_start) + &
           f*(this % shells(i_shell) % cross_section(i_grid+1-i_start) - &
           this % shells(i_shell) % cross_section(i_grid-i_start)))
    end do

    ! Calculate microscopic pair production cross section
    xs % pair_production = exp(&
         this % pair_production_total(i_grid) + f*(&
         this % pair_production_total(i_grid+1) - &
         this % pair_production_total(i_grid)))

    ! Calculate microscopic total cross section
    xs % total = xs % coherent + xs % incoherent + xs % photoelectric + &
         xs % pair_production

    xs % last_E = E

  end subroutine photon_calculate_xs

!===============================================================================
! FREE_MEMORY_PHOTON deallocates/resets global variables in this module
!===============================================================================

  subroutine free_memory_photon()
    ! Deallocate photon cross section data
    if (allocated(elements)) deallocate(elements)
    if (allocated(compton_profile_pz)) deallocate(compton_profile_pz)
    n_elements = 0
    call element_dict % clear()

    ! Clear TTB-related arrays
    if (allocated(ttb_e_grid)) deallocate(ttb_e_grid)
    if (allocated(ttb)) deallocate(ttb)
  end subroutine free_memory_photon

  function micro_photon_xs_ptr() result(ptr) bind(C)
    type(C_PTR) :: ptr
    if (size(micro_photon_xs) > 0) then
      ptr = C_LOC(micro_photon_xs(1))
    else
      ptr = C_NULL_PTR
    end if
  end function

end module photon_header
