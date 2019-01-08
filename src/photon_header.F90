module photon_header

  use, intrinsic :: ISO_C_BINDING

  use algorithm,        only: binary_search
  use constants
  use dict_header,      only: DictCharInt
  use hdf5_interface
  use settings

  real(8), allocatable :: compton_profile_pz(:)
  real(8), allocatable :: ttb_e_grid(:) ! energy T of incident electron
  real(8), allocatable :: ttb_k_grid(:) ! reduced energy W/T of emitted photon

  type PhotonInteraction
    character(3) :: name  ! atomic symbol, e.g. 'Zr'
    integer      :: Z     ! atomic number

    ! Microscopic cross sections
    real(8), allocatable :: energy(:)

    ! Stopping power data
    real(8) :: I ! mean excitation energy
    real(8), allocatable :: stopping_power_collision(:)
    real(8), allocatable :: stopping_power_radiative(:)

    ! Bremsstrahlung scaled DCS
    real(8), allocatable :: dcs(:,:)

  contains
    procedure :: from_hdf5 => photon_from_hdf5
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

    integer          :: i
    integer(HID_T)   :: rgroup
    integer(HID_T)   :: dset_id
    integer(HSIZE_T) :: dims(1), dims2(2)
    integer          :: n_energy
    integer          :: n_k
    integer          :: n_e
    real(8)          :: f
    real(8)          :: y
    real(8), allocatable :: electron_energy(:)
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

  end subroutine photon_from_hdf5

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
