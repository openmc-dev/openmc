module tally_filter_energy

  use, intrinsic :: ISO_C_BINDING

  use hdf5, only: HID_T

  use algorithm,           only: binary_search
  use constants
  use error
  use hdf5_interface
  use mgxs_header,         only: num_energy_groups, rev_energy_bins
  use particle_header,     only: Particle
  use settings,            only: run_CE
  use string,              only: to_str
  use tally_filter_header
  use xml_interface

  implicit none
  private
  public :: openmc_energy_filter_get_bins
  public :: openmc_energy_filter_set_bins

!===============================================================================
! ENERGYFILTER bins the incident neutron energy.
!===============================================================================

  type, public, extends(TallyFilter) :: EnergyFilter
    real(8), allocatable :: bins(:)

    ! True if transport group number can be used directly to get bin number
    logical :: matches_transport_groups = .false.
  contains
    procedure :: from_xml => from_xml_energy
    procedure :: get_all_bins => get_all_bins_energy
    procedure :: to_statepoint => to_statepoint_energy
    procedure :: text_label => text_label_energy
  end type EnergyFilter

!===============================================================================
! ENERGYOUTFILTER bins the outgoing neutron energy.  Only scattering events use
! the get_all_bins functionality.  Nu-fission tallies manually iterate over the
! filter bins.
!===============================================================================

  type, public, extends(EnergyFilter) :: EnergyoutFilter
  contains
    ! Inherit from_xml from EnergyFilter
    procedure :: get_all_bins => get_all_bins_energyout
    procedure :: to_statepoint => to_statepoint_energyout
    procedure :: text_label => text_label_energyout
  end type EnergyoutFilter

contains

!===============================================================================
! EnergyFilter methods
!===============================================================================

  subroutine from_xml_energy(this, node)
    class(EnergyFilter), intent(inout) :: this
    type(XMLNode), intent(in) :: node

    integer :: n

    n = node_word_count(node, "bins")

    ! Allocate and store bins
    this % n_bins = n - 1
    allocate(this % bins(n))
    call get_node_array(node, "bins", this % bins)

    ! We can save tallying time if we know that the tally bins match
    ! the energy group structure.  In that case, the matching bin
    ! index is simply the group (after flipping for the different
    ! ordering of the library and tallying systems).
    if (.not. run_CE) then
      if (n == num_energy_groups + 1) then
        if (all(this % bins == rev_energy_bins)) &
             then
          this % matches_transport_groups = .true.
        end if
      end if
    end if
  end subroutine from_xml_energy

  subroutine get_all_bins_energy(this, p, estimator, match)
    class(EnergyFilter), intent(in)  :: this
    type(Particle),      intent(in)  :: p
    integer,             intent(in)  :: estimator
    type(TallyFilterMatch),   intent(inout) :: match

    integer :: n
    integer :: bin
    real(8) :: E

    n = this % n_bins

    if (p % g /= NONE .and. this % matches_transport_groups) then
      if (estimator == ESTIMATOR_TRACKLENGTH) then
        call match % bins % push_back(num_energy_groups - p % g + 1)
        call match % weights % push_back(ONE)
      else
        call match % bins % push_back(num_energy_groups - p % last_g + 1)
        call match % weights % push_back(ONE)
      end if

    else
      ! Pre-collision energy of particle
      E = p % last_E

      ! Search to find incoming energy bin.
      bin = binary_search(this % bins, n + 1, E)
      if (bin /= NO_BIN_FOUND) then
        call match % bins % push_back(bin)
        call match % weights % push_back(ONE)
      end if
    end if
  end subroutine get_all_bins_energy

  subroutine to_statepoint_energy(this, filter_group)
    class(EnergyFilter), intent(in) :: this
    integer(HID_T),      intent(in) :: filter_group

    call write_dataset(filter_group, "type", "energy")
    call write_dataset(filter_group, "n_bins", this % n_bins)
    call write_dataset(filter_group, "bins", this % bins)
  end subroutine to_statepoint_energy

  function text_label_energy(this, bin) result(label)
    class(EnergyFilter), intent(in) :: this
    integer,             intent(in) :: bin
    character(MAX_LINE_LEN)         :: label

    real(8) :: E0, E1

    E0 = this % bins(bin)
    E1 = this % bins(bin + 1)
    label = "Incoming Energy [" // trim(to_str(E0)) // ", " &
         // trim(to_str(E1)) // ")"
  end function text_label_energy

!===============================================================================
! EnergyoutFilter methods
!===============================================================================

  subroutine get_all_bins_energyout(this, p, estimator, match)
    class(EnergyoutFilter), intent(in)  :: this
    type(Particle),         intent(in)  :: p
    integer,                intent(in)  :: estimator
    type(TallyFilterMatch),      intent(inout) :: match

    integer :: n
    integer :: bin

    n = this % n_bins

    if (p % g /= NONE .and. this % matches_transport_groups) then
      ! Tallies are ordered in increasing groups, group indices
      ! however are the opposite, so switch
      call match % bins % push_back(num_energy_groups - p % g + 1)
      call match % weights % push_back(ONE)

    else

      ! Search to find incoming energy bin.
      bin = binary_search(this % bins, n + 1, p % E)
      if (bin /= NO_BIN_FOUND) then
        call match % bins % push_back(bin)
        call match % weights % push_back(ONE)
      end if
    end if
  end subroutine get_all_bins_energyout

  subroutine to_statepoint_energyout(this, filter_group)
    class(EnergyoutFilter), intent(in) :: this
    integer(HID_T),      intent(in) :: filter_group

    call write_dataset(filter_group, "type", "energyout")
    call write_dataset(filter_group, "n_bins", this % n_bins)
    call write_dataset(filter_group, "bins", this % bins)
  end subroutine to_statepoint_energyout

  function text_label_energyout(this, bin) result(label)
    class(EnergyoutFilter), intent(in) :: this
    integer,             intent(in) :: bin
    character(MAX_LINE_LEN)         :: label

    real(8) :: E0, E1

    E0 = this % bins(bin)
    E1 = this % bins(bin + 1)
    label = "Outgoing Energy [" // trim(to_str(E0)) // ", " &
         // trim(to_str(E1)) // ")"
  end function text_label_energyout

!===============================================================================
!                               C API FUNCTIONS
!===============================================================================

  function openmc_energy_filter_get_bins(index, energies, n) result(err) bind(C)
    ! Return the bounding energies for an energy filter
    integer(C_INT32_T), value :: index
    type(C_PTR), intent(out) :: energies
    integer(C_INT32_T), intent(out) :: n
    integer(C_INT) :: err

    if (index >= 1 .and. index <= n_filters) then
      if (allocated(filters(index) % obj)) then
        select type (f => filters(index) % obj)
        type is (EnergyFilter)
          energies = C_LOC(f % bins)
          n = size(f % bins)
          err = 0
        class default
          err = E_INVALID_TYPE
          call set_errmsg("Tried to get energy bins on a non-energy filter.")
        end select
      else
        err = E_ALLOCATE
        call set_errmsg("Filter type has not been set yet.")
      end if
    else
      err = E_OUT_OF_BOUNDS
      call set_errmsg("Index in filters array out of bounds.")
    end if
  end function openmc_energy_filter_get_bins


  function openmc_energy_filter_set_bins(index, n, energies) result(err) bind(C)
    ! Set the bounding energies for an energy filter
    integer(C_INT32_T), value, intent(in) :: index
    integer(C_INT32_T), value, intent(in) :: n
    real(C_DOUBLE), intent(in) :: energies(n)
    integer(C_INT) :: err

    err = 0
    if (index >= 1 .and. index <= n_filters) then
      if (allocated(filters(index) % obj)) then
        select type (f => filters(index) % obj)
        type is (EnergyFilter)
          f % n_bins = n - 1
          if (allocated(f % bins)) deallocate(f % bins)
          allocate(f % bins(n))
          f % bins(:) = energies
        class default
          err = E_INVALID_TYPE
          call set_errmsg("Tried to get energy bins on a non-energy filter.")
        end select
      else
        err = E_ALLOCATE
        call set_errmsg("Filter type has not been set yet.")
      end if
    else
      err = E_OUT_OF_BOUNDS
      call set_errmsg("Index in filters array out of bounds.")
    end if
  end function openmc_energy_filter_set_bins

end module tally_filter_energy
