module tally_filter_energy

  use, intrinsic :: ISO_C_BINDING

  use algorithm,           only: binary_search
  use constants
  use error
  use hdf5_interface
  use mgxs_interface,      only: num_energy_groups, rev_energy_bins
  use particle_header,     only: Particle
  use settings,            only: run_CE
  use string,              only: to_str
  use tally_filter_cpp
  use xml_interface

  implicit none
  private
  public :: openmc_energy_filter_get_bins
  public :: openmc_energy_filter_set_bins

!===============================================================================
! ENERGYFILTER bins the incident neutron energy.
!===============================================================================

  type, public, extends(CppTallyFilter) :: EnergyFilter
    real(8), allocatable :: bins(:)

    ! True if transport group number can be used directly to get bin number
    logical :: matches_transport_groups = .false.
  contains
    procedure :: from_xml => from_xml_energy
  end type EnergyFilter

!===============================================================================
! ENERGYOUTFILTER bins the outgoing neutron energy.  Only scattering events use
! the get_all_bins functionality.  Nu-fission tallies manually iterate over the
! filter bins.
!===============================================================================

  type, public, extends(EnergyFilter) :: EnergyoutFilter
  end type EnergyoutFilter

contains

!===============================================================================
! EnergyFilter methods
!===============================================================================

  subroutine from_xml_energy(this, node)
    class(EnergyFilter), intent(inout) :: this
    type(XMLNode), intent(in) :: node

    integer :: n

    call this % from_xml_cpp_inner(node)

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

!===============================================================================
!                               C API FUNCTIONS
!===============================================================================

  function openmc_energy_filter_get_bins(index, energies, n) result(err) bind(C)
    ! Return the bounding energies for an energy filter
    integer(C_INT32_T), value :: index
    type(C_PTR), intent(out) :: energies
    integer(C_INT32_T), intent(out) :: n
    integer(C_INT) :: err

    err = verify_filter(index)
    if (err == 0) then
      select type (f => filters(index) % obj)
      type is (EnergyFilter)
        energies = C_LOC(f % bins)
        n = size(f % bins)
        err = 0
      type is (EnergyoutFilter)
        energies = C_LOC(f % bins)
        n = size(f % bins)
        err = 0
        class default
        err = E_INVALID_TYPE
        call set_errmsg("Tried to get energy bins on a non-energy filter.")
      end select
    end if
  end function openmc_energy_filter_get_bins


  function openmc_energy_filter_set_bins(index, n, energies) result(err) bind(C)
    ! Set the bounding energies for an energy filter
    integer(C_INT32_T), value, intent(in) :: index
    integer(C_INT32_T), value, intent(in) :: n
    real(C_DOUBLE), intent(in) :: energies(n)
    integer(C_INT) :: err

    err = verify_filter(index)
    if (err == 0) then
      select type (f => filters(index) % obj)
      type is (EnergyFilter)
        f % n_bins = n - 1
        if (allocated(f % bins)) deallocate(f % bins)
        allocate(f % bins(n))
        f % bins(:) = energies
      type is (EnergyoutFilter)
        f % n_bins = n - 1
        if (allocated(f % bins)) deallocate(f % bins)
        allocate(f % bins(n))
        f % bins(:) = energies
        class default
        err = E_INVALID_TYPE
        call set_errmsg("Tried to get energy bins on a non-energy filter.")
      end select
    end if
  end function openmc_energy_filter_set_bins

end module tally_filter_energy
