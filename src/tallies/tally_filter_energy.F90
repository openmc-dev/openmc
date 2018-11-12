module tally_filter_energy

  use, intrinsic :: ISO_C_BINDING

  use tally_filter_header
  use xml_interface

  implicit none
  private
  public :: openmc_energy_filter_get_bins
  public :: openmc_energy_filter_set_bins

  interface
    function openmc_energy_filter_get_bins(index, energies, n) result(err) bind(C)
      import C_INT32_T, C_PTR, C_INT
      integer(C_INT32_T), value :: index
      type(C_PTR), intent(out) :: energies
      integer(C_INT32_T), intent(out) :: n
      integer(C_INT) :: err
    end function
    function openmc_energy_filter_set_bins(index, n, energies) result(err) bind(C)
      import C_INT32_T, C_DOUBLE, C_INT
      integer(C_INT32_T), value, intent(in) :: index
      integer(C_INT32_T), value, intent(in) :: n
      real(C_DOUBLE), intent(in) :: energies(n)
      integer(C_INT) :: err
    end function
  end interface

!===============================================================================
! ENERGYFILTER bins the incident neutron energy.
!===============================================================================

  type, public, extends(TallyFilter) :: EnergyFilter
    ! True if transport group number can be used directly to get bin number
    logical :: matches_transport_groups = .false.
  contains
    procedure :: from_xml => from_xml_energy
    procedure :: search
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

    interface
      function energy_filter_matches_transport_groups(filt) result(matches) &
           bind(C)
        import C_PTR, C_BOOL
        type(C_PTR), value :: filt
        logical(C_BOOL)    :: matches
      end function
    end interface

    call this % from_xml_cpp(node)
    this % n_bins = this % n_bins_cpp()
    this % matches_transport_groups = &
         energy_filter_matches_transport_groups(this % ptr)
  end subroutine from_xml_energy

  function search(this, val) result(bin)
    class(EnergyFilter), intent(in) :: this
    real(8),             intent(in) :: val
    integer                         :: bin
    interface
      function energy_filter_search(filt, val) result(bin) bind(C)
        import C_PTR, C_DOUBLE, C_INT
        type(C_PTR), value    :: filt
        real(C_DOUBLE), value :: val
        integer(C_INT)        :: bin
      end function
    end interface
    bin = energy_filter_search(this % ptr, val)
  end function search

end module tally_filter_energy
