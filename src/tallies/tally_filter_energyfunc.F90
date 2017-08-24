module tally_filter_energyfunc

  use, intrinsic :: ISO_C_BINDING

  use hdf5

  use algorithm,          only: binary_search
  use constants
  use error,              only: fatal_error
  use hdf5_interface
  use particle_header,    only: Particle
  use settings,           only: run_CE
  use string,             only: to_str
  use tally_filter_header
  use xml_interface

  implicit none
  private

!===============================================================================
! EnergyFunctionFilter multiplies tally scores by an arbitrary function of
! incident energy described by a piecewise linear-linear interpolation.
!===============================================================================

  type, public, extends(TallyFilter) :: EnergyFunctionFilter
    real(8), allocatable :: energy(:)
    real(8), allocatable :: y(:)

  contains
    procedure :: from_xml
    procedure :: get_all_bins => get_all_bins_energyfunction
    procedure :: to_statepoint => to_statepoint_energyfunction
    procedure :: text_label => text_label_energyfunction
  end type EnergyFunctionFilter

contains

  subroutine from_xml(this, node)
    class(EnergyFunctionFilter), intent(inout) :: this
    type(XMLNode), intent(in) :: node

    integer :: n

    this % n_bins = 1
    ! Make sure this is continuous-energy mode.
    if (.not. run_CE) then
      call fatal_error("EnergyFunction filters are only supported for &
           &continuous-energy transport calculations")
    end if

    ! Allocate and store energy grid.
    if (.not. check_for_node(node, "energy")) then
      call fatal_error("Energy grid not specified for EnergyFunction &
           &filter.")
    end if
    n = node_word_count(node, "energy")
    allocate(this % energy(n))
    call get_node_array(node, "energy", this % energy)

    ! Allocate and store interpolant values.
    if (.not. check_for_node(node, "y")) then
      call fatal_error("y values not specified for EnergyFunction &
           &filter.")
    end if
    n = node_word_count(node, "y")
    allocate(this % y(n))
    call get_node_array(node, "y", this % y)
  end subroutine from_xml

  subroutine get_all_bins_energyfunction(this, p, estimator, match)
    class(EnergyFunctionFilter), intent(in)  :: this
    type(Particle),              intent(in)  :: p
    integer,                     intent(in)  :: estimator
    type(TallyFilterMatch),      intent(inout) :: match

    integer :: n, indx
    real(8) :: E, f, weight

    select type(this)
    type is (EnergyFunctionFilter)
      n = size(this % energy)

      ! Get pre-collision energy of particle
      E = p % last_E

      ! Search to find incoming energy bin.
      indx = binary_search(this % energy, n, E)

      ! Compute an interpolation factor between nearest bins.
      f = (E - this % energy(indx)) &
           / (this % energy(indx+1) - this % energy(indx))

      ! Interpolate on the lin-lin grid.
      call match % bins % push_back(1)
      weight = (ONE - f) * this % y(indx) + f * this % y(indx+1)
      call match % weights % push_back(weight)
    end select
  end subroutine get_all_bins_energyfunction

  subroutine to_statepoint_energyfunction(this, filter_group)
    class(EnergyFunctionFilter), intent(in) :: this
    integer(HID_T),              intent(in) :: filter_group

    select type(this)
    type is (EnergyFunctionFilter)
      call write_dataset(filter_group, "type", "energyfunction")
      call write_dataset(filter_group, "energy", this % energy)
      call write_dataset(filter_group, "y", this % y)
    end select
  end subroutine to_statepoint_energyfunction

  function text_label_energyfunction(this, bin) result(label)
    class(EnergyFunctionFilter), intent(in) :: this
    integer,                     intent(in) :: bin
    character(MAX_LINE_LEN)                 :: label

    select type(this)
    type is (EnergyFunctionFilter)
      write(label, FMT="(A, ES8.1, A, ES8.1, A, ES8.1, A, ES8.1, A)") &
           "Energy Function f([", this % energy(1), ", ..., ", &
           this % energy(size(this % energy)), "]) = [", this % y(1), &
           ", ..., ", this % y(size(this % y)), "]"
    end select
  end function text_label_energyfunction

end module tally_filter_energyfunc
