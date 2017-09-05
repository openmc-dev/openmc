module tally_filter_legendre

  use, intrinsic :: ISO_C_BINDING

  use hdf5, only: HID_T

  use constants
  use error
  use hdf5_interface
  use math,                only: calc_pn
  use particle_header,     only: Particle
  use string,              only: to_str
  use tally_filter_header
  use xml_interface

  implicit none
  private

!===============================================================================
! LEGENDREFILTER gives Legendre moments of the change in scattering angle
!===============================================================================

  type, public, extends(TallyFilter) :: LegendreFilter
    integer :: order
  contains
    procedure :: from_xml
    procedure :: get_all_bins
    procedure :: to_statepoint
    procedure :: text_label
  end type LegendreFilter

contains

!===============================================================================
! LegendreFilter methods
!===============================================================================

  subroutine from_xml(this, node)
    class(LegendreFilter), intent(inout) :: this
    type(XMLNode), intent(in) :: node

    ! Get specified order
    call get_node_value(node, "order", this % order)
    this % n_bins = this % order + 1
  end subroutine from_xml

  subroutine get_all_bins(this, p, estimator, match)
    class(LegendreFilter), intent(in)  :: this
    type(Particle),      intent(in)  :: p
    integer,             intent(in)  :: estimator
    type(TallyFilterMatch),   intent(inout) :: match

    integer :: i
    real(8) :: wgt

    ! TODO: Use recursive formula to calculate higher orders
    do i = 0, this % order
      wgt = calc_pn(i, p % mu)
      call match % bins % push_back(i + 1)
      call match % weights % push_back(wgt)
    end do
  end subroutine get_all_bins

  subroutine to_statepoint(this, filter_group)
    class(LegendreFilter), intent(in) :: this
    integer(HID_T),      intent(in) :: filter_group

    call write_dataset(filter_group, "type", "legendre")
    call write_dataset(filter_group, "n_bins", this % n_bins)
    call write_dataset(filter_group, "order", this % order)
  end subroutine to_statepoint

  function text_label(this, bin) result(label)
    class(LegendreFilter), intent(in) :: this
    integer,             intent(in) :: bin
    character(MAX_LINE_LEN)         :: label

    label = "Legendre expansion, P" // trim(to_str(bin - 1))
  end function text_label

!===============================================================================
!                               C API FUNCTIONS
!===============================================================================


end module tally_filter_legendre
