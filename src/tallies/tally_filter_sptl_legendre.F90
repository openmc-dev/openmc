module tally_filter_sptl_legendre

  use, intrinsic :: ISO_C_BINDING

  use hdf5, only: HID_T

  use constants
  use error
  use hdf5_interface
  use math,                only: calc_pn
  use particle_header,     only: Particle
  use string,              only: to_str, to_lower
  use tally_filter_header
  use xml_interface

  implicit none
  private

  integer, parameter :: AXIS_X = 1
  integer, parameter :: AXIS_Y = 2
  integer, parameter :: AXIS_Z = 3

!===============================================================================
! SpatialLEGENDREFILTER gives Legendre moments
!===============================================================================

  type, public, extends(TallyFilter) :: SpatialLegendreFilter
    integer :: order
    integer :: axis
    real(8) :: min
    real(8) :: max
  contains
    procedure :: from_xml
    procedure :: get_all_bins
    procedure :: to_statepoint
    procedure :: text_label
  end type SpatialLegendreFilter

contains

!===============================================================================
! SpatialLegendreFilter methods
!===============================================================================

  subroutine from_xml(this, node)
    class(SpatialLegendreFilter), intent(inout) :: this
    type(XMLNode), intent(in) :: node

    character(MAX_WORD_LEN) :: axis

    ! Get attributes from XML
    call get_node_value(node, "order", this % order)
    call get_node_value(node, "axis", axis)
    select case (to_lower(axis))
    case ('x')
      this % axis = AXIS_X
    case ('y')
      this % axis = AXIS_Y
    case ('z')
      this % axis = AXIS_Z
    end select
    call get_node_value(node, "min", this % min)
    call get_node_value(node, "max", this % max)

    this % n_bins = this % order + 1
  end subroutine from_xml

  subroutine get_all_bins(this, p, estimator, match)
    class(SpatialLegendreFilter), intent(in)  :: this
    type(Particle),      intent(in)  :: p
    integer,             intent(in)  :: estimator
    type(TallyFilterMatch),   intent(inout) :: match

    integer :: i
    real(8) :: wgt
    real(8) :: x       ! Position on specified axis
    real(8) :: x_norm  ! Normalized position

    x = p % coord(1) % xyz(this % axis)
    if (this % min <= x .and. x <= this % max) then
      ! Calculate normalized position between min and max
      x_norm = TWO*(x - this % min)/(this % max - this % min) - ONE

      ! TODO: Use recursive formula to calculate higher orders
      do i = 0, this % order
        wgt = calc_pn(i, x_norm)
        call match % bins % push_back(i + 1)
        call match % weights % push_back(wgt)
      end do
    end if
  end subroutine get_all_bins

  subroutine to_statepoint(this, filter_group)
    class(SpatialLegendreFilter), intent(in) :: this
    integer(HID_T),      intent(in) :: filter_group

    call write_dataset(filter_group, "type", "spatiallegendre")
    call write_dataset(filter_group, "n_bins", this % n_bins)
    call write_dataset(filter_group, "order", this % order)
    select case (this % axis)
    case (AXIS_X)
      call write_dataset(filter_group, 'axis', 'x')
    case (AXIS_Y)
      call write_dataset(filter_group, 'axis', 'y')
    case (AXIS_Z)
      call write_dataset(filter_group, 'axis', 'z')
    end select
    call write_dataset(filter_group, 'min', this % min)
    call write_dataset(filter_group, 'max', this % max)
  end subroutine to_statepoint

  function text_label(this, bin) result(label)
    class(SpatialLegendreFilter), intent(in) :: this
    integer,             intent(in) :: bin
    character(MAX_LINE_LEN)         :: label

    label = "Legendre expansion, P" // trim(to_str(bin - 1))
  end function text_label

!===============================================================================
!                               C API FUNCTIONS
!===============================================================================


end module tally_filter_sptl_legendre
