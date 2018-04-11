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
  public :: openmc_spatial_legendre_filter_get_order
  public :: openmc_spatial_legendre_filter_get_params
  public :: openmc_spatial_legendre_filter_set_order
  public :: openmc_spatial_legendre_filter_set_params

  integer, parameter :: AXIS_X = 1
  integer, parameter :: AXIS_Y = 2
  integer, parameter :: AXIS_Z = 3

!===============================================================================
! SPATIALLEGENDREFILTER gives Legendre moments of the particle's normalized
! position along an axis
!===============================================================================

  type, public, extends(TallyFilter) :: SpatialLegendreFilter
    integer(C_INT) :: order
    integer(C_INT) :: axis
    real(C_DOUBLE) :: min
    real(C_DOUBLE) :: max
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
    character(kind=C_CHAR) :: axis

    call write_dataset(filter_group, "type", "spatiallegendre")
    call write_dataset(filter_group, "n_bins", this % n_bins)
    call write_dataset(filter_group, "order", this % order)
    select case (this % axis)
    case (AXIS_X)
      axis = 'x'
    case (AXIS_Y)
      axis = 'y'
    case (AXIS_Z)
      axis = 'z'
    end select
    call write_dataset(filter_group, 'axis', axis)
    call write_dataset(filter_group, 'min', this % min)
    call write_dataset(filter_group, 'max', this % max)
  end subroutine to_statepoint

  function text_label(this, bin) result(label)
    class(SpatialLegendreFilter), intent(in) :: this
    integer,             intent(in) :: bin
    character(MAX_LINE_LEN)         :: label

    character(1) :: axis

    select case (this % axis)
    case (AXIS_X)
      axis = 'x'
    case (AXIS_Y)
      axis = 'y'
    case (AXIS_Z)
      axis = 'z'
    end select
    label = "Legendre expansion, " // axis // " axis, P" // &
         trim(to_str(bin - 1))
  end function text_label

!===============================================================================
!                               C API FUNCTIONS
!===============================================================================

  function openmc_spatial_legendre_filter_get_order(index, order) result(err) bind(C)
    ! Get the order of an expansion filter
    integer(C_INT32_T), value       :: index
    integer(C_INT),     intent(out) :: order
    integer(C_INT) :: err

    err = verify_filter(index)
    if (err == 0) then
      select type (f => filters(index) % obj)
      type is (SpatialLegendreFilter)
        order = f % order
      class default
        err = E_INVALID_TYPE
        call set_errmsg("Not a spatial Legendre filter.")
      end select
    end if
  end function openmc_spatial_legendre_filter_get_order


  function openmc_spatial_legendre_filter_set_order(index, order) result(err) bind(C)
    ! Set the order of an expansion filter
    integer(C_INT32_T), value :: index
    integer(C_INT),     value :: order
    integer(C_INT) :: err

    err = verify_filter(index)
    if (err == 0) then
      select type (f => filters(index) % obj)
      type is (SpatialLegendreFilter)
        f % order = order
        f % n_bins = order + 1
      class default
        err = E_INVALID_TYPE
        call set_errmsg("Not a spatial Legendre filter.")
      end select
    end if
  end function openmc_spatial_legendre_filter_set_order


  function openmc_spatial_legendre_filter_get_params(index, axis, min, max) &
       result(err) bind(C)
    ! Get the parameters for a spatial Legendre filter
    integer(C_INT32_T), value :: index
    integer(C_INT), intent(out) :: axis
    real(C_DOUBLE), intent(out) :: min
    real(C_DOUBLE), intent(out) :: max
    integer(C_INT) :: err

    err = verify_filter(index)
    if (err == 0) then
      select type(f => filters(index) % obj)
      type is (SpatialLegendreFilter)
        axis = f % axis
        min = f % min
        max = f % max
      class default
        err = E_INVALID_TYPE
        call set_errmsg("Not a spatial Legendre filter.")
      end select
    end if
  end function openmc_spatial_legendre_filter_get_params


  function openmc_spatial_legendre_filter_set_params(index, axis, min, max) &
       result(err) bind(C)
    ! Set the parameters for a spatial Legendre filter
    integer(C_INT32_T), value :: index
    integer(C_INT), intent(in), optional :: axis
    real(C_DOUBLE), intent(in), optional :: min
    real(C_DOUBLE), intent(in), optional :: max
    integer(C_INT) :: err

    err = verify_filter(index)
    if (err == 0) then
      select type(f => filters(index) % obj)
      type is (SpatialLegendreFilter)
        if (present(axis)) f % axis = axis
        if (present(min)) f % min = min
        if (present(max)) f % max = max
      class default
        err = E_INVALID_TYPE
        call set_errmsg("Not a spatial Legendre filter.")
      end select
    end if
  end function openmc_spatial_legendre_filter_set_params

end module tally_filter_sptl_legendre
