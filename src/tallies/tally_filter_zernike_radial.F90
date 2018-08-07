module tally_filter_zernike_radial

  use, intrinsic :: ISO_C_BINDING

  use constants
  use error
  use hdf5_interface
  use math,                only: calc_zn_rad
  use particle_header,     only: Particle
  use string,              only: to_str
  use tally_filter_header
  use xml_interface

  implicit none
  private

!===============================================================================
! ZERNIKERADIALFILTER gives even order radial Zernike polynomial moments of a 
! particle's position
!===============================================================================

  type, public, extends(TallyFilter) :: ZernikeRadialFilter
    integer :: order
    real(8) :: x
    real(8) :: y
    real(8) :: r
  contains
    procedure :: from_xml
    procedure :: get_all_bins
    procedure :: to_statepoint
    procedure :: text_label
  end type ZernikeRadialFilter

contains

!===============================================================================
! ZernikeRadialFilter methods
!===============================================================================

  subroutine from_xml(this, node)
    class(ZernikeRadialFilter), intent(inout) :: this
    type(XMLNode), intent(in) :: node

    integer :: n

    ! Get center of cylinder and radius
    call get_node_value(node, "x", this % x)
    call get_node_value(node, "y", this % y)
    call get_node_value(node, "r", this % r)

    ! Get specified order
    call get_node_value(node, "order", n)
    this % order = n
    this % n_bins = n/2 + 1 
  end subroutine from_xml

  subroutine get_all_bins(this, p, estimator, match)
    class(ZernikeRadialFilter), intent(in)  :: this
    type(Particle),      intent(in)  :: p
    integer,             intent(in)  :: estimator
    type(TallyFilterMatch),   intent(inout) :: match

    integer :: i
    real(8) :: x, y, r
    real(C_DOUBLE) :: zn_rad(this % n_bins)

    ! Determine normalized (r,theta) positions
    x = p % coord(1) % xyz(1) - this % x
    y = p % coord(1) % xyz(2) - this % y
    r = sqrt(x*x + y*y)/this % r
    if (r <= 1) then

      ! Get moments for even order Zernike polynomial orders 0..n
      call calc_zn_rad(this % order, r, zn_rad)

      do i = 1, this % n_bins
        call match % bins % push_back(i)
        call match % weights % push_back(zn_rad(i))
      end do
    endif
  end subroutine get_all_bins

  subroutine to_statepoint(this, filter_group)
    class(ZernikeRadialFilter), intent(in) :: this
    integer(HID_T),      intent(in) :: filter_group

    call write_dataset(filter_group, "type", "zernikeradial")
    call write_dataset(filter_group, "n_bins", this % n_bins)
    call write_dataset(filter_group, "order", this % order)
    call write_dataset(filter_group, "x", this % x)
    call write_dataset(filter_group, "y", this % y)
    call write_dataset(filter_group, "r", this % r)
  end subroutine to_statepoint

  function text_label(this, bin) result(label)
    class(ZernikeRadialFilter), intent(in) :: this
    integer,             intent(in) :: bin
    character(MAX_LINE_LEN)         :: label

    label = "Zernike expansion, Z" // trim(to_str(2*bin)) // ",0"
  end function text_label

!===============================================================================
!                               C API FUNCTIONS
!===============================================================================

  function openmc_zernike_radial_filter_get_order(index, order) result(err) bind(C)
    ! Get the order of an expansion filter
    integer(C_INT32_T), value       :: index
    integer(C_INT),     intent(out) :: order
    integer(C_INT) :: err

    err = verify_filter(index)
    if (err == 0) then
      select type (f => filters(index) % obj)
      type is (ZernikeRadialFilter)
        order = f % order
      class default
        err = E_INVALID_TYPE
        call set_errmsg("Not a Zernike Radial filter.")
      end select
    end if
  end function openmc_zernike_radial_filter_get_order


  function openmc_zernike_radial_filter_get_params(index, x, y, r) result(err) bind(C)
    ! Get the Zernike filter parameters
    integer(C_INT32_T), value :: index
    real(C_DOUBLE), intent(out) :: x
    real(C_DOUBLE), intent(out) :: y
    real(C_DOUBLE), intent(out) :: r
    integer(C_INT) :: err

    err = verify_filter(index)
    if (err == 0) then
      select type (f => filters(index) % obj)
      type is (ZernikeRadialFilter)
        x = f % x
        y = f % y
        r = f % r
      class default
        err = E_INVALID_TYPE
        call set_errmsg("Not a Zernike Radial filter.")
      end select
    end if
  end function openmc_zernike_radial_filter_get_params


  function openmc_zernike_radial_filter_set_order(index, order) result(err) bind(C)
    ! Set the order of an expansion filter
    integer(C_INT32_T), value :: index
    integer(C_INT),     value :: order
    integer(C_INT) :: err

    err = verify_filter(index)
    if (err == 0) then
      select type (f => filters(index) % obj)
      type is (ZernikeRadialFilter)
        f % order = order
        f % n_bins = order/2 + 1
      class default
        err = E_INVALID_TYPE
        call set_errmsg("Not a Zernike Radial filter.")
      end select
    end if
  end function openmc_zernike_radial_filter_set_order


  function openmc_zernike_radial_filter_set_params(index, x, y, r) result(err) bind(C)
    ! Set the Zernike filter parameters
    integer(C_INT32_T), value :: index
    real(C_DOUBLE), intent(in), optional :: x
    real(C_DOUBLE), intent(in), optional :: y
    real(C_DOUBLE), intent(in), optional :: r
    integer(C_INT) :: err

    err = verify_filter(index)
    if (err == 0) then
      select type (f => filters(index) % obj)
      type is (ZernikeRadialFilter)
        if (present(x)) f % x = x
        if (present(y)) f % y = y
        if (present(r)) f % r = r
      class default
        err = E_INVALID_TYPE
        call set_errmsg("Not a Zernike Radial filter.")
      end select
    end if
  end function openmc_zernike_radial_filter_set_params

end module tally_filter_zernike_radial
