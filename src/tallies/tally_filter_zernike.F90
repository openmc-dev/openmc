module tally_filter_zernike

  use, intrinsic :: ISO_C_BINDING

  use constants
  use error
  use hdf5_interface
  use math,                only: calc_zn, calc_zn_rad
  use particle_header,     only: Particle
  use string,              only: to_str
  use tally_filter_header
  use xml_interface

  implicit none
  private
  public :: openmc_zernike_filter_get_order
  public :: openmc_zernike_filter_get_params
  public :: openmc_zernike_filter_set_order
  public :: openmc_zernike_filter_set_params

!===============================================================================
! ZERNIKEFILTER gives Zernike polynomial moments of a particle's position
!===============================================================================

  type, public, extends(TallyFilter) :: ZernikeFilter
    integer :: order
    real(8) :: x
    real(8) :: y
    real(8) :: r
  contains
    procedure :: calc_n_bins => calc_n_bins_zn
    procedure :: from_xml => from_xml_zn
    procedure :: get_all_bins => get_all_bins_zn
    procedure :: to_statepoint => to_statepoint_zn
    procedure :: text_label => text_label_zn
  end type ZernikeFilter

!===============================================================================
! ZERNIKERADIALFILTER gives even order radial Zernike polynomial moments of a
! particle's position
!===============================================================================

  type, public, extends(ZernikeFilter) :: ZernikeRadialFilter
  contains
    procedure :: calc_n_bins => calc_n_bins_zn_rad
    ! Inherit from_xml from ZernikeFilter
    procedure :: get_all_bins => get_all_bins_zn_rad
    procedure :: to_statepoint => to_statepoint_zn_rad
    procedure :: text_label => text_label_zn_rad
  end type ZernikeRadialFilter

contains

!===============================================================================
! ZernikeFilter methods
!===============================================================================

  function calc_n_bins_zn(this) result(n_bins)
    class(ZernikeFilter), intent(in) :: this
    integer :: n
    integer :: n_bins

    n = this % order
    n_bins = ((n+1) * (n+2))/2
  end function calc_n_bins_zn

  subroutine from_xml_zn(this, node)
    class(ZernikeFilter), intent(inout) :: this
    type(XMLNode), intent(in) :: node

    integer :: n

    ! Get center of cylinder and radius
    call get_node_value(node, "x", this % x)
    call get_node_value(node, "y", this % y)
    call get_node_value(node, "r", this % r)

    ! Get specified order
    call get_node_value(node, "order", n)
    this % order = n
    this % n_bins = this % calc_n_bins()
  end subroutine from_xml_zn

  subroutine get_all_bins_zn(this, p, estimator, match)
    class(ZernikeFilter), intent(in)  :: this
    type(Particle),      intent(in)  :: p
    integer,             intent(in)  :: estimator
    type(TallyFilterMatch),   intent(inout) :: match

    integer :: i
    real(8) :: x, y, r, theta
    real(C_DOUBLE) :: zn(this % n_bins)

    ! Determine normalized (r,theta) positions
    x = p % coord(1) % xyz(1) - this % x
    y = p % coord(1) % xyz(2) - this % y
    r = sqrt(x*x + y*y)/this % r
    if (r <= ONE) then
      theta = atan2(y, x)
      
      ! Get moments for Zernike polynomial orders 0..n
      call calc_zn(this % order, r, theta, zn)

      do i = 1, this % n_bins
        call match % bins % push_back(i)
        call match % weights % push_back(zn(i))
      end do
    endif
  end subroutine get_all_bins_zn

  subroutine to_statepoint_zn(this, filter_group)
    class(ZernikeFilter), intent(in) :: this
    integer(HID_T),      intent(in) :: filter_group

    call write_dataset(filter_group, "type", "zernike")
    call write_dataset(filter_group, "n_bins", this % n_bins)
    call write_dataset(filter_group, "order", this % order)
    call write_dataset(filter_group, "x", this % x)
    call write_dataset(filter_group, "y", this % y)
    call write_dataset(filter_group, "r", this % r)
  end subroutine to_statepoint_zn

  function text_label_zn(this, bin) result(label)
    class(ZernikeFilter), intent(in) :: this
    integer,             intent(in) :: bin
    character(MAX_LINE_LEN)         :: label

    integer :: n, m
    integer :: first, last

    do n = 0, this % order
      last = (n + 1)*(n + 2)/2
      if (bin <= last) then
        first = last - n
        m = -n + (bin - first)*2
        label = "Zernike expansion, Z" // trim(to_str(n)) // "," &
             // trim(to_str(m))
        exit
      end if
    end do
  end function text_label_zn

!===============================================================================
! ZernikeRadialFilter methods
!===============================================================================

  function calc_n_bins_zn_rad(this) result(num_n_bins)
    class(ZernikeRadialFilter), intent(in) :: this
    integer :: n
    integer :: num_n_bins

    n = this % order
    num_n_bins = n/2 + 1
  end function calc_n_bins_zn_rad

  subroutine get_all_bins_zn_rad(this, p, estimator, match)
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
    if (r <= ONE) then

      ! Get moments for even order Zernike polynomial orders 0..n
      call calc_zn_rad(this % order, r, zn_rad)

      do i = 1, this % n_bins
        call match % bins % push_back(i)
        call match % weights % push_back(zn_rad(i))
      end do
    endif
  end subroutine get_all_bins_zn_rad

  subroutine to_statepoint_zn_rad(this, filter_group)
    class(ZernikeRadialFilter), intent(in) :: this
    integer(HID_T),      intent(in) :: filter_group

    call write_dataset(filter_group, "type", "zernikeradial")
    call write_dataset(filter_group, "n_bins", this % n_bins)
    call write_dataset(filter_group, "order", this % order)
    call write_dataset(filter_group, "x", this % x)
    call write_dataset(filter_group, "y", this % y)
    call write_dataset(filter_group, "r", this % r)
  end subroutine to_statepoint_zn_rad

  function text_label_zn_rad(this, bin) result(label)
    class(ZernikeRadialFilter), intent(in) :: this
    integer,             intent(in) :: bin
    character(MAX_LINE_LEN)         :: label

    label = "Zernike expansion, Z" // trim(to_str(2*(bin-1))) // ",0"
  end function text_label_zn_rad

!===============================================================================
!                               C API FUNCTIONS
!===============================================================================

  function openmc_zernike_filter_get_order(index, order) result(err) bind(C)
    ! Get the order of an expansion filter
    integer(C_INT32_T), value       :: index
    integer(C_INT),     intent(out) :: order
    integer(C_INT) :: err

    err = verify_filter(index)
    if (err == 0) then
      select type (f => filters(index) % obj)
      type is (ZernikeFilter)
        order = f % order
      type is (ZernikeRadialFilter)
        order = f % order
      class default
        err = E_INVALID_TYPE
        call set_errmsg("Not a Zernike filter.")
      end select
    end if
  end function openmc_zernike_filter_get_order


  function openmc_zernike_filter_get_params(index, x, y, r) result(err) bind(C)
    ! Get the Zernike filter parameters
    integer(C_INT32_T), value :: index
    real(C_DOUBLE), intent(out) :: x
    real(C_DOUBLE), intent(out) :: y
    real(C_DOUBLE), intent(out) :: r
    integer(C_INT) :: err

    err = verify_filter(index)
    if (err == 0) then
      select type (f => filters(index) % obj)
      class is (ZernikeFilter)
        x = f % x
        y = f % y
        r = f % r
      class default
        err = E_INVALID_TYPE
        call set_errmsg("Not a Zernike filter.")
      end select
    end if
  end function openmc_zernike_filter_get_params


  function openmc_zernike_filter_set_order(index, order) result(err) bind(C)
    ! Set the order of an expansion filter
    integer(C_INT32_T), value :: index
    integer(C_INT),     value :: order
    integer(C_INT) :: err

    err = verify_filter(index)
    if (err == 0) then
      select type (f => filters(index) % obj)
      type is (ZernikeFilter)
        f % order = order
        f % n_bins = ((order + 1)*(order + 2))/2
      type is (ZernikeRadialFilter)
        f % order = order
        f % n_bins = order/2 + 1
      class default
        err = E_INVALID_TYPE
        call set_errmsg("Not a Zernike filter.")
      end select
    end if
  end function openmc_zernike_filter_set_order


  function openmc_zernike_filter_set_params(index, x, y, r) result(err) bind(C)
    ! Set the Zernike filter parameters
    integer(C_INT32_T), value :: index
    real(C_DOUBLE), intent(in), optional :: x
    real(C_DOUBLE), intent(in), optional :: y
    real(C_DOUBLE), intent(in), optional :: r
    integer(C_INT) :: err

    err = verify_filter(index)
    if (err == 0) then
      select type (f => filters(index) % obj)
      class is (ZernikeFilter)
        if (present(x)) f % x = x
        if (present(y)) f % y = y
        if (present(r)) f % r = r
      class default
        err = E_INVALID_TYPE
        call set_errmsg("Not a Zernike filter.")
      end select
    end if
  end function openmc_zernike_filter_set_params

end module tally_filter_zernike
