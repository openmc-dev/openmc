module tally_filter_zernike

  use, intrinsic :: ISO_C_BINDING

  use hdf5, only: HID_T

  use constants
  use error
  use hdf5_interface
  use particle_header,     only: Particle
  use string,              only: to_str
  use tally_filter_header
  use xml_interface

  implicit none
  private

!===============================================================================
! ZERNIKEFILTER gives Zernike polynomial moments of a particle's position
!===============================================================================

  type, public, extends(TallyFilter) :: ZernikeFilter
    integer :: order
    real(8) :: x
    real(8) :: y
    real(8) :: r
  contains
    procedure :: from_xml
    procedure :: get_all_bins
    procedure :: to_statepoint
    procedure :: text_label
  end type ZernikeFilter

contains

!===============================================================================
! ZernikeFilter methods
!===============================================================================

  subroutine from_xml(this, node)
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
    this % n_bins = (n + 1)*(n + 2)/2
  end subroutine from_xml

  subroutine get_all_bins(this, p, estimator, match)
    class(ZernikeFilter), intent(in)  :: this
    type(Particle),      intent(in)  :: p
    integer,             intent(in)  :: estimator
    type(TallyFilterMatch),   intent(inout) :: match

    integer :: n
    integer :: j
    integer :: i
    real(8) :: wgt
    real(8) :: tmp(this % order + 1)
    real(8) :: x, y, r, theta

    ! Determine normalized (r,theta) positions
    x = p % coord(1) % xyz(1) - this % x
    y = p % coord(1) % xyz(2) - this % y
    r = sqrt(x*x + y*y)/this % r
    theta = atan2(y, x)

    i = 0
    do n = 0, this % order
      ! Get moments for n-th order Zernike polynomial
      tmp(1:n+1) = calc_zn_scaled(n, r, theta)

      ! Indicate matching bins/weights
      do j = 1, n + 1
        call match % bins % push_back(i + j)
        call match % weights % push_back(tmp(j))
      end do
      i = i + n + 1
    end do
  end subroutine get_all_bins

  subroutine to_statepoint(this, filter_group)
    class(ZernikeFilter), intent(in) :: this
    integer(HID_T),      intent(in) :: filter_group

    call write_dataset(filter_group, "type", "zernike")
    call write_dataset(filter_group, "n_bins", this % n_bins)
    call write_dataset(filter_group, "order", this % order)
    call write_dataset(filter_group, "x", this % x)
    call write_dataset(filter_group, "y", this % y)
    call write_dataset(filter_group, "r", this % r)
  end subroutine to_statepoint

  function text_label(this, bin) result(label)
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
        label = "Zernike expansion, Y" // trim(to_str(n)) // "," &
             // trim(to_str(m))
        exit
      end if
    end do
  end function text_label


!===============================================================================
! CALC_ZN calculates the n-th order Zernike polynomial moment for a given angle
! (rho, theta) location in the unit disk.
!===============================================================================

  pure function calc_zn(n, rho, phi) result(zn)

    integer, intent(in) :: n           ! Order requested
    real(8), intent(in) :: rho         ! Radial location in the unit disk
    real(8), intent(in) :: phi         ! Theta (radians) location in the unit disk
    real(8)             :: zn(n + 1)   ! The resultant Z_n(uvw)

    ! n == radial degree
    ! m == azimuthal frequency

    select case(n)
    case(0)
      ! n = 0, m = 0
      zn(1) = ( ( 1.00 ) ) &
           * ( 1.000000000000 )
    case(1)
      ! n = 1, m = -1
      zn(1) = ( ( 1.00 ) * rho ) &
           * ( 2.000000000000 ) * sin(1.00 * phi)
      ! n = 1, m = 1
      zn(2) = ( ( 1.00 ) * rho ) &
           * ( 2.000000000000 ) * cos(1.00 * phi)
    case(2)
      ! n = 2, m = -2
      zn(1) = ( ( 1.00 ) *rho * rho ) &
           * ( 2.449489742783 ) * sin(2.00 * phi)
      ! n = 2, m = 0
      zn(2) = ( ( -1.00 ) + &
           ( 2.00 ) *rho * rho ) &
           * ( 1.732050807569 )
      ! n = 2, m = 2
      zn(3) = ( ( 1.00 ) *rho * rho ) &
           * ( 2.449489742783 ) * cos(2.00 * phi)
    case(3)
      ! n = 3, m = -3
      zn(1) = ( ( 1.00 ) *rho *rho * rho ) &
           * ( 2.828427124746 ) * sin(3.00 * phi)
      ! n = 3, m = -1
      zn(2) = ( ( -2.00 ) * rho + &
           ( 3.00 ) *rho *rho * rho ) &
           * ( 2.828427124746 ) * sin(3.00 * phi)
      ! n = 3, m = 1
      zn(3) = ( ( -2.00 ) * rho + &
           ( 3.00 ) *rho *rho * rho ) &
           * ( 2.828427124746 ) * cos(3.00 * phi)
      ! n = 3, m = 3
      zn(4) = ( ( 1.00 ) *rho *rho * rho ) &
           * ( 2.828427124746 ) * cos(3.00 * phi)
    case(4)
      ! n = 4, m = -4
      zn(1) = ( ( 1.00 ) *rho *rho *rho * rho ) &
           * ( 3.162277660168 ) * sin(4.00 * phi)
      ! n = 4, m = -2
      zn(2) = ( ( -3.00 ) *rho * rho + &
           ( 4.00 ) *rho *rho *rho * rho ) &
           * ( 3.162277660168 ) * sin(4.00 * phi)
      ! n = 4, m = 0
      zn(3) = ( ( 1.00 ) + &
           ( -6.00 ) *rho * rho + &
           ( 6.00 ) *rho *rho *rho * rho ) &
           * ( 2.236067977500 )
      ! n = 4, m = 2
      zn(4) = ( ( -3.00 ) *rho * rho + &
           ( 4.00 ) *rho *rho *rho * rho ) &
           * ( 3.162277660168 ) * cos(4.00 * phi)
      ! n = 4, m = 4
      zn(5) = ( ( 1.00 ) *rho *rho *rho * rho ) &
           * ( 3.162277660168 ) * cos(4.00 * phi)
    case(5)
      ! n = 5, m = -5
      zn(1) = ( ( 1.00 ) *rho *rho *rho *rho * rho ) &
           * ( 3.464101615138 ) * sin(5.00 * phi)
      ! n = 5, m = -3
      zn(2) = ( ( -4.00 ) *rho *rho * rho + &
           ( 5.00 ) *rho *rho *rho *rho * rho ) &
           * ( 3.464101615138 ) * sin(5.00 * phi)
      ! n = 5, m = -1
      zn(3) = ( ( 3.00 ) * rho + &
           ( -12.00 ) *rho *rho * rho + &
           ( 10.00 ) *rho *rho *rho *rho * rho ) &
           * ( 3.464101615138 ) * sin(5.00 * phi)
      ! n = 5, m = 1
      zn(4) = ( ( 3.00 ) * rho + &
           ( -12.00 ) *rho *rho * rho + &
           ( 10.00 ) *rho *rho *rho *rho * rho ) &
           * ( 3.464101615138 ) * cos(5.00 * phi)
      ! n = 5, m = 3
      zn(5) = ( ( -4.00 ) *rho *rho * rho + &
           ( 5.00 ) *rho *rho *rho *rho * rho ) &
           * ( 3.464101615138 ) * cos(5.00 * phi)
      ! n = 5, m = 5
      zn(6) = ( ( 1.00 ) *rho *rho *rho *rho * rho ) &
           * ( 3.464101615138 ) * cos(5.00 * phi)
    case(6)
      ! n = 6, m = -6
      zn(1) = ( ( 1.00 ) *rho *rho *rho *rho *rho * rho ) &
           * ( 3.741657386774 ) * sin(6.00 * phi)
      ! n = 6, m = -4
      zn(2) = ( ( -5.00 ) *rho *rho *rho * rho + &
           ( 6.00 ) *rho *rho *rho *rho *rho * rho ) &
           * ( 3.741657386774 ) * sin(6.00 * phi)
      ! n = 6, m = -2
      zn(3) = ( ( 6.00 ) *rho * rho + &
           ( -20.00 ) *rho *rho *rho * rho + &
           ( 15.00 ) *rho *rho *rho *rho *rho * rho ) &
           * ( 3.741657386774 ) * sin(6.00 * phi)
      ! n = 6, m = 0
      zn(4) = ( ( -1.00 ) + &
           ( 12.00 ) *rho * rho + &
           ( -30.00 ) *rho *rho *rho * rho + &
           ( 20.00 ) *rho *rho *rho *rho *rho * rho ) &
           * ( 2.645751311065 )
      ! n = 6, m = 2
      zn(5) = ( ( 6.00 ) *rho * rho + &
           ( -20.00 ) *rho *rho *rho * rho + &
           ( 15.00 ) *rho *rho *rho *rho *rho * rho ) &
           * ( 3.741657386774 ) * cos(6.00 * phi)
      ! n = 6, m = 4
      zn(6) = ( ( -5.00 ) *rho *rho *rho * rho + &
           ( 6.00 ) *rho *rho *rho *rho *rho * rho ) &
           * ( 3.741657386774 ) * cos(6.00 * phi)
      ! n = 6, m = 6
      zn(7) = ( ( 1.00 ) *rho *rho *rho *rho *rho * rho ) &
           * ( 3.741657386774 ) * cos(6.00 * phi)
    case(7)
      ! n = 7, m = -7
      zn(1) = ( ( 1.00 ) *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.000000000000 ) * sin(7.00 * phi)
      ! n = 7, m = -5
      zn(2) = ( ( -6.00 ) *rho *rho *rho *rho * rho + &
           ( 7.00 ) *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.000000000000 ) * sin(7.00 * phi)
      ! n = 7, m = -3
      zn(3) = ( ( 10.00 ) *rho *rho * rho + &
           ( -30.00 ) *rho *rho *rho *rho * rho + &
           ( 21.00 ) *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.000000000000 ) * sin(7.00 * phi)
      ! n = 7, m = -1
      zn(4) = ( ( -4.00 ) * rho + &
           ( 30.00 ) *rho *rho * rho + &
           ( -60.00 ) *rho *rho *rho *rho * rho + &
           ( 35.00 ) *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.000000000000 ) * sin(7.00 * phi)
      ! n = 7, m = 1
      zn(5) = ( ( -4.00 ) * rho + &
           ( 30.00 ) *rho *rho * rho + &
           ( -60.00 ) *rho *rho *rho *rho * rho + &
           ( 35.00 ) *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.000000000000 ) * cos(7.00 * phi)
      ! n = 7, m = 3
      zn(6) = ( ( 10.00 ) *rho *rho * rho + &
           ( -30.00 ) *rho *rho *rho *rho * rho + &
           ( 21.00 ) *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.000000000000 ) * cos(7.00 * phi)
      ! n = 7, m = 5
      zn(7) = ( ( -6.00 ) *rho *rho *rho *rho * rho + &
           ( 7.00 ) *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.000000000000 ) * cos(7.00 * phi)
      ! n = 7, m = 7
      zn(8) = ( ( 1.00 ) *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.000000000000 ) * cos(7.00 * phi)
    case(8)
      ! n = 8, m = -8
      zn(1) = ( ( 1.00 ) *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.242640687119 ) * sin(8.00 * phi)
      ! n = 8, m = -6
      zn(2) = ( ( -7.00 ) *rho *rho *rho *rho *rho * rho + &
           ( 8.00 ) *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.242640687119 ) * sin(8.00 * phi)
      ! n = 8, m = -4
      zn(3) = ( ( 15.00 ) *rho *rho *rho * rho + &
           ( -42.00 ) *rho *rho *rho *rho *rho * rho + &
           ( 28.00 ) *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.242640687119 ) * sin(8.00 * phi)
      ! n = 8, m = -2
      zn(4) = ( ( -10.00 ) *rho * rho + &
           ( 60.00 ) *rho *rho *rho * rho + &
           ( -105.00 ) *rho *rho *rho *rho *rho * rho + &
           ( 56.00 ) *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.242640687119 ) * sin(8.00 * phi)
      ! n = 8, m = 0
      zn(5) = ( ( 1.00 ) + &
           ( -20.00 ) *rho * rho + &
           ( 90.00 ) *rho *rho *rho * rho + &
           ( -140.00 ) *rho *rho *rho *rho *rho * rho + &
           ( 70.00 ) *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 3.000000000000 )
      ! n = 8, m = 2
      zn(6) = ( ( -10.00 ) *rho * rho + &
           ( 60.00 ) *rho *rho *rho * rho + &
           ( -105.00 ) *rho *rho *rho *rho *rho * rho + &
           ( 56.00 ) *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.242640687119 ) * cos(8.00 * phi)
      ! n = 8, m = 4
      zn(7) = ( ( 15.00 ) *rho *rho *rho * rho + &
           ( -42.00 ) *rho *rho *rho *rho *rho * rho + &
           ( 28.00 ) *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.242640687119 ) * cos(8.00 * phi)
      ! n = 8, m = 6
      zn(8) = ( ( -7.00 ) *rho *rho *rho *rho *rho * rho + &
           ( 8.00 ) *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.242640687119 ) * cos(8.00 * phi)
      ! n = 8, m = 8
      zn(9) = ( ( 1.00 ) *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.242640687119 ) * cos(8.00 * phi)
    case(9)
      ! n = 9, m = -9
      zn(1) = ( ( 1.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.472135955000 ) * sin(9.00 * phi)
      ! n = 9, m = -7
      zn(2) = ( ( -8.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( 9.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.472135955000 ) * sin(9.00 * phi)
      ! n = 9, m = -5
      zn(3) = ( ( 21.00 ) *rho *rho *rho *rho * rho + &
           ( -56.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( 36.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.472135955000 ) * sin(9.00 * phi)
      ! n = 9, m = -3
      zn(4) = ( ( -20.00 ) *rho *rho * rho + &
           ( 105.00 ) *rho *rho *rho *rho * rho + &
           ( -168.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( 84.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.472135955000 ) * sin(9.00 * phi)
      ! n = 9, m = -1
      zn(5) = ( ( 5.00 ) * rho + &
           ( -60.00 ) *rho *rho * rho + &
           ( 210.00 ) *rho *rho *rho *rho * rho + &
           ( -280.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( 126.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.472135955000 ) * sin(9.00 * phi)
      ! n = 9, m = 1
      zn(6) = ( ( 5.00 ) * rho + &
           ( -60.00 ) *rho *rho * rho + &
           ( 210.00 ) *rho *rho *rho *rho * rho + &
           ( -280.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( 126.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.472135955000 ) * cos(9.00 * phi)
      ! n = 9, m = 3
      zn(7) = ( ( -20.00 ) *rho *rho * rho + &
           ( 105.00 ) *rho *rho *rho *rho * rho + &
           ( -168.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( 84.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.472135955000 ) * cos(9.00 * phi)
      ! n = 9, m = 5
      zn(8) = ( ( 21.00 ) *rho *rho *rho *rho * rho + &
           ( -56.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( 36.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.472135955000 ) * cos(9.00 * phi)
      ! n = 9, m = 7
      zn(9) = ( ( -8.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( 9.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.472135955000 ) * cos(9.00 * phi)
      ! n = 9, m = 9
      zn(10) = ( ( 1.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.472135955000 ) * cos(9.00 * phi)
    case(10)
      ! n = 10, m = -10
      zn(1) = ( ( 1.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.690415759823 ) * sin(10.00 * phi)
      ! n = 10, m = -8
      zn(2) = ( ( -9.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 10.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.690415759823 ) * sin(10.00 * phi)
      ! n = 10, m = -6
      zn(3) = ( ( 28.00 ) *rho *rho *rho *rho *rho * rho + &
           ( -72.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 45.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.690415759823 ) * sin(10.00 * phi)
      ! n = 10, m = -4
      zn(4) = ( ( -35.00 ) *rho *rho *rho * rho + &
           ( 168.00 ) *rho *rho *rho *rho *rho * rho + &
           ( -252.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 120.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.690415759823 ) * sin(10.00 * phi)
      ! n = 10, m = -2
      zn(5) = ( ( 15.00 ) *rho * rho + &
           ( -140.00 ) *rho *rho *rho * rho + &
           ( 420.00 ) *rho *rho *rho *rho *rho * rho + &
           ( -504.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 210.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.690415759823 ) * sin(10.00 * phi)
      ! n = 10, m = 0
      zn(6) = ( ( -1.00 ) + &
           ( 30.00 ) *rho * rho + &
           ( -210.00 ) *rho *rho *rho * rho + &
           ( 560.00 ) *rho *rho *rho *rho *rho * rho + &
           ( -630.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 252.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 3.316624790355 )
      ! n = 10, m = 2
      zn(7) = ( ( 15.00 ) *rho * rho + &
           ( -140.00 ) *rho *rho *rho * rho + &
           ( 420.00 ) *rho *rho *rho *rho *rho * rho + &
           ( -504.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 210.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.690415759823 ) * cos(10.00 * phi)
      ! n = 10, m = 4
      zn(8) = ( ( -35.00 ) *rho *rho *rho * rho + &
           ( 168.00 ) *rho *rho *rho *rho *rho * rho + &
           ( -252.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 120.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.690415759823 ) * cos(10.00 * phi)
      ! n = 10, m = 6
      zn(9) = ( ( 28.00 ) *rho *rho *rho *rho *rho * rho + &
           ( -72.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 45.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.690415759823 ) * cos(10.00 * phi)
      ! n = 10, m = 8
      zn(10) = ( ( -9.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 10.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.690415759823 ) * cos(10.00 * phi)
      ! n = 10, m = 10
      zn(11) = ( ( 1.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.690415759823 ) * cos(10.00 * phi)
    case(11)
      ! n = 11, m = -11
      zn(1) = ( ( 1.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.898979485566 ) * sin(11.00 * phi)
      ! n = 11, m = -9
      zn(2) = ( ( -10.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 11.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.898979485566 ) * sin(11.00 * phi)
      ! n = 11, m = -7
      zn(3) = ( ( 36.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( -90.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 55.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.898979485566 ) * sin(11.00 * phi)
      ! n = 11, m = -5
      zn(4) = ( ( -56.00 ) *rho *rho *rho *rho * rho + &
           ( 252.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( -360.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 165.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.898979485566 ) * sin(11.00 * phi)
      ! n = 11, m = -3
      zn(5) = ( ( 35.00 ) *rho *rho * rho + &
           ( -280.00 ) *rho *rho *rho *rho * rho + &
           ( 756.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( -840.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 330.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.898979485566 ) * sin(11.00 * phi)
      ! n = 11, m = -1
      zn(6) = ( ( -6.00 ) * rho + &
           ( 105.00 ) *rho *rho * rho + &
           ( -560.00 ) *rho *rho *rho *rho * rho + &
           ( 1260.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( -1260.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 462.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.898979485566 ) * sin(11.00 * phi)
      ! n = 11, m = 1
      zn(7) = ( ( -6.00 ) * rho + &
           ( 105.00 ) *rho *rho * rho + &
           ( -560.00 ) *rho *rho *rho *rho * rho + &
           ( 1260.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( -1260.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 462.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.898979485566 ) * cos(11.00 * phi)
      ! n = 11, m = 3
      zn(8) = ( ( 35.00 ) *rho *rho * rho + &
           ( -280.00 ) *rho *rho *rho *rho * rho + &
           ( 756.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( -840.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 330.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.898979485566 ) * cos(11.00 * phi)
      ! n = 11, m = 5
      zn(9) = ( ( -56.00 ) *rho *rho *rho *rho * rho + &
           ( 252.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( -360.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 165.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.898979485566 ) * cos(11.00 * phi)
      ! n = 11, m = 7
      zn(10) = ( ( 36.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( -90.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 55.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.898979485566 ) * cos(11.00 * phi)
      ! n = 11, m = 9
      zn(11) = ( ( -10.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 11.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.898979485566 ) * cos(11.00 * phi)
      ! n = 11, m = 11
      zn(12) = ( ( 1.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.898979485566 ) * cos(11.00 * phi)
    case(12)
      ! n = 12, m = -12
      zn(1) = ( ( 1.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.099019513593 ) * sin(12.00 * phi)
      ! n = 12, m = -10
      zn(2) = ( ( -11.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 12.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.099019513593 ) * sin(12.00 * phi)
      ! n = 12, m = -8
      zn(3) = ( ( 45.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -110.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 66.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.099019513593 ) * sin(12.00 * phi)
      ! n = 12, m = -6
      zn(4) = ( ( -84.00 ) *rho *rho *rho *rho *rho * rho + &
           ( 360.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -495.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 220.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.099019513593 ) * sin(12.00 * phi)
      ! n = 12, m = -4
      zn(5) = ( ( 70.00 ) *rho *rho *rho * rho + &
           ( -504.00 ) *rho *rho *rho *rho *rho * rho + &
           ( 1260.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -1320.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 495.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.099019513593 ) * sin(12.00 * phi)
      ! n = 12, m = -2
      zn(6) = ( ( -21.00 ) *rho * rho + &
           ( 280.00 ) *rho *rho *rho * rho + &
           ( -1260.00 ) *rho *rho *rho *rho *rho * rho + &
           ( 2520.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -2310.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 792.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.099019513593 ) * sin(12.00 * phi)
      ! n = 12, m = 0
      zn(7) = ( ( 1.00 ) + &
           ( -42.00 ) *rho * rho + &
           ( 420.00 ) *rho *rho *rho * rho + &
           ( -1680.00 ) *rho *rho *rho *rho *rho * rho + &
           ( 3150.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -2772.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 924.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 3.605551275464 )
      ! n = 12, m = 2
      zn(8) = ( ( -21.00 ) *rho * rho + &
           ( 280.00 ) *rho *rho *rho * rho + &
           ( -1260.00 ) *rho *rho *rho *rho *rho * rho + &
           ( 2520.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -2310.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 792.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.099019513593 ) * cos(12.00 * phi)
      ! n = 12, m = 4
      zn(9) = ( ( 70.00 ) *rho *rho *rho * rho + &
           ( -504.00 ) *rho *rho *rho *rho *rho * rho + &
           ( 1260.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -1320.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 495.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.099019513593 ) * cos(12.00 * phi)
      ! n = 12, m = 6
      zn(10) = ( ( -84.00 ) *rho *rho *rho *rho *rho * rho + &
           ( 360.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -495.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 220.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.099019513593 ) * cos(12.00 * phi)
      ! n = 12, m = 8
      zn(11) = ( ( 45.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -110.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 66.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.099019513593 ) * cos(12.00 * phi)
      ! n = 12, m = 10
      zn(12) = ( ( -11.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 12.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.099019513593 ) * cos(12.00 * phi)
      ! n = 12, m = 12
      zn(13) = ( ( 1.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.099019513593 ) * cos(12.00 * phi)
    case(13)
      ! n = 13, m = -13
      zn(1) = ( ( 1.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.291502622129 ) * sin(13.00 * phi)
      ! n = 13, m = -11
      zn(2) = ( ( -12.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 13.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.291502622129 ) * sin(13.00 * phi)
      ! n = 13, m = -9
      zn(3) = ( ( 55.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -132.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 78.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.291502622129 ) * sin(13.00 * phi)
      ! n = 13, m = -7
      zn(4) = ( ( -120.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( 495.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -660.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 286.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.291502622129 ) * sin(13.00 * phi)
      ! n = 13, m = -5
      zn(5) = ( ( 126.00 ) *rho *rho *rho *rho * rho + &
           ( -840.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( 1980.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -1980.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 715.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.291502622129 ) * sin(13.00 * phi)
      ! n = 13, m = -3
      zn(6) = ( ( -56.00 ) *rho *rho * rho + &
           ( 630.00 ) *rho *rho *rho *rho * rho + &
           ( -2520.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( 4620.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -3960.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 1287.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.291502622129 ) * sin(13.00 * phi)
      ! n = 13, m = -1
      zn(7) = ( ( 7.00 ) * rho + &
           ( -168.00 ) *rho *rho * rho + &
           ( 1260.00 ) *rho *rho *rho *rho * rho + &
           ( -4200.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( 6930.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -5544.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 1716.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.291502622129 ) * sin(13.00 * phi)
      ! n = 13, m = 1
      zn(8) = ( ( 7.00 ) * rho + &
           ( -168.00 ) *rho *rho * rho + &
           ( 1260.00 ) *rho *rho *rho *rho * rho + &
           ( -4200.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( 6930.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -5544.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 1716.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.291502622129 ) * cos(13.00 * phi)
      ! n = 13, m = 3
      zn(9) = ( ( -56.00 ) *rho *rho * rho + &
           ( 630.00 ) *rho *rho *rho *rho * rho + &
           ( -2520.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( 4620.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -3960.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 1287.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.291502622129 ) * cos(13.00 * phi)
      ! n = 13, m = 5
      zn(10) = ( ( 126.00 ) *rho *rho *rho *rho * rho + &
           ( -840.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( 1980.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -1980.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 715.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.291502622129 ) * cos(13.00 * phi)
      ! n = 13, m = 7
      zn(11) = ( ( -120.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( 495.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -660.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 286.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.291502622129 ) * cos(13.00 * phi)
      ! n = 13, m = 9
      zn(12) = ( ( 55.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -132.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 78.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.291502622129 ) * cos(13.00 * phi)
      ! n = 13, m = 11
      zn(13) = ( ( -12.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 13.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.291502622129 ) * cos(13.00 * phi)
      ! n = 13, m = 13
      zn(14) = ( ( 1.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.291502622129 ) * cos(13.00 * phi)
    case(14)
      ! n = 14, m = -14
      zn(1) = ( ( 1.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.477225575052 ) * sin(14.00 * phi)
      ! n = 14, m = -12
      zn(2) = ( ( -13.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 14.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.477225575052 ) * sin(14.00 * phi)
      ! n = 14, m = -10
      zn(3) = ( ( 66.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -156.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 91.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.477225575052 ) * sin(14.00 * phi)
      ! n = 14, m = -8
      zn(4) = ( ( -165.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 660.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -858.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 364.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.477225575052 ) * sin(14.00 * phi)
      ! n = 14, m = -6
      zn(5) = ( ( 210.00 ) *rho *rho *rho *rho *rho * rho + &
           ( -1320.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 2970.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -2860.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 1001.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.477225575052 ) * sin(14.00 * phi)
      ! n = 14, m = -4
      zn(6) = ( ( -126.00 ) *rho *rho *rho * rho + &
           ( 1260.00 ) *rho *rho *rho *rho *rho * rho + &
           ( -4620.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 7920.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -6435.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 2002.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.477225575052 ) * sin(14.00 * phi)
      ! n = 14, m = -2
      zn(7) = ( ( 28.00 ) *rho * rho + &
           ( -504.00 ) *rho *rho *rho * rho + &
           ( 3150.00 ) *rho *rho *rho *rho *rho * rho + &
           ( -9240.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 13860.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -10296.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 3003.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.477225575052 ) * sin(14.00 * phi)
      ! n = 14, m = 0
      zn(8) = ( ( -1.00 ) + &
           ( 56.00 ) *rho * rho + &
           ( -756.00 ) *rho *rho *rho * rho + &
           ( 4200.00 ) *rho *rho *rho *rho *rho * rho + &
           ( -11550.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 16632.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -12012.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 3432.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 3.872983346207 )
      ! n = 14, m = 2
      zn(9) = ( ( 28.00 ) *rho * rho + &
           ( -504.00 ) *rho *rho *rho * rho + &
           ( 3150.00 ) *rho *rho *rho *rho *rho * rho + &
           ( -9240.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 13860.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -10296.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 3003.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.477225575052 ) * cos(14.00 * phi)
      ! n = 14, m = 4
      zn(10) = ( ( -126.00 ) *rho *rho *rho * rho + &
           ( 1260.00 ) *rho *rho *rho *rho *rho * rho + &
           ( -4620.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 7920.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -6435.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 2002.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.477225575052 ) * cos(14.00 * phi)
      ! n = 14, m = 6
      zn(11) = ( ( 210.00 ) *rho *rho *rho *rho *rho * rho + &
           ( -1320.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 2970.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -2860.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 1001.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.477225575052 ) * cos(14.00 * phi)
      ! n = 14, m = 8
      zn(12) = ( ( -165.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 660.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -858.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 364.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.477225575052 ) * cos(14.00 * phi)
      ! n = 14, m = 10
      zn(13) = ( ( 66.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -156.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 91.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.477225575052 ) * cos(14.00 * phi)
      ! n = 14, m = 12
      zn(14) = ( ( -13.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 14.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.477225575052 ) * cos(14.00 * phi)
      ! n = 14, m = 14
      zn(15) = ( ( 1.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.477225575052 ) * cos(14.00 * phi)
    case(15)
      ! n = 15, m = -15
      zn(1) = ( ( 1.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.656854249492 ) * sin(15.00 * phi)
      ! n = 15, m = -13
      zn(2) = ( ( -14.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 15.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.656854249492 ) * sin(15.00 * phi)
      ! n = 15, m = -11
      zn(3) = ( ( 78.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -182.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 105.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.656854249492 ) * sin(15.00 * phi)
      ! n = 15, m = -9
      zn(4) = ( ( -220.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 858.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -1092.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 455.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.656854249492 ) * sin(15.00 * phi)
      ! n = 15, m = -7
      zn(5) = ( ( 330.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( -1980.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 4290.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -4004.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 1365.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.656854249492 ) * sin(15.00 * phi)
      ! n = 15, m = -5
      zn(6) = ( ( -252.00 ) *rho *rho *rho *rho * rho + &
           ( 2310.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( -7920.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 12870.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -10010.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 3003.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.656854249492 ) * sin(15.00 * phi)
      ! n = 15, m = -3
      zn(7) = ( ( 84.00 ) *rho *rho * rho + &
           ( -1260.00 ) *rho *rho *rho *rho * rho + &
           ( 6930.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( -18480.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 25740.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -18018.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 5005.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.656854249492 ) * sin(15.00 * phi)
      ! n = 15, m = -1
      zn(8) = ( ( -8.00 ) * rho + &
           ( 252.00 ) *rho *rho * rho + &
           ( -2520.00 ) *rho *rho *rho *rho * rho + &
           ( 11550.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( -27720.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 36036.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -24024.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 6435.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.656854249492 ) * sin(15.00 * phi)
      ! n = 15, m = 1
      zn(9) = ( ( -8.00 ) * rho + &
           ( 252.00 ) *rho *rho * rho + &
           ( -2520.00 ) *rho *rho *rho *rho * rho + &
           ( 11550.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( -27720.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 36036.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -24024.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 6435.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.656854249492 ) * cos(15.00 * phi)
      ! n = 15, m = 3
      zn(10) = ( ( 84.00 ) *rho *rho * rho + &
           ( -1260.00 ) *rho *rho *rho *rho * rho + &
           ( 6930.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( -18480.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 25740.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -18018.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 5005.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.656854249492 ) * cos(15.00 * phi)
      ! n = 15, m = 5
      zn(11) = ( ( -252.00 ) *rho *rho *rho *rho * rho + &
           ( 2310.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( -7920.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 12870.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -10010.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 3003.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.656854249492 ) * cos(15.00 * phi)
      ! n = 15, m = 7
      zn(12) = ( ( 330.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( -1980.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 4290.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -4004.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 1365.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.656854249492 ) * cos(15.00 * phi)
      ! n = 15, m = 9
      zn(13) = ( ( -220.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 858.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -1092.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 455.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.656854249492 ) * cos(15.00 * phi)
      ! n = 15, m = 11
      zn(14) = ( ( 78.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -182.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 105.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.656854249492 ) * cos(15.00 * phi)
      ! n = 15, m = 13
      zn(15) = ( ( -14.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 15.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.656854249492 ) * cos(15.00 * phi)
      ! n = 15, m = 15
      zn(16) = ( ( 1.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.656854249492 ) * cos(15.00 * phi)
    case(16)
      ! n = 16, m = -16
      zn(1) = ( ( 1.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.830951894845 ) * sin(16.00 * phi)
      ! n = 16, m = -14
      zn(2) = ( ( -15.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 16.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.830951894845 ) * sin(16.00 * phi)
      ! n = 16, m = -12
      zn(3) = ( ( 91.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -210.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 120.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.830951894845 ) * sin(16.00 * phi)
      ! n = 16, m = -10
      zn(4) = ( ( -286.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 1092.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -1365.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 560.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.830951894845 ) * sin(16.00 * phi)
      ! n = 16, m = -8
      zn(5) = ( ( 495.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -2860.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 6006.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -5460.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 1820.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.830951894845 ) * sin(16.00 * phi)
      ! n = 16, m = -6
      zn(6) = ( ( -462.00 ) *rho *rho *rho *rho *rho * rho + &
           ( 3960.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -12870.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 20020.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -15015.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 4368.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.830951894845 ) * sin(16.00 * phi)
      ! n = 16, m = -4
      zn(7) = ( ( 210.00 ) *rho *rho *rho * rho + &
           ( -2772.00 ) *rho *rho *rho *rho *rho * rho + &
           ( 13860.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -34320.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 45045.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -30030.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 8008.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.830951894845 ) * sin(16.00 * phi)
      ! n = 16, m = -2
      zn(8) = ( ( -36.00 ) *rho * rho + &
           ( 840.00 ) *rho *rho *rho * rho + &
           ( -6930.00 ) *rho *rho *rho *rho *rho * rho + &
           ( 27720.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -60060.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 72072.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -45045.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 11440.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.830951894845 ) * sin(16.00 * phi)
      ! n = 16, m = 0
      zn(9) = ( ( 1.00 ) + &
           ( -72.00 ) *rho * rho + &
           ( 1260.00 ) *rho *rho *rho * rho + &
           ( -9240.00 ) *rho *rho *rho *rho *rho * rho + &
           ( 34650.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -72072.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 84084.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -51480.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 12870.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.123105625618 )
      ! n = 16, m = 2
      zn(10) = ( ( -36.00 ) *rho * rho + &
           ( 840.00 ) *rho *rho *rho * rho + &
           ( -6930.00 ) *rho *rho *rho *rho *rho * rho + &
           ( 27720.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -60060.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 72072.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -45045.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 11440.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.830951894845 ) * cos(16.00 * phi)
      ! n = 16, m = 4
      zn(11) = ( ( 210.00 ) *rho *rho *rho * rho + &
           ( -2772.00 ) *rho *rho *rho *rho *rho * rho + &
           ( 13860.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -34320.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 45045.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -30030.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 8008.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.830951894845 ) * cos(16.00 * phi)
      ! n = 16, m = 6
      zn(12) = ( ( -462.00 ) *rho *rho *rho *rho *rho * rho + &
           ( 3960.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -12870.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 20020.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -15015.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 4368.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.830951894845 ) * cos(16.00 * phi)
      ! n = 16, m = 8
      zn(13) = ( ( 495.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -2860.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 6006.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -5460.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 1820.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.830951894845 ) * cos(16.00 * phi)
      ! n = 16, m = 10
      zn(14) = ( ( -286.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 1092.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -1365.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 560.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.830951894845 ) * cos(16.00 * phi)
      ! n = 16, m = 12
      zn(15) = ( ( 91.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -210.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 120.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.830951894845 ) * cos(16.00 * phi)
      ! n = 16, m = 14
      zn(16) = ( ( -15.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 16.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.830951894845 ) * cos(16.00 * phi)
      ! n = 16, m = 16
      zn(17) = ( ( 1.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 5.830951894845 ) * cos(16.00 * phi)
    case(17)
      ! n = 17, m = -17
      zn(1) = ( ( 1.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 6.000000000000 ) * sin(17.00 * phi)
      ! n = 17, m = -15
      zn(2) = ( ( -16.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 17.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 6.000000000000 ) * sin(17.00 * phi)
      ! n = 17, m = -13
      zn(3) = ( ( 105.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -240.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 136.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 6.000000000000 ) * sin(17.00 * phi)
      ! n = 17, m = -11
      zn(4) = ( ( -364.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 1365.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -1680.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 680.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 6.000000000000 ) * sin(17.00 * phi)
      ! n = 17, m = -9
      zn(5) = ( ( 715.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -4004.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 8190.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -7280.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 2380.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 6.000000000000 ) * sin(17.00 * phi)
      ! n = 17, m = -7
      zn(6) = ( ( -792.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( 6435.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -20020.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 30030.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -21840.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 6188.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 6.000000000000 ) * sin(17.00 * phi)
      ! n = 17, m = -5
      zn(7) = ( ( 462.00 ) *rho *rho *rho *rho * rho + &
           ( -5544.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( 25740.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -60060.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 75075.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -48048.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 12376.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 6.000000000000 ) * sin(17.00 * phi)
      ! n = 17, m = -3
      zn(8) = ( ( -120.00 ) *rho *rho * rho + &
           ( 2310.00 ) *rho *rho *rho *rho * rho + &
           ( -16632.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( 60060.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -120120.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 135135.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -80080.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 19448.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 6.000000000000 ) * sin(17.00 * phi)
      ! n = 17, m = -1
      zn(9) = ( ( 9.00 ) * rho + &
           ( -360.00 ) *rho *rho * rho + &
           ( 4620.00 ) *rho *rho *rho *rho * rho + &
           ( -27720.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( 90090.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -168168.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 180180.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -102960.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 24310.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 6.000000000000 ) * sin(17.00 * phi)
      ! n = 17, m = 1
      zn(10) = ( ( 9.00 ) * rho + &
           ( -360.00 ) *rho *rho * rho + &
           ( 4620.00 ) *rho *rho *rho *rho * rho + &
           ( -27720.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( 90090.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -168168.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 180180.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -102960.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 24310.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 6.000000000000 ) * cos(17.00 * phi)
      ! n = 17, m = 3
      zn(11) = ( ( -120.00 ) *rho *rho * rho + &
           ( 2310.00 ) *rho *rho *rho *rho * rho + &
           ( -16632.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( 60060.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -120120.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 135135.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -80080.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 19448.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 6.000000000000 ) * cos(17.00 * phi)
      ! n = 17, m = 5
      zn(12) = ( ( 462.00 ) *rho *rho *rho *rho * rho + &
           ( -5544.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( 25740.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -60060.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 75075.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -48048.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 12376.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 6.000000000000 ) * cos(17.00 * phi)
      ! n = 17, m = 7
      zn(13) = ( ( -792.00 ) *rho *rho *rho *rho *rho *rho * rho + &
           ( 6435.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -20020.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 30030.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -21840.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 6188.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 6.000000000000 ) * cos(17.00 * phi)
      ! n = 17, m = 9
      zn(14) = ( ( 715.00 ) *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -4004.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 8190.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -7280.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 2380.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 6.000000000000 ) * cos(17.00 * phi)
      ! n = 17, m = 11
      zn(15) = ( ( -364.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 1365.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -1680.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 680.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 6.000000000000 ) * cos(17.00 * phi)
      ! n = 17, m = 13
      zn(16) = ( ( 105.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -240.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 136.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 6.000000000000 ) * cos(17.00 * phi)
      ! n = 17, m = 15
      zn(17) = ( ( -16.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 17.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 6.000000000000 ) * cos(17.00 * phi)
      ! n = 17, m = 17
      zn(18) = ( ( 1.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 6.000000000000 ) * cos(17.00 * phi)
    case(18)
      ! n = 18, m = -18
      zn(1) = ( ( 1.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 6.164414002969 ) * sin(18.00 * phi)
      ! n = 18, m = -16
      zn(2) = ( ( -17.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 18.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 6.164414002969 ) * sin(18.00 * phi)
      ! n = 18, m = -14
      zn(3) = ( ( 120.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -272.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 153.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 6.164414002969 ) * sin(18.00 * phi)
      ! n = 18, m = -12
      zn(4) = ( ( -455.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 1680.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -2040.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 816.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 6.164414002969 ) * sin(18.00 * phi)
      ! n = 18, m = -10
      zn(5) = ( ( 1001.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -5460.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 10920.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -9520.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 3060.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 6.164414002969 ) * sin(18.00 * phi)
      ! n = 18, m = -8
      zn(6) = ( ( -1287.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 10010.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -30030.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 43680.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -30940.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 8568.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 6.164414002969 ) * sin(18.00 * phi)
      ! n = 18, m = -6
      zn(7) = ( ( 924.00 ) *rho *rho *rho *rho *rho * rho + &
           ( -10296.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 45045.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -100100.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 120120.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -74256.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 18564.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 6.164414002969 ) * sin(18.00 * phi)
      ! n = 18, m = -4
      zn(8) = ( ( -330.00 ) *rho *rho *rho * rho + &
           ( 5544.00 ) *rho *rho *rho *rho *rho * rho + &
           ( -36036.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 120120.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -225225.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 240240.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -136136.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 31824.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 6.164414002969 ) * sin(18.00 * phi)
      ! n = 18, m = -2
      zn(9) = ( ( 45.00 ) *rho * rho + &
           ( -1320.00 ) *rho *rho *rho * rho + &
           ( 13860.00 ) *rho *rho *rho *rho *rho * rho + &
           ( -72072.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 210210.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -360360.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 360360.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -194480.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 43758.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 6.164414002969 ) * sin(18.00 * phi)
      ! n = 18, m = 0
      zn(10) = ( ( -1.00 ) + &
           ( 90.00 ) *rho * rho + &
           ( -1980.00 ) *rho *rho *rho * rho + &
           ( 18480.00 ) *rho *rho *rho *rho *rho * rho + &
           ( -90090.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 252252.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -420420.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 411840.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -218790.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 48620.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 4.358898943541 )
      ! n = 18, m = 2
      zn(11) = ( ( 45.00 ) *rho * rho + &
           ( -1320.00 ) *rho *rho *rho * rho + &
           ( 13860.00 ) *rho *rho *rho *rho *rho * rho + &
           ( -72072.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 210210.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -360360.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 360360.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -194480.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 43758.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 6.164414002969 ) * cos(18.00 * phi)
      ! n = 18, m = 4
      zn(12) = ( ( -330.00 ) *rho *rho *rho * rho + &
           ( 5544.00 ) *rho *rho *rho *rho *rho * rho + &
           ( -36036.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 120120.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -225225.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 240240.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -136136.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 31824.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 6.164414002969 ) * cos(18.00 * phi)
      ! n = 18, m = 6
      zn(13) = ( ( 924.00 ) *rho *rho *rho *rho *rho * rho + &
           ( -10296.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 45045.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -100100.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 120120.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -74256.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 18564.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 6.164414002969 ) * cos(18.00 * phi)
      ! n = 18, m = 8
      zn(14) = ( ( -1287.00 ) *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 10010.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -30030.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 43680.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -30940.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 8568.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 6.164414002969 ) * cos(18.00 * phi)
      ! n = 18, m = 10
      zn(15) = ( ( 1001.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -5460.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 10920.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -9520.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 3060.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 6.164414002969 ) * cos(18.00 * phi)
      ! n = 18, m = 12
      zn(16) = ( ( -455.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 1680.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -2040.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 816.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 6.164414002969 ) * cos(18.00 * phi)
      ! n = 18, m = 14
      zn(17) = ( ( 120.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( -272.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 153.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 6.164414002969 ) * cos(18.00 * phi)
      ! n = 18, m = 16
      zn(18) = ( ( -17.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho + &
           ( 18.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 6.164414002969 ) * cos(18.00 * phi)
      ! n = 18, m = 18
      zn(19) = ( ( 1.00 ) *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho *rho * rho ) &
           * ( 6.164414002969 ) * cos(18.00 * phi)
    case default
      zn = ONE
    end select

  end function calc_zn

!===============================================================================
! CALC_ZN_SCALED calculates the n-th order Zernike polynomial moment for a given
! angle (rho, theta) location in the unit disk, scaled correctly for orthogonal
! integration.
!===============================================================================

  pure function calc_zn_scaled(n, rho, phi) result(zn)

    integer, intent(in) :: n           ! Order requested
    real(8), intent(in) :: rho         ! Radial location in the unit disk
    real(8), intent(in) :: phi         ! Theta (radians) location in the unit disk
    real(8)             :: zn(n + 1)   ! The resultant Z_n(uvw)

    zn = calc_zn(n, rho, phi) / SQRT_PI

  end function calc_zn_scaled

!===============================================================================
!                               C API FUNCTIONS
!===============================================================================


end module tally_filter_zernike
