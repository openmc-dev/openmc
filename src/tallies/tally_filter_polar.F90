module tally_filter_polar

  use, intrinsic :: ISO_C_BINDING

  use hdf5

  use algorithm,          only: binary_search
  use constants
  use error,              only: fatal_error
  use hdf5_interface
  use particle_header,    only: Particle
  use string,             only: to_str
  use tally_filter_header
  use xml_interface

  implicit none
  private

!===============================================================================
! POLARFILTER bins the incident neutron polar angle (relative to the global
! z-axis).
!===============================================================================

  type, public, extends(TallyFilter) :: PolarFilter
    real(8), allocatable :: bins(:)
  contains
    procedure :: from_xml
    procedure :: get_all_bins => get_all_bins_polar
    procedure :: to_statepoint => to_statepoint_polar
    procedure :: text_label => text_label_polar
  end type PolarFilter

contains

  subroutine from_xml(this, node)
    class(PolarFilter), intent(inout) :: this
    type(XMLNode), intent(in) :: node

    integer :: i
    integer :: n_angle
    integer :: n
    real(8) :: d_angle

    n = node_word_count(node, "bins")

    ! Allocate and store bins
    this % n_bins = n - 1
    allocate(this % bins(n))
    call get_node_array(node, "bins", this % bins)

    ! Allow a user to input a lone number which will mean that you
    ! subdivide [0,pi] evenly with the input being the number of bins
    if (n == 1) then
      n_angle = int(this % bins(1))
      if (n_angle > 1) then
        this % n_bins = n_angle
        d_angle = PI / real(n_angle,8)
        deallocate(this % bins)
        allocate(this % bins(n_angle + 1))
        do i = 1, n_angle
          this % bins(i) = (i - 1) * d_angle
        end do
        this % bins(n_angle + 1) = PI
      else
        call fatal_error("Number of bins for polar filter must be&
             & greater than 1.")
      end if
    end if
  end subroutine from_xml

  subroutine get_all_bins_polar(this, p, estimator, match)
    class(PolarFilter), intent(in)  :: this
    type(Particle),     intent(in)  :: p
    integer,            intent(in)  :: estimator
    type(TallyFilterMatch),  intent(inout) :: match

    integer :: n
    integer :: bin
    real(8) :: theta

    n = this % n_bins

    ! Make sure the correct direction vector is used.
    if (estimator == ESTIMATOR_TRACKLENGTH) then
      theta = acos(p % coord(1) % uvw(3))
    else
      theta = acos(p % last_uvw(3))
    end if

    ! Search to find polar angle bin.
    bin = binary_search(this % bins, n + 1, theta)
    if (bin /= NO_BIN_FOUND) then
      call match % bins % push_back(bin)
      call match % weights % push_back(ONE)
    end if
  end subroutine get_all_bins_polar

  subroutine to_statepoint_polar(this, filter_group)
    class(PolarFilter), intent(in) :: this
    integer(HID_T),     intent(in) :: filter_group

    call write_dataset(filter_group, "type", "polar")
    call write_dataset(filter_group, "n_bins", this % n_bins)
    call write_dataset(filter_group, "bins", this % bins)
  end subroutine to_statepoint_polar

  function text_label_polar(this, bin) result(label)
    class(PolarFilter), intent(in) :: this
    integer,            intent(in) :: bin
    character(MAX_LINE_LEN)        :: label

    real(8) :: E0, E1

    E0 = this % bins(bin)
    E1 = this % bins(bin + 1)
    label = "Polar Angle [" // trim(to_str(E0)) // ", " // trim(to_str(E1)) &
         // ")"
  end function text_label_polar

end module tally_filter_polar
