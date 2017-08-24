module tally_filter_azimuthal

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
! AZIMUTHALFILTER bins the incident neutron azimuthal angle (relative to the
! global xy-plane).
!===============================================================================

  type, public, extends(TallyFilter) :: AzimuthalFilter
    real(8), allocatable :: bins(:)
  contains
    procedure :: from_xml
    procedure :: get_all_bins => get_all_bins_azimuthal
    procedure :: to_statepoint => to_statepoint_azimuthal
    procedure :: text_label => text_label_azimuthal
  end type AzimuthalFilter

contains

  subroutine from_xml(this, node)
    class(AzimuthalFilter), intent(inout) :: this
    type(XMLNode), intent(in) :: node

    integer :: i
    integer :: n
    integer :: n_angle
    real(8) :: d_angle

    n = node_word_count(node, "bins")

    ! Allocate and store bins
    this % n_bins = n - 1
    allocate(this % bins(n))
    call get_node_array(node, "bins", this % bins)

    ! Allow a user to input a lone number which will mean that you
    ! subdivide [-pi,pi) evenly with the input being the number of
    ! bins
    if (n == 1) then
      n_angle = int(this % bins(1))
      if (n_angle > 1) then
        this % n_bins = n_angle
        d_angle = TWO * PI / n_angle
        deallocate(this % bins)
        allocate(this % bins(n_angle + 1))
        do i = 1, n_angle
          this % bins(i) = -PI + (i - 1) * d_angle
        end do
        this % bins(n_angle + 1) = PI
      else
        call fatal_error("Number of bins for azimuthal filter must be&
             & greater than 1.")
      end if
    end if
  end subroutine from_xml

  subroutine get_all_bins_azimuthal(this, p, estimator, match)
    class(AzimuthalFilter), intent(in)  :: this
    type(Particle),         intent(in)  :: p
    integer,                intent(in)  :: estimator
    type(TallyFilterMatch),      intent(inout) :: match

    integer :: n
    integer :: bin
    real(8) :: phi

    n = this % n_bins

    ! Make sure the correct direction vector is used.
    if (estimator == ESTIMATOR_TRACKLENGTH) then
      phi = atan2(p % coord(1) % uvw(2), p % coord(1) % uvw(1))
    else
      phi = atan2(p % last_uvw(2), p % last_uvw(1))
    end if

    ! Search to find azimuthal angle bin.
    bin = binary_search(this % bins, n + 1, phi)
    if (bin /= NO_BIN_FOUND) then
      call match % bins % push_back(bin)
      call match % weights % push_back(ONE)
    end if

  end subroutine get_all_bins_azimuthal

  subroutine to_statepoint_azimuthal(this, filter_group)
    class(AzimuthalFilter), intent(in) :: this
    integer(HID_T),         intent(in) :: filter_group

    call write_dataset(filter_group, "type", "azimuthal")
    call write_dataset(filter_group, "n_bins", this % n_bins)
    call write_dataset(filter_group, "bins", this % bins)
  end subroutine to_statepoint_azimuthal

  function text_label_azimuthal(this, bin) result(label)
    class(AzimuthalFilter), intent(in) :: this
    integer,                intent(in) :: bin
    character(MAX_LINE_LEN)            :: label

    real(8) :: E0, E1

    E0 = this % bins(bin)
    E1 = this % bins(bin + 1)
    label = "Azimuthal Angle [" // trim(to_str(E0)) // ", " &
         // trim(to_str(E1)) // ")"
  end function text_label_azimuthal

end module tally_filter_azimuthal
