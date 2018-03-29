module tally_filter_surface

  use, intrinsic :: ISO_C_BINDING

  use hdf5

  use constants,          only: ONE, MAX_LINE_LEN
  use dict_header,        only: EMPTY
  use error,              only: fatal_error
  use hdf5_interface
  use surface_header
  use particle_header,    only: Particle
  use string,             only: to_str
  use tally_filter_header
  use xml_interface

  implicit none
  private

!===============================================================================
! SURFACEFILTER specifies which surface particles are crossing
!===============================================================================

  type, public, extends(TallyFilter) :: SurfaceFilter
    integer, allocatable :: surfaces(:)

    ! True if this filter is used for surface currents
    logical              :: current = .false.
  contains
    procedure :: from_xml
    procedure :: get_all_bins => get_all_bins_surface
    procedure :: to_statepoint => to_statepoint_surface
    procedure :: text_label => text_label_surface
    procedure :: initialize => initialize_surface
  end type SurfaceFilter

contains

  subroutine from_xml(this, node)
    class(SurfaceFilter), intent(inout) :: this
    type(XMLNode), intent(in) :: node

    integer :: n

    n = node_word_count(node, "bins")

    ! Allocate and store bins
    this % n_bins = n
    allocate(this % surfaces(n))
    call get_node_array(node, "bins", this % surfaces)
  end subroutine from_xml

  subroutine get_all_bins_surface(this, p, estimator, match)
    class(SurfaceFilter), intent(in)  :: this
    type(Particle),       intent(in)  :: p
    integer,              intent(in)  :: estimator
    type(TallyFilterMatch),    intent(inout) :: match

    integer :: i

      do i = 1, this % n_bins
        if (abs(p % surface) == this % surfaces(i)) then
          call match % bins % push_back(i)
          if (p % surface < 0) then
            call match % weights % push_back(-ONE)
          else
            call match % weights % push_back(ONE)
          end if
          exit
        end if
      end do

  end subroutine get_all_bins_surface

  subroutine to_statepoint_surface(this, filter_group)
    class(SurfaceFilter), intent(in) :: this
    integer(HID_T),       intent(in) :: filter_group

    call write_dataset(filter_group, "type", "surface")
    call write_dataset(filter_group, "n_bins", this % n_bins)
    call write_dataset(filter_group, "bins", this % surfaces)
  end subroutine to_statepoint_surface

  subroutine initialize_surface(this)
    class(SurfaceFilter), intent(inout) :: this

    integer :: i, id
    integer :: val

    ! Convert ids to indices.
    do i = 1, this % n_bins
      id = this % surfaces(i)
      val = surface_dict % get(id)
      if (val /= EMPTY) then
        this % surfaces(i) = val
      else
        call fatal_error("Could not find surface " // trim(to_str(id)) &
             &// " specified on tally filter.")
      end if
    end do
  end subroutine initialize_surface

  function text_label_surface(this, bin) result(label)
    class(SurfaceFilter), intent(in) :: this
    integer,              intent(in) :: bin
    character(MAX_LINE_LEN)          :: label

    label = "Surface " // to_str(surfaces(this % surfaces(bin)) % id())
  end function text_label_surface

end module tally_filter_surface
