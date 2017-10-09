module tally_filter_universe

  use, intrinsic :: ISO_C_BINDING

  use hdf5

  use constants,          only: ONE, MAX_LINE_LEN
  use dict_header,        only: EMPTY
  use error,              only: fatal_error
  use hdf5_interface
  use geometry_header
  use particle_header,    only: Particle
  use string,             only: to_str
  use tally_filter_header
  use xml_interface

  implicit none
  private

!===============================================================================
! UNIVERSEFILTER specifies which geometric universes tally events reside in.
!===============================================================================

  type, public, extends(TallyFilter) :: UniverseFilter
    integer, allocatable :: universes(:)
    type(DictIntInt)     :: map
  contains
    procedure :: from_xml
    procedure :: get_all_bins => get_all_bins_universe
    procedure :: to_statepoint => to_statepoint_universe
    procedure :: text_label => text_label_universe
    procedure :: initialize => initialize_universe
  end type UniverseFilter

contains

  subroutine from_xml(this, node)
    class(UniverseFilter), intent(inout) :: this
    type(XMLNode), intent(in) :: node

    integer :: n

    n = node_word_count(node, "bins")

    ! Allocate and store bins
    this % n_bins = n
    allocate(this % universes(n))
    call get_node_array(node, "bins", this % universes)
  end subroutine from_xml

  subroutine get_all_bins_universe(this, p, estimator, match)
    class(UniverseFilter), intent(in)  :: this
    type(Particle),        intent(in)  :: p
    integer,               intent(in)  :: estimator
    type(TallyFilterMatch),     intent(inout) :: match

    integer :: i
    integer :: val

    ! Iterate over coordinate levels to see which universes match
    do i = 1, p % n_coord
      val = this % map % get(p % coord(i) % universe)
      if (val /= EMPTY) then
        call match % bins % push_back(val)
        call match % weights % push_back(ONE)
      end if
    end do

  end subroutine get_all_bins_universe

  subroutine to_statepoint_universe(this, filter_group)
    class(UniverseFilter), intent(in) :: this
    integer(HID_T),        intent(in) :: filter_group

    integer :: i
    integer, allocatable :: universe_ids(:)

    call write_dataset(filter_group, "type", "universe")
    call write_dataset(filter_group, "n_bins", this % n_bins)

    allocate(universe_ids(size(this % universes)))
    do i = 1, size(this % universes)
      universe_ids(i) = universes(this % universes(i)) % id
    end do
    call write_dataset(filter_group, "bins", universe_ids)
  end subroutine to_statepoint_universe

  subroutine initialize_universe(this)
    class(UniverseFilter), intent(inout) :: this

    integer :: i, id
    integer :: val

    ! Convert ids to indices.
    do i = 1, this % n_bins
      id = this % universes(i)
      val = universe_dict % get(id)
      if (val /= EMPTY) then
        this % universes(i) = val
      else
        call fatal_error("Could not find universe " // trim(to_str(id)) &
             &// " specified on a tally filter.")
      end if
    end do

    ! Generate mapping from universe indices to filter bins.
    do i = 1, this % n_bins
      call this % map % set(this % universes(i), i)
    end do
  end subroutine initialize_universe

  function text_label_universe(this, bin) result(label)
    class(UniverseFilter), intent(in) :: this
    integer,               intent(in) :: bin
    character(MAX_LINE_LEN)           :: label

    label = "Universe " // to_str(universes(this % universes(bin)) % id)
  end function text_label_universe

end module tally_filter_universe
