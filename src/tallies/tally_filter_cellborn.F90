module tally_filter_cellborn

  use, intrinsic :: ISO_C_BINDING

  use hdf5

  use constants,          only: ONE, MAX_LINE_LEN
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
! CELLBORNFILTER specifies which cell the particle was born in.
!===============================================================================

  type, public, extends(TallyFilter) :: CellbornFilter
    integer, allocatable :: cells(:)
    type(DictIntInt)     :: map
  contains
    procedure :: from_xml
    procedure :: get_all_bins => get_all_bins_cellborn
    procedure :: to_statepoint => to_statepoint_cellborn
    procedure :: text_label => text_label_cellborn
    procedure :: initialize => initialize_cellborn
  end type CellbornFilter

contains

  subroutine from_xml(this, node)
    class(CellbornFilter), intent(inout) :: this
    type(XMLNode), intent(in) :: node

    integer :: n

    ! Determine number of bins
    n = node_word_count(node, "bins")

    ! Allocate and store bins
    this % n_bins = n
    allocate(this % cells(n))
    call get_node_array(node, "bins", this % cells)
  end subroutine from_xml

  subroutine get_all_bins_cellborn(this, p, estimator, match)
    class(CellbornFilter), intent(in)  :: this
    type(Particle),        intent(in)  :: p
    integer,               intent(in)  :: estimator
    type(TallyFilterMatch),     intent(inout) :: match

      if (this % map % has_key(p % cell_born)) then
        call match % bins % push_back(this % map % get_key(p % cell_born))
        call match % weights % push_back(ONE)
      end if

  end subroutine get_all_bins_cellborn

  subroutine to_statepoint_cellborn(this, filter_group)
    class(CellbornFilter), intent(in) :: this
    integer(HID_T),        intent(in) :: filter_group

    integer :: i
    integer, allocatable :: cell_ids(:)

    call write_dataset(filter_group, "type", "cellborn")
    call write_dataset(filter_group, "n_bins", this % n_bins)
    allocate(cell_ids(size(this % cells)))
    do i = 1, size(this % cells)
      cell_ids(i) = cells(this % cells(i)) % id
    end do
    call write_dataset(filter_group, "bins", cell_ids)
  end subroutine to_statepoint_cellborn

  subroutine initialize_cellborn(this)
    class(CellbornFilter), intent(inout) :: this

    integer :: i, id

    ! Convert ids to indices.
    do i = 1, this % n_bins
      id = this % cells(i)
      if (cell_dict % has_key(id)) then
        this % cells(i) = cell_dict % get_key(id)
      else
        call fatal_error("Could not find cell " // trim(to_str(id)) &
             &// " specified on tally filter.")
      end if
    end do

    ! Generate mapping from cell indices to filter bins.
    do i = 1, this % n_bins
      call this % map % add_key(this % cells(i), i)
    end do
  end subroutine initialize_cellborn

  function text_label_cellborn(this, bin) result(label)
    class(CellbornFilter), intent(in) :: this
    integer,               intent(in) :: bin
    character(MAX_LINE_LEN)           :: label

    label = "Birth Cell " // to_str(cells(this % cells(bin)) % id)
  end function text_label_cellborn

end module tally_filter_cellborn
