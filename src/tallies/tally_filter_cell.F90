module tally_filter_cell

  use, intrinsic :: ISO_C_BINDING

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
  public :: openmc_cell_filter_get_bins

!===============================================================================
! CELLFILTER specifies which geometric cells tally events reside in.
!===============================================================================

  type, public, extends(CppTallyFilter) :: CellFilter
    integer(C_INT32_T), allocatable :: cells(:)
    type(DictIntInt)     :: map
  contains
    procedure :: from_xml
    procedure :: get_all_bins => get_all_bins_cell
    procedure :: to_statepoint => to_statepoint_cell
    procedure :: text_label => text_label_cell
    procedure :: initialize => initialize_cell
  end type CellFilter

contains

  subroutine from_xml(this, node)
    class(CellFilter), intent(inout) :: this
    type(XMLNode), intent(in) :: node

    integer :: n

    call this % from_xml_c(node)

    ! Determine how many bins were given
    n = node_word_count(node, "bins")

    ! Allocate and store bins
    this % n_bins = n
    allocate(this % cells(n))
    call get_node_array(node, "bins", this % cells)
  end subroutine from_xml

  subroutine get_all_bins_cell(this, p, estimator, match)
    class(CellFilter), intent(in)  :: this
    type(Particle),    intent(in)  :: p
    integer,           intent(in)  :: estimator
    type(TallyFilterMatch), intent(inout) :: match

    call this % get_all_bins_c(p, estimator, match)

  end subroutine get_all_bins_cell

  subroutine to_statepoint_cell(this, filter_group)
    class(CellFilter), intent(in) :: this
    integer(HID_T),    intent(in) :: filter_group

    integer :: i
    integer, allocatable :: cell_ids(:)

    call write_dataset(filter_group, "type", "cell")
    call write_dataset(filter_group, "n_bins", this % n_bins)

    allocate(cell_ids(size(this % cells)))
    do i = 1, size(this % cells)
      cell_ids(i) = cells(this % cells(i)) % id()
    end do
    call write_dataset(filter_group, "bins", cell_ids)
  end subroutine to_statepoint_cell

  subroutine initialize_cell(this)
    class(CellFilter), intent(inout) :: this

    integer :: i, id
    integer :: val

    call this % initialize_c()

    ! Convert ids to indices.
    do i = 1, this % n_bins
      id = this % cells(i)
      val = cell_dict % get(id)
      if (val /= EMPTY) then
        this % cells(i) = val
      else
        call fatal_error("Could not find cell " // trim(to_str(id)) &
             &// " specified on tally filter.")
      end if
    end do

    ! Generate mapping from cell indices to filter bins.
    do i = 1, this % n_bins
      call this % map % set(this % cells(i), i)
    end do
  end subroutine initialize_cell

  function text_label_cell(this, bin) result(label)
    class(CellFilter), intent(in) :: this
    integer,           intent(in) :: bin
    character(MAX_LINE_LEN)       :: label

    label = "Cell " // to_str(cells(this % cells(bin)) % id())
  end function text_label_cell

!===============================================================================
!                               C API FUNCTIONS
!===============================================================================

  function openmc_cell_filter_get_bins(index, cells, n) result(err) bind(C)
    ! Return the cells associated with a cell filter
    integer(C_INT32_T), value :: index
    type(C_PTR), intent(out) :: cells
    integer(C_INT32_T), intent(out) :: n
    integer(C_INT) :: err

    err = verify_filter(index)
    if (err == 0) then
       select type (f => filters(index) % obj)
       type is (CellFilter)
          cells = C_LOC(f % cells)
          n = size(f % cells)
       class default
          err = E_INVALID_TYPE
          call set_errmsg("Tried to get cells from a non-cell filter.")
       end select
    end if
  end function openmc_cell_filter_get_bins

end module tally_filter_cell
