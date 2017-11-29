module tally_filter_cellfrom

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
  use tally_filter_cell
  use xml_interface

  implicit none
  private

!===============================================================================
! CELLFROMFILTER specifies which geometric cells particles exit when crossing a
! surface.
!===============================================================================

  type, public, extends(CellFilter) :: CellFromFilter
  contains
    ! Inherit from_xml from CellFilter
    procedure :: get_all_bins => get_all_bins_cell_from
    procedure :: to_statepoint => to_statepoint_cell_from
    procedure :: text_label => text_label_cell_from
  end type CellFromFilter

contains

  subroutine get_all_bins_cell_from(this, p, estimator, match)
    class(CellFromFilter), intent(in)  :: this
    type(Particle),    intent(in)  :: p
    integer,           intent(in)  :: estimator
    type(TallyFilterMatch), intent(inout) :: match

    integer :: i
    integer :: val

    ! Starting one coordinate level deeper, find the next bin.
    do i = 1, p % last_n_coord
      val = this % map % get(p % last_cell(i))
      if (val /= EMPTY) then
        call match % bins % push_back(val)
        call match % weights % push_back(ONE)
        exit
      end if
    end do

  end subroutine get_all_bins_cell_from

  subroutine to_statepoint_cell_from(this, filter_group)
    class(CellFromFilter), intent(in) :: this
    integer(HID_T),        intent(in) :: filter_group

    integer :: i
    integer, allocatable :: cell_ids(:)

    call write_dataset(filter_group, "type", "cellfrom")
    call write_dataset(filter_group, "n_bins", this % n_bins)

    allocate(cell_ids(size(this % cells)))
    do i = 1, size(this % cells)
      cell_ids(i) = cells(this % cells(i)) % id
    end do
    call write_dataset(filter_group, "bins", cell_ids)
  end subroutine to_statepoint_cell_from

  function text_label_cell_from(this, bin) result(label)
    class(CellFromFilter), intent(in) :: this
    integer,               intent(in) :: bin
    character(MAX_LINE_LEN)           :: label

    label = "Cell from " // to_str(cells(this % cells(bin)) % id)
  end function text_label_cell_from

end module tally_filter_cellfrom
