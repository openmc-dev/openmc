module tally_filter_cellfrom

  use, intrinsic :: ISO_C_BINDING

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

    call this % get_all_bins_c(p, estimator, match)

  end subroutine get_all_bins_cell_from

  subroutine to_statepoint_cell_from(this, filter_group)
    class(CellFromFilter), intent(in) :: this
    integer(HID_T),        intent(in) :: filter_group

    call this % to_statepoint_c(filter_group)

  end subroutine to_statepoint_cell_from

  function text_label_cell_from(this, bin) result(label)
    class(CellFromFilter), intent(in) :: this
    integer,               intent(in) :: bin
    character(MAX_LINE_LEN)           :: label

    label = this % text_label_c(bin)
  end function text_label_cell_from

end module tally_filter_cellfrom
