module tally_filter_distribcell

  use, intrinsic :: ISO_C_BINDING

  use constants
  use dict_header,     only: EMPTY
  use error
  use geometry_header
  use hdf5_interface
  use particle_header, only: Particle
  use string,          only: to_str
  use tally_filter_header
  use xml_interface

  implicit none
  private

!===============================================================================
! DISTRIBCELLFILTER specifies which distributed geometric cells tally events
! reside in.
!===============================================================================

  type, public, extends(TallyFilter) :: DistribcellFilter
    integer :: cell
  contains
    procedure :: from_xml
    procedure :: get_all_bins => get_all_bins_distribcell
    procedure :: to_statepoint => to_statepoint_distribcell
    procedure :: text_label => text_label_distribcell
    procedure :: initialize => initialize_distribcell
  end type DistribcellFilter

contains

  subroutine from_xml(this, node)
    class(DistribcellFilter), intent(inout) :: this
    type(XMLNode), intent(in) :: node

    integer :: n

    n = node_word_count(node, "bins")
    if (n /= 1) call fatal_error("Only one cell can be &
         &specified per distribcell filter.")

    ! Store bins
    call get_node_value(node, "bins", this % cell)
  end subroutine from_xml

  subroutine get_all_bins_distribcell(this, p, estimator, match)
    class(DistribcellFilter), intent(in)  :: this
    type(Particle),           intent(in)  :: p
    integer,                  intent(in)  :: estimator
    type(TallyFilterMatch),   intent(inout) :: match

    integer :: distribcell_index, offset, i

      distribcell_index = cells(this % cell) % distribcell_index()
      offset = 0
      do i = 1, p % n_coord
        if (cells(p % coord(i) % cell) % type() == FILL_UNIVERSE) then
          offset = offset + cells(p % coord(i) % cell) &
               % offset(distribcell_index-1)
        elseif (cells(p % coord(i) % cell) % type() == FILL_LATTICE) then
          if (lattices(p % coord(i + 1) % lattice) % obj &
               % are_valid_indices([&
               p % coord(i + 1) % lattice_x, &
               p % coord(i + 1) % lattice_y, &
               p % coord(i + 1) % lattice_z])) then
            offset = offset + lattices(p % coord(i + 1) % lattice) % obj &
                 % offset(distribcell_index - 1, &
                 [p % coord(i + 1) % lattice_x, &
                 p % coord(i + 1) % lattice_y, &
                 p % coord(i + 1) % lattice_z])
          end if
        end if
        if (this % cell == p % coord(i) % cell) then
          call match % bins % push_back(offset + 1)
          call match % weights % push_back(ONE)
          return
        end if
      end do
  end subroutine get_all_bins_distribcell

  subroutine to_statepoint_distribcell(this, filter_group)
    class(DistribcellFilter), intent(in) :: this
    integer(HID_T),           intent(in) :: filter_group

    call write_dataset(filter_group, "type", "distribcell")
    call write_dataset(filter_group, "n_bins", this % n_bins)
    call write_dataset(filter_group, "bins", cells(this % cell) % id())
  end subroutine to_statepoint_distribcell

  subroutine initialize_distribcell(this)
    class(DistribcellFilter), intent(inout) :: this

    integer :: id
    integer :: val

    ! Convert id to index.
    id = this % cell
    val = cell_dict % get(id)
    if (val /= EMPTY) then
      this % cell = val
      this % n_bins = cells(this % cell) % n_instances()
    else
      call fatal_error("Could not find cell " // trim(to_str(id)) &
           &// " specified on tally filter.")
    end if
  end subroutine initialize_distribcell

  function text_label_distribcell(this, bin) result(label)
    class(DistribcellFilter), intent(in) :: this
    integer,                  intent(in) :: bin
    character(MAX_LINE_LEN)              :: label

    label = ''
    call find_offset(this % cell, bin-1, label)
    label = "Distributed Cell " // label
  end function text_label_distribcell

!===============================================================================
! FIND_OFFSET (for distribcell) uses a given map number, a target cell ID, and
! a target offset to build a string which is the path from the base universe to
! the target cell with the given offset
!===============================================================================

  subroutine find_offset(i_cell, target_offset, path)
    integer, intent(in) :: i_cell         ! The target cell index
    integer, intent(in) :: target_offset  ! Target offset
    character(*), intent(inout) :: path   ! Path to offset

    integer :: map                  ! Index in maps vector
    integer :: i                    ! Index over cells

    integer(C_INT) :: path_len
    character(kind=C_CHAR), allocatable, target :: path_c(:)

    interface
      function distribcell_path_len(target_cell, map, target_offset, root_univ)&
           bind(C) result(len)
        import C_INT32_T, C_INT
        integer(C_INT32_T), intent(in), value :: target_cell, map, &
                                                 target_offset, root_univ
        integer(C_INT) :: len
      end function distribcell_path_len

      subroutine distribcell_path(target_cell, map, target_offset, root_univ, &
                                  path) bind(C)
        import C_INT32_T, C_CHAR
        integer(C_INT32_T), intent(in), value :: target_cell, map, &
                                                 target_offset, root_univ
        character(kind=C_CHAR), intent(out)   :: path(*)
      end subroutine distribcell_path
    end interface

    ! Get the distribcell index for this cell
    map = cells(i_cell) % distribcell_index()

    path_len = distribcell_path_len(i_cell-1, map-1, target_offset, &
                                    root_universe)
    allocate(path_c(path_len))
    call distribcell_path(i_cell-1, map-1, target_offset, root_universe, &
                          path_c)
    do i = 1, min(path_len, MAX_LINE_LEN)
      if (path_c(i) == C_NULL_CHAR) exit
      path(i:i) = path_c(i)
    end do
  end subroutine find_offset

end module tally_filter_distribcell
