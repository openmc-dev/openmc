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
  public :: find_offset

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

      distribcell_index = cells(this % cell) % distribcell_index
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
                 [p % coord(i + 1) % lattice_x - 1, &
                 p % coord(i + 1) % lattice_y - 1, &
                 p % coord(i + 1) % lattice_z - 1])
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

    integer                 :: offset
    type(Universe), pointer :: univ

    univ => universes(root_universe)
    offset = 0
    label = ''
    call find_offset(this % cell, univ, bin-1, offset, label)
    label = "Distributed Cell " // label
  end function text_label_distribcell

!===============================================================================
! FIND_OFFSET (for distribcell) uses a given map number, a target cell ID, and
! a target offset to build a string which is the path from the base universe to
! the target cell with the given offset
!===============================================================================

  recursive subroutine find_offset(i_cell, univ, target_offset, offset, path)

    integer, intent(in) :: i_cell         ! The target cell index
    type(Universe), intent(in) :: univ  ! Universe to begin search
    integer, intent(in) :: target_offset        ! Target offset
    integer, intent(inout) :: offset    ! Current offset
    character(*), intent(inout) :: path ! Path to offset

    integer :: map                  ! Index in maps vector
    integer :: i                    ! Index over cells
    integer :: k, l, m              ! Indices in lattice
    integer :: n_x, n_y, n_z        ! Lattice cell array dimensions
    integer :: temp_offset          ! Looped sum of offsets
    integer :: i_univ               ! index in universes array
    type(Cell), pointer :: c           ! Pointer to current cell
    type(Universe), pointer :: next_univ  ! Next universe to loop through
    class(Lattice), pointer :: lat        ! Pointer to current lattice

    ! Get the distribcell index for this cell
    map = cells(i_cell) % distribcell_index

    ! Write to the geometry stack
    i_univ = universe_dict % get(univ % id)
    if (i_univ == root_universe) then
      path = trim(path) // "u" // to_str(univ%id)
    else
      path = trim(path) // "->u" // to_str(univ%id)
    end if

    ! Look through all cells in this universe
    do i = 1, size(univ % cells)
      ! If the cell matches the goal and the offset matches final, write to the
      ! geometry stack
      if (univ % cells(i) == i_cell .and. offset == target_offset) then
        c => cells(univ % cells(i))
        path = trim(path) // "->c" // to_str(c % id())
        return
      end if
    end do

    ! Find the fill cell or lattice cell that contains the target offset.
    do i = size(univ % cells), 2, -1
      c => cells(univ % cells(i))

      ! Material cells don't contain other cells so ignore them.
      if (c % type() /= FILL_MATERIAL) then

        if (c % type() == FILL_UNIVERSE) then
          temp_offset = offset + c % offset(map-1)
        else
          ! Get the offset of the first lattice location
          lat => lattices(c % fill() + 1) % obj
          temp_offset = offset + lat % offset(map-1, [0, 0, 0])
        end if

        ! The desired cell is the first cell that gives an offset smaller or
        ! equal to the target offset.
        if (temp_offset <= target_offset) exit
      end if
    end do

    ! Add the cell to the path string.
    c => cells(univ % cells(i))
    path = trim(path) // "->c" // to_str(c%id())

    if (c % type() == FILL_UNIVERSE) then
      ! Recurse into the fill cell.
      offset = c % offset(map-1) + offset
      next_univ => universes(c % fill() + 1)
      call find_offset(i_cell, next_univ, target_offset, offset, path)

    elseif (c % type() == FILL_LATTICE) then
      ! Recurse into the lattice cell.
      lat => lattices(c % fill() + 1) % obj
      select type (lat)

      type is (RectLattice)
        path = trim(path) // "->l" // to_str(lat%id())

        n_x = lat % n_cells(1)
        n_y = lat % n_cells(2)
        n_z = lat % n_cells(3)

        do m = n_z, 1, -1
          do l = n_y, 1, -1
            do k = n_x, 1, -1
              temp_offset = offset + lat % offset(map-1, [k-1, l-1, m-1])
              if (temp_offset <= target_offset) then
                next_univ => universes(lat % get([k-1, l-1, m-1])+1)
                if (lat % is_3d) then
                  path = trim(path) // "(" // trim(to_str(k-1)) // &
                       "," // trim(to_str(l-1)) // "," // &
                       trim(to_str(m-1)) // ")"
                else
                  path = trim(path) // "(" // trim(to_str(k-1)) // &
                       "," // trim(to_str(l-1)) // ")"
                end if
                offset = temp_offset
                call find_offset(i_cell, next_univ, target_offset, offset, path)
                return
              end if
            end do
          end do
        end do

      type is (HexLattice)
        path = trim(path) // "->l" // to_str(lat%id())

        n_z = lat % n_axial
        n_y = 2 * lat % n_rings - 1
        n_x = 2 * lat % n_rings - 1

        do m = n_z, 1, -1
          do l = n_y, 1, -1
            do k = n_x, 1, -1
              if (k + l < lat % n_rings + 1) then
                cycle
              else if (k + l > 3*lat % n_rings - 1) then
                cycle
              end if
              temp_offset = offset + lat % offset(map-1, [k-1, l-1, m-1])
              if (temp_offset <= target_offset) then
                next_univ => universes(lat % get([k-1, l-1, m-1])+1)
                if (lat % is_3d) then
                  path = trim(path) // "(" // &
                       trim(to_str(k - lat % n_rings)) // "," // &
                       trim(to_str(l - lat % n_rings)) // "," // &
                       trim(to_str(m - 1)) // ")"
                else
                  path = trim(path) // "(" // &
                       trim(to_str(k - lat % n_rings)) // "," // &
                       trim(to_str(l - lat % n_rings)) // ")"
                end if
                offset = temp_offset
                call find_offset(i_cell, next_univ, target_offset, offset, path)
                return
              end if
            end do
          end do
        end do

      end select

    end if
  end subroutine find_offset

end module tally_filter_distribcell
