module tally_filter_distribcell

  use, intrinsic :: ISO_C_BINDING

  use hdf5, only: HID_T

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
        if (cells(p % coord(i) % cell) % type == FILL_UNIVERSE) then
          offset = offset + cells(p % coord(i) % cell) % &
               offset(distribcell_index)
        elseif (cells(p % coord(i) % cell) % type == FILL_LATTICE) then
          if (lattices(p % coord(i + 1) % lattice) % obj &
               % are_valid_indices([&
               p % coord(i + 1) % lattice_x, &
               p % coord(i + 1) % lattice_y, &
               p % coord(i + 1) % lattice_z])) then
            offset = offset + lattices(p % coord(i + 1) % lattice) % obj % &
                 offset(distribcell_index, &
                 p % coord(i + 1) % lattice_x, &
                 p % coord(i + 1) % lattice_y, &
                 p % coord(i + 1) % lattice_z)
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
    call write_dataset(filter_group, "bins", cells(this % cell) % id)
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
      this % n_bins = cells(this % cell) % instances
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
    integer :: i, j                 ! Index over cells
    integer :: k, l, m              ! Indices in lattice
    integer :: old_k, old_l, old_m  ! Previous indices in lattice
    integer :: n_x, n_y, n_z        ! Lattice cell array dimensions
    integer :: n                    ! Number of cells to search
    integer :: cell_index           ! Index in cells array
    integer :: lat_offset           ! Offset from lattice
    integer :: temp_offset          ! Looped sum of offsets
    integer :: i_univ               ! index in universes array
    logical :: this_cell = .false.  ! Advance in this cell?
    logical :: later_cell = .false. ! Fill cells after this one?
    type(Cell), pointer :: c           ! Pointer to current cell
    type(Universe), pointer :: next_univ  ! Next universe to loop through
    class(Lattice), pointer :: lat        ! Pointer to current lattice

    ! Get the distribcell index for this cell
    map = cells(i_cell) % distribcell_index

    n = size(univ % cells)

    ! Write to the geometry stack
    i_univ = universe_dict % get(univ % id)
    if (i_univ == root_universe) then
      path = trim(path) // "u" // to_str(univ%id)
    else
      path = trim(path) // "->u" // to_str(univ%id)
    end if

    ! Look through all cells in this universe
    do i = 1, n
      ! If the cell matches the goal and the offset matches final, write to the
      ! geometry stack
      if (univ % cells(i) == i_cell .and. offset == target_offset) then
        c => cells(univ % cells(i))
        path = trim(path) // "->c" // to_str(c % id)
        return
      end if
    end do

    ! Find the fill cell or lattice cell that we need to enter
    do i = 1, n

      later_cell = .false.

      c => cells(univ % cells(i))

      this_cell = .false.

      ! If we got here, we still think the target is in this universe
      ! or further down, but it's not this exact cell.
      ! Compare offset to next cell to see if we should enter this cell
      if (i /= n) then

        do j = i+1, n

          c => cells(univ % cells(j))

          ! Skip normal cells which do not have offsets
          if (c % type == FILL_MATERIAL) cycle

          ! Break loop once we've found the next cell with an offset
          exit
        end do

        ! Ensure we didn't just end the loop by iteration
        if (c % type /= FILL_MATERIAL) then

          ! There are more cells in this universe that it could be in
          later_cell = .true.

          ! Two cases, lattice or fill cell
          if (c % type == FILL_UNIVERSE) then
            temp_offset = c % offset(map)

          ! Get the offset of the first lattice location
          else
            lat => lattices(c % fill) % obj
            temp_offset = lat % offset(map, 1, 1, 1)
          end if

          ! If the final offset is in the range of offset - temp_offset+offset
          ! then the goal is in this cell
          if (target_offset < temp_offset + offset) then
            this_cell = .true.
          end if
        end if
      end if

      if (n == 1 .and. c % type /= FILL_MATERIAL) then
        this_cell = .true.
      end if

      if (.not. later_cell) then
        this_cell = .true.
      end if

      ! Get pointer to THIS cell because target must be in this cell
      if (this_cell) then

        cell_index = univ % cells(i)
        c => cells(cell_index)

        path = trim(path) // "->c" // to_str(c%id)

        ! ====================================================================
        ! CELL CONTAINS LOWER UNIVERSE, RECURSIVELY FIND CELL
        if (c % type == FILL_UNIVERSE) then

          ! Enter this cell to update the current offset
          offset = c % offset(map) + offset

          next_univ => universes(c % fill)
          call find_offset(i_cell, next_univ, target_offset, offset, path)
          return

        ! ====================================================================
        ! CELL CONTAINS LATTICE, RECURSIVELY FIND CELL
        elseif (c % type == FILL_LATTICE) then

          ! Set current lattice
          lat => lattices(c % fill) % obj

          select type (lat)

          ! ==================================================================
          ! RECTANGULAR LATTICES
          type is (RectLattice)

            ! Write to the geometry stack
            path = trim(path) // "->l" // to_str(lat%id)

            n_x = lat % n_cells(1)
            n_y = lat % n_cells(2)
            n_z = lat % n_cells(3)
            old_m = 1
            old_l = 1
            old_k = 1

            ! Loop over lattice coordinates
            do k = 1, n_x
              do l = 1, n_y
                do m = 1, n_z

                  if (target_offset >= lat % offset(map, k, l, m) + offset) then
                    if (k == n_x .and. l == n_y .and. m == n_z) then
                      ! This is last lattice cell, so target must be here
                      lat_offset = lat % offset(map, k, l, m)
                      offset = offset + lat_offset
                      next_univ => universes(lat % universes(k, l, m))
                      if (lat % is_3d) then
                        path = trim(path) // "(" // trim(to_str(k-1)) // &
                             "," // trim(to_str(l-1)) // "," // &
                             trim(to_str(m-1)) // ")"
                      else
                        path = trim(path) // "(" // trim(to_str(k-1)) // &
                             "," // trim(to_str(l-1)) // ")"
                      end if
                      call find_offset(i_cell, next_univ, target_offset, offset, path)
                      return
                    else
                      old_m = m
                      old_l = l
                      old_k = k
                      cycle
                    end if
                  else
                    ! Target is at this lattice position
                    lat_offset = lat % offset(map, old_k, old_l, old_m)
                    offset = offset + lat_offset
                    next_univ => universes(lat % universes(old_k, old_l, old_m))
                    if (lat % is_3d) then
                      path = trim(path) // "(" // trim(to_str(old_k-1)) // &
                           "," // trim(to_str(old_l-1)) // "," // &
                           trim(to_str(old_m-1)) // ")"
                    else
                      path = trim(path) // "(" // trim(to_str(old_k-1)) // &
                           "," // trim(to_str(old_l-1)) // ")"
                    end if
                    call find_offset(i_cell, next_univ, target_offset, offset, path)
                    return
                  end if

                end do
              end do
            end do

          ! ==================================================================
          ! HEXAGONAL LATTICES
          type is (HexLattice)

            ! Write to the geometry stack
            path = trim(path) // "->l" // to_str(lat%id)

            n_z = lat % n_axial
            n_y = 2 * lat % n_rings - 1
            n_x = 2 * lat % n_rings - 1
            old_m = 1
            old_l = 1
            old_k = 1

            ! Loop over lattice coordinates
            do m = 1, n_z
              do l = 1, n_y
                do k = 1, n_x

                  ! This array position is never used
                  if (k + l < lat % n_rings + 1) then
                    cycle
                  ! This array position is never used
                  else if (k + l > 3*lat % n_rings - 1) then
                    cycle
                  end if

                  if (target_offset >= lat % offset(map, k, l, m) + offset) then
                    if (k == lat % n_rings .and. l == n_y .and. m == n_z) then
                      ! This is last lattice cell, so target must be here
                      lat_offset = lat % offset(map, k, l, m)
                      offset = offset + lat_offset
                      next_univ => universes(lat % universes(k, l, m))
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
                      call find_offset(i_cell, next_univ, target_offset, offset, path)
                      return
                    else
                      old_m = m
                      old_l = l
                      old_k = k
                      cycle
                    end if
                  else
                    ! Target is at this lattice position
                    lat_offset = lat % offset(map, old_k, old_l, old_m)
                    offset = offset + lat_offset
                    next_univ => universes(lat % universes(old_k, old_l, old_m))
                    if (lat % is_3d) then
                      path = trim(path) // "(" // &
                           trim(to_str(old_k - lat % n_rings)) // "," // &
                           trim(to_str(old_l - lat % n_rings)) // "," // &
                           trim(to_str(old_m - 1)) // ")"
                    else
                      path = trim(path) // "(" // &
                           trim(to_str(old_k - lat % n_rings)) // "," // &
                           trim(to_str(old_l - lat % n_rings)) // ")"
                    end if
                    call find_offset(i_cell, next_univ, target_offset, offset, path)
                    return
                  end if

                end do
              end do
            end do

          end select

        end if
      end if
    end do
  end subroutine find_offset

end module tally_filter_distribcell
