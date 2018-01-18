module geometry

  use constants
  use error,                  only: fatal_error, warning, write_message
  use geometry_header
  use particle_header,        only: LocalCoord, Particle
  use simulation_header
  use settings
  use surface_header
  use stl_vector,             only: VectorInt
  use string,                 only: to_str

  implicit none

contains

!===============================================================================
! CELL_CONTAINS determines if a cell contains the particle at a given
! location. The bounds of the cell are detemined by a logical expression
! involving surface half-spaces. At initialization, the expression was converted
! to RPN notation.
!
! The function is split into two cases, one for simple cells (those involving
! only the intersection of half-spaces) and one for complex cells. Simple cells
! can be evaluated with short circuit evaluation, i.e., as soon as we know that
! one half-space is not satisfied, we can exit. This provides a performance
! benefit for the common case. In complex_cell_contains, we evaluate the RPN
! expression using a stack, similar to how a RPN calculator would work.
!===============================================================================

  pure function cell_contains(c, p) result(in_cell)
    type(Cell), intent(in) :: c
    type(Particle), intent(in) :: p
    logical :: in_cell

    if (c%simple) then
      in_cell = simple_cell_contains(c, p)
    else
      in_cell = complex_cell_contains(c, p)
    end if
  end function cell_contains

  pure function simple_cell_contains(c, p) result(in_cell)
    type(Cell), intent(in) :: c
    type(Particle), intent(in) :: p
    logical :: in_cell

    integer :: i
    integer :: token
    logical :: actual_sense    ! sense of particle wrt surface

    in_cell = .true.
    do i = 1, size(c%rpn)
      token = c%rpn(i)
      if (token < OP_UNION) then
        ! If the token is not an operator, evaluate the sense of particle with
        ! respect to the surface and see if the token matches the sense. If the
        ! particle's surface attribute is set and matches the token, that
        ! overrides the determination based on sense().
        if (token == p%surface) then
          cycle
        elseif (-token == p%surface) then
          in_cell = .false.
          exit
        else
          actual_sense = surfaces(abs(token))%obj%sense(&
               p%coord(p%n_coord)%xyz, p%coord(p%n_coord)%uvw)
          if (actual_sense .neqv. (token > 0)) then
            in_cell = .false.
            exit
          end if
        end if
      end if
    end do
  end function simple_cell_contains

  pure function complex_cell_contains(c, p) result(in_cell)
    type(Cell), intent(in) :: c
    type(Particle), intent(in) :: p
    logical :: in_cell

    integer :: i
    integer :: token
    integer :: i_stack
    logical :: actual_sense    ! sense of particle wrt surface
    logical :: stack(size(c%rpn))

    i_stack = 0
    do i = 1, size(c%rpn)
      token = c%rpn(i)

      ! If the token is a binary operator (intersection/union), apply it to
      ! the last two items on the stack. If the token is a unary operator
      ! (complement), apply it to the last item on the stack.
      select case (token)
      case (OP_UNION)
        stack(i_stack - 1) = stack(i_stack - 1) .or. stack(i_stack)
        i_stack = i_stack - 1
      case (OP_INTERSECTION)
        stack(i_stack - 1) = stack(i_stack - 1) .and. stack(i_stack)
        i_stack = i_stack - 1
      case (OP_COMPLEMENT)
        stack(i_stack) = .not. stack(i_stack)
      case default
        ! If the token is not an operator, evaluate the sense of particle with
        ! respect to the surface and see if the token matches the sense. If the
        ! particle's surface attribute is set and matches the token, that
        ! overrides the determination based on sense().
        i_stack = i_stack + 1
        if (token == p%surface) then
          stack(i_stack) = .true.
        elseif (-token == p%surface) then
          stack(i_stack) = .false.
        else
          actual_sense = surfaces(abs(token))%obj%sense(&
               p%coord(p%n_coord)%xyz, p%coord(p%n_coord)%uvw)
          stack(i_stack) = (actual_sense .eqv. (token > 0))
        end if
      end select

    end do

    if (i_stack == 1) then
      ! The one remaining logical on the stack indicates whether the particle is
      ! in the cell.
      in_cell = stack(i_stack)
    else
      ! This case occurs if there is no region specification since i_stack will
      ! still be zero.
      in_cell = .true.
    end if
  end function complex_cell_contains

!===============================================================================
! CHECK_CELL_OVERLAP checks for overlapping cells at the current particle's
! position using cell_contains and the LocalCoord's built up by find_cell
!===============================================================================

  subroutine check_cell_overlap(p)

    type(Particle), intent(inout) :: p

    integer :: i                       ! cell loop index on a level
    integer :: j                       ! coordinate level index
    integer :: n_coord                 ! saved number of coordinate levels
    integer :: n                       ! number of cells to search on a level
    integer :: index_cell              ! index in cells array
    type(Cell),       pointer :: c     ! pointer to cell
    type(Universe),   pointer :: univ  ! universe to search in

    ! loop through each coordinate level
    n_coord = p % n_coord
    do j = 1, n_coord
      p % n_coord = j
      univ => universes(p % coord(j) % universe)
      n = size(univ % cells)

      ! loop through each cell on this level
      do i = 1, n
        index_cell = univ % cells(i)
        c => cells(index_cell)

        if (cell_contains(c, p)) then
          ! the particle should only be contained in one cell per level
          if (index_cell /= p % coord(j) % cell) then
            call fatal_error("Overlapping cells detected: " &
                 &// trim(to_str(cells(index_cell) % id)) // ", " &
                 &// trim(to_str(cells(p % coord(j) % cell) % id)) &
                 &// " on universe " // trim(to_str(univ % id)))
          end if

          overlap_check_cnt(index_cell) = overlap_check_cnt(index_cell) + 1

        end if

      end do
    end do

  end subroutine check_cell_overlap

!===============================================================================
! FIND_CELL determines what cell a source particle is in within a particular
! universe. If the base universe is passed, the particle should be found as long
! as it's within the geometry
!===============================================================================

  recursive subroutine find_cell(p, found, search_cells)

    type(Particle), intent(inout) :: p
    logical,        intent(inout) :: found
    integer,        optional      :: search_cells(:)
    integer :: i                    ! index over cells
    integer :: j, k                 ! coordinate level index
    integer :: offset               ! instance # of a distributed cell
    integer :: distribcell_index
    integer :: i_xyz(3)             ! indices in lattice
    integer :: n                    ! number of cells to search
    integer :: i_cell               ! index in cells array
    integer :: i_universe           ! index in universes array
    logical :: use_search_cells     ! use cells provided as argument

    do j = p % n_coord + 1, MAX_COORD
      call p % coord(j) % reset()
    end do
    j = p % n_coord

    ! set size of list to search
    i_universe = p % coord(j) % universe
    if (present(search_cells)) then
      use_search_cells = .true.
      n = size(search_cells)
    else
      use_search_cells = .false.
      n = size(universes(i_universe) % cells)
    end if

    found = .false.
    CELL_LOOP: do i = 1, n
      ! select cells based on whether we are searching a universe or a provided
      ! list of cells (this would be for lists of neighbor cells)
      if (use_search_cells) then
        i_cell = search_cells(i)
        ! check to make sure search cell is in same universe
        if (cells(i_cell) % universe /= i_universe) cycle
      else
        i_cell = universes(i_universe) % cells(i)
      end if

      ! Move on to the next cell if the particle is not inside this cell
      if (cell_contains(cells(i_cell), p)) then
        ! Set cell on this level
        p % coord(j) % cell = i_cell

        ! Show cell information on trace
        if (verbosity >= 10 .or. trace) then
          call write_message("    Entering cell " // trim(to_str(&
               cells(i_cell) % id)))
        end if

        found = .true.
        exit
      end if
    end do CELL_LOOP

    if (found) then
      associate(c => cells(i_cell))
        CELL_TYPE: if (c % type == FILL_MATERIAL) then
          ! ======================================================================
          ! AT LOWEST UNIVERSE, TERMINATE SEARCH

          ! Save previous material and temperature
          p % last_material = p % material
          p % last_sqrtkT = p % sqrtkT

          ! Get distributed offset
          if (size(c % material) > 1 .or. size(c % sqrtkT) > 1) then
            ! Distributed instances of this cell have different
            ! materials/temperatures. Determine which instance this is for
            ! assigning the matching material/temperature.
            distribcell_index = c % distribcell_index
            offset = 0
            do k = 1, p % n_coord
              if (cells(p % coord(k) % cell) % type == FILL_UNIVERSE) then
                offset = offset + cells(p % coord(k) % cell) % &
                     offset(distribcell_index)
              elseif (cells(p % coord(k) % cell) % type == FILL_LATTICE) then
                if (lattices(p % coord(k + 1) % lattice) % obj &
                     % are_valid_indices([&
                     p % coord(k + 1) % lattice_x, &
                     p % coord(k + 1) % lattice_y, &
                     p % coord(k + 1) % lattice_z])) then
                  offset = offset + lattices(p % coord(k + 1) % lattice) % obj % &
                       offset(distribcell_index, &
                       p % coord(k + 1) % lattice_x, &
                       p % coord(k + 1) % lattice_y, &
                       p % coord(k + 1) % lattice_z)
                end if
              end if
            end do

            ! Keep track of which instance of the cell the particle is in
            p % cell_instance = offset + 1
          else
            p % cell_instance = 1
          end if

          ! Save the material
          if (size(c % material) > 1) then
            p % material = c % material(offset + 1)
          else
            p % material = c % material(1)
          end if

          ! Save the temperature
          if (size(c % sqrtkT) > 1) then
            p % sqrtkT = c % sqrtkT(offset + 1)
          else
            p % sqrtkT = c % sqrtkT(1)
          end if

        elseif (c % type == FILL_UNIVERSE) then CELL_TYPE
          ! ======================================================================
          ! CELL CONTAINS LOWER UNIVERSE, RECURSIVELY FIND CELL

          ! Store lower level coordinates
          p % coord(j + 1) % xyz = p % coord(j) % xyz
          p % coord(j + 1) % uvw = p % coord(j) % uvw

          ! Move particle to next level and set universe
          j = j + 1
          p % n_coord = j
          p % coord(j) % universe = c % fill

          ! Apply translation
          if (allocated(c % translation)) then
            p % coord(j) % xyz = p % coord(j) % xyz - c % translation
          end if

          ! Apply rotation
          if (allocated(c % rotation_matrix)) then
            p % coord(j) % xyz = matmul(c % rotation_matrix, p % coord(j) % xyz)
            p % coord(j) % uvw = matmul(c % rotation_matrix, p % coord(j) % uvw)
            p % coord(j) % rotated = .true.
          end if

          call find_cell(p, found)
          j = p % n_coord

        elseif (c % type == FILL_LATTICE) then CELL_TYPE
          ! ======================================================================
          ! CELL CONTAINS LATTICE, RECURSIVELY FIND CELL

          associate (lat => lattices(c % fill) % obj)
            ! Determine lattice indices
            i_xyz = lat % get_indices(p % coord(j) % xyz + TINY_BIT * p % coord(j) % uvw)

            ! Store lower level coordinates
            p % coord(j + 1) % xyz = lat % get_local_xyz(p % coord(j) % xyz, i_xyz)
            p % coord(j + 1) % uvw = p % coord(j) % uvw

            ! set particle lattice indices
            p % coord(j + 1) % lattice   = c % fill
            p % coord(j + 1) % lattice_x = i_xyz(1)
            p % coord(j + 1) % lattice_y = i_xyz(2)
            p % coord(j + 1) % lattice_z = i_xyz(3)

            ! Set the next lowest coordinate level.
            if (lat % are_valid_indices(i_xyz)) then
              ! Particle is inside the lattice.
              p % coord(j + 1) % universe = &
                   lat % universes(i_xyz(1), i_xyz(2), i_xyz(3))

            else
              ! Particle is outside the lattice.
              if (lat % outer == NO_OUTER_UNIVERSE) then
                call p % mark_as_lost("Particle " // trim(to_str(p %id)) &
                     // " is outside lattice " // trim(to_str(lat % id)) &
                     // " but the lattice has no defined outer universe.")
                return
              else
                p % coord(j + 1) % universe = lat % outer
              end if
            end if
          end associate

          ! Move particle to next level and search for the lower cells.
          j = j + 1
          p % n_coord = j

          call find_cell(p, found)
          j = p % n_coord

        end if CELL_TYPE
      end associate
    end if

  end subroutine find_cell

!===============================================================================
! CROSS_LATTICE moves a particle into a new lattice element
!===============================================================================

  subroutine cross_lattice(p, lattice_translation)

    type(Particle), intent(inout) :: p
    integer,        intent(in)    :: lattice_translation(3)
    integer :: j
    integer :: i_xyz(3)       ! indices in lattice
    logical :: found          ! particle found in cell?
    class(Lattice),   pointer :: lat

    j = p % n_coord
    lat => lattices(p % coord(j) % lattice) % obj

    if (verbosity >= 10 .or. trace) then
      call write_message("    Crossing lattice " // trim(to_str(lat % id)) &
           &// ". Current position (" // trim(to_str(p % coord(j) % lattice_x)) &
           &// "," // trim(to_str(p % coord(j) % lattice_y)) // "," &
           &// trim(to_str(p % coord(j) % lattice_z)) // ")")
    end if

    ! Set the lattice indices.
    p % coord(j) % lattice_x = p % coord(j) % lattice_x + lattice_translation(1)
    p % coord(j) % lattice_y = p % coord(j) % lattice_y + lattice_translation(2)
    p % coord(j) % lattice_z = p % coord(j) % lattice_z + lattice_translation(3)
    i_xyz(1) = p % coord(j) % lattice_x
    i_xyz(2) = p % coord(j) % lattice_y
    i_xyz(3) = p % coord(j) % lattice_z

    ! Set the new coordinate position.
    p % coord(j) % xyz = lat % get_local_xyz(p % coord(j - 1) % xyz, i_xyz)

    OUTSIDE_LAT: if (.not. lat % are_valid_indices(i_xyz)) then
      ! The particle is outside the lattice.  Search for it from base coord
      p % n_coord = 1
      call find_cell(p, found)
      if (.not. found) then
        if (p % alive) then ! Particle may have been killed in find_cell
          call p % mark_as_lost("Could not locate particle " &
               // trim(to_str(p % id)) // " after crossing a lattice boundary.")
          return
        end if
      end if

    else OUTSIDE_LAT

      ! Find cell in next lattice element
      p % coord(j) % universe = lat % universes(i_xyz(1), i_xyz(2), i_xyz(3))

      call find_cell(p, found)
      if (.not. found) then
        ! In some circumstances, a particle crossing the corner of a cell may
        ! not be able to be found in the next universe. In this scenario we cut
        ! off all lower-level coordinates and search from universe zero

        ! Remove lower coordinates
        p % n_coord = 1

        ! Search for particle
        call find_cell(p, found)
        if (.not. found) then
          call p % mark_as_lost("Could not locate particle " // &
               trim(to_str(p % id)) // " after crossing a lattice boundary.")
          return
        end if
      end if
    end if OUTSIDE_LAT

  end subroutine cross_lattice

!===============================================================================
! DISTANCE_TO_BOUNDARY calculates the distance to the nearest boundary for a
! particle 'p' traveling in a certain direction. For a cell in a subuniverse
! that has a parent cell, also include the surfaces of the edge of the universe.
!===============================================================================

  subroutine distance_to_boundary(p, dist, surface_crossed, lattice_translation, &
       next_level)
    type(Particle), intent(inout) :: p
    real(8),        intent(out)   :: dist
    integer,        intent(out)   :: surface_crossed
    integer,        intent(out)   :: lattice_translation(3)
    integer,        intent(out)   :: next_level

    integer :: i                  ! index for surface in cell
    integer :: j
    integer :: index_surf         ! index in surfaces array (with sign)
    integer :: i_xyz(3)           ! lattice indices
    integer :: level_surf_cross   ! surface crossed on current level
    integer :: level_lat_trans(3) ! lattice translation on current level
    real(8) :: x,y,z              ! particle coordinates
    real(8) :: xyz_t(3)           ! local particle coordinates
    real(8) :: beta, gama         ! skewed particle coordiantes
    real(8) :: u,v,w              ! particle directions
    real(8) :: beta_dir           ! skewed particle direction
    real(8) :: gama_dir           ! skewed particle direction
    real(8) :: edge               ! distance to oncoming edge
    real(8) :: d                  ! evaluated distance
    real(8) :: d_lat              ! distance to lattice boundary
    real(8) :: d_surf             ! distance to surface
    real(8) :: x0,y0,z0           ! coefficients for surface
    real(8) :: xyz_cross(3)       ! coordinates at projected surface crossing
    logical :: coincident         ! is particle on surface?
    type(Cell),       pointer :: c
    class(Surface),   pointer :: surf
    class(Lattice),   pointer :: lat

    ! inialize distance to infinity (huge)
    dist = INFINITY
    d_lat = INFINITY
    d_surf = INFINITY
    lattice_translation(:) = [0, 0, 0]

    next_level = 0

    ! Loop over each universe level
    LEVEL_LOOP: do j = 1, p % n_coord

      ! get pointer to cell on this level
      c => cells(p % coord(j) % cell)

      ! copy directional cosines
      u = p % coord(j) % uvw(1)
      v = p % coord(j) % uvw(2)
      w = p % coord(j) % uvw(3)

      ! =======================================================================
      ! FIND MINIMUM DISTANCE TO SURFACE IN THIS CELL

      SURFACE_LOOP: do i = 1, size(c % region)
        index_surf = c % region(i)
        coincident = (index_surf == p % surface)

        ! ignore this token if it corresponds to an operator rather than a
        ! region.
        index_surf = abs(index_surf)
        if (index_surf >= OP_UNION) cycle

        ! Calculate distance to surface
        surf => surfaces(index_surf) % obj
        d = surf % distance(p % coord(j) % xyz, p % coord(j) % uvw, coincident)

        ! Check if calculated distance is new minimum
        if (d < d_surf) then
          if (abs(d - d_surf)/d_surf >= FP_PRECISION) then
            d_surf = d
            level_surf_cross = -c % region(i)
          end if
        end if
      end do SURFACE_LOOP

      ! =======================================================================
      ! FIND MINIMUM DISTANCE TO LATTICE SURFACES

      LAT_COORD: if (p % coord(j) % lattice /= NONE) then
        lat => lattices(p % coord(j) % lattice) % obj

        LAT_TYPE: select type(lat)

        type is (RectLattice)
          ! copy local coordinates
          x = p % coord(j) % xyz(1)
          y = p % coord(j) % xyz(2)
          z = p % coord(j) % xyz(3)

          ! determine oncoming edge
          x0 = sign(lat % pitch(1) * HALF, u)
          y0 = sign(lat % pitch(2) * HALF, v)

          ! left and right sides
          if (abs(x - x0) < FP_PRECISION) then
            d = INFINITY
          elseif (u == ZERO) then
            d = INFINITY
          else
            d = (x0 - x)/u
          end if

          d_lat = d
          if (u > 0) then
            level_lat_trans(:) = [1, 0, 0]
          else
            level_lat_trans(:) = [-1, 0, 0]
          end if

          ! front and back sides
          if (abs(y - y0) < FP_PRECISION) then
            d = INFINITY
          elseif (v == ZERO) then
            d = INFINITY
          else
            d = (y0 - y)/v
          end if

          if (d < d_lat) then
            d_lat = d
            if (v > 0) then
              level_lat_trans(:) = [0, 1, 0]
            else
              level_lat_trans(:) = [0, -1, 0]
            end if
          end if

          if (lat % is_3d) then
            z0 = sign(lat % pitch(3) * HALF, w)

            ! top and bottom sides
            if (abs(z - z0) < FP_PRECISION) then
              d = INFINITY
            elseif (w == ZERO) then
              d = INFINITY
            else
              d = (z0 - z)/w
            end if

            if (d < d_lat) then
              d_lat = d
              if (w > 0) then
                level_lat_trans(:) = [0, 0, 1]
              else
                level_lat_trans(:) = [0, 0, -1]
              end if
            end if
          end if

        type is (HexLattice) LAT_TYPE
          ! Copy local coordinates.
          z = p % coord(j) % xyz(3)
          i_xyz(1) = p % coord(j) % lattice_x
          i_xyz(2) = p % coord(j) % lattice_y
          i_xyz(3) = p % coord(j) % lattice_z

          ! Compute velocities along the hexagonal axes.
          beta_dir = u*sqrt(THREE)/TWO + v/TWO
          gama_dir = u*sqrt(THREE)/TWO - v/TWO

          ! Note that hexagonal lattice distance calculations are performed
          ! using the particle's coordinates relative to the neighbor lattice
          ! cells, not relative to the particle's current cell.  This is done
          ! because there is significant disagreement between neighboring cells
          ! on where the lattice boundary is due to the worse finite precision
          ! of hex lattices.

          ! Upper right and lower left sides.
          edge = -sign(lat % pitch(1)/TWO, beta_dir)  ! Oncoming edge
          if (beta_dir > ZERO) then
            xyz_t = lat % get_local_xyz(p % coord(j - 1) % xyz, i_xyz+[1, 0, 0])
          else
            xyz_t = lat % get_local_xyz(p % coord(j - 1) % xyz, i_xyz+[-1, 0, 0])
          end if
          beta = xyz_t(1)*sqrt(THREE)/TWO + xyz_t(2)/TWO
          if (abs(beta - edge) < FP_PRECISION) then
            d = INFINITY
          else if (beta_dir == ZERO) then
            d = INFINITY
          else
            d = (edge - beta)/beta_dir
          end if

          d_lat = d
          if (beta_dir > 0) then
            level_lat_trans(:) = [1, 0, 0]
          else
            level_lat_trans(:) = [-1, 0, 0]
          end if

          ! Lower right and upper left sides.
          edge = -sign(lat % pitch(1)/TWO, gama_dir)  ! Oncoming edge
          if (gama_dir > ZERO) then
            xyz_t = lat % get_local_xyz(p % coord(j - 1) % xyz, i_xyz+[1, -1, 0])
          else
            xyz_t = lat % get_local_xyz(p % coord(j - 1) % xyz, i_xyz+[-1, 1, 0])
          end if
          gama = xyz_t(1)*sqrt(THREE)/TWO - xyz_t(2)/TWO
          if (abs(gama - edge) < FP_PRECISION) then
            d = INFINITY
          else if (gama_dir == ZERO) then
            d = INFINITY
          else
            d = (edge - gama)/gama_dir
          end if

          if (d < d_lat) then
            d_lat = d
            if (gama_dir > 0) then
              level_lat_trans(:) = [1, -1, 0]
            else
              level_lat_trans(:) = [-1, 1, 0]
            end if
          end if

          ! Upper and lower sides.
          edge = -sign(lat % pitch(1)/TWO, v)  ! Oncoming edge
          if (v > ZERO) then
            xyz_t = lat % get_local_xyz(p % coord(j - 1) % xyz, i_xyz+[0, 1, 0])
          else
            xyz_t = lat % get_local_xyz(p % coord(j - 1) % xyz, i_xyz+[0, -1, 0])
          end if
          if (abs(xyz_t(2) - edge) < FP_PRECISION) then
            d = INFINITY
          else if (v == ZERO) then
            d = INFINITY
          else
            d = (edge - xyz_t(2))/v
          end if

          if (d < d_lat) then
            d_lat = d
            if (v > 0) then
              level_lat_trans(:) = [0, 1, 0]
            else
              level_lat_trans(:) = [0, -1, 0]
            end if
          end if

          ! Top and bottom sides.
          if (lat % is_3d) then
            z0 = sign(lat % pitch(2) * HALF, w)

            if (abs(z - z0) < FP_PRECISION) then
              d = INFINITY
            elseif (w == ZERO) then
              d = INFINITY
            else
              d = (z0 - z)/w
            end if

            if (d < d_lat) then
              d_lat = d
              if (w > 0) then
                level_lat_trans(:) = [0, 0, 1]
              else
                level_lat_trans(:) = [0, 0, -1]
              end if
            end if
          end if
        end select LAT_TYPE

        if (d_lat < ZERO) then
          call p % mark_as_lost("Particle " // trim(to_str(p % id)) &
               //" had a negative distance to a lattice boundary. d = " &
               //trim(to_str(d_lat)))
        end if
      end if LAT_COORD

      ! If the boundary on this lattice level is coincident with a boundary on
      ! a higher level then we need to make sure that the higher level boundary
      ! is selected.  This logic must include consideration of floating point
      ! precision.
      if (d_surf < d_lat) then
        if ((dist - d_surf)/dist >= FP_REL_PRECISION) then
          dist = d_surf

          ! If the cell is not simple, it is possible that both the negative and
          ! positive half-space were given in the region specification. Thus, we
          ! have to explicitly check which half-space the particle would be
          ! traveling into if the surface is crossed
          if (.not. c % simple) then
            xyz_cross(:) = p % coord(j) % xyz + d_surf*p % coord(j) % uvw
            surf => surfaces(abs(level_surf_cross)) % obj
            if (dot_product(p % coord(j) % uvw, &
                 surf % normal(xyz_cross)) > ZERO) then
              surface_crossed = abs(level_surf_cross)
            else
              surface_crossed = -abs(level_surf_cross)
            end if
          else
            surface_crossed = level_surf_cross
          end if

          lattice_translation(:) = [0, 0, 0]
          next_level = j
        end if
      else
        if ((dist - d_lat)/dist >= FP_REL_PRECISION) then
          dist = d_lat
          surface_crossed = NONE
          lattice_translation(:) = level_lat_trans
          next_level = j
        end if
      end if

    end do LEVEL_LOOP

  end subroutine distance_to_boundary

!===============================================================================
! NEIGHBOR_LISTS builds a list of neighboring cells to each surface to speed up
! searches when a cell boundary is crossed.
!===============================================================================

  subroutine neighbor_lists()

    integer :: i  ! index in cells/surfaces array
    integer :: j  ! index in region specification
    integer :: k  ! surface half-space spec
    integer :: n  ! size of vector
    type(VectorInt), allocatable :: neighbor_pos(:)
    type(VectorInt), allocatable :: neighbor_neg(:)

    call write_message("Building neighboring cells lists for each surface...", &
         6)

    allocate(neighbor_pos(n_surfaces))
    allocate(neighbor_neg(n_surfaces))

    do i = 1, n_cells
      do j = 1, size(cells(i)%region)
        ! Get token from region specification and skip any tokens that
        ! correspond to operators rather than regions
        k = cells(i)%region(j)
        if (abs(k) >= OP_UNION) cycle

        ! Add this cell ID to neighbor list for k-th surface
        if (k > 0) then
          call neighbor_pos(abs(k))%push_back(i)
        else
          call neighbor_neg(abs(k))%push_back(i)
        end if
      end do
    end do

    do i = 1, n_surfaces
      ! Copy positive neighbors to Surface instance
      n = neighbor_pos(i)%size()
      if (n > 0) then
        allocate(surfaces(i)%obj%neighbor_pos(n))
        surfaces(i)%obj%neighbor_pos(:) = neighbor_pos(i)%data(1:n)
      end if

      ! Copy negative neighbors to Surface instance
      n = neighbor_neg(i)%size()
      if (n > 0) then
        allocate(surfaces(i)%obj%neighbor_neg(n))
        surfaces(i)%obj%neighbor_neg(:) = neighbor_neg(i)%data(1:n)
      end if
    end do

  end subroutine neighbor_lists

!===============================================================================
! CALC_OFFSETS calculates and stores the offsets in all fill cells. This
! routine is called once upon initialization.
!===============================================================================

  subroutine calc_offsets(univ_id, map, univ, counts, found)

    integer, intent(in)        :: univ_id         ! target universe ID
    integer, intent(in)        :: map          ! map index in vector of maps
    type(Universe), intent(in) :: univ         ! universe searching in
    integer, intent(inout)     :: counts(:,:)  ! target count
    logical, intent(inout)     :: found(:,:)   ! target found

    integer :: i                          ! index over cells
    integer :: j, k, m                    ! indices in lattice
    integer :: n                          ! number of cells to search
    integer :: offset                     ! total offset for a given cell
    integer :: cell_index                 ! index in cells array
    type(Cell),     pointer :: c          ! pointer to current cell
    type(Universe), pointer :: next_univ  ! next universe to cycle through
    class(Lattice), pointer :: lat        ! pointer to current lattice

    n = size(univ % cells)
    offset = 0

    do i = 1, n

      cell_index = univ % cells(i)

      ! get pointer to cell
      c => cells(cell_index)

      ! ====================================================================
      ! AT LOWEST UNIVERSE, TERMINATE SEARCH
      if (c % type == FILL_MATERIAL) then

      ! ====================================================================
      ! CELL CONTAINS LOWER UNIVERSE, RECURSIVELY FIND CELL
      elseif (c % type == FILL_UNIVERSE) then
        ! Set offset for the cell on this level
        c % offset(map) = offset

        ! Count contents of this cell
        next_univ => universes(c % fill)
        offset = offset + count_target(next_univ, counts, found, univ_id, map)

        ! Move into the next universe
        next_univ => universes(c % fill)
        c => cells(cell_index)

      ! ====================================================================
      ! CELL CONTAINS LATTICE, RECURSIVELY FIND CELL
      elseif (c % type == FILL_LATTICE) then

        ! Set current lattice
        lat => lattices(c % fill) % obj

        select type (lat)

        type is (RectLattice)

          ! Loop over lattice coordinates
          do j = 1, lat % n_cells(1)
            do k = 1, lat % n_cells(2)
              do m = 1, lat % n_cells(3)
                lat % offset(map, j, k, m) = offset
                next_univ => universes(lat % universes(j, k, m))
                offset = offset + &
                     count_target(next_univ, counts, found, univ_id, map)
              end do
            end do
          end do

        type is (HexLattice)

          ! Loop over lattice coordinates
          do m = 1, lat % n_axial
            do k = 1, 2*lat % n_rings - 1
              do j = 1, 2*lat % n_rings - 1
                ! This array location is never used
                if (j + k < lat % n_rings + 1) then
                  cycle
                ! This array location is never used
                else if (j + k > 3*lat % n_rings - 1) then
                  cycle
                else
                  lat % offset(map, j, k, m) = offset
                  next_univ => universes(lat % universes(j, k, m))
                  offset = offset + &
                       count_target(next_univ, counts, found, univ_id, map)
                end if
              end do
            end do
          end do
        end select

      end if
    end do

  end subroutine calc_offsets

!===============================================================================
! COUNT_TARGET recursively totals the numbers of occurances of a given
! universe ID beginning with the universe given.
!===============================================================================

  recursive function count_target(univ, counts, found, univ_id, map) result(count)

    type(Universe), intent(in) :: univ         ! universe to search through
    integer, intent(inout)     :: counts(:,:)  ! target count
    logical, intent(inout)     :: found(:,:)   ! target found
    integer, intent(in)        :: univ_id         ! target universe ID
    integer, intent(in)        :: map          ! current map

    integer :: i                           ! index over cells
    integer :: j, k, m                     ! indices in lattice
    integer :: n                           ! number of cells to search
    integer :: cell_index                  ! index in cells array
    integer :: count                       ! number of times target located
    type(Cell),     pointer :: c           ! pointer to current cell
    type(Universe), pointer :: next_univ   ! next univ to loop through
    class(Lattice), pointer :: lat         ! pointer to current lattice

    ! Don't research places already checked
    if (found(universe_dict % get(univ % id), map)) then
      count = counts(universe_dict % get(univ % id), map)
      return
    end if

    ! If this is the target, it can't contain itself.
    ! Count = 1, then quit
    if (univ % id == univ_id) then
      count = 1
      counts(universe_dict % get(univ % id), map) = 1
      found(universe_dict % get(univ % id), map) = .true.
      return
    end if

    count = 0
    n = size(univ % cells)

    do i = 1, n

      cell_index = univ % cells(i)

      ! get pointer to cell
      c => cells(cell_index)

      ! ====================================================================
      ! AT LOWEST UNIVERSE, TERMINATE SEARCH
      if (c % type == FILL_MATERIAL) then

      ! ====================================================================
      ! CELL CONTAINS LOWER UNIVERSE, RECURSIVELY FIND CELL
      elseif (c % type == FILL_UNIVERSE) then

        next_univ => universes(c % fill)

        ! Found target - stop since target cannot contain itself
        if (next_univ % id == univ_id) then
          count = count + 1
          return
        end if

        count = count + count_target(next_univ, counts, found, univ_id, map)
        c => cells(cell_index)

      ! ====================================================================
      ! CELL CONTAINS LATTICE, RECURSIVELY FIND CELL
      elseif (c % type == FILL_LATTICE) then

        ! Set current lattice
        lat => lattices(c % fill) % obj

        select type (lat)

        type is (RectLattice)

          ! Loop over lattice coordinates
          do j = 1, lat % n_cells(1)
            do k = 1, lat % n_cells(2)
              do m = 1, lat % n_cells(3)
                next_univ => universes(lat % universes(j, k, m))

                ! Found target - stop since target cannot contain itself
                if (next_univ % id == univ_id) then
                  count = count + 1
                  cycle
                end if

                count = count + &
                     count_target(next_univ, counts, found, univ_id, map)

              end do
            end do
          end do

          type is (HexLattice)

            ! Loop over lattice coordinates
            do m = 1, lat % n_axial
              do k = 1, 2*lat % n_rings - 1
                do j = 1, 2*lat % n_rings - 1
                  ! This array location is never used
                  if (j + k < lat % n_rings + 1) then
                    cycle
                  ! This array location is never used
                  else if (j + k > 3*lat % n_rings - 1) then
                    cycle
                  else
                    next_univ => universes(lat % universes(j, k, m))

                    ! Found target - stop since target cannot contain itself
                    if (next_univ % id == univ_id) then
                      count = count + 1
                      cycle
                    end if

                    count = count + &
                         count_target(next_univ, counts, found, univ_id, map)
                  end if
                end do
              end do
            end do

          end select

      end if
    end do

    counts(universe_dict % get(univ % id), map) = count
    found(universe_dict % get(univ % id), map) = .true.

  end function count_target

!===============================================================================
! COUNT_INSTANCE recursively totals the number of occurrences of all cells
! beginning with the universe given.
!===============================================================================

  recursive subroutine count_instance(univ)

    type(Universe), intent(in) :: univ  ! universe to search through

    integer :: i        ! index over cells
    integer :: j, k, m  ! indices in lattice
    integer :: n        ! number of cells to search

    n = size(univ % cells)

    do i = 1, n
      associate (c => cells(univ % cells(i)))
        c % instances = c % instances + 1

        ! ====================================================================
        ! CELL CONTAINS LOWER UNIVERSE, RECURSIVELY FIND CELL
        if (c % type == FILL_UNIVERSE) then

          call count_instance(universes(c % fill))

        ! ====================================================================
        ! CELL CONTAINS LATTICE, RECURSIVELY FIND CELL
        elseif (c % type == FILL_LATTICE) then

          ! Set current lattice
          associate (lat => lattices(c % fill) % obj)

            select type (lat)
            type is (RectLattice)

              ! Loop over lattice coordinates
              do j = 1, lat % n_cells(1)
                do k = 1, lat % n_cells(2)
                  do m = 1, lat % n_cells(3)
                    call count_instance(universes(lat % universes(j, k, m)))
                  end do
                end do
              end do

            type is (HexLattice)

              ! Loop over lattice coordinates
              do m = 1, lat % n_axial
                do k = 1, 2*lat % n_rings - 1
                  do j = 1, 2*lat % n_rings - 1
                    ! This array location is never used
                    if (j + k < lat % n_rings + 1) then
                      cycle
                      ! This array location is never used
                    else if (j + k > 3*lat % n_rings - 1) then
                      cycle
                    else
                      call count_instance(universes(lat % universes(j, k, m)))
                    end if
                  end do
                end do
              end do

            end select
          end associate
        end if
      end associate
    end do

  end subroutine count_instance

!===============================================================================
! MAXIMUM_LEVELS determines the maximum number of nested coordinate levels in
! the geometry
!===============================================================================

  recursive function maximum_levels(univ) result(levels)

    type(Universe), intent(in) :: univ  ! universe to search through
    integer :: levels                   ! maximum number of levels for this universe

    integer :: i                          ! index over cells
    integer :: j, k, m                    ! indices in lattice
    integer :: levels_below               ! max levels below this universe
    type(Cell),     pointer :: c          ! pointer to current cell
    type(Universe), pointer :: next_univ  ! next universe to loop through
    class(Lattice), pointer :: lat        ! pointer to current lattice

    levels_below = 0
    do i = 1, size(univ % cells)
      c => cells(univ % cells(i))

      ! ====================================================================
      ! CELL CONTAINS LOWER UNIVERSE, RECURSIVELY FIND CELL
      if (c % type == FILL_UNIVERSE) then

        next_univ => universes(c % fill)
        levels_below = max(levels_below, maximum_levels(next_univ))

      ! ====================================================================
      ! CELL CONTAINS LATTICE, RECURSIVELY FIND CELL
      elseif (c % type == FILL_LATTICE) then

        ! Set current lattice
        lat => lattices(c % fill) % obj

        select type (lat)

        type is (RectLattice)

          ! Loop over lattice coordinates
          do j = 1, lat % n_cells(1)
            do k = 1, lat % n_cells(2)
              do m = 1, lat % n_cells(3)
                next_univ => universes(lat % universes(j, k, m))
                levels_below = max(levels_below, maximum_levels(next_univ))
              end do
            end do
          end do

        type is (HexLattice)

          ! Loop over lattice coordinates
          do m = 1, lat % n_axial
            do k = 1, 2*lat % n_rings - 1
              do j = 1, 2*lat % n_rings - 1
                ! This array location is never used
                if (j + k < lat % n_rings + 1) then
                  cycle
                ! This array location is never used
                else if (j + k > 3*lat % n_rings - 1) then
                  cycle
                else
                  next_univ => universes(lat % universes(j, k, m))
                  levels_below = max(levels_below, maximum_levels(next_univ))
                end if
              end do
            end do
          end do

        end select

      end if
    end do

    levels = 1 + levels_below

  end function maximum_levels

end module geometry
