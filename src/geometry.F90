module geometry

  use constants
  use error,                  only: fatal_error, warning
  use geometry_header,        only: Cell, Universe, Lattice, &
                                    &RectLattice, HexLattice
  use global
  use output,                 only: write_message
  use particle_header,        only: LocalCoord, Particle
  use particle_restart_write, only: write_particle_restart
  use surface_header
  use stl_vector,             only: VectorInt
  use string,                 only: to_str
  use tally,                  only: score_surface_current

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
      n = univ % n_cells

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
    integer :: index_cell           ! index in cells array
    logical :: use_search_cells     ! use cells provided as argument
    type(Cell),     pointer :: c    ! pointer to cell
    class(Lattice), pointer :: lat  ! pointer to lattice
    type(Universe), pointer :: univ ! universe to search in

    do j = p % n_coord + 1, MAX_COORD
      call p % coord(j) % reset()
    end do
    j = p % n_coord

    ! set size of list to search
    if (present(search_cells)) then
      use_search_cells = .true.
      n = size(search_cells)
    else
      use_search_cells = .false.
      univ => universes(p % coord(j) % universe)
      n = univ % n_cells
    end if

    CELL_LOOP: do i = 1, n
      ! select cells based on whether we are searching a universe or a provided
      ! list of cells (this would be for lists of neighbor cells)
      if (use_search_cells) then
        index_cell = search_cells(i)
        ! check to make sure search cell is in same universe
        if (cells(index_cell) % universe /= p % coord(j) % universe) cycle
      else
        index_cell = univ % cells(i)
      end if

      ! get pointer to cell
      c => cells(index_cell)

      ! Move on to the next cell if the particle is not inside this cell
      if (.not. cell_contains(c, p)) cycle

      ! Set cell on this level
      p % coord(j) % cell = index_cell

      ! Show cell information on trace
      if (verbosity >= 10 .or. trace) then
        call write_message("    Entering cell " // trim(to_str(c % id)))
      end if

      CELL_TYPE: if (c % type == CELL_NORMAL) then
        ! ======================================================================
        ! AT LOWEST UNIVERSE, TERMINATE SEARCH

        ! Set the particle material
        p % last_material = p % material
        if (size(c % material) == 1) then
          ! Only one material for this cell; assign that one to the particle.
          p % material = c % material(1)
        else
          ! Distributed instances of this cell have different materials.
          ! Determine which instance this is and assign the matching material.
          distribcell_index = c % distribcell_index
          offset = 0
          do k = 1, p % n_coord
            if (cells(p % coord(k) % cell) % type == CELL_FILL) then
              offset = offset + cells(p % coord(k) % cell) % &
                   offset(distribcell_index)
            elseif (cells(p % coord(k) % cell) % type == CELL_LATTICE) then
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
          p % material = c % material(offset + 1)
        end if

      elseif (c % type == CELL_FILL) then CELL_TYPE
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
        if (.not. found) exit

      elseif (c % type == CELL_LATTICE) then CELL_TYPE
        ! ======================================================================
        ! CELL CONTAINS LATTICE, RECURSIVELY FIND CELL

        ! Set current lattice
        lat => lattices(c % fill) % obj

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
            call handle_lost_particle(p, "Particle " // trim(to_str(p %id)) &
                 // " is outside lattice " // trim(to_str(lat % id)) &
                 // " but the lattice has no defined outer universe.")
            return
          else
            p % coord(j + 1) % universe = lat % outer
          end if
        end if

        ! Move particle to next level and search for the lower cells.
        j = j + 1
        p % n_coord = j

        call find_cell(p, found)
        j = p % n_coord
        if (.not. found) exit

      end if CELL_TYPE

      ! Found cell so we can return
      found = .true.
      return
    end do CELL_LOOP

    found = .false.

  end subroutine find_cell

!===============================================================================
! CROSS_SURFACE handles all surface crossings, whether the particle leaks out of
! the geometry, is reflected, or crosses into a new lattice or cell
!===============================================================================

  subroutine cross_surface(p, last_cell)
    type(Particle), intent(inout) :: p
    integer,        intent(in)    :: last_cell  ! last cell particle was in

    real(8) :: u         ! x-component of direction
    real(8) :: v         ! y-component of direction
    real(8) :: w         ! z-component of direction
    real(8) :: norm      ! "norm" of surface normal
    real(8) :: xyz(3)    ! Saved global coordinate
    integer :: i_surface ! index in surfaces
    logical :: found     ! particle found in universe?
    class(Surface), pointer :: surf

    i_surface = abs(p % surface)
    surf => surfaces(i_surface)%obj
    if (verbosity >= 10 .or. trace) then
      call write_message("    Crossing surface " // trim(to_str(surf % id)))
    end if

    if (surf % bc == BC_VACUUM .and. (run_mode /= MODE_PLOTTING)) then
      ! =======================================================================
      ! PARTICLE LEAKS OUT OF PROBLEM

      ! Kill particle
      p % alive = .false.

      ! Score any surface current tallies -- note that the particle is moved
      ! forward slightly so that if the mesh boundary is on the surface, it is
      ! still processed

      if (active_current_tallies % size() > 0) then
        ! TODO: Find a better solution to score surface currents than
        ! physically moving the particle forward slightly

        p % coord(1) % xyz = p % coord(1) % xyz + TINY_BIT * p % coord(1) % uvw
        call score_surface_current(p)
      end if

      ! Score to global leakage tally
      if (tallies_on) then
        global_tally_leakage = global_tally_leakage + p % wgt
      end if

      ! Display message
      if (verbosity >= 10 .or. trace) then
        call write_message("    Leaked out of surface " &
             &// trim(to_str(surf % id)))
      end if
      return

    elseif (surf % bc == BC_REFLECT .and. (run_mode /= MODE_PLOTTING)) then
      ! =======================================================================
      ! PARTICLE REFLECTS FROM SURFACE

      ! Do not handle reflective boundary conditions on lower universes
      if (p % n_coord /= 1) then
        call handle_lost_particle(p, "Cannot reflect particle " &
             &// trim(to_str(p % id)) // " off surface in a lower universe.")
        return
      end if

      ! Score surface currents since reflection causes the direction of the
      ! particle to change -- artificially move the particle slightly back in
      ! case the surface crossing is coincident with a mesh boundary

      if (active_current_tallies % size() > 0) then
        xyz = p % coord(1) % xyz
        p % coord(1) % xyz = p % coord(1) % xyz - TINY_BIT * p % coord(1) % uvw
        call score_surface_current(p)
        p % coord(1) % xyz = xyz
      end if

      ! Reflect particle off surface
      call surf%reflect(p%coord(1)%xyz, p%coord(1)%uvw)

      ! Make sure new particle direction is normalized
      u = p%coord(1)%uvw(1)
      v = p%coord(1)%uvw(2)
      w = p%coord(1)%uvw(3)
      norm = sqrt(u*u + v*v + w*w)
      p%coord(1)%uvw(:) = [u, v, w] / norm

      ! Reassign particle's cell and surface
      p % coord(1) % cell = last_cell
      p % surface = -p % surface

      ! If a reflective surface is coincident with a lattice or universe
      ! boundary, it is necessary to redetermine the particle's coordinates in
      ! the lower universes.

      p % n_coord = 1
      call find_cell(p, found)
      if (.not. found) then
        call handle_lost_particle(p, "Couldn't find particle after reflecting&
             & from surface " // trim(to_str(surf%id)) // ".")
        return
      end if

      ! Set previous coordinate going slightly past surface crossing
      p % last_xyz = p % coord(1) % xyz + TINY_BIT * p % coord(1) % uvw

      ! Diagnostic message
      if (verbosity >= 10 .or. trace) then
        call write_message("    Reflected from surface " &
             &// trim(to_str(surf%id)))
      end if
      return
    elseif (surf % bc == BC_PERIODIC .and. run_mode /= MODE_PLOTTING) then
      ! =======================================================================
      ! PERIODIC BOUNDARY

      ! Do not handle periodic boundary conditions on lower universes
      if (p % n_coord /= 1) then
        call handle_lost_particle(p, "Cannot transfer particle " &
             // trim(to_str(p % id)) // " across surface in a lower universe.&
             & Boundary conditions must be applied to universe 0.")
        return
      end if

      ! Score surface currents since reflection causes the direction of the
      ! particle to change -- artificially move the particle slightly back in
      ! case the surface crossing is coincident with a mesh boundary

      if (active_current_tallies % size() > 0) then
        xyz = p % coord(1) % xyz
        p % coord(1) % xyz = p % coord(1) % xyz - TINY_BIT * p % coord(1) % uvw
        call score_surface_current(p)
        p % coord(1) % xyz = xyz
      end if

      select type (surf)
      type is (SurfaceXPlane)
        select type (opposite => surfaces(surf % i_periodic) % obj)
        type is (SurfaceXPlane)
          p % coord(1) % xyz(1) = opposite % x0
        end select

      type is (SurfaceYPlane)
        select type (opposite => surfaces(surf % i_periodic) % obj)
        type is (SurfaceYPlane)
          p % coord(1) % xyz(2) = opposite % y0
        end select

      type is (SurfaceZPlane)
        select type (opposite => surfaces(surf % i_periodic) % obj)
        type is (SurfaceZPlane)
          p % coord(1) % xyz(3) = opposite % z0
        end select
      end select

      ! Reassign particle's surface
      p % surface = sign(surf % i_periodic, p % surface)

      ! Figure out what cell particle is in now
      p % n_coord = 1
      call find_cell(p, found)
      if (.not. found) then
        call handle_lost_particle(p, "Couldn't find particle after hitting &
             &periodic boundary on surface " // trim(to_str(surf%id)) // ".")
        return
      end if

      ! Set previous coordinate going slightly past surface crossing
      p % last_xyz = p % coord(1) % xyz + TINY_BIT * p % coord(1) % uvw

      ! Diagnostic message
      if (verbosity >= 10 .or. trace) then
        call write_message("    Hit periodic boundary on surface " &
             // trim(to_str(surf%id)))
      end if
      return
    end if

    ! ==========================================================================
    ! SEARCH NEIGHBOR LISTS FOR NEXT CELL

    if (p % surface > 0 .and. allocated(surf%neighbor_pos)) then
      ! If coming from negative side of surface, search all the neighboring
      ! cells on the positive side

      call find_cell(p, found, surf%neighbor_pos)
      if (found) return

    elseif (p % surface < 0  .and. allocated(surf%neighbor_neg)) then
      ! If coming from positive side of surface, search all the neighboring
      ! cells on the negative side

      call find_cell(p, found, surf%neighbor_neg)
      if (found) return

    end if

    ! ==========================================================================
    ! COULDN'T FIND PARTICLE IN NEIGHBORING CELLS, SEARCH ALL CELLS

    ! Remove lower coordinate levels and assignment of surface
    p % surface = NONE
    p % n_coord = 1
    call find_cell(p, found)

    if (run_mode /= MODE_PLOTTING .and. (.not. found)) then
      ! If a cell is still not found, there are two possible causes: 1) there is
      ! a void in the model, and 2) the particle hit a surface at a tangent. If
      ! the particle is really traveling tangent to a surface, if we move it
      ! forward a tiny bit it should fix the problem.

      p % n_coord = 1
      p % coord(1) % xyz = p % coord(1) % xyz + TINY_BIT * p % coord(1) % uvw
      call find_cell(p, found)

      ! Couldn't find next cell anywhere! This probably means there is an actual
      ! undefined region in the geometry.

      if (.not. found) then
        call handle_lost_particle(p, "After particle " // trim(to_str(p % id)) &
             // " crossed surface " // trim(to_str(surf%id)) &
             // " it could not be located in any cell and it did not leak.")
        return
      end if
    end if

  end subroutine cross_surface

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
          call handle_lost_particle(p, "Could not locate particle " &
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
          call handle_lost_particle(p, "Could not locate particle " &
               // trim(to_str(p % id)) &
               // " after crossing a lattice boundary.")
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
          call handle_lost_particle(p, "Particle " // trim(to_str(p % id)) &
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
         4)

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
! HANDLE_LOST_PARTICLE
!===============================================================================

  subroutine handle_lost_particle(p, message)

    type(Particle), intent(inout) :: p
    character(*)                  :: message

    integer(8) :: tot_n_particles

    ! Print warning and write lost particle file
    call warning(message)
    call write_particle_restart(p)

    ! Increment number of lost particles
    p % alive = .false.
!$omp atomic
    n_lost_particles = n_lost_particles + 1

    ! Count the total number of simulated particles (on this processor)
    tot_n_particles = n_batches * gen_per_batch * work

    ! Abort the simulation if the maximum number of lost particles has been
    ! reached
    if (n_lost_particles >= MAX_LOST_PARTICLES .and. &
         n_lost_particles >= REL_MAX_LOST_PARTICLES * tot_n_particles) then
      call fatal_error("Maximum number of lost particles has been reached.")
    end if

  end subroutine handle_lost_particle

!===============================================================================
! CALC_OFFSETS calculates and stores the offsets in all fill cells. This
! routine is called once upon initialization.
!===============================================================================

  subroutine calc_offsets(goal, map, univ, counts, found)

    integer, intent(in)        :: goal         ! target universe ID
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

    n = univ % n_cells
    offset = 0

    do i = 1, n

      cell_index = univ % cells(i)

      ! get pointer to cell
      c => cells(cell_index)

      ! ====================================================================
      ! AT LOWEST UNIVERSE, TERMINATE SEARCH
      if (c % type == CELL_NORMAL) then

      ! ====================================================================
      ! CELL CONTAINS LOWER UNIVERSE, RECURSIVELY FIND CELL
      elseif (c % type == CELL_FILL) then
        ! Set offset for the cell on this level
        c % offset(map) = offset

        ! Count contents of this cell
        next_univ => universes(c % fill)
        offset = offset + count_target(next_univ, counts, found, goal, map)

        ! Move into the next universe
        next_univ => universes(c % fill)
        c => cells(cell_index)

      ! ====================================================================
      ! CELL CONTAINS LATTICE, RECURSIVELY FIND CELL
      elseif (c % type == CELL_LATTICE) then

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
                     count_target(next_univ, counts, found, goal, map)
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
                       count_target(next_univ, counts, found, goal, map)
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

  recursive function count_target(univ, counts, found, goal, map) result(count)

    type(Universe), intent(in) :: univ         ! universe to search through
    integer, intent(inout)     :: counts(:,:)  ! target count
    logical, intent(inout)     :: found(:,:)   ! target found
    integer, intent(in)        :: goal         ! target universe ID
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
    if (found(universe_dict % get_key(univ % id), map)) then
      count = counts(universe_dict % get_key(univ % id), map)
      return
    end if

    ! If this is the target, it can't contain itself.
    ! Count = 1, then quit
    if (univ % id == goal) then
      count = 1
      counts(universe_dict % get_key(univ % id), map) = 1
      found(universe_dict % get_key(univ % id), map) = .true.
      return
    end if

    count = 0
    n = univ % n_cells

    do i = 1, n

      cell_index = univ % cells(i)

      ! get pointer to cell
      c => cells(cell_index)

      ! ====================================================================
      ! AT LOWEST UNIVERSE, TERMINATE SEARCH
      if (c % type == CELL_NORMAL) then

      ! ====================================================================
      ! CELL CONTAINS LOWER UNIVERSE, RECURSIVELY FIND CELL
      elseif (c % type == CELL_FILL) then

        next_univ => universes(c % fill)

        ! Found target - stop since target cannot contain itself
        if (next_univ % id == goal) then
          count = count + 1
          return
        end if

        count = count + count_target(next_univ, counts, found, goal, map)
        c => cells(cell_index)

      ! ====================================================================
      ! CELL CONTAINS LATTICE, RECURSIVELY FIND CELL
      elseif (c % type == CELL_LATTICE) then

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
                if (next_univ % id == goal) then
                  count = count + 1
                  cycle
                end if

                count = count + &
                     count_target(next_univ, counts, found, goal, map)

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
                    if (next_univ % id == goal) then
                      count = count + 1
                      cycle
                    end if

                    count = count + &
                         count_target(next_univ, counts, found, goal, map)
                  end if
                end do
              end do
            end do

          end select

      end if
    end do

    counts(universe_dict % get_key(univ % id), map) = count
    found(universe_dict % get_key(univ % id), map) = .true.

  end function count_target

!===============================================================================
! COUNT_INSTANCE recursively totals the number of occurrences of all cells
! beginning with the universe given.
!===============================================================================

  recursive subroutine count_instance(univ)

    type(Universe), intent(in) :: univ  ! universe to search through

    integer :: i                          ! index over cells
    integer :: j, k, m                    ! indices in lattice
    integer :: n                          ! number of cells to search
    integer :: cell_index                 ! index in cells array
    type(Cell),     pointer :: c          ! pointer to current cell
    type(Universe), pointer :: next_univ  ! next universe to loop through
    class(Lattice), pointer :: lat        ! pointer to current lattice

    n = univ % n_cells

    do i = 1, n

      cell_index = univ % cells(i)

      ! get pointer to cell
      c => cells(cell_index)
      c % instances = c % instances + 1

      ! ====================================================================
      ! AT LOWEST UNIVERSE, TERMINATE SEARCH
      if (c % type == CELL_NORMAL) then

      ! ====================================================================
      ! CELL CONTAINS LOWER UNIVERSE, RECURSIVELY FIND CELL
      elseif (c % type == CELL_FILL) then

        next_univ => universes(c % fill)

        call count_instance(next_univ)
        c => cells(cell_index)

      ! ====================================================================
      ! CELL CONTAINS LATTICE, RECURSIVELY FIND CELL
      elseif (c % type == CELL_LATTICE) then

        ! Set current lattice
        lat => lattices(c % fill) % obj

        select type (lat)

        type is (RectLattice)

          ! Loop over lattice coordinates
          do j = 1, lat % n_cells(1)
            do k = 1, lat % n_cells(2)
              do m = 1, lat % n_cells(3)
                next_univ => universes(lat % universes(j, k, m))
                call count_instance(next_univ)
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
                  call count_instance(next_univ)
                end if
              end do
            end do
          end do

        end select

      end if
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
    do i = 1, univ % n_cells
      c => cells(univ % cells(i))

      ! ====================================================================
      ! CELL CONTAINS LOWER UNIVERSE, RECURSIVELY FIND CELL
      if (c % type == CELL_FILL) then

        next_univ => universes(c % fill)
        levels_below = max(levels_below, maximum_levels(next_univ))

      ! ====================================================================
      ! CELL CONTAINS LATTICE, RECURSIVELY FIND CELL
      elseif (c % type == CELL_LATTICE) then

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
