module geometry

  use constants
  use error,                  only: fatal_error, warning
  use geometry_header,        only: Cell, Surface, Universe, Lattice, &
                                    &RectLattice, HexLattice
  use global
  use output,                 only: write_message
  use particle_header,        only: LocalCoord, deallocate_coord, Particle
  use particle_restart_write, only: write_particle_restart
  use string,                 only: to_str
  use tally,                  only: score_surface_current

  implicit none

contains

!===============================================================================
! SIMPLE_CELL_CONTAINS determines whether a given the current coordinates of the
! particle are inside a cell defined as the intersection of a series of surfaces
!===============================================================================

  function simple_cell_contains(c, p) result(in_cell)

    type(Cell),     pointer       :: c
    type(Particle), intent(inout) :: p
    logical                       :: in_cell

    integer :: i               ! index of surfaces in cell
    integer :: i_surface       ! index in surfaces array (with sign)
    logical :: specified_sense ! specified sense of surface in list
    logical :: actual_sense    ! sense of particle wrt surface
    type(Surface), pointer :: s

    SURFACE_LOOP: do i = 1, c % n_surfaces
      ! Lookup surface
      i_surface = c % surfaces(i)

      ! Check if the particle is currently on the specified surface
      if (i_surface == p % surface) then
        ! Particle is heading into the cell
        cycle
      elseif (i_surface == -p % surface) then
        ! Particle is heading out of the cell
        in_cell = .false.
        return
      end if

      ! Determine the specified sense of the surface in the cell and the actual
      ! sense of the particle with respect to the surface
      s => surfaces(abs(i_surface))
      actual_sense = sense(p, s)
      specified_sense = (c % surfaces(i) > 0)

      ! Compare sense of point to specified sense
      if (actual_sense .neqv. specified_sense) then
        in_cell = .false.
        return
      end if
    end do SURFACE_LOOP

    ! If we've reached here, then the sense matched on every surface or there
    ! are no surfaces.
    in_cell = .true.

  end function simple_cell_contains


!===============================================================================
! CHECK_CELL_OVERLAP checks for overlapping cells at the current particle's
! position using simple_cell_contains and the LocalCoord's built up by find_cell
!===============================================================================

  subroutine check_cell_overlap(p)

    type(Particle), intent(inout) :: p

    integer :: i                       ! cell loop index on a level
    integer :: n                       ! number of cells to search on a level
    integer :: index_cell              ! index in cells array
    type(Cell),       pointer :: c     ! pointer to cell
    type(Universe),   pointer :: univ  ! universe to search in
    type(LocalCoord), pointer :: coord ! particle coordinate to search on

    coord => p % coord0

    ! loop through each coordinate level
    do while (associated(coord))

      p % coord => coord

      univ => universes(coord % universe)
      n = univ % n_cells

      ! loop through each cell on this level
      do i = 1, n
        index_cell = univ % cells(i)
        c => cells(index_cell)

        if (simple_cell_contains(c, p)) then
          ! the particle should only be contained in one cell per level
          if (index_cell /= coord % cell) then
            call fatal_error("Overlapping cells detected: " &
                 &// trim(to_str(cells(index_cell) % id)) // ", " &
                 &// trim(to_str(cells(coord % cell) % id)) &
                 &// " on universe " // trim(to_str(univ % id)))
          end if

          overlap_check_cnt(index_cell) = overlap_check_cnt(index_cell) + 1

        end if

      end do

      coord => coord % next

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
    integer :: i_xyz(3)             ! indices in lattice
    integer :: n                    ! number of cells to search
    integer :: index_cell           ! index in cells array
    logical :: use_search_cells     ! use cells provided as argument
    type(Cell),     pointer :: c    ! pointer to cell
    class(Lattice), pointer :: lat  ! pointer to lattice
    type(Universe), pointer :: univ ! universe to search in

    ! Remove coordinates for any lower levels
    call deallocate_coord(p % coord % next)

    ! set size of list to search
    if (present(search_cells)) then
      use_search_cells = .true.
      n = size(search_cells)
    else
      use_search_cells = .false.
      univ => universes(p % coord % universe)
      n = univ % n_cells
    end if

    CELL_LOOP: do i = 1, n
      ! select cells based on whether we are searching a universe or a provided
      ! list of cells (this would be for lists of neighbor cells)
      if (use_search_cells) then
        index_cell = search_cells(i)
        ! check to make sure search cell is in same universe
        if (cells(index_cell) % universe /= p % coord % universe) cycle
      else
        index_cell = univ % cells(i)
      end if

      ! get pointer to cell
      c => cells(index_cell)

      ! Move on to the next cell if the particle is not inside this cell
      if (.not. simple_cell_contains(c, p)) cycle

      ! Set cell on this level
      p % coord % cell = index_cell

      ! Show cell information on trace
      if (verbosity >= 10 .or. trace) then
        call write_message("    Entering cell " // trim(to_str(c % id)))
      end if

      CELL_TYPE: if (c % type == CELL_NORMAL) then
        ! ======================================================================
        ! AT LOWEST UNIVERSE, TERMINATE SEARCH

        ! set material
        p % last_material = p % material
        p % material = c % material

      elseif (c % type == CELL_FILL) then CELL_TYPE
        ! ======================================================================
        ! CELL CONTAINS LOWER UNIVERSE, RECURSIVELY FIND CELL

        ! Create new level of coordinates
        allocate(p % coord % next)
        p % coord % next % xyz = p % coord % xyz
        p % coord % next % uvw = p % coord % uvw

        ! Move particle to next level and set universe
        p % coord => p % coord % next
        p % coord % universe = c % fill

        ! Apply translation
        if (allocated(c % translation)) then
          p % coord % xyz = p % coord % xyz - c % translation
        end if

        ! Apply rotation
        if (allocated(c % rotation_matrix)) then
          p % coord % xyz = matmul(c % rotation_matrix, p % coord % xyz)
          p % coord % uvw = matmul(c % rotation_matrix, p % coord % uvw)
          p % coord % rotated = .true.
        end if

        call find_cell(p, found)
        if (.not. found) exit

      elseif (c % type == CELL_LATTICE) then CELL_TYPE
        ! ======================================================================
        ! CELL CONTAINS LATTICE, RECURSIVELY FIND CELL

        ! Set current lattice
        lat => lattices(c % fill) % obj

        ! Determine lattice indices
        i_xyz = lat % get_indices(p % coord % xyz + TINY_BIT * p % coord % uvw)

        ! Create new level of coordinates
        allocate(p % coord % next)
        p % coord % next % xyz = lat % get_local_xyz(p % coord % xyz, i_xyz)
        p % coord % next % uvw = p % coord % uvw

        ! set particle lattice indices
        p % coord % next% lattice   = c % fill
        p % coord % next% lattice_x = i_xyz(1)
        p % coord % next% lattice_y = i_xyz(2)
        p % coord % next% lattice_z = i_xyz(3)

        ! Set the next lowest coordinate level.
        if (lat % are_valid_indices(i_xyz)) then
          ! Particle is inside the lattice.
          p % coord % next % universe = &
               &lat % universes(i_xyz(1), i_xyz(2), i_xyz(3))

        else
          ! Particle is outside the lattice.
          if (lat % outer == NO_OUTER_UNIVERSE) then
            call fatal_error("A particle is outside latttice " &
                 &// trim(to_str(lat % id)) // " but the lattice has no &
                 &defined outer universe.")
          else
            p % coord % next % universe = lat % outer
          end if
        end if

        ! Move particle to next level and search for the lower cells.
        p % coord => p % coord % next

        call find_cell(p, found)
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

    real(8) :: x         ! x-x0 for sphere
    real(8) :: y         ! y-y0 for sphere
    real(8) :: z         ! z-z0 for sphere
    real(8) :: R         ! radius of sphere
    real(8) :: u         ! x-component of direction
    real(8) :: v         ! y-component of direction
    real(8) :: w         ! z-component of direction
    real(8) :: n1        ! x-component of surface normal
    real(8) :: n2        ! y-component of surface normal
    real(8) :: n3        ! z-component of surface normal
    real(8) :: dot_prod  ! dot product of direction and normal
    real(8) :: norm      ! "norm" of surface normal
    integer :: i_surface ! index in surfaces
    logical :: found     ! particle found in universe?
    type(Surface), pointer :: surf

    i_surface = abs(p % surface)
    surf => surfaces(i_surface)
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

        p % coord0 % xyz = p % coord0 % xyz + TINY_BIT * p % coord0 % uvw
        call score_surface_current(p)
      end if

      ! Score to global leakage tally
      if (tallies_on) then
!$omp atomic
        global_tallies(LEAKAGE) % value = &
           global_tallies(LEAKAGE) % value + p % wgt
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
      if (.not. associated(p % coord, p % coord0)) then
        call handle_lost_particle(p, "Cannot reflect particle " &
             &// trim(to_str(p % id)) // " off surface in a lower universe.")
        return
      end if

      ! Score surface currents since reflection causes the direction of the
      ! particle to change -- artificially move the particle slightly back in
      ! case the surface crossing in coincident with a mesh boundary

      if (active_current_tallies % size() > 0) then
        p % coord0 % xyz = p % coord0 % xyz - TINY_BIT * p % coord0 % uvw
        call score_surface_current(p)
        p % coord0 % xyz = p % coord0 % xyz + TINY_BIT * p % coord0 % uvw
      end if

      ! Copy particle's direction cosines
      u = p % coord0 % uvw(1)
      v = p % coord0 % uvw(2)
      w = p % coord0 % uvw(3)

      select case (surf%type)
      case (SURF_PX)
        u = -u

      case (SURF_PY)
        v = -v

      case (SURF_PZ)
        w = -w

      case (SURF_PLANE)
        ! Find surface coefficients and norm of vector normal to surface
        n1 = surf % coeffs(1)
        n2 = surf % coeffs(2)
        n3 = surf % coeffs(3)
        norm = n1*n1 + n2*n2 + n3*n3
        dot_prod = u*n1 + v*n2 + w*n3

        ! Reflect direction according to normal
        u = u - 2*dot_prod*n1/norm
        v = v - 2*dot_prod*n2/norm
        w = w - 2*dot_prod*n3/norm

      case (SURF_CYL_X)
        ! Find y-y0, z-z0 and dot product of direction and surface normal
        y = p % coord0 % xyz(2) - surf % coeffs(1)
        z = p % coord0 % xyz(3) - surf % coeffs(2)
        R = surf % coeffs(3)
        dot_prod = v*y + w*z

        ! Reflect direction according to normal
        v = v - 2*dot_prod*y/(R*R)
        w = w - 2*dot_prod*z/(R*R)

      case (SURF_CYL_Y)
        ! Find x-x0, z-z0 and dot product of direction and surface normal
        x = p % coord0 % xyz(1) - surf % coeffs(1)
        z = p % coord0 % xyz(3) - surf % coeffs(2)
        R = surf % coeffs(3)
        dot_prod = u*x + w*z

        ! Reflect direction according to normal
        u = u - 2*dot_prod*x/(R*R)
        w = w - 2*dot_prod*z/(R*R)

      case (SURF_CYL_Z)
        ! Find x-x0, y-y0 and dot product of direction and surface normal
        x = p % coord0 % xyz(1) - surf % coeffs(1)
        y = p % coord0 % xyz(2) - surf % coeffs(2)
        R = surf % coeffs(3)
        dot_prod = u*x + v*y

        ! Reflect direction according to normal
        u = u - 2*dot_prod*x/(R*R)
        v = v - 2*dot_prod*y/(R*R)

      case (SURF_SPHERE)
        ! Find x-x0, y-y0, z-z0 and dot product of direction and surface
        ! normal
        x = p % coord0 % xyz(1) - surf % coeffs(1)
        y = p % coord0 % xyz(2) - surf % coeffs(2)
        z = p % coord0 % xyz(3) - surf % coeffs(3)
        R = surf % coeffs(4)
        dot_prod = u*x + v*y + w*z

        ! Reflect direction according to normal
        u = u - 2*dot_prod*x/(R*R)
        v = v - 2*dot_prod*y/(R*R)
        w = w - 2*dot_prod*z/(R*R)

      case (SURF_CONE_X)
        ! Find x-x0, y-y0, z-z0 and dot product of direction and surface
        ! normal
        x = p % coord0 % xyz(1) - surf % coeffs(1)
        y = p % coord0 % xyz(2) - surf % coeffs(2)
        z = p % coord0 % xyz(3) - surf % coeffs(3)
        R = surf % coeffs(4)
        dot_prod = (v*y + w*z - R*u*x)/((R + ONE)*R*x*x)

        ! Reflect direction according to normal
        u = u + 2*dot_prod*R*x
        v = v - 2*dot_prod*y
        w = w - 2*dot_prod*z

      case (SURF_CONE_Y)
        ! Find x-x0, y-y0, z-z0 and dot product of direction and surface
        ! normal
        x = p % coord0 % xyz(1) - surf % coeffs(1)
        y = p % coord0 % xyz(2) - surf % coeffs(2)
        z = p % coord0 % xyz(3) - surf % coeffs(3)
        R = surf % coeffs(4)
        dot_prod = (u*x + w*z - R*v*y)/((R + ONE)*R*y*y)

        ! Reflect direction according to normal
        u = u - 2*dot_prod*x
        v = v + 2*dot_prod*R*y
        w = w - 2*dot_prod*z

      case (SURF_CONE_Z)
        ! Find x-x0, y-y0, z-z0 and dot product of direction and surface
        ! normal
        x = p % coord0 % xyz(1) - surf % coeffs(1)
        y = p % coord0 % xyz(2) - surf % coeffs(2)
        z = p % coord0 % xyz(3) - surf % coeffs(3)
        R = surf % coeffs(4)
        dot_prod = (u*x + v*y - R*w*z)/((R + ONE)*R*z*z)

        ! Reflect direction according to normal
        u = u - 2*dot_prod*x
        v = v - 2*dot_prod*y
        w = w + 2*dot_prod*R*z

      case default
        call fatal_error("Reflection not supported for surface " &
             &// trim(to_str(surf % id)))
      end select

      ! Set new particle direction
      norm = sqrt(u*u + v*v + w*w)
      p % coord0 % uvw = [u, v, w] / norm

      ! Reassign particle's cell and surface
      p % coord0 % cell = last_cell
      p % surface = -p % surface

      ! If a reflective surface is coincident with a lattice or universe
      ! boundary, it is necessary to redetermine the particle's coordinates in
      ! the lower universes.

      if (associated(p % coord0 % next)) then
        call deallocate_coord(p % coord0 % next)
        call find_cell(p, found)
        if (.not. found) then
          call handle_lost_particle(p, "Couldn't find particle after reflecting&
               & from surface.")
          return
        end if
      end if

      ! Set previous coordinate going slightly past surface crossing
      p % last_xyz = p % coord0 % xyz + TINY_BIT * p % coord0 % uvw

      ! Diagnostic message
      if (verbosity >= 10 .or. trace) then
        call write_message("    Reflected from surface " &
             &// trim(to_str(surf%id)))
      end if
      return
    end if

    ! ==========================================================================
    ! SEARCH NEIGHBOR LISTS FOR NEXT CELL

    if (p % surface > 0 .and. allocated(surf % neighbor_pos)) then
      ! If coming from negative side of surface, search all the neighboring
      ! cells on the positive side

      call find_cell(p, found, surf % neighbor_pos)
      if (found) return

    elseif (p % surface < 0  .and. allocated(surf % neighbor_neg)) then
      ! If coming from positive side of surface, search all the neighboring
      ! cells on the negative side

      call find_cell(p, found, surf % neighbor_neg)
      if (found) return

    end if

    ! ==========================================================================
    ! COULDN'T FIND PARTICLE IN NEIGHBORING CELLS, SEARCH ALL CELLS

    ! Remove lower coordinate levels and assignment of surface
    p % surface = NONE
    p % coord => p % coord0
    call deallocate_coord(p % coord % next)
    call find_cell(p, found)

    if (run_mode /= MODE_PLOTTING .and. (.not. found)) then
      ! If a cell is still not found, there are two possible causes: 1) there is
      ! a void in the model, and 2) the particle hit a surface at a tangent. If
      ! the particle is really traveling tangent to a surface, if we move it
      ! forward a tiny bit it should fix the problem.

      p % coord => p % coord0
      call deallocate_coord(p % coord % next)
      p % coord % xyz = p % coord % xyz + TINY_BIT * p % coord % uvw
      call find_cell(p, found)

      ! Couldn't find next cell anywhere! This probably means there is an actual
      ! undefined region in the geometry.

      if (.not. found) then
        call handle_lost_particle(p, "After particle " // trim(to_str(p % id)) &
             &// " crossed surface " // trim(to_str(surfaces(i_surface) % id)) &
             &// " it could not be located in any cell and it did not leak.")
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
    integer :: i_xyz(3)       ! indices in lattice
    logical :: found          ! particle found in cell?
    class(Lattice),   pointer :: lat
    type(LocalCoord), pointer :: parent_coord

    lat => lattices(p % coord % lattice) % obj

    if (verbosity >= 10 .or. trace) then
      call write_message("    Crossing lattice " // trim(to_str(lat % id)) &
           &// ". Current position (" // trim(to_str(p % coord % lattice_x)) &
           &// "," // trim(to_str(p % coord % lattice_y)) // "," &
           &// trim(to_str(p % coord % lattice_z)) // ")")
    end if

    ! Find the coordiante level just above the current one.
    parent_coord => p % coord0
    do while(.not. associated(parent_coord % next, p % coord))
      parent_coord => parent_coord % next
    end do

    ! Set the lattice indices.
    p % coord % lattice_x = p % coord % lattice_x + lattice_translation(1)
    p % coord % lattice_y = p % coord % lattice_y + lattice_translation(2)
    p % coord % lattice_z = p % coord % lattice_z + lattice_translation(3)
    i_xyz(1) = p % coord % lattice_x
    i_xyz(2) = p % coord % lattice_y
    i_xyz(3) = p % coord % lattice_z

    ! Set the new coordinate position.
    p % coord % xyz = lat % get_local_xyz(parent_coord % xyz, i_xyz)

    OUTSIDE_LAT: if (.not. lat % are_valid_indices(i_xyz)) then
      ! The particle is outside the lattice.  Search for it from coord0.
      call deallocate_coord(p % coord0 % next)
      p % coord => p % coord0
      call find_cell(p, found)
      if (.not. found) then
        call handle_lost_particle(p, "Could not locate particle " &
             &// trim(to_str(p % id)) // " after crossing a lattice boundary.")
        return
      end if

    else OUTSIDE_LAT

      ! Find cell in next lattice element
      p % coord % universe = lat % universes(i_xyz(1), i_xyz(2), i_xyz(3))

      call find_cell(p, found)
      if (.not. found) then
        ! In some circumstances, a particle crossing the corner of a cell may
        ! not be able to be found in the next universe. In this scenario we cut
        ! off all lower-level coordinates and search from universe zero

        ! Remove lower coordinates
        call deallocate_coord(p % coord0 % next)
        p % coord => p % coord0

        ! Search for particle
        call find_cell(p, found)
        if (.not. found) then
          call handle_lost_particle(p, "Could not locate particle " &
               &// trim(to_str(p % id)) &
               &// " after crossing a lattice boundary.")
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

  subroutine distance_to_boundary(p, dist, surface_crossed, lattice_translation)

    type(Particle), intent(inout) :: p
    real(8),        intent(out)   :: dist
    integer,        intent(out)   :: surface_crossed
    integer,        intent(out)   :: lattice_translation(3)

    integer :: i                  ! index for surface in cell
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
    real(8) :: r                  ! radius for quadratic surfaces
    real(8) :: tmp                ! dot product of surface normal with direction
    real(8) :: a,b,c,k            ! quadratic equation coefficients
    real(8) :: quad               ! discriminant of quadratic equation
    logical :: on_surface         ! is particle on surface?
    type(Cell),       pointer :: cl
    type(Surface),    pointer :: surf
    class(Lattice),   pointer :: lat
    type(LocalCoord), pointer :: coord
    type(LocalCoord), pointer :: final_coord
    type(LocalCoord), pointer :: parent_coord

    ! inialize distance to infinity (huge)
    dist = INFINITY
    d_lat = INFINITY
    d_surf = INFINITY
    lattice_translation(:) = [0, 0, 0]
    nullify(final_coord)

    ! Get pointer to top-level coordinates
    coord => p % coord0

    ! Loop over each universe level
    LEVEL_LOOP: do while(associated(coord))

      ! get pointer to cell on this level
      cl => cells(coord % cell)

      ! copy directional cosines
      u = coord % uvw(1)
      v = coord % uvw(2)
      w = coord % uvw(3)

      ! =======================================================================
      ! FIND MINIMUM DISTANCE TO SURFACE IN THIS CELL

      SURFACE_LOOP: do i = 1, cl % n_surfaces

        ! copy local coordinates of particle
        x = coord % xyz(1)
        y = coord % xyz(2)
        z = coord % xyz(3)

        ! check for coincident surface -- note that we can't skip the
        ! calculation in general because a particle could be on one side of a
        ! cylinder and still have a positive distance to the other

        index_surf = cl % surfaces(i)
        if (index_surf == p % surface) then
          on_surface = .true.
        else
          on_surface = .false.
        end if

        ! check for operators
        index_surf = abs(index_surf)
        if (index_surf >= OP_DIFFERENCE) cycle

        ! get pointer to surface
        surf => surfaces(index_surf)

        ! TODO: Can probably combines a lot of the cases to reduce repetition
        ! since the algorithm is the same for (x-plane, y-plane, z-plane),
        ! (x-cylinder, y-cylinder, z-cylinder), etc.

        select case (surf % type)
        case (SURF_PX)
          if (on_surface .or. u == ZERO) then
            d = INFINITY
          else
            x0 = surf % coeffs(1)
            d = (x0 - x)/u
            if (d < ZERO) d = INFINITY
          end if

        case (SURF_PY)
          if (on_surface .or. v == ZERO) then
            d = INFINITY
          else
            y0 = surf % coeffs(1)
            d = (y0 - y)/v
            if (d < ZERO) d = INFINITY
          end if

        case (SURF_PZ)
          if (on_surface .or. w == ZERO) then
            d = INFINITY
          else
            z0 = surf % coeffs(1)
            d = (z0 - z)/w
            if (d < ZERO) d = INFINITY
          end if

        case (SURF_PLANE)
          A = surf % coeffs(1)
          B = surf % coeffs(2)
          C = surf % coeffs(3)
          D = surf % coeffs(4)

          tmp = A*u + B*v + C*w
          if (on_surface .or. tmp == ZERO) then
            d = INFINITY
          else
            d = -(A*x + B*y + C*z - D)/tmp
            if (d < ZERO) d = INFINITY
          end if

        case (SURF_CYL_X)
          a = ONE - u*u  ! v^2 + w^2
          if (a == ZERO) then
            d = INFINITY
          else
            y0 = surf % coeffs(1)
            z0 = surf % coeffs(2)
            r = surf % coeffs(3)

            y = y - y0
            z = z - z0
            k = y*v + z*w
            c = y*y + z*z - r*r
            quad = k*k - a*c

            if (quad < ZERO) then
              ! no intersection with cylinder

              d = INFINITY

            elseif (on_surface) then
              ! particle is on the cylinder, thus one distance is
              ! positive/negative and the other is zero. The sign of k
              ! determines if we are facing in or out

              if (k >= ZERO) then
                d = INFINITY
              else
                d = (-k + sqrt(quad))/a
              end if

            elseif (c < ZERO) then
              ! particle is inside the cylinder, thus one distance must be
              ! negative and one must be positive. The positive distance
              ! will be the one with negative sign on sqrt(quad)

              d = (-k + sqrt(quad))/a

            else
              ! particle is outside the cylinder, thus both distances are
              ! either positive or negative. If positive, the smaller
              ! distance is the one with positive sign on sqrt(quad)

              d = (-k - sqrt(quad))/a
              if (d < ZERO) d = INFINITY

            end if
          end if

        case (SURF_CYL_Y)
          a = ONE - v*v  ! u^2 + w^2
          if (a == ZERO) then
            d = INFINITY
          else
            x0 = surf % coeffs(1)
            z0 = surf % coeffs(2)
            r = surf % coeffs(3)

            x = x - x0
            z = z - z0
            k = x*u + z*w
            c = x*x + z*z - r*r
            quad = k*k - a*c

            if (quad < ZERO) then
              ! no intersection with cylinder

              d = INFINITY

            elseif (on_surface) then
              ! particle is on the cylinder, thus one distance is
              ! positive/negative and the other is zero. The sign of k
              ! determines if we are facing in or out

              if (k >= ZERO) then
                d = INFINITY
              else
                d = (-k + sqrt(quad))/a
              end if

            elseif (c < ZERO) then
              ! particle is inside the cylinder, thus one distance must be
              ! negative and one must be positive. The positive distance
              ! will be the one with negative sign on sqrt(quad)

              d = (-k + sqrt(quad))/a

            else
              ! particle is outside the cylinder, thus both distances are
              ! either positive or negative. If positive, the smaller
              ! distance is the one with positive sign on sqrt(quad)

              d = (-k - sqrt(quad))/a
              if (d < ZERO) d = INFINITY

            end if
          end if

        case (SURF_CYL_Z)
          a = ONE - w*w  ! u^2 + v^2
          if (a == ZERO) then
            d = INFINITY
          else
            x0 = surf % coeffs(1)
            y0 = surf % coeffs(2)
            r = surf % coeffs(3)

            x = x - x0
            y = y - y0
            k = x*u + y*v
            c = x*x + y*y - r*r
            quad = k*k - a*c

            if (quad < ZERO) then
              ! no intersection with cylinder

              d = INFINITY

            elseif (on_surface) then
              ! particle is on the cylinder, thus one distance is
              ! positive/negative and the other is zero. The sign of k
              ! determines if we are facing in or out

              if (k >= ZERO) then
                d = INFINITY
              else
                d = (-k + sqrt(quad))/a
              end if

            elseif (c < ZERO) then
              ! particle is inside the cylinder, thus one distance must be
              ! negative and one must be positive. The positive distance
              ! will be the one with negative sign on sqrt(quad)

              d = (-k + sqrt(quad))/a

            else
              ! particle is outside the cylinder, thus both distances are
              ! either positive or negative. If positive, the smaller
              ! distance is the one with positive sign on sqrt(quad)

              d = (-k - sqrt(quad))/a
              if (d <= ZERO) d = INFINITY

            end if
          end if

        case (SURF_SPHERE)
          x0 = surf % coeffs(1)
          y0 = surf % coeffs(2)
          z0 = surf % coeffs(3)
          r = surf % coeffs(4)

          x = x - x0
          y = y - y0
          z = z - z0
          k = x*u + y*v + z*w
          c = x*x + y*y + z*z - r*r
          quad = k*k - c

          if (quad < ZERO) then
            ! no intersection with sphere

            d = INFINITY

          elseif (on_surface) then
            ! particle is on the sphere, thus one distance is
            ! positive/negative and the other is zero. The sign of k
            ! determines if we are facing in or out

            if (k >= ZERO) then
              d = INFINITY
            else
              d = -k + sqrt(quad)
            end if

          elseif (c < ZERO) then
            ! particle is inside the sphere, thus one distance must be
            ! negative and one must be positive. The positive distance will
            ! be the one with negative sign on sqrt(quad)

            d = -k + sqrt(quad)

          else
            ! particle is outside the sphere, thus both distances are either
            ! positive or negative. If positive, the smaller distance is the
            ! one with positive sign on sqrt(quad)

            d = -k - sqrt(quad)
            if (d < ZERO) d = INFINITY

          end if

        case (SURF_CONE_X)
          x0 = surf % coeffs(1)
          y0 = surf % coeffs(2)
          z0 = surf % coeffs(3)
          r = surf % coeffs(4)

          x = x - x0
          y = y - y0
          z = z - z0
          a = v*v + w*w - r*u*u
          k = y*v + z*w - r*x*u
          c = y*y + z*z - r*x*x
          quad = k*k - a*c

          if (quad < ZERO) then
            ! no intersection with cone

            d = INFINITY

          elseif (on_surface) then
            ! particle is on the cone, thus one distance is positive/negative
            ! and the other is zero. The sign of k determines which distance is
            ! zero and which is not.

            if (k >= ZERO) then
              d = (-k - sqrt(quad))/a
            else
              d = (-k + sqrt(quad))/a
            end if

          else
            ! calculate both solutions to the quadratic
            quad = sqrt(quad)
            d = (-k - quad)/a
            b = (-k + quad)/a

            ! determine the smallest positive solution
            if (d < ZERO) then
              if (b > ZERO) then
                d = b
              end if
            else
              if (b > ZERO) d = min(d, b)
            end if
          end if

          ! If the distance was negative, set boundary distance to infinity
          if (d <= ZERO) d = INFINITY

        case (SURF_CONE_Y)
          x0 = surf % coeffs(1)
          y0 = surf % coeffs(2)
          z0 = surf % coeffs(3)
          r = surf % coeffs(4)

          x = x - x0
          y = y - y0
          z = z - z0
          a = u*u + w*w - r*v*v
          k = x*u + z*w - r*y*v
          c = x*x + z*z - r*y*y
          quad = k*k - a*c

          if (quad < ZERO) then
            ! no intersection with cone

            d = INFINITY

          elseif (on_surface) then
            ! particle is on the cone, thus one distance is positive/negative
            ! and the other is zero. The sign of k determines which distance is
            ! zero and which is not.

            if (k >= ZERO) then
              d = (-k - sqrt(quad))/a
            else
              d = (-k + sqrt(quad))/a
            end if

          else
            ! calculate both solutions to the quadratic
            quad = sqrt(quad)
            d = (-k - quad)/a
            b = (-k + quad)/a

            ! determine the smallest positive solution
            if (d < ZERO) then
              if (b > ZERO) then
                d = b
              end if
            else
              if (b > ZERO) d = min(d, b)
            end if
          end if

          ! If the distance was negative, set boundary distance to infinity
          if (d <= ZERO) d = INFINITY

        case (SURF_CONE_Z)
          x0 = surf % coeffs(1)
          y0 = surf % coeffs(2)
          z0 = surf % coeffs(3)
          r = surf % coeffs(4)

          x = x - x0
          y = y - y0
          z = z - z0
          a = u*u + v*v - r*w*w
          k = x*u + y*v - r*z*w
          c = x*x + y*y - r*z*z
          quad = k*k - a*c

          if (quad < ZERO) then
            ! no intersection with cone

            d = INFINITY

          elseif (on_surface) then
            ! particle is on the cone, thus one distance is positive/negative
            ! and the other is zero. The sign of k determines which distance is
            ! zero and which is not.

            if (k >= ZERO) then
              d = (-k - sqrt(quad))/a
            else
              d = (-k + sqrt(quad))/a
            end if

          else
            ! calculate both solutions to the quadratic
            quad = sqrt(quad)
            d = (-k - quad)/a
            b = (-k + quad)/a

            ! determine the smallest positive solution
            if (d < ZERO) then
              if (b > ZERO) then
                d = b
              end if
            else
              if (b > ZERO) d = min(d, b)
            end if
          end if

          ! If the distance was negative, set boundary distance to infinity
          if (d <= ZERO) d = INFINITY

        end select

        ! Check is calculated distance is new minimum
        if (d < d_surf) then
          if (abs(d - d_surf)/d_surf >= FP_PRECISION) then
            d_surf = d
            level_surf_cross = -cl % surfaces(i)
          end if
        end if

      end do SURFACE_LOOP

      ! =======================================================================
      ! FIND MINIMUM DISTANCE TO LATTICE SURFACES

      LAT_COORD: if (coord % lattice /= NONE) then
        lat => lattices(coord % lattice) % obj

        LAT_TYPE: select type(lat)

        type is (RectLattice)
          ! copy local coordinates
          x = coord % xyz(1)
          y = coord % xyz(2)
          z = coord % xyz(3)

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
          z = coord % xyz(3)
          i_xyz(1) = coord % lattice_x
          i_xyz(2) = coord % lattice_y
          i_xyz(3) = coord % lattice_z
          parent_coord => p % coord0
          do while(.not. associated(parent_coord % next, coord))
            parent_coord => parent_coord % next
          end do

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
            xyz_t = lat % get_local_xyz(parent_coord % xyz, i_xyz+[1, 0, 0])
          else
            xyz_t = lat % get_local_xyz(parent_coord % xyz, i_xyz+[-1, 0, 0])
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
            xyz_t = lat % get_local_xyz(parent_coord % xyz, i_xyz+[1, -1, 0])
          else
            xyz_t = lat % get_local_xyz(parent_coord % xyz, i_xyz+[-1, 1, 0])
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
            xyz_t = lat % get_local_xyz(parent_coord % xyz, i_xyz+[0, 1, 0])
          else
            xyz_t = lat % get_local_xyz(parent_coord % xyz, i_xyz+[0, -1, 0])
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
          surface_crossed = level_surf_cross
          lattice_translation(:) = [0, 0, 0]
          final_coord => coord
        end if
      else
        if ((dist - d_lat)/dist >= FP_REL_PRECISION) then
          dist = d_lat
          surface_crossed = None
          lattice_translation(:) = level_lat_trans
          final_coord => coord
        end if
      end if

      coord => coord % next

    end do LEVEL_LOOP

    ! Move particle to appropriate coordinate level
    if (associated(final_coord)) p % coord => final_coord

  end subroutine distance_to_boundary

!===============================================================================
! SENSE determines whether a point is on the 'positive' or 'negative' side of a
! surface. This routine is crucial for determining what cell a particular point
! is in.
!===============================================================================

  recursive function sense(p, surf) result(s)

    type(Particle), intent(inout) :: p
    type(Surface),  pointer       :: surf   ! surface
    logical                       :: s      ! sense of particle

    real(8) :: x,y,z    ! coordinates of particle
    real(8) :: func     ! surface function evaluated at point
    real(8) :: A        ! coefficient on x for plane
    real(8) :: B        ! coefficient on y for plane
    real(8) :: C        ! coefficient on z for plane
    real(8) :: D        ! coefficient for plane
    real(8) :: x0,y0,z0 ! coefficients for quadratic surfaces / box
    real(8) :: r        ! radius for quadratic surfaces

    x = p % coord % xyz(1)
    y = p % coord % xyz(2)
    z = p % coord % xyz(3)

    select case (surf % type)
    case (SURF_PX)
      x0 = surf % coeffs(1)
      func = x - x0

    case (SURF_PY)
      y0 = surf % coeffs(1)
      func = y - y0

    case (SURF_PZ)
      z0 = surf % coeffs(1)
      func = z - z0

    case (SURF_PLANE)
      A = surf % coeffs(1)
      B = surf % coeffs(2)
      C = surf % coeffs(3)
      D = surf % coeffs(4)
      func = A*x + B*y + C*z - D

    case (SURF_CYL_X)
      y0 = surf % coeffs(1)
      z0 = surf % coeffs(2)
      r = surf % coeffs(3)
      y = y - y0
      z = z - z0
      func = y*y + z*z - r*r

    case (SURF_CYL_Y)
      x0 = surf % coeffs(1)
      z0 = surf % coeffs(2)
      r = surf % coeffs(3)
      x = x - x0
      z = z - z0
      func = x*x + z*z - r*r

    case (SURF_CYL_Z)
      x0 = surf % coeffs(1)
      y0 = surf % coeffs(2)
      r = surf % coeffs(3)
      x = x - x0
      y = y - y0
      func = x*x + y*y - r*r

    case (SURF_SPHERE)
      x0 = surf % coeffs(1)
      y0 = surf % coeffs(2)
      z0 = surf % coeffs(3)
      r = surf % coeffs(4)
      x = x - x0
      y = y - y0
      z = z - z0
      func = x*x + y*y + z*z - r*r

    case (SURF_CONE_X)
      x0 = surf % coeffs(1)
      y0 = surf % coeffs(2)
      z0 = surf % coeffs(3)
      r = surf % coeffs(4)
      x = x - x0
      y = y - y0
      z = z - z0
      func = y*y + z*z - r*x*x

    case (SURF_CONE_Y)
      x0 = surf % coeffs(1)
      y0 = surf % coeffs(2)
      z0 = surf % coeffs(3)
      r = surf % coeffs(4)
      x = x - x0
      y = y - y0
      z = z - z0
      func = x*x + z*z - r*y*y

    case (SURF_CONE_Z)
      x0 = surf % coeffs(1)
      y0 = surf % coeffs(2)
      z0 = surf % coeffs(3)
      r = surf % coeffs(4)
      x = x - x0
      y = y - y0
      z = z - z0
      func = x*x + y*y - r*z*z

    end select

    ! Check which side of surface the point is on
    if (abs(func) < FP_COINCIDENT) then
      ! Particle may be coincident with this surface. Artifically move the
      ! particle forward a tiny bit.
      p % coord % xyz = p % coord % xyz + TINY_BIT * p % coord % uvw
      s = sense(p, surf)
    elseif (func > 0) then
      s = .true.
    else
      s = .false.
    end if

  end function sense

!===============================================================================
! NEIGHBOR_LISTS builds a list of neighboring cells to each surface to speed up
! searches when a cell boundary is crossed.
!===============================================================================

  subroutine neighbor_lists()

    integer :: i          ! index in cells/surfaces array
    integer :: j          ! index of surface in cell
    integer :: i_surface  ! index in count arrays
    integer, allocatable :: count_positive(:) ! # of cells on positive side
    integer, allocatable :: count_negative(:) ! # of cells on negative side
    logical :: positive   ! positive side specified in surface list
    type(Cell),    pointer  :: c
    type(Surface), pointer  :: surf

    call write_message("Building neighboring cells lists for each surface...", &
         &4)

    allocate(count_positive(n_surfaces))
    allocate(count_negative(n_surfaces))
    count_positive = 0
    count_negative = 0

    do i = 1, n_cells
      c => cells(i)

      ! loop over each surface specification
      do j = 1, c % n_surfaces
        i_surface = c % surfaces(j)
        positive = (i_surface > 0)
        i_surface = abs(i_surface)
        if (positive) then
          count_positive(i_surface) = count_positive(i_surface) + 1
        else
          count_negative(i_surface) = count_negative(i_surface) + 1
        end if
      end do
    end do

    ! allocate neighbor lists for each surface
    do i = 1, n_surfaces
      surf => surfaces(i)
      if (count_positive(i) > 0) then
        allocate(surf%neighbor_pos(count_positive(i)))
      end if
      if (count_negative(i) > 0) then
        allocate(surf%neighbor_neg(count_negative(i)))
      end if
    end do

    count_positive = 0
    count_negative = 0

    ! loop over all cells
    do i = 1, n_cells
      c => cells(i)

      ! loop over each surface specification
      do j = 1, c % n_surfaces
        i_surface = c % surfaces(j)
        positive = (i_surface > 0)
        i_surface = abs(i_surface)

        surf => surfaces(i_surface)
        if (positive) then
          count_positive(i_surface) = count_positive(i_surface) + 1
          surf%neighbor_pos(count_positive(i_surface)) = i
        else
          count_negative(i_surface) = count_negative(i_surface) + 1
          surf%neighbor_neg(count_negative(i_surface)) = i
        end if
      end do
    end do

    deallocate(count_positive)
    deallocate(count_negative)

  end subroutine neighbor_lists

!===============================================================================
! HANDLE_LOST_PARTICLE
!===============================================================================

  subroutine handle_lost_particle(p, message)

    type(Particle), intent(inout) :: p
    character(*)                  :: message

    ! Print warning and write lost particle file
    call warning(message)
    call write_particle_restart(p)

    ! Increment number of lost particles
    p % alive = .false.
!$omp atomic
    n_lost_particles = n_lost_particles + 1

    ! Abort the simulation if the maximum number of lost particles has been
    ! reached
    if (n_lost_particles == MAX_LOST_PARTICLES) then
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


end module geometry
