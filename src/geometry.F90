module geometry

  use constants
  use error,                  only: fatal_error, warning
  use geometry_header,        only: Cell, Surface, Universe, Lattice
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
    type(Surface), pointer, save :: s => null()
!$omp threadprivate(s)

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

    ! If we've reached here, then the sense matched on every surface
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
    type(Cell),       pointer, save :: c => null()     ! pointer to cell
    type(Universe),   pointer, save :: univ => null()  ! universe to search in
    type(LocalCoord), pointer, save :: coord => null() ! particle coordinate to search on
!$omp threadprivate(c, univ, coord)

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
            message = "Overlapping cells detected: " //               &
                      trim(to_str(cells(index_cell) % id)) // ", " // &
                      trim(to_str(cells(coord % cell) % id)) //       &
                      " on universe " // trim(to_str(univ % id))
            call fatal_error()
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

    integer :: i                       ! index over cells
    integer :: i_xyz(3)                ! indices in lattice
    integer :: n                       ! number of cells to search
    integer :: index_cell              ! index in cells array
    real(8) :: xyz(3)                  ! temporary location
    logical :: use_search_cells        ! use cells provided as argument
    type(Cell),     pointer, save :: c => null()    ! pointer to cell
    type(Lattice),  pointer, save :: lat => null()  ! pointer to lattice
    type(Universe), pointer, save :: univ => null() ! universe to search in
!$omp threadprivate(c, lat, univ)

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
        message = "    Entering cell " // trim(to_str(c % id))
        call write_message()
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
        if (allocated(c % rotation)) then
          p % coord % xyz = matmul(c % rotation, p % coord % xyz)
          p % coord % uvw = matmul(c % rotation, p % coord % uvw)
          p % coord % rotated = .true.
        end if

        call find_cell(p, found)
        if (.not. found) exit

      elseif (c % type == CELL_LATTICE) then CELL_TYPE
        ! ======================================================================
        ! CELL CONTAINS LATTICE, RECURSIVELY FIND CELL

        ! Set current lattice
        lat => lattices(c % fill)

        ! Determine lattice indices
        i_xyz = get_lat_indices(lat, &
            &p % coord % xyz + TINY_BIT * p % coord % uvw, p % coord % uvw)

        ! Create new level of coordinates
        allocate(p % coord % next)
        p % coord % next % xyz = get_lat_trans(lat, p % coord % xyz, i_xyz)
        p % coord % next % uvw = p % coord % uvw

        ! set particle lattice indices
        p % coord % next% lattice   = c % fill
        p % coord % next% lattice_x = i_xyz(1)
        p % coord % next% lattice_y = i_xyz(2)
        p % coord % next% lattice_z = i_xyz(3)

        if (is_valid_lat_index(lat, i_xyz)) then
          p % coord % next % universe = &
              &lat % universes(i_xyz(1), i_xyz(2), i_xyz(3))

          ! Move particle to next level and search for the lower cells.
          p % coord => p % coord % next
          call find_cell(p, found)
          if (.not. found) exit

        else
          ! We're outside the lattice which means the next coordinates will
          ! have the same universe as the current coordinates, but we still
          ! need the lattice coordinates so distance_to_boundary will
          ! correctly handle the particle if it reenters the lattice.
          p % last_material = p % material
          p % material = lat % outside
          p % coord % next % universe = p % coord % universe
          p % coord % next % cell = index_cell

          ! Move particle to next level.
          p % coord => p % coord % next
        end if

      end if CELL_TYPE

      ! Found cell so we can return
      found = .true.
      return
    end do CELL_LOOP

    found = .false.

  end subroutine find_cell

!===============================================================================
! GET_LAT_TRANS returns the translated xyz coordinates in the specified lattice
! cell for a particle currently at the given xyz coordinates.
!===============================================================================

  function get_lat_trans(lat, xyz, i_xyz) result(xyz_t)
    type(Lattice), pointer, intent(in)  :: lat
    real(8)               , intent(in)  :: xyz(3)
    integer               , intent(in)  :: i_xyz(3)

    real(8) :: xyz_t(3)

    integer :: n_rings

    if (lat % type == LATTICE_RECT) then
      xyz_t(1) = xyz(1) - (lat % lower_left(1) + &
          &(i_xyz(1) - 0.5_8)*lat % width(1))
      xyz_t(2) = xyz(2) - (lat % lower_left(2) + &
          &(i_xyz(2) - 0.5_8)*lat % width(2))
      if (lat % n_dimension == 3) then
        xyz_t(3) = xyz(3) - (lat % lower_left(3) + &
            &(i_xyz(3) - 0.5_8)*lat % width(3))
      else
        xyz_t(3) = xyz(3)
      end if
 
    else if (lat % type == LATTICE_HEX) then
      n_rings = lat % dimension(1)
 
      xyz_t(1) = xyz(1) - (lat % lower_left(1) + &
          &sqrt(3.0_8) * (i_xyz(1) - n_rings) * lat % width(1))
      xyz_t(2) = xyz(2) - (lat % lower_left(2) + &
          &(2 * (i_xyz(2) - n_rings)) * lat % width(1) + &
          &(i_xyz(1) - n_rings) * lat % width(1))
      if (lat % n_dimension == 3) then
        xyz_t(3) = xyz(3) - (lat % lower_left(3) + &
            &(i_xyz(3) - 0.5_8) * lat % width(3))
      else
        xyz_t(3) = xyz(3)
      end if

    end if

  end function get_lat_trans

!===============================================================================
! IS_VALID_LAT_INDEX returns .true. if the given lattice indices fit within
! the bounds of the given lattice.  Returns false otherwise.
!===============================================================================

  function is_valid_lat_index(lat, i_xyz) result(is_valid)
    type(Lattice), pointer, intent(in)  :: lat
    integer               , intent(in)  :: i_xyz(3)
    logical                             :: is_valid

    integer  :: n_rings, n_x, n_y, n_z

    if (lat % type == LATTICE_RECT) then
      n_x = lat % dimension(1)
      n_y = lat % dimension(2)
      if (lat % n_dimension == 3) then
        n_z = lat % dimension(3)
      else
        n_z = 1
      end if
 
      is_valid = i_xyz(1) > 0 .and. i_xyz(1) <= n_x .and. &
                &i_xyz(2) > 0 .and. i_xyz(2) <= n_y .and. &
                &i_xyz(3) > 0 .and. i_xyz(3) <= n_z

    else if (lat % type == LATTICE_HEX) then
      n_rings = lat % dimension(1)

      if (lat % n_dimension == 2) then
        n_z = lat % dimension(2)
      else
        n_z = 1
      end if

      is_valid = (i_xyz(1) > 0 .and. i_xyz(1) < 2*n_rings .and. &
                 &i_xyz(2) > 0 .and. i_xyz(2) < 2*n_rings .and. &
                 &i_xyz(1) + i_xyz(2) > n_rings .and. &
                 &i_xyz(1) + i_xyz(2) < 3*n_rings .and. &
                 &i_xyz(3) > 0 .and. i_xyz(3) < n_z + 1)

    end if

  end function is_valid_lat_index

!===============================================================================
! GET_LAT_INDICES returns the indices in a lattice for the given xyz
! coordinates.
!===============================================================================

  function get_lat_indices(lat, xyz, uvw) result(i_xyz)
    type(Lattice), pointer, intent(in) :: lat
    real(8)               , intent(in) :: xyz(3)
    real(8)               , intent(in) :: uvw(3)

    integer :: i_xyz(3)

    real(8) :: alpha, beta_dir
    real(8) :: xyz_t(3)
    real(8) :: dists(4)
    integer :: n_rings
    integer :: i, j, k, loc(1)
 
    if (lat % type == LATTICE_RECT) then
      ! Find approximate indices using ceiling division.
      i_xyz(1) = ceiling((xyz(1) - lat % lower_left(1))/lat % width(1))
      i_xyz(2) = ceiling((xyz(2) - lat % lower_left(2))/lat % width(2))
      if (lat % n_dimension == 3) then
        i_xyz(3) = ceiling((xyz(3) - lat % lower_left(3))/lat % width(3))
      else
        i_xyz(3) = 1
      end if

      ! Adjust indices if particle is very close to a boundary and moving toward
      ! that boundary.
      !xyz_t = get_lat_trans(lat, xyz, i_xyz)
      !if (abs(xyz_t(1) - 0.5_8*lat % width(1)) &
      !    &< FP_REL_PRECISION*lat % width(1) .and. uvw(1) > 0.0) then
      !  i_xyz = i_xyz + (/1, 0, 0/)
      !else if (abs(xyz_t(1) + 0.5_8*lat % width(1)) &
      !    &< FP_REL_PRECISION*lat % width(1) .and. uvw(1) < 0.0) then
      !  i_xyz = i_xyz + (/-1, 0, 0/)
      !else if (abs(xyz_t(2) - 0.5_8*lat % width(2)) &
      !    &< FP_REL_PRECISION*lat % width(2) .and. uvw(2) > 0.0) then
      !  i_xyz = i_xyz + (/0, 1, 0/)
      !else if (abs(xyz_t(2) + 0.5_8*lat % width(2)) &
      !    &< FP_REL_PRECISION*lat % width(2) .and. uvw(2) < 0.0) then
      !  i_xyz = i_xyz + (/0, -1, 0/)
      !end if
      !if (lat % n_dimension == 3) then
      !  if (abs(xyz_t(3) - 0.5_8*lat % width(3)) &
      !      &< FP_REL_PRECISION*lat % width(3) .and. uvw(3) > 0.0) then
      !    i_xyz = i_xyz + (/0, 0, 1/)
      !  else if (abs(xyz_t(3) + 0.5_8*lat % width(3)) &
      !      &< FP_REL_PRECISION*lat % width(3) .and. uvw(3) < 0.0) then
      !    i_xyz = i_xyz + (/0, 0, -1/)
      !  end if
      !end if

    else if (lat % type == LATTICE_HEX) then
      ! Index z direction.
      if (lat % n_dimension == 2) then
        i_xyz(3) = ceiling((xyz(3) - lat % lower_left(3)) / lat % width(2))
      else
        i_xyz(3) = 1
      end if

      ! Convert coordinates into skewed bases.  The (x, alpha) basis is
      ! used to find the index of the particle coordinates to within 4
      ! cells.
      alpha = xyz(2) - xyz(1) / sqrt(3.0_8)
      i_xyz(1) = floor(xyz(1) / (sqrt(3.0_8) * lat % width(1)))
      i_xyz(2) = floor(alpha / (2.0_8 * lat % width(1)))

      ! Add offset to indices (the center cell is (i_x, i_alpha) = (0, 0)
      ! but the array is offset so that the indices never go below 1).
      n_rings = lat % dimension(1)
      i_xyz(1) = i_xyz(1) + n_rings
      i_xyz(2) = i_xyz(2) + n_rings

      ! Calculate the (squared) distance between the particle and the centers of
      ! the four possible cells.  Regular hexagonal tiles form a centroidal
      ! Voronoi tessellation so the particle should be in the hexagonal cell
      ! that it is closest to the center of.  This method is used over a
      ! remainder method (such as the rectangular one) becasue it provides
      ! better finite precision performance.  Squared distances are used becasue
      ! they are more computationally efficient than normal distances.
      k = 1
      do i=0,1
        do j=0,1
          xyz_t = get_lat_trans(lat, xyz, i_xyz + (/j, i, 0/))
          dists(k) = xyz_t(1)**2 + xyz_t(2)**2
          k = k + 1
        end do
      end do

      ! Now the lattice indices will be adjusted to one of the four possible
      ! lattice cells.  First we will check the special cases where the particle
      ! is in one cell but is very close to the cell that it is going towards.
      loc = minloc(dists)
      beta_dir = uvw(1)*sqrt(3.0_8)/2.0_8 + uvw(2)/2.0_8

      !! Check 1-2-3 corner.
      !if (abs(dists(1) - dists(2)) < FP_REL_PRECISION*dists(1) .and. &
      !    &abs(dists(3) - dists(1)) < FP_REL_PRECISION*dists(1)) then
      !  if (beta_dir > 0.0 .and. beta_dir > uvw(2)) then
      !    i_xyz = i_xyz + (/1, 0, 0/)
      !  else if (uvw(2) > 0.0) then
      !    i_xyz = i_xyz + (/0, 1, 0/)
      !  end if

      !! Check 2-3-4 corner.
      !else if (abs(dists(2) - dists(3)) < FP_REL_PRECISION*dists(2) .and. &
      !    &abs(dists(2) - dists(4)) < FP_REL_PRECISION*dists(2)) then
      !  if (beta_dir - uvw(2) > 0.0 .and. beta_dir - uvw(2) > beta_dir) then
      !    i_xyz = i_xyz + (/1, 0, 0/)
      !  else if (beta_dir > 0.0) then
      !    i_xyz = i_xyz + (/1, 1, 0/)
      !  else
      !    i_xyz = i_xyz + (/0, 1, 0/)
      !  end if

      !! Check 1-2 border.
      !else if (abs(dists(1) - dists(2)) < FP_REL_PRECISION*dists(1) .and. &
      !    &(loc(1) == 1 .or. loc(1) == 2)) then
      !  if (uvw(1)*sqrt(3.0_8)/2.0_8 + uvw(2)/2.0_8 > 0.0) then
      !    i_xyz = i_xyz + (/1, 0, 0/)
      !  end if

      !! Check 1-3 border.
      !else if (abs(dists(1) - dists(3)) < FP_REL_PRECISION*dists(1) .and. &
      !    &(loc(1) == 1 .or. loc(1) == 3)) then
      !  if (uvw(2) > 0.0) then
      !    i_xyz = i_xyz + (/0, 1, 0/)
      !  end if

      !! Check 2-3 border.
      !else if (abs(dists(2) - dists(3)) < FP_REL_PRECISION*dists(2) .and. &
      !    &(loc(1) == 2 .or. loc(1) == 3)) then
      !  if (uvw(1)*sqrt(3.0_8)/2.0_8 - uvw(2)/2.0_8 > 0.0) then
      !    i_xyz = i_xyz + (/1, 0, 0/)
      !  else
      !    i_xyz = i_xyz + (/0, 1, 0/)
      !  end if

      !! Check 2-4 border.
      !else if (abs(dists(2) - dists(4)) < FP_REL_PRECISION*dists(2) .and. &
      !    &(loc(1) == 2 .or. loc(1) == 4)) then
      !  if (uvw(2) > 0.0) then
      !    i_xyz = i_xyz + (/1, 1, 0/)
      !  else
      !    i_xyz = i_xyz + (/1, 0, 0/)
      !  end if

      !! Check 3-4 border.
      !else if (abs(dists(3) - dists(4)) < FP_REL_PRECISION*dists(3) .and. &
      !    &(loc(1) == 3 .or. loc(1) == 4)) then
      !  if (uvw(1)*sqrt(3.0_8)/2.0_8 + uvw(2)/2.0_8 > 0.0) then
      !    i_xyz = i_xyz + (/1, 1, 0/)
      !  else
      !    i_xyz = i_xyz + (/0, 1, 0/)
      !  end if

      ! The particle is not close to any of the lattice borders or corners.  It
      ! is in the cell that we would expect from the distance calculation.
      if (loc(1) == 2) then
        i_xyz = i_xyz + (/1, 0, 0/)

      else if (loc(1) == 3) then
        i_xyz = i_xyz + (/0, 1, 0/)

      else if (loc(1) == 4) then
        i_xyz = i_xyz + (/1, 1, 0/)

      end if
 
    end if

  end function get_lat_indices

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
    type(Surface), pointer, save :: surf => null()
!$omp threadprivate(surf)

    i_surface = abs(p % surface)
    surf => surfaces(i_surface)
    if (verbosity >= 10 .or. trace) then
      message = "    Crossing surface " // trim(to_str(surf % id))
      call write_message()
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
!$omp critical
        global_tallies(LEAKAGE) % value = &
           global_tallies(LEAKAGE) % value + p % wgt
!$omp end critical
      end if

      ! Display message
      if (verbosity >= 10 .or. trace) then
        message = "    Leaked out of surface " // trim(to_str(surf % id))
        call write_message()
      end if
      return

    elseif (surf % bc == BC_REFLECT .and. (run_mode /= MODE_PLOTTING)) then
      ! =======================================================================
      ! PARTICLE REFLECTS FROM SURFACE

      ! Do not handle reflective boundary conditions on lower universes
      if (.not. associated(p % coord, p % coord0)) then
        message = "Cannot reflect particle " // trim(to_str(p % id)) // &
             " off surface in a lower universe."
        call handle_lost_particle(p)
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
        message = "Reflection not supported for surface " // &
             trim(to_str(surf % id))
        call fatal_error()
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
          message = "Couldn't find particle after reflecting from surface."
          call handle_lost_particle(p)
          return
        end if
      end if

      ! Set previous coordinate going slightly past surface crossing
      p % last_xyz = p % coord0 % xyz + TINY_BIT * p % coord0 % uvw

      ! Diagnostic message
      if (verbosity >= 10 .or. trace) then
        message = "    Reflected from surface " // trim(to_str(surf%id))
        call write_message()
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
        message = "After particle " // trim(to_str(p % id)) // " crossed surface " &
             // trim(to_str(surfaces(i_surface) % id)) // " it could not be &
             &located in any cell and it did not leak."
        call handle_lost_particle(p)
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
    real(8) :: x0, y0, z0     ! half width of lattice element
    real(8) :: half_pitch     ! half pitch of a hexagonal lattice
    logical :: found          ! particle found in cell?
    type(Lattice), pointer, save :: lat => null()
    type(LocalCoord), pointer :: parent_coord
!$omp threadprivate(lat)

    lat => lattices(p % coord % lattice)

    if (verbosity >= 10 .or. trace) then
      message = "    Crossing lattice " // trim(to_str(lat % id)) // &
           ". Current position (" // trim(to_str(p % coord % lattice_x)) &
           // "," // trim(to_str(p % coord % lattice_y)) // "," // &
           trim(to_str(p % coord % lattice_z)) // "). " // "Translation (" // &
           trim(to_str(lattice_translation(1))) // "," // &
           trim(to_str(lattice_translation(2))) // "," // &
           trim(to_str(lattice_translation(3))) // ")"
      call write_message()
    end if

    ! Find the coordiante just above the current one.
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

    ! Set the coordinate position.
    p % coord % xyz = get_lat_trans(lat, parent_coord % xyz, i_xyz)

    OUTSIDE_LAT: if (.not. is_valid_lat_index(lat, i_xyz)) then
      ! The particle is outside the lattice.  Search for it from coord0.
      call deallocate_coord(p % coord0 % next)
      p % coord => p % coord0
      call find_cell(p, found)
      if (.not. found) then
        message = "Could not locate particle " // trim(to_str(p % id)) // &
             " after crossing a lattice boundary."
        call handle_lost_particle(p)
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
          message = "Could not locate particle " // trim(to_str(p % id)) // &
               " after crossing a lattice boundary."
          call handle_lost_particle(p)
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

    integer :: i                 ! index for surface in cell
    integer :: index_surf        ! index in surfaces array (with sign)
    real(8) :: x,y,z             ! particle coordinates
    real(8) :: alpha, beta, gama ! skewed particle coordiantes
    real(8) :: u,v,w             ! particle directions
    real(8) :: alpha_dir         ! skewed particle direction
    real(8) :: beta_dir          ! skewed particle direction
    real(8) :: gama_dir          ! skewed particle direction
    real(8) :: edge              ! distance to oncoming edge
    real(8) :: d                 ! evaluated distance
    real(8) :: x0,y0,z0          ! coefficients for surface
    real(8) :: r                 ! radius for quadratic surfaces
    real(8) :: tmp               ! dot product of surface normal with direction
    real(8) :: a,b,c,k           ! quadratic equation coefficients
    real(8) :: quad              ! discriminant of quadratic equation
    logical :: on_surface        ! is particle on surface?
    type(Cell),       pointer, save :: cl => null()
    type(Surface),    pointer, save :: surf => null()
    type(Lattice),    pointer, save :: lat => null()
    type(LocalCoord), pointer, save :: coord => null()
    type(LocalCoord), pointer, save :: final_coord => null()

    integer :: i_xyz(3)
    real(8) :: xyz_t(3)
    type(LocalCoord), pointer :: parent_coord
    real(8) :: d_lat
!$omp threadprivate(cl, surf, lat, coord, final_coord)

    ! inialize distance to infinity (huge)
    dist = INFINITY
    lattice_translation = (/0, 0, 0/)
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
            d = -(A*x + B*y + C*w - D)/tmp
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
        if (d < dist) then
          if (abs(d - dist)/dist >= FP_PRECISION) then
            dist = d
            surface_crossed = -cl % surfaces(i)
            lattice_translation = (/0, 0, 0/)
            final_coord => coord
          end if
        end if

      end do SURFACE_LOOP

      ! =======================================================================
      ! FIND MINIMUM DISTANCE TO LATTICE SURFACES

      LAT_COORD: if (coord % lattice /= NONE) then
        lat => lattices(coord % lattice)
        LAT_TYPE: if (lat % type == LATTICE_RECT) then
          ! copy local coordinates
          x = coord % xyz(1)
          y = coord % xyz(2)
          z = coord % xyz(3)

          ! determine oncoming edge
          x0 = sign(lat % width(1) * 0.5_8, u)
          y0 = sign(lat % width(2) * 0.5_8, v)

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
            lattice_translation = (/ 1, 0, 0 /)
          else
            lattice_translation = (/ -1, 0, 0 /)
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
              lattice_translation = (/ 0, 1, 0 /)
            else
              lattice_translation = (/ 0, -1, 0 /)
            end if
          end if

          if (lat % n_dimension == 3) then
            z0 = -sign(lat % width(3) * 0.5_8, w)

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
                lattice_translation = (/ 0, 0, 1 /)
              else
                lattice_translation = (/ 0, 0, -1 /)
              end if
            end if
          end if

        elseif (lat % type == LATTICE_HEX) then LAT_TYPE
          ! Copy local coordinates.
          x = coord % xyz(1)
          y = coord % xyz(2)
          z = coord % xyz(3)
          i_xyz(1) = coord % lattice_x
          i_xyz(2) = coord % lattice_y
          i_xyz(3) = coord % lattice_z
          parent_coord => p % coord0
          do while(.not. associated(parent_coord % next, coord))
            parent_coord => parent_coord % next
          end do

          ! Compute skewed coordinates and velocities.
          beta = y/2.0_8 + sqrt(3.0_8)*x/2.0_8
          beta_dir = v/2.0_8 + sqrt(3.0_8)*u/2.0_8
          gama = -y/2.0_8 + sqrt(3.0_8)*x/2.0_8
          gama_dir = -v/2.0_8 + sqrt(3.0_8)*u/2.0_8

          ! Upper right and lower left sides.
          edge = -sign(lat % width(1), beta_dir)  ! Oncoming edge
          if (beta_dir > 0.0) then
            xyz_t = get_lat_trans(lat, parent_coord % xyz, i_xyz + (/1, 0, 0/))
          else
            xyz_t = get_lat_trans(lat, parent_coord % xyz, i_xyz + (/-1, 0, 0/))
          end if
          beta = xyz_t(1)*sqrt(3.0_8)/2.0_8 + xyz_t(2)/2.0_8
          if (abs(beta - edge) < FP_PRECISION) then
            d = INFINITY
          else if (beta_dir == ZERO) then
            d = INFINITY
          else
            d = (edge - beta)/beta_dir
          end if

          d_lat = d
          if (beta_dir > 0) then
            lattice_translation = (/1, 0, 0/)
          else
            lattice_translation = (/-1, 0, 0/)
          end if

          ! Lower right and upper left sides.
          edge = -sign(lat % width(1), gama_dir)  ! Oncoming edge
          if (gama_dir > 0.0) then
            xyz_t = get_lat_trans(lat, parent_coord % xyz, i_xyz + (/1, -1, 0/))
          else
            xyz_t = get_lat_trans(lat, parent_coord % xyz, i_xyz + (/-1, 1, 0/))
          end if
          gama = xyz_t(1)*sqrt(3.0_8)/2.0_8 - xyz_t(2)/2.0_8
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
              lattice_translation = (/1, -1, 0/)
            else
              lattice_translation = (/-1, 1, 0/)
            end if
          end if

          ! Upper and lower sides.
          edge = -sign(lat % width(1), v)  ! Oncoming edge
          if (v > 0.0) then
            xyz_t = get_lat_trans(lat, parent_coord % xyz, i_xyz + (/0, 1, 0/))
          else
            xyz_t = get_lat_trans(lat, parent_coord % xyz, i_xyz + (/0, -1, 0/))
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
              lattice_translation = (/0, 1, 0/)
            else
              lattice_translation = (/0, -1, 0/)
            end if
          end if

          ! Top and bottom sides.
          if (lat % n_dimension == 2) then
            z0 = sign(lat % width(3) * 0.5_8, w)

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
                lattice_translation = (/ 0, 0, 1 /)
              else
                lattice_translation = (/ 0, 0, -1 /)
              end if
            end if
          end if
        end if LAT_TYPE

        ! If the lattice boundary is coincident with the parent cell boundary,
        ! we need to make sure that the lattice is not selected. This is
        ! complicated by the fact that floating point may determine that one
        ! is closer than the other (can't check direct equality). Thus, the
        ! logic here checks whether the relative difference is within floating
        ! point precision.
        if ((dist - d_lat)/dist >= FP_REL_PRECISION) then
          dist = d_lat
          final_coord => coord
        else
          lattice_translation = (/0, 0, 0/)
        end if

        if (d_lat < 0.0) then
          message = "Particle " // trim(to_str(p % id)) // &
              &" had a negative distance to a lattice boundary. d = " // &
              &trim(to_str(d_lat))
          call handle_lost_particle(p)
        end if
      end if LAT_COORD

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

    message = "Building neighboring cells lists for each surface..."
    call write_message(4)

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

  subroutine handle_lost_particle(p)

    type(Particle), intent(inout) :: p

    ! Print warning and write lost particle file
    call warning(force = .true.)
    call write_particle_restart(p)

    ! Increment number of lost particles
    p % alive = .false.
!$omp critical
    n_lost_particles = n_lost_particles + 1
!$omp end critical

    ! Abort the simulation if the maximum number of lost particles has been
    ! reached
    if (n_lost_particles == MAX_LOST_PARTICLES) then
      message = "Maximum number of lost particles has been reached."
      call fatal_error()
    end if

  end subroutine handle_lost_particle

end module geometry
