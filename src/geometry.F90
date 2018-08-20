module geometry

  use constants
  use error,                  only: fatal_error, warning, write_message
  use geometry_header
  use particle_header
  use simulation_header
  use settings
  use surface_header
  use stl_vector,             only: VectorInt
  use string,                 only: to_str

  use, intrinsic :: ISO_C_BINDING

  implicit none

  interface
    function cell_contains_c(cell_ptr, xyz, uvw, on_surface) &
         bind(C, name="cell_contains") result(in_cell)
      import C_PTR, C_DOUBLE, C_INT32_T, C_BOOL
      type(C_PTR),        intent(in), value :: cell_ptr
      real(C_DOUBLE),     intent(in)        :: xyz(3)
      real(C_DOUBLE),     intent(in)        :: uvw(3)
      integer(C_INT32_T), intent(in), value :: on_surface
      logical(C_BOOL)                       :: in_cell
    end function cell_contains_c

    function count_universe_instances(search_univ, target_univ_id) bind(C) &
         result(count)
      import C_INT32_T, C_INT
      integer(C_INT32_T), intent(in), value :: search_univ
      integer(C_INT32_T), intent(in), value :: target_univ_id
      integer(C_INT)                        :: count
    end function
  end interface

contains

  function cell_contains(c, p) result(in_cell)
    type(Cell), intent(in) :: c
    type(Particle), intent(in) :: p
    logical :: in_cell
    in_cell = cell_contains_c(c%ptr, p%coord(p%n_coord)%xyz, &
                              p%coord(p%n_coord)%uvw, p%surface)
  end function cell_contains

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
                 &// trim(to_str(cells(index_cell) % id())) // ", " &
                 &// trim(to_str(cells(p % coord(j) % cell) % id())) &
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
      call reset_coord(p % coord(j))
    end do
    j = p % n_coord

    ! Determine universe (if not yet set, use root universe)
    i_universe = p % coord(j) % universe
    if (i_universe == C_NONE) then
      p % coord(j) % universe = root_universe
      i_universe = root_universe
    end if

    ! set size of list to search
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
        if (cells(i_cell) % universe() /= i_universe) cycle
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
               cells(i_cell) % id())))
        end if

        found = .true.
        exit
      end if
    end do CELL_LOOP

    if (found) then
      associate(c => cells(i_cell))
        CELL_TYPE: if (c % type() == FILL_MATERIAL) then
          ! ======================================================================
          ! AT LOWEST UNIVERSE, TERMINATE SEARCH

          ! Save previous material and temperature
          p % last_material = p % material
          p % last_sqrtkT = p % sqrtkT

          ! Get distributed offset
          if (c % material_size() > 1 .or. size(c % sqrtkT) > 1) then
            ! Distributed instances of this cell have different
            ! materials/temperatures. Determine which instance this is for
            ! assigning the matching material/temperature.
            distribcell_index = c % distribcell_index
            offset = 0
            do k = 1, p % n_coord
              if (cells(p % coord(k) % cell) % type() == FILL_UNIVERSE) then
                offset = offset + cells(p % coord(k) % cell) &
                     % offset(distribcell_index-1)
              elseif (cells(p % coord(k) % cell) % type() == FILL_LATTICE) then
                if (lattices(p % coord(k + 1) % lattice) % obj &
                     % are_valid_indices([&
                     p % coord(k + 1) % lattice_x, &
                     p % coord(k + 1) % lattice_y, &
                     p % coord(k + 1) % lattice_z])) then
                  offset = offset + lattices(p % coord(k + 1) % lattice) % obj &
                       % offset(distribcell_index - 1, &
                       [p % coord(k + 1) % lattice_x, &
                       p % coord(k + 1) % lattice_y, &
                       p % coord(k + 1) % lattice_z])
                end if
              end if
            end do

            ! Keep track of which instance of the cell the particle is in
            p % cell_instance = offset + 1
          else
            p % cell_instance = 1
          end if

          ! Save the material
          if (c % material_size() > 1) then
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

        elseif (c % type() == FILL_UNIVERSE) then CELL_TYPE
          ! ======================================================================
          ! CELL CONTAINS LOWER UNIVERSE, RECURSIVELY FIND CELL

          ! Store lower level coordinates
          p % coord(j + 1) % xyz = p % coord(j) % xyz
          p % coord(j + 1) % uvw = p % coord(j) % uvw

          ! Move particle to next level and set universe
          j = j + 1
          p % n_coord = j
          p % coord(j) % universe = c % fill() + 1

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

        elseif (c % type() == FILL_LATTICE) then CELL_TYPE
          ! ======================================================================
          ! CELL CONTAINS LATTICE, RECURSIVELY FIND CELL

          associate (lat => lattices(c % fill() + 1) % obj)
            ! Determine lattice indices
            i_xyz = lat % get_indices(p % coord(j) % xyz + TINY_BIT * p % coord(j) % uvw)

            ! Store lower level coordinates
            p % coord(j + 1) % xyz = lat % get_local_xyz(p % coord(j) % xyz, i_xyz)
            p % coord(j + 1) % uvw = p % coord(j) % uvw

            ! set particle lattice indices
            p % coord(j + 1) % lattice   = c % fill() + 1
            p % coord(j + 1) % lattice_x = i_xyz(1)
            p % coord(j + 1) % lattice_y = i_xyz(2)
            p % coord(j + 1) % lattice_z = i_xyz(3)

            ! Set the next lowest coordinate level.
            if (lat % are_valid_indices(i_xyz)) then
              ! Particle is inside the lattice.
              p % coord(j + 1) % universe = &
                   lat % get([i_xyz(1), i_xyz(2), i_xyz(3)]) + 1

            else
              ! Particle is outside the lattice.
              if (lat % outer() == NO_OUTER_UNIVERSE) then
                call warning("Particle " // trim(to_str(p %id)) &
                     // " is outside lattice " // trim(to_str(lat % id())) &
                     // " but the lattice has no defined outer universe.")
                found = .false.
                return
              else
                p % coord(j + 1) % universe = lat % outer() + 1
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
      call write_message("    Crossing lattice " // trim(to_str(lat % id())) &
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
          call particle_mark_as_lost(p, "Could not locate particle " &
               // trim(to_str(p % id)) // " after crossing a lattice boundary.")
          return
        end if
      end if

    else OUTSIDE_LAT

      ! Find cell in next lattice element
      p % coord(j) % universe = &
           lat % get([i_xyz(1), i_xyz(2), i_xyz(3)]) + 1

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
          call particle_mark_as_lost(p, "Could not locate particle " // &
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

    integer :: j
    integer :: i_xyz(3)           ! lattice indices
    integer :: level_surf_cross   ! surface crossed on current level
    integer :: level_lat_trans(3) ! lattice translation on current level
    real(8) :: xyz_t(3)           ! local particle coordinates
    real(8) :: d_lat              ! distance to lattice boundary
    real(8) :: d_surf             ! distance to surface
    real(8) :: xyz_cross(3)       ! coordinates at projected surface crossing
    real(8) :: surf_uvw(3)        ! surface normal direction
    type(Cell),       pointer :: c
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

      ! =======================================================================
      ! FIND MINIMUM DISTANCE TO SURFACE IN THIS CELL

      call c % distance(p % coord(j) % xyz, p % coord(j) % uvw, p % surface, &
                        d_surf, level_surf_cross)

      ! =======================================================================
      ! FIND MINIMUM DISTANCE TO LATTICE SURFACES

      LAT_COORD: if (p % coord(j) % lattice /= NONE) then
        lat => lattices(p % coord(j) % lattice) % obj

        i_xyz(1) = p % coord(j) % lattice_x
        i_xyz(2) = p % coord(j) % lattice_y
        i_xyz(3) = p % coord(j) % lattice_z

        LAT_TYPE: select type(lat)

        type is (RectLattice)
          call lat % distance(p % coord(j) % xyz, p % coord(j) % uvw, &
                              i_xyz, d_lat, level_lat_trans)

        type is (HexLattice) LAT_TYPE
          xyz_t(1) = p % coord(j-1) % xyz(1)
          xyz_t(2) = p % coord(j-1) % xyz(2)
          xyz_t(3) = p % coord(j) % xyz(3)
          call lat % distance(xyz_t, p % coord(j) % uvw, &
                              i_xyz, d_lat, level_lat_trans)
        end select LAT_TYPE

        if (d_lat < ZERO) then
          call particle_mark_as_lost(p, "Particle " // trim(to_str(p % id)) &
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
          if (.not. c % simple()) then
            xyz_cross(:) = p % coord(j) % xyz + d_surf*p % coord(j) % uvw
            call surfaces(abs(level_surf_cross)) % normal(xyz_cross, surf_uvw)
            if (dot_product(p % coord(j) % uvw, surf_uvw) > ZERO) then
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
        allocate(surfaces(i)%neighbor_pos(n))
        surfaces(i)%neighbor_pos(:) = neighbor_pos(i)%data(1:n)
      end if

      ! Copy negative neighbors to Surface instance
      n = neighbor_neg(i)%size()
      if (n > 0) then
        allocate(surfaces(i)%neighbor_neg(n))
        surfaces(i)%neighbor_neg(:) = neighbor_neg(i)%data(1:n)
      end if
    end do

  end subroutine neighbor_lists

end module geometry
