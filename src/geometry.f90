module geometry

  use constants
  use datatypes,       only: dict_get_key
  use error,           only: fatal_error
  use geometry_header, only: Cell, Surface, Universe, Lattice
  use global
  use output,          only: write_message
  use particle_header, only: Particle, LocalCoord, deallocate_coord
  use string,          only: int_to_str
  use tally,           only: score_surface_current

  implicit none
     
contains

!===============================================================================
! CELL_CONTAINS determines whether a given point is inside a cell
!===============================================================================

  function cell_contains(c, p) result(in_cell)

    type(Cell),     pointer :: c
    type(Particle), pointer :: p
    logical                 :: in_cell

    integer, allocatable :: expression(:) ! copy of surfaces list
    integer :: specified_sense ! specified sense of surface in list
    integer :: actual_sense    ! sense of particle wrt surface
    integer :: n_surf          ! number of surfaces in cell
    integer :: i               ! index of surfaces in cell
    integer :: surf_num        ! index in surfaces array (with sign)
    integer :: current_surface ! current surface of particle (with sign)
    type(Surface), pointer  :: surf => null()

    current_surface = p % surface

    n_surf = size(c % surfaces)
    allocate(expression(n_surf))
    expression = c % surfaces
    do i = 1, n_surf

       ! Don't change logical operator
       if (expression(i) >= OP_DIFFERENCE) then
          cycle
       end if

       ! Lookup surface
       surf_num = expression(i)
       surf => surfaces(abs(surf_num))

       ! Check if the particle is currently on the specified surface
       if (surf_num == current_surface) then
          ! particle is on specified surface heading into cell
          expression(i) = 1
          cycle
       elseif (surf_num == -current_surface) then
          ! particle is on specified surface, but heading other direction
          expression(i) = 0
          cycle
       end if

       ! Compare sense of point to specified sense
       specified_sense = sign(1,expression(i))
       actual_sense = sense(surf, p % coord % xyz)
       if (actual_sense == specified_sense) then
          expression(i) = 1
       else
          expression(i) = 0
       end if

    end do

    ! TODO: Need to replace this with a 'lgeval' like subroutine that can
    ! actually test expressions with unions and parentheses
    if (all(expression == 1)) then
       in_cell = .true.
    else
       in_cell = .false.
    end if

    ! Free up memory from expression
    deallocate(expression)

  end function cell_contains

!===============================================================================
! FIND_CELL determines what cell a source particle is in within a particular
! universe. If the base universe is passed, the particle should be found as long
! as it's within the geometry
!===============================================================================

  recursive subroutine find_cell(p, found, search_cells)

    type(Particle), pointer  :: p
    logical, intent(inout)   :: found
    integer, optional        :: search_cells(:)

    integer :: i                    ! index over cells
    integer :: x                    ! x-index for lattice
    integer :: y                    ! y-index for lattice
    integer :: n                    ! number of cells to search
    integer :: index_cell           ! index in cells array
    real(8) :: xyz(3)               ! temporary location
    logical :: use_search_cells     ! use cells provided as argument
    type(Cell),     pointer :: c    ! pointer to cell
    type(Lattice),  pointer :: lat  ! pointer to lattice
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

    do i = 1, n
       ! select cells based on whether we are searching a universe or a provided
       ! list of cells (this would be for lists of neighbor cells)
       if (use_search_cells) then
          index_cell = search_cells(i) 
       else
          index_cell = univ % cells(i)
       end if

       ! get pointer to cell
       c => cells(index_cell)

       if (cell_contains(c, p)) then
          ! Set cell on this level
          p % coord % cell = index_cell

          if (c % type == CELL_NORMAL) then
             ! =================================================================
             ! AT LOWEST UNIVERSE, TERMINATE SEARCH

             ! set material
             p % material = c % material

          elseif (c % type == CELL_FILL) then
             ! =================================================================
             ! CELL CONTAINS LOWER UNIVERSE, RECURSIVELY FIND CELL

             ! Create new level of coordinates
             p % in_lower_universe = .true.
             allocate(p % coord % next)

             ! Move particle to next level and set universe
             p % coord => p % coord % next
             p % coord % universe = c % fill

             call find_cell(p, found)
             if (.not. found) exit

          elseif (c % type == CELL_LATTICE) then
             ! =================================================================
             ! CELL CONTAINS LATTICE, RECURSIVELY FIND CELL

             ! Set current lattice
             lat => lattices(c % fill)

             ! determine universe based on lattice position
             xyz = p % coord % xyz + TINY_BIT * p % coord % uvw
             x = ceiling((xyz(1) - lat % x0)/lat % width_x)
             y = ceiling((xyz(2) - lat % y0)/lat % width_y)

             ! Create new level of coordinates
             p % in_lower_universe = .true.
             allocate(p % coord % next)
             
             ! adjust local position of particle
             p % coord % next % xyz(1) = p % coord % xyz(1) - &
                  (lat%x0 + (x-0.5_8)*lat%width_x)
             p % coord % next % xyz(2) = p % coord % xyz(2) - &
                  (lat%y0 + (y-0.5_8)*lat%width_y)
             p % coord % next % xyz(3) = p % coord % xyz(3)
             p % coord % next % uvw    = p % coord % uvw

             ! Move particle to next level
             p % coord => p % coord % next
             
             ! set particle lattice indices
             p % coord % lattice   = c % fill
             p % coord % lattice_x = x
             p % coord % lattice_y = y
             p % coord % universe  = lat % element(x,y)

             call find_cell(p, found)
             if (.not. found) exit
          end if
          
          ! Found cell so we can return
          found = .true.
          return
       end if
    end do

    found = .false.

  end subroutine find_cell

!===============================================================================
! CROSS_SURFACE handles all surface crossings, whether the particle leaks out of
! the geometry, is reflected, or crosses into a new lattice or cell
!===============================================================================

  subroutine cross_surface(p, last_cell)

    type(Particle), pointer :: p
    integer, intent(in)     :: last_cell  ! last cell particle was in

    real(8)                 :: x          ! x-x0 for sphere
    real(8)                 :: y          ! y-y0 for sphere
    real(8)                 :: z          ! z-z0 for sphere
    real(8)                 :: R          ! radius of sphere
    real(8)                 :: u          ! x-component of direction
    real(8)                 :: v          ! y-component of direction
    real(8)                 :: w          ! z-component of direction
    real(8)                 :: n1         ! x-component of surface normal
    real(8)                 :: n2         ! y-component of surface normal
    real(8)                 :: n3         ! z-component of surface normal
    real(8)                 :: dot_prod   ! dot product of direction and normal
    real(8)                 :: norm       ! "norm" of surface normal
    logical                 :: found      ! particle found in universe?
    type(Surface),  pointer :: surf => null()

    surf => surfaces(abs(p % surface))
    if (verbosity >= 10 .or. trace) then
       message = "    Crossing surface " // trim(int_to_str(surf % id))
       call write_message()
    end if

    if (surf % bc == BC_VACUUM .and. (.not. plotting)) then
       ! =======================================================================
       ! PARTICLE LEAKS OUT OF PROBLEM

       ! Kill particle
       p % alive = .false.

       ! Score any surface current tallies -- note that the particle is moved
       ! forward slightly so that if the mesh boundary is on the surface, it is
       ! still processed

       ! TODO: Find a better solution to score surface currents than physically
       ! moving the particle forward slightly

       p % coord0 % xyz = p % coord0 % xyz + TINY_BIT * p % coord0 % uvw
       call score_surface_current(p)

       ! Display message
       if (verbosity >= 10 .or. trace) then
          message = "    Leaked out of surface " // trim(int_to_str(surf % id))
          call write_message()
       end if
       return

    elseif (surf % bc == BC_REFLECT .and. (.not. plotting)) then
       ! =======================================================================
       ! PARTICLE REFLECTS FROM SURFACE

       ! Do not handle reflective boundary conditions on lower universes
       if (p % in_lower_universe) then
          message = "Cannot reflect particle off surface in a lower universe."
          call fatal_error()
       end if

       ! Score surface currents since reflection causes the direction of the
       ! particle to change -- artificially move the particle slightly back in
       ! case the surface crossing in coincident with a mesh boundary
       p % coord0 % xyz = p % coord0 % xyz - TINY_BIT * p % coord0 % uvw
       call score_surface_current(p)
       p % coord0 % xyz = p % coord0 % xyz + TINY_BIT * p % coord0 % uvw

       ! Copy particle's direction cosines
       u = p % coord0 % uvw(1)
       v = p % coord0 % uvw(2)
       w = p % coord0 % uvw(3)

       select case (surf%type)
       case (SURF_PX)
          p % coord0 % uvw = (/ -u, v, w /)
       case (SURF_PY)
          p % coord0 % uvw = (/ u, -v, w /)
       case (SURF_PZ)
          p % coord0 % uvw = (/ u, v, -w /)
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

          ! Set vector
          p % coord0 % uvw = (/ u, v, w /)
       case (SURF_CYL_X)
          ! Find y-y0, z-z0 and dot product of direction and surface normal
          y = p % coord0 % xyz(2) - surf % coeffs(1)
          z = p % coord0 % xyz(3) - surf % coeffs(2)
          R = surf % coeffs(3)
          dot_prod = v*y + w*z

          ! Reflect direction according to normal
          v = v - 2*dot_prod*y/(R*R)
          w = w - 2*dot_prod*z/(R*R)

          ! Set vector
          p % coord0 % uvw = (/ u, v, w /)
       case (SURF_CYL_Y)
          ! Find x-x0, z-z0 and dot product of direction and surface normal
          x = p % coord0 % xyz(1) - surf % coeffs(1)
          z = p % coord0 % xyz(3) - surf % coeffs(2)
          R = surf % coeffs(3)
          dot_prod = u*x + w*z

          ! Reflect direction according to normal
          u = u - 2*dot_prod*x/(R*R)
          w = w - 2*dot_prod*z/(R*R)

          ! Set vector
          p % coord0 % uvw = (/ u, v, w /)
       case (SURF_CYL_Z)
          ! Find x-x0, y-y0 and dot product of direction and surface normal
          x = p % coord0 % xyz(1) - surf % coeffs(1)
          y = p % coord0 % xyz(2) - surf % coeffs(2)
          R = surf % coeffs(3)
          dot_prod = u*x + v*y

          ! Reflect direction according to normal
          u = u - 2*dot_prod*x/(R*R)
          v = v - 2*dot_prod*y/(R*R)

          ! Set vector
          p % coord0 % uvw = (/ u, v, w /)
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

          ! Set vector
          p % coord0 % uvw = (/ u, v, w /)
       case default
          message = "Reflection not supported for surface " // &
               trim(int_to_str(surf % id))
          call fatal_error()
       end select

       ! Reassign particle's cell and surface
       p % coord0 % cell = last_cell
       p % surface = -p % surface

       ! Set previous coordinate going slightly past surface crossing
       p % last_xyz = p % coord0 % xyz + TINY_BIT * p % coord0 % uvw

       ! Diagnostic message
       if (verbosity >= 10 .or. trace) then
          message = "    Reflected from surface " // trim(int_to_str(surf%id))
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

    call find_cell(p, found)

    ! Couldn't find next cell anywhere!
    if ((.not. found) .and. (.not. plotting)) then
       message = "After particle crossed surface " // trim(int_to_str(p%surface)) &
            // ", it could not be located in any cell and it did not leak."
       call fatal_error()
    end if
       
  end subroutine cross_surface

!===============================================================================
! CROSS_LATTICE moves a particle into a new lattice element
!===============================================================================

  subroutine cross_lattice(p)

    type(Particle), pointer :: p

    integer        :: i_x      ! x index in lattice
    integer        :: i_y      ! y index in lattice
    real(8)        :: d_left   ! distance to left side
    real(8)        :: d_right  ! distance to right side
    real(8)        :: d_bottom ! distance to bottom side
    real(8)        :: d_top    ! distance to top side
    real(8)        :: dist     ! shortest distance
    real(8)        :: x        ! x coordinate in local lattice element
    real(8)        :: y        ! y coordinate in local lattice element
    real(8)        :: z        ! z coordinate in local lattice element
    real(8)        :: u        ! cosine of angle with x axis
    real(8)        :: v        ! cosine of angle with y axis
    real(8)        :: x0       ! half the width of lattice element
    real(8)        :: y0       ! half the height of lattice element
    logical        :: found    ! particle found in cell?
    type(Lattice),  pointer :: lat => null()

    lat => lattices(p % coord % lattice)

    if (verbosity >= 10 .or. trace) then
       message = "    Crossing lattice " // int_to_str(lat % id)
       call write_message()
    end if

    u = p % coord % uvw(1)
    v = p % coord % uvw(2)

    if (lat % type == LATTICE_RECT) then
       x = p % coord % xyz(1)
       y = p % coord % xyz(2)
       z = p % coord % xyz(3)
       x0 = lat % width_x * 0.5_8
       y0 = lat % width_y * 0.5_8
       
       dist = INFINITY

       ! left and right sides
       if (u == ZERO) then
          d_left = INFINITY
          d_right = INFINITY
       elseif (u > 0) then
          d_left = INFINITY
          d_right = (x0 - x)/u
       else
          d_left = -(x0 + x)/u
          d_right = INFINITY
       end if

       ! top and bottom sides
       if (v == ZERO) then
          d_bottom = INFINITY
          d_top = INFINITY
       elseif (v > 0) then
          d_bottom = INFINITY
          d_top = (y0 - y)/v
       else
          d_bottom = -(y0 + y)/v
          d_top = INFINITY
       end if

       dist = min(d_left, d_right, d_top, d_bottom)
       if (dist == d_left) then
          ! Move particle to left element
          p % coord % lattice_x = p % coord % lattice_x - 1
          p % coord % xyz(1) = x0
          
       elseif (dist == d_right) then
          ! Move particle to right element
          p % coord % lattice_x = p % coord % lattice_x + 1
          p % coord % xyz(1) = -x0

       elseif (dist == d_bottom) then
          ! Move particle to bottom element
          p % coord % lattice_y = p % coord % lattice_y - 1
          p % coord % xyz(2) = y0

       elseif (dist == d_top) then
          ! Move particle to top element
          p % coord % lattice_y = p % coord % lattice_y + 1
          p % coord % xyz(2) = -y0

       end if
    elseif (lat % type == LATTICE_HEX) then
       ! TODO: Add hex lattice support
    end if
    
    ! Check to make sure still in lattice
    i_x = p % coord % lattice_x
    i_y = p % coord % lattice_y
    if (i_x < 1 .or. i_x > lat % n_x) then
       message = "Reached edge of lattice " // trim(int_to_str(lat % id)) // &
            " at position (" // trim(int_to_str(i_x)) // "," // &
            trim(int_to_str(i_y)) // ")."
       call fatal_error()
    elseif (i_y < 1 .or. i_y > lat % n_y) then
       message = "Reached edge of lattice " // trim(int_to_str(lat % id)) // &
            " at position (" // trim(int_to_str(i_x)) // "," // &
            trim(int_to_str(i_y)) // ")."
       call fatal_error()
    end if

    ! Find universe for next lattice element
    p % coord % universe = lat % element(i_x, i_y)

    ! Find cell in next lattice element
    call find_cell(p, found)
    if (.not. found) then
       message = "Could not locate particle in universe: " // &
            int_to_str(universes(p % coord % universe) % id)
       call fatal_error()
    end if

  end subroutine cross_lattice

!===============================================================================
! DISTANCE_TO_BOUNDARY calculates the distance to the nearest boundary for a
! particle 'p' traveling in a certain direction. For a cell in a subuniverse
! that has a parent cell, also include the surfaces of the edge of the universe.
!===============================================================================

  subroutine distance_to_boundary(p, dist, surface_crossed, lattice_crossed)

    type(Particle), pointer     :: p
    real(8),        intent(out) :: dist
    integer,        intent(out) :: surface_crossed
    logical,        intent(out) :: lattice_crossed

    integer :: i            ! index for surface in cell
    integer :: index_surf   ! index in surfaces array (with sign)
    real(8) :: x,y,z        ! particle coordinates
    real(8) :: u,v,w        ! particle directions
    real(8) :: d            ! evaluated distance
    real(8) :: x0,y0,z0     ! coefficients for surface
    real(8) :: r            ! radius for quadratic surfaces
    real(8) :: tmp          ! dot product of surface normal with direction
    real(8) :: a,b,c,k      ! quadratic equation coefficients
    real(8) :: quad         ! discriminant of quadratic equation
    logical :: on_surface   ! is particle on surface?
    type(Cell),       pointer :: cl => null()
    type(Surface),    pointer :: surf => null()
    type(Lattice),    pointer :: lat => null()
    type(LocalCoord), pointer :: coord => null()
    type(LocalCoord), pointer :: final_coord => null()

    ! inialize distance to infinity (huge)
    dist = INFINITY
    lattice_crossed = .false.
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

          case (SURF_GQ)
             message = "Surface distance not yet implement for general quadratic."
             call fatal_error()

          end select

          ! Check is calculated distance is new minimum
          if (d < dist) then
             dist = d
             surface_crossed = -cl % surfaces(i)
             lattice_crossed = .false.
             final_coord => coord
          end if

       end do SURFACE_LOOP

       ! =======================================================================
       ! FIND MINIMUM DISTANCE TO LATTICE SURFACES

       if (coord % lattice /= NONE) then
          lat => lattices(coord % lattice)
          if (lat % type == LATTICE_RECT) then
             ! copy local coordinates
             x = coord % xyz(1)
             y = coord % xyz(2)
             z = coord % xyz(3)

             ! determine oncoming edge
             x0 = lat % width_x * 0.5_8
             y0 = lat % width_y * 0.5_8

             ! left and right sides
             if (u == ZERO) then
                d = INFINITY
             elseif (u > 0) then
                d = (x0 - x)/u
             else
                d = -(x + x0)/u
             end if

             ! If the lattice boundary is coincident with the parent cell boundary,
             ! we need to make sure that the lattice is not selected. This is
             ! complicated by the fact that floating point may determine that one
             ! is closer than the other (can't check direct equality). Thus, the
             ! logic here checks whether the relative difference is within floating
             ! point precision.

             if (d < dist) then 
                if (abs(d - dist)/dist >= FP_PRECISION) then
                   dist = d
                   lattice_crossed = .true.
                   final_coord => coord
                end if
             end if

             ! top and bottom sides
             if (v == ZERO) then
                d = INFINITY
             elseif (v > 0) then
                d = (y0 - y)/v
             else
                d = -(y + y0)/v
             end if

             if (d < dist) then
                if (abs(d - dist)/dist >= FP_PRECISION) then
                   dist = d
                   lattice_crossed = .true.
                   final_coord => coord
                end if
             end if

          elseif (lat % type == LATTICE_HEX) then
             ! TODO: Add hex lattice support
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

  function sense(surf, xyz) result(s)

    type(Surface), pointer    :: surf   ! surface
    real(8),       intent(in) :: xyz(3) ! coordinates of particle
    integer                   :: s      ! sense of particle

    real(8) :: x,y,z    ! coordinates of particle
    real(8) :: func     ! surface function evaluated at point
    real(8) :: A        ! coefficient on x**2 term in GQ
    real(8) :: B        ! coefficient on y**2 term in GQ
    real(8) :: C        ! coefficient on z**2 term in GQ
    real(8) :: D        ! coefficient on x*y term in GQ
    real(8) :: E        ! coefficient on y*z term in GQ
    real(8) :: F        ! coefficient on x*z term in GQ
    real(8) :: G        ! coefficient on x term in GQ
    real(8) :: H        ! coefficient on y term in GQ
    real(8) :: I        ! coefficient on z term in GQ
    real(8) :: J        ! coefficient on constant term in GQ
    real(8) :: x0,y0,z0 ! coefficients for quadratic surfaces / box
    real(8) :: r        ! radius for quadratic surfaces
    real(8) :: x1,y1,z1 ! upper-right corner of box

    x = xyz(1)
    y = xyz(2)
    z = xyz(3)

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

    case (SURF_BOX_X)
       y0 = surf % coeffs(1)
       z0 = surf % coeffs(2)
       y1 = surf % coeffs(3)
       z1 = surf % coeffs(4)
       if (y >= y0 .and. y < y1 .and. z >= z0 .and. z < z1) then
          s = SENSE_NEGATIVE
       else
          s = SENSE_POSITIVE
       end if
       return

    case (SURF_BOX_Y)
       x0 = surf % coeffs(1)
       z0 = surf % coeffs(2)
       x1 = surf % coeffs(3)
       z1 = surf % coeffs(4)
       if (x >= x0 .and. x < x1 .and. z >= z0 .and. z < z1) then
          s = SENSE_NEGATIVE
       else
          s = SENSE_POSITIVE
       end if
       return

    case (SURF_BOX_Z)
       x0 = surf % coeffs(1)
       y0 = surf % coeffs(2)
       x1 = surf % coeffs(3)
       y1 = surf % coeffs(4)
       if (x >= x0 .and. x < x1 .and. y >= y0 .and. y < y1) then
          s = SENSE_NEGATIVE
       else
          s = SENSE_POSITIVE
       end if
       return

    case (SURF_BOX)
       x0 = surf % coeffs(1)
       y0 = surf % coeffs(2)
       z0 = surf % coeffs(3)
       x1 = surf % coeffs(4)
       y1 = surf % coeffs(5)
       z1 = surf % coeffs(6)
       if (x >= x0 .and. x < x1 .and. y >= y0 .and. y < y1 .and. & 
            & z >= z0 .and. z < z1) then
          s = SENSE_NEGATIVE
       else
          s = SENSE_POSITIVE
       end if
       return

    case (SURF_GQ)
       A = surf % coeffs(1)
       B = surf % coeffs(2)
       C = surf % coeffs(3)
       D = surf % coeffs(4)
       E = surf % coeffs(5)
       F = surf % coeffs(6)
       G = surf % coeffs(7)
       H = surf % coeffs(8)
       I = surf % coeffs(9)
       J = surf % coeffs(10)
       func = A*x*x + B*y*y + C*z*z + D*x*y + E*y*z + F*x*z + G*x &
            & + H*y + I*z + J

    end select

    ! Check which side of surface the point is on
    if (func > 0) then
       s = SENSE_POSITIVE
    else
       s = SENSE_NEGATIVE
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

end module geometry
