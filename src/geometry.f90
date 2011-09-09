module geometry

  use constants
  use datatypes,       only: dict_get_key
  use error,           only: fatal_error
  use geometry_header, only: Cell, Surface, Universe, Lattice
  use global
  use output,          only: message
  use particle_header, only: Particle
  use string,          only: int_to_str

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
    integer :: n_surfaces      ! number of surfaces in cell
    integer :: i               ! index of surfaces in cell
    integer :: surf_num        ! index in surfaces array (with sign)
    integer :: current_surface ! current surface of particle (with sign)
    type(Surface), pointer  :: surf => null()

    current_surface = p%surface

    n_surfaces = size(c%surfaces)
    allocate(expression(n_surfaces))
    expression = c%surfaces
    do i = 1, n_surfaces

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
       actual_sense = sense(surf, p%xyz_local)
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

  recursive subroutine find_cell(univ, p, found)

    type(Universe), pointer :: univ       ! universe to search in
    type(Particle), pointer :: p          ! pointer to particle
    logical,  intent(inout) :: found      ! particle found?

    character(MAX_LINE_LEN) :: msg        ! error message
    integer                 :: i          ! index over cells
    integer                 :: x, y
    type(Cell),     pointer :: c          ! pointer to cell
    type(Lattice),  pointer :: lat        ! pointer to lattice
    type(Universe), pointer :: lower_univ ! if particle is in lower
                                          ! universe, use this pointer
                                          ! to call recursively

    found = .false.

    ! determine what region in
    do i = 1, univ % n_cells
       c => cells(univ % cells(i))
       
       if (cell_contains(c, p)) then
          ! If this cell contains a universe or lattice, search for the particle
          ! in that universe/lattice
          if (c % type == CELL_NORMAL) then
             ! set current pointers
             found = .true.         
             cCell => c
             cMaterial => materials(cCell%material)
             cUniverse => univ
             
             ! set particle attributes
             p%cell = univ % cells(i)
             p%universe = dict_get_key(universe_dict, univ % uid)
             exit
          elseif (c % type == CELL_FILL) then
             lower_univ => universes(c % fill)
             call find_cell(lower_univ, p, found)
             if (found) then
                exit
             else
                msg = "Could not locate particle in universe: "
                call fatal_error(msg)
             end if
          elseif (c % type == CELL_LATTICE) then
             ! Set current lattice
             lat => lattices(c % fill)
             cLattice => lat
             p % lattice = c % fill

             ! determine universe based on lattice position
             x = ceiling((p%xyz(1) - lat%x0)/lat%width_x)
             y = ceiling((p%xyz(2) - lat%y0)/lat%width_y)
             lower_univ => universes(lat % element(x,y))

             ! adjust local position of particle
             p%xyz_local(1) = p%xyz(1) - (lat%x0 + (x-0.5)*lat%width_x)
             p%xyz_local(2) = p%xyz(2) - (lat%y0 + (y-0.5)*lat%width_y)
             p%xyz_local(3) = p%xyz(3)
             
             ! set particle lattice indices
             p % index_x = x
             p % index_y = y

             call find_cell(lower_univ, p, found)
             if (found) then
                exit
             else
                msg = "Could not locate particle in lattice: " & 
                     & // int_to_str(lat % uid)
                call fatal_error(msg)
             end if
          end if
       end if

    end do

  end subroutine find_cell

!===============================================================================
! CROSS_SURFACE moves a particle into a new cell
!===============================================================================

  subroutine cross_surface(p, last_cell)

    type(Particle), pointer :: p
    integer, intent(in)     :: last_cell  ! last cell particle was in

    integer                 :: i          ! index of neighbors
    integer                 :: index_cell ! index in cells array
    real(8)                 :: u          ! x-component of direction
    real(8)                 :: v          ! y-component of direction
    real(8)                 :: w          ! z-component of direction
    real(8)                 :: n1         ! x-component of surface normal
    real(8)                 :: n2         ! y-component of surface normal
    real(8)                 :: n3         ! z-component of surface normal
    real(8)                 :: norm       ! "norm" of surface normal
    logical                 :: found      ! particle found in universe?
    character(MAX_LINE_LEN) :: msg        ! output/error message?
    type(Surface),  pointer :: surf
    type(Cell),     pointer :: c
    type(Universe), pointer :: lower_univ => null()

    surf => surfaces(abs(p%surface))
    if (verbosity >= 10) then
       msg = "    Crossing surface " // trim(int_to_str(surf%uid))
       call message(msg)
    end if

    if (surf%bc == BC_VACUUM) then
       ! =======================================================================
       ! PARTICLE LEAKS OUT OF PROBLEM

       p%alive = .false.
       if (verbosity >= 10) then
          msg = "    Leaked out of surface " // trim(int_to_str(surf%uid))
          call message(msg)
       end if
       return

    elseif (surf%bc == BC_REFLECT) then
       ! =======================================================================
       ! PARTICLE REFLECTS FROM SURFACE

       ! Copy particle's direction cosines
       u = p % uvw(1)
       v = p % uvw(2)
       w = p % uvw(3)

       select case (surf%type)
       case (SURF_PX)
          p % uvw = (/ -u, v, w /)
       case (SURF_PY)
          p % uvw = (/ u, -v, w /)
       case (SURF_PZ)
          p % uvw = (/ u, v, -w /)
       case (SURF_PLANE)
          ! Find surface coefficients and norm of vector normal to surface
          n1 = surf % coeffs(1)
          n2 = surf % coeffs(2)
          n3 = surf % coeffs(3)
          norm = n1*n1 + n2*n2 + n3*n3

          ! Reflect direction according to normal
          u = u - 2*u*n1*n1/norm
          v = v - 2*v*n2*n2/norm
          w = w - 2*w*n3*n3/norm

          ! Set vector
          p % uvw = (/ u, v, w /)
       case default
          msg = "Reflection not supported for surface " // &
               trim(int_to_str(surf % uid))
          call fatal_error(msg)
       end select

       ! Reassign particle's cell and surface
       p % cell = last_cell
       p % surface = -p % surface

       ! Diagnostic message
       if (verbosity >= 10) then
          msg = "    Reflected from surface " // trim(int_to_str(surf%uid))
          call message(msg)
       end if
       return
    end if

    if (p%surface > 0 .and. allocated(surf%neighbor_pos)) then
       ! If coming from negative side of surface, search all the neighboring
       ! cells on the positive side
       do i = 1, size(surf%neighbor_pos)
          index_cell = surf%neighbor_pos(i)
          c => cells(index_cell)
          if (cell_contains(c, p)) then
             if (c % type == CELL_FILL) then
                lower_univ => universes(c % fill)
                call find_cell(lower_univ, p, found)
                if (.not. found) then
                   msg = "Could not locate particle in universe: "
                   call fatal_error(msg)
                end if
             else
                ! set current pointers
                p%cell = index_cell
                cCell => c
                cMaterial => materials(cCell%material)
             end if
             return
          end if
       end do
    elseif (p%surface < 0  .and. allocated(surf%neighbor_neg)) then
       ! If coming from positive side of surface, search all the neighboring
       ! cells on the negative side
       do i = 1, size(surf%neighbor_neg)
          index_cell = surf%neighbor_neg(i)
          c => cells(index_cell)
          if (cell_contains(c, p)) then
             if (c % type == CELL_FILL) then
                lower_univ => universes(c % fill)
                call find_cell(lower_univ, p, found)
                if (.not. found) then
                   msg = "Could not locate particle in universe: "
                   call fatal_error(msg)
                end if
             else
                ! set current pointers
                p%cell = index_cell
                cCell => c
                cMaterial => materials(cCell%material)
             end if
             return
          end if
       end do
    end if

    ! Couldn't find particle in neighboring cells, search through all cells
    do i = 1, size(cells)
       c => cells(i)
       if (cell_contains(c, p)) then
          p%cell = i
          cCell => c
          cMaterial => materials(cCell%material)
          return
       end if
    end do

    ! Couldn't find next cell anywhere!
    msg = "After particle crossed surface " // trim(int_to_str(p%surface)) // &
         & ", it could not be located in any cell and it did not leak."
    call fatal_error(msg)
       
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
    character(MAX_LINE_LEN) :: msg  ! output/error message
    type(Lattice),  pointer :: lat
    type(Universe), pointer :: univ

    if (verbosity >= 10) then
       msg = "    Crossing lattice"
       call message(msg)
    end if

    lat => lattices(p % lattice)

    u = p % uvw(1)
    v = p % uvw(2)

    if (lat % type == LATTICE_RECT) then
       x = p % xyz_local(1)
       y = p % xyz_local(2)
       z = p % xyz_local(3)
       x0 = lat % width_x * 0.5
       y0 = lat % width_y * 0.5
       
       dist = INFINITY

       ! left and right sides
       if (u > 0) then
          d_left = INFINITY
          d_right = (x0 - x)/u
       else
          d_left = -(x0 + x)/u
          d_right = INFINITY
       end if

       ! top and bottom sides
       if (v > 0) then
          d_bottom = INFINITY
          d_top = (y0 - y)/v
       else
          d_bottom = -(y0 + y)/v
          d_top = INFINITY
       end if

       dist = min(d_left, d_right, d_top, d_bottom)
       if (dist == d_left) then
          ! Move particle to left element
          p % index_x = p % index_x - 1
          p % xyz_local(1) = x0
          
       elseif (dist == d_right) then
          ! Move particle to right element
          p % index_x = p % index_x + 1
          p % xyz_local(1) = -x0

       elseif (dist == d_bottom) then
          ! Move particle to bottom element
          p % index_y = p % index_y - 1
          p % xyz_local(2) = y0

       elseif (dist == d_top) then
          ! Move particle to top element
          p % index_y = p % index_y + 1
          p % xyz_local(2) = -y0

       end if
    elseif (lat % type == LATTICE_HEX) then
       ! TODO: Add hex lattice support
    end if
    
    ! Check to make sure still in lattice
    i_x = p % index_x
    i_y = p % index_y
    if (i_x < 1 .or. i_x > lat % n_x) then
       msg = "Reached edge of lattice."
       call fatal_error(msg)
    elseif (i_y < 1 .or. i_y > lat % n_y) then
       msg = "Reached edge of lattice."
       call fatal_error(msg)
    end if

    ! Find universe for next lattice element
    univ => universes(lat % element(i_x,i_y))

    ! Find cell in next lattice element
    call find_cell(univ, p, found)
    if (.not. found) then
       msg = "Could not locate particle in universe: "
       call fatal_error(msg)
    end if

  end subroutine cross_lattice

!===============================================================================
! DIST_TO_BOUNDARY calculates the distance to the nearest boundary for a
! particle 'p' traveling in a certain direction. For a cell in a subuniverse
! that has a parent cell, also include the surfaces of the edge of the universe.
!===============================================================================

  subroutine dist_to_boundary(p, dist, surf, in_lattice)

    type(Particle), pointer     :: p
    real(8),        intent(out) :: dist
    integer,        intent(out) :: surf
    logical,        intent(out) :: in_lattice

    integer, allocatable :: expression(:) ! copy of surface list
    integer :: i            ! index for surface in cell
    integer :: n_surfaces   ! total number of surfaces to check
    integer :: n1           ! number of surfaces in current cell
    integer :: n2           ! number of surfaces in parent cell
    integer :: index_surf   ! index in surfaces array (with sign)
    integer :: current_surf ! current surface 
    real(8) :: x,y,z        ! particle coordinates
    real(8) :: u,v,w        ! particle directions
    real(8) :: d            ! evaluated distance
    real(8) :: x0,y0,z0     ! coefficients for surface
    real(8) :: r            ! radius for quadratic surfaces
    real(8) :: tmp          ! dot product of surface normal with direction
    real(8) :: a,b,c,k      ! quadratic equation coefficients
    real(8) :: quad         ! discriminant of quadratic equation
    logical :: on_surface   ! is particle on surface?
    character(MAX_LINE_LEN) :: msg   ! output/error message
    type(Cell),    pointer  :: cell_p => null()
    type(Cell),    pointer  :: parent_p => null()
    type(Surface), pointer  :: surf_p => null()
    type(Lattice), pointer  :: lat => null()

    cell_p => cells(p%cell)

    current_surf = p%surface

    ! determine number of surfaces to check
    n1 = cell_p % n_surfaces
    n2 = 0
    if (cell_p % parent > 0) then
       parent_p => cells(cell_p % parent)
       n2 = parent_p % n_surfaces
    end if
    n_surfaces = n1 + n2

    ! allocate space and assign expression
    allocate(expression(n_surfaces))
    expression(1:n1) = cell_p % surfaces
    if (cell_p % parent > 0) then
       expression(n1+1:n1+n2) = parent_p % surfaces
    end if

    u = p % uvw(1)
    v = p % uvw(2)
    w = p % uvw(3)

    ! loop over all surfaces
    dist = INFINITY
    do i = 1, n_surfaces
       if (i <= n1) then
          ! in local cell, so use xyz_local
          x = p % xyz_local(1)
          y = p % xyz_local(2)
          z = p % xyz_local(3)
       else
          ! in parent cell, so use xyz
          x = p % xyz(1)
          y = p % xyz(2)
          z = p % xyz(3)
       end if

       ! check for coincident surface -- note that we can't skip the calculation
       ! in general because a particle could be on one side of a cylinder and
       ! still have a positive distance to the other

       index_surf = expression(i)
       if (index_surf == current_surf) then
          on_surface = .true.
       else
          on_surface = .false.
       end if

       ! check for operators
       index_surf = abs(index_surf)
       if (index_surf >= OP_DIFFERENCE) cycle

       surf_p => surfaces(index_surf)
       select case (surf_p%type)
       case (SURF_PX)
          if (on_surface .or. u == ZERO) then
             d = INFINITY
          else
             x0 = surf_p % coeffs(1)
             d = (x0 - x)/u
             if (d < ZERO) d = INFINITY
          end if
          
       case (SURF_PY)
          if (on_surface .or. v == ZERO) then
             d = INFINITY
          else
             y0 = surf_p % coeffs(1)
             d = (y0 - y)/v
             if (d < ZERO) d = INFINITY
          end if
          
       case (SURF_PZ)
          if (on_surface .or. w == ZERO) then
             d = INFINITY
          else
             z0 = surf_p % coeffs(1)
             d = (z0 - z)/w
             if (d < ZERO) d = INFINITY
          end if
          
       case (SURF_PLANE)
          A = surf_p % coeffs(1)
          B = surf_p % coeffs(2)
          C = surf_p % coeffs(3)
          D = surf_p % coeffs(4)
          
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
             y0 = surf_p % coeffs(1)
             z0 = surf_p % coeffs(2)
             r = surf_p % coeffs(3)

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
                ! negative and one must be positive. The positive distance will
                ! be the one with negative sign on sqrt(quad)

                d = (-k + sqrt(quad))/a

             else
                ! particle is outside the cylinder, thus both distances are
                ! either positive or negative. If positive, the smaller distance
                ! is the one with positive sign on sqrt(quad)

                d = (-k - sqrt(quad))/a
                if (d < ZERO) d = INFINITY

             end if
          end if

       case (SURF_CYL_Y)
          a = ONE - v*v  ! u^2 + w^2
          if (a == ZERO) then
             d = INFINITY
          else
             x0 = surf_p % coeffs(1)
             z0 = surf_p % coeffs(2)
             r = surf_p % coeffs(3)

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
                ! negative and one must be positive. The positive distance will
                ! be the one with negative sign on sqrt(quad)

                d = (-k + sqrt(quad))/a

             else
                ! particle is outside the cylinder, thus both distances are
                ! either positive or negative. If positive, the smaller distance
                ! is the one with positive sign on sqrt(quad)

                d = (-k - sqrt(quad))/a
                if (d < ZERO) d = INFINITY

             end if
          end if

       case (SURF_CYL_Z)
          a = ONE - w*w  ! u^2 + v^2
          if (a == ZERO) then
             d = INFINITY
          else
             x0 = surf_p % coeffs(1)
             y0 = surf_p % coeffs(2)
             r = surf_p % coeffs(3)

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
                ! negative and one must be positive. The positive distance will
                ! be the one with negative sign on sqrt(quad)

                d = (-k + sqrt(quad))/a

             else
                ! particle is outside the cylinder, thus both distances are
                ! either positive or negative. If positive, the smaller distance
                ! is the one with positive sign on sqrt(quad)

                d = (-k - sqrt(quad))/a
                if (d <= ZERO) d = INFINITY

             end if
          end if

       case (SURF_SPHERE)
          x0 = surf_p % coeffs(1)
          y0 = surf_p % coeffs(2)
          z0 = surf_p % coeffs(3)
          r = surf_p % coeffs(4)

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
             ! particle is on the sphere, thus one distance is positive/negative
             ! and the other is zero. The sign of k determines if we are facing
             ! in or out
             
             if (k >= ZERO) then
                d = INFINITY
             else
                d = -k + sqrt(quad)
             end if

          elseif (c < ZERO) then
             ! particle is inside the sphere, thus one distance must be negative
             ! and one must be positive. The positive distance will be the one
             ! with negative sign on sqrt(quad)

             d = -k + sqrt(quad)

          else
             ! particle is outside the sphere, thus both distances are either
             ! positive or negative. If positive, the smaller distance is the
             ! one with positive sign on sqrt(quad)

             d = -k - sqrt(quad)
             if (d < ZERO) d = INFINITY

          end if

       case (SURF_GQ)
          msg = "Surface distance not yet implement for general quadratic."
          call fatal_error(msg)

       end select

       ! Check is calculated distance is new minimum
       if (d < dist) then
          dist = d
          surf = -expression(i)
       end if

    end do

    ! Check lattice surfaces
    in_lattice = .false.
    if (p % lattice > 0) then
       lat => lattices(p % lattice)
       if (lat % type == LATTICE_RECT) then
          x = p % xyz_local(1)
          y = p % xyz_local(2)
          z = p % xyz_local(3)
          x0 = lat % width_x * 0.5
          y0 = lat % width_y * 0.5
          
          
          ! left and right sides
          if (u > 0) then
             d = (x0 - x)/u
          else
             d = -(x + x0)/u
          end if
          if (d < dist) then
             dist = d
             in_lattice = .true.
          end if

          ! top and bottom sides
          if (v > 0) then
             d = (y0 - y)/v
          else
             d = -(y + y0)/v
          end if
          if (d < dist) then
             dist = d
             in_lattice = .true.
          end if

       elseif (lat % type == LATTICE_HEX) then
          ! TODO: Add hex lattice support
       end if
    end if

    ! deallocate expression
    deallocate(expression)

  end subroutine dist_to_boundary

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
    integer :: index      ! index in count arrays
    integer, allocatable :: count_positive(:) ! # of cells on positive side
    integer, allocatable :: count_negative(:) ! # of cells on negative side
    logical :: positive   ! positive side specified in surface list
    character(MAX_LINE_LEN) :: msg ! output/error message
    type(Cell),    pointer  :: c
    type(Surface), pointer  :: surf

    msg = "Building neighboring cells lists for each surface..."
    call message(msg, 4)

    allocate(count_positive(n_surfaces))
    allocate(count_negative(n_surfaces))
    count_positive = 0
    count_negative = 0

    do i = 1, n_cells
       c => cells(i)

       ! loop over each surface specification
       do j = 1, c % n_surfaces
          index = c % surfaces(j)
          positive = (index > 0)
          index = abs(index)
          if (positive) then
             count_positive(index) = count_positive(index) + 1
          else
             count_negative(index) = count_negative(index) + 1
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
          index = c % surfaces(j)
          positive = (index > 0)
          index = abs(index)

          surf => surfaces(index)
          if (positive) then
             count_positive(index) = count_positive(index) + 1
             surf%neighbor_pos(count_positive(index)) = i
          else
             count_negative(index) = count_negative(index) + 1
             surf%neighbor_neg(count_negative(index)) = i
          end if
       end do
    end do

    deallocate(count_positive)
    deallocate(count_negative)

  end subroutine neighbor_lists

end module geometry
