module geometry

  use global
  use types,  only: Cell, Surface
  use output, only: error, message
  use data_structures, only: dict_get_key

  implicit none
     
contains

!=====================================================================
! CELL_CONTAINS determines whether a given point is inside a cell
!=====================================================================

  function cell_contains(c, p) result(in_cell)

    type(Cell),     pointer :: c
    type(Particle), pointer :: p
    logical                 :: in_cell

    integer, allocatable :: expression(:)
    integer :: specified_sense
    integer :: actual_sense
    integer :: n_surfaces
    integer :: i, j
    integer :: surf_num
    integer :: current_surface
    type(Surface), pointer :: surf => null()
    character(250) :: msg

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
          ! particle is on specified surface, but heading other
          ! direction
          expression(i) = 0
          cycle
       end if

       ! Compare sense of point to specified sense
       specified_sense = sign(1,expression(i))
       call sense(surf, p%xyz, actual_sense)
       if (actual_sense == specified_sense) then
          expression(i) = 1
       else
          expression(i) = 0
       end if

    end do

    ! TODO: Need to replace this with a 'lgeval' like subroutine that
    ! can actually test expressions with unions and parentheses
    if (all(expression == 1)) then
       in_cell = .true.
    else
       in_cell = .false.
    end if

    ! Free up memory from expression
    deallocate(expression)

  end function cell_contains

!=====================================================================
! FIND_CELL determines what cell a source particle is in within a
! particular universe. If the base universe is passed, the particle
! should be found as long as it's within the geometry
!=====================================================================

  recursive subroutine find_cell(univ, p, found)

    type(Universe), pointer :: univ       ! universe to search in
    type(Particle), pointer :: p          ! pointer to particle
    logical,  intent(inout) :: found      ! particle found?

    character(250) :: msg                 ! error message
    integer        :: i                   ! index over cells
    type(Cell),     pointer :: c          ! pointer to cell
    type(Universe), pointer :: lower_univ ! if particle is in lower
                                          ! universe, use this pointer
                                          ! to call recursively

    found = .false.

    ! determine what region in
    do i = 1, univ % n_cells
       c => cells(univ % cells(i))
       
       if (cell_contains(c, p)) then
          p%cell = dict_get_key(cell_dict, c % uid)

          ! If this cell contains a universe of lattice, search for
          ! the particle in that universe/lattice
          if (c % fill > 0) then
             lower_univ => universes(c % fill)
             call find_cell(lower_univ, p, found)
             if (found) then
                exit
             else
                msg = "Could not locate particle in universe: "
                call error(msg)
             end if
          else
             ! set current pointers
             found = .true.         
             cCell => c
             cMaterial => materials(cCell%material)
             exit
          end if
       end if

    end do

  end subroutine find_cell

!=====================================================================
! CROSS_BOUNDARY moves a particle into a new cell
!=====================================================================

  subroutine cross_boundary(p)

    type(Particle), pointer :: p

    type(Surface), pointer :: surf
    type(Cell),    pointer :: c
    integer :: i
    integer :: index_cell
    character(250) :: msg
    logical :: found
    type(Universe), pointer :: lower_univ => null()

    surf => surfaces(abs(p%surface))
    if (verbosity >= 10) then
       msg = "    Crossing surface " // trim(int_to_str(surf%uid))
       call message(msg, 10)
    end if

    ! check for leakage
    if (surf%bc == BC_VACUUM) then
       p%alive = .false.
       if (verbosity >= 10) then
          msg = "    Leaked out of surface " // trim(int_to_str(surf%uid))
          call message(msg, 10)
       end if
       return
    end if

    if (p%surface > 0 .and. allocated(surf%neighbor_pos)) then
       ! If coming from negative side of surface, search all the
       ! neighboring cells on the positive side
       do i = 1, size(surf%neighbor_pos)
          index_cell = surf%neighbor_pos(i)
          c => cells(index_cell)
          if (cell_contains(c, p)) then
             if (c % fill > 0) then
                lower_univ => universes(c % fill)
                call find_cell(lower_univ, p, found)
                if (.not. found) then
                   msg = "Could not locate particle in universe: "
                   call error(msg)
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
       ! If coming from positive side of surface, search all the
       ! neighboring cells on the negative side
       do i = 1, size(surf%neighbor_neg)
          index_cell = surf%neighbor_neg(i)
          c => cells(index_cell)
          if (cell_contains(c, p)) then
             if (c % fill > 0) then
                lower_univ => universes(c % fill)
                call find_cell(lower_univ, p, found)
                if (.not. found) then
                   msg = "Could not locate particle in universe: "
                   call error(msg)
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

    ! Couldn't find particle in neighboring cells, search through all
    ! cells
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
    call error(msg)
       
  end subroutine cross_boundary

!=====================================================================
! DIST_TO_BOUNDARY calculates the distance to the nearest boundary of
! the cell 'cl' for a particle 'p' traveling in a certain
! direction. For a cell in a subuniverse that has a parent cell, also
! include the surfaces of the edge of the universe.
!=====================================================================

  subroutine dist_to_boundary(p, dist, surf, other_cell)

    type(Particle),    pointer     :: p
    real(8),           intent(out) :: dist
    integer,           intent(out) :: surf
    integer, optional, intent(in)  :: other_cell

    type(Cell),    pointer :: cell_p => null()
    type(Cell),    pointer :: parent_p => null()
    type(Surface), pointer :: surf_p => null()
    integer :: i
    integer :: n_surfaces, n1, n2
    integer, allocatable :: expression(:)
    integer :: index_surf
    integer :: current_surf
    real(8) :: x,y,z,u,v,w
    real(8) :: d
    real(8) :: x0,y0,z0,r
    real(8) :: tmp
    real(8) :: a,b,c,k
    real(8) :: quad
    character(250) :: msg

    if (present(other_cell)) then
       cell_p => cells(other_cell)
    else
       cell_p => cells(p%cell)
    end if

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

    ! loop over all surfaces
    dist = INFINITY
    do i = 1, n_surfaces
       x = p % xyz(1)
       y = p % xyz(2)
       z = p % xyz(3)
       u = p % uvw(1)
       v = p % uvw(2)
       w = p % uvw(3)

       ! check for coincident surface
       index_surf = expression(i)
       if (index_surf == current_surf) cycle

       ! check for operators
       index_surf = abs(index_surf)
       if (index_surf >= OP_DIFFERENCE) cycle

       surf_p => surfaces(index_surf)
       select case (surf_p%type)
       case (SURF_PX)
          if (u == ZERO) then
             d = INFINITY
          else
             x0 = surf_p % coeffs(1)
             d = (x0 - x)/u
             if (d < ZERO) d = INFINITY
          end if
          
       case (SURF_PY)
          if (v == ZERO) then
             d = INFINITY
          else
             y0 = surf_p % coeffs(1)
             d = (y0 - y)/v
             if (d < ZERO) d = INFINITY
          end if
          
       case (SURF_PZ)
          if (w == ZERO) then
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
          if (tmp == ZERO) then
             d = INFINITY
          else
             d = -(A*x + B*y + C*w - D)/tmp
             if (d < ZERO) d = INFINITY
          end if

       case (SURF_CYL_X)
          a = ONE - u**2  ! v^2 + w^2
          if (a == ZERO) then
             d = INFINITY
          else
             y0 = surf_p % coeffs(1)
             z0 = surf_p % coeffs(2)
             r = surf_p % coeffs(3)

             y = y - y0
             z = z - z0
             k = y*v + z*w
             c = y**2 + z**2 - r**2
             quad = k**2 - a*c

             if (c < ZERO) then
                ! particle is inside the cylinder, thus one distance
                ! must be negative and one must be positive. The
                ! positive distance will be the one with negative sign
                ! on sqrt(quad)

                d = -(k - sqrt(quad))/a

             else
                ! particle is outside the cylinder, thus both
                ! distances are either positive or negative. If
                ! positive, the smaller distance is the one with
                ! positive sign on sqrt(quad)

                d = -(k + sqrt(quad))/a
                if (d < ZERO) d = INFINITY

             end if
          end if

       case (SURF_CYL_Y)
          a = ONE - v**2  ! u^2 + w^2
          if (a == ZERO) then
             d = INFINITY
          else
             x0 = surf_p % coeffs(1)
             z0 = surf_p % coeffs(2)
             r = surf_p % coeffs(3)

             x = x - x0
             z = z - z0
             k = x*u + z*w
             c = x**2 + z**2 - r**2
             quad = k**2 - a*c

             if (c < ZERO) then
                ! particle is inside the cylinder, thus one distance
                ! must be negative and one must be positive. The
                ! positive distance will be the one with negative sign
                ! on sqrt(quad)

                d = -(k - sqrt(quad))/a

             else
                ! particle is outside the cylinder, thus both
                ! distances are either positive or negative. If
                ! positive, the smaller distance is the one with
                ! positive sign on sqrt(quad)

                d = -(k + sqrt(quad))/a
                if (d < ZERO) d = INFINITY

             end if
          end if

       case (SURF_CYL_Z)
          a = ONE - w**2  ! u^2 + v^2
          if (a == ZERO) then
             d = INFINITY
          else
             x0 = surf_p % coeffs(1)
             y0 = surf_p % coeffs(2)
             r = surf_p % coeffs(3)

             x = x - x0
             y = y - y0
             k = x*u + y*v
             c = x**2 + y**2 - r**2
             quad = k**2 - a*c
             
             if (quad < ZERO) then
                ! no intersection with cylinder

                d = INFINITY 

             elseif (c < ZERO) then
                ! particle is inside the cylinder, thus one distance
                ! must be negative and one must be positive. The
                ! positive distance will be the one with negative sign
                ! on sqrt(quad)

                d = -(k - sqrt(quad))/a

             else
                ! particle is outside the cylinder, thus both
                ! distances are either positive or negative. If
                ! positive, the smaller distance is the one with
                ! positive sign on sqrt(quad)

                d = -(k + sqrt(quad))/a
                if (d < ZERO) d = INFINITY

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
          c = x**2 + y**2 + z**2 - r**2
          quad = k**2 - c

          if (quad < ZERO) then
             ! no intersection with sphere

             d = INFINITY 

          elseif (c < ZERO) then
             ! particle is inside the sphere, thus one distance
             ! must be negative and one must be positive. The
             ! positive distance will be the one with negative sign
             ! on sqrt(quad)

             d = -(k - sqrt(quad))

          else
             ! particle is outside the sphere, thus both
             ! distances are either positive or negative. If
             ! positive, the smaller distance is the one with
             ! positive sign on sqrt(quad)

             d = -(k + sqrt(quad))
             if (d < ZERO) d = INFINITY

          end if

       case (SURF_GQ)
          msg = "Surface distance not yet implement for general quadratic."
          call error(msg)

       end select

       ! Check is calculated distance is new minimum
       if (d < dist) then
          dist = d
          surf = -expression(i)
       end if

    end do

    ! deallocate expression
    deallocate(expression)

  end subroutine dist_to_boundary

!=====================================================================
! SENSE determines whether a point is on the 'positive' or 'negative'
! side of a surface. This routine is crucial for determining what cell
! a particular point is in.
!=====================================================================

  subroutine sense(surf, xyz, s)

    type(Surface), intent(in) :: surf
    real(8), intent(in) :: xyz(3)
    integer, intent(out) :: s

    real(8) :: x,y,z
    real(8) :: func
    real(8) :: A,B,C,D,E,F,G,H,I,J
    real(8) :: x0, y0, z0, r
    real(8) :: x1, y1, z1

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
       func = (y-y0)**2 + (z-z0)**2 - r**2

    case (SURF_CYL_Y)
       x0 = surf % coeffs(1)
       z0 = surf % coeffs(2)
       r = surf % coeffs(3)
       func = (x-x0)**2 + (z-z0)**2 - r**2

    case (SURF_CYL_Z)
       x0 = surf % coeffs(1)
       y0 = surf % coeffs(2)
       r = surf % coeffs(3)
       func = (x-x0)**2 + (y-y0)**2 - r**2

    case (SURF_SPHERE)
       x0 = surf % coeffs(1)
       y0 = surf % coeffs(2)
       z0 = surf % coeffs(3)
       r = surf % coeffs(4)
       func = (x-x0)**2 + (y-y0)**2 + (z-z0)**2 - r**2

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
       func = A*x**2 + B*y**2 + C*z**2 + D*x*y + E*y*z + F*x*z + G*x &
            & + H*y + I*z + J

    end select

    ! Check which side of surface the point is on
    if (func > 0) then
       s = SENSE_POSITIVE
    else
       s = SENSE_NEGATIVE
    end if

  end subroutine sense

!=====================================================================
! NEIGHBOR_LISTS builds a list of neighboring cells to each surface to
! speed up searches when a cell boundary is crossed.
!=====================================================================

  subroutine neighbor_lists()

    type(Cell),    pointer :: c
    type(Surface), pointer :: surf
    integer :: i, j
    integer :: index
    integer :: surf_num
    logical :: positive
    character(250) :: msg

    integer, allocatable :: count_positive(:)
    integer, allocatable :: count_negative(:)

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

  end subroutine neighbor_lists

end module geometry
