module geometry

  use global
  use types,  only: Cell, Surface
  use output, only: error

  implicit none
     
contains

!=====================================================================
! CELL_CONTAINS determines whether a given point is inside a cell
!=====================================================================

  function cell_contains( c, xyz )

    type(Cell), intent(in)  :: c
    real(8),    intent(in)  :: xyz(3)
    logical :: cell_contains

    integer, allocatable :: expression(:)
    integer :: specified_sense
    integer :: actual_sense
    integer :: n_boundaries
    integer :: i, j
    integer :: surf_num
    type(Surface), pointer :: surf => null()
    logical :: surface_found
    character(250) :: msg

    n_boundaries = size(c%boundary_list)
    allocate( expression(n_boundaries) )
    expression = c%boundary_list
    do i = 1, n_boundaries
       surface_found = .false.

       ! Don't change logical operator
       if ( expression(i) >= OP_DIFFERENCE ) then
          cycle
       end if

       ! Lookup surface
       surf_num = abs(expression(i))
       ! TODO: replace this loop with a hash since this lookup is O(N)
       do j = 1, n_surfaces
          surf => surfaces(j)
          if ( surf%uid == surf_num ) then
             surface_found = .true.
             exit
          end if
       end do

       ! Report error if can't find specified surface
       if ( .not. surface_found ) then
          deallocate( expression )
          msg = "Count not find surface " // trim(int_to_str(surf_num)) // & 
               & " specified on cell " // int_to_str(c%uid)
          call error( msg )
       end if

       ! Compare sense of point to specified sense
       specified_sense = sign(1,expression(i))
       call sense( surf, xyz, actual_sense )
       if ( actual_sense == specified_sense ) then
          expression(i) = 1
       else
          expression(i) = 0
       end if

    end do

    ! TODO: Need to replace this with a 'lgeval' like subroutine that
    ! can actually test expressions with unions and parentheses
    if ( all( expression == 1 ) ) then
       cell_contains = .true.
    else
       cell_contains = .false.
    end if

    ! Free up memory from expression
    deallocate( expression )

  end function cell_contains

!=====================================================================
! FIND_CELL determines what cell a source neutron is in
!=====================================================================

  subroutine find_cell( neut )

    type(Neutron), pointer, intent(inout) :: neut

    logical :: found_cell
    character(250) :: msg
    integer :: i

    ! determine what region in
    do i = 1, n_cells
       if ( cell_contains(cells(i), neut%xyz) ) then
          neut%cell = i
          found_cell = .true.
          exit
       end if
    end do

    ! if neutron couldn't be located, print error
    if ( .not. found_cell ) then
       write(msg, 100) "Could not locate cell for neutron at: ", neut%xyz
100    format (A,3ES11.3)
       call error( msg )
    end if

  end subroutine find_cell

!=====================================================================
! CROSS_BOUNDARY moves a neutron into a new cell
!=====================================================================

  subroutine cross_boundary( neut )

    type(Neutron), pointer, intent(in) :: neut

    type(Surface), pointer :: surf
    type(Cell), pointer :: c
    integer :: i
    integer :: index_cell
    character(100) :: msg

    surf => surfaces(abs(neut%surface))

    ! check for leakage
    if ( surf%bc == BC_VACUUM ) then
       neut%alive = .false.
       return
    end if

    if ( neut%surface < 0 ) then
       ! If coming from negative side of surface, search all the
       ! neighboring cells on the positive side
       do i = 1, size(surf%neighbor_pos(:))
          index_cell = surf%neighbor_pos(i)
          c = cells(index_cell)
          if ( cell_contains(c, neut%xyz) ) then
             neut%cell = index_cell
             return
          end if
       end do
    elseif ( neut%surface > 0 ) then
       ! If coming from positive side of surface, search all the
       ! neighboring cells on the negative side
       do i = 1, size(surf%neighbor_neg(:))
          index_cell = surf%neighbor_neg(i)
          c = cells(index_cell)
          if ( cell_contains(c, neut%xyz) ) then
             neut%cell = index_cell
             return
          end if
       end do
    end if

    ! Couldn't find particle in neighboring cells, search through all
    ! cells
    do i = 1, size(cells)
       c = cells(i)
       if ( cell_contains(c, neut%xyz) ) then
          neut%cell = i
          return
       end if
    end do

    ! Couldn't find next cell anywhere!
    msg = "After neutron crossed surface " // int_to_str(neut%surface) // &
         & ", it could not be located in any cell and it did not leak."
    call error( msg )
       
  end subroutine cross_boundary

!=====================================================================
! DIST_TO_BOUNDARY calculates the distance to the nearest boundary of
! the cell 'cl' for a particle 'neut' traveling in a certain
! direction.
!=====================================================================

  subroutine dist_to_boundary( neut, dist, surf, other_cell )

    type(Neutron),     intent(in)  :: neut
    real(8),           intent(out) :: dist
    integer,           intent(out) :: surf
    integer, optional, intent(in)  :: other_cell

    type(Cell),    pointer :: cell_p => null()
    type(Surface), pointer :: surf_p => null()
    integer :: i
    integer :: n_boundaries
    integer, allocatable :: expression(:)
    integer :: index_surf
    real(8) :: x,y,z,u,v,w
    real(8) :: d
    real(8) :: x0,y0,z0,r
    real(8) :: tmp
    real(8) :: a,b,c,k
    real(8) :: quad
    character(250) :: msg

    if ( present(other_cell) ) then
       cell_p => cells(other_cell)
    else
       cell_p => cells(neut%cell)
    end if

    x = neut%xyz(1)
    y = neut%xyz(2)
    z = neut%xyz(3)
    u = neut%uvw(1)
    v = neut%uvw(2)
    z = neut%uvw(3)

    dist = INFINITY
    n_boundaries = size(cell_p%boundary_list)
    allocate( expression(n_boundaries) )
    expression = cell_p%boundary_list
    do i = 1, n_boundaries
       index_surf = abs(expression(i))
       if ( index_surf >= OP_DIFFERENCE ) cycle

       surf_p => surfaces(index_surf)
       select case ( surf_p%type )
       case ( SURF_PX )
          if ( u == 0.0 ) then
             d = INFINITY
          else
             x0 = surf_p%coeffs(1)
             d = (x0 - x)/u
             if ( d < 0 ) d = INFINITY
          end if
          
       case ( SURF_PY )
          if ( v == 0.0 ) then
             d = INFINITY
          else
             y0 = surf_p%coeffs(1)
             d = (y0 - y)/v
             if ( d < 0 ) d = INFINITY
          end if
          
       case ( SURF_PZ )
          if ( w == 0.0 ) then
             d = INFINITY
          else
             z0 = surf_p%coeffs(1)
             d = (z0 - z)/w
             if ( d < 0.0 ) d = INFINITY
          end if
          
       case ( SURF_PLANE )
          A = surf_p%coeffs(1)
          B = surf_p%coeffs(2)
          C = surf_p%coeffs(3)
          D = surf_p%coeffs(4)
          
          tmp = A*u + B*v + C*w
          if ( tmp == 0.0 ) then
             d = INFINITY
          else
             d = -(A*x + B*y + C*w - D)/tmp
             if ( d < 0.0 ) d = INFINITY
          end if

       case ( SURF_CYL_X )
          a = 1.0 - u**2  ! v^2 + w^2
          if ( a == 0.0 ) then
             d = INFINITY
          else
             y0 = surf_p%coeffs(1)
             z0 = surf_p%coeffs(2)
             r = surf_p%coeffs(3)

             y = y - y0
             z = z - z0
             k = y*v + z*w
             c = y**2 + z**2 - r**2
             quad = k**2 - a*c

             if ( c < 0 ) then
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
                if ( d < 0 ) d = INFINITY

             end if
          end if

       case ( SURF_CYL_Y )
          a = 1.0 - v**2  ! u^2 + w^2
          if ( a == 0.0 ) then
             d = INFINITY
          else
             x0 = surf_p%coeffs(1)
             z0 = surf_p%coeffs(2)
             r = surf_p%coeffs(3)

             x = x - x0
             z = z - z0
             k = x*u + z*w
             c = x**2 + z**2 - r**2
             quad = k**2 - a*c

             if ( c < 0 ) then
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
                if ( d < 0 ) d = INFINITY

             end if
          end if

       case ( SURF_CYL_Z )
          a = 1.0 - w**2  ! u^2 + v^2
          if ( a == 0.0 ) then
             d = INFINITY
          else
             x0 = surf_p%coeffs(1)
             y0 = surf_p%coeffs(2)
             r = surf_p%coeffs(3)

             x = x - x0
             y = y - y0
             k = x*u + y*v
             c = x**2 + y**2 - r**2
             quad = k**2 - a*c
             
             if ( quad < 0 ) then
                ! no intersection with cylinder

                d = INFINITY 

             elseif ( c < 0 ) then
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
                if ( d < 0 ) d = INFINITY

             end if
          end if

       case ( SURF_SPHERE )
          x0 = surf_p%coeffs(1)
          y0 = surf_p%coeffs(2)
          z0 = surf_p%coeffs(3)
          r = surf_p%coeffs(4)

          x = x - x0
          y = y - y0
          z = z - z0
          k = x*u + y*v + z*w
          c = x**2 + y**2 + z**2 - r**2
          quad = k**2 - c

          if ( quad < 0 ) then
             ! no intersection with sphere

             d = INFINITY 

          elseif ( c < 0 ) then
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
             if ( d < 0 ) d = INFINITY

          end if

       case ( SURF_GQ )
          msg = "Surface distance not yet implement for general quadratic."
          call error( msg )

       end select

       ! Check is calculated distance is new minimum
       if ( d < dist ) then
          dist = d
          surf = expression(i)
       end if

    end do

    ! deallocate expression
    deallocate( expression )

  end subroutine dist_to_boundary

!=====================================================================
! SENSE determines whether a point is on the 'positive' or 'negative'
! side of a surface. This routine is crucial for determining what cell
! a particular point is in.
!=====================================================================

  subroutine sense( surf, xyz, s )

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

    select case ( surf%type )
    case ( SURF_PX )
       x0 = surf%coeffs(1)
       func = x - x0
       
    case ( SURF_PY )
       y0 = surf%coeffs(1)
       func = y - y0
       
    case ( SURF_PZ )
       z0 = surf%coeffs(1)
       func = z - z0
       
    case ( SURF_PLANE )
       A = surf%coeffs(1)
       B = surf%coeffs(2)
       C = surf%coeffs(3)
       D = surf%coeffs(4)
       func = A*x + B*y + C*z - D

    case ( SURF_CYL_X )
       y0 = surf%coeffs(1)
       z0 = surf%coeffs(2)
       r = surf%coeffs(3)
       func = (y-y0)**2 + (z-z0)**2 - r**2

    case ( SURF_CYL_Y )
       x0 = surf%coeffs(1)
       z0 = surf%coeffs(2)
       r = surf%coeffs(3)
       func = (x-x0)**2 + (z-z0)**2 - r**2

    case ( SURF_CYL_Z )
       x0 = surf%coeffs(1)
       y0 = surf%coeffs(2)
       r = surf%coeffs(3)
       func = (x-x0)**2 + (y-y0)**2 - r**2

    case ( SURF_SPHERE )
       x0 = surf%coeffs(1)
       y0 = surf%coeffs(2)
       z0 = surf%coeffs(3)
       r = surf%coeffs(4)
       func = (x-x0)**2 + (y-y0)**2 + (z-z0)**2 - r**2

    case ( SURF_BOX_X )
       y0 = surf%coeffs(1)
       z0 = surf%coeffs(2)
       y1 = surf%coeffs(3)
       z1 = surf%coeffs(4)
       if ( y >= y0 .and. y < y1 .and. z >= z0 .and. z < z1 ) then
          s = SENSE_NEGATIVE
       else
          s = SENSE_POSITIVE
       end if
       return

    case ( SURF_BOX_Y )
       x0 = surf%coeffs(1)
       z0 = surf%coeffs(2)
       x1 = surf%coeffs(3)
       z1 = surf%coeffs(4)
       if ( x >= x0 .and. x < x1 .and. z >= z0 .and. z < z1 ) then
          s = SENSE_NEGATIVE
       else
          s = SENSE_POSITIVE
       end if
       return

    case ( SURF_BOX_Z )
       x0 = surf%coeffs(1)
       y0 = surf%coeffs(2)
       x1 = surf%coeffs(3)
       y1 = surf%coeffs(4)
       if ( x >= x0 .and. x < x1 .and. y >= y0 .and. y < y1 ) then
          s = SENSE_NEGATIVE
       else
          s = SENSE_POSITIVE
       end if
       return

    case ( SURF_BOX )
       x0 = surf%coeffs(1)
       y0 = surf%coeffs(2)
       z0 = surf%coeffs(3)
       x1 = surf%coeffs(4)
       y1 = surf%coeffs(5)
       z1 = surf%coeffs(6)
       if ( x >= x0 .and. x < x1 .and. y >= y0 .and. y < y1 .and. & 
            & z >= z0 .and. z < z1 ) then
          s = SENSE_NEGATIVE
       else
          s = SENSE_POSITIVE
       end if
       return

    case ( SURF_GQ )
       func = A*x**2 + B*y**2 + C*z**2 + D*x*y + E*y*z + F*x*z + G*x &
            & + H*y + I*z + J

    end select

    ! Check which side of surface the point is on
    if ( func > 0 ) then
       s = SENSE_POSITIVE
    else
       s = SENSE_NEGATIVE
    end if

  end subroutine sense

end module geometry
