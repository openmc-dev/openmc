module m_contours
#ifndef DUMMYLIB
!  use m_common_format
! change to by GT to test in FoX4.0.3
   use fox_common
   use m_common_error

  implicit none
  private

  type richline
    real, pointer :: x(:) => null()
    real, pointer :: y(:) => null()
    integer, pointer :: contains_lt(:) => null() ! used for polygons
    integer, pointer :: contains_eq(:) => null() ! used for polygons
    integer, pointer :: contains_gt(:) => null() ! used for polygons
    integer, pointer :: borders_lt(:) => null() ! used for polygons
    integer, pointer :: borders_eq(:) => null() ! used for polygons
    integer, pointer :: borders_gt(:) => null() ! used for polygons
    integer :: mountain = 0 ! used for polygons - is this a valley (1) or a mountain (2)?
    logical :: boundary = .false. ! Does part of this line run along a boundary?
    integer :: se = 0           ! used temporarily for unterminated lines
    integer :: ee = 0           ! used temporarily for unterminated lines
    logical :: sedone = .false. ! used temporarily for unterminated lines
    logical :: eedone = .false. ! used temporarily for unterminated lines
    logical :: uphill           ! when this line terminates at an edge, is that edge going clockwise uphill?
    real :: area
  end type richline

  type line
    real, pointer :: x(:) => null()
    real, pointer :: y(:) => null()
  end type line

  type lineList
    type(line), pointer :: list(:) => null()
  end type lineList
  
  type richlineList
    type(richline), pointer :: list(:) => null()
  end type richlineList

  type contourLines
    type(richlineList), pointer :: polys(:) => null()
    type(richlineList), pointer :: lines(:) => null()
    type(richline) :: current
    integer :: current_se, current_ee
    real, pointer :: contours(:)
  end type contourLines

  type polygon
    type(line) :: outerBoundary
    type(lineList) :: innerBoundaries
  end type polygon

  type polygonList
    type(polygon), pointer :: list(:) => null()
  end type polygonList

  type contourObject
    type(lineList), pointer :: lines(:) => null() ! List of lines at each contour level, 0th is boundary box
    type(polygonList), pointer :: polys(:) => null() ! List of polygons for each between-contour-layer.
    real, pointer :: contours(:) => null()
    real, pointer :: zmax => null()
  end type contourObject

  public :: line
  public :: lineList
  public :: polygon
  public :: polygonList
  public :: contourObject

  public :: destroy_contourObject
  public :: make_contours_on_simplest_grid  
  public :: make_contours_on_a_simple_grid
  public :: make_contours_on_a_complex_grid

  logical :: debug = .false.

contains

! EXPLANATION
! This module provides a set of routines to calculate contour lines & contour polygons from a set of 2D data on
! a putatively rectangular grid (ie the grid can be scaled later, but all interpolation herein is strictly
! linear and univariate)

! There are only two points of entry:
! make_contours_on_a_simple_grid - which makes the contours and scales them onto a simple rectangular grid
! make_contours_on_a_complex_grid - which does bivariate linear interpolation onto a topologically rectangular grid.

! The product of either of these is a contourObject, as described above.
! This contains:
!   a real array of the contour values
!   a list of lists of lines, one list for each contour value, and a 0th list containing the boundary box.
!   a list of lists of polygons, one list for each between-contour-values, one above, and one below.
!   each polygon contains an outerBoundary, and a list of innerBoundaries, each of which are lines
!   each line is two real arrays, one of x values, one of y.

! The algorithm must give us:
!  A set of contour lines
!  A set of polygons (potentially within inner and outer boundaries) and knowledge
!   of whether the contained points are greater or smaller than the associated contour level.

! The algorithm is as follows:
!  1: use TOMS 531 to generate contour lines.
!      This will produce some whole polygons, and some lines (terminating at boundaries.
!  2: If there are no lines, then add the box perimeter as an additional polygon.
!  3: For all lines terminating at boundaries, walk round the box, constructing the smallest 
!      possible polygon consisting of line segments and boundary segments. Note that this polygon
!      references only contours at its own level, or the level above.
!  4: Search over all contour lines at a given level (both those which are polygons as found
!      in step 1, and joined up lines from step 2) and see if any other polygons - from the same
!      contour level, or one above or one below - are within them.
!     If they are, we can tell whether the area between them contains points higher or lower
!      than the given contour value, from simple logic.
!  5: Check over all polygons which contain boundary segments - which other boundary polygons
!      do they share an edge with (these must be either at the same contour level, or one above).
!  6: There are probably some polygons whose status we don't know about - but we can work this
!      out from the fact that they border polygons whose status we *do* know about.
!  7: Finally - it is possible that only one contour level is visible. If so, then working out
!      polygon status by logic above fails - we need to test a point inside the polygon. Do so,
!      and propagate this information.
!  8: Copy all the information from the temporary structure we've been using into
!      a contourObject (duplicating all innerBoundaries).
!  9: Scale all points as appropriate.



  subroutine init_line(thisLine)
    type(richline), intent(inout) :: thisLine

    allocate(thisLine%x(0))
    allocate(thisLine%y(0))
  end subroutine init_line

  subroutine destroy_contourLines(cp)
    type(contourLines), intent(inout) :: cp
    integer :: i, j

    do i = 0, size(cp%contours)
      do j = 1, size(cp%polys(i)%list)
        deallocate(cp%polys(i)%list(j)%x)
        deallocate(cp%polys(i)%list(j)%y)
        deallocate(cp%polys(i)%list(j)%contains_lt)
        deallocate(cp%polys(i)%list(j)%contains_eq)
        deallocate(cp%polys(i)%list(j)%contains_gt)
        deallocate(cp%polys(i)%list(j)%borders_lt)
        deallocate(cp%polys(i)%list(j)%borders_eq)
        deallocate(cp%polys(i)%list(j)%borders_gt)
      enddo
      do j = 1, size(cp%lines(i)%list)
        deallocate(cp%lines(i)%list(j)%x)
        deallocate(cp%lines(i)%list(j)%y)
      enddo
      deallocate(cp%polys(i)%list)
      deallocate(cp%lines(i)%list)
    enddo
    deallocate(cp%polys)
    deallocate(cp%lines)
    deallocate(cp%contours)
  end subroutine destroy_contourLines

  subroutine init_contourObject(o, nx, ny, cv, zmax)
    type(contourObject), intent(inout) :: o
    integer, intent(in) :: nx, ny
    real, intent(in) :: cv(:)
    real, intent(in), optional :: zmax

    integer :: i, k, ncv

    ncv = size(cv)

    allocate(o%contours(ncv))
    o%contours = cv
    allocate(o%lines(0:ncv))
    allocate(o%polys(0:ncv))
    ! Put in bounding box:
    allocate(o%lines(0)%list(1))
    allocate(o%lines(0)%list(1)%x(2*nx+2*ny-3))
    allocate(o%lines(0)%list(1)%y(2*nx+2*ny-3))
    o%lines(0)%list(1)%x(1) = 1.0
    o%lines(0)%list(1)%y(1) = 1.0
    k = 2
    do i = 2, ny
      o%lines(0)%list(1)%x(k) = 1.0
      o%lines(0)%list(1)%y(k) = i
      k = k + 1
    enddo
    do i = 2, nx
      o%lines(0)%list(1)%x(k) = i
      o%lines(0)%list(1)%y(k) = ny
      k = k + 1
    enddo
    do i = ny - 1, 1, -1
      o%lines(0)%list(1)%x(k) = nx
      o%lines(0)%list(1)%y(k) = i
      k = k + 1
    enddo
    do i = nx - 1, 1, -1
      o%lines(0)%list(1)%x(k) = i
      o%lines(0)%list(1)%y(k) = 1.0
      k = k + 1
    enddo
    ! Everything else
    allocate(o%polys(0)%list(0))
    do i = 1, size(cv)
      allocate(o%lines(i)%list(0))
      allocate(o%polys(i)%list(0))
    enddo

    if (present(zmax)) then
      allocate(o%zmax)
      o%zmax = zmax
    endif

  end subroutine init_contourObject

  subroutine destroy_contourObject(o)
    type(contourObject), intent(inout) :: o
    integer :: i, j, k
    
    do i = 0, size(o%contours)
      do j = 1, size(o%lines(i)%list)
        deallocate(o%lines(i)%list(j)%x)
        deallocate(o%lines(i)%list(j)%y)
      enddo
      deallocate(o%lines(i)%list)
      do j = 1, size(o%polys(i)%list)
        deallocate(o%polys(i)%list(j)%outerBoundary%x)
        deallocate(o%polys(i)%list(j)%outerBoundary%y)
        do k = 1, size(o%polys(i)%list(j)%innerBoundaries%list)
          deallocate(o%polys(i)%list(j)%innerBoundaries%list(k)%x)
          deallocate(o%polys(i)%list(j)%innerBoundaries%list(k)%y)
        enddo
        deallocate(o%polys(i)%list(j)%innerBoundaries%list)
      enddo
      deallocate(o%polys(i)%list)
    enddo
    deallocate(o%lines)
    deallocate(o%polys)
    deallocate(o%contours)
    if (associated(o%zmax)) deallocate(o%zmax)
    
  end subroutine destroy_contourObject

  subroutine add_point(thisLine, x, y)
    type(richline), intent(inout) :: thisLine
    real, intent(in) :: x, y

    real, pointer :: t(:)

    ! Don't add the same point twice. Sometimes it wants to.
    if (size(thisLine%x)>0) then
      if (abs(thisLine%x(size(thisLine%x))-x)<1e-5 &
        .and.abs(thisLine%y(size(thisLine%y))-y)<1e-5) return
    endif
    
    t => thisLine%x
    allocate(thisLine%x(size(t)+1))
    thisLine%x(:size(t)) = t
    thisLine%x(size(t)+1) = x
    deallocate(t)
    t => thisLine%y
    allocate(thisLine%y(size(t)+1))
    thisLine%y(:size(t)) = t
    thisLine%y(size(t)+1) = y
    deallocate(t)

  end subroutine add_point

  subroutine add_line(plot, ci, se, ee, uphill)
    type(contourLines), intent(inout) :: plot
    integer, intent(in) :: ci, se, ee
    logical, intent(in) :: uphill

    type(richline), pointer :: tempLines(:)
    integer :: i


    ! Don't add line if there's only one point!
    if (size(plot%current%x)==1) then
      deallocate(plot%current%x)
      deallocate(plot%current%y)
      return
    endif

    tempLines => plot%lines(ci)%list
    allocate(plot%lines(ci)%list(size(tempLines)+1))
    do i = 1, size(tempLines)
      plot%lines(ci)%list(i)%x => tempLines(i)%x
      plot%lines(ci)%list(i)%y => tempLines(i)%y
      plot%lines(ci)%list(i)%se = tempLines(i)%se
      plot%lines(ci)%list(i)%ee = tempLines(i)%ee
      plot%lines(ci)%list(i)%uphill = tempLines(i)%uphill
    enddo
    deallocate(tempLines)
    plot%lines(ci)%list(i)%x => plot%current%x
    plot%lines(ci)%list(i)%y => plot%current%y
    plot%lines(ci)%list(i)%se = se
    plot%lines(ci)%list(i)%ee = ee
    plot%lines(ci)%list(i)%uphill = uphill

    nullify(plot%current%x)
    nullify(plot%current%y)
  end subroutine add_line

  subroutine shave_polygon(l)
    type(richline) :: l

    real :: x1, x2, x3, y1, y2, y3
    real, pointer :: x(:), y(:), tempx(:), tempy(:)
    integer :: i, shrink, ix1, ix2, ix3

    x => l%x; y=> l%y
    i = 1
    shrink = 1
    ! print*, "STARTING LENGTH ", size(x)
    do while (.false.)!(shrink<size(x)+3.and.size(x)>2)
      ix1 = i
      ix2 = mod(i, size(x))+1
      ix3 = mod(i+1, size(x))+1
      x1 = x(ix1); y1 = y(ix1)
      x2 = x(ix2); y2 = y(ix2)
      x3 = x(ix3); y3 = y(ix3)
      if (x2==x3.and.y2==y3) then
        shrink = 1
      elseif (x1==x3.and.y1==y3) then
        shrink = 1
        i = modulo(i-1, size(x))
      else
        i = i + 1
        if (i>size(x)) i = 1
        shrink = shrink + 1
      endif
      ! print*, shrink, size(x)
      if (shrink==1) then
        !print*, "SHRINKING"
        tempx => x
        tempy => y
        allocate(x(size(tempx)-1))
        allocate(y(size(tempx)-1))
        if (ix3==1) then
          ! Lose one item from the end of the array
          x = tempx(:ix1)
          y = tempy(:ix1)
        elseif (ix3==2) then
          ! Lose one item from the beginning of the array
          x = tempx(2:)
          y = tempy(2:)
        else
          ! Lose an item from the middle
          x(:ix1) = tempx(:ix1)
          x(ix2:) = tempx(ix3:)
          y(:ix1) = tempy(:ix1)
          y(ix2:) = tempy(ix3:)
        endif
        deallocate(tempx)
        deallocate(tempy)
      endif
    enddo
    ! print*, "ENDING LENGTH ", size(x)
    ! Stupid KML requires last & first coordinate the same, so restore that:
    tempx => x
    tempy => y
    allocate(x(size(tempx)+1))
    allocate(y(size(tempx)+1))
    x(:size(tempx)) = tempx
    y(:size(tempx)) = tempy
    x(size(x)) = x(1)
    y(size(x)) = y(1)
    deallocate(tempx)
    deallocate(tempy)

    l%x => x; l%y => y

  end subroutine shave_polygon

  subroutine shave_polygons(plot)
    type(contourLines), intent(inout) :: plot

    integer :: i, j

    do i = lbound(plot%polys, 1), ubound(plot%polys, 1)
      do j = 1, size(plot%polys(i)%list)
        call shave_polygon(plot%polys(i)%list(j))
      enddo
      ! FIXME and remove if necessary
    enddo
  end subroutine shave_polygons


  subroutine add_polygon(plot, ci)
    type(contourLines), intent(inout) :: plot
    integer, intent(in) :: ci

    type(richline), pointer :: tempLines(:)
    integer :: i
    real :: area

    ! Check polygon for wrong bits
    call shave_polygon(plot%current)

    ! Don't record empty polygons:
    if (size(plot%current%x)<4) then
      deallocate(plot%current%x)
      deallocate(plot%current%y)
      return
    endif
    
    ! Record area of polygon
    area = 0.0
    do i = 1, size(plot%current%x)-1
      area = area + plot%current%x(i)*plot%current%y(i+1) - plot%current%x(i+1)*plot%current%y(i)
    enddo
    area = area + plot%current%x(i)*plot%current%y(1) - plot%current%x(1)*plot%current%y(i)

    tempLines => plot%polys(ci)%list
    allocate(plot%polys(ci)%list(size(tempLines)+1))
    do i = 1, size(tempLines)
      plot%polys(ci)%list(i)%x => tempLines(i)%x
      plot%polys(ci)%list(i)%y => tempLines(i)%y
      plot%polys(ci)%list(i)%mountain = tempLines(i)%mountain
      plot%polys(ci)%list(i)%boundary = tempLines(i)%boundary
      plot%polys(ci)%list(i)%se = tempLines(i)%se
      plot%polys(ci)%list(i)%area = tempLines(i)%area
    enddo
    deallocate(tempLines)
    plot%polys(ci)%list(i)%x => plot%current%x
    plot%polys(ci)%list(i)%y => plot%current%y
    plot%polys(ci)%list(i)%area = area
    nullify(plot%current%x)
    nullify(plot%current%y)

  end subroutine add_polygon

  subroutine concatenate_lines(l1, l2, forward)
    type(richline), intent(inout) :: l1
    type(richline), intent(in) :: l2
    logical, intent(in) :: forward

    type(line) :: t
    
    t%x => l1%x
    t%y => l1%y
    
    allocate(l1%x(size(t%x)+size(l2%x)))
    allocate(l1%y(size(t%y)+size(l2%y)))
    l1%x(:size(t%x)) = t%x
    l1%y(:size(t%y)) = t%y
    if (forward) then
      l1%x(size(t%x)+1:) = l2%x
      l1%y(size(t%y)+1:) = l2%y
    else
      l1%x(size(t%x)+1:) = l2%x(size(l2%x):1:-1)
      l1%y(size(t%y)+1:) = l2%y(size(l2%y):1:-1)
    endif
    deallocate(t%x)
    deallocate(t%y)

    l1%uphill = l2%uphill.eqv.forward

  end subroutine concatenate_lines
      
  subroutine simpleScalePoints(x, y, thisX, thisY)
    real, intent(in) :: x(:), y(:)
    real, intent(inout) :: thisX(:), thisY(:)
    !simple linear interpolation on a rectilinear grid

    integer :: i, j
    real :: t

    do i = 1, size(thisX)
      t = thisX(i)
      if (int(t)==size(x)) then
        thisX(i) = x(size(x))
      else
        thisX(i) = x(int(t)) + (x(int(t)+1) - x(int(t))) * (t - floor(t))
      endif
      t = thisY(i)
      if (int(t)==size(y)) then
        thisY(i) = y(size(y))
      else
        thisY(i) = y(int(t)) + (y(int(t)+1) - y(int(t))) * (t - floor(t))
      endif
    end do

  end subroutine simpleScalePoints

  subroutine complexScalePoints(x, y, thisX, thisY)
    real, intent(in) :: x(:,:), y(:,:)
    real, intent(inout) :: thisX(:), thisY(:)
    ! simplest possible rescaling: linear interpolation & find where they cross.

    integer :: i, j, tx1, tx2, ty1, ty2
    real :: tx, ty, txfrac, tyfrac 
    real, dimension(2) :: bottomPoint, topPoint, bottomToTop, scaledT

    do i = 1, size(thisX)
      tx = thisX(i)
      ty = thisY(i)
      ! Mark out the square in rectilinear space within which this point exists: from (tx1, ty1) to (tx2, ty2)
      tx1 = int(tx)
      ty1 = int(ty)
      if (int(tx)==size(x,1)) then
        tx2 = tx1
      else
        tx2 = tx1 + 1
      endif
      if (int(ty)==size(y,2)) then
        ty2 = ty1
      else
        ty2 = ty1 + 1
      endif
      txfrac = tx - tx1
      tyfrac = ty - ty1
      
      ! The true coordinate box is from (x(tx1,ty1),y(tx1,ty1)) up to (x(tx2,ty2),y(tx2,ty2))
      ! Find the point along the bottom and top lines of that box which is txfrac of the way along:
      bottomPoint = (/ x(tx1,ty1) + txfrac * (x(tx2,ty1) - x(tx1,ty1)), &
        y(tx1,ty1) + txfrac * (y(tx2,ty1) - y(tx1,ty1)) /)
      topPoint = (/ x(tx1,ty2) + txfrac * (x(tx2,ty2) - x(tx1,ty2)), &
        y(tx1,ty2) + txfrac * (y(tx2,ty2) - y(tx1,ty2)) /)
      ! Draw a line between these two points, and go tyfrac of the way up the line:
      bottomToTop = topPoint - bottomPoint
      scaledT = bottomPoint + tyfrac * bottomToTop
    
      thisX(i) = scaledT(1)
      thisY(i) = scaledT(2)
    enddo

  end subroutine complexScalePoints

  subroutine scale_contours_simple(cp, x, y)
    type(contourObject), intent(inout) :: cp
    real, intent(in) :: x(:), y(:)

    integer :: i, j, k
    real, pointer :: thisX(:), thisY(:)
    do i = 0, size(cp%contours)
      do j = 1, size(cp%lines(i)%list)
        thisX => cp%lines(i)%list(j)%x
        thisY => cp%lines(i)%list(j)%y
        call simpleScalePoints(x, y, thisX, thisY)
      enddo
      do j = 1, size(cp%polys(i)%list)
        thisX => cp%polys(i)%list(j)%outerBoundary%x
        thisY => cp%polys(i)%list(j)%outerBoundary%y
        call simplescalePoints(x, y, thisX, thisY)
        do k = 1, size(cp%polys(i)%list(j)%innerBoundaries%list)
          thisX => cp%polys(i)%list(j)%innerBoundaries%list(k)%x
          thisY => cp%polys(i)%list(j)%innerBoundaries%list(k)%y
          call simpleScalePoints(x, y, thisX, thisY)
        enddo
      enddo
    enddo
  end subroutine scale_contours_simple

  subroutine scale_contours_complex(cp, x, y)
    type(contourObject), intent(inout) :: cp
    real, intent(in) :: x(:,:), y(:,:)

    integer :: i, j, k
    real, pointer :: thisX(:), thisY(:)
    do i = 0, size(cp%contours)
      do j = 1, size(cp%lines(i)%list)
        thisX => cp%lines(i)%list(j)%x
        thisY => cp%lines(i)%list(j)%y
        call complexScalePoints(x, y, thisX, thisY)
      enddo
      do j = 1, size(cp%polys(i)%list)
        thisX => cp%polys(i)%list(j)%outerBoundary%x
        thisY => cp%polys(i)%list(j)%outerBoundary%y
        call complexScalePoints(x, y, thisX, thisY)
        do k = 1, size(cp%polys(i)%list(j)%innerBoundaries%list)
          thisX => cp%polys(i)%list(j)%innerBoundaries%list(k)%x
          thisY => cp%polys(i)%list(j)%innerBoundaries%list(k)%y
          call complexScalePoints(x, y, thisX, thisY)
        enddo
      enddo
    enddo
  end subroutine scale_contours_complex

  function make_contours_on_simplest_grid(e, w, n, s, z, contour_values, ncv, ignore_gt) result(o)
    real, intent(in) :: e, w, n, s, z(:,:)
    real, intent(in), optional :: contour_values(:)
    integer, intent(in), optional :: ncv
    real, intent(in), optional :: ignore_gt
    type(contourObject) :: o
    type(contourLines) :: cp

    integer :: i, ncv_
    real :: x(size(z,1)), y(size(z,2))
    real :: zmin, zmax, zinc
    real, allocatable :: contour_values_(:)

    if (present(contour_values)) then
      if (present(ncv)) then
        call FoX_error('cannot specify ncv and contour_values at the same time')
      endif
      do i = 2, size(contour_values)
        if (contour_values(i)<=contour_values(i-1)) then
          call FoX_error('contour values must be monotonically increasing')
        endif
      enddo
      if (present(ignore_gt)) then
        if (contour_values(size(contour_values))>ignore_gt) then
          call FoX_error('cannot have contours above the maximum value')
        endif
        allocate(contour_values_(size(contour_values)+1))
        contour_values_(:size(contour_values)) = contour_values
        contour_values_(size(contour_values_)) = contour_values(size(contour_values)) + &
          0.999*(ignore_gt-contour_values(size(contour_values)))
      else
        allocate(contour_values_(size(contour_values)))
        contour_values_ = contour_values
      endif
    else
      if (present(ncv)) then
        ncv_ = ncv
      else
        ncv_ = 5
      endif
      zmin = minval(z)
      if (present(ignore_gt)) then
        allocate(contour_values_(ncv_+1))
        zmax = minval(abs(z-ignore_gt))
      else
        allocate(contour_values_(ncv_))
        zmax = maxval(z)
      endif
      zinc = (zmax-zmin)/(ncv_+1)
      do i = 1, ncv_
        contour_values_(i) = zmin + i * zinc 
      enddo
      if (present(ignore_gt)) then
        contour_values_(i) = zmin + (i-0.001)*zinc
      endif
    endif

    ! Note that we need extra poly lists - to surround areas which are less than the minimum contour, or larger than the greatest
    call gcontr(z, contour_values_, o, ignoreval=ignore_gt)

    do i = 1, size(z,1)
      x(i) = w + (i-1)*(e-w)/(size(z,1)-1)
    enddo
    do i = 1, size(z,2)
      y(i) = n - (i-1)*(n-s)/(size(z,2)-1)
    enddo

    call scale_contours_simple(o, x, y)
  end function make_contours_on_simplest_grid

  function make_contours_on_a_simple_grid(x, y, z, contour_values, ncv, ignore_gt) result(o)
    real, intent(in) :: x(:), y(:), z(:,:)
    real, intent(in), optional :: contour_values(:)
    integer, intent(in), optional :: ncv
    real, intent(in), optional :: ignore_gt
    type(contourObject) :: o
    type(contourLines) :: cp

    integer :: i, ncv_
    real :: zmin, zmax, zinc
    real, allocatable :: contour_values_(:)

    if (size(x)/=size(z,1)) then
      call FoX_error('wrong number of x coordinates')
    elseif (size(y)/=size(z,2)) then
      call FoX_error('wrong number of x coordinates')
    endif

    if (present(contour_values)) then
      if (present(ncv)) then
        call FoX_error('cannot specify ncv and contour_values at the same time')
      endif
      do i = 2, size(contour_values)
        if (contour_values(i)<=contour_values(i-1)) then
          call FoX_error('contour values must be monotonically increasing')
        endif
      enddo
      if (present(ignore_gt)) then
        if (contour_values(size(contour_values))>ignore_gt) then
          call FoX_error('cannot have contours above the maximum value')
        endif
        allocate(contour_values_(size(contour_values)+1))
        contour_values_(:size(contour_values)) = contour_values
        contour_values_(size(contour_values_)) = contour_values(size(contour_values)) + &
          0.999*(ignore_gt-contour_values(size(contour_values)))
      else
        allocate(contour_values_(size(contour_values)))
        contour_values_ = contour_values
      endif
    else
      if (present(ncv)) then
        ncv_ = ncv
      else
        ncv_ = 5
      endif
      zmin = minval(z)
      if (present(ignore_gt)) then
        allocate(contour_values_(ncv_+1))
        zmax = minval(abs(z-ignore_gt))
      else
        allocate(contour_values_(ncv_))
        zmax = maxval(z)
      endif
      zinc = (zmax-zmin)/(ncv_+1)
      do i = 1, ncv_
        contour_values_(i) = zmin + i * zinc 
      enddo
      if (present(ignore_gt)) then
        contour_values_(i) = zmin + (i-0.001)*zinc
      endif
    endif

    ! Note that we need extra poly lists - to surround areas which are less than the minimum contour, or larger than the greatest
    call gcontr(z, contour_values_, o, ignoreval=ignore_gt)

    call scale_contours_simple(o, x, y)
  end function make_contours_on_a_simple_grid

  function make_contours_on_a_complex_grid(x, y, z, contour_values, ncv, ignore_gt) result(o)
    real, intent(in) :: x(:,:), y(:,:), z(:,:)
    real, intent(in), optional :: contour_values(:)
    integer, intent(in), optional :: ncv
    real, intent(in), optional :: ignore_gt
    type(contourObject) :: o
    type(contourLines) :: cp

    integer :: i, ncv_
    real :: zmin, zmax, zinc
    real, allocatable :: contour_values_(:)

    if (size(x,1)/=size(z,1).and.size(y,1)/=size(z,2)) then
      call FoX_error('wrong number of coordinates in the x direction')
    elseif (size(x,2)/=size(z,2).and.size(y,2)/=size(z,2)) then
      call FoX_error('wrong number of coordinates in the y direction')
    endif

    if (present(contour_values)) then
      if (present(ncv)) then
        call FoX_error('cannot specify ncv and contour_values at the same time')
      endif
      do i = 2, size(contour_values)
        if (contour_values(i)<=contour_values(i-1)) then
          call FoX_error('contour values must be monotonically increasing')
        endif
      enddo
      if (present(ignore_gt)) then
        if (contour_values(size(contour_values))>ignore_gt) then
          call FoX_error('cannot have contours above the maximum value')
        endif
        allocate(contour_values_(size(contour_values)+1))
        contour_values_(:size(contour_values)) = contour_values
        contour_values_(size(contour_values_)) = ignore_gt
      else
        allocate(contour_values_(size(contour_values)))
        contour_values_ = contour_values
      endif
    else
      if (present(ncv)) then
        ncv_ = ncv
      else
        ncv_ = 5
      endif
      zmin = minval(z)
      if (present(ignore_gt)) then
        allocate(contour_values_(ncv_+1))
        zmax = minval(abs(z-ignore_gt))
      else
        allocate(contour_values_(ncv_))
        zmax = maxval(z)
      endif
      zinc = (zmax-zmin)/(ncv_+1)
      do i = 1, ncv_
        contour_values_(i) = zmin + i * zinc 
      enddo
      if (present(ignore_gt)) then
        contour_values_(i) = zmin + (i-0.001)*zinc
      endif
    endif

    ! Note that we need extra poly lists - to surround areas which are less than the minimum contour, or larger than the greatest
    call gcontr(z, contour_values_, o, ignoreval=ignore_gt)

    call scale_contours_complex(o, x, y)
  end function make_contours_on_a_complex_grid

  subroutine draw(x, y, nx, ny, icv, iflag, cp, uphill)
    real, intent(in) :: x, y
    integer, intent(in) :: nx, ny
    integer, intent(in) :: icv
    integer, intent(in) :: iflag
    type(contourLines), intent(inout) :: cp
    logical, intent(in) :: uphill

    select case (iflag)
    case (1)
      ! carry on along the way ...
      call add_point(cp%current, x, y)
    case (2)
      if (x==1.0) then
        cp%current_se = 1
      else if (y==ny) then
        cp%current_se = 2
      else if (x==nx) then
        cp%current_se = 3
      else if (y==1) then
        cp%current_se = 4
      else
        cp%current_se = 0
      endif
      call init_line(cp%current)
      call add_point(cp%current, x, y)
    case (3)
      cp%current_se = 0
      call init_line(cp%current)
      call add_point(cp%current, x, y)
    case (4)
      call add_point(cp%current, x, y)
      if (cp%current_se==0) then
        call add_polygon(cp, icv)
      else
        if (x==1) then
          cp%current_ee = 1
        else if (y==ny) then
          cp%current_ee = 2
        else if (x==nx) then
          cp%current_ee = 3
        else if (y==1) then
          cp%current_ee = 4
        else
          cp%current_ee = 0
        endif
        call add_line(cp, icv, cp%current_se, cp%current_ee, uphill)
      endif
    case (5) 
      call add_point(cp%current, x, y)
      call add_polygon(cp, icv)
      ! If we were going uphill, and clockwise ...
      if (uphill.eqv.(cp%polys(icv)%list(size(cp%polys(icv)%list))%area<0.0)) then
        cp%polys(icv)%list(size(cp%polys(icv)%list))%mountain = 2
      else
        cp%polys(icv)%list(size(cp%polys(icv)%list))%mountain = 1
      endif
    end select

  end subroutine draw

  SUBROUTINE GCONTR(Z, CV, o, ignoreval)

    !     THIS SUBROUTINE DRAWS A CONTOUR THROUGH EQUAL VALUES OF AN ARRAY.
    !
    !     *****     FORMAL ARGUMENTS     ***********************************
    !
    !     Z IS THE ARRAY FOR WHICH CONTOURS ARE TO BE DRAWN.  THE ELEMENTS
    !     OF Z ARE ASSUMED TO LIE UPON THE NODES OF A TOPOLOGICALLY
    !     RECTANGULAR COORDINATE SYSTEM - E.G. CARTESIAN, POLAR (EXCEPT
    !     THE ORIGIN), ETC.
    !
    !     NRZ IS THE NUMBER OF ROWS DECLARED FOR Z IN THE CALLING PROGRAM.
    !
    !     NX IS THE LIMIT FOR THE FIRST SUBSCRIPT OF Z.
    !
    !     NY IS THE LIMIT FOR THE SECOND SUBSCRIPT OF Z.
    !
    !     CV ARE THE VALUES OF THE CONTOURS TO BE DRAWN.
    !
    !     NCV IS THE NUMBER OF CONTOUR VALUES IN CV.
    !
    !     ZMAX IS THE MAXIMUM VALUE OF Z FOR CONSIDERATION.  A VALUE OF
    !     Z(I,J) GREATER THAN ZMAX IS A SIGNAL THAT THAT POINT AND THE
    !     GRID LINE SEGMENTS RADIATING FROM THAT POINT TO IT'S NEIGHBORS
    !     ARE TO BE EXCLUDED FROM CONTOURING.
    !

    !
    !     DRAW IS A USER-PROVIDED SUBROUTINE USED TO DRAW CONTOURS.
    !     THE CALLING SEQUENCE FOR DRAW IS:
    !
    !         CALL DRAW (X,Y,IFLAG)
    !         LET NX = INTEGER PART OF X, FX = FRACTIONAL PART OF X.
    !         THEN X SHOULD BE INTERPRETED SUCH THAT INCREASES IN NX
    !         CORRESPOND TO INCREASES IN THE FIRST SUBSCRIPT OF Z, AND
    !         FX IS THE FRACTIONAL DISTANCE FROM THE ABSCISSA CORRESPONDING
    !         TO NX TO THE ABSCISSA CORRESPONDING TO NX+1,
    !         AND Y SHOULD BE INTERPRETED SIMILARLY FOR THE SECOND
    !         SUBSCRIPT OF Z.
    !         THE LOW-ORDER DIGIT OF IFLAG WILL HAVE ONE OF THE VALUES:
    !             1 - CONTINUE A CONTOUR,
    !             2 - START A CONTOUR AT A BOUNDARY,
    !             3 - START A CONTOUR NOT AT A BOUNDARY,
    !             4 - FINISH A CONTOUR AT A BOUNDARY,
    !             5 - FINISH A CLOSED CONTOUR (NOT AT A BOUNDARY).
    !                 NOTE THAT REQUESTS 1, 4 AND 5 ARE FOR PEN-DOWN
    !                 MOVES, AND THAT REQUESTS 2 AND 3 ARE FOR PEN-UP
    !                 MOVES.
    !             6 - SET X AND Y TO THE APPROXIMATE 'PEN' POSITION, USING
    !                 THE NOTATION DISCUSSED ABOVE.  THIS CALL MAY BE
    !                 IGNORED, THE RESULT BEING THAT THE 'PEN' POSITION
    !                 IS TAKEN TO CORRESPOND TO Z(1,1).
    !         IFLAG/10 IS THE CONTOUR NUMBER.
    !
    !     *****     EXTERNAL SUBPROGRAMS     *******************************
    !
    !     DRAW IS THE USER-SUPPLIED LINE DRAWING SUBPROGRAM DESCRIBED ABOVE.
    !     DRAW MAY BE SENSITIVE TO THE HOST COMPUTER AND TO THE PLOT DEVICE.
    !
    !     ******************************************************************
    !
    REAL, intent(in) :: Z(:,:)
    REAL, intent(in) :: CV(:)
    type(contourObject), intent(out) :: o
    real, intent(in), optional :: ignoreval
    ! If we encounter any values greater than ignoreval, we ignore them, and do
    ! not draw contours through the region
    ! This lets us do contour plots over incomplete data sets - set all incomplete
    ! data to be greater than ignoreval

    real :: x, y
    integer :: i, j, imin, imax, jmin, jmax
    real :: dmax, cval
    integer :: ibkey
    integer :: icv, icur, iflag, ii, iedge, ix, jcur, jump, jj, ks, k, l, idir, nxidir, ni, m, n
    real :: zz, z1, z2

    integer :: nx, ny, ncv
    integer, allocatable :: tempList(:)

    integer :: se, ee

    INTEGER L1(4), L2(4), IJ(2)
    !
    !     L1 AND L2 CONTAIN LIMITS USED DURING THE SPIRAL SEARCH FOR THE
    !     BEGINNING OF A CONTOUR.
    !     IJ STORES SUBCRIPTS USED DURING THE SPIRAL SEARCH.
    !
    INTEGER I1(2), I2(2), I3(6)
    !
    !     I1, I2 AND I3 ARE USED FOR SUBSCRIPT COMPUTATIONS DURING THE
    !     EXAMINATION OF LINES FROM Z(I,J) TO IT'S NEIGHBORS.
    !
    REAL XINT(4)
    !
    !     XINT IS USED TO MARK INTERSECTIONS OF THE CONTOUR UNDER
    !     CONSIDERATION WITH THE EDGES OF THE CELL BEING EXAMINED.
    !
    REAL XY(2)
    !
    !     XY IS USED TO COMPUTE COORDINATES FOR THE DRAW SUBROUTINE.
    !

    logical :: contourTest(2, size(z,1), size(z,2), size(cv))
    ! 2 for vertical/horizontal
    ! 1 for each point (though we don't need the last surely?)
    ! 1 for each contour value
    logical :: uphill, uphills(4)
    ! need to keep track of boundaries as well:
    integer :: boundary1, boundary2
    real :: highpoint(2), firstpoint(2)

    type(contourLines) :: cp

    logical jump_to ! To avoid the assigned goto deleted from 
                    ! Fortran 95 - if .false. goto 100
                    !              if .true. goto 280...

    EQUIVALENCE (L2(1),IMAX), (L2(2),JMAX), (L2(3),IMIN), (L2(4),JMIN)
    EQUIVALENCE (IJ(1),I), (IJ(2),J)
    EQUIVALENCE (XY(1),X), (XY(2),Y)

    !
    !  Initialize contour plot
    nx = size(z, 1)
    ny = size(z, 2)
    ncv = size(cv)
    L1 = (/nx, ny, -1, -1/)
    I1 = (/1, 0/)
    I2 = (/1,-1/)
    I3 = (/1,0,0,1,1,0/)

    allocate(cp%contours(ncv))
    cp%contours = cv
    allocate(cp%polys(0:ncv))
    allocate(cp%lines(0:ncv))
    do i = 0, ncv
      allocate(cp%polys(i)%list(0))
      allocate(cp%lines(i)%list(0))
    enddo

    !if (present(ignoreval)) then
    !  dmax = ignoreval
    !else
      DMAX = huge(z)
    !endif
    !
    !     SET THE CURRENT PEN POSITION.  THE DEFAULT POSITION CORRESPONDS
    !     TO Z(1,1).
    !
    X = 1.0
    Y = 1.0
    ICUR = 1
    JCUR = 1
    !
    !     CLEAR THE BITMAP
    !
    contourTest = .false.
    !
    !     SEARCH ALONG A RECTANGULAR SPIRAL PATH FOR A LINE SEGMENT HAVING
    !     THE FOLLOWING PROPERTIES:
    !          1.  THE END POINTS ARE NOT EXCLUDED,
    !          2.  NO MARK HAS BEEN RECORDED FOR THE SEGMENT,
    !          3.  THE VALUES OF Z AT THE ENDS OF THE SEGMENT ARE SUCH THAT
    !              ONE Z IS LESS THAN THE CURRENT CONTOUR VALUE, AND THE
    !              OTHER IS GREATER THAN OR EQUAL TO THE CURRENT CONTOUR
    !              VALUE.
    !
    !     SEARCH ALL BOUNDARIES FIRST, THEN SEARCH INTERIOR LINE SEGMENTS.
    !     NOTE THAT THE INTERIOR LINE SEGMENTS NEAR EXCLUDED POINTS MAY BE
    !     BOUNDARIES.
    !
    IBKEY = 0
10  I = ICUR
    J = JCUR
20  IMAX = I
    IMIN = -I
    JMAX = J
    JMIN = -J
    IDIR = 0
    !     DIRECTION ZERO IS +I, 1 IS +J, 2 IS -I, 3 IS -J.
30  NXIDIR = IDIR + 1
    K = NXIDIR
    IF (NXIDIR.GT.3) NXIDIR = 0
40  I = IABS(I)
    J = IABS(J)
    IF (Z(I,J).GT.DMAX) GO TO 140
    L = 1
    !     L=1 MEANS HORIZONTAL LINE, L=2 MEANS VERTICAL LINE.
50  IF (IJ(L).GE.L1(L)) GO TO 130 
    II = I + I1(L)
    JJ = J + I1(3-L)
    IF (Z(II,JJ).GT.DMAX) GO TO 130
    ! ASSIGN 100 TO JUMP
    jump_to = .false.
    !     THE NEXT 15 STATEMENTS (OR SO) DETECT BOUNDARIES.
60  IX = 1
    IF (IJ(3-L).EQ.1) GO TO 80
    II = I - I1(3-L)
    JJ = J - I1(L)
    IF (Z(II,JJ).GT.DMAX) GO TO 70
    II = I + I2(L)
    JJ = J + I2(3-L)
    IF (Z(II,JJ).LT.DMAX) IX = 0
70  IF (IJ(3-L).GE.L1(3-L)) GO TO 90
80  II = I + I1(3-L)
    JJ = J + I1(L)
    IF (Z(II,JJ).GT.DMAX) GO TO 90
    IF (Z(I+1,J+1).LT.DMAX) then ! GO TO JUMP, (100, 280)
      if (jump_to) then
        go to 280
      else 
        go to 100
      endif
    endif
90  IX = IX + 2
    !GO TO JUMP, (100, 280)
    if (jump_to) then
      go to 280
    else 
      go to 100
    endif
100 IF (IX.EQ.3) GO TO 130
    IF (IX+IBKEY.EQ.0) GO TO 130
    !     NOW DETERMINE WHETHER THE LINE SEGMENT IS CROSSED BY THE CONTOUR.
    II = I + I1(L)
    JJ = J + I1(3-L)
    ! The two points of interest are:
    Z1 = Z(I,J)
    Z2 = Z(II,JJ)
    DO 120 ICV = 1, NCV
      if (contourTest(l, i, j, icv)) goto 120 ! We've checked this edge before for this contour
      IF (CV(ICV).LE.MIN(Z1,Z2)) GO TO 110 ! Nope, not this one.
      IF (CV(ICV).LE.MAX(Z1,Z2)) then
        GO TO 190 ! Yes, the contour crosses here.
      endif
110   continue
      contourTest(l, i, j, icv) = .true.
120 continue
130 L = L + 1
    IF (L.LE.2) GO TO 50 ! search vertically instead of horizontally
140 L = MOD(IDIR,2) + 1 
    IJ(L) = ISIGN(IJ(L),L1(K))
    !
    !     LINES FROM Z(I,J) TO Z(I+1,J) AND Z(I,J+1) ARE NOT SATISFACTORY.
    !     CONTINUE THE SPIRAL.
    !
150 IF (IJ(L).GE.L1(K)) GO TO 170
    IJ(L) = IJ(L) + 1
    IF (IJ(L).GT.L2(K)) GO TO 160
    GO TO 40
160 L2(K) = IJ(L)
    IDIR = NXIDIR
    GO TO 30
170 IF (IDIR.EQ.NXIDIR) GO TO 180
    NXIDIR = NXIDIR + 1
    IJ(L) = L1(K)
    K = NXIDIR
    L = 3 - L
    IJ(L) = L2(K)
    IF (NXIDIR.GT.3) NXIDIR = 0
    GO TO 150
180 IF (IBKEY.NE.0) goto 300
    IBKEY = 1
    GO TO 10
    !
    !     AN ACCEPTABLE LINE SEGMENT HAS BEEN FOUND.
    !     FOLLOW THE CONTOUR UNTIL IT EITHER HITS A BOUNDARY OR CLOSES.
    !
190 continue ! At this point, z1 and z2 hold the two points of interest
    IEDGE = L ! horizontal or vertical?
    CVAL = CV(ICV) ! contour value of interest
    IF (IX.NE.1) IEDGE = IEDGE + 2 ! ... ?
    IFLAG = 2 + IBKEY ! We've only just found one, so this is the start of a line at either a boundary or inside. (according to ibkey)
    XINT(IEDGE) = (CVAL-Z1)/(Z2-Z1) ! Linear interpolation to find fraction in the y-direction
200 XY(L) = FLOAT(IJ(L)) + XINT(IEDGE) ! X (or Y, according to L) is the value at current ij +- the interpolate along the line
    XY(3-L) = FLOAT(IJ(3-L)) ! and the other one is constant of course
    contourTest(l, i, j, icv) = .true.
    call draw(x, y, nx, ny, icv, iflag, cp, uphill)
    IF (IFLAG.LT.4)  then ! we need to carry on with the present line ...
      GO TO 210
    endif
    ! Otherwise we've just finished a line, so we go back to the beginning:
    ICUR = I
    JCUR = J
    GO TO 20
    !
    !     CONTINUE A CONTOUR.  THE EDGES ARE NUMBERED CLOCKWISE WITH
    !     THE BOTTOM EDGE BEING EDGE NUMBER ONE.
    !
210 NI = 1
    IF (IEDGE.LT.3) GO TO 220
    I = I - I3(IEDGE)
    J = J - I3(IEDGE+2)
220 DO K=1,4
      IF (K.EQ.IEDGE) cycle
      II = I + I3(K)
      JJ = J + I3(K+1)
      Z1 = Z(II,JJ)
      II = I + I3(K+1)
      JJ = J + I3(K+2)
      Z2 = Z(II,JJ)
      IF (CVAL.LE.MIN(Z1,Z2)) cycle
      IF (CVAL.GT.MAX(Z1,Z2)) cycle
      uphills(k) = z2>z1
      IF (K.EQ.1) GO TO 230
      IF (K.NE.4) GO TO 240
230   ZZ = Z1
      Z1 = Z2
      Z2 = ZZ
240   XINT(K) = (CVAL-Z1)/(Z2-Z1)
      NI = NI + 1
      KS = K
    enddo
    IF (NI.EQ.2) GO TO 260
    !
    !     THE CONTOUR CROSSES ALL FOUR EDGES OF THE CELL BEING EXAMINED.
    !     CHOOSE THE LINES TOP-TO-LEFT AND BOTTOM-TO-RIGHT IF THE
    !     INTERPOLATION POINT ON THE TOP EDGE IS LESS THAN THE INTERPOLATION
    !     POINT ON THE BOTTOM EDGE.  OTHERWISE, CHOOSE THE OTHER PAIR.  THIS
    !     METHOD PRODUCES THE SAME RESULTS IF THE AXES ARE REVERSED.  THE
    !     CONTOUR MAY CLOSE AT ANY EDGE, BUT MUST NOT CROSS ITSELF INSIDE
    !     ANY CELL.
    !
    KS = 5 - IEDGE
    IF (XINT(3).LT.XINT(1)) GO TO 260
    KS = 3 - IEDGE
    IF (KS.LE.0) KS = KS + 4
    !
    !     DETERMINE WHETHER THE CONTOUR WILL CLOSE OR RUN INTO A BOUNDARY
    !     AT EDGE KS OF THE CURRENT CELL.
    !
260 L = KS
    IFLAG = 1
    ! ASSIGN 280 TO JUMP
    jump_to = .true.
    IF (KS.LT.3) GO TO 270
    I = I + I3(KS)
    J = J + I3(KS+2)
    L = KS - 2
270 continue
    if (.not.contourTest(l, i, j, icv)) goto 60
    IFLAG = 5
    GO TO 290
280 IF (IX.NE.0) IFLAG = 4
290 IEDGE = KS + 2
    IF (IEDGE.GT.4) IEDGE = IEDGE - 4
    XINT(IEDGE) = XINT(KS)
    uphill = uphills(ks)
    GO TO 200

300 continue
    call init_contourObject(o, nx, ny, cv, ignoreval)

    call addLinesToContourObject(cp, o)

    if (size(cp%lines)>1) then
      call joinUpLines(cp, nx, ny)
    else
      ! We need to add the box perimeter
      ! If we only add the corners, then Google Earth seemsto get a bit unhappy,
      ! and doesn't join them up nicely ...
      call init_line(cp%current)
      do i = 2, ny
        call add_point(cp%current, 1.0, real(i))
      enddo
      do i = 2, nx
        call add_point(cp%current, real(i), real(ny))
      enddo
      do i = ny - 1, 1, -1
        call add_point(cp%current, real(nx), real(i))
      enddo
      do i = nx - 1, 1, -1
        call add_point(cp%current, real(i), 1.0)
      enddo
      ! Which level should this be at? Well, whatever level (1,1) should be at.
      k = whichLevel(z(1,1), cv)
      call add_polygon(cp, whichLevel(z(1,1), cv))
    endif
            
    call polygonsInPolygons(cp, z)

    !call shave_polygons(cp)

    call addPolygonsToContourObject(cp, o)

    call destroy_contourLines(cp)

  end subroutine gcontr

  function whichLevel(z, cv) result(level)
    real, intent(in) :: z
    real, intent(in) :: cv(:)
    integer :: level
    
    integer :: i
    if (z<cv(1)) then
      level = 0
      return
    endif
    do i = 1, size(cv)-1
      if (z>cv(i) .and. z<cv(i+1)) then
        level = i
        return
      endif
    enddo
    level = i + 1
  end function whichLevel
  
  function checkPointOnPolygon(l, x, y) result(p)
    type(richline), intent(in) :: l
    real, intent(in) :: x, y
    logical :: p
    integer :: i, j

    p = .false.
    do i = 1, size(l%x)
      if (abs(l%x(i)-x)<0.01 .and. abs(l%y(i)-y)<0.01) then
        p = .true.
        return
!        print*, 'point nearly on polygon:'
!        print*, 'point', x,y
!        print*, 'polygon', l%x(i), l%y(i)
      endif
    enddo

  end function checkPointOnPolygon

  function pointInPolygon(l, x, y) result(p)
    type(richline), intent(in) :: l
    real, intent(in) :: x, y
    logical :: p
    integer :: i, j
    
    p = .false.
    j = 1
    do i = 1, size(l%x)
      j = j + 1
      if (j>size(l%x)) j = 1
      if ((l%y(i) < y .and. l%y(j) >= y) .or. (l%y(j)<y .and. l%y(i) >= y)) then
        if ( l%x(i) &
          +(y-l%y(i))/(l%y(j)-l%y(i))*(l%x(j)-l%x(i)) &
          < x ) then
          p = .not. p
        endif
      endif
    enddo
  end function pointInPolygon

  function nextLineOnEdge(list, current, edge, minLine, forward, lastlinelength, lastse, whichhill) &
    result(minGap)
    type(richline), intent(in) :: list(:)
    type(richline), intent(in) :: current
    integer, intent(in) :: edge ! which edge are we on?
    integer, intent(inout) :: minLine ! on input, current j so we know not to find the same line (-1 if a different level)
                                      ! on output, the j of the line that has been found.
    logical, intent(inout) :: forward !on entrance, which way round is the current line.
                                      ! on exit, which way round is the line we found?
    integer, intent(in) :: lastlinelength ! what is the length of the last line segment we added to current?
    integer, intent(in) :: lastse ! on which edge started the last line segment we added to current?
    logical, intent(in) :: whichhill ! are we searching in the correct (up|down)hill direction?
    real :: minGap

    real :: p, px, py, pnew, pold, pnewnew, area, minArea
    integer :: i, j, ix, ix2, c_j, b, c, thisse
    logical :: c_f

    c_j = minLine
    c_f = forward

    ! First work out where we are starting from:
    px = current%x(size(current%x))
    py = current%y(size(current%x))
    if (edge==1.or.edge==3) then
      p = py
    elseif (edge==2.or.edge==4) then
      p = px
    endif

    minGap = -1.0
    minArea = 1.0
    minLine = 0
    ! Now loop over all lines in the list we are looking at, and check (with setpnew) how
    ! close they are to the current endpoint.
    do j = 1, size(list)
      ! We are allowed to find only the end of the line in the direction we did not start with!
      ! This means that if a) we are comparing lines on the same level
      !                    b) the line index is the same (these are both checked by j==c_j)
      !                    c) if this line goes forward the initial point of the line is the same as our finish point,
      !                       if this line goes backward, the final point of the line is the same as our finish point.
      !                 then it is the line we are coming from, we should not use it.
      if (debug) print*,"DEBUG going to check ", j, c_j, list(j)%se==edge, list(j)%ee==edge, j==c_j, c_f
      if (debug) print*, "DEBUG starts at ", list(j)%x(1), list(j)%y(1)
      if (debug) print*, "DEBUG ends at ", list(j)%x(size(list(j)%x)), list(j)%y(size(list(j)%x))
      ! For each line, check its starting point:
      if (list(j)%se==edge.and..not.(j==c_j &
        .and.list(j)%x(1)==px.and.list(j)%y(1)==py)) &
        then
        ix = 1
        ix2 = size(list(j)%x)
        call setpnew(.true.)
      endif
      ! And it's ending point.
      if (list(j)%ee==edge.and..not.(j==c_j &
        .and.list(j)%x(size(list(j)%x))==px.and.list(j)%y(size(list(j)%x))==py)) &
        then
        ix = size(list(j)%x)
        ix2 = 1
        call setpnew(.false.)
      endif
    enddo

  contains
    subroutine setpnew(whichWay)
      logical, intent(in) :: whichWay
      integer :: lastix
      lastix = size(current%x)-lastlinelength+1
      select case(edge)
      case(1)
        pold = current%y(lastix)
        pnew = list(j)%y(ix)
        pnewnew = list(j)%y(ix2)
      case(2)
        pold = current%x(lastix)
        pnew = list(j)%x(ix)
        pnewnew = list(j)%x(ix2)
      case(3)
        pold = 2.0*p - current%y(lastix)
        pnew = 2.0*p - list(j)%y(ix)
        pnewnew = 2.0*p - list(j)%y(ix2)
      case(4)
        pold = 2.0*p - current%x(lastix)
        pnew = 2.0*p - list(j)%x(ix)
        pnewnew = 2.0*p - list(j)%x(ix2)
      end select
      ! Next if statement does the following:
      ! Where this line hits the edge at the current endpoint, then the distance
      ! (pnew-p) is zero. However, this is only the closest line for our purposes
      ! if the area between the two lines is > 0; ie the line heads off pointing
      ! further along our tour round the perimeter.
      ! Indeed, in case there is more than one line present for which this
      ! is the case, we need to pick the line for which the area is the
      ! smallest positive value available.
      if (debug) print*, "DEBUG pp", p, pnew
      if (pnew>p .and. ((pnew-p)<minGap .or. minLine==0)) then
        ! Also, we want to record this as minimum only if - we have had no minimum
        ! before (so minLine==0) or it is a lower gap than previously recorded.
        if (debug) print*, "DEBUG best so far"
        minGap = pnew - p
        minLine = j
        forward = whichWay
      elseif (pnew==p) then
        ! compare the starting and finishing edges ...
        if (whichWay) then
          thisse = list(j)%ee
        else
          thisse = list(j)%se
        endif
        b = modulo(thisse - lastse, 4)
        c = modulo(edge - lastse, 4)
        if (debug) print*, "DEBUG bc", lastse, thisse, edge, b, c, p, pold, pnewnew
        if (b==0.and.c>0) then
          if (debug) print*, "DEBUG linetype 1"
          ! we started on a different side, and the candidate line will finish
          ! on the same, other, side.
          ! To know if we want it we must calculate the area of the polygon
          ! inscribed by the two lines. If it is negative, then it goes the
          ! right way round.
          area = 0.0
          do i = size(current%x)-lastlinelength+1, size(current%x)-1
            area = area + current%x(i)*current%y(i+1)-current%x(i+1)*current%y(i)
            if (debug) print*, "DEBUG a", current%x(i), current%y(i), area
          enddo
          if (whichWay) then
            area = area + current%x(i)*list(j)%y(1) - list(j)%x(1)*current%y(i)
            if (debug) print*, "DEBUG a", current%x(i), current%y(i), area
            do i = 1, size(list(j)%x)-1
              area = area + list(j)%x(i)*list(j)%y(i+1) - list(j)%x(i+1)*list(j)%y(i)
              if (debug) print*, "DEBUG a", list(j)%x(i), list(j)%y(i), area
            enddo
          else
            if (debug) print*, "DEBUG a", current%x(i), current%y(i), area
            area = area + current%x(i)*list(j)%y(size(list(j)%x)) - list(j)%x(size(list(j)%x))*current%y(i)
            do i = size(list(j)%x), 2, -1
              area = area + list(j)%x(i)*list(j)%y(i-1) - list(j)%x(i-1)*list(j)%y(i)
              if (debug) print*, "DEBUG a", list(j)%x(i), list(j)%y(i), area
            enddo
          endif
          area = area + list(j)%x(i)*current%y(size(current%x)-lastlinelength+1) &
            - current%x(size(current%x)-lastlinelength+1)*list(j)%y(i)
          if (debug) print*, "DEBUG area", current%x(size(current%x)-lastlinelength+1), &
            current%y(size(current%x)-lastlinelength+1), area
        elseif (b<c) then
          if (debug) print*, "DEBUG linetype 2"
          area = 1.0
          ! Then this line is not the one we want, mark its area positive
        elseif (b>c) then
          if (debug) print*, "DEBUG linetype 3"
          ! Then this is a genuine candidate line, mark its area negative
          area = -1.0
          ! FIXME we might need to calculate this area correctly in some cases ...
        elseif (c==0) then
          if (debug) print*, "DEBUG linetype 4"
          ! Last line started on this edge, and the candidate line finishes there
          ! too. Are they correctly ordered ...
          if ((pold<p .and. p<pnewnew) &
            .or. (pnewnew<=pold .and. pold<p) &
            .or. (p<pnewnew .and. pnewnew<=pold)) then
            area = -1.0
          else
            area = 1.0
          endif
        else
          if (debug) print*, "DEBUG linetype 5"
          ! Last line started on another edge, but the candidate line will finish on
          ! this edge. We only want it if it goes clockwise ...
          if (pnewnew>p) then
            area = -1.0
          else
            area = 1.0
          endif
        endif
        if (debug) print*, "DEBUG same line", area, c_j, whichhill
        if ((abs(area)<1e-3 .and. whichhill) &
          !These are essentially the same line, but on different levels
          ! (if it were on the same level, it could only be the same line
          ! & that is excluded above)
          ! then we only want this if we are changing levels in the correct direction.
          .or. area<0.0) &
          ! The two lines are separated in the right direction
          then
          minGap = 0.0
          if (debug) print*, "DEBUG checked area", area, minArea
          if (minArea==1.0 .or. area>minArea) then
            if (debug) print*, "DEBUG best so far"
            minArea = area
            minLine = j
            forward = whichWay
          endif
        endif
      endif
    end subroutine setpnew

  end function nextLineOnEdge

  subroutine joinUpLines(cp, nx, ny)
    type(contourLines), intent(inout) :: cp
    integer, intent(in) :: nx, ny

    type(richline), pointer :: l
    integer :: i, ii, j, jj, ix, jjj

    integer :: s, ce, once, lastlinelength, lastse
    integer :: minLine, minLine1, minLine2, minLevel
    logical :: forward, forward1, forward2, allOneLevel
    real :: x, y, minGap1, minGap2

    logical :: infinite1, infinite2, infinite3, uphill

    ! This is the algorithm:
    ! For each line on each contour:
    ! Is (se|ee)done set? If so then this side of the line is already part of a polygon started on the level below.
    ! Start a new polygon with the points of the current line.
    ! Start at the end of the line. Follow the edge of the box counterclockwise
    !    until we hit the end of another line. (finding the line is complex - see nextLineOnEdge above)
    ! Concatenate the new line onto our polygon. Set its (se|ee)done variable so we know we've counted it.
    ! Continue round until we hit the other end of the line we started with.
    ! Add the completed polygon to the list of polygons
    ! Repeat, starting from the *end* of the line instead.

    do i = 1, size(cp%contours)
      ! At each contour level:
      j = 1
      allOneLevel = .true.
      do j = 1, size(cp%lines(i)%list)
        if (cp%lines(i)%list(j)%se==0) then
          cycle
        elseif (cp%lines(i)%list(j)%ee==0) then
          cycle
        endif
        do s = 1, 2 ! 1 = start at beginning of line, 2 = start at end
          !(we need to do two sides of each line)
          ! For each line discovered at this level ...
          ! Start at one end of the line
          l => cp%lines(i)%list(j)
!          debug = (i==3.and.j==1.and.s==1)
          if (debug) then
            if (s==1) then
              if (debug) print*, "DEBUG from ",  l%x(1), l%y(1), " to ",l%x(size(l%x)), l%y(size(l%y)), l%uphill
            elseif (s==2) then
              if (debug) print*, "DEBUG from ", l%x(size(l%x)), l%y(size(l%y)), " to ", l%x(1), l%y(1), l%uphill
            endif
          endif
          if (s==1) then
            if (l%eedone) cycle
            allocate(cp%current%x(size(l%x)))
            allocate(cp%current%y(size(l%y)))
            cp%current%x = l%x
            cp%current%y = l%y
            cp%current%uphill = l%uphill
            ce = l%ee
            l%eedone = .true.
            lastse = l%se
          elseif (s==2) then
            if (l%sedone) cycle
            allocate(cp%current%x(size(l%x)))
            allocate(cp%current%y(size(l%y)))
            cp%current%x = l%x(size(l%x):1:-1)
            cp%current%y = l%y(size(l%y):1:-1)
            cp%current%uphill = .not.l%uphill
            ce = l%se
            l%sedone = .true.
            lastse = l%ee
          endif
          if (debug) print*, "DEBUG uphill", cp%current%uphill
          if (debug) then 
          open(file="line."//str(i)//"."//str(j)//"."//str(s)//".dat", unit=8, status="unknown")
            do jjj = 1, size(cp%current%x)
              write(8, *) cp%current%x(jjj), cp%current%y(jjj)
            enddo
          endif
          once = 1
          forward = (s==1)
          minLevel = i
          minLine = j
          lastlinelength = size(l%x)
          do
            if (ce==5) ce = 1
            if (once==5) then
              ! We may have to check the side we are starting on twice.
              call FoX_error("GCONTR INTERNAL ERROR GONE ROUND BOX")
            endif
            x = cp%current%x(size(cp%current%x)) ! Get coordinates at the end of the line we are building
            y = cp%current%y(size(cp%current%y))
            if (minLevel==i) then ! we are on the current level.
              minLine1 = minLine
              minLine2 = 0 
              if (debug) print*, "DEBUG checking from lower level"
            else
              minLine1 = 0 
              minLine2 = minLine
              if (debug) print*, "DEBUG checking from upper level"
            endif
            forward1 = forward
            forward2 = forward
            if (debug) print*, "DEBUG checking lower level", (minLevel>i.and..not.cp%current%uphill)
            minGap1 = nextLineOnEdge(cp%lines(i)%list, cp%current, ce, minLine1, forward1, lastlinelength, lastse, &
              (minLevel>i.and..not.cp%current%uphill))
            ii = i + 1
            minGap2 = -1.0
            if (ii<=size(cp%contours)) then
              if (debug) print*, "DEBUG checking upper level", (minLevel==i.and.cp%current%uphill)
               minGap2 = nextLineOnEdge(cp%lines(ii)%list, cp%current, ce, minLine2, forward2, lastlinelength, lastse, &
                (minLevel==i.and.cp%current%uphill))
            endif
            if (debug) print*, "DEBUG next line (",s,")", minLine1, minLine2, minGap1, minGap2, forward1, forward2, forward
            if (minLine1==0.and.minLine2==0) then
              ! Didn't find anything
              if (debug) print*, "DEBUG round corner"
              minLine = 0
            elseif (minLine2==0 & ! There is no candidate line for the next level up
              .or. (minLine1>0.and.(minGap1<minGap2))) & ! There are two candidate lines, but the one on this level is closer
              then ! two lines are equidistant, stay at current level
              ! next line is on the current level
              if (debug) print*, "DEBUG to lower level"
              minLine = minLine1
              minLevel = i
              forward = forward1
            elseif (minLine1==0 & ! There is no candidate line for this level
              .or. (minLine2>0.and.(minGap2<minGap1))) & ! There are two candidate lines, but the one for the next level is closer
              then ! two lines are equidistant, stay at current level
              ! next line is on the next level
              if (debug) print*, "DEBUG to upper level"
              minLine = minLine2
              minLevel = ii
              forward = forward2
              allOneLevel = .false.
            elseif (minGap2==minGap1) then
              !there are two equidistant candidate lines; but in all such cases we want that from the lower (current) level
              if (debug) print*, "DEBUG same to lower level"
              minLine = minLine1
              minLevel = i
              forward = forward1
            endif
            if (minLine==0) then
              ! Nothing on this side. Add corner to polygon & check next side ...
              select case(ce)
              case(1)
                call add_point(cp%current, 1.0, real(ny))
                if (debug) write(8,*) 1.0, real(ny)
              case(2)
                call add_point(cp%current, real(nx), real(ny))
                if (debug) write(8,*) real(nx), real(ny)
              case(3)
                call add_point(cp%current, real(nx), 1.0)
                if (debug) write(8,*) real(nx), 1.0
              case(4)
                call add_point(cp%current, 1.0, 1.0)
                if (debug) write(8,*) 1.0, 1.0
              end select
              ce = ce + 1
              once = once + 1
            else if (minLevel==i.and.minLine==j) then
              ! We've just found the other end of this line, so this polygon is now complete.
              if (debug) print*, "DEBUG found ok"
              exit
            else
              ! Add the points from the new line we've found to the end of the old line
              l => cp%lines(minLevel)%list(minLine)
              if (forward) then
                ce = l%ee
                if (l%eedone) then
                  call FoX_error("GCONTR INTERNAL ERROR: Adding line twice")
                endif
                if (debug) then
                  print*, "DEBUG now from ", l%x(1), l%y(1), " to ", l%x(size(l%x)), l%y(size(l%x)), l%uphill
                  do jjj = 1, size(l%x)
                    write(8,*) l%x(jjj), l%y(jjj)
                  enddo
                endif
                l%eedone = .true.
                lastse = l%se
              else
                ce = l%se
                if (l%sedone) then
                  call FoX_error("GCONTR INTERNAL ERROR: Adding line twice")
                endif
                if (debug) then
                  print*, "DEBUG now from ", l%x(size(l%x)), l%y(size(l%y)), " to ", l%x(1), l%y(1), l%uphill
                  do jjj = size(l%x), 1, -1
                    write(8,*) l%x(jjj), l%y(jjj)
                  enddo
                endif
                l%sedone = .true.
                lastse = l%ee
              endif
              call concatenate_lines(cp%current, l, forward)
              lastlinelength = size(l%x)
            end if
          enddo
          uphill = cp%current%uphill
          call add_polygon(cp, i)
          cp%polys(i)%list(size(cp%polys(i)%list))%boundary = .true.
          if (.not.allOneLevel) then
            ! This is a mountain, since it starts on level i but also goes up to level i+1
            cp%polys(i)%list(size(cp%polys(i)%list))%mountain = 2
          else
            ! Then we need to remember the starting edge of this polygon to test later:
            if (s==1) then
              cp%polys(i)%list(size(cp%polys(i)%list))%se = cp%lines(i)%list(j)%se
            else
              cp%polys(i)%list(size(cp%polys(i)%list))%se = cp%lines(i)%list(j)%ee
            endif
            if (uphill) then ! the created polygon is a mountain
              cp%polys(i)%list(size(cp%polys(i)%list))%mountain = 2
            else
              cp%polys(i)%list(size(cp%polys(i)%list))%mountain = 1
            endif
          endif
          if (debug) close(8)
        enddo 
      enddo
    enddo

  end subroutine joinUpLines

  subroutine polygonsInPolygons(cp, z)
    type(contourLines), intent(inout) :: cp
    real, dimension(:,:), intent(in) :: z

    integer :: i, j, ii, jj, k, kk
    integer, pointer :: tempList(:)
    real :: x, y, testX, testY, z1, z2, a

    do i = 1, size(cp%contours)
      ! Go over all contour levels
      do j = 1, size(cp%polys(i)%list)
        ! For each level, go over all contours
        allocate(cp%polys(i)%list(j)%contains_lt(0))
        allocate(cp%polys(i)%list(j)%contains_eq(0))
        allocate(cp%polys(i)%list(j)%contains_gt(0))
        allocate(cp%polys(i)%list(j)%borders_lt(0))
        allocate(cp%polys(i)%list(j)%borders_eq(0))
        allocate(cp%polys(i)%list(j)%borders_gt(0))
        if (cp%polys(i)%list(j)%mountain==0) then
          call FoX_error("SHOULDNT BE HERE")
        endif
        a = cp%polys(i)%list(j)%area
        if (a==0.0) cycle
        ! If this polygon has zero area, it's not going to contain anything (by convention)

        ii = i - 1
        ! Find the contour level *below* the current one.
        if (ii>0) then 
          ! Does part of this polygon lie on the boundary?
          do jj = 1, size(cp%polys(ii)%list)
            if (cp%polys(ii)%list(jj)%boundary) cycle ! a polygon lying on the boundary will never be within another.
            if (cp%polys(ii)%list(jj)%area > a) cycle ! A bigger polygon cant be inside a smaller.
            ! Go over all contours at this level - to find out if they are contained
            ! by the previous one.
            ! First, find a point that is definitely not lying on the current polygon
            do kk = 1, size(cp%polys(ii)%list(jj)%x)
              x = cp%polys(ii)%list(jj)%x(kk)
              y = cp%polys(ii)%list(jj)%y(kk)
              if (.not.checkPointOnPolygon(cp%polys(i)%list(j), x, y)) exit
              ! FIXME this is a horribly inefficient loop. Need proper sorting really.
            enddo
            if (kk==size(cp%polys(ii)%list(jj)%x).and.(cp%polys(i)%list(j)%mountain==2)) cycle
            ! If we left the above loop with all points on the polygon, then we
            ! only include this if we are in a valley
            if (pointInPolygon(cp%polys(i)%list(j), x, y)) then
              ! This is definitively a valley, and so is the contained polygon.
              cp%polys(i)%list(j)%mountain = 1
              cp%polys(ii)%list(jj)%mountain = 1
              tempList => cp%polys(i)%list(j)%contains_lt
              allocate(cp%polys(i)%list(j)%contains_lt(size(tempList)+1))
              cp%polys(i)%list(j)%contains_lt(:size(tempList)) = tempList
              cp%polys(i)%list(j)%contains_lt(size(tempList)+1) = jj
              deallocate(tempList)
            endif
          end do
        endif

        ii = i + 1
        if (ii<=size(cp%contours)) then
          do jj = 1, size(cp%polys(ii)%list)
            if (cp%polys(ii)%list(jj)%boundary) cycle ! a polygon lying on the boundary will never be within another.
            if (cp%polys(ii)%list(jj)%area > a) cycle ! A bigger polygon cant be inside a smaller.
            ! Go over all contours at this level - to find out if they are contained
            ! by the previous one.
            ! First, find a point that is definitely not lying on the current polygon
            do kk = 1, size(cp%polys(ii)%list(jj)%x)
              x = cp%polys(ii)%list(jj)%x(kk)
              y = cp%polys(ii)%list(jj)%y(kk)
              if (.not.checkPointOnPolygon(cp%polys(i)%list(j), x, y)) exit
              ! FIXME this is a horribly inefficient loop. Need proper sorting really.
            enddo
            if (kk==size(cp%polys(ii)%list(jj)%x).and.(cp%polys(i)%list(j)%mountain==1)) cycle
            ! If we left the above loop with all points on the polygon, then we
            ! only include this if we are in a mountain
            if (pointInPolygon(cp%polys(i)%list(j), x, y)) then
              cp%polys(i)%list(j)%mountain = 2
              cp%polys(ii)%list(jj)%mountain = 2
              tempList => cp%polys(i)%list(j)%contains_gt
              allocate(cp%polys(i)%list(j)%contains_gt(size(tempList)+1))
              cp%polys(i)%list(j)%contains_gt(:size(tempList)) = tempList
              cp%polys(i)%list(j)%contains_gt(size(tempList)+1) = jj
              deallocate(tempList)
            endif
          end do
        endif

        ! Check all other polygons at the same level - are any contained?
        ii = i
        do jj = 1, size(cp%polys(ii)%list)
          if (cp%polys(ii)%list(jj)%boundary) cycle ! a polygon lying on the boundary will never be within another.
          if (cp%polys(ii)%list(jj)%area > a) cycle ! A bigger polygon cant be inside a smaller.
          if (jj==j) cycle ! A polygon can't be inside itself
          ! Go over all contours at this level - to find out if they are contained
          ! by the previous one.
          ! First, find a point that is definitely not lying on the current polygon
          do kk = 1, size(cp%polys(ii)%list(jj)%x)
            x = cp%polys(ii)%list(jj)%x(kk)
            y = cp%polys(ii)%list(jj)%y(kk)
            if (.not.checkPointOnPolygon(cp%polys(i)%list(j), x, y)) exit
            ! FIXME this is a horribly inefficient loop. Need proper sorting really.
          enddo
          ! It will never be the case that we exit the loop above normally
          if (pointInPolygon(cp%polys(i)%list(j), x, y)) then
            ! propagate mountain/valleyness between two the two polygons, whichever way we can ...
            if (cp%polys(i)%list(j)%mountain==1) then
              cp%polys(ii)%list(jj)%mountain = 2
            elseif (cp%polys(i)%list(j)%mountain==2) then
              cp%polys(ii)%list(jj)%mountain = 1
            elseif (cp%polys(ii)%list(jj)%mountain==1) then
              cp%polys(i)%list(j)%mountain = 2
            elseif (cp%polys(ii)%list(jj)%mountain==2) then
              cp%polys(i)%list(j)%mountain = 1
            endif
            tempList => cp%polys(i)%list(j)%contains_eq
            allocate(cp%polys(i)%list(j)%contains_eq(size(tempList)+1))
            cp%polys(i)%list(j)%contains_eq(:size(tempList)) = tempList
            cp%polys(i)%list(j)%contains_eq(size(tempList)+1) = jj
            deallocate(tempList)
          endif
        end do
      enddo
    enddo
!!$
!!$    ! We also need to know which (boundary-adhering) polygons are next to which other (b-a) polygons.
!!$    ! In principle, we can get this information from the joinUpLines subroutine, but keeping
!!$    ! track of the indices is an arse. Instead, this brute force attack will work, since it relies
!!$    ! on floating point comparisons ... it looks inefficient, but shouldn't be too bad on normal datasets though.
!!$
!!$    do i = 1, size(cp%contours)
!!$      do j = 1, size(cp%polys(i)%list)
!!$        if (cp%polys(i)%list(j)%boundary) then
!!$          l1: do jj = j+1, size(cp%polys(i)%list)
!!$            if (cp%polys(i)%list(jj)%boundary) then
!!$              do k = 1, size(cp%polys(i)%list(j)%x)
!!$                do kk = 1, size(cp%polys(i)%list(jj)%x)
!!$                  if (cp%polys(i)%list(j)%x(k)==cp%polys(i)%list(jj)%x(kk)) then
!!$                    if (cp%polys(i)%list(j)%y(k)==cp%polys(i)%list(jj)%y(kk)) then
!!$                      tempList => cp%polys(i)%list(j)%borders_eq
!!$                      allocate(cp%polys(i)%list(j)%borders_eq(size(tempList)+1))
!!$                      cp%polys(i)%list(j)%borders_eq(:size(tempList)) = tempList
!!$                      cp%polys(i)%list(j)%borders_eq(size(tempList)+1) = jj
!!$                      deallocate(tempList)
!!$                      tempList => cp%polys(i)%list(jj)%borders_eq
!!$                      allocate(cp%polys(i)%list(jj)%borders_eq(size(tempList)+1))
!!$                      cp%polys(i)%list(jj)%borders_eq(:size(tempList)) = tempList
!!$                      cp%polys(i)%list(jj)%borders_eq(size(tempList)+1) = j
!!$                      deallocate(tempList)
!!$                      cycle l1
!!$                    endif
!!$                  endif
!!$                enddo
!!$              enddo
!!$            endif
!!$          enddo l1
!!$          ii = i + 1
!!$          if (ii<=size(cp%contours)) then
!!$            l2: do jj = 1, size(cp%polys(ii)%list)
!!$              if (cp%polys(ii)%list(jj)%boundary) then
!!$                do k = 1, size(cp%polys(i)%list(j)%x)
!!$                  do kk = 1, size(cp%polys(ii)%list(jj)%x)
!!$                    if (cp%polys(i)%list(j)%x(k)==cp%polys(ii)%list(jj)%x(kk)) then
!!$                      if (cp%polys(i)%list(j)%y(k)==cp%polys(ii)%list(jj)%y(kk)) then
!!$                        tempList => cp%polys(i)%list(j)%borders_gt
!!$                        allocate(cp%polys(i)%list(j)%borders_gt(size(tempList)+1))
!!$                        cp%polys(i)%list(j)%borders_gt(:size(tempList)) = tempList
!!$                        cp%polys(i)%list(j)%borders_gt(size(tempList)+1) = jj
!!$                        deallocate(tempList)
!!$                        tempList => cp%polys(ii)%list(jj)%borders_lt
!!$                        allocate(cp%polys(ii)%list(jj)%borders_lt(size(tempList)+1))
!!$                        cp%polys(ii)%list(jj)%borders_lt(:size(tempList)) = tempList
!!$                        cp%polys(ii)%list(jj)%borders_lt(size(tempList)+1) = j
!!$                        deallocate(tempList)
!!$                        cycle l2
!!$                      endif
!!$                    endif
!!$                  enddo
!!$                enddo
!!$              endif
!!$            enddo l2
!!$          endif
!!$        endif
!!$      enddo
!!$    enddo
!!$
!!$    ! Calculating which polygons are *above* their contour (mountains) and which are *below* (valleys)
!!$    ! Rules: if contains_gt has any in it, it is a mountain, and so are the contained polygons,
!!$    !           and any contains_eq are valleys and so are any polygons contained within that.
!!$    !        if contains_lt has any in it, it is a valley, and so are the contained polygons,
!!$    !           and any contains_eq are mountains (and so are any polygons contained within that.
!!$    ! Those rules are captured by the code above. What is not captured ...
!!$    !        if the polygon contains no other, nor is itself contained, then 
!!$    !           it can only be that part of this polygon is a boundary. Test a point on the boundary.
!!$    !   The first point in any such polygon lies on the boundary, as does the last point
!!$    !     (by construction, due to the way we walked the boundaries in JoinUpLines.
!!$    !   Therefore - test the value halfway between the first and last points -- it will lie
!!$    !     on the boundary, and have the characteristics of the inside of the polygon.
!!$    !     If it is greater than the contour, this is a mountain, if less, a valley.
!!$
!!$    do i = 1, size(cp%contours)
!!$      do j = 1, size(cp%polys(i)%list)
!!$        if (cp%polys(i)%list(j)%boundary.and.cp%polys(i)%list(j)%mountain/=0) then
!!$         ! if (cp%polys(i)%list(j)%mountain==1) then
!!$         !   do k = 1, size(cp%polys(i)%list(j)%borders_eq)
!!$         !     jj = cp%polys(i)%list(j)%borders_eq(k)
!!$         !     cp%polys(i)%list(jj)%mountain = 2
!!$         !   enddo
!!$         !   do k = 1, size(cp%polys(i)%list(j)%borders_lt)
!!$         !     jj = cp%polys(i)%list(j)%borders_lt(k)
!!$         !     cp%polys(i-1)%list(jj)%mountain = 1
!!$         !   enddo
!!$         !   do k = 1, size(cp%polys(i)%list(j)%borders_gt)
!!$         !     jj = cp%polys(i)%list(j)%borders_gt(k)
!!$         !     cp%polys(i+1)%list(jj)%mountain = 2
!!$         !   enddo
!!$          if (cp%polys(i)%list(j)%mountain==2) then
!!$            do k = 1, size(cp%polys(i)%list(j)%borders_eq)
!!$              jj = cp%polys(i)%list(j)%borders_eq(k)
!!$              cp%polys(i)%list(jj)%mountain = 1
!!$            enddo
!!$            !do k = 1, size(cp%polys(i)%list(j)%borders_lt)
!!$            !  jj = cp%polys(i)%list(j)%borders_lt(k)
!!$            !  cp%polys(i-1)%list(jj)%mountain = 1
!!$            !enddo
!!$            do k = 1, size(cp%polys(i)%list(j)%borders_gt)
!!$              jj = cp%polys(i)%list(j)%borders_gt(k)
!!$              cp%polys(i+1)%list(jj)%mountain = 2
!!$            enddo
!!$          endif
!!$          do k = 1, size(cp%polys(i)%list(j)%contains_eq)
!!$            jj = cp%polys(i)%list(j)%contains_eq(k)
!!$            if (cp%polys(i)%list(j)%mountain==1) then
!!$              cp%polys(i)%list(jj)%mountain = 2
!!$            else
!!$              cp%polys(i)%list(jj)%mountain = 1
!!$            endif
!!$          enddo
!!$        endif
!!$      enddo
!!$    enddo
!!$
!!$    ! After all that, in 99% of cases we have assigned everywhere. However, it is possible
!!$    ! in regions where only one contour line is present, that we may have some unassigned
!!$    ! regions left over ...
!!$
!!$    do i = 1, size(cp%contours)
!!$      do j = 1, size(cp%polys(i)%list)
!!$        if (cp%polys(i)%list(j)%mountain==0) then
!!$          !the mountain-ness has not yet been set;
!!$          ! Find a grid point between the first & last points - this is guaranteed to be on 
!!$          ! the boundary, and to be a grid point with the correct magnitude ...
!!$          ! This will *fail* if the first and last points are not separated by a grid point - 
!!$          ! but I don't think such polygons can be generated by the algorithm above.
!!$          ! Find the next grid point after 1, in the direction of the last one:
!!$          if (cp%polys(i)%list(j)%x(1) < cp%polys(i)%list(j)%x(size(cp%polys(i)%list(j)%x))) then
!!$            testX = ceiling(cp%polys(i)%list(j)%x(1))
!!$          else
!!$            testX = floor(cp%polys(i)%list(j)%x(1))
!!$          endif
!!$          if (cp%polys(i)%list(j)%y(1) < cp%polys(i)%list(j)%y(size(cp%polys(i)%list(j)%y))) then
!!$            testY = ceiling(cp%polys(i)%list(j)%y(1))
!!$          else
!!$            testY = floor(cp%polys(i)%list(j)%y(1))
!!$          endif
!!$          select case(cp%polys(i)%list(j)%se)
!!$          case(1)
!!$            z1 = z(1, int(testY))
!!$          case(2)
!!$            z1 = z(int(testX), size(z,2))
!!$          case(3)
!!$            z1 = z(size(z,1), int(testY))
!!$          case(4)
!!$            z1 = z(int(testX), 1)
!!$          end select
!!$          if (z1>cp%contours(i)) then
!!$            cp%polys(i)%list(j)%mountain = 1
!!$          else
!!$            cp%polys(i)%list(j)%mountain = 2
!!$          endif
!!$          ! There might be single polygons within here which we need to notify (there can only be 
!!$          ! a single layer - if there were a stack we would have notified them already.
!!$          do k = 1, size(cp%polys(i)%list(j)%contains_eq)
!!$            jj = cp%polys(i)%list(j)%contains_eq(k)
!!$            if (cp%polys(i)%list(j)%mountain==1) then
!!$              cp%polys(i)%list(jj)%mountain = 2
!!$            else
!!$              cp%polys(i)%list(jj)%mountain = 1
!!$            endif
!!$          enddo
!!$        endif
!!$      enddo
!!$    enddo

  end subroutine polygonsInPolygons
  
  subroutine addLineToLineList(l, ls)
    type(richline), intent(in) :: l
    type(lineList), intent(inout) :: ls

    integer :: i
    type(line), pointer :: t(:)

    t => ls%list
    allocate(ls%list(size(t)+1))
    do i = 1, size(t)
      ls%list(i)%x => t(i)%x
      ls%list(i)%y => t(i)%y
    enddo
    deallocate(t)
    allocate(ls%list(i)%x(size(l%x)))
    allocate(ls%list(i)%y(size(l%y)))
    ls%list(i)%x = l%x
    ls%list(i)%y = l%y

  end subroutine addLineToLineList

  subroutine addPolygonToPolygonList(cp, c, n, ls)
    type(contourLines), intent(in) :: cp
    integer, intent(in) :: c, n
    type(polygonList), intent(inout) :: ls

    integer :: i, j, j2, k
    type(polygon), pointer :: t(:)

    t => ls%list
    allocate(ls%list(size(t)+1))
    do i = 1, size(t)
      ls%list(i)%outerBoundary%x => t(i)%outerBoundary%x
      ls%list(i)%outerBoundary%y => t(i)%outerBoundary%y
      ls%list(i)%innerBoundaries%list => t(i)%innerBoundaries%list
    enddo
    deallocate(t)
    ! Now add new one:
    allocate(ls%list(i)%outerBoundary%x(size(cp%polys(c)%list(n)%x)))
    allocate(ls%list(i)%outerBoundary%y(size(cp%polys(c)%list(n)%y)))
    ls%list(i)%outerBoundary%x = cp%polys(c)%list(n)%x
    ls%list(i)%outerBoundary%y = cp%polys(c)%list(n)%y
    ! And add all the relevant innerBoundaries:
    allocate(ls%list(i)%innerBoundaries%list( &
      size(cp%polys(c)%list(n)%contains_lt)&
      +size(cp%polys(c)%list(n)%contains_eq) &
      +size(cp%polys(c)%list(n)%contains_gt)))
    do j = 1, size(cp%polys(c)%list(n)%contains_lt)
      k = cp%polys(c)%list(n)%contains_lt(j)
      allocate(ls%list(i)%innerBoundaries%list(j)%x(size(cp%polys(c-1)%list(k)%x)))
      allocate(ls%list(i)%innerBoundaries%list(j)%y(size(cp%polys(c-1)%list(k)%y)))
      ls%list(i)%innerBoundaries%list(j)%x = cp%polys(c-1)%list(k)%x
      ls%list(i)%innerBoundaries%list(j)%y = cp%polys(c-1)%list(k)%y
    enddo
    do j = 1, size(cp%polys(c)%list(n)%contains_eq)
      k = cp%polys(c)%list(n)%contains_eq(j)
      j2 = size(cp%polys(c)%list(n)%contains_lt) + j
      allocate(ls%list(i)%innerBoundaries%list(j2)%x(size(cp%polys(c)%list(k)%x)))
      allocate(ls%list(i)%innerBoundaries%list(j2)%y(size(cp%polys(c)%list(k)%y)))
      ls%list(i)%innerBoundaries%list(j2)%x = cp%polys(c)%list(k)%x
      ls%list(i)%innerBoundaries%list(j2)%y = cp%polys(c)%list(k)%y
    enddo
    do j = 1, size(cp%polys(c)%list(n)%contains_gt)
      k = cp%polys(c)%list(n)%contains_gt(j)
      j2 = size(cp%polys(c)%list(n)%contains_lt) + size(cp%polys(c)%list(n)%contains_eq) + j
      allocate(ls%list(i)%innerBoundaries%list(j2)%x(size(cp%polys(c+1)%list(k)%x)))
      allocate(ls%list(i)%innerBoundaries%list(j2)%y(size(cp%polys(c+1)%list(k)%y)))
      ls%list(i)%innerBoundaries%list(j2)%x = cp%polys(c+1)%list(k)%x
      ls%list(i)%innerBoundaries%list(j2)%y = cp%polys(c+1)%list(k)%y
    enddo
    
  end subroutine addPolygonToPolygonList

  subroutine addLinesToContourObject(cp, o)
    type(contourLines), intent(in) :: cp
    type(contourObject), intent(inout) :: o
    
    integer :: i, j, k


    if (debug) open(file="alllines.dat", unit=9, status='replace', action='write')
    do i = lbound(cp%contours, 1), ubound(cp%contours, 1)
      if (debug) write(9, *) "LEVEL", i
        do j = 1, size(cp%lines(i)%list)
          if (debug) then
            write(9, *) "LINE", size(cp%lines(i)%list(j)%x)
            do k = 1, size(cp%lines(i)%list(j)%x)
              write(9, *) cp%lines(i)%list(j)%x(k), cp%lines(i)%list(j)%y(k)
            enddo
          endif
        call addLineToLineList(cp%lines(i)%list(j), o%lines(i))
      enddo
      do j = 1, size(cp%polys(i)%list)
        if (debug) then
          write(9, *) "LINE", size(cp%polys(i)%list(j)%x)+1
          do k = 1, size(cp%polys(i)%list(j)%x)
            write(9, *) cp%polys(i)%list(j)%x(k), cp%polys(i)%list(j)%y(k)
          enddo
          write(9, *) cp%polys(i)%list(j)%x(1), cp%polys(i)%list(j)%y(1)
        endif
        call addLineToLineList(cp%polys(i)%list(j), o%lines(i))
      enddo
    enddo
    if (debug) close(9)

  end subroutine addLinesToContourObject

  subroutine addPolygonsToContourObject(cp, o)
    type(contourLines), intent(in) :: cp
    type(contourObject), intent(inout) :: o

    integer :: i, j, ii

    do i = 1, size(cp%contours)
      do j = 1, size(cp%polys(i)%list)
        if (cp%polys(i)%list(j)%mountain==1) then
          ii = i - 1
        elseif (cp%polys(i)%list(j)%mountain==2) then
          ii = i
        else
          call FoX_error('GCONTR INTERNAL ERROR: POLYGON CALCULATING ERROR!!! at ')
        endif
        ! print*, "AREA", i, j, cp%polys(i)%list(j)%area
        if (cp%polys(i)%list(j)%area==0.0) cycle
        call addPolygonToPolygonList(cp, i, j, o%polys(ii))
      enddo
    enddo
  end subroutine addPolygonsToContourObject
#endif
end module m_contours
