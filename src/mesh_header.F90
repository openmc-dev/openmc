module mesh_header

  use constants, only: NO_BIN_FOUND

  implicit none

!===============================================================================
! STRUCTUREDMESH represents a tessellation of n-dimensional Euclidean space by
! congruent squares or cubes
!===============================================================================

  type RegularMesh
    integer :: id                          ! user-specified id
    integer :: type                        ! rectangular, hexagonal
    integer :: n_dimension                 ! rank of mesh
    real(8) :: volume_frac                 ! volume fraction of each cell
    integer, allocatable :: dimension(:)   ! number of cells in each direction
    real(8), allocatable :: lower_left(:)  ! lower-left corner of mesh
    real(8), allocatable :: upper_right(:) ! upper-right corner of mesh
    real(8), allocatable :: width(:)       ! width of each mesh cell
  contains
    procedure :: get_bin => regular_get_bin
    procedure :: get_indices => regular_get_indices
    procedure :: get_bin_from_indices => regular_get_bin_from_indices
    procedure :: get_indices_from_bin => regular_get_indices_from_bin
    procedure :: intersects => regular_intersects
  end type RegularMesh

contains

!===============================================================================
! GET_MESH_BIN determines the tally bin for a particle in a structured mesh
!===============================================================================

  pure subroutine regular_get_bin(this, xyz, bin)
    class(RegularMesh), intent(in) :: this
    real(8), intent(in)           :: xyz(:) ! coordinates
    integer, intent(out)          :: bin    ! tally bin

    integer :: n       ! size of mesh
    integer :: d       ! mesh dimension index
    integer :: ijk(3)  ! indices in mesh
    logical :: in_mesh ! was given coordinate in mesh at all?

    ! Get number of dimensions
    n = this % n_dimension

    ! Loop over the dimensions of the mesh
    do d = 1, n

      ! Check for cases where particle is outside of mesh
      if (xyz(d) < this % lower_left(d)) then
        bin = NO_BIN_FOUND
        return
      elseif (xyz(d) > this % upper_right(d)) then
        bin = NO_BIN_FOUND
        return
      end if
    end do

    ! Determine indices
    call this % get_indices(xyz, ijk, in_mesh)

    ! Convert indices to bin
    if (in_mesh) then
      bin = this % get_bin_from_indices(ijk)
    else
      bin = NO_BIN_FOUND
    end if

  end subroutine regular_get_bin

!===============================================================================
! GET_MESH_INDICES determines the indices of a particle in a structured mesh
!===============================================================================

  pure subroutine regular_get_indices(this, xyz, ijk, in_mesh)
    class(RegularMesh), intent(in) :: this
    real(8), intent(in)           :: xyz(:)  ! coordinates to check
    integer, intent(out)          :: ijk(:)  ! indices in mesh
    logical, intent(out)          :: in_mesh ! were given coords in mesh?

    ! Find particle in mesh
    ijk(:this % n_dimension) = ceiling((xyz(:this % n_dimension) - &
         this % lower_left)/this % width)

    ! Determine if particle is in mesh
    if (any(ijk(:this % n_dimension) < 1) .or. &
         any(ijk(:this % n_dimension) > this % dimension)) then
      in_mesh = .false.
    else
      in_mesh = .true.
    end if

  end subroutine regular_get_indices

!===============================================================================
! MESH_INDICES_TO_BIN maps (i), (i,j), or (i,j,k) indices to a single bin number
! for use in a TallyObject results array
!===============================================================================

  pure function regular_get_bin_from_indices(this, ijk) result(bin)
    class(RegularMesh), intent(in) :: this
    integer, intent(in)           :: ijk(:)
    integer                       :: bin

    if (this % n_dimension == 1) then
      bin = ijk(1)
    elseif (this % n_dimension == 2) then
      bin = (ijk(2) - 1) * this % dimension(1) + ijk(1)
    elseif (this % n_dimension == 3) then
      bin = ((ijk(3) - 1) * this % dimension(2) + (ijk(2) - 1)) &
           * this % dimension(1) + ijk(1)
    end if

  end function regular_get_bin_from_indices

!===============================================================================
! BIN_TO_MESH_INDICES maps a single mesh bin from a TallyObject results array to
! (i), (i,j), or (i,j,k) indices
!===============================================================================

  pure subroutine regular_get_indices_from_bin(this, bin, ijk)
    class(RegularMesh), intent(in) :: this
    integer, intent(in)           :: bin
    integer, intent(out)          :: ijk(:)

    if (this % n_dimension == 1) then
      ijk(1) = bin
    else if (this % n_dimension == 2) then
      ijk(1) = mod(bin - 1, this % dimension(1)) + 1
      ijk(2) = (bin - 1)/this % dimension(1) + 1
    else if (this % n_dimension == 3) then
      ijk(1) = mod(bin - 1, this % dimension(1)) + 1
      ijk(2) = mod(bin - 1, this % dimension(1) * this % dimension(2)) &
           / this % dimension(1) + 1
      ijk(3) = (bin - 1)/(this % dimension(1) * this % dimension(2)) + 1
    end if

  end subroutine regular_get_indices_from_bin

!===============================================================================
! MESH_INTERSECTS determines if a line between xyz0 and xyz1 intersects the
! outer boundary of the given mesh. This is important for determining whether a
! track will score to a mesh tally.
!===============================================================================

  pure function regular_intersects(this, xyz0, xyz1) result(intersects)
    class(RegularMesh), intent(in) :: this
    real(8), intent(in) :: xyz0(1)
    real(8), intent(in) :: xyz1(1)
    logical :: intersects

    select case(this % n_dimension)
    case (1)
      intersects = mesh_intersects_1d(this, xyz0, xyz1)
    case (2)
      intersects = mesh_intersects_2d(this, xyz0, xyz1)
    case (3)
      intersects = mesh_intersects_3d(this, xyz0, xyz1)
    end select
  end function regular_intersects

  pure function mesh_intersects_1d(m, xyz0, xyz1) result(intersects)
    type(RegularMesh), intent(in) :: m
    real(8), intent(in) :: xyz0(:)
    real(8), intent(in) :: xyz1(:)
    logical :: intersects

    real(8) :: x0   ! track start point
    real(8) :: x1   ! track end point
    real(8) :: xm0 ! lower-left coordinates of mesh
    real(8) :: xm1 ! upper-right coordinates of mesh

    ! Copy coordinates of starting point
    x0 = xyz0(1)

    ! Copy coordinates of ending point
    x1 = xyz1(1)

    ! Copy coordinates of mesh lower_left
    xm0 = m % lower_left(1)

    ! Copy coordinates of mesh upper_right
    xm1 = m % upper_right(1)

    ! Set default value for intersects
    intersects = .false.

    ! Check if line intersects left surface
    if ((x0 < xm0 .and. x1 > xm0) .or. (x0 > xm0 .and. x1 < xm0)) then
      intersects = .true.
      return
    end if

    ! Check if line intersects right surface
    if ((x0 < xm1 .and. x1 > xm1) .or. (x0 > xm1 .and. x1 < xm1)) then
      intersects = .true.
      return
    end if

  end function mesh_intersects_1d

  pure function mesh_intersects_2d(m, xyz0, xyz1) result(intersects)
    type(RegularMesh), intent(in) :: m
    real(8), intent(in) :: xyz0(:)
    real(8), intent(in) :: xyz1(:)
    logical :: intersects

    real(8) :: x0, y0   ! track start point
    real(8) :: x1, y1   ! track end point
    real(8) :: xi, yi   ! track intersection point with mesh
    real(8) :: xm0, ym0 ! lower-left coordinates of mesh
    real(8) :: xm1, ym1 ! upper-right coordinates of mesh

    ! Copy coordinates of starting point
    x0 = xyz0(1)
    y0 = xyz0(2)

    ! Copy coordinates of ending point
    x1 = xyz1(1)
    y1 = xyz1(2)

    ! Copy coordinates of mesh lower_left
    xm0 = m % lower_left(1)
    ym0 = m % lower_left(2)

    ! Copy coordinates of mesh upper_right
    xm1 = m % upper_right(1)
    ym1 = m % upper_right(2)

    ! Set default value for intersects
    intersects = .false.

    ! Check if line intersects left surface -- calculate the intersection point
    ! y
    if ((x0 < xm0 .and. x1 > xm0) .or. (x0 > xm0 .and. x1 < xm0)) then
      yi = y0 + (xm0 - x0) * (y1 - y0) / (x1 - x0)
      if (yi >= ym0 .and. yi < ym1) then
        intersects = .true.
        return
      end if
    end if

    ! Check if line intersects back surface -- calculate the intersection point
    ! x
    if ((y0 < ym0 .and. y1 > ym0) .or. (y0 > ym0 .and. y1 < ym0)) then
      xi = x0 + (ym0 - y0) * (x1 - x0) / (y1 - y0)
      if (xi >= xm0 .and. xi < xm1) then
        intersects = .true.
        return
      end if
    end if

    ! Check if line intersects right surface -- calculate the intersection
    ! point y
    if ((x0 < xm1 .and. x1 > xm1) .or. (x0 > xm1 .and. x1 < xm1)) then
      yi = y0 + (xm1 - x0) * (y1 - y0) / (x1 - x0)
      if (yi >= ym0 .and. yi < ym1) then
        intersects = .true.
        return
      end if
    end if

    ! Check if line intersects front surface -- calculate the intersection point
    ! x
    if ((y0 < ym1 .and. y1 > ym1) .or. (y0 > ym1 .and. y1 < ym1)) then
      xi = x0 + (ym1 - y0) * (x1 - x0) / (y1 - y0)
      if (xi >= xm0 .and. xi < xm1) then
        intersects = .true.
        return
      end if
    end if

  end function mesh_intersects_2d

  pure function mesh_intersects_3d(m, xyz0, xyz1) result(intersects)
    type(RegularMesh), intent(in) :: m
    real(8), intent(in) :: xyz0(:)
    real(8), intent(in) :: xyz1(:)
    logical :: intersects

    real(8) :: x0, y0, z0    ! track start point
    real(8) :: x1, y1, z1    ! track end point
    real(8) :: xi, yi, zi    ! track intersection point with mesh
    real(8) :: xm0, ym0, zm0 ! lower-left coordinates of mesh
    real(8) :: xm1, ym1, zm1 ! upper-right coordinates of mesh

    ! Copy coordinates of starting point
    x0 = xyz0(1)
    y0 = xyz0(2)
    z0 = xyz0(3)

    ! Copy coordinates of ending point
    x1 = xyz1(1)
    y1 = xyz1(2)
    z1 = xyz1(3)

    ! Copy coordinates of mesh lower_left
    xm0 = m % lower_left(1)
    ym0 = m % lower_left(2)
    zm0 = m % lower_left(3)

    ! Copy coordinates of mesh upper_right
    xm1 = m % upper_right(1)
    ym1 = m % upper_right(2)
    zm1 = m % upper_right(3)

    ! Set default value for intersects
    intersects = .false.

    ! Check if line intersects left surface -- calculate the intersection point
    ! (y,z)
    if ((x0 < xm0 .and. x1 > xm0) .or. (x0 > xm0 .and. x1 < xm0)) then
      yi = y0 + (xm0 - x0) * (y1 - y0) / (x1 - x0)
      zi = z0 + (xm0 - x0) * (z1 - z0) / (x1 - x0)
      if (yi >= ym0 .and. yi < ym1 .and. zi >= zm0 .and. zi < zm1) then
        intersects = .true.
        return
      end if
    end if

    ! Check if line intersects back surface -- calculate the intersection point
    ! (x,z)
    if ((y0 < ym0 .and. y1 > ym0) .or. (y0 > ym0 .and. y1 < ym0)) then
      xi = x0 + (ym0 - y0) * (x1 - x0) / (y1 - y0)
      zi = z0 + (ym0 - y0) * (z1 - z0) / (y1 - y0)
      if (xi >= xm0 .and. xi < xm1 .and. zi >= zm0 .and. zi < zm1) then
        intersects = .true.
        return
      end if
    end if

    ! Check if line intersects bottom surface -- calculate the intersection
    ! point (x,y)
    if ((z0 < zm0 .and. z1 > zm0) .or. (z0 > zm0 .and. z1 < zm0)) then
      xi = x0 + (zm0 - z0) * (x1 - x0) / (z1 - z0)
      yi = y0 + (zm0 - z0) * (y1 - y0) / (z1 - z0)
      if (xi >= xm0 .and. xi < xm1 .and. yi >= ym0 .and. yi < ym1) then
        intersects = .true.
        return
      end if
    end if

    ! Check if line intersects right surface -- calculate the intersection point
    ! (y,z)
    if ((x0 < xm1 .and. x1 > xm1) .or. (x0 > xm1 .and. x1 < xm1)) then
      yi = y0 + (xm1 - x0) * (y1 - y0) / (x1 - x0)
      zi = z0 + (xm1 - x0) * (z1 - z0) / (x1 - x0)
      if (yi >= ym0 .and. yi < ym1 .and. zi >= zm0 .and. zi < zm1) then
        intersects = .true.
        return
      end if
    end if

    ! Check if line intersects front surface -- calculate the intersection point
    ! (x,z)
    if ((y0 < ym1 .and. y1 > ym1) .or. (y0 > ym1 .and. y1 < ym1)) then
      xi = x0 + (ym1 - y0) * (x1 - x0) / (y1 - y0)
      zi = z0 + (ym1 - y0) * (z1 - z0) / (y1 - y0)
      if (xi >= xm0 .and. xi < xm1 .and. zi >= zm0 .and. zi < zm1) then
        intersects = .true.
        return
      end if
    end if

    ! Check if line intersects top surface -- calculate the intersection point
    ! (x,y)
    if ((z0 < zm1 .and. z1 > zm1) .or. (z0 > zm1 .and. z1 < zm1)) then
      xi = x0 + (zm1 - z0) * (x1 - x0) / (z1 - z0)
      yi = y0 + (zm1 - z0) * (y1 - y0) / (z1 - z0)
      if (xi >= xm0 .and. xi < xm1 .and. yi >= ym0 .and. yi < ym1) then
        intersects = .true.
        return
      end if
    end if

  end function mesh_intersects_3d

end module mesh_header
