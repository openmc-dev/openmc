module mesh_header

  use, intrinsic :: ISO_C_BINDING

  use hdf5

  use constants
  use dict_header, only: DictIntInt
  use error
  use hdf5_interface
  use string, only: to_str, to_lower
  use xml_interface

  implicit none
  private
  public :: free_memory_mesh
  public :: openmc_extend_meshes
  public :: openmc_get_mesh_index
  public :: openmc_mesh_get_id
  public :: openmc_mesh_get_dimension
  public :: openmc_mesh_get_params
  public :: openmc_mesh_set_id
  public :: openmc_mesh_set_dimension
  public :: openmc_mesh_set_params

!===============================================================================
! STRUCTUREDMESH represents a tessellation of n-dimensional Euclidean space by
! congruent squares or cubes
!===============================================================================

  type, public :: RegularMesh
    integer :: id = -1                     ! user-specified id
    integer :: type = MESH_REGULAR         ! rectangular, hexagonal
    integer(C_INT) :: n_dimension          ! rank of mesh
    real(8) :: volume_frac                 ! volume fraction of each cell
    integer(C_INT), allocatable :: dimension(:)   ! number of cells in each direction
    real(C_DOUBLE), allocatable :: lower_left(:)  ! lower-left corner of mesh
    real(C_DOUBLE), allocatable :: upper_right(:) ! upper-right corner of mesh
    real(C_DOUBLE), allocatable :: width(:)       ! width of each mesh cell
  contains
    procedure :: from_xml => regular_from_xml
    procedure :: get_bin => regular_get_bin
    procedure :: get_indices => regular_get_indices
    procedure :: get_bin_from_indices => regular_get_bin_from_indices
    procedure :: get_indices_from_bin => regular_get_indices_from_bin
    procedure :: intersects => regular_intersects
    procedure :: to_hdf5 => regular_to_hdf5
  end type RegularMesh

  integer(C_INT32_T), public, bind(C) :: n_meshes = 0 ! # of structured meshes

  type(RegularMesh), public, allocatable, target :: meshes(:)

  ! Dictionary that maps user IDs to indices in 'meshes'
  type(DictIntInt), public :: mesh_dict

contains

  subroutine regular_from_xml(this, node)
    class(RegularMesh), intent(inout) :: this
    type(XMLNode), intent(in) :: node

    integer :: n
    character(MAX_LINE_LEN) :: temp_str

    ! Copy mesh id
    if (check_for_node(node, "id")) then
      call get_node_value(node, "id", this % id)

      ! Check to make sure 'id' hasn't been used
      if (mesh_dict % has(this % id)) then
        call fatal_error("Two or more meshes use the same unique ID: " &
             // to_str(this % id))
      end if
    end if

    ! Read mesh type
    if (check_for_node(node, "type")) then
      call get_node_value(node, "type", temp_str)
      select case (to_lower(temp_str))
      case ('rect', 'rectangle', 'rectangular')
        call warning("Mesh type '" // trim(temp_str) // "' is deprecated. &
             &Please use 'regular' instead.")
        this % type = MESH_REGULAR
      case ('regular')
        this % type = MESH_REGULAR
      case default
        call fatal_error("Invalid mesh type: " // trim(temp_str))
      end select
    else
      this % type = MESH_REGULAR
    end if

    ! Determine number of dimensions for mesh
    if (check_for_node(node, "dimension")) then
      n = node_word_count(node, "dimension")
      if (n /= 1 .and. n /= 2 .and. n /= 3) then
        call fatal_error("Mesh must be one, two, or three dimensions.")
      end if
      this % n_dimension = n

      ! Allocate attribute arrays
      allocate(this % dimension(n))

      ! Check that dimensions are all greater than zero
      call get_node_array(node, "dimension", this % dimension)
      if (any(this % dimension <= 0)) then
        call fatal_error("All entries on the <dimension> element for a tally &
             &mesh must be positive.")
      end if
    end if

    ! Check for lower-left coordinates
    if (check_for_node(node, "lower_left")) then
      n = node_word_count(node, "lower_left")
      allocate(this % lower_left(n))

      ! Read mesh lower-left corner location
      call get_node_array(node, "lower_left", this % lower_left)
    else
      call fatal_error("Must specify <lower_left> on a mesh.")
    end if

    if (check_for_node(node, "width")) then
      ! Make sure both upper-right or width were specified
      if (check_for_node(node, "upper_right")) then
        call fatal_error("Cannot specify both <upper_right> and <width> on a &
             &mesh.")
      end if

      n = node_word_count(node, "width")
      allocate(this % width(n))
      allocate(this % upper_right(n))

      ! Check to ensure width has same dimensions
      if (n /= size(this % lower_left)) then
        call fatal_error("Number of entries on <width> must be the same as &
             &the number of entries on <lower_left>.")
      end if

      ! Check for negative widths
      call get_node_array(node, "width", this % width)
      if (any(this % width < ZERO)) then
        call fatal_error("Cannot have a negative <width> on a tally mesh.")
      end if

      ! Set width and upper right coordinate
      this % upper_right = this % lower_left + this % dimension * this % width

    elseif (check_for_node(node, "upper_right")) then
      n = node_word_count(node, "upper_right")
      allocate(this % upper_right(n))
      allocate(this % width(n))

      ! Check to ensure width has same dimensions
      if (n /= size(this % lower_left)) then
        call fatal_error("Number of entries on <upper_right> must be the &
             &same as the number of entries on <lower_left>.")
      end if

      ! Check that upper-right is above lower-left
      call get_node_array(node, "upper_right", this % upper_right)
      if (any(this % upper_right < this % lower_left)) then
        call fatal_error("The <upper_right> coordinates must be greater than &
             &the <lower_left> coordinates on a tally mesh.")
      end if

      ! Set width and upper right coordinate
      this % width = (this % upper_right - this % lower_left) / this % dimension
    else
      call fatal_error("Must specify either <upper_right> and <width> on a &
           &mesh.")
    end if

    if (allocated(this % dimension)) then
      if (size(this % dimension) /= size(this % lower_left)) then
        call fatal_error("Number of entries on <lower_left> must be the same &
             &as the number of entries on <dimension>.")
      end if

      ! Set volume fraction
      this % volume_frac = ONE/real(product(this % dimension),8)
    end if

  end subroutine regular_from_xml

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
    real(8), intent(in) :: xyz0(:)
    real(8), intent(in) :: xyz1(:)
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

!===============================================================================
! TO_HDF5 writes the mesh data to an HDF5 group
!===============================================================================

  subroutine regular_to_hdf5(this, group)
    class(RegularMesh), intent(in) :: this
    integer(HID_T), intent(in) :: group

    integer(HID_T) :: mesh_group

    mesh_group = create_group(group, "mesh " // trim(to_str(this % id)))

    call write_dataset(mesh_group, "type", "regular")
    call write_dataset(mesh_group, "dimension", this % dimension)
    call write_dataset(mesh_group, "lower_left", this % lower_left)
    call write_dataset(mesh_group, "upper_right", this % upper_right)
    call write_dataset(mesh_group, "width", this % width)

    call close_group(mesh_group)
  end subroutine regular_to_hdf5

!===============================================================================
! FREE_MEMORY_MESH deallocates global arrays defined in this module
!===============================================================================

  subroutine free_memory_mesh()
    n_meshes = 0
    if (allocated(meshes)) deallocate(meshes)
    call mesh_dict % clear()
  end subroutine free_memory_mesh

!===============================================================================
!                               C API FUNCTIONS
!===============================================================================

  function openmc_extend_meshes(n, index_start, index_end) result(err) bind(C)
    ! Extend the meshes array by n elements
    integer(C_INT32_T), value, intent(in) :: n
    integer(C_INT32_T), optional, intent(out) :: index_start
    integer(C_INT32_T), optional, intent(out) :: index_end
    integer(C_INT) :: err

    type(RegularMesh), allocatable :: temp(:) ! temporary meshes array

    if (n_meshes == 0) then
      ! Allocate meshes array
      allocate(meshes(n))
    else
      ! Allocate meshes array with increased size
      allocate(temp(n_meshes + n))

      ! Copy original meshes to temporary array
      temp(1:n_meshes) = meshes

      ! Move allocation from temporary array
      call move_alloc(FROM=temp, TO=meshes)
    end if

    ! Return indices in meshes array
    if (present(index_start)) index_start = n_meshes + 1
    if (present(index_end)) index_end = n_meshes + n
    n_meshes = n_meshes + n

    err = 0
  end function openmc_extend_meshes


  function openmc_get_mesh_index(id, index) result(err) bind(C)
    ! Return the index in the meshes array of a mesh with a given ID
    integer(C_INT32_T), value :: id
    integer(C_INT32_T), intent(out) :: index
    integer(C_INT) :: err

    if (allocated(meshes)) then
      if (mesh_dict % has(id)) then
        index = mesh_dict % get(id)
        err = 0
      else
        err = E_INVALID_ID
        call set_errmsg("No mesh exists with ID=" // trim(to_str(id)) // ".")
      end if
    else
      err = E_ALLOCATE
      call set_errmsg("Memory has not been allocated for meshes.")
    end if
  end function openmc_get_mesh_index


  function openmc_mesh_get_id(index, id) result(err) bind(C)
    ! Return the ID of a mesh
    integer(C_INT32_T), value       :: index
    integer(C_INT32_T), intent(out) :: id
    integer(C_INT) :: err

    if (index >= 1 .and. index <= size(meshes)) then
      id = meshes(index) % id
      err = 0
    else
      err = E_OUT_OF_BOUNDS
      call set_errmsg("Index in meshes array is out of bounds.")
    end if
  end function openmc_mesh_get_id


  function openmc_mesh_set_id(index, id) result(err) bind(C)
    ! Set the ID of a mesh
    integer(C_INT32_T), value, intent(in) :: index
    integer(C_INT32_T), value, intent(in) :: id
    integer(C_INT) :: err

    if (index >= 1 .and. index <= n_meshes) then
      meshes(index) % id = id
      call mesh_dict % set(id, index)
      err = 0
    else
      err = E_OUT_OF_BOUNDS
      call set_errmsg("Index in meshes array is out of bounds.")
    end if
  end function openmc_mesh_set_id


  function openmc_mesh_get_dimension(index, dims, n) result(err) bind(C)
    ! Get the dimension of a mesh
    integer(C_INT32_T), value, intent(in) :: index
    type(C_PTR),               intent(out) :: dims
    integer(C_INT),            intent(out) :: n
    integer(C_INT) :: err

    if (index >= 1 .and. index <= n_meshes) then
      dims = C_LOC(meshes(index) % dimension)
      n = meshes(index) % n_dimension
      err = 0
    else
      err = E_OUT_OF_BOUNDS
      call set_errmsg("Index in meshes array is out of bounds.")
    end if
  end function openmc_mesh_get_dimension


  function openmc_mesh_set_dimension(index, n, dims) result(err) bind(C)
    ! Set the dimension of a mesh
    integer(C_INT32_T), value, intent(in) :: index
    integer(C_INT),     value, intent(in) :: n
    integer(C_INT),            intent(in) :: dims(n)
    integer(C_INT) :: err

    if (index >= 1 .and. index <= n_meshes) then
      associate (m => meshes(index))
        if (allocated(m % dimension)) deallocate (m % dimension)
        if (allocated(m % lower_left)) deallocate (m % lower_left)
        if (allocated(m % upper_right)) deallocate (m % upper_right)
        if (allocated(m % width)) deallocate (m % width)

        m % n_dimension = n
        allocate(m % dimension(n))
        allocate(m % lower_left(n))
        allocate(m % upper_right(n))
        allocate(m % width(n))

        ! Copy dimension
        m % dimension(:) = dims
      end associate
      err = 0
    else
      err = E_OUT_OF_BOUNDS
      call set_errmsg("Index in meshes array is out of bounds.")
    end if
  end function openmc_mesh_set_dimension


  function openmc_mesh_get_params(index, ll, ur, width, n) result(err) bind(C)
    ! Get the mesh parameters
    integer(C_INT32_T), value, intent(in) :: index
    type(C_PTR), intent(out) :: ll
    type(C_PTR), intent(out) :: ur
    type(C_PTR), intent(out) :: width
    integer(C_INT), intent(out) :: n
    integer(C_INT) :: err

    err = 0
    if (index >= 1 .and. index <= n_meshes) then
      associate (m => meshes(index))
        if (allocated(m % lower_left)) then
          ll = C_LOC(m % lower_left(1))
          ur = C_LOC(m % upper_right(1))
          width = C_LOC(m % width(1))
          n = m % n_dimension
        else
          err = E_ALLOCATE
          call set_errmsg("Mesh parameters have not been set.")
        end if
      end associate
    else
      err = E_OUT_OF_BOUNDS
      call set_errmsg("Index in meshes array is out of bounds.")
    end if
  end function openmc_mesh_get_params


  function openmc_mesh_set_params(index, n, ll, ur, width) result(err) bind(C)
    ! Set the mesh parameters
    integer(C_INT32_T), value, intent(in) :: index
    integer(C_INT),     value, intent(in) :: n
    real(C_DOUBLE), intent(in), optional :: ll(n)
    real(C_DOUBLE), intent(in), optional :: ur(n)
    real(C_DOUBLE), intent(in), optional :: width(n)
    integer(C_INT) :: err

    err = 0
    if (index >= 1 .and. index <= n_meshes) then
      associate (m => meshes(index))
        if (allocated(m % lower_left)) deallocate (m % lower_left)
        if (allocated(m % upper_right)) deallocate (m % upper_right)
        if (allocated(m % width)) deallocate (m % width)

        allocate(m % lower_left(n))
        allocate(m % upper_right(n))
        allocate(m % width(n))

        if (present(ll) .and. present(ur)) then
          m % lower_left(:) = ll
          m % upper_right(:) = ur
          m % width(:) = (ur - ll) / m % dimension
        elseif (present(ll) .and. present(width)) then
          m % lower_left(:) = ll
          m % width(:) = width
          m % upper_right(:) = ll + width * m % dimension
        elseif (present(ur) .and. present(width)) then
          m % upper_right(:) = ur
          m % width(:) = width
          m % lower_left(:) = ur - width * m % dimension
        else
          err = E_INVALID_ARGUMENT
          call set_errmsg("At least two parameters must be specified.")
        end if
      end associate
    else
      err = E_OUT_OF_BOUNDS
      call set_errmsg("Index in meshes array is out of bounds.")
    end if
  end function openmc_mesh_set_params

end module mesh_header
