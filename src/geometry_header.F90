module geometry_header

  use, intrinsic :: ISO_C_BINDING

  use algorithm,       only: find
  use constants,       only: HALF, TWO, THREE, INFINITY, K_BOLTZMANN, &
                             MATERIAL_VOID, NONE
  use dict_header,     only: DictCharInt, DictIntInt
  use material_header, only: Material, materials, material_dict, n_materials
  use nuclide_header
  use sab_header
  use stl_vector,      only: VectorReal
  use string,          only: to_lower

  implicit none

!===============================================================================
! UNIVERSE defines a geometry that fills all phase space
!===============================================================================

  type Universe
    integer :: id                     ! Unique ID
    integer :: type                   ! Type
    integer, allocatable :: cells(:)  ! List of cells within
    real(8) :: x0                     ! Translation in x-coordinate
    real(8) :: y0                     ! Translation in y-coordinate
    real(8) :: z0                     ! Translation in z-coordinate
  end type Universe

!===============================================================================
! LATTICE abstract type for ordered array of universes.
!===============================================================================

  type, abstract :: Lattice
    integer              :: id               ! Universe number for lattice
    character(len=104)   :: name = ""        ! User-defined name
    real(8), allocatable :: pitch(:)         ! Pitch along each axis
    integer, allocatable :: universes(:,:,:) ! Specified universes
    integer              :: outside          ! Material to fill area outside
    integer              :: outer            ! universe to tile outside the lat
    logical              :: is_3d            ! Lattice has cells on z axis
    integer, allocatable :: offset(:,:,:,:)  ! Distribcell offsets
  contains
    procedure(lattice_are_valid_indices_), deferred :: are_valid_indices
    procedure(lattice_get_indices_),       deferred :: get_indices
    procedure(lattice_get_local_xyz_),     deferred :: get_local_xyz
  end type Lattice

  abstract interface

!===============================================================================
! ARE_VALID_INDICES returns .true. if the given lattice indices fit within the
! bounds of the lattice.  Returns false otherwise.

    function lattice_are_valid_indices_(this, i_xyz) result(is_valid)
      import Lattice
      class(Lattice), intent(in) :: this
      integer,        intent(in) :: i_xyz(3)
      logical                    :: is_valid
    end function lattice_are_valid_indices_

!===============================================================================
! GET_INDICES returns the indices in a lattice for the given global xyz.

    function lattice_get_indices_(this, global_xyz) result(i_xyz)
      import Lattice
      class(Lattice), intent(in) :: this
      real(8),        intent(in) :: global_xyz(3)
      integer                    :: i_xyz(3)
    end function lattice_get_indices_

!===============================================================================
! GET_LOCAL_XYZ returns the translated local version of the given global xyz.

    function lattice_get_local_xyz_(this, global_xyz, i_xyz) result(local_xyz)
      import Lattice
      class(Lattice), intent(in) :: this
      real(8),        intent(in) :: global_xyz(3)
      integer,        intent(in) :: i_xyz(3)
      real(8)                    :: local_xyz(3)
    end function lattice_get_local_xyz_
  end interface

!===============================================================================
! RECTLATTICE extends LATTICE for rectilinear arrays.
!===============================================================================

  type, extends(Lattice) :: RectLattice
    integer              :: n_cells(3)     ! Number of cells along each axis
    real(8), allocatable :: lower_left(:)  ! Global lower-left corner of lat

    contains

    procedure :: are_valid_indices => valid_inds_rect
    procedure :: get_indices => get_inds_rect
    procedure :: get_local_xyz => get_local_rect
  end type RectLattice

!===============================================================================
! HEXLATTICE extends LATTICE for hexagonal (sometimes called triangular) arrays.
!===============================================================================

  type, extends(Lattice) :: HexLattice
    integer              :: n_rings   ! Number of radial ring cell positoins
    integer              :: n_axial   ! Number of axial cell positions
    real(8), allocatable :: center(:) ! Global center of lattice

    contains

    procedure :: are_valid_indices => valid_inds_hex
    procedure :: get_indices => get_inds_hex
    procedure :: get_local_xyz => get_local_hex
  end type HexLattice

!===============================================================================
! LATTICECONTAINER pointer array for storing lattices
!===============================================================================

  type LatticeContainer
    class(Lattice), allocatable :: obj
  end type LatticeContainer

!===============================================================================
! CELL defines a closed volume by its bounding surfaces
!===============================================================================

  type Cell
    integer :: id                          ! Unique ID
    character(len=104) :: name = ""        ! User-defined name
    integer :: type                        ! Type of cell (normal, universe,
                                           !  lattice)
    integer :: universe                    ! universe # this cell is in
    integer :: fill                        ! universe # filling this cell
    integer :: instances                   ! number of instances of this cell in
                                           !  the geom
    integer, allocatable :: material(:)    ! Material within cell.  Multiple
                                           !  materials for distribcell
                                           !  instances.  0 signifies a universe
    integer, allocatable :: offset(:)      ! Distribcell offset for tally
                                           !  counter
    integer, allocatable :: region(:)      ! Definition of spatial region as
                                           !  Boolean expression of half-spaces
    integer, allocatable :: rpn(:)         ! Reverse Polish notation for region
                                           !  expression
    logical :: simple                      ! Is the region simple (intersections
                                           !  only)
    integer :: distribcell_index           ! Index corresponding to this cell in
                                           !  distribcell arrays
    real(8), allocatable :: sqrtkT(:)      ! Square root of k_Boltzmann *
                                           !  temperature in eV.  Multiple for
                                           !  distribcell

    ! Rotation matrix and translation vector
    real(8), allocatable :: translation(:)
    real(8), allocatable :: rotation(:)
    real(8), allocatable :: rotation_matrix(:,:)
  end type Cell

  ! array index of the root universe
  integer :: root_universe = -1

  integer(C_INT32_T), bind(C) :: n_cells     ! # of cells
  integer(C_INT32_T), bind(C) :: n_universes ! # of universes
  integer(C_INT32_T), bind(C) :: n_lattices  ! # of lattices

  type(Cell),             allocatable, target :: cells(:)
  type(Universe),         allocatable, target :: universes(:)
  type(LatticeContainer), allocatable, target :: lattices(:)

  ! Dictionaries which map user IDs to indices in the global arrays
  type(DictIntInt) :: cell_dict
  type(DictIntInt) :: universe_dict
  type(DictIntInt) :: lattice_dict

contains

!===============================================================================

  function valid_inds_rect(this, i_xyz) result(is_valid)
    class(RectLattice), intent(in) :: this
    integer,            intent(in) :: i_xyz(3)
    logical                        :: is_valid

    is_valid = all(i_xyz > 0 .and. i_xyz <= this % n_cells)
  end function valid_inds_rect

!===============================================================================

  function valid_inds_hex(this, i_xyz) result(is_valid)
    class(HexLattice), intent(in) :: this
    integer,           intent(in) :: i_xyz(3)
    logical                       :: is_valid

    is_valid = (all(i_xyz > 0) .and. &
               &i_xyz(1) < 2*this % n_rings .and. &
               &i_xyz(2) < 2*this % n_rings .and. &
               &i_xyz(1) + i_xyz(2) > this % n_rings .and. &
               &i_xyz(1) + i_xyz(2) < 3*this % n_rings .and. &
               &i_xyz(3) <= this % n_axial)
  end function valid_inds_hex

!===============================================================================

  function get_inds_rect(this, global_xyz) result(i_xyz)
    class(RectLattice), intent(in) :: this
    real(8),            intent(in) :: global_xyz(3)
    integer                        :: i_xyz(3)

    real(8) :: xyz(3)  ! global_xyz alias

    xyz = global_xyz

    i_xyz(1) = ceiling((xyz(1) - this % lower_left(1))/this % pitch(1))
    i_xyz(2) = ceiling((xyz(2) - this % lower_left(2))/this % pitch(2))
    if (this % is_3d) then
      i_xyz(3) = ceiling((xyz(3) - this % lower_left(3))/this % pitch(3))
    else
      i_xyz(3) = 1
    end if
  end function get_inds_rect

!===============================================================================

  function get_inds_hex(this, global_xyz) result(i_xyz)
    class(HexLattice), intent(in) :: this
    real(8),           intent(in) :: global_xyz(3)
    integer                       :: i_xyz(3)

    real(8) :: xyz(3)    ! global xyz relative to the center
    real(8) :: alpha     ! Skewed coord axis
    real(8) :: xyz_t(3)  ! Local xyz
    real(8) :: d, d_min  ! Squared distance from cell centers
    integer :: i, j, k   ! Iterators
    integer :: k_min     ! Minimum distance index

    xyz(1) = global_xyz(1) - this % center(1)
    xyz(2) = global_xyz(2) - this % center(2)

    ! Index z direction.
    if (this % is_3d) then
      xyz(3) = global_xyz(3) - this % center(3)
      i_xyz(3) = ceiling(xyz(3)/this % pitch(2) + HALF*this % n_axial)
    else
      xyz(3) = global_xyz(3)
      i_xyz(3) = 1
    end if

    ! Convert coordinates into skewed bases.  The (x, alpha) basis is used to
    ! find the index of the global coordinates to within 4 cells.
    alpha = xyz(2) - xyz(1) / sqrt(THREE)
    i_xyz(1) = floor(xyz(1) / (sqrt(THREE) / TWO * this % pitch(1)))
    i_xyz(2) = floor(alpha / this % pitch(1))

    ! Add offset to indices (the center cell is (i_x, i_alpha) = (0, 0) but
    ! the array is offset so that the indices never go below 1).
    i_xyz(1) = i_xyz(1) + this % n_rings
    i_xyz(2) = i_xyz(2) + this % n_rings

    ! Calculate the (squared) distance between the particle and the centers of
    ! the four possible cells.  Regular hexagonal tiles form a centroidal
    ! Voronoi tessellation so the global xyz should be in the hexagonal cell
    ! that it is closest to the center of.  This method is used over a
    ! method that uses the remainders of the floor divisions above because it
    ! provides better finite precision performance.  Squared distances are
    ! used becasue they are more computationally efficient than normal
    ! distances.
    k = 1
    d_min = INFINITY
    do i = 0, 1
      do j = 0, 1
        xyz_t = this % get_local_xyz(global_xyz, i_xyz + [j, i, 0])
        d = xyz_t(1)**2 + xyz_t(2)**2
        if (d < d_min) then
          d_min = d
          k_min = k
        end if
        k = k + 1
      end do
    end do

    ! Select the minimum squared distance which corresponds to the cell the
    ! coordinates are in.
    if (k_min == 2) then
      i_xyz(1) = i_xyz(1) + 1
    else if (k_min == 3) then
      i_xyz(2) = i_xyz(2) + 1
    else if (k_min == 4) then
      i_xyz(1) = i_xyz(1) + 1
      i_xyz(2) = i_xyz(2) + 1
    end if
  end function get_inds_hex

!===============================================================================

  function get_local_rect(this, global_xyz, i_xyz) result(local_xyz)
    class(RectLattice), intent(in) :: this
    real(8),            intent(in) :: global_xyz(3)
    integer,            intent(in) :: i_xyz(3)
    real(8)                        :: local_xyz(3)

    real(8) :: xyz(3)  ! global_xyz alias

    xyz = global_xyz

    local_xyz(1) = xyz(1) - (this % lower_left(1) + &
         (i_xyz(1) - HALF)*this % pitch(1))
    local_xyz(2) = xyz(2) - (this % lower_left(2) + &
         (i_xyz(2) - HALF)*this % pitch(2))
    if (this % is_3d) then
      local_xyz(3) = xyz(3) - (this % lower_left(3) + &
           (i_xyz(3) - HALF)*this % pitch(3))
    else
      local_xyz(3) = xyz(3)
    end if
  end function get_local_rect

!===============================================================================

  function get_local_hex(this, global_xyz, i_xyz) result(local_xyz)
    class(HexLattice), intent(in) :: this
    real(8),           intent(in) :: global_xyz(3)
    integer,           intent(in) :: i_xyz(3)
    real(8)                       :: local_xyz(3)

    real(8) :: xyz(3)  ! global_xyz alias

    xyz = global_xyz

    ! x_l = x_g - (center + pitch_x*cos(30)*index_x)
    local_xyz(1) = xyz(1) - (this % center(1) + &
         sqrt(THREE) / TWO * (i_xyz(1) - this % n_rings) * this % pitch(1))
    ! y_l = y_g - (center + pitch_x*index_x + pitch_y*sin(30)*index_y)
    local_xyz(2) = xyz(2) - (this % center(2) + &
         (i_xyz(2) - this % n_rings) * this % pitch(1) + &
         (i_xyz(1) - this % n_rings) * this % pitch(1) / TWO)
    if (this % is_3d) then
      local_xyz(3) = xyz(3) - this % center(3) &
           + (HALF*this % n_axial - i_xyz(3) + HALF) * this % pitch(2)
    else
      local_xyz(3) = xyz(3)
    end if
  end function get_local_hex

!===============================================================================
! GET_TEMPERATURES returns a list of temperatures that each nuclide/S(a,b) table
! appears at in the model. Later, this list is used to determine the actual
! temperatures to read (which may be different if interpolation is used)
!===============================================================================

  subroutine get_temperatures(nuc_temps, sab_temps)
    type(VectorReal),            allocatable, intent(out) :: nuc_temps(:)
    type(VectorReal),  optional, allocatable, intent(out) :: sab_temps(:)

    integer :: i, j, k
    integer :: i_nuclide    ! index in nuclides array
    integer :: i_sab        ! index in S(a,b) array
    integer :: i_material
    real(8) :: temperature  ! temperature in Kelvin

    allocate(nuc_temps(n_nuclides))
    if (present(sab_temps)) allocate(sab_temps(n_sab_tables))

    do i = 1, size(cells)
      do j = 1, size(cells(i) % material)
        ! Skip any non-material cells and void materials
        if (cells(i) % material(j) == NONE .or. &
             cells(i) % material(j) == MATERIAL_VOID) cycle

        ! Get temperature of cell (rounding to nearest integer)
        if (size(cells(i) % sqrtkT) > 1) then
          temperature = cells(i) % sqrtkT(j)**2 / K_BOLTZMANN
        else
          temperature = cells(i) % sqrtkT(1)**2 / K_BOLTZMANN
        end if

        i_material = cells(i) % material(j)

        associate (mat => materials(i_material))
          NUC_NAMES_LOOP: do k = 1, size(mat % names)
            ! Get index in nuc_temps array
            i_nuclide = nuclide_dict % get(to_lower(mat % names(k)))

            ! Add temperature if it hasn't already been added
            if (find(nuc_temps(i_nuclide), temperature) == -1) then
              call nuc_temps(i_nuclide) % push_back(temperature)
            end if
          end do NUC_NAMES_LOOP

          if (present(sab_temps) .and. mat % n_sab > 0) then
            SAB_NAMES_LOOP: do k = 1, size(mat % sab_names)
              ! Get index in nuc_temps array
              i_sab = sab_dict % get(to_lower(mat % sab_names(k)))

              ! Add temperature if it hasn't already been added
              if (find(sab_temps(i_sab), temperature) == -1) then
                call sab_temps(i_sab) % push_back(temperature)
              end if
            end do SAB_NAMES_LOOP
          end if
        end associate
      end do
    end do

  end subroutine get_temperatures

!===============================================================================
! FREE_MEMORY_GEOMETRY deallocates global arrays defined in this module
!===============================================================================

  subroutine free_memory_geometry()

    n_cells = 0
    n_universes = 0
    n_lattices = 0

    if (allocated(cells)) deallocate(cells)
    if (allocated(universes)) deallocate(universes)
    if (allocated(lattices)) deallocate(lattices)

    call cell_dict % clear()
    call universe_dict % clear()
    call lattice_dict % clear()

  end subroutine free_memory_geometry

!===============================================================================
!                               C API FUNCTIONS
!===============================================================================

  function openmc_extend_cells(n, index_start, index_end) result(err) bind(C)
    ! Extend the cells array by n elements
    integer(C_INT32_T), value, intent(in) :: n
    integer(C_INT32_T), optional, intent(out) :: index_start
    integer(C_INT32_T), optional, intent(out) :: index_end
    integer(C_INT) :: err

    type(Cell), allocatable :: temp(:) ! temporary cells array

    if (n_cells == 0) then
      ! Allocate cells array
      allocate(cells(n))
    else
      ! Allocate cells array with increased size
      allocate(temp(n_cells + n))

      ! Copy original cells to temporary array
      temp(1:n_cells) = cells

      ! Move allocation from temporary array
      call move_alloc(FROM=temp, TO=cells)
    end if

    ! Return indices in cells array
    if (present(index_start)) index_start = n_cells + 1
    if (present(index_end)) index_end = n_cells + n
    n_cells = n_cells + n

    err = 0
  end function openmc_extend_cells


  function openmc_get_cell_index(id, index) result(err) bind(C)
    ! Return the index in the cells array of a cell with a given ID
    integer(C_INT32_T), value :: id
    integer(C_INT32_T), intent(out) :: index
    integer(C_INT) :: err

    if (allocated(cells)) then
      if (cell_dict % has(id)) then
        index = cell_dict % get(id)
        err = 0
      else
        err = E_INVALID_ID
        call set_errmsg("No cell exists with ID=" // trim(to_str(id)) // ".")
      end if
    else
      err = E_ALLOCATE
      call set_errmsg("Memory has not been allocated for cells.")
    end if
  end function openmc_get_cell_index


  function openmc_cell_get_fill(index, type, indices, n) result(err) bind(C)
    integer(C_INT32_T), value, intent(in) :: index
    integer(C_INT), intent(out) :: type
    integer(C_INT32_T), intent(out) :: n
    type(C_PTR), intent(out) :: indices
    integer(C_INT) :: err

    err = 0
    if (index >= 1 .and. index <= size(cells)) then
      associate (c => cells(index))
        type = c % type
        select case (type)
        case (FILL_MATERIAL)
          n = size(c % material)
          indices = C_LOC(c % material(1))
        case (FILL_UNIVERSE, FILL_LATTICE)
          n = 1
          indices = C_LOC(c % fill)
        end select
      end associate
    else
      err = E_OUT_OF_BOUNDS
      call set_errmsg("Index in cells array is out of bounds.")
    end if
  end function openmc_cell_get_fill


  function openmc_cell_get_id(index, id) result(err) bind(C)
    ! Return the ID of a cell
    integer(C_INT32_T), value       :: index
    integer(C_INT32_T), intent(out) :: id
    integer(C_INT) :: err

    if (index >= 1 .and. index <= size(cells)) then
      id = cells(index) % id
      err = 0
    else
      err = E_OUT_OF_BOUNDS
      call set_errmsg("Index in cells array is out of bounds.")
    end if
  end function openmc_cell_get_id


  function openmc_cell_set_fill(index, type, n, indices) result(err) bind(C)
    ! Set the fill for a fill
    integer(C_INT32_T), value, intent(in) :: index    ! index in cells
    integer(C_INT), value, intent(in)     :: type
    integer(c_INT32_T), value, intent(in) :: n
    integer(C_INT32_T), intent(in) :: indices(n)
    integer(C_INT) :: err

    integer :: i, j

    err = 0
    if (index >= 1 .and. index <= size(cells)) then
      associate (c => cells(index))
        select case (type)
        case (FILL_MATERIAL)
          if (allocated(c % material)) deallocate(c % material)
          allocate(c % material(n))

          c % type = FILL_MATERIAL
          do i = 1, n
            j = indices(i)
            if (j == 0) then
              c % material(i) = MATERIAL_VOID
            else
              if (j >= 1 .and. j <= n_materials) then
                c % material(i) = j
              else
                err = E_OUT_OF_BOUNDS
                call set_errmsg("Index " // trim(to_str(j)) // " in the &
                     &materials array is out of bounds.")
              end if
            end if
          end do
        case (FILL_UNIVERSE)
          c % type = FILL_UNIVERSE
        case (FILL_LATTICE)
          c % type = FILL_LATTICE
        end select
      end associate
    else
      err = E_OUT_OF_BOUNDS
      call set_errmsg("Index in cells array is out of bounds.")
    end if

  end function openmc_cell_set_fill


  function openmc_cell_set_id(index, id) result(err) bind(C)
    ! Set the ID of a cell
    integer(C_INT32_T), value, intent(in) :: index
    integer(C_INT32_T), value, intent(in) :: id
    integer(C_INT) :: err

    if (index >= 1 .and. index <= n_cells) then
      cells(index) % id = id
      call cell_dict % set(id, index)
      err = 0
    else
      err = E_OUT_OF_BOUNDS
      call set_errmsg("Index in cells array is out of bounds.")
    end if
  end function openmc_cell_set_id


  function openmc_cell_set_temperature(index, T, instance) result(err) bind(C)
    ! Set the temperature of a cell
    integer(C_INT32_T), value, intent(in)    :: index    ! index in cells
    real(C_DOUBLE), value, intent(in)        :: T        ! temperature
    integer(C_INT32_T), optional, intent(in) :: instance ! cell instance

    integer(C_INT) :: err     ! error code
    integer :: j              ! looping variable
    integer :: n              ! number of cell instances
    integer :: material_index ! material index in materials array
    integer :: num_nuclides   ! num nuclides in material
    integer :: nuclide_index  ! index of nuclide in nuclides array
    real(8) :: min_temp       ! min common-denominator avail temp
    real(8) :: max_temp       ! max common-denominator avail temp
    real(8) :: temp           ! actual temp we'll assign
    logical :: outside_low    ! lower than available data
    logical :: outside_high   ! higher than available data

    outside_low = .false.
    outside_high = .false.

    err = E_UNASSIGNED

    if (index >= 1 .and. index <= size(cells)) then

      ! error if the cell is filled with another universe
      if (cells(index) % fill /= NONE) then
        err = E_GEOMETRY
        call set_errmsg("Cannot set temperature on a cell filled &
             &with a universe.")
      else
        ! find which material is associated with this cell (material_index
        ! is the index into the materials array)
        if (present(instance)) then
          material_index = cells(index) % material(instance)
        else
          material_index = cells(index) % material(1)
        end if

        ! number of nuclides associated with this material
        num_nuclides = size(materials(material_index) % nuclide)

        min_temp = ZERO
        max_temp = INFINITY

        do j = 1, num_nuclides
          nuclide_index = materials(material_index) % nuclide(j)
          min_temp = max(min_temp, minval(nuclides(nuclide_index) % kTs))
          max_temp = min(max_temp, maxval(nuclides(nuclide_index) % kTs))
        end do

        ! adjust the temperature to be within bounds if necessary
        if (K_BOLTZMANN * T < min_temp) then
          outside_low = .true.
          temp = min_temp / K_BOLTZMANN
        else if (K_BOLTZMANN * T > max_temp) then
          outside_high = .true.
          temp = max_temp / K_BOLTZMANN
        else
          temp = T
        end if

        associate (c => cells(index))
          if (allocated(c % sqrtkT)) then
            n = size(c % sqrtkT)
            if (present(instance) .and. n > 1) then
              if (instance >= 0 .and. instance < n) then
                c % sqrtkT(instance + 1) = sqrt(K_BOLTZMANN * temp)
                err = 0
              end if
            else
              c % sqrtkT(:) = sqrt(K_BOLTZMANN * temp)
              err = 0
            end if
          end if
        end associate

        ! Assign error codes for outside of temperature bounds provided the
        ! temperature was changed correctly. This needs to be done after
        ! changing the temperature based on the logical structure above.
        if (err == 0) then
          if (outside_low) then
            err = E_WARNING
            call set_errmsg("Nuclear data has not been loaded beyond lower &
                 &bound of T=" // trim(to_str(T)) // " K.")
          else if (outside_high) then
            err = E_WARNING
            call set_errmsg("Nuclear data has not been loaded beyond upper &
                 &bound of T=" // trim(to_str(T)) // " K.")
          end if
        end if

      end if

    else
      err = E_OUT_OF_BOUNDS
      call set_errmsg("Index in cells array is out of bounds.")
    end if
  end function openmc_cell_set_temperature

end module geometry_header
