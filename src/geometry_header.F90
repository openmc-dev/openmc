module geometry_header

  use, intrinsic :: ISO_C_BINDING

  use algorithm,       only: find
  use constants,       only: HALF, TWO, THREE, INFINITY, K_BOLTZMANN, &
                             MATERIAL_VOID
  use dict_header,     only: DictCharInt, DictIntInt
  use hdf5_interface,  only: HID_T
  use material_header, only: Material, materials, material_dict, n_materials
  use nuclide_header
  use sab_header
  use stl_vector,      only: VectorReal
  use string,          only: to_lower

  implicit none

  interface
    function cell_pointer(cell_ind) bind(C) result(ptr)
      import C_PTR, C_INT32_T
      integer(C_INT32_T), intent(in), value :: cell_ind
      type(C_PTR)                           :: ptr
    end function cell_pointer

    function cell_id_c(cell_ptr) bind(C, name='cell_id') result(id)
      import C_PTR, C_INT32_T
      type(C_PTR), intent(in), value :: cell_ptr
      integer(C_INT32_T)             :: id
    end function cell_id_c

    subroutine cell_set_id_c(cell_ptr, id) bind(C, name='cell_set_id')
      import C_PTR, C_INT32_T
      type(C_PTR),        intent(in), value :: cell_ptr
      integer(C_INT32_T), intent(in), value :: id
    end subroutine cell_set_id_c

    function cell_type_c(cell_ptr) bind(C, name='cell_type') result(type)
      import C_PTR, C_INT
      type(C_PTR), intent(in), value :: cell_ptr
      integer(C_INT)                 :: type
    end function cell_type_c

    function cell_universe_c(cell_ptr) bind(C, name='cell_universe') &
         result(universe)
      import C_PTR, C_INT32_T
      type(C_PTR), intent(in), value :: cell_ptr
      integer(C_INT32_T)             :: universe
    end function cell_universe_c

    function cell_fill_c(cell_ptr) bind(C, name="cell_fill") result(fill)
      import C_PTR, C_INT32_T
      type(C_PTR), intent(in), value :: cell_ptr
      integer(C_INT32_T)             :: fill
    end function cell_fill_c

    function cell_n_instances_c(cell_ptr) bind(C, name='cell_n_instances') &
         result(n_instances)
      import C_PTR, C_INT32_T
      type(C_PTR), intent(in), value :: cell_ptr
      integer(C_INT32_T)             :: n_instances
    end function cell_n_instances_c

    function cell_material_size_c(cell_ptr) bind(C, name='cell_material_size') &
         result(n)
      import C_PTR, C_INT
      type(C_PTR), intent(in), value :: cell_ptr
      integer(C_INT)                 :: n
    end function cell_material_size_c

    function cell_material_c(cell_ptr, i) bind(C, name='cell_material') &
         result(mat)
      import C_PTR, C_INT, C_INT32_T
      type(C_PTR),    intent(in), value :: cell_ptr
      integer(C_INT), intent(in), value :: i
      integer(C_INT32_T)                :: mat
    end function cell_material_c

    function cell_simple_c(cell_ptr) bind(C, name='cell_simple') result(simple)
      import C_PTR, C_BOOL
      type(C_PTR), intent(in), value :: cell_ptr
      logical(C_BOOL)                :: simple
    end function cell_simple_c

    subroutine cell_distance_c(cell_ptr, xyz, uvw, on_surface, min_dist, &
                               i_surf) bind(C, name="cell_distance")
      import C_PTR, C_INT32_T, C_DOUBLE
      type(C_PTR),        intent(in), value :: cell_ptr
      real(C_DOUBLE),     intent(in)        :: xyz(3)
      real(C_DOUBLE),     intent(in)        :: uvw(3)
      integer(C_INT32_T), intent(in), value :: on_surface
      real(C_DOUBLE),     intent(out)       :: min_dist
      integer(C_INT32_T), intent(out)       :: i_surf
    end subroutine cell_distance_c

    function cell_offset_c(cell_ptr, map) bind(C, name="cell_offset") &
         result(offset)
      import C_PTR, C_INT, C_INT32_T
      type(C_PTR),    intent(in), value :: cell_ptr
      integer(C_INT), intent(in), value :: map
      integer(C_INT32_T)                :: offset
    end function cell_offset_c

    subroutine cell_to_hdf5_c(cell_ptr, group) bind(C, name='cell_to_hdf5')
      import HID_T, C_PTR
      type(C_PTR),    intent(in), value :: cell_ptr
      integer(HID_T), intent(in), value :: group
    end subroutine cell_to_hdf5_c

    function lattice_pointer(lat_ind) bind(C) result(ptr)
      import C_PTR, C_INT32_T
      integer(C_INT32_T), intent(in), value :: lat_ind
      type(C_PTR)                           :: ptr
    end function lattice_pointer

    function lattice_id_c(lat_ptr) bind(C, name='lattice_id') result(id)
      import C_PTR, C_INT32_T
      type(C_PTR), intent(in), value :: lat_ptr
      integer(C_INT32_T)             :: id
    end function lattice_id_c

    function lattice_are_valid_indices_c(lat_ptr, i_xyz) &
         bind(C, name='lattice_are_valid_indices') result (is_valid)
      import C_PTR, C_INT, C_BOOL
      type(C_PTR),    intent(in), value :: lat_ptr
      integer(C_INT), intent(in)        :: i_xyz(3)
      logical(C_BOOL)                   :: is_valid
    end function lattice_are_valid_indices_c

    subroutine lattice_distance_c(lat_ptr, xyz, uvw, i_xyz, d, lattice_trans) &
         bind(C, name='lattice_distance')
      import C_PTR, C_INT, C_DOUBLE
      type(C_PTR),    intent(in), value :: lat_ptr
      real(C_DOUBLE), intent(in)        :: xyz(3)
      real(C_DOUBLE), intent(in)        :: uvw(3)
      integer(C_INT), intent(in)        :: i_xyz(3)
      real(C_DOUBLE), intent(out)       :: d
      integer(C_INT), intent(out)       :: lattice_trans(3)
    end subroutine lattice_distance_c

    subroutine lattice_get_indices_c(lat_ptr, xyz, i_xyz) &
         bind(C, name='lattice_get_indices')
      import C_PTR, C_INT, C_DOUBLE
      type(C_PTR),    intent(in), value :: lat_ptr
      real(C_DOUBLE), intent(in)        :: xyz(3)
      integer(C_INT), intent(out)       :: i_xyz(3)
    end subroutine lattice_get_indices_c

    subroutine lattice_get_local_xyz_c(lat_ptr, global_xyz, i_xyz, local_xyz) &
         bind(C, name='lattice_get_local_xyz')
      import C_PTR, C_INT, C_DOUBLE
      type(C_PTR),    intent(in), value :: lat_ptr
      real(C_DOUBLE), intent(in)        :: global_xyz(3)
      integer(C_INT), intent(in)        :: i_xyz(3)
      real(C_DOUBLE), intent(out)       :: local_xyz(3)
    end subroutine lattice_get_local_xyz_c

    subroutine lattice_to_hdf5_c(lat_ptr, group) bind(C, name='lattice_to_hdf5')
      import HID_T, C_PTR
      type(C_PTR),    intent(in), value :: lat_ptr
      integer(HID_T), intent(in), value :: group
    end subroutine lattice_to_hdf5_c

    function lattice_offset_c(lat_ptr, map, i_xyz) &
         bind(C, name='lattice_offset') result(offset)
      import C_PTR, C_INT, C_INT32_T
      type(C_PTR),    intent(in), value :: lat_ptr
      integer(C_INT), intent(in), value :: map
      integer(C_INT), intent(in)        :: i_xyz(3)
      integer(C_INT32_T)                :: offset
    end function lattice_offset_c

    function lattice_outer_c(lat_ptr) bind(C, name='lattice_outer') &
         result(outer)
      import C_PTR, C_INT32_T
      type(C_PTR), intent(in), value :: lat_ptr
      integer(C_INT32_T)             :: outer
    end function lattice_outer_c

    function lattice_universe_c(lat_ptr, i_xyz) &
         bind(C, name='lattice_universe') result(univ)
      import C_PTR, C_INT32_t, C_INT
      type(C_PTR),    intent(in), value :: lat_ptr
      integer(C_INT), intent(in)        :: i_xyz(3)
      integer(C_INT32_T)                :: univ
    end function lattice_universe_c

    subroutine extend_cells_c(n) bind(C)
      import C_INT32_t
      integer(C_INT32_T), intent(in), value :: n
    end subroutine extend_cells_c
  end interface

!===============================================================================
! UNIVERSE defines a geometry that fills all phase space
!===============================================================================

  type Universe
    integer :: id                     ! Unique ID
    integer, allocatable :: cells(:)  ! List of cells within
    real(8) :: x0                     ! Translation in x-coordinate
    real(8) :: y0                     ! Translation in y-coordinate
    real(8) :: z0                     ! Translation in z-coordinate
  end type Universe

!===============================================================================
! LATTICE abstract type for ordered array of universes.
!===============================================================================

  type, abstract :: Lattice
    type(C_PTR) :: ptr
  contains
    procedure :: id => lattice_id
    procedure :: are_valid_indices => lattice_are_valid_indices
    procedure :: distance => lattice_distance
    procedure :: get => lattice_get
    procedure :: get_indices => lattice_get_indices
    procedure :: get_local_xyz => lattice_get_local_xyz
    procedure :: offset => lattice_offset
    procedure :: outer => lattice_outer
    procedure :: to_hdf5 => lattice_to_hdf5
  end type Lattice

!===============================================================================
! RECTLATTICE extends LATTICE for rectilinear arrays.
!===============================================================================

  type, extends(Lattice) :: RectLattice
  end type RectLattice

!===============================================================================
! HEXLATTICE extends LATTICE for hexagonal (sometimes called triangular) arrays.
!===============================================================================

  type, extends(Lattice) :: HexLattice
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
    type(C_PTR) :: ptr

    integer, allocatable :: region(:)      ! Definition of spatial region as
                                           !  Boolean expression of half-spaces
    integer :: distribcell_index           ! Index corresponding to this cell in
                                           !  distribcell arrays
    real(8), allocatable :: sqrtkT(:)      ! Square root of k_Boltzmann *
                                           !  temperature in eV.  Multiple for
                                           !  distribcell

    ! Rotation matrix and translation vector
    real(8), allocatable :: translation(:)
    real(8), allocatable :: rotation(:)
    real(8), allocatable :: rotation_matrix(:,:)

  contains

    procedure :: id => cell_id
    procedure :: set_id => cell_set_id
    procedure :: type => cell_type
    procedure :: universe => cell_universe
    procedure :: fill => cell_fill
    procedure :: n_instances => cell_n_instances
    procedure :: material_size => cell_material_size
    procedure :: material => cell_material
    procedure :: simple => cell_simple
    procedure :: distance => cell_distance
    procedure :: offset => cell_offset
    procedure :: to_hdf5 => cell_to_hdf5

  end type Cell

  ! array index of the root universe
  integer(C_INT), bind(C, name='openmc_root_universe') :: root_universe = -1

  integer(C_INT32_T), bind(C) :: n_cells     ! # of cells
  integer(C_INT32_T), bind(C) :: n_universes ! # of universes

  type(Cell),             allocatable, target :: cells(:)
  type(Universe),         allocatable, target :: universes(:)
  type(LatticeContainer), allocatable, target :: lattices(:)

  ! Dictionaries which map user IDs to indices in the global arrays
  type(DictIntInt) :: cell_dict
  type(DictIntInt) :: universe_dict
  type(DictIntInt) :: lattice_dict

contains

  function lattice_id(this) result(id)
    class(Lattice), intent(in) :: this
    integer(C_INT32_T)         :: id
    id = lattice_id_c(this % ptr)
  end function lattice_id

  function lattice_are_valid_indices(this, i_xyz) result (is_valid)
    class(Lattice), intent(in) :: this
    integer(C_INT), intent(in) :: i_xyz(3)
    logical(C_BOOL)            :: is_valid
    is_valid = lattice_are_valid_indices_c(this % ptr, i_xyz)
  end function lattice_are_valid_indices

  subroutine lattice_distance(this, xyz, uvw, i_xyz, d, lattice_trans)
    class(Lattice), intent(in)  :: this
    real(C_DOUBLE), intent(in)  :: xyz(3)
    real(C_DOUBLE), intent(in)  :: uvw(3)
    integer(C_INT), intent(in)  :: i_xyz(3)
    real(C_DOUBLE), intent(out) :: d
    integer(C_INT), intent(out) :: lattice_trans(3)
    call lattice_distance_c(this % ptr, xyz, uvw, i_xyz, d, lattice_trans)
  end subroutine lattice_distance

  function lattice_get(this, i_xyz) result(univ)
    class(Lattice), intent(in)  :: this
    integer(C_INT), intent(in)  :: i_xyz(3)
    integer(C_INT32_T)          :: univ
    univ = lattice_universe_c(this % ptr, i_xyz)
  end function lattice_get

  function lattice_get_indices(this, xyz) result(i_xyz)
    class(Lattice), intent(in)  :: this
    real(C_DOUBLE), intent(in)  :: xyz(3)
    integer(C_INT)              :: i_xyz(3)
    call lattice_get_indices_c(this % ptr, xyz, i_xyz)
  end function lattice_get_indices

  function lattice_get_local_xyz(this, global_xyz, i_xyz) &
       result(local_xyz)
    class(Lattice), intent(in) :: this
    real(C_DOUBLE), intent(in) :: global_xyz(3)
    integer(C_INT), intent(in) :: i_xyz(3)
    real(C_DOUBLE)             :: local_xyz(3)
    call lattice_get_local_xyz_c(this % ptr, global_xyz, i_xyz, local_xyz)
  end function lattice_get_local_xyz

  function lattice_offset(this, map, i_xyz) result(offset)
    class(Lattice), intent(in) :: this
    integer(C_INT), intent(in) :: map
    integer(C_INT), intent(in) :: i_xyz(3)
    integer(C_INT32_T)         :: offset
    offset = lattice_offset_c(this % ptr, map, i_xyz)
  end function lattice_offset

  function lattice_outer(this) result(outer)
    class(Lattice), intent(in) :: this
    integer(C_INT32_T)         :: outer
    outer = lattice_outer_c(this % ptr)
  end function lattice_outer

  subroutine lattice_to_hdf5(this, group)
    class(Lattice), intent(in) :: this
    integer(HID_T), intent(in) :: group
    call lattice_to_hdf5_c(this % ptr, group)
  end subroutine lattice_to_hdf5

!===============================================================================

  function cell_id(this) result(id)
    class(Cell), intent(in) :: this
    integer(C_INT32_T)      :: id
    id = cell_id_c(this % ptr)
  end function cell_id

  subroutine cell_set_id(this, id)
    class(Cell),        intent(in) :: this
    integer(C_INT32_T), intent(in) :: id
    call cell_set_id_c(this % ptr, id)
  end subroutine cell_set_id

  function cell_type(this) result(type)
    class(Cell), intent(in) :: this
    integer(C_INT)          :: type
    type = cell_type_c(this % ptr)
  end function cell_type

  function cell_universe(this) result(universe)
    class(Cell), intent(in) :: this
    integer(C_INT32_T)      :: universe
    universe = cell_universe_c(this % ptr)
  end function cell_universe

  function cell_fill(this) result(fill)
    class(Cell), intent(in) :: this
    integer(C_INT32_T)      :: fill
    fill = cell_fill_c(this % ptr)
  end function cell_fill

  function cell_n_instances(this) result(n_instances)
    class(Cell), intent(in) :: this
    integer(C_INT32_T)      :: n_instances
    n_instances = cell_n_instances_c(this % ptr)
  end function cell_n_instances

  function cell_material_size(this) result(n)
    class(Cell), intent(in) :: this
    integer(C_INT)          :: n
    n = cell_material_size_c(this % ptr)
  end function cell_material_size

  function cell_material(this, i) result(mat)
    class(Cell), intent(in) :: this
    integer,     intent(in) :: i
    integer(C_INT32_T)      :: mat
    mat = cell_material_c(this % ptr, i)
  end function cell_material

  function cell_simple(this) result(simple)
    class(Cell), intent(in) :: this
    logical(C_BOOL)         :: simple
    simple = cell_simple_c(this % ptr)
  end function cell_simple

  subroutine cell_distance(this, xyz, uvw, on_surface, min_dist, i_surf)
    class(Cell),        intent(in)  :: this
    real(C_DOUBLE),     intent(in)  :: xyz(3)
    real(C_DOUBLE),     intent(in)  :: uvw(3)
    integer(C_INT32_T), intent(in)  :: on_surface
    real(C_DOUBLE),     intent(out) :: min_dist
    integer(C_INT32_T), intent(out) :: i_surf
    call cell_distance_c(this % ptr, xyz, uvw, on_surface, min_dist, i_surf)
  end subroutine cell_distance

  function cell_offset(this, map) result(offset)
    class(Cell), intent(in)    :: this
    integer(C_INT), intent(in) :: map
    integer(C_INT32_T)         :: offset
    offset = cell_offset_c(this % ptr, map)
  end function cell_offset

  subroutine cell_to_hdf5(this, group)
    class(Cell),    intent(in) :: this
    integer(HID_T), intent(in) :: group
    call cell_to_hdf5_c(this % ptr, group)
  end subroutine cell_to_hdf5

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
      ! Skip non-material cells.
      if (cells(i) % fill() /= C_NONE) cycle

      do j = 1, cells(i) % material_size()
        ! Skip void materials
        if (cells(i) % material(j) == MATERIAL_VOID) cycle

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
    interface
      subroutine free_memory_geometry_c() bind(C)
      end subroutine free_memory_geometry_c
    end interface

    call free_memory_geometry_c()

    n_cells = 0
    n_universes = 0

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
    integer(C_INT32_T) :: i
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

    ! Extend the C++ cells array and get pointers to the C++ objects
    call extend_cells_c(n)
    do i = n_cells - n, n_cells
      cells(i) % ptr = cell_pointer(i - 1)
    end do

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


  function openmc_cell_get_id(index, id) result(err) bind(C)
    ! Return the ID of a cell
    integer(C_INT32_T), value       :: index
    integer(C_INT32_T), intent(out) :: id
    integer(C_INT) :: err

    if (index >= 1 .and. index <= size(cells)) then
      id = cells(index) % id()
      err = 0
    else
      err = E_OUT_OF_BOUNDS
      call set_errmsg("Index in cells array is out of bounds.")
    end if
  end function openmc_cell_get_id


  function openmc_cell_set_id(index, id) result(err) bind(C)
    ! Set the ID of a cell
    integer(C_INT32_T), value, intent(in) :: index
    integer(C_INT32_T), value, intent(in) :: id
    integer(C_INT) :: err

    if (index >= 1 .and. index <= n_cells) then
      call cells(index) % set_id(id)
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
      if (cells(index) % fill() /= C_NONE) then
        err = E_GEOMETRY
        call set_errmsg("Cannot set temperature on a cell filled &
             &with a universe.")
      else
        ! find which material is associated with this cell (material_index
        ! is the index into the materials array)
        if (present(instance)) then
          material_index = cells(index) % material(instance + 1)
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
