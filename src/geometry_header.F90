module geometry_header

  use algorithm,       only: find
  use constants,       only: HALF, TWO, THREE, INFINITY, K_BOLTZMANN, &
                             MATERIAL_VOID, NONE
  use dict_header,     only: DictCharInt, DictIntInt
  use material_header, only: Material
  use stl_vector,      only: VectorReal
  use string,          only: to_lower

  implicit none

!===============================================================================
! UNIVERSE defines a geometry that fills all phase space
!===============================================================================

  type Universe
    integer :: id                     ! Unique ID
    integer :: type                   ! Type
    integer :: n_cells                ! # of cells within
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

  ! array index of universe 0
  integer :: BASE_UNIVERSE

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

  subroutine get_temperatures(cells, materials, material_dict, nuclide_dict, &
                              n_nucs, nuc_temps, sab_dict, n_sabs, sab_temps)
    type(Cell),                  allocatable, intent(in)  :: cells(:)
    type(Material),              allocatable, intent(in)  :: materials(:)
    type(DictIntInt),                         intent(in)  :: material_dict
    type(DictCharInt),                        intent(in)  :: nuclide_dict
    integer,                                  intent(in)  :: n_nucs
    type(VectorReal),            allocatable, intent(out) :: nuc_temps(:)
    type(DictCharInt), optional,              intent(in)  :: sab_dict
    integer,           optional,              intent(in)  :: n_sabs
    type(VectorReal),  optional, allocatable, intent(out) :: sab_temps(:)

    integer :: i, j, k
    integer :: i_nuclide    ! index in nuclides array
    integer :: i_sab        ! index in S(a,b) array
    integer :: i_material
    real(8) :: temperature  ! temperature in Kelvin

    allocate(nuc_temps(n_nucs))
    if (present(n_sabs) .and. present(sab_temps)) allocate(sab_temps(n_sabs))

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

        i_material = material_dict % get_key(cells(i) % material(j))
        associate (mat => materials(i_material))
          NUC_NAMES_LOOP: do k = 1, size(mat % names)
            ! Get index in nuc_temps array
            i_nuclide = nuclide_dict % get_key(to_lower(mat % names(k)))

            ! Add temperature if it hasn't already been added
            if (find(nuc_temps(i_nuclide), temperature) == -1) then
              call nuc_temps(i_nuclide) % push_back(temperature)
            end if
          end do NUC_NAMES_LOOP

          if (present(sab_temps) .and. present(sab_dict) .and. &
               mat % n_sab > 0) then
            SAB_NAMES_LOOP: do k = 1, size(mat % sab_names)
              ! Get index in nuc_temps array
              i_sab = sab_dict % get_key(to_lower(mat % sab_names(k)))

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

end module geometry_header
