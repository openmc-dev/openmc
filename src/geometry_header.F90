module geometry_header

  implicit none

!===============================================================================
! UNIVERSE defines a geometry that fills all phase space
!===============================================================================

  type Universe
     integer :: id                    ! Unique ID
     integer :: type                  ! Type
     integer :: n_cells               ! # of cells within
     integer, allocatable :: cells(:) ! List of cells within
     real(8) :: x0                    ! Translation in x-coordinate
     real(8) :: y0                    ! Translation in y-coordinate
     real(8) :: z0                    ! Translation in z-coordinate
  end type Universe

!===============================================================================
! LATTICE is an ordered array of elements (either rectangular, hexagonal, or
! triangular)
!===============================================================================

  type, abstract :: Lattice
    integer :: id          ! Universe number for lattice
    real(8), allocatable :: pitch(:)         ! Pitch along each axis
    integer, allocatable :: universes(:,:,:) ! Specified universes
    integer              :: outside          ! Material to fill area outside
    logical              :: is_3d            ! Lattice has cells on z axis
  end type Lattice

  type, extends(Lattice) :: RectLattice
    integer               :: n_cells(3)     ! Number of cells along each axis
    real(8), allocatable  :: lower_left(:)  ! Global lower-left corner of lat
  end type RectLattice

  type, extends(Lattice) :: HexLattice
    integer               :: n_rings   ! Number of radial ring cell positoins
    integer               :: n_axial   ! Number of axial cell positions
    real(8), allocatable  :: center(:) ! Global center of lattice
  end type HexLattice

  type LatticeContainer
    class(Lattice), allocatable :: obj
  end type LatticeContainer

!===============================================================================
! SURFACE type defines a first- or second-order surface that can be used to
! construct closed volumes (cells)
!===============================================================================

  type Surface
     integer :: id                     ! Unique ID
     integer :: type                   ! Type of surface
     real(8), allocatable :: coeffs(:) ! Definition of surface
     integer, allocatable :: & 
          neighbor_pos(:), &           ! List of cells on positive side
          neighbor_neg(:)              ! List of cells on negative side
     integer :: bc                     ! Boundary condition
  end type Surface

!===============================================================================
! CELL defines a closed volume by its bounding surfaces
!===============================================================================

  type Cell
     integer :: id         ! Unique ID
     integer :: type       ! Type of cell (normal, universe, lattice)
     integer :: universe   ! universe # this cell is in
     integer :: fill       ! universe # filling this cell
     integer :: material   ! Material within cell (0 for universe)
     integer :: n_surfaces ! Number of surfaces within
     integer, allocatable :: & 
          & surfaces(:)    ! List of surfaces bounding cell -- note that
                           ! parentheses, union, etc operators will be listed
                           ! here too

     ! Rotation matrix and translation vector
     real(8), allocatable :: rotation(:,:)
     real(8), allocatable :: translation(:)
  end type Cell

  ! array index of universe 0
  integer :: BASE_UNIVERSE

end module geometry_header
