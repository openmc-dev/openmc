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

  type Lattice
     integer :: id          ! Universe number for lattice
     integer :: type        ! Type of lattice (rectangular, hex, etc)
     integer :: level       ! Level of lattice
     integer :: n_dimension ! Number of dimensions
     integer, allocatable :: dimension(:)     ! number of cells in each direction
     real(8), allocatable :: lower_left(:)    ! lower-left corner of lattice
     real(8), allocatable :: width(:)         ! width of each lattice cell
     integer, allocatable :: universes(:,:,:) ! specified universes
     integer              :: outside          ! material to fill area outside
  end type Lattice

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
