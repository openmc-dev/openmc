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
     integer :: id       ! Universe number for lattice
     integer :: type     ! Type of lattice (rectangular, hex, etc)
     integer :: level    ! Level of lattice
     integer :: n_x      ! number of lattice cells in x-direction
     integer :: n_y      ! number of lattice cells in y-direction
     real(8) :: x0       ! x-coordinate of lattice origin
     real(8) :: y0       ! y-coordinate of lattice origin
     real(8) :: width_x  ! width of lattice cell 
     real(8) :: width_y  ! width of lattice cell
     integer, allocatable :: element(:,:) ! specified universes
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
  end type Cell

  ! array index of universe 0
  integer :: BASE_UNIVERSE

end module geometry_header
