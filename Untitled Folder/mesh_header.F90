module mesh_header

  implicit none

!===============================================================================
! STRUCTUREDMESH represents a tessellation of n-dimensional Euclidean space by
! congruent squares or cubes
!===============================================================================

  type StructuredMesh
    integer :: id                          ! user-specified id
    integer :: type                        ! rectangular, hexagonal
    integer :: n_dimension                 ! rank of mesh
    real(8) :: volume_frac                 ! volume fraction of each cell
    integer, allocatable :: dimension(:)   ! number of cells in each direction
    real(8), allocatable :: lower_left(:)  ! lower-left corner of mesh
    real(8), allocatable :: upper_right(:) ! upper-right corner of mesh
    real(8), allocatable :: width(:)       ! width of each mesh cell
  end type StructuredMesh

end module mesh_header
