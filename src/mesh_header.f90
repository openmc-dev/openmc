module mesh_header

  implicit none

!===============================================================================
! STRUCTUREDMESH represents a tesslation of n-dimensional Euclidean space by
! congruent squares or cubes
!===============================================================================

  type StructuredMesh
     integer :: uid
     integer :: type
     integer :: n_dimension
     integer, allocatable :: dimension(:)
     integer, allocatable :: origin(:)
     integer, allocatable :: width(:)
  end type StructuredMesh

end module mesh_header
