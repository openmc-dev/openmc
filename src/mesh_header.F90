module mesh_header

  implicit none

!===============================================================================
! STRUCTUREDMESH represents a tesslation of n-dimensional Euclidean space by
! congruent squares or cubes
!===============================================================================

  type StructuredMesh
     integer :: id
     integer :: type
     integer :: n_dimension
     integer, allocatable :: dimension(:)
     real(8), allocatable :: origin(:)
     real(8), allocatable :: width(:)
  end type StructuredMesh

end module mesh_header
