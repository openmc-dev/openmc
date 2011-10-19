program test_mesh_cross

  use mesh, only: surface_crossings
  use mesh_header, only: StructuredMesh
  use particle_header, only: Particle

  implicit none

  integer, parameter   :: n = 3
  integer, parameter   :: n_x = 5
  integer, parameter   :: n_y = 5
  integer, parameter   :: n_z = 5

  integer :: i, j, k
  real(8) :: diff(3)
  type(StructuredMesh), pointer :: m => null()
  type(Particle),       pointer :: p => null()

  ! Set up mesh
  allocate(m)
  m % n_dimension = n
  allocate(m % dimension(n))
  allocate(m % origin(n))
  allocate(m % width(n))

  m % dimension = (/ 5, 5, 5 /)
  m % origin = (/ 0.0, 0.0, 0.0 /)
  m % width = (/ 10.0, 10.0, 10.0 /)

  ! Set up particle
  allocate(p)
  p % last_xyz = (/ 9.3, 0.2, 0.5 /)
  p % xyz = (/ 47.5, 0.2, 0.2 /)
  diff = p % xyz - p % last_xyz
  p % uvw = diff / norm2(diff)

  call surface_crossings(m, p)

end program test_mesh_cross
