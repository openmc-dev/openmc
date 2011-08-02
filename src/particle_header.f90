module particle_header

  implicit none

!===============================================================================
! PARTICLE describes the state of a particle being transported through the
! geometry
!===============================================================================

  type Particle
    integer(8) :: uid          ! Unique ID
    integer    :: type         ! Particle type (n, p, e, etc)
    real(8)    :: xyz(3)       ! location
    real(8)    :: xyz_local(3) ! local location (after transformations)
    real(8)    :: uvw(3)       ! directional cosines
    real(8)    :: wgt          ! particle weight
    real(8)    :: E            ! energy
    integer    :: IE           ! index on energy grid
    real(8)    :: interp       ! interpolation factor for energy grid
    integer    :: cell         ! current cell
    integer    :: universe     ! current universe
    integer    :: lattice      ! current lattice
    integer    :: surface      ! current surface
    integer    :: index_x      ! lattice index for x direction
    integer    :: index_y      ! lattice index for y direction
    logical    :: alive        ! is particle alive?
    integer    :: n_coll       ! # of collisions
  end type Particle

end module particle_header
