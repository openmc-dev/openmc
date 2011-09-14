module particle_header

  implicit none

!===============================================================================
! PARTICLE describes the state of a particle being transported through the
! geometry
!===============================================================================

  type Particle
     ! Basic data
     integer(8) :: uid           ! Unique ID
     integer    :: type          ! Particle type (n, p, e, etc)

     ! Physical data
     real(8)    :: xyz(3)        ! location
     real(8)    :: xyz_local(3)  ! local location (after transformations)
     real(8)    :: uvw(3)        ! directional cosines
     real(8)    :: wgt           ! particle weight
     real(8)    :: E             ! energy
     logical    :: alive         ! is particle alive?

     ! Energy grid data
     integer    :: IE            ! index on energy grid
     real(8)    :: interp        ! interpolation factor for energy grid

     ! Indices for various arrays
     integer    :: cell          ! index for current cell
     integer    :: universe      ! index for current universe
     integer    :: lattice       ! index for current lattice
     integer    :: surface       ! index for current surface
     integer    :: material      ! index for current material
     integer    :: last_material ! index for last material
     integer    :: index_x       ! lattice index for x direction
     integer    :: index_y       ! lattice index for y direction

     ! Statistical data
     integer    :: n_collision   ! # of collisions

  end type Particle

end module particle_header
