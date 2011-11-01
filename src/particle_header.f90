module particle_header

  use constants, only: NEUTRON, ONE

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
     real(8)    :: mu            ! angle of scatter
     logical    :: alive         ! is particle alive?

     ! Pre-collision physical data
     real(8)    :: last_xyz(3)   ! previous coordinates
     real(8)    :: last_wgt      ! last particle weight
     real(8)    :: last_E        ! last energy

     ! Post-collision physical data
     integer    :: n_bank        ! number of fission sites banked

     ! Energy grid data
     integer    :: IE            ! index on energy grid
     real(8)    :: interp        ! interpolation factor for energy grid

     ! Indices for various arrays
     integer    :: cell          ! index for current cell
     integer    :: cell_born     ! index for cell particle was born in
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

contains

!===============================================================================
! INITIALIZE_PARTICLE sets default attributes for a particle from the source
! bank
!===============================================================================

  subroutine initialize_particle(p)

    type(Particle), pointer :: p

    ! TODO: if information on the cell, lattice, universe, and material is
    ! passed through the fission bank to the source bank, no lookup would be
    ! needed at the beginning of a cycle

    p % type          = NEUTRON
    p % wgt           = ONE
    p % last_wgt      = ONE
    p % alive         = .true.
    p % n_bank        = 0
    p % cell          = 0
    p % cell_born     = 0
    p % universe      = 0
    p % lattice       = 0
    p % surface       = 0
    p % material      = 0
    p % last_material = 0
    p % index_x       = 0
    p % index_y       = 0
    p % n_collision   = 0

  end subroutine initialize_particle

end module particle_header
