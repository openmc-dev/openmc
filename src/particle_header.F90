module particle_header

  use constants, only: NEUTRON, ONE, NONE, ZERO

  implicit none

!===============================================================================
! LOCALCOORD describes the location of a particle local to a single
! universe. When the geometry consists of nested universes, a particle will have
! a list of coordinates in each level
!===============================================================================

  type LocalCoord
     ! Indices in various arrays for this level
     integer :: cell      = NONE
     integer :: universe  = NONE
     integer :: lattice   = NONE
     integer :: lattice_x = NONE
     integer :: lattice_y = NONE

     ! Particle position and direction for this level
     real(8) :: xyz(3)
     real(8) :: uvw(3)

     ! Pointer to next (more local) set of coordinates
     type(LocalCoord), pointer :: next => null()
  end type LocalCoord
     
!===============================================================================
! PARTICLE describes the state of a particle being transported through the
! geometry
!===============================================================================

  type Particle
     ! Basic data
     integer(8) :: id            ! Unique ID
     integer    :: type          ! Particle type (n, p, e, etc)

     ! Particle coordinates
     type(LocalCoord), pointer :: coord0 => null() ! coordinates on universe 0
     type(LocalCoord), pointer :: coord  => null() ! coordinates on lowest universe

     ! Other physical data
     real(8)    :: wgt           ! particle weight
     real(8)    :: E             ! energy
     real(8)    :: mu            ! angle of scatter
     logical    :: alive         ! is particle alive?

     ! Pre-collision physical data
     real(8)    :: last_xyz(3)   ! previous coordinates
     real(8)    :: last_wgt      ! last particle weight
     real(8)    :: last_E        ! last energy

     ! What event last took place
     integer    :: event         ! scatter, absorption, fission
     integer    :: event_nuclide ! index in nuclides array

     ! Post-collision physical data
     integer    :: n_bank        ! number of fission sites banked

     ! Energy grid data
     integer    :: IE            ! index on energy grid
     real(8)    :: interp        ! interpolation factor for energy grid

     ! Indices for various arrays
     integer    :: surface       ! index for surface particle is on
     integer    :: cell_born     ! index for cell particle was born in
     integer    :: material      ! index for current material
     integer    :: last_material ! index for last material

     ! Statistical data
     integer    :: n_collision   ! # of collisions

  end type Particle

contains

!===============================================================================
! DEALLOCATE_COORD removes all levels of coordinates below a given level. This
! is used in distance_to_boundary when the particle moves from a lower universe
! to a higher universe since the data for the lower one is not needed anymore.
!===============================================================================

  recursive subroutine deallocate_coord(coord)

    type(LocalCoord), pointer :: coord

    if (associated(coord)) then 
       ! recursively deallocate lower coordinates
       if (associated(coord % next)) call deallocate_coord(coord%next)

       ! deallocate original coordinate
       deallocate(coord)
    end if

  end subroutine deallocate_coord

end module particle_header
