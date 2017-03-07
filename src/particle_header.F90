module particle_header

  use constants,       only: NEUTRON, ONE, NONE, ZERO
  use geometry_header, only: BASE_UNIVERSE

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
    integer :: lattice_z = NONE

    ! Particle position and direction for this level
    real(8) :: xyz(3)
    real(8) :: uvw(3)

    ! Is this level rotated?
    logical :: rotated = .false.

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
    real(8)    :: last_uvw(3)   ! previous direction coordinates
    real(8)    :: last_wgt      ! pre-collision particle weight
    real(8)    :: last_E        ! pre-collision energy
    real(8)    :: absorb_wgt    ! weight absorbed for survival biasing

    ! What event last took place
    logical    :: fission       ! did the particle cause implicit fission
    integer    :: event         ! scatter, absorption
    integer    :: event_nuclide ! index in nuclides array
    integer    :: event_MT      ! reaction MT

    ! Post-collision physical data
    integer    :: n_bank        ! number of fission sites banked
    real(8)    :: wgt_bank      ! weight of fission sites banked

    ! Indices for various arrays
    integer    :: surface       ! index for surface particle is on
    integer    :: cell_born     ! index for cell particle was born in
    integer    :: material      ! index for current material
    integer    :: last_material ! index for last material

    ! Statistical data
    integer    :: n_collision   ! # of collisions

    ! Track output
    logical    :: write_track = .false.

  contains
    procedure :: initialize => initialize_particle
    procedure :: clear => clear_particle
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

!===============================================================================
! INITIALIZE_PARTICLE sets default attributes for a particle from the source
! bank
!===============================================================================

  subroutine initialize_particle(this)

    class(Particle) :: this

    ! Clear coordinate lists
    call this % clear()

    ! Set particle to neutron that's alive
    this % type  = NEUTRON
    this % alive = .true.

    ! clear attributes
    this % surface       = NONE
    this % cell_born     = NONE
    this % material      = NONE
    this % last_material = NONE
    this % wgt           = ONE
    this % last_wgt      = ONE
    this % absorb_wgt    = ZERO
    this % n_bank        = 0
    this % wgt_bank      = ZERO
    this % n_collision   = 0
    this % fission       = .false.

    ! Set up base level coordinates
    allocate(this % coord0)
    this % coord0 % universe = BASE_UNIVERSE
    this % coord             => this % coord0

  end subroutine initialize_particle

!===============================================================================
! CLEAR_PARTICLE
!===============================================================================

  subroutine clear_particle(this)

    class(Particle) :: this

    ! remove any coordinate levels
    call deallocate_coord(this % coord0)

    ! Make sure coord pointer is nullified
    nullify(this % coord)

  end subroutine clear_particle

end module particle_header
