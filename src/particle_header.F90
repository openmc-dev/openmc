module particle_header

  use bank_header,     only: Bank
  use constants,       only: NEUTRON, ONE, NONE, ZERO, MAX_SECONDARY, &
                             MAX_DELAYED_GROUPS
! asking for help: circular error, compiling failed
!XX  use error,           only: fatal_error
  use geometry_header, only: BASE_UNIVERSE
  use random_lcg_header, only: N_STREAMS

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
  contains
    procedure :: reset => reset_coord
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
    integer          :: n_coord          ! number of current coordinates
    type(LocalCoord) :: coord(MAX_COORD) ! coordinates for all levels

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
    integer    :: delayed_group ! delayed group

    ! Post-collision physical data
    integer    :: n_bank        ! number of fission sites banked
    real(8)    :: wgt_bank      ! weight of fission sites banked
    integer    :: n_delayed_bank(MAX_DELAYED_GROUPS) ! number of delayed fission
                                                     ! sites banked

    ! Indices for various arrays
    integer    :: surface       ! index for surface particle is on
    integer    :: cell_born     ! index for cell particle was born in
    integer    :: material      ! index for current material
    integer    :: last_material ! index for last material

    ! Statistical data
    integer    :: n_collision   ! # of collisions

    ! Track output
    logical    :: write_track = .false.

    ! Secondary particles created
    integer    :: n_secondary = 0
    type(Bank) :: secondary_bank(MAX_SECONDARY)
    ! Was this particle just created?
    logical    :: new_particle = .true.

    ! Data needed to restart a particle stored in a bank after changing domains
    real(8)    :: stored_xyz(3)
    real(8)    :: stored_uvw(3)
    real(8)    :: stored_distance ! sampled distance to go after changing domain
    real(8)    :: fly_dd_distance ! accumulated distance to domain boundary
    integer(8) :: prn_seed(N_STREAMS) ! the  next random number seed
    integer(8) :: xs_seed(N_STREAMS)  ! the previously seed used for xs gen

    ! Domain information
    integer    :: dd_meshbin    ! DD meshbin the particle is to be run in next
    integer    :: outscatter_destination ! Which domain to transmit particle to

  contains
    procedure :: initialize => initialize_particle
    procedure :: clear => clear_particle
    procedure :: initialize_from_source
    procedure :: create_secondary
  end type Particle

!===============================================================================
! PARTICLEBUFFER is a copy of the particle data structure for sending to other
! processes via MPI - it's needed to get proper alignment of fields, and it
! provides a way to cut down on some of the info we need to send
!===============================================================================

  type ParticleBuffer

    sequence

    integer(8) :: id
    integer(8) :: prn_seed(N_STREAMS)
    integer(8) :: xs_seed(N_STREAMS)
    real(8)    :: wgt
    real(8)    :: E
    real(8)    :: mu
    real(8)    :: last_wgt
    real(8)    :: last_E
    real(8)    :: absorb_wgt
    real(8)    :: wgt_bank
    real(8)    :: stored_distance
    real(8)    :: fly_dd_distance
    real(8)    :: last_xyz(3)
    real(8)    :: stored_xyz(3)
    real(8)    :: stored_uvw(3)
    integer    :: type
    integer    :: event
    integer    :: event_nuclide
    integer    :: event_MT
    integer    :: n_bank
    integer    :: surface
    integer    :: cell_born
    integer    :: n_collision
    integer    :: material
    integer    :: last_material

  end type ParticleBuffer

contains

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
    this % surface           = NONE
    this % cell_born         = NONE
    this % material          = NONE
    this % last_material     = NONE
    this % wgt               = ONE
    this % last_wgt          = ONE
    this % absorb_wgt        = ZERO
    this % n_bank            = 0
    this % wgt_bank          = ZERO
    this % n_collision       = 0
    this % fission           = .false.
    this % delayed_group     = 0
    this % n_delayed_bank(:) = 0
    this % write_track       = .false.
    this % new_particle      = .true.
    this % stored_distance   = ZERO
    this % fly_dd_distance   = ZERO

    ! Set up base level coordinates
    this % coord(1) % universe = BASE_UNIVERSE
    this % n_coord = 1

  end subroutine initialize_particle

!===============================================================================
! PARTICLE_TO_BUFFER copies particle information to the ParticleBuffer data type
! for sending via MPI
!===============================================================================

  subroutine particle_to_buffer(part, buf)

    type(Particle),       intent(in)  :: part
    type(ParticleBuffer), intent(out) :: buf

    buf % id              = part % id
    buf % type            = part % type
    buf % wgt             = part % wgt
    buf % E               = part % E
    buf % mu              = part % mu
    buf % last_xyz        = part % last_xyz
    buf % last_wgt        = part % last_wgt
    buf % last_E          = part % last_E
    buf % absorb_wgt      = part % absorb_wgt
    buf % event           = part % event
    buf % event_nuclide   = part % event_nuclide
    buf % event_MT        = part % event_MT
    buf % n_bank          = part % n_bank
    buf % wgt_bank        = part % wgt_bank
    buf % surface         = part % surface
    buf % cell_born       = part % cell_born
    buf % n_collision     = part % n_collision
    buf % stored_xyz      = part % stored_xyz
    buf % stored_uvw      = part % stored_uvw
    buf % material        = part % material
    buf % last_material   = part % last_material
    buf % prn_seed        = part % prn_seed
    buf % xs_seed         = part % xs_seed
    buf % stored_distance = part % stored_distance
    buf % fly_dd_distance = part % fly_dd_distance

  end subroutine particle_to_buffer

!===============================================================================
! BUFFER_TO_PARTICLE copies particle information out of ParticleBuffer type
! received via MPI back into a particle data type
!===============================================================================

  subroutine buffer_to_particle(buf, part)

    type(ParticleBuffer), intent(in)  :: buf
    type(Particle),       intent(out) :: part

    part % id              = buf % id
    part % type            = buf % type
    part % wgt             = buf % wgt
    part % E               = buf % E
    part % mu              = buf % mu
    part % alive           = .true.
    part % last_xyz        = buf % last_xyz
    part % last_wgt        = buf % last_wgt
    part % last_E          = buf % last_E
    part % absorb_wgt      = buf % absorb_wgt
    part % event           = buf % event
    part % event_nuclide   = buf % event_nuclide
    part % event_MT        = buf % event_MT
    part % n_bank          = buf % n_bank
    part % wgt_bank        = buf % wgt_bank
    part % surface         = buf % surface
    part % cell_born       = buf % cell_born
    part % n_collision     = buf % n_collision
    part % new_particle    = .false.
    part % stored_xyz      = buf % stored_xyz
    part % stored_uvw      = buf % stored_uvw
    part % material        = buf % material
    part % last_material   = buf % last_material
    part % prn_seed        = buf % prn_seed
    part % xs_seed         = buf % xs_seed
    part % stored_distance = buf % stored_distance
    part % fly_dd_distance = buf % fly_dd_distance

  end subroutine buffer_to_particle

!===============================================================================
! CLEAR_PARTICLE resets all coordinate levels for the particle
!===============================================================================

  subroutine clear_particle(this)

    class(Particle) :: this
    integer :: i

    ! remove any coordinate levels
    do i = 1, MAX_COORD
      call this % coord(i) % reset()
    end do

  end subroutine clear_particle

!===============================================================================
! RESET_COORD clears data from a single coordinate level
!===============================================================================

  elemental subroutine reset_coord(this)
    class(LocalCoord), intent(inout) :: this

    this % cell = NONE
    this % universe = NONE
    this % lattice = NONE
    this % lattice_x = NONE
    this % lattice_y = NONE
    this % lattice_z = NONE
    this % rotated = .false.

  end subroutine reset_coord

!===============================================================================
! INITIALIZE_FROM_SOURCE initializes a particle from data stored in a source
! site. The source site may have been produced from an external source, from
! fission, or simply as a secondary particle.
!===============================================================================

  subroutine initialize_from_source(this, src)
    class(Particle), intent(inout) :: this
    type(Bank),      intent(in)    :: src

    ! set defaults
    call this % initialize()

    ! copy attributes from source bank site
    this % wgt            = src % wgt
    this % last_wgt       = src % wgt
    this % coord(1) % xyz = src % xyz
    this % coord(1) % uvw = src % uvw
    this % last_xyz       = src % xyz
    this % last_uvw       = src % uvw
    this % E              = src % E
    this % last_E         = src % E
    this % prn_seed       = src % prn_seed
    this % xs_seed        = src % prn_seed

  end subroutine initialize_from_source

!===============================================================================
! CREATE_SECONDARY stores the current phase space attributes of the particle in
! the secondary bank and increments the number of sites in the secondary bank.
!===============================================================================

  subroutine create_secondary(this, uvw, type)
    class(Particle), intent(inout) :: this
    real(8),         intent(in)    :: uvw(3)
    integer,         intent(in)    :: type

    integer :: n

    ! Check to make sure that the hard-limit on secondary particles is not
    ! exceeded.
!asking for help: circular error, compiling failed
!XX    if (this % n_secondary == MAX_SECONDARY) then
!XX      call fatal_error("Too many secondary particles created.")
!XX    end if

    n = this % n_secondary + 1
    this % secondary_bank(n) % wgt    = this % wgt
    this % secondary_bank(n) % xyz(:) = this % coord(1) % xyz
    this % secondary_bank(n) % uvw(:) = uvw
    this % secondary_bank(n) % E      = this % E
    this % n_secondary = n

  end subroutine create_secondary

end module particle_header
