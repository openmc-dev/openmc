module particle_header

  use hdf5, only: HID_T

  use bank_header,     only: Bank, source_bank
  use constants
  use error,           only: fatal_error, warning
  use geometry_header, only: root_universe
  use hdf5_interface
  use settings
  use simulation_header
  use string,          only: to_str

  implicit none

  private

!===============================================================================
! LOCALCOORD describes the location of a particle local to a single
! universe. When the geometry consists of nested universes, a particle will have
! a list of coordinates in each level
!===============================================================================

  type, public :: LocalCoord

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

  type, public :: Particle
    ! Basic data
    integer(8) :: id            ! Unique ID
    integer    :: type          ! Particle type (n, p, e, etc)

    ! Particle coordinates
    integer          :: n_coord          ! number of current coordinates
    integer          :: cell_instance    ! offset for distributed properties
    type(LocalCoord) :: coord(MAX_COORD) ! coordinates for all levels

    ! Particle coordinates before crossing a surface
    integer :: last_n_coord         ! number of current coordinates
    integer :: last_cell(MAX_COORD) ! coordinates for all levels

    ! Energy Data
    real(8)    :: E      ! post-collision energy
    real(8)    :: last_E ! pre-collision energy
    integer    :: g      ! post-collision energy group (MG only)
    integer    :: last_g ! pre-collision energy group (MG only)

    ! Other physical data
    real(8)    :: wgt           ! particle weight
    real(8)    :: mu            ! angle of scatter
    logical    :: alive         ! is particle alive?

    ! Pre-collision physical data
    real(8)    :: last_xyz_current(3) ! coordinates of the last collision or
                                      !  reflective/periodic surface crossing
                                      !  for current tallies
    real(8)    :: last_xyz(3)         ! previous coordinates
    real(8)    :: last_uvw(3)         ! previous direction coordinates
    real(8)    :: last_wgt            ! pre-collision particle weight
    real(8)    :: absorb_wgt          ! weight absorbed for survival biasing

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

    ! Temperature of the current cell
    real(8)    :: sqrtkT        ! sqrt(k_Boltzmann * temperature) in eV
    real(8)    :: last_sqrtKT   ! last temperature

    ! Statistical data
    integer    :: n_collision   ! # of collisions

    ! Track output
    logical    :: write_track = .false.

    ! Secondary particles created
    integer(8) :: n_secondary = 0
    type(Bank) :: secondary_bank(MAX_SECONDARY)

  contains
    procedure :: clear
    procedure :: create_secondary
    procedure :: initialize
    procedure :: initialize_from_source
    procedure :: mark_as_lost
    procedure :: write_restart
  end type Particle

contains

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
! CLEAR_PARTICLE resets all coordinate levels for the particle
!===============================================================================

  subroutine clear(this)
    class(Particle) :: this

    integer :: i

    ! remove any coordinate levels
    do i = 1, MAX_COORD
      call this % coord(i) % reset()
    end do
  end subroutine clear

!===============================================================================
! CREATE_SECONDARY stores the current phase space attributes of the particle in
! the secondary bank and increments the number of sites in the secondary bank.
!===============================================================================

  subroutine create_secondary(this, uvw, type, run_CE)
    class(Particle), intent(inout) :: this
    real(8),         intent(in)    :: uvw(3)
    integer,         intent(in)    :: type
    logical,         intent(in)    :: run_CE

    integer(8) :: n

    ! Check to make sure that the hard-limit on secondary particles is not
    ! exceeded.
    if (this % n_secondary == MAX_SECONDARY) then
      call fatal_error("Too many secondary particles created.")
    end if

    n = this % n_secondary + 1
    this % secondary_bank(n) % wgt    = this % wgt
    this % secondary_bank(n) % xyz(:) = this % coord(1) % xyz
    this % secondary_bank(n) % uvw(:) = uvw
    this % n_secondary = n
    this % secondary_bank(this % n_secondary) % E = this % E
    if (.not. run_CE) then
      this % secondary_bank(this % n_secondary) % E = real(this % g, 8)
    end if

  end subroutine create_secondary

!===============================================================================
! INITIALIZE sets default attributes for a particle from the source bank
!===============================================================================

  subroutine initialize(this)

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
    this % last_sqrtkT       = NONE
    this % wgt               = ONE
    this % last_wgt          = ONE
    this % absorb_wgt        = ZERO
    this % n_bank            = 0
    this % wgt_bank          = ZERO
    this % sqrtkT            = ERROR_REAL
    this % n_collision       = 0
    this % fission           = .false.
    this % delayed_group     = 0
    this % n_delayed_bank(:) = 0
    this % g = NONE

    ! Set up base level coordinates
    this % coord(1) % universe = root_universe
    this % n_coord = 1
    this % last_n_coord = 1

  end subroutine initialize

!===============================================================================
! INITIALIZE_FROM_SOURCE initializes a particle from data stored in a source
! site. The source site may have been produced from an external source, from
! fission, or simply as a secondary particle.
!===============================================================================

  subroutine initialize_from_source(this, src, run_CE, energy_bin_avg)
    class(Particle), intent(inout)   :: this
    type(Bank),      intent(in)      :: src
    logical,         intent(in)      :: run_CE
    real(8), allocatable, intent(in) :: energy_bin_avg(:)

    ! set defaults
    call this % initialize()

    ! copy attributes from source bank site
    this % wgt              = src % wgt
    this % last_wgt         = src % wgt
    this % coord(1) % xyz   = src % xyz
    this % coord(1) % uvw   = src % uvw
    this % last_xyz_current = src % xyz
    this % last_xyz         = src % xyz
    this % last_uvw         = src % uvw
    if (run_CE) then
      this % E              = src % E
      this % g              = NONE
    else
      this % g              = int(src % E)
      this % last_g         = int(src % E)
      this % E              = energy_bin_avg(this % g)
    end if
    this % last_E           = this % E

  end subroutine initialize_from_source

!===============================================================================
! MARK_AS_LOST
!===============================================================================

  subroutine mark_as_lost(this, message)
    class(Particle), intent(inout) :: this
    character(*)                   :: message

    integer(8) :: tot_n_particles

    ! Print warning and write lost particle file
    call warning(message)
    call this % write_restart()

    ! Increment number of lost particles
    this % alive = .false.
!$omp atomic
    n_lost_particles = n_lost_particles + 1

    ! Count the total number of simulated particles (on this processor)
    tot_n_particles = current_batch * gen_per_batch * work

    ! Abort the simulation if the maximum number of lost particles has been
    ! reached
    if (n_lost_particles >= MAX_LOST_PARTICLES .and. &
         n_lost_particles >= REL_MAX_LOST_PARTICLES * tot_n_particles) then
      call fatal_error("Maximum number of lost particles has been reached.")
    end if

  end subroutine mark_as_lost

!===============================================================================
! WRITE_RESTART creates a particle restart file
!===============================================================================

  subroutine write_restart(this)
    class(Particle), intent(in) :: this

    integer(HID_T) :: file_id
    character(MAX_FILE_LEN) :: filename

    ! Dont write another restart file if in particle restart mode
    if (run_mode == MODE_PARTICLE) return

    ! Set up file name
    filename = trim(path_output) // 'particle_' // trim(to_str(current_batch)) &
         // '_' // trim(to_str(this % id)) // '.h5'

!$omp critical (WriteParticleRestart)
    ! Create file
    file_id = file_create(filename)

    associate (src => source_bank(current_work))
      ! Write filetype and version info
      call write_attribute(file_id, 'filetype', 'particle restart')
      call write_attribute(file_id, 'version', VERSION_PARTICLE_RESTART)
      call write_attribute(file_id, "openmc_version", VERSION)
#ifdef GIT_SHA1
      call write_attribute(file_id, "git_sha1", GIT_SHA1)
#endif

      ! Write data to file
      call write_dataset(file_id, 'current_batch', current_batch)
      call write_dataset(file_id, 'generations_per_batch', gen_per_batch)
      call write_dataset(file_id, 'current_generation', current_gen)
      call write_dataset(file_id, 'n_particles', n_particles)
      select case(run_mode)
      case (MODE_FIXEDSOURCE)
        call write_dataset(file_id, 'run_mode', 'fixed source')
      case (MODE_EIGENVALUE)
        call write_dataset(file_id, 'run_mode', 'eigenvalue')
      case (MODE_PARTICLE)
        call write_dataset(file_id, 'run_mode', 'particle restart')
      end select
      call write_dataset(file_id, 'id', this % id)
      call write_dataset(file_id, 'weight', src % wgt)
      call write_dataset(file_id, 'energy', src % E)
      call write_dataset(file_id, 'xyz', src % xyz)
      call write_dataset(file_id, 'uvw', src % uvw)
    end associate

    ! Close file
    call file_close(file_id)
!$omp end critical (WriteParticleRestart)

  end subroutine write_restart

end module particle_header
