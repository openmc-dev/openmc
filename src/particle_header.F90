module particle_header

  use, intrinsic :: ISO_C_BINDING

  use bank_header,     only: Bank
  use constants
  use string,          only: to_c_string

  implicit none

!===============================================================================
! LOCALCOORD describes the location of a particle local to a single
! universe. When the geometry consists of nested universes, a particle will have
! a list of coordinates in each level
!===============================================================================

  type, bind(C) :: LocalCoord
    ! Indices in various arrays for this level
    integer(C_INT) :: cell      = NONE
    integer(C_INT) :: universe  = NONE
    integer(C_INT) :: lattice   = NONE
    integer(C_INT) :: lattice_x = NONE
    integer(C_INT) :: lattice_y = NONE
    integer(C_INT) :: lattice_z = NONE

    ! Particle position and direction for this level
    real(C_DOUBLE) :: xyz(3)
    real(C_DOUBLE) :: uvw(3)

    ! Is this level rotated?
    logical(C_BOOL) :: rotated = .false.
  end type LocalCoord

!===============================================================================
! PARTICLE describes the state of a particle being transported through the
! geometry
!===============================================================================

  type, bind(C) :: Particle
    ! Basic data
    integer(C_INT64_T) :: id            ! Unique ID
    integer(C_INT)     :: type          ! Particle type (n, p, e, etc)

    ! Particle coordinates
    integer(C_INT)   :: n_coord          ! number of current coordinates
    integer(C_INT)   :: cell_instance    ! offset for distributed properties
    type(LocalCoord) :: coord(MAX_COORD) ! coordinates for all levels

    ! Particle coordinates before crossing a surface
    integer(C_INT) :: last_n_coord         ! number of current coordinates
    integer(C_INT) :: last_cell(MAX_COORD) ! coordinates for all levels

    ! Energy Data
    real(C_DOUBLE) :: E      ! post-collision energy
    real(C_DOUBLE) :: last_E ! pre-collision energy
    integer(C_INT) :: g      ! post-collision energy group (MG only)
    integer(C_INT) :: last_g ! pre-collision energy group (MG only)

    ! Other physical data
    real(C_DOUBLE) :: wgt           ! particle weight
    real(C_DOUBLE) :: mu            ! angle of scatter
    logical(C_BOOL) :: alive         ! is particle alive?

    ! Pre-collision physical data
    real(C_DOUBLE) :: last_xyz_current(3) ! coordinates of the last collision or
                                          ! reflective/periodic surface crossing
                                          ! for current tallies
    real(C_DOUBLE) :: last_xyz(3)         ! previous coordinates
    real(C_DOUBLE) :: last_uvw(3)         ! previous direction coordinates
    real(C_DOUBLE) :: last_wgt            ! pre-collision particle weight
    real(C_DOUBLE) :: absorb_wgt          ! weight absorbed for survival biasing

    ! What event last took place
    logical(C_BOOL) :: fission       ! did the particle cause implicit fission
    integer(C_INT) :: event         ! scatter, absorption
    integer(C_INT) :: event_nuclide ! index in nuclides array
    integer(C_INT) :: event_MT      ! reaction MT
    integer(C_INT) :: delayed_group ! delayed group

    ! Post-collision physical data
    integer(C_INT) :: n_bank        ! number of fission sites banked
    real(C_DOUBLE) :: wgt_bank      ! weight of fission sites banked
    integer(C_INT) :: n_delayed_bank(MAX_DELAYED_GROUPS) ! number of delayed fission
                                                     ! sites banked

    ! Indices for various arrays
    integer(C_INT) :: surface       ! index for surface particle is on
    integer(C_INT) :: cell_born     ! index for cell particle was born in
    integer(C_INT) :: material      ! index for current material
    integer(C_INT) :: last_material ! index for last material

    ! Temperature of the current cell
    real(C_DOUBLE) :: sqrtkT        ! sqrt(k_Boltzmann * temperature) in eV
    real(C_DOUBLE) :: last_sqrtKT   ! last temperature

    ! Statistical data
    integer(C_INT) :: n_collision   ! # of collisions

    ! Track output
    logical(C_BOOL) :: write_track = .false.

    ! Secondary particles created
    integer(C_INT64_T) :: n_secondary = 0
    type(Bank)         :: secondary_bank(MAX_SECONDARY)
  end type Particle

  interface
    subroutine reset_coord(c) bind(C)
      import LocalCoord
      type(LocalCoord), intent(inout) :: c
    end subroutine reset_coord

    subroutine particle_clear(p) bind(C)
      import Particle
      type(Particle), intent(inout) :: p
    end subroutine particle_clear

    subroutine particle_create_secondary(p, uvw, E, type, run_CE) bind(C)
      import Particle, C_DOUBLE, C_INT, C_BOOL
      type(Particle), intent(inout) :: p
      real(C_DOUBLE), intent(in)    :: uvw(3)
      real(C_DOUBLE), value         :: E
      integer(C_INT), value         :: type
      logical(C_BOOL), value        :: run_CE
    end subroutine particle_create_secondary

    subroutine particle_initialize(p) bind(C)
      import Particle
      type(Particle), intent(inout) :: p
    end subroutine particle_initialize

    subroutine particle_from_source(p, src) bind(C)
      import Particle, Bank, C_DOUBLE
      type(Particle), intent(inout) :: p
      type(Bank),     intent(in)    :: src
    end subroutine particle_from_source

    subroutine particle_mark_as_lost_c(p, message) bind(C, name='particle_mark_as_lost')
      import Particle, C_CHAR
      type(Particle), intent(in) :: p
      character(kind=C_CHAR), intent(in) :: message(*)
    end subroutine particle_mark_as_lost_c

    subroutine particle_write_restart(p) bind(C)
      import Particle
      type(Particle), intent(in) :: p
    end subroutine particle_write_restart
  end interface

contains

  subroutine particle_mark_as_lost(this, message)
    type(Particle), intent(inout) :: this
    character(*)                  :: message

    call particle_mark_as_lost_c(this, to_c_string(message))
  end subroutine particle_mark_as_lost

end module particle_header
