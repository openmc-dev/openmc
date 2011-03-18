module source

  use global
  use mcnp_random, only: rang, RN_init_particle

  implicit none

  integer :: source_index

contains

!=====================================================================
! INIT_SOURCE initializes particles in the source bank
!=====================================================================

  subroutine init_source()

    integer :: i,j
    real(8) :: r(3)
    real(8) :: phi ! azimuthal angle
    real(8) :: mu  ! cosine of polar angle
    real(8) :: p_min(3), p_max(3)
    type(Particle), pointer :: p => null()

    ! Allocate fission and source banks
    allocate( source_bank(n_particles)    )
    allocate( fission_bank(3*n_particles) )

    ! Check external source type
    if ( external_source%type == SRC_BOX ) then
       p_min = external_source%values(1:3)
       p_max = external_source%values(4:6)
    end if

    ! Initialize first cycle source bank
    do i = 1, n_particles
       call RN_init_particle(int(i,8))
       p => source_bank(i)

       ! position
       r = (/ (rang(), j = 1,3) /)
       p % uid = i
       p % xyz = p_min + r*(p_max - p_min)

       ! angle
       phi = TWO*PI*rang()
       mu = TWO*rang() - ONE
       p % uvw(1) = mu
       p % uvw(2) = sqrt(ONE - mu*mu) * cos(phi)
       p % uvw(3) = sqrt(ONE - mu*mu) * sin(phi)

       ! set defaults
       p % type = NEUTRON
       p % cell = 0
       p % surface = 0
       p % universe = 0
       p % wgt = ONE
       p % alive = .true.

    end do

    ! Reset source index
    source_index = 0

  end subroutine init_source

!=====================================================================
! GET_SOURCE_PARTICLE returns the next source particle 
!=====================================================================

  function get_source_particle()

    type(Particle), pointer :: get_source_particle

    ! increment index
    source_index = source_index + 1

    ! if at end of bank, return null pointer
    if (source_index > size(source_bank)) then
       get_source_particle => null()
       return
    end if

    ! point to next source particle
    get_source_particle => source_bank(source_index)

  end function get_source_particle

!=====================================================================
! ADD_BANK_SITES
!=====================================================================

  subroutine add_bank_sites(p, table, n)

    type(Particle),      pointer    :: p
    type(AceContinuous), pointer    :: table
    integer,             intent(in) :: n

    integer :: i

    if (n == 0 .or. n_bank == 3*n_particles) return
    do i = n_bank + 1, min(n_bank + n, 3*n_particles)
      ! Copy particle data
      fission_bank(i)%uid = p%uid
      fission_bank(i)%xyz = p%xyz

      ! TODO: Sample angle and energy from secondary distribution
    end do
    n_bank = min(n_bank + n, 3*n_particles)

  end subroutine add_bank_sites

end module source
