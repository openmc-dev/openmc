module source

  use global
  use mcnp_random, only: rang, RN_init_particle
  use output, only: message
  use physics, only: watt_spectrum

  implicit none

  integer :: source_index

contains

!=====================================================================
! INIT_SOURCE initializes particles in the source bank
!=====================================================================

  subroutine init_source()

    type(Particle), pointer :: p => null()
    integer :: i,j,k
    integer :: maxwork    ! maxinum # of particles per processor
    real(8) :: r(3)       ! sampled coordinates
    real(8) :: phi        ! azimuthal angle
    real(8) :: mu         ! cosine of polar angle
    real(8) :: p_min(3)   ! minimum coordinates of source
    real(8) :: p_max(3)   ! maximum coordinates of source
    character(250) :: msg ! error message

    msg = 'Initializing source particles...'
    call message(msg, 6)

    ! Allocate fission and source banks
    allocate(source_bank(n_particles))
    allocate(fission_bank(3*n_particles))

    ! Determine maximum amount of particles to simulate on each
    ! processor
    maxwork = ceiling(real(n_particles)/n_procs)

    ! Check external source type
    if (external_source%type == SRC_BOX) then
       p_min = external_source%values(1:3)
       p_max = external_source%values(4:6)
    end if

    ! Initialize first cycle source bank
    do i = 0, n_procs - 1
       if (rank == i) then
          ! UID's of first and last source particles
          bank_first = i*maxwork + 1
          bank_last  = min((i+1)*maxwork, n_particles)

          ! number of particles for this processor
          work = bank_last - bank_first + 1

          do j = bank_first, bank_last
             p => source_bank(j - bank_first + 1)

             ! initialize random number seed
             call RN_init_particle(int(j,8))

             ! position
             r = (/ (rang(), k = 1,3) /)
             p % uid = j
             p % xyz = p_min + r*(p_max - p_min)

             ! angle
             phi = TWO*PI*rang()
             mu = TWO*rang() - ONE
             p % uvw(1) = mu
             p % uvw(2) = sqrt(ONE - mu*mu) * cos(phi)
             p % uvw(3) = sqrt(ONE - mu*mu) * sin(phi)

             ! set defaults
             p % type     = NEUTRON
             p % cell     = 0
             p % surface  = 0
             p % universe = 0
             p % wgt      = ONE
             p % alive    = .true.

             ! sample energy from Watt fission energy spectrum for U-235
             p % E = watt_spectrum(0.988_8, 2.249_8)
          end do
       end if
    end do

    ! Reset source index
    source_index = 0

  end subroutine init_source

!=====================================================================
! GET_SOURCE_PARTICLE returns the next source particle 
!=====================================================================

  function get_source_particle() result(p)

    type(Particle), pointer :: p

    ! increment index
    source_index = source_index + 1

    ! if at end of bank, return null pointer
    if (source_index > work) then
       p => null()
       return
    end if

    ! point to next source particle
    p => source_bank(source_index)

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

!=====================================================================
! COPY_FROM_BANK
!=====================================================================

  subroutine copy_from_bank(temp_bank, index, n_sites)

    type(Bank), intent(in) :: temp_bank(n_sites)
    integer,    intent(in) :: index    ! starting index in source_bank
    integer,    intent(in) :: n_sites  ! # of bank sites to copy

    integer :: i ! index in temp_bank
    integer :: j ! index in source_bank
    
    do i = 1, n_sites
       j = index + i - 1
       source_bank(j) % uid  = temp_bank(i) % uid
       source_bank(j) % xyz  = temp_bank(i) % xyz
       source_bank(j) % uvw  = temp_bank(i) % uvw
       source_bank(j) % E    = temp_bank(i) % E

       ! set defaults
       source_bank(j) % type     = NEUTRON
       source_bank(j) % cell     = 0
       source_bank(j) % surface  = 0
       source_bank(j) % universe = 0
       source_bank(j) % wgt      = ONE
       source_bank(j) % alive    = .true.
    end do

  end subroutine copy_from_bank

end module source
