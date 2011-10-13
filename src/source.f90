module source

  use bank_header,          only: Bank
  use constants,            only: ONE, MAX_LINE_LEN
  use cross_section_header, only: Nuclide
  use global
  use mcnp_random,          only: rang, RN_init_particle
  use output,               only: message
  use particle_header,      only: Particle, initialize_particle
  use physics,              only: watt_spectrum

  implicit none

contains

!===============================================================================
! INITIALIZE_SOURCE initializes particles in the source bank
!===============================================================================

  subroutine initialize_source()

    type(Particle), pointer :: p => null()
    integer    :: i          ! loop index over processors
    integer(8) :: j          ! loop index over bank sites
    integer    :: k          ! dummy loop index
    integer(8) :: maxwork    ! maxinum # of particles per processor
    real(8)    :: r(3)       ! sampled coordinates
    real(8)    :: phi        ! azimuthal angle
    real(8)    :: mu         ! cosine of polar angle
    real(8)    :: E          ! outgoing energy
    real(8)    :: p_min(3)   ! minimum coordinates of source
    real(8)    :: p_max(3)   ! maximum coordinates of source
    character(MAX_LINE_LEN) :: msg    ! error message

    msg = 'Initializing source particles...'
    call message(msg, 6)

    ! Determine maximum amount of particles to simulate on each processor
    maxwork = ceiling(real(n_particles)/n_procs,8)

    ! Allocate fission and source banks
    allocate(source_bank(maxwork))
    allocate(fission_bank(3*maxwork))

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

             ! sample position
             r = (/ (rang(), k = 1,3) /)
             p % uid = j
             p % xyz = p_min + r*(p_max - p_min)
             p % xyz_local = p % xyz

             ! sample angle
             phi = TWO*PI*rang()
             mu = TWO*rang() - ONE
             p % uvw(1) = mu
             p % uvw(2) = sqrt(ONE - mu*mu) * cos(phi)
             p % uvw(3) = sqrt(ONE - mu*mu) * sin(phi)

             ! set defaults
             call initialize_particle(p)

             ! sample energy from Watt fission energy spectrum for U-235
             do
                E = watt_spectrum(0.988_8, 2.249_8)
                ! resample if energy is >= 20 MeV
                if (E < 20) exit
             end do

             ! set particle energy
             p % E = E
             p % last_E = E
          end do
       end if
    end do

    ! Reset source index
    source_index = 0_8

  end subroutine initialize_source

!===============================================================================
! GET_SOURCE_PARTICLE returns the next source particle 
!===============================================================================

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

    ! set uid
    p % uid = bank_first + source_index - 1

  end function get_source_particle

end module source
