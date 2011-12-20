module source

  use bank_header,          only: Bank
  use constants,            only: ONE, MAX_LINE_LEN
  use cross_section_header, only: Nuclide
  use error,                only: fatal_error
  use global
  use output,               only: write_message
  use particle_header,      only: Particle, initialize_particle
  use physics,              only: watt_spectrum
  use random_lcg,           only: prn, set_particle_seed
  use string,               only: int_to_str

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
    integer    :: alloc_err  ! allocation error code
    integer(8) :: bytes      ! size of fission/source bank
    real(8)    :: r(3)       ! sampled coordinates
    real(8)    :: phi        ! azimuthal angle
    real(8)    :: mu         ! cosine of polar angle
    real(8)    :: E          ! outgoing energy
    real(8)    :: p_min(3)   ! minimum coordinates of source
    real(8)    :: p_max(3)   ! maximum coordinates of source
    type(Bank) :: bank_obj

    message = "Initializing source particles..."
    call write_message(6)

    ! Determine maximum amount of particles to simulate on each processor
    maxwork = ceiling(real(n_particles)/n_procs,8)

    ! Allocate source bank
    allocate(source_bank(maxwork), STAT=alloc_err)
    if (alloc_err /= 0) then
#ifndef NO_F2008
       bytes = maxwork * storage_size(bank_obj) / 8
#else
       bytes = maxwork * 64 / 8
#endif
       message = "Could not allocate source bank. Attempted to allocate " &
            // trim(int_to_str(bytes)) // " bytes."
       call fatal_error()
    end if

    ! Allocate fission bank
    allocate(fission_bank(3*maxwork), STAT=alloc_err)
    if (alloc_err /= 0) then
#ifndef NO_F2008
       bytes = 3 * maxwork * storage_size(bank_obj) / 8
#else
       bytes = 3 * maxwork * 64 / 8
#endif
       message = "Could not allocate fission bank. Attempted to allocate " &
            // trim(int_to_str(bytes)) // " bytes."
       call fatal_error()
    end if

    ! Check external source type
    if (external_source%type == SRC_BOX) then
       p_min = external_source%values(1:3)
       p_max = external_source%values(4:6)
    else
       message = "Unsupported external source type: " // &
            int_to_str(external_source%type)
       call fatal_error()
    end if

    ! Initialize first cycle source bank
    do i = 0, n_procs - 1
       if (rank == i) then
          ! ID's of first and last source particles
          bank_first = i*maxwork + 1
          bank_last  = min((i+1)*maxwork, n_particles)

          ! number of particles for this processor
          work = bank_last - bank_first + 1

          do j = bank_first, bank_last
             p => source_bank(j - bank_first + 1)

             ! set defaults
             call initialize_particle(p)

             ! initialize random number seed
             call set_particle_seed(int(j,8))

             ! sample position
             r = (/ (prn(), k = 1,3) /)
             p % id = j
             p % coord0 % xyz = p_min + r*(p_max - p_min)
             p % last_xyz = p % coord0 % xyz

             ! sample angle
             phi = TWO*PI*prn()
             mu = TWO*prn() - ONE
             p % coord0 % uvw(1) = mu
             p % coord0 % uvw(2) = sqrt(ONE - mu*mu) * cos(phi)
             p % coord0 % uvw(3) = sqrt(ONE - mu*mu) * sin(phi)

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

    ! set id
    p % id = bank_first + source_index - 1

  end function get_source_particle

end module source
