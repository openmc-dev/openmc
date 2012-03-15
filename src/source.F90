module source

  use bank_header,     only: Bank
  use constants,       only: ONE
  use error,           only: fatal_error
  use geometry_header, only: BASE_UNIVERSE
  use global
  use output,          only: write_message
  use particle_header, only: deallocate_coord
  use physics,         only: watt_spectrum
  use random_lcg,      only: prn, set_particle_seed
  use string,          only: to_str

  implicit none

contains

!===============================================================================
! INITIALIZE_SOURCE initializes particles in the source bank
!===============================================================================

  subroutine initialize_source()

    integer    :: i          ! loop index over processors
    integer(8) :: j          ! loop index over bank sites
    integer    :: k          ! dummy loop index
    integer(8) :: bytes      ! size of fission/source bank
    integer(8) :: id         ! particle id
    integer    :: alloc_err  ! allocation error code
    real(8)    :: r(3)       ! sampled coordinates
    real(8)    :: phi        ! azimuthal angle
    real(8)    :: mu         ! cosine of polar angle
    real(8)    :: E          ! outgoing energy
    real(8)    :: p_min(3)   ! minimum coordinates of source
    real(8)    :: p_max(3)   ! maximum coordinates of source
#ifndef NO_F2008
    type(Bank) :: bank_obj
#endif

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
            // trim(to_str(bytes)) // " bytes."
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
            // trim(to_str(bytes)) // " bytes."
       call fatal_error()
    end if

    ! Initialize first cycle source bank
    do i = 0, n_procs - 1
       if (rank == i) then
          ! ID's of first and last source particles
          bank_first = rank*maxwork + 1
          bank_last  = min((rank+1)*maxwork, n_particles)

          ! number of particles for this processor
          work = bank_last - bank_first + 1

          do j = 1, work
             id = bank_first + j - 1
             source_bank(j) % id = id

             ! initialize random number seed
             call set_particle_seed(id)

             ! sample position from external source
             select case (external_source % type)
             case (SRC_BOX)
                p_min = external_source % values(1:3)
                p_max = external_source % values(4:6)
                r = (/ (prn(), k = 1,3) /)
                source_bank(j) % xyz = p_min + r*(p_max - p_min)
             end select

             ! sample angle
             phi = TWO*PI*prn()
             mu = TWO*prn() - ONE
             source_bank(j) % uvw(1) = mu
             source_bank(j) % uvw(2) = sqrt(ONE - mu*mu) * cos(phi)
             source_bank(j) % uvw(3) = sqrt(ONE - mu*mu) * sin(phi)

             ! sample energy from Watt fission energy spectrum for U-235
             do
                E = watt_spectrum(0.988_8, 2.249_8)
                ! resample if energy is >= 20 MeV
                if (E < 20) exit
             end do

             ! set particle energy
             source_bank(j) % E = E
          end do
       end if
    end do

  end subroutine initialize_source

!===============================================================================
! GET_SOURCE_PARTICLE returns the next source particle 
!===============================================================================

  subroutine get_source_particle(index_source)

    integer(8), intent(in) :: index_source

    integer(8) :: particle_seed  ! unique index for particle
    type(Bank), pointer :: src => null()

    ! set defaults
    call initialize_particle()

    ! point to next source particle
    src => source_bank(index_source)

    ! copy attributes from source bank site
    p % coord % xyz = src % xyz
    p % coord % uvw = src % uvw
    p % last_xyz    = src % xyz
    p % E           = src % E
    p % last_E      = src % E

    ! set identifier for particle
    p % id = bank_first + index_source - 1

    ! set random number seed
    particle_seed = ((current_batch - 1)*gen_per_batch + & 
         current_gen - 1)*n_particles + p % id
    call set_particle_seed(particle_seed)
          
    ! set particle trace
    trace = .false.
    if (current_batch == trace_batch .and. current_gen == trace_gen .and. &
         p % id == trace_particle) trace = .true.

  end subroutine get_source_particle

!===============================================================================
! INITIALIZE_PARTICLE sets default attributes for a particle from the source
! bank
!===============================================================================

  subroutine initialize_particle()

    ! Set particle to neutron that's alive
    p % type  = NEUTRON
    p % alive = .true.

    ! clear attributes
    p % surface       = NONE
    p % cell_born     = NONE
    p % material      = NONE
    p % last_material = NONE
    p % wgt           = ONE
    p % last_wgt      = ONE
    p % n_bank        = 0
    p % n_collision   = 0

    ! remove any original coordinates
    call deallocate_coord(p % coord0)
    
    ! Set up base level coordinates
    allocate(p % coord0)
    p % coord0 % universe = BASE_UNIVERSE
    p % coord             => p % coord0

  end subroutine initialize_particle

end module source
