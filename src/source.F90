module source

  use bank_header,     only: Bank
  use constants
  use error,           only: fatal_error
  use geometry_header, only: BASE_UNIVERSE
  use global
  use output,          only: write_message
  use particle_header, only: deallocate_coord
  use physics,         only: watt_spectrum
  use random_lcg,      only: prn, set_particle_seed
  use string,          only: to_str

#ifdef MPI
  use mpi
#endif

  implicit none

contains

!===============================================================================
! ALLOCATE_BANKS allocates memory for the fission and source banks
!===============================================================================

  subroutine allocate_banks()

    integer(8) :: bytes      ! size of fission/source bank
    integer    :: alloc_err  ! allocation error code
#ifndef NO_F2008
    type(Bank) :: bank_obj
#endif

    ! Determine maximum amount of particles to simulate on each processor
    maxwork = ceiling(real(n_particles)/n_procs,8)

    ! ID's of first and last source particles
    bank_first = rank*maxwork + 1
    bank_last  = min((rank+1)*maxwork, n_particles)

    ! number of particles for this processor
    work = bank_last - bank_first + 1

    ! Allocate source bank
    allocate(source_bank(maxwork), STAT=alloc_err)

    ! Check for allocation errors 
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

    ! Check for allocation errors 
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

  end subroutine allocate_banks

!===============================================================================
! INITIALIZE_SOURCE initializes particles in the source bank
!===============================================================================

  subroutine initialize_source()

    integer(8) :: i          ! loop index over bank sites
    integer    :: j          ! dummy loop index
    integer(8) :: id         ! particle id
    real(8)    :: r(3)       ! sampled coordinates
    real(8)    :: phi        ! azimuthal angle
    real(8)    :: mu         ! cosine of polar angle
    real(8)    :: E          ! outgoing energy
    real(8)    :: p_min(3)   ! minimum coordinates of source
    real(8)    :: p_max(3)   ! maximum coordinates of source

    message = "Initializing source particles..."
    call write_message(6)

    if (external_source % type == SRC_FILE) then
       ! Read the source from a binary file instead of sampling from some
       ! assumed source distribution

       call read_source_binary()

    else
       ! Generation source sites from specified distribution in user input
       do i = 1, work
          id = bank_first + i - 1
          source_bank(i) % id = id

          ! Set weight to one
          source_bank(i) % wgt = ONE

          ! initialize random number seed
          call set_particle_seed(id)

          ! sample position from external source
          select case (external_source % type)
          case (SRC_BOX)
             p_min = external_source % values(1:3)
             p_max = external_source % values(4:6)
             r = (/ (prn(), j = 1,3) /)
             source_bank(i) % xyz = p_min + r*(p_max - p_min)
          case (SRC_POINT)
             source_bank(i) % xyz = external_source % values
          end select

          ! sample angle
          phi = TWO*PI*prn()
          mu = TWO*prn() - ONE
          source_bank(i) % uvw(1) = mu
          source_bank(i) % uvw(2) = sqrt(ONE - mu*mu) * cos(phi)
          source_bank(i) % uvw(3) = sqrt(ONE - mu*mu) * sin(phi)
          
          ! sample energy from Watt fission energy spectrum for U-235
          do
             E = watt_spectrum(0.988_8, 2.249_8)
             ! resample if energy is >= 20 MeV
             if (E < 20) exit
          end do

          ! set particle energy
          source_bank(i) % E = E
       end do
    end if
 
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
    p % wgt         = src % wgt
    p % last_wgt    = src % wgt
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

    ! Add paricle's starting weight to count for normalizing tallies later
    total_weight = total_weight + src % wgt

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
    p % wgt_bank      = ZERO
    p % n_collision   = 0

    ! remove any original coordinates
    call deallocate_coord(p % coord0)
    
    ! Set up base level coordinates
    allocate(p % coord0)
    p % coord0 % universe = BASE_UNIVERSE
    p % coord             => p % coord0

  end subroutine initialize_particle

!===============================================================================
! WRITE_SOURCE writes out the final source distribution to a binary file that
! can be used as a starting source in a new simulation
!===============================================================================

  subroutine write_source_binary()
    
#ifdef MPI
    integer                  :: fh     ! file handle
    integer(MPI_OFFSET_KIND) :: offset ! offset in memory (0=beginning of file)

    ! ==========================================================================
    ! PARALLEL I/O USING MPI-2 ROUTINES

    ! Open binary source file for reading
    call MPI_FILE_OPEN(MPI_COMM_WORLD, 'source.binary', MPI_MODE_CREATE + &
         MPI_MODE_WRONLY, MPI_INFO_NULL, fh, mpi_err)

    if (master) then
       offset = 0
       call MPI_FILE_WRITE_AT(fh, offset, n_particles, 1, MPI_INTEGER8, &
            MPI_STATUS_IGNORE, mpi_err)
    end if

    ! Set proper offset for source data on this processor
    offset = 8*(1 + rank*maxwork*9)

    ! Write all source sites
    call MPI_FILE_WRITE_AT(fh, offset, source_bank(1), work, MPI_BANK, &
         MPI_STATUS_IGNORE, mpi_err)

    ! Close binary source file
    call MPI_FILE_CLOSE(fh, mpi_err)

#else
    ! ==========================================================================
    ! SERIAL I/O USING FORTRAN INTRINSIC ROUTINES

    ! Open binary source file for writing
    open(UNIT=UNIT_SOURCE, FILE='source.binary', STATUS='replace', &
         ACCESS='stream')

    ! Write the number of particles
    write(UNIT=UNIT_SOURCE) n_particles

    ! Write information from the source bank
    write(UNIT=UNIT_SOURCE) source_bank(1:work)

    ! Close binary source file
    close(UNIT=UNIT_SOURCE)
#endif

  end subroutine write_source_binary

!===============================================================================
! READ_SOURCE_BINARY reads a source distribution from a source.binary file and
! initializes the source bank
!===============================================================================

  subroutine read_source_binary()

    integer(8) :: n_sites  ! number of sites in binary file
#ifdef MPI
    integer    :: fh       ! file handle
    integer(MPI_OFFSET_KIND) :: offset ! offset in memory (0=beginning of file)
#else
    integer    :: i        ! loop over repeating sites
    integer    :: n_repeat ! number of times to repeat a site
#endif

#ifdef MPI
    ! ==========================================================================
    ! PARALLEL I/O USING MPI-2 ROUTINES

    ! Open binary source file for reading
    call MPI_FILE_OPEN(MPI_COMM_WORLD, 'source.binary', MPI_MODE_RDONLY, &
         MPI_INFO_NULL, fh, mpi_err)

    ! Read number of source sites in file
    offset = 0
    call MPI_FILE_READ_AT(fh, offset, n_sites, 1, MPI_INTEGER8, &
         MPI_STATUS_IGNORE, mpi_err)

    if (n_particles > n_sites) then
       message = "No support yet for source files of smaller size than &
            &specified number of particles per generation."
       call fatal_error()
    else
       ! Set proper offset for source data on this processor
       offset = 8*(1 + rank*maxwork*9)

       ! Read all source sites
       call MPI_FILE_READ_AT(fh, offset, source_bank(1), work, MPI_BANK, &
            MPI_STATUS_IGNORE, mpi_err)

       ! Close binary source file
       call MPI_FILE_CLOSE(fh, mpi_err)
    end if

#else
    ! ==========================================================================
    ! SERIAL I/O USING FORTRAN INTRINSIC ROUTINES

    ! Open binary source file for reading
    open(UNIT=UNIT_SOURCE, FILE='source.binary', STATUS='old', &
         ACCESS='stream')

    ! Read number of source sites in file
    read(UNIT=UNIT_SOURCE) n_sites

    if (n_particles > n_sites) then
       ! The size of the source file is smaller than the number of particles we
       ! need. Thus, read all sites and then duplicate sites as necessary.

       read(UNIT=UNIT_SOURCE) source_bank(1:n_sites)

       ! Let's say we have 300 sites and we need to fill in 1000. This do loop
       ! will fill in sites 301 - 900.

       n_repeat = int(n_particles / n_sites)
       do i = 1, n_repeat - 1
          source_bank(i*n_sites + 1:(i+1)*n_sites) = &
               source_bank((i-1)*n_sites + 1:i*n_sites)
       end do

       ! This final statement would fill sites 901 - 1000 in the above example.

       source_bank(n_repeat*n_sites + 1:n_particles) = &
            source_bank(1:n_particles - n_repeat * n_sites)

    else
       ! The size of the source file is bigger than or equal to the number of
       ! particles we need for one generation. Thus, we can just read as many
       ! sites as we need.

       read(UNIT=UNIT_SOURCE) source_bank(1:n_particles)

    end if

    ! Close binary source file
    close(UNIT=UNIT_SOURCE)
#endif

  end subroutine read_source_binary

end module source
