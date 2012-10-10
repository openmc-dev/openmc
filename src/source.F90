module source

  use bank_header,     only: Bank
  use constants
  use error,           only: fatal_error
  use geometry_header, only: BASE_UNIVERSE
  use global
  use output,          only: write_message
  use particle_header, only: deallocate_coord
  use physics,         only: maxwell_spectrum, watt_spectrum
  use random_lcg,      only: prn, set_particle_seed
  use string,          only: to_str

#ifdef MPI
  use mpi
#endif

  implicit none

contains

!===============================================================================
! INITIALIZE_SOURCE initializes particles in the source bank
!===============================================================================

  subroutine initialize_source()

    integer(8) :: i          ! loop index over bank sites
    integer(8) :: id         ! particle id

    type(Bank), pointer :: src => null() ! source bank site

    message = "Initializing source particles..."
    call write_message(6)

    if (path_source /= '') then
       ! Read the source from a binary file instead of sampling from some
       ! assumed source distribution

       call read_source_binary()

    else
       ! Generation source sites from specified distribution in user input
       do i = 1, work
          ! Get pointer to source bank site
          src => source_bank(i)

          ! initialize random number seed
          id = bank_first + i - 1
          call set_particle_seed(id)

          ! sample external source distribution
          call sample_external_source(src)
       end do
    end if
 
  end subroutine initialize_source

!===============================================================================
! SAMPLE_EXTERNAL_SOURCE
!===============================================================================

  subroutine sample_external_source(site)

    type(Bank), pointer :: site ! source site

    integer :: i          ! dummy loop index
    real(8) :: r(3)       ! sampled coordinates
    real(8) :: phi        ! azimuthal angle
    real(8) :: mu         ! cosine of polar angle
    real(8) :: p_min(3)   ! minimum coordinates of source
    real(8) :: p_max(3)   ! maximum coordinates of source
    real(8) :: a          ! Arbitrary parameter 'a'
    real(8) :: b          ! Arbitrary parameter 'b'

    ! Set weight to one by default
    site % wgt = ONE

    ! Sample position
    select case (external_source % type_space)
    case (SRC_SPACE_BOX)
       ! Coordinates sampled uniformly over a box
       p_min = external_source % params_space(1:3)
       p_max = external_source % params_space(4:6)
       r = (/ (prn(), i = 1,3) /)
       site % xyz = p_min + r*(p_max - p_min)

    case (SRC_SPACE_POINT)
       ! Point source
       site % xyz = external_source % params_space

    end select
    
    ! Sample angle
    select case (external_source % type_angle)
    case (SRC_ANGLE_ISOTROPIC)
       ! Sample isotropic distribution
       phi = TWO*PI*prn()
       mu = TWO*prn() - ONE
       site % uvw(1) = mu
       site % uvw(2) = sqrt(ONE - mu*mu) * cos(phi)
       site % uvw(3) = sqrt(ONE - mu*mu) * sin(phi)

    case (SRC_ANGLE_MONO)
       ! Monodirectional source
       site % uvw = external_source % params_angle

    case default
       message = "No angle distribution specified for external source!"
       call fatal_error()
    end select
          
    ! Sample energy distribution
    select case (external_source % type_energy)
    case (SRC_ENERGY_MONO)
       ! Monoenergtic source
       site % E = external_source % params_energy(1)

    case (SRC_ENERGY_MAXWELL)
       a = external_source % params_energy(1)
       do
          ! Sample Maxwellian fission spectrum
          site % E = maxwell_spectrum(a)

          ! resample if energy is >= 20 MeV
          if (site % E < 20) exit
       end do

    case (SRC_ENERGY_WATT)
       a = external_source % params_energy(1)
       b = external_source % params_energy(2)
       do
          ! Sample Watt fission spectrum
          site % E = watt_spectrum(a, b)

          ! resample if energy is >= 20 MeV
          if (site % E < 20) exit
       end do

    case default
       message = "No energy distribution specified for external source!"
       call fatal_error()
    end select

  end subroutine sample_external_source

!===============================================================================
! GET_SOURCE_PARTICLE returns the next source particle 
!===============================================================================

  subroutine get_source_particle(index_source)

    integer(8), intent(in) :: index_source

    integer(8) :: particle_seed  ! unique index for particle
    type(Bank), pointer :: src => null()

    ! set defaults
    call initialize_particle()

    ! Copy attributes from source to particle
    src => source_bank(index_source)
    call copy_source_attributes(src)

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
! COPY_SOURCE_ATTRIBUTES
!===============================================================================

  subroutine copy_source_attributes(src)

    type(Bank), pointer :: src

    ! copy attributes from source bank site
    p % wgt         = src % wgt
    p % last_wgt    = src % wgt
    p % coord % xyz = src % xyz
    p % coord % uvw = src % uvw
    p % last_xyz    = src % xyz
    p % E           = src % E
    p % last_E      = src % E

  end subroutine copy_source_attributes

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
    call MPI_FILE_OPEN(MPI_COMM_WORLD, path_source, MPI_MODE_CREATE + &
         MPI_MODE_WRONLY, MPI_INFO_NULL, fh, mpi_err)

    if (master) then
       offset = 0
       call MPI_FILE_WRITE_AT(fh, offset, n_particles, 1, MPI_INTEGER8, &
            MPI_STATUS_IGNORE, mpi_err)
    end if

    ! Set proper offset for source data on this processor
    offset = 8*(1 + rank*maxwork*8)

    ! Write all source sites
    call MPI_FILE_WRITE_AT(fh, offset, source_bank(1), work, MPI_BANK, &
         MPI_STATUS_IGNORE, mpi_err)

    ! Close binary source file
    call MPI_FILE_CLOSE(fh, mpi_err)

#else
    ! ==========================================================================
    ! SERIAL I/O USING FORTRAN INTRINSIC ROUTINES

    ! Open binary source file for writing
    open(UNIT=UNIT_SOURCE, FILE=path_source, STATUS='replace', &
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

    integer    :: i        ! loop over repeating sites
    integer(8) :: n_sites  ! number of sites in binary file
    integer    :: n_repeat ! number of times to repeat a site
#ifdef MPI
    integer    :: fh       ! file handle
    integer(MPI_OFFSET_KIND) :: offset ! offset in memory (0=beginning of file)
    integer    :: n_read   ! number of sites to read on a single process
#endif

#ifdef MPI
    ! ==========================================================================
    ! PARALLEL I/O USING MPI-2 ROUTINES

    ! Open binary source file for reading
    call MPI_FILE_OPEN(MPI_COMM_WORLD, path_source, MPI_MODE_RDONLY, &
         MPI_INFO_NULL, fh, mpi_err)

    ! Read number of source sites in file
    offset = 0
    call MPI_FILE_READ_AT(fh, offset, n_sites, 1, MPI_INTEGER8, &
         MPI_STATUS_IGNORE, mpi_err)

    if (n_particles > n_sites) then
       ! Determine number of sites to read and offset
       if (rank <= mod(n_sites,n_procs) - 1) then
          n_read = n_sites/n_procs + 1
          offset = 8*(1 + rank*n_read*8)
       else
          n_read = n_sites/n_procs
          offset = 8*(1 + (rank*n_read + mod(n_sites,n_procs))*8)
       end if

       ! Read source sites
       call MPI_FILE_READ_AT(fh, offset, source_bank(1), n_read, MPI_BANK, &
            MPI_STATUS_IGNORE, mpi_err)

       ! Let's say we have 30 sites and we need to fill in 200. This do loop
       ! will fill in sites 31 - 180.

       n_repeat = int(work / n_read)
       do i = 1, n_repeat - 1
          source_bank(i*n_read + 1:(i+1)*n_read) = &
               source_bank((i-1)*n_read + 1:i*n_read)
       end do

       ! This final statement would fill sites 181 - 200 in the above example.

       if (mod(work, n_repeat*n_read) > 0) then
          source_bank(n_repeat*n_read + 1:work) = &
               source_bank(1:work - n_repeat * n_read)
       end if
       
    else
       ! Set proper offset for source data on this processor
       offset = 8*(1 + rank*maxwork*8)

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
    open(UNIT=UNIT_SOURCE, FILE=path_source, STATUS='old', &
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
