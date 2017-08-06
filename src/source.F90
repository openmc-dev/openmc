module source

  use hdf5, only: HID_T
#ifdef MPI
  use message_passing
#endif

  use algorithm,        only: binary_search
  use bank_header,      only: Bank
  use constants
  use distribution_univariate, only: Discrete
  use distribution_multivariate, only: SpatialBox
  use error,            only: fatal_error
  use geometry,         only: find_cell
  use global
  use hdf5_interface,   only: file_create, file_open, file_close, read_dataset
  use message_passing,  only: rank
  use output,           only: write_message
  use particle_header,  only: Particle
  use random_lcg,       only: prn, set_particle_seed, prn_set_stream
  use string,           only: to_str
  use math
  use state_point,      only: read_source_bank, write_source_bank

  implicit none

contains

!===============================================================================
! INITIALIZE_SOURCE initializes particles in the source bank
!===============================================================================

  subroutine initialize_source()

    integer(8) :: i          ! loop index over bank sites
    integer(8) :: id         ! particle id
    integer(HID_T) :: file_id
    character(MAX_WORD_LEN) :: filetype
    character(MAX_FILE_LEN) :: filename
    type(Bank), pointer :: src ! source bank site

    call write_message("Initializing source particles...", 5)

    if (path_source /= '') then
      ! Read the source from a binary file instead of sampling from some
      ! assumed source distribution

      call write_message('Reading source file from ' // trim(path_source) &
           &// '...', 6)

      ! Open the binary file
      file_id = file_open(path_source, 'r', parallel=.true.)

      ! Read the file type
      call read_dataset(filetype, file_id, "filetype")

      ! Check to make sure this is a source file
      if (filetype /= 'source') then
        call fatal_error("Specified starting source file not a source file &
             &type.")
      end if

      ! Read in the source bank
      call read_source_bank(file_id)

      ! Close file
      call file_close(file_id)

    else
      ! Generation source sites from specified distribution in user input
      do i = 1, work
        ! Get pointer to source bank site
        src => source_bank(i)

        ! initialize random number seed
        id = total_gen*n_particles + work_index(rank) + i
        call set_particle_seed(id)

        ! sample external source distribution
        call sample_external_source(src)
      end do
    end if

    ! Write out initial source
    if (write_initial_source) then
      call write_message('Writing out initial source...', 5)
      filename = trim(path_output) // 'initial_source.h5'
      file_id = file_create(filename, parallel=.true.)
      call write_source_bank(file_id)
      call file_close(file_id)
    end if

  end subroutine initialize_source

!===============================================================================
! SAMPLE_EXTERNAL_SOURCE samples the user-specified external source and stores
! the position, angle, and energy in a Bank type.
!===============================================================================

  subroutine sample_external_source(site)
    type(Bank), intent(inout) :: site ! source site

    integer :: i          ! dummy loop index
    integer :: n_source   ! number of source distributions
    real(8) :: c          ! cumulative frequency
    real(8) :: r(3)       ! sampled coordinates
    logical :: found      ! Does the source particle exist within geometry?
    type(Particle) :: p   ! Temporary particle for using find_cell
    integer, save :: n_accept = 0  ! Number of samples accepted
    integer, save :: n_reject = 0  ! Number of samples rejected

    ! Set weight to one by default
    site % wgt = ONE

    ! Set the random number generator to the source stream.
    call prn_set_stream(STREAM_SOURCE)

    ! Sample from among multiple source distributions
    n_source = size(external_source)
    if (n_source > 1) then
      r(1) = prn()*sum(external_source(:) % strength)
      c = ZERO
      do i = 1, n_source
        c = c + external_source(i) % strength
        if (r(1) < c) exit
      end do
    else
      i = 1
    end if

    ! Repeat sampling source location until a good site has been found
    found = .false.
    do while (.not. found)
      ! Set particle defaults
      call p % initialize()

      ! Sample spatial distribution
      site % xyz(:) = external_source(i) % space % sample()

      ! Fill p with needed data
      p % coord(1) % xyz(:) = site % xyz
      p % coord(1) % uvw(:) = [ ONE, ZERO, ZERO ]

      ! Now search to see if location exists in geometry
      call find_cell(p, found)

      ! Check if spatial site is in fissionable material
      if (found) then
        select type (space => external_source(i) % space)
        type is (SpatialBox)
          if (space % only_fissionable) then
            if (p % material == MATERIAL_VOID) then
              found = .false.
            elseif (.not. materials(p % material) % fissionable) then
              found = .false.
            end if
          end if
        end select
      end if

      ! Check for rejection
      if (.not. found) then
        n_reject = n_reject + 1
        if (n_reject >= EXTSRC_REJECT_THRESHOLD .and. &
             real(n_accept, 8)/n_reject <= EXTSRC_REJECT_FRACTION) then
          call fatal_error("More than 95% of external source sites sampled &
               &were rejected. Please check your external source definition.")
        end if
      end if
    end do

    ! Increment number of accepted samples
    n_accept = n_accept + 1

    call p % clear()

    ! Sample angle
    site % uvw(:) = external_source(i) % angle % sample()

    ! Check for monoenergetic source above maximum neutron energy
    select type (energy => external_source(i) % energy)
    type is (Discrete)
      if (any(energy % x >= energy_max_neutron)) then
        call fatal_error("Source energy above range of energies of at least &
             &one cross section table")
      end if
    end select

    do
      ! Sample energy spectrum
      site % E = external_source(i) % energy % sample()

      ! resample if energy is greater than maximum neutron energy
      if (site % E < energy_max_neutron) exit
    end do

    ! Set delayed group
    site % delayed_group = 0

    ! If running in MG, convert site % E to group
    if (.not. run_CE) then
      if (site % E <= energy_bins(1)) then
        site % E = real(1, 8)
      else if (site % E > energy_bins(num_energy_groups + 1)) then
        site % E = real(num_energy_groups, 8)
      else
        site % E = real(binary_search(energy_bins, num_energy_groups + 1, &
             site % E), 8)
      end if
    end if

    ! Set the random number generator back to the tracking stream.
    call prn_set_stream(STREAM_TRACKING)

  end subroutine sample_external_source

end module source
