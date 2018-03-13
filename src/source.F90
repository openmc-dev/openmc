module source

  use hdf5, only: HID_T
#ifdef OPENMC_MPI
  use message_passing
#endif

  use algorithm,        only: binary_search
  use bank_header,      only: Bank, source_bank
  use constants
  use distribution_univariate, only: Discrete
  use distribution_multivariate, only: SpatialBox
  use error,            only: fatal_error
  use geometry,         only: find_cell
  use hdf5_interface,   only: file_create, file_open, file_close, read_dataset
  use math
  use message_passing,  only: rank
  use mgxs_header,      only: rev_energy_bins, num_energy_groups
  use output,           only: write_message
  use particle_header,  only: Particle
  use random_lcg,       only: prn, set_particle_seed, prn_set_stream
  use settings
  use simulation_header
  use source_header,    only: external_source
  use string,           only: to_str
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
    real(8) :: xi

    ! Set the random number generator to the source stream.
    call prn_set_stream(STREAM_SOURCE)

    ! Sample from among multiple source distributions
    n_source = size(external_source)
    if (n_source > 1) then
      xi = prn()*sum(external_source(:) % strength)
      c = ZERO
      do i = 1, n_source
        c = c + external_source(i) % strength
        if (xi < c) exit
      end do
    else
      i = 1
    end if

    ! Sample source site from i-th source distribution
    site = external_source(i) % sample()

    ! If running in MG, convert site % E to group
    if (.not. run_CE) then
      site % E = real(binary_search(rev_energy_bins, num_energy_groups + 1, &
           site % E), 8)
      site % E = num_energy_groups + 1 - site % E
    end if

    ! Set the random number generator back to the tracking stream.
    call prn_set_stream(STREAM_TRACKING)

  end subroutine sample_external_source

end module source
