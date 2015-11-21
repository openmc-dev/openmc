module source

  use bank_header,      only: Bank
  use constants
  use error,            only: fatal_error
  use geometry,         only: find_cell
  use geometry_header,  only: BASE_UNIVERSE
  use global
  use hdf5_interface,   only: file_create, file_open, file_close, read_dataset
  use math,             only: maxwell_spectrum, watt_spectrum
  use output,           only: write_message
  use particle_header,  only: Particle
  use random_lcg,       only: prn, set_particle_seed, prn_set_stream
  use state_point,      only: read_source_bank, write_source_bank
  use string,           only: to_str

#ifdef MPI
  use message_passing
#endif

  use hdf5, only: HID_T

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

    call write_message("Initializing source particles...", 6)

    if (path_source /= '') then
      ! Read the source from a binary file instead of sampling from some
      ! assumed source distribution

      call write_message('Reading source file from ' // trim(path_source) &
           &// '...', 6)

      ! Open the binary file
      file_id = file_open(path_source, 'r', parallel=.true.)

      ! Read the file type
      call read_dataset(file_id, "filetype", filetype)

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
        id = work_index(rank) + i
        call set_particle_seed(id)

        ! sample external source distribution
        call sample_external_source(src)
      end do
    end if

    ! Write out initial source
    if (write_initial_source) then
      call write_message('Writing out initial source...', 1)
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
    real(8) :: r(3)       ! sampled coordinates
    real(8) :: phi        ! azimuthal angle
    real(8) :: mu         ! cosine of polar angle
    real(8) :: p_min(3)   ! minimum coordinates of source
    real(8) :: p_max(3)   ! maximum coordinates of source
    real(8) :: a          ! Arbitrary parameter 'a'
    real(8) :: b          ! Arbitrary parameter 'b'
    logical :: found      ! Does the source particle exist within geometry?
    type(Particle) :: p   ! Temporary particle for using find_cell
    integer, save :: num_resamples = 0 ! Number of resamples encountered

    ! Set weight to one by default
    site%wgt = ONE

    ! Set the random number generator to the source stream.
    call prn_set_stream(STREAM_SOURCE)

    ! Sample position
    select case (external_source%type_space)
    case (SRC_SPACE_BOX)
      ! Set particle defaults
      call p%initialize()
      ! Repeat sampling source location until a good site has been found
      found = .false.
      do while (.not.found)
        ! Coordinates sampled uniformly over a box
        p_min = external_source%params_space(1:3)
        p_max = external_source%params_space(4:6)
        r = (/ (prn(), i = 1,3) /)
        site%xyz = p_min + r*(p_max - p_min)

        ! Fill p with needed data
        p%coord(1)%xyz = site%xyz
        p%coord(1)%uvw = [ ONE, ZERO, ZERO ]

        ! Now search to see if location exists in geometry
        call find_cell(p, found)
        if (.not. found) then
          num_resamples = num_resamples + 1
          if (num_resamples == MAX_EXTSRC_RESAMPLES) then
            call fatal_error("Maximum number of external source spatial &
                 &resamples reached!")
          end if
        end if
      end do
      call p%clear()

    case (SRC_SPACE_FISSION)
      ! Repeat sampling source location until a good site has been found
      found = .false.
      do while (.not.found)
        ! Set particle defaults
        call p%initialize()

        ! Coordinates sampled uniformly over a box
        p_min = external_source%params_space(1:3)
        p_max = external_source%params_space(4:6)
        r = (/ (prn(), i = 1,3) /)
        site%xyz = p_min + r*(p_max - p_min)

        ! Fill p with needed data
        p%coord(1)%xyz = site%xyz
        p%coord(1)%uvw = [ ONE, ZERO, ZERO ]

        ! Now search to see if location exists in geometry
        call find_cell(p, found)
        if (.not. found) then
          num_resamples = num_resamples + 1
          if (num_resamples == MAX_EXTSRC_RESAMPLES) then
            call fatal_error("Maximum number of external source spatial &
                 &resamples reached!")
          end if
          cycle
        end if
        if (p%material == MATERIAL_VOID) then
          found = .false.
          cycle
        end if
        if (.not. materials(p%material)%fissionable) found = .false.
      end do
      call p%clear()

    case (SRC_SPACE_POINT)
      ! Point source
      site%xyz = external_source%params_space

    end select

    ! Sample angle
    select case (external_source%type_angle)
    case (SRC_ANGLE_ISOTROPIC)
      ! Sample isotropic distribution
      phi = TWO*PI*prn()
      mu = TWO*prn() - ONE
      site%uvw(1) = mu
      site%uvw(2) = sqrt(ONE - mu*mu) * cos(phi)
      site%uvw(3) = sqrt(ONE - mu*mu) * sin(phi)

    case (SRC_ANGLE_MONO)
      ! Monodirectional source
      site%uvw = external_source%params_angle

    case default
      call fatal_error("No angle distribution specified for external source!")
    end select

    ! Sample energy distribution
    select case (external_source%type_energy)
    case (SRC_ENERGY_MONO)
      ! Monoenergtic source
      site%E = external_source%params_energy(1)
      if (site%E >= energy_max_neutron) then
        call fatal_error("Source energy above range of energies of at least &
             &one cross section table")
      end if

    case (SRC_ENERGY_MAXWELL)
      a = external_source%params_energy(1)
      do
        ! Sample Maxwellian fission spectrum
        site%E = maxwell_spectrum(a)

        ! resample if energy is greater than maximum neutron energy
        if (site%E < energy_max_neutron) exit
      end do

    case (SRC_ENERGY_WATT)
      a = external_source%params_energy(1)
      b = external_source%params_energy(2)
      do
        ! Sample Watt fission spectrum
        site%E = watt_spectrum(a, b)

        ! resample if energy is greater than maximum neutron energy
        if (site%E < energy_max_neutron) exit
      end do

    case default
      call fatal_error("No energy distribution specified for external source!")
    end select

    ! Set the random number generator back to the tracking stream.
    call prn_set_stream(STREAM_TRACKING)

  end subroutine sample_external_source

end module source
