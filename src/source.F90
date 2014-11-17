module source

  use bank_header,        only: Bank
  use constants
  use dd_comm,            only: distribute_source
  use error,              only: fatal_error
  use geometry,           only: find_cell
  use geometry_header,    only: BASE_UNIVERSE
  use global
  use math,               only: maxwell_spectrum, watt_spectrum
  use output,             only: write_message
  use output_interface,   only: BinaryOutput
  use particle_header,    only: Particle
  use random_lcg,         only: prn, set_particle_seed, prn_set_stream, prn_seed
  use random_lcg_header,  only: STREAM_TRACKING, STREAM_SOURCE
  use string,             only: to_str

#ifdef HDF5
  use hdf5_interface
#endif
#ifdef MPI
  use mpi
#endif

  implicit none

contains

!===============================================================================
! INITIALIZE_SOURCE initializes particles in the source bank
!===============================================================================

  subroutine initialize_source()

    character(MAX_FILE_LEN) :: filename
    integer(8) :: i          ! loop index over bank sites
    integer(8) :: id         ! particle id
    integer(4) :: itmp       ! temporary integer
    type(Bank), pointer :: src => null() ! source bank site
    type(BinaryOutput) :: sp ! statepoint/source binary file

    call write_message("Initializing source particles...", 6)

    if (path_source /= '') then
      ! Read the source from a binary file instead of sampling from some
      ! assumed source distribution

      if (dd_run) then
        call fatal_error("Reading source from binary file not " // &
                         "implementedfor domain decomposition.")
      end if

      call write_message('Reading source file from ' // trim(path_source) &
           &// '...', 6)

      ! Open the binary file
      call sp % file_open(path_source, 'r', serial = .false.)

      ! Read the file type
      call sp % read_data(itmp, "filetype")

      ! Check to make sure this is a source file
      if (itmp /= FILETYPE_SOURCE) then
        call fatal_error("Specified starting source file not a source file &
             &type.")
      end if

      ! Read in the source bank
      call read_source_bank(sp)

      ! Close file
      call sp % file_close()

    else
      ! Generate source sites from specified distribution in user input
      
      do i = 1, work
      
        ! Get pointer to source bank site
        src => source_bank(i)

        ! Initialize random number seed
        id = work_index(rank) + i
        call set_particle_seed(id)

        ! Sample external source distribution
        call sample_external_source(src)
        
        ! Store the LCG seed for the particle
        src % prn_seed = prn_seed
        
      end do

      if (dd_run) then

        ! Send source sites to the processes they'll be transported on
        call distribute_source(domain_decomp)
        
      end if
      
    end if

    ! Write out initial source
    if (write_initial_source) then
      call write_message('Writing out initial source guess...', 1)
#ifdef HDF5
      filename = trim(path_output) // 'initial_source.h5'
#else
      filename = trim(path_output) // 'initial_source.binary'
#endif
      call sp % file_create(filename, serial = .false.)
      call write_source_bank(sp)
      call sp % file_close()
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
    logical :: found      ! Does the source particle exist within geometry?
    type(Particle) :: p   ! Temporary particle for using find_cell
    integer, save :: num_resamples = 0 ! Number of resamples encountered

    ! Set weight to one by default
    site % wgt = ONE

    ! Set the random number generator to the source stream.
    call prn_set_stream(STREAM_SOURCE)

    ! Sample position
    select case (external_source % type_space)
    case (SRC_SPACE_BOX)
      ! Set particle defaults
      call p % initialize()
      ! Repeat sampling source location until a good site has been found
      found = .false.
      do while (.not.found)
        ! Coordinates sampled uniformly over a box
        p_min = external_source % params_space(1:3)
        p_max = external_source % params_space(4:6)
        r = (/ (prn(), i = 1,3) /)
        site % xyz = p_min + r*(p_max - p_min)

        ! Fill p with needed data
        p % coord0 % xyz = site % xyz
        p % coord0 % uvw = [ ONE, ZERO, ZERO ]

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
      call p % clear()

    case (SRC_SPACE_FISSION)
      ! Repeat sampling source location until a good site has been found
      found = .false.
      do while (.not.found)
        ! Set particle defaults
        call p % initialize()

        ! Coordinates sampled uniformly over a box
        p_min = external_source % params_space(1:3)
        p_max = external_source % params_space(4:6)
        r = (/ (prn(), i = 1,3) /)
        site % xyz = p_min + r*(p_max - p_min)

        ! Fill p with needed data
        p % coord0 % xyz = site % xyz
        p % coord0 % uvw = [ ONE, ZERO, ZERO ]

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
        if (p % material == MATERIAL_VOID) then
          found = .false.
          cycle
        end if
        if (.not. materials(p % material) % fissionable) found = .false.
      end do
      call p % clear()

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
      call fatal_error("No angle distribution specified for external source!")
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
      call fatal_error("No energy distribution specified for external source!")
    end select

    ! Set the random number generator back to the tracking stream.
    call prn_set_stream(STREAM_TRACKING)

  end subroutine sample_external_source

!===============================================================================
! GET_SOURCE_PARTICLE returns the next source particle
!===============================================================================

  subroutine get_source_particle(p, index_source)

    type(Particle), intent(inout) :: p
    integer(8),     intent(in)    :: index_source

    integer(8) :: particle_seed  ! unique index for particle
    integer :: i
    type(Bank), pointer, save :: src => null()
!$omp threadprivate(src)

    if (.not. dd_run) then

      ! set defaults
      call p % initialize()

      ! Copy attributes from source to particle
      src => source_bank(index_source)
      call copy_source_attributes(p, src)
    
      ! set identifier for particle
      p % id = work_index(rank) + index_source
      
      ! set random number seed
      particle_seed = (overall_gen - 1)*n_particles + p % id
      call set_particle_seed(particle_seed)
    
    else

      if (current_stage == 1) then

        ! On the first stage, we set the particle from the source bank like for
        ! non-DD runs

        ! set defaults
        call p % initialize()

        ! Copy attributes from source to particle
        src => source_bank(index_source)
        call copy_source_attributes(p, src)

      else

        ! On all other stages, we're running a particle in the buffer
        ! transferred from other domains, so all we need to do is set the coord
        call p % clear()
        allocate(p % coord0)
        p % coord0 % universe = BASE_UNIVERSE
        p % coord             => p % coord0
        p % coord0 % xyz =  p % stored_xyz
        p % coord0 % uvw =  p % stored_uvw

        ! Re-initialize the stored outscatter_destination
        p % outscatter_destination = NO_OUTSCATTER

      end if

      ! The order of the source bank is scrambled for domain-decomposed runs,
      ! so unless we want to do a distributed sort over all n_partcles, we're
      ! not going to have particle id numbers that are reproducible with
      ! different numbers of processes.
      p % id = index_source

      ! Results will still be reproducible for DD runs with varying numbers of
      ! domains/proceses since we're carrying around the random number streams
      prn_seed = p % prn_seed

      ! Note: tally results will not match the non-DD version of the same run
      ! (unless you run in DD mode with one domain).

    end if

    ! set particle trace
    trace = .false.
    if (current_batch == trace_batch .and. current_gen == trace_gen .and. &
         p % id == trace_particle) trace = .true.
    
    ! Set particle track.
    p % write_track = .false.
    if (write_all_tracks) then
      p % write_track = .true.
    else if (allocated(track_identifiers)) then
      do i=1, size(track_identifiers(1,:))
        if (current_batch == track_identifiers(1,i) .and. &
             &current_gen == track_identifiers(2,i) .and. &
             &p % id == track_identifiers(3,i)) then
          p % write_track = .true.
          exit
        end if
      end do
    end if

  end subroutine get_source_particle

!===============================================================================
! COPY_SOURCE_ATTRIBUTES
!===============================================================================

  subroutine copy_source_attributes(p, src)

    type(Particle), intent(inout) :: p
    type(Bank),     pointer       :: src

    ! copy attributes from source bank site
    p % wgt         = src % wgt
    p % last_wgt    = src % wgt
    p % coord % xyz = src % xyz
    p % coord % uvw = src % uvw
    p % last_xyz    = src % xyz
    p % last_uvw    = src % uvw
    p % E           = src % E
    p % last_E      = src % E
    p % prn_seed    = src % prn_seed
    p % xs_seed     = src % prn_seed

  end subroutine copy_source_attributes

!===============================================================================
! WRITE_SOURCE_BANK writes OpenMC source_bank data
!===============================================================================

  subroutine write_source_bank(sp)

    type(BinaryOutput), intent(inout) :: sp

#ifdef MPI
    integer                  :: mpiio_err
# ifndef HDF5
    integer(MPI_OFFSET_KIND) :: offset           ! offset of data
    integer                  :: size_offset_kind ! the data offset kind
    integer                  :: size_bank        ! size of bank to write
# endif
# ifdef HDF5
    integer(8)               :: offset(1)        ! source data offset
# endif
#endif

#ifdef HDF5
# ifdef MPI

    ! Set size of total dataspace for all procs and rank
    dims1(1) = n_particles
    hdf5_rank = 1

    ! Create that dataspace
    call h5screate_simple_f(hdf5_rank, dims1, dspace, hdf5_err)

    ! Create the dataset for that dataspace
    call h5dcreate_f(sp % hdf5_fh, "source_bank", hdf5_bank_t, dspace, dset, hdf5_err)

    ! Close the dataspace
    call h5sclose_f(dspace, hdf5_err)

    ! Create another data space but for each proc individually
    dims1(1) = work
    call h5screate_simple_f(hdf5_rank, dims1, memspace, hdf5_err)

    ! Get the individual local proc dataspace
    call h5dget_space_f(dset, dspace, hdf5_err)

    ! Select hyperslab for this dataspace
    offset(1) = work_index(rank)
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, offset, dims1, hdf5_err)

    ! Set up the property list for parallel writing
    call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)
    call h5pset_dxpl_mpio_f(plist, H5FD_MPIO_COLLECTIVE_F, hdf5_err)

    ! Set up pointer to data
    f_ptr = c_loc(source_bank(1))

    ! Write data to file in parallel
    call h5dwrite_f(dset, hdf5_bank_t, f_ptr, hdf5_err, &
         file_space_id = dspace, mem_space_id = memspace, &
         xfer_prp = plist)

    ! Close all ids
    call h5sclose_f(dspace, hdf5_err)
    call h5sclose_f(memspace, hdf5_err)
    call h5dclose_f(dset, hdf5_err)
    call h5pclose_f(plist, hdf5_err)

# else

    ! Set size
    dims1(1) = work
    hdf5_rank = 1

    ! Create dataspace
    call h5screate_simple_f(hdf5_rank, dims1, dspace, hdf5_err)

    ! Create dataset
    call h5dcreate_f(sp % hdf5_fh, "source_bank", hdf5_bank_t, &
         dspace, dset, hdf5_err)

    ! Set up pointer to data
    f_ptr = c_loc(source_bank(1))

    ! Write dataset to file
    call h5dwrite_f(dset, hdf5_bank_t, f_ptr, hdf5_err)

    ! Close all ids
    call h5dclose_f(dset, hdf5_err)
    call h5sclose_f(dspace, hdf5_err)

# endif 

#elif MPI

    ! Get current offset for master 
    if (master) call MPI_FILE_GET_POSITION(sp % unit_fh, offset, mpiio_err)

    ! Determine offset on master process and broadcast to all processors
    call MPI_SIZEOF(offset, size_offset_kind, mpi_err)
    select case (size_offset_kind)
    case (4)
      call MPI_BCAST(offset, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)
    case (8)
      call MPI_BCAST(offset, 1, MPI_INTEGER8, 0, MPI_COMM_WORLD, mpi_err)
    end select

    ! Set the proper offset for source data on this processor
    call MPI_TYPE_SIZE(MPI_BANK, size_bank, mpi_err)
    offset = offset + size_bank*work_index(rank)

    ! Write all source sites
    call MPI_FILE_WRITE_AT(sp % unit_fh, offset, source_bank(1), work, MPI_BANK, &
         MPI_STATUS_IGNORE, mpiio_err)

#else

    ! Write out source sites
    write(sp % unit_fh) source_bank

#endif

  end subroutine write_source_bank

!===============================================================================
! READ_SOURCE_BANK reads OpenMC source_bank data
!===============================================================================

  subroutine read_source_bank(sp)

    type(BinaryOutput), intent(inout) :: sp

#ifdef MPI
    integer                  :: mpiio_err
# ifndef HDF5
    integer(MPI_OFFSET_KIND) :: offset           ! offset of data
    integer                  :: size_bank        ! size of bank to read
# endif
# ifdef HDF5
    integer(8)               :: offset(1)        ! offset of data
# endif
#endif

#ifdef HDF5
# ifdef MPI

    ! Set size of total dataspace for all procs and rank
    dims1(1) = n_particles
    hdf5_rank = 1

    ! Open the dataset
    call h5dopen_f(sp % hdf5_fh, "source_bank", dset, hdf5_err)

    ! Create another data space but for each proc individually
    dims1(1) = work
    call h5screate_simple_f(hdf5_rank, dims1, memspace, hdf5_err)

    ! Get the individual local proc dataspace
    call h5dget_space_f(dset, dspace, hdf5_err)

    ! Select hyperslab for this dataspace
    offset(1) = work_index(rank)
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, offset, dims1, hdf5_err)

    ! Set up the property list for parallel writing
    call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)
    call h5pset_dxpl_mpio_f(plist, H5FD_MPIO_COLLECTIVE_F, hdf5_err)

    ! Set up pointer to data
    f_ptr = c_loc(source_bank(1))

    ! Read data from file in parallel
    call h5dread_f(dset, hdf5_bank_t, f_ptr, hdf5_err, &
         file_space_id = dspace, mem_space_id = memspace, &
         xfer_prp = plist)

    ! Close all ids
    call h5sclose_f(dspace, hdf5_err)
    call h5sclose_f(memspace, hdf5_err)
    call h5dclose_f(dset, hdf5_err)
    call h5pclose_f(plist, hdf5_err)

# else

    ! Open dataset
    call h5dopen_f(sp % hdf5_fh, "source_bank", dset, hdf5_err)

    ! Set up pointer to data
    f_ptr = c_loc(source_bank(1))

    ! Read dataset from file
    call h5dread_f(dset, hdf5_bank_t, f_ptr, hdf5_err)

    ! Close all ids
    call h5dclose_f(dset, hdf5_err)

# endif 

#elif MPI

    ! Go to the end of the file to set file pointer
    offset = 0
    call MPI_FILE_SEEK(sp % unit_fh, offset, MPI_SEEK_END, &
         mpiio_err)

    ! Get current offset (will be at EOF) 
    call MPI_FILE_GET_POSITION(sp % unit_fh, offset, mpiio_err)

    ! Get the size of the source bank on all procs
    call MPI_TYPE_SIZE(MPI_BANK, size_bank, mpi_err)

    ! Calculate offset where the source bank will begin
    offset = offset - n_particles*size_bank

    ! Set the proper offset for source data on this processor
    offset = offset + size_bank*work_index(rank)

    ! Write all source sites
    call MPI_FILE_READ_AT(sp % unit_fh, offset, source_bank(1), work, MPI_BANK, &
         MPI_STATUS_IGNORE, mpiio_err)

#else

    ! Write out source sites
    read(sp % unit_fh) source_bank

#endif

  end subroutine read_source_bank

end module source
