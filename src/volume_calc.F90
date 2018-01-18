module volume_calc

  use, intrinsic :: ISO_C_BINDING

  use hdf5, only: HID_T
#ifdef _OPENMP
  use omp_lib
#endif

  use constants
  use error,        only: write_message
  use geometry,     only: find_cell
  use geometry_header, only: universes, cells
  use hdf5_interface, only: file_create, file_close, write_attribute, &
       create_group, close_group, write_dataset
  use output,       only: header, time_stamp
  use material_header, only: materials
  use message_passing
  use nuclide_header, only: nuclides
  use particle_header, only: Particle
  use random_lcg,   only: prn, prn_set_stream, set_particle_seed
  use settings,     only: path_output
  use stl_vector,   only: VectorInt, VectorReal
  use string,       only: to_str
  use timer_header, only: Timer
  use volume_header

  implicit none
  private

  public :: openmc_calculate_volumes

contains

!===============================================================================
! OPENMC_CALCULATE_VOLUMES runs each of the stochastic volume calculations that
! the user has specified and writes results to HDF5 files
!===============================================================================

  subroutine openmc_calculate_volumes() bind(C)
    integer :: i, j
    integer :: n
    real(8), allocatable :: volume(:,:)  ! volume mean/stdev in each domain
    character(10) :: domain_type
    character(MAX_FILE_LEN) :: filename  ! filename for HDF5 file
    type(Timer) :: time_volume           ! timer for volume calculation
    type(VectorInt), allocatable :: nuclide_vec(:) ! indices in nuclides array
    type(VectorReal), allocatable :: atoms_vec(:) ! total # of atoms of each nuclide
    type(VectorReal), allocatable :: uncertainty_vec(:) ! uncertainty of total # of atoms

    if (master) then
      call header("STOCHASTIC VOLUME CALCULATION", 3)
      call time_volume % start()
    end if

    do i = 1, size(volume_calcs)
      n = size(volume_calcs(i) % domain_id)
      allocate(nuclide_vec(n))
      allocate(atoms_vec(n), uncertainty_vec(n))
      allocate(volume(2,n))

      if (master) then
        call write_message("Running volume calculation " // trim(to_str(i)) &
             // "...", 4)
      end if

      call get_volume(volume_calcs(i), volume, nuclide_vec, atoms_vec, &
           uncertainty_vec)

      if (master) then
        select case (volume_calcs(i) % domain_type)
        case (FILTER_CELL)
          domain_type = '  Cell'
        case (FILTER_MATERIAL)
          domain_type = '  Material'
        case (FILTER_UNIVERSE)
          domain_type = '  Universe'
        end select

        ! Display domain volumes
        do j = 1, size(volume_calcs(i) % domain_id)
          call write_message(trim(domain_type) // " " // trim(to_str(&
               volume_calcs(i) % domain_id(j))) // ": " // trim(to_str(&
               volume(1,j))) // " +/- " // trim(to_str(volume(2,j))) // &
               " cm^3", 4)
        end do
        call write_message("", 4)

        filename = trim(path_output) // 'volume_' // trim(to_str(i)) // '.h5'
        call write_volume(volume_calcs(i), filename, volume, nuclide_vec, &
             atoms_vec, uncertainty_vec)
      end if

      deallocate(nuclide_vec, atoms_vec, uncertainty_vec, volume)
    end do

    ! Show elapsed time
    if (master) then
      call time_volume % stop()
      call write_message("Elapsed time: " // trim(to_str(time_volume % &
           get_value())) // " s", 6)
    end if
  end subroutine openmc_calculate_volumes

!===============================================================================
! GET_VOLUME stochastically determines the volume of a set of domains along with
! the average number densities of nuclides within the domain
!===============================================================================

  subroutine get_volume(this, volume, nuclide_vec, atoms_vec, uncertainty_vec)
    type(VolumeCalculation), intent(in) :: this
    real(8),          intent(out) :: volume(:,:)     ! volume mean/stdev in each domain
    type(VectorInt),  intent(out) :: nuclide_vec(:)  ! indices in nuclides array
    type(VectorReal), intent(out) :: atoms_vec(:)    ! total # of atoms of each nuclide
    type(VectorReal), intent(out) :: uncertainty_vec(:) ! uncertainty of total # of atoms

    ! Variables that are private to each thread
    integer(8) :: i
    integer :: j, k
    integer :: i_domain      ! index in domain_id array
    integer :: i_material  ! index in materials array
    integer :: level       ! local coordinate level
    integer :: n_mat(size(this % domain_id))  ! Number of materials for each domain
    integer, allocatable :: indices(:,:) ! List of material indices for each domain
    integer, allocatable :: hits(:,:)    ! Number of hits for each material in each domain
    logical :: found_cell
    type(Particle) :: p

    ! Shared variables
    integer :: i_start, i_end  ! Starting/ending sample for each process
    type(VectorInt) :: master_indices(size(this % domain_id))
    type(VectorInt) :: master_hits(size(this % domain_id))

    ! Variables used outside of parallel region
    integer :: i_nuclide   ! index in nuclides array
    integer :: total_hits  ! total hits for a single domain (summed over materials)
    integer :: min_samples ! minimum number of samples per process
    integer :: remainder   ! leftover samples from uneven divide
#ifdef MPI
    integer :: mpi_err ! MPI error code
    integer :: m  ! index over materials
    integer :: n  ! number of materials
    integer, allocatable :: data(:) ! array used to send number of hits
#endif
    real(8) :: f              ! fraction of hits
    real(8) :: var_f          ! variance of fraction of hits
    real(8) :: volume_sample  ! total volume of sampled region
    real(8) :: atoms(2, size(nuclides))

    ! Divide work over MPI processes
    min_samples = this % samples / n_procs
    remainder = mod(this % samples, n_procs)
    if (rank < remainder) then
      i_start = (min_samples + 1)*rank
      i_end = i_start + min_samples
    else
      i_start = (min_samples + 1)*remainder + (rank - remainder)*min_samples
      i_end = i_start + min_samples - 1
    end if

    call p % initialize()

!$omp parallel private(i, j, k, i_domain, i_material, level, found_cell, &
!$omp&                 indices, hits, n_mat) firstprivate(p)

    ! Create space for material indices and number of hits for each
    allocate(indices(size(this % domain_id), 8))
    allocate(hits(size(this % domain_id), 8))
    n_mat(:) = 0

    call prn_set_stream(STREAM_VOLUME)

    ! ==========================================================================
    ! SAMPLES LOCATIONS AND COUNT HITS

!$omp do
    SAMPLE_LOOP: do i = i_start, i_end
      call set_particle_seed(i)

      p % n_coord = 1
      p % coord(1) % xyz(1) = this % lower_left(1) + prn()*(&
           this % upper_right(1) - this % lower_left(1))
      p % coord(1) % xyz(2) = this % lower_left(2) + prn()*(&
           this % upper_right(2) - this % lower_left(2))
      p % coord(1) % xyz(3) = this % lower_left(3) + prn()*(&
           this % upper_right(3) - this % lower_left(3))
      p % coord(1) % uvw(:) = [HALF, HALF, HALF]

      ! If this location is not in the geometry at all, move on to the next
      ! block
      call find_cell(p, found_cell)
      if (.not. found_cell) cycle

      if (this % domain_type == FILTER_MATERIAL) then
        i_material = p % material
        do i_domain = 1, size(this % domain_id)
          if (i_material == materials(i_domain) % id) then
            call check_hit(i_domain, i_material, indices, hits, n_mat)
          end if
        end do

      elseif (this % domain_type == FILTER_CELL) THEN
        do level = 1, p % n_coord
          do i_domain = 1, size(this % domain_id)
            if (cells(p % coord(level) % cell) % id == this % domain_id(i_domain)) then
              i_material = p % material
              call check_hit(i_domain, i_material, indices, hits, n_mat)
            end if
          end do
        end do

      elseif (this % domain_type == FILTER_UNIVERSE) then
        do level = 1, p % n_coord
          do i_domain = 1, size(this % domain_id)
            if (universes(p % coord(level) % universe) % id == &
                 this % domain_id(i_domain)) then
              i_material = p % material
              call check_hit(i_domain, i_material, indices, hits, n_mat)
            end if
          end do
        end do

      end if
    end do SAMPLE_LOOP
!$omp end do

    ! ==========================================================================
    ! REDUCE HITS ONTO MASTER THREAD

    ! At this point, each thread has its own pair of index/hits lists and we now
    ! need to reduce them. OpenMP is not nearly smart enough to do this on its
    ! own, so we have to manually reduce them.

#ifdef _OPENMP
!$omp do ordered schedule(static)
    THREAD_LOOP: do i = 1, omp_get_num_threads()
!$omp ordered
      do i_domain = 1, size(this % domain_id)
        INDEX_LOOP: do j = 1, n_mat(i_domain)
          ! Check if this material has been added to the master list and if so,
          ! accumulate the number of hits
          do k = 1, master_indices(i_domain) % size()
            if (indices(i_domain, j) == master_indices(i_domain) % data(k)) then
              master_hits(i_domain) % data(k) = &
                   master_hits(i_domain) % data(k) + hits(i_domain, j)
              cycle INDEX_LOOP
            end if
          end do

          ! If we made it here, this means the material hasn't yet been added to
          ! the master list, so add an entry to both the master indices and master
          ! hits lists
          call master_indices(i_domain) % push_back(indices(i_domain, j))
          call master_hits(i_domain) % push_back(hits(i_domain, j))
        end do INDEX_LOOP
      end do
!$omp end ordered
    end do THREAD_LOOP
!$omp end do
#else
    do i_domain = 1, size(this % domain_id)
      do j = 1, n_mat(i_domain)
        call master_indices(i_domain) % push_back(indices(i_domain, j))
        call master_hits(i_domain) % push_back(hits(i_domain, j))
      end do
    end do
#endif

    call prn_set_stream(STREAM_TRACKING)
!$omp end parallel

    ! ==========================================================================
    ! REDUCE HITS ONTO MASTER PROCESS

    volume_sample = product(this % upper_right - this % lower_left)

    do i_domain = 1, size(this % domain_id)
      atoms(:, :) = ZERO
      total_hits = 0

      if (master) then
#ifdef MPI
        do j = 1, n_procs - 1
          call MPI_RECV(n, 1, MPI_INTEGER, j, 0, mpi_intracomm, &
               MPI_STATUS_IGNORE, mpi_err)

          allocate(data(2*n))
          call MPI_RECV(data, 2*n, MPI_INTEGER, j, 1, mpi_intracomm, &
               MPI_STATUS_IGNORE, mpi_err)
          do k = 0, n - 1
            do m = 1, master_indices(i_domain) % size()
              if (data(2*k + 1) == master_indices(i_domain) % data(m)) then
                master_hits(i_domain) % data(m) = master_hits(i_domain) % data(m) + &
                     data(2*k + 2)
              end if
            end do
          end do
          deallocate(data)
        end do
#endif

        do j = 1, master_indices(i_domain) % size()
          total_hits = total_hits + master_hits(i_domain) % data(j)
          f = real(master_hits(i_domain) % data(j), 8) / this % samples
          var_f = f*(ONE - f) / this % samples

          i_material = master_indices(i_domain) % data(j)
          if (i_material == MATERIAL_VOID) cycle

          associate (mat => materials(i_material))
            do k = 1, size(mat % nuclide)
              ! Accumulate nuclide density
              i_nuclide = mat % nuclide(k)
              atoms(1, i_nuclide) = atoms(1, i_nuclide) + &
                   mat % atom_density(k) * f
              atoms(2, i_nuclide) = atoms(2, i_nuclide) + &
                   mat % atom_density(k)**2 * var_f
            end do
          end associate
        end do

        ! Determine volume
        volume(1, i_domain) = real(total_hits, 8) / this % samples * volume_sample
        volume(2, i_domain) = sqrt(volume(1, i_domain) * (volume_sample - &
             volume(1, i_domain)) / this % samples)

        ! Determine total number of atoms. At this point, we have values in
        ! atoms/b-cm. To get to atoms we multiple by 10^24 V.
        do j = 1, size(atoms, 2)
          atoms(1, j) = 1.0e24_8 * volume_sample * atoms(1, j)
          atoms(2, j) = 1.0e24_8 * volume_sample * sqrt(atoms(2, j))
        end do

        ! Convert full arrays to vectors
        do j = 1, size(nuclides)
          if (atoms(1, j) > ZERO) then
            call nuclide_vec(i_domain) % push_back(j)
            call atoms_vec(i_domain) % push_back(atoms(1, j))
            call uncertainty_vec(i_domain) % push_back(atoms(2, j))
          end if
        end do

      else
#ifdef MPI
        n = master_indices(i_domain) % size()
        allocate(data(2*n))
        do k = 0, n - 1
          data(2*k + 1) = master_indices(i_domain) % data(k + 1)
          data(2*k + 2) = master_hits(i_domain) % data(k + 1)
        end do

        call MPI_SEND(n, 1, MPI_INTEGER, 0, 0, mpi_intracomm, mpi_err)
        call MPI_SEND(data, 2*n, MPI_INTEGER, 0, 1, mpi_intracomm, mpi_err)
        deallocate(data)
#endif
      end if
    end do

  contains

    !===========================================================================
    ! CHECK_HIT is an internal subroutine that checks for whether a material has
    ! already been hit for a given domain. If not, it increases the list size by
    ! one (taking care of re-allocation if needed).
    !===========================================================================

    subroutine check_hit(i_domain, i_material, indices, hits, n_mat)
      integer :: i_domain
      integer :: i_material
      integer, allocatable :: indices(:,:)
      integer, allocatable :: hits(:,:)
      integer :: n_mat(:)

      integer, allocatable :: temp(:,:)
      logical :: already_hit
      integer :: j, k, nm

      ! Check if we've already had a hit in this material and if so,
      ! simply add one
      already_hit = .false.
      nm = n_mat(i_domain)
      do j = 1, nm
        if (indices(i_domain, j) == i_material) then
          hits(i_domain, j) = hits(i_domain, j) + 1
          already_hit = .true.
        end if
      end do

      if (.not. already_hit) then
        ! If we make it here, that means we haven't yet had a hit in this
        ! material. First check if the indices and hits arrays are large enough
        ! and if not, double them.
        if (nm == size(indices, 2)) then
          k = 2*size(indices, 2)
          allocate(temp(size(this % domain_id), k))
          temp(:, 1:nm) = indices(:, 1:nm)
          call move_alloc(FROM=temp, TO=indices)

          allocate(temp(size(this % domain_id), k))
          temp(:, 1:nm) = hits(:, 1:nm)
          call move_alloc(FROM=temp, TO=hits)
        end if

        ! Add an entry to both the indices list and the hits list
        n_mat(i_domain) = n_mat(i_domain) + 1
        indices(i_domain, n_mat(i_domain)) = i_material
        hits(i_domain, n_mat(i_domain)) = 1
      end if
    end subroutine check_hit

  end subroutine get_volume

!===============================================================================
! WRITE_VOLUME writes the results of a single stochastic volume calculation to
! an HDF5 file
!===============================================================================

  subroutine write_volume(this, filename, volume, nuclide_vec, atoms_vec, &
       uncertainty_vec)
    type(VolumeCalculation), intent(in) :: this
    character(*),     intent(in) :: filename       ! filename for HDF5 file
    real(8),          intent(in) :: volume(:,:)    ! volume mean/stdev in each domain
    type(VectorInt),  intent(in) :: nuclide_vec(:) ! indices in nuclides array
    type(VectorReal), intent(in) :: atoms_vec(:)   ! total # of atoms of each nuclide
    type(VectorReal), intent(in) :: uncertainty_vec(:) ! uncertainty of total # of atoms

    integer :: i, j
    integer :: n
    integer(HID_T) :: file_id
    integer(HID_T) :: group_id
    real(8), allocatable :: atom_data(:,:) ! mean/stdev of total # of atoms for
                                           ! each nuclide
    character(MAX_WORD_LEN), allocatable :: nucnames(:) ! names of nuclides

    ! Create HDF5 file
    file_id = file_create(filename)

    ! Write header info
    call write_attribute(file_id, "filetype", "volume")
    call write_attribute(file_id, "version", VERSION_VOLUME)
    call write_attribute(file_id, "openmc_version", VERSION)
#ifdef GIT_SHA1
    call write_attribute(file_id, "git_sha1", GIT_SHA1)
#endif

    ! Write current date and time
    call write_attribute(file_id, "date_and_time", time_stamp())

    ! Write basic metadata
    select case (this % domain_type)
    case (FILTER_CELL)
      call write_attribute(file_id, "domain_type", "cell")
    case (FILTER_MATERIAL)
      call write_attribute(file_id, "domain_type", "material")
    case (FILTER_UNIVERSE)
      call write_attribute(file_id, "domain_type", "universe")
    end select
    call write_attribute(file_id, "samples", this % samples)
    call write_attribute(file_id, "lower_left", this % lower_left)
    call write_attribute(file_id, "upper_right", this % upper_right)

    do i = 1, size(this % domain_id)
      group_id = create_group(file_id, "domain_" // trim(to_str(&
           this % domain_id(i))))

      ! Write volume for domain
      call write_dataset(group_id, "volume", volume(:, i))

      ! Create array of nuclide names from the vector
      n = nuclide_vec(i) % size()
      if (n > 0) then
        allocate(nucnames(n))
        do j = 1, n
          nucnames(j) = nuclides(nuclide_vec(i) % data(j)) % name
        end do

        ! Create array of total # of atoms with uncertainty for each nuclide
        allocate(atom_data(2, n))
        atom_data(1, :) = atoms_vec(i) % data(1:n)
        atom_data(2, :) = uncertainty_vec(i) % data(1:n)

        ! Write results
        call write_dataset(group_id, "nuclides", nucnames)
        call write_dataset(group_id, "atoms", atom_data)

        deallocate(nucnames)
        deallocate(atom_data)
      end if

      call close_group(group_id)
    end do
    call file_close(file_id)
  end subroutine write_volume

end module volume_calc
