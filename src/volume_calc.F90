module volume_calc

  use hdf5, only: HID_T
#ifdef _OPENMP
  use omp_lib
#endif

  use constants
  use geometry,     only: find_cell
  use global
  use hdf5_interface, only: file_create, file_close, write_attribute, &
       create_group, close_group, write_dataset
  use output,       only: write_message, header
  use message_passing
  use particle_header, only: Particle
  use random_lcg,   only: prn, prn_set_stream, set_particle_seed
  use stl_vector,   only: VectorInt, VectorReal
  use timer_header, only: Timer
  use volume_header

  implicit none
  private

  public :: run_volume_calculations

contains

!===============================================================================
! RUN_VOLUME_CALCULATIONS runs each of the stochastic volume calculations that
! the user has specified and writes results to HDF5 files
!===============================================================================

  subroutine run_volume_calculations()
    integer :: i, j
    integer :: n
    real(8), allocatable :: volume(:,:)  ! volume mean/stdev in each cell
    character(MAX_FILE_LEN) :: filename  ! filename for HDF5 file
    type(Timer) :: time_volume           ! timer for volume calculation
    type(VectorInt), allocatable :: nuclide_vec(:) ! indices in nuclides array
    type(VectorReal), allocatable :: atoms_vec(:) ! total # of atoms of each nuclide
    type(VectorReal), allocatable :: uncertainty_vec(:) ! uncertainty of total # of atoms

    if (master) then
      call header("STOCHASTIC VOLUME CALCULATION", level=1)
      call time_volume % start()
    end if

    do i = 1, size(volume_calcs)
      n = size(volume_calcs(i) % cell_id)
      allocate(nuclide_vec(n))
      allocate(atoms_vec(n), uncertainty_vec(n))
      allocate(volume(2,n))

      if (master) then
        call write_message("Running volume calculation " // trim(to_str(i)) &
             // "...")
      end if

      call get_volume(volume_calcs(i), volume, nuclide_vec, atoms_vec, &
           uncertainty_vec)

      if (master) then
        ! Display cell volumes
        do j = 1, size(volume_calcs(i) % cell_id)
          call write_message("  Cell " // trim(to_str(volume_calcs(i) % &
               cell_id(j))) // ": " // trim(to_str(volume(1,j))) // " +/- " // &
               trim(to_str(volume(2,j))) // " cm^3")
        end do
        call write_message("")

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
           get_value())) // " s")
    end if
  end subroutine run_volume_calculations

!===============================================================================
! GET_VOLUME stochastically determines the volume of a set of cells along with
! the average number densities of nuclides within the cell
!===============================================================================

  subroutine get_volume(this, volume, nuclide_vec, atoms_vec, uncertainty_vec)
    type(VolumeCalculation), intent(in) :: this
    real(8),          intent(out) :: volume(:,:)     ! volume mean/stdev in each cell
    type(VectorInt),  intent(out) :: nuclide_vec(:)  ! indices in nuclides array
    type(VectorReal), intent(out) :: atoms_vec(:)    ! total # of atoms of each nuclide
    type(VectorReal), intent(out) :: uncertainty_vec(:) ! uncertainty of total # of atoms

    ! Variables that are private to each thread
    integer(8) :: i
    integer :: j, k
    integer :: i_cell      ! index in cell_id array
    integer :: i_material  ! index in materials array
    integer :: level       ! local coordinate level
    logical :: found_cell
    type(VectorInt) :: indices(size(this % cell_id)) ! List of material indices
    type(VectorInt) :: hits(size(this % cell_id))    ! Number of hits for each material
    type(Particle) :: p

    ! Shared variables
    integer :: i_start, i_end  ! Starting/ending sample for each process
    type(VectorInt) :: master_indices(size(this % cell_id))
    type(VectorInt) :: master_hits(size(this % cell_id))

    ! Variables used outside of parallel region
    integer :: i_nuclide    ! index in nuclides array
    integer :: total_hits   ! total hits for a single cell (summed over materials)
    integer :: min_samples ! minimum number of samples per process
    integer :: remainder        ! leftover samples from uneven divide
#ifdef MPI
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

!$omp parallel private(i, j, k, i_cell, i_material, level, found_cell) &
!$omp&         firstprivate(p, indices, hits)

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

      ! Determine if point is within desired cell
      LEVEL_LOOP: do level = 1, p % n_coord
        CELL_CHECK_LOOP: do i_cell = 1, size(this % cell_id)
          if (cells(p % coord(level) % cell) % id == this % cell_id(i_cell)) then

            ! Determine what material this is
            i_material = p % material

            ! Check if we've already had a hit in this material and if so,
            ! simply add one
            do j = 1, indices(i_cell) % size()
              if (indices(i_cell) % data(j) == i_material) then
                hits(i_cell) % data(j) = hits(i_cell) % data(j) + 1
                cycle CELL_CHECK_LOOP
              end if
            end do

            ! If we make it here, that means we haven't yet had a hit in this
            ! material. Add an entry to both the indices list and the hits list
            call indices(i_cell) % push_back(i_material)
            call hits(i_cell) % push_back(1)
          end if
        end do CELL_CHECK_LOOP
      end do LEVEL_LOOP
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
      do i_cell = 1, size(this % cell_id)
        INDEX_LOOP: do j = 1, indices(i_cell) % size()
          ! Check if this material has been added to the master list and if so,
          ! accumulate the number of hits
          do k = 1, master_indices(i_cell) % size()
            if (indices(i_cell) % data(j) == master_indices(i_cell) % data(k)) then
              master_hits(i_cell) % data(k) = &
                   master_hits(i_cell) % data(k) + hits(i_cell) % data(j)
              cycle INDEX_LOOP
            end if
          end do

          ! If we made it here, this means the material hasn't yet been added to
          ! the master list, so add an entry to both the master indices and master
          ! hits lists
          call master_indices(i_cell) % push_back(indices(i_cell) % data(j))
          call master_hits(i_cell) % push_back(hits(i_cell) % data(j))
        end do INDEX_LOOP
      end do
!$omp end ordered
    end do THREAD_LOOP
!$omp end do
#else
    master_indices = indices
    master_hits = hits
#endif

    call prn_set_stream(STREAM_TRACKING)
!$omp end parallel

    ! ==========================================================================
    ! REDUCE HITS ONTO MASTER PROCESS

    volume_sample = product(this % upper_right - this % lower_left)

    do i_cell = 1, size(this % cell_id)
      atoms(:, :) = ZERO
      total_hits = 0

      if (master) then
#ifdef MPI
        do j = 1, n_procs - 1
          call MPI_RECV(n, 1, MPI_INTEGER, j, 0, MPI_COMM_WORLD, &
               MPI_STATUS_IGNORE, mpi_err)

          allocate(data(2*n))
          call MPI_RECV(data, 2*n, MPI_INTEGER, j, 1, MPI_COMM_WORLD, &
               MPI_STATUS_IGNORE, mpi_err)
          do k = 0, n - 1
            do m = 1, master_indices(i_cell) % size()
              if (data(2*k + 1) == master_indices(i_cell) % data(m)) then
                master_hits(i_cell) % data(m) = master_hits(i_cell) % data(m) + &
                     data(2*k + 2)
              end if
            end do
          end do
          deallocate(data)
        end do
#endif

        do j = 1, master_indices(i_cell) % size()
          total_hits = total_hits + master_hits(i_cell) % data(j)
          f = real(master_hits(i_cell) % data(j), 8) / this % samples
          var_f = f*(ONE - f) / this % samples

          i_material = master_indices(i_cell) % data(j)
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
        volume(1, i_cell) = real(total_hits, 8) / this % samples * volume_sample
        volume(2, i_cell) = sqrt(volume(1, i_cell) * (volume_sample - &
             volume(1, i_cell)) / this % samples)

        ! Determine total number of atoms. At this point, we have values in
        ! atoms/b-cm. To get to atoms we multiple by 10^24 V.
        do j = 1, size(atoms, 2)
          atoms(1, j) = 1.0e24_8 * volume_sample * atoms(1, j)
          atoms(2, j) = 1.0e24_8 * volume_sample * sqrt(atoms(2, j))
        end do

        ! Convert full arrays to vectors
        do j = 1, size(nuclides)
          if (atoms(1, j) > ZERO) then
            call nuclide_vec(i_cell) % push_back(j)
            call atoms_vec(i_cell) % push_back(atoms(1, j))
            call uncertainty_vec(i_cell) % push_back(atoms(2, j))
          end if
        end do

      else
#ifdef MPI
        n = master_indices(i_cell) % size()
        allocate(data(2*n))
        do k = 0, n - 1
          data(2*k + 1) = master_indices(i_cell) % data(k + 1)
          data(2*k + 2) = master_hits(i_cell) % data(k + 1)
        end do

        call MPI_SEND(n, 1, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, mpi_err)
        call MPI_SEND(data, 2*n, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, mpi_err)
        deallocate(data)
#endif
      end if
    end do

  end subroutine get_volume

!===============================================================================
! WRITE_VOLUME writes the results of a single stochastic volume calculation to
! an HDF5 file
!===============================================================================

  subroutine write_volume(this, filename, volume, nuclide_vec, atoms_vec, &
       uncertainty_vec)
    type(VolumeCalculation), intent(in) :: this
    character(*),     intent(in) :: filename       ! filename for HDF5 file
    real(8),          intent(in) :: volume(:,:)    ! volume mean/stdev in each cell
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

    ! Write basic metadata
    call write_attribute(file_id, "samples", this % samples)
    call write_attribute(file_id, "lower_left", this % lower_left)
    call write_attribute(file_id, "upper_right", this % upper_right)

    do i = 1, size(this % cell_id)
      group_id = create_group(file_id, "cell_" // trim(to_str(&
           this % cell_id(i))))

      ! Write volume for cell
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
