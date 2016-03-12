module eigenvalue

#ifdef MPI
  use message_passing
#endif

  use constants,   only: ZERO
  use error,       only: fatal_error, warning
  use global
  use math,        only: t_percentile
  use mesh,        only: count_bank_sites
  use mesh_header, only: RegularMesh
  use random_lcg,  only: prn, set_particle_seed, advance_prn_seed
  use search,      only: binary_search
  use string,      only: to_str

  implicit none

  real(8) :: keff_generation ! Single-generation k on each processor
  real(8) :: k_sum(2) = ZERO ! Used to reduce sum and sum_sq

contains

!===============================================================================
! SYNCHRONIZE_BANK samples source sites from the fission sites that were
! accumulated during the generation. This routine is what allows this Monte
! Carlo to scale to large numbers of processors where other codes cannot.
!===============================================================================

  subroutine synchronize_bank()

    integer    :: i            ! loop indices
    integer    :: j            ! loop indices
    integer(8) :: start        ! starting index in global bank
    integer(8) :: finish       ! ending index in global bank
    integer(8) :: total        ! total sites in global fission bank
    integer(8) :: index_temp   ! index in temporary source bank
    integer(8) :: sites_needed ! # of sites to be sampled
    real(8)    :: p_sample     ! probability of sampling a site
    type(Bank), save, allocatable :: &
         & temp_sites(:)       ! local array of extra sites on each node

#ifdef MPI
    integer(8) :: n            ! number of sites to send/recv
    integer    :: neighbor     ! processor to send/recv data from
#ifdef MPIF08
    type(MPI_Request) :: request(20)
#else
    integer    :: request(20)  ! communication request for send/recving sites
#endif
    integer    :: n_request    ! number of communication requests
    integer(8) :: index_local  ! index in local source bank
    integer(8), save, allocatable :: &
         & bank_position(:)    ! starting positions in global source bank
#endif

    ! In order to properly understand the fission bank algorithm, you need to
    ! think of the fission and source bank as being one global array divided
    ! over multiple processors. At the start, each processor has a random amount
    ! of fission bank sites -- each processor needs to know the total number of
    ! sites in order to figure out the probability for selecting
    ! sites. Furthermore, each proc also needs to know where in the 'global'
    ! fission bank its own sites starts in order to ensure reproducibility by
    ! skipping ahead to the proper seed.

#ifdef MPI
    start = 0_8
    call MPI_EXSCAN(n_bank, start, 1, MPI_INTEGER8, MPI_SUM, &
         MPI_COMM_WORLD, mpi_err)

    ! While we would expect the value of start on rank 0 to be 0, the MPI
    ! standard says that the receive buffer on rank 0 is undefined and not
    ! significant
    if (rank == 0) start = 0_8

    finish = start + n_bank
    total = finish
    call MPI_BCAST(total, 1, MPI_INTEGER8, n_procs - 1, &
         MPI_COMM_WORLD, mpi_err)

#else
    start  = 0_8
    finish = n_bank
    total  = n_bank
#endif

    ! If there are not that many particles per generation, it's possible that no
    ! fission sites were created at all on a single processor. Rather than add
    ! extra logic to treat this circumstance, we really want to ensure the user
    ! runs enough particles to avoid this in the first place.

    if (n_bank == 0) then
      call fatal_error("No fission sites banked on processor " // to_str(rank))
    end if

    ! Make sure all processors start at the same point for random sampling. Then
    ! skip ahead in the sequence using the starting index in the 'global'
    ! fission bank for each processor.

    call set_particle_seed(int((current_batch - 1)*gen_per_batch + &
         current_gen,8))
    call advance_prn_seed(start)

    ! Determine how many fission sites we need to sample from the source bank
    ! and the probability for selecting a site.

    if (total < n_particles) then
      sites_needed = mod(n_particles,total)
    else
      sites_needed = n_particles
    end if
    p_sample = real(sites_needed,8)/real(total,8)

    call time_bank_sample % start()

    ! ==========================================================================
    ! SAMPLE N_PARTICLES FROM FISSION BANK AND PLACE IN TEMP_SITES

    ! Allocate temporary source bank
    index_temp = 0_8
    if (.not. allocated(temp_sites)) allocate(temp_sites(3*work))

    do i = 1, int(n_bank,4)

      ! If there are less than n_particles particles banked, automatically add
      ! int(n_particles/total) sites to temp_sites. For example, if you need
      ! 1000 and 300 were banked, this would add 3 source sites per banked site
      ! and the remaining 100 would be randomly sampled.
      if (total < n_particles) then
        do j = 1, int(n_particles/total)
          index_temp = index_temp + 1
          temp_sites(index_temp) = fission_bank(i)
        end do
      end if

      ! Randomly sample sites needed
      if (prn() < p_sample) then
        index_temp = index_temp + 1
        temp_sites(index_temp) = fission_bank(i)
      end if
    end do

    ! At this point, the sampling of source sites is done and now we need to
    ! figure out where to send source sites. Since it is possible that one
    ! processor's share of the source bank spans more than just the immediate
    ! neighboring processors, we have to perform an ALLGATHER to determine the
    ! indices for all processors

#ifdef MPI
    ! First do an exclusive scan to get the starting indices for
    start = 0_8
    call MPI_EXSCAN(index_temp, start, 1, MPI_INTEGER8, MPI_SUM, &
         MPI_COMM_WORLD, mpi_err)
    finish = start + index_temp

    ! Allocate space for bank_position if this hasn't been done yet
    if (.not. allocated(bank_position)) allocate(bank_position(n_procs))
    call MPI_ALLGATHER(start, 1, MPI_INTEGER8, bank_position, 1, &
         MPI_INTEGER8, MPI_COMM_WORLD, mpi_err)
#else
    start  = 0_8
    finish = index_temp
#endif

    ! Now that the sampling is complete, we need to ensure that we have exactly
    ! n_particles source sites. The way this is done in a reproducible manner is
    ! to adjust only the source sites on the last processor.

    if (rank == n_procs - 1) then
      if (finish > n_particles) then
        ! If we have extra sites sampled, we will simply discard the extra
        ! ones on the last processor
        index_temp = n_particles - start

      elseif (finish < n_particles) then
        ! If we have too few sites, repeat sites from the very end of the
        ! fission bank
        sites_needed = n_particles - finish
        do i = 1, int(sites_needed,4)
          index_temp = index_temp + 1
          temp_sites(index_temp) = fission_bank(n_bank - sites_needed + i)
        end do
      end if

      ! the last processor should not be sending sites to right
      finish = work_index(rank + 1)
    end if

    call time_bank_sample % stop()
    call time_bank_sendrecv % start()

#ifdef MPI
    ! ==========================================================================
    ! SEND BANK SITES TO NEIGHBORS

    index_local = 1
    n_request = 0

    if (start < n_particles) then
      ! Determine the index of the processor which has the first part of the
      ! source_bank for the local processor
      neighbor = binary_search(work_index, n_procs + 1, start) - 1

      SEND_SITES: do while (start < finish)
        ! Determine the number of sites to send
        n = min(work_index(neighbor + 1), finish) - start

        ! Initiate an asynchronous send of source sites to the neighboring
        ! process
        if (neighbor /= rank) then
          n_request = n_request + 1
          call MPI_ISEND(temp_sites(index_local), int(n), MPI_BANK, neighbor, &
               rank, MPI_COMM_WORLD, request(n_request), mpi_err)
        end if

        ! Increment all indices
        start       = start       + n
        index_local = index_local + n
        neighbor    = neighbor    + 1

        ! Check for sites out of bounds -- this only happens in the rare
        ! circumstance that a processor close to the end has so many sites that
        ! it would exceed the bank on the last processor
        if (neighbor > n_procs - 1) exit
      end do SEND_SITES
    end if

    ! ==========================================================================
    ! RECEIVE BANK SITES FROM NEIGHBORS OR TEMPORARY BANK

    start = work_index(rank)
    index_local = 1

    ! Determine what process has the source sites that will need to be stored at
    ! the beginning of this processor's source bank.

    if (start >= bank_position(n_procs)) then
      neighbor = n_procs - 1
    else
      neighbor = binary_search(bank_position, n_procs, start) - 1
    end if

    RECV_SITES: do while (start < work_index(rank + 1))
      ! Determine how many sites need to be received
      if (neighbor == n_procs - 1) then
        n = work_index(rank + 1) - start
      else
        n = min(bank_position(neighbor + 2), work_index(rank + 1)) - start
      end if

      if (neighbor /= rank) then
        ! If the source sites are not on this processor, initiate an
        ! asynchronous receive for the source sites

        n_request = n_request + 1
        call MPI_IRECV(source_bank(index_local), int(n), MPI_BANK, &
             neighbor, neighbor, MPI_COMM_WORLD, request(n_request), mpi_err)

      else
        ! If the source sites are on this procesor, we can simply copy them
        ! from the temp_sites bank

        index_temp = start - bank_position(rank+1) + 1
        source_bank(index_local:index_local+n-1) = &
             temp_sites(index_temp:index_temp+n-1)
      end if

      ! Increment all indices
      start       = start       + n
      index_local = index_local + n
      neighbor    = neighbor    + 1
    end do RECV_SITES

    ! Since we initiated a series of asynchronous ISENDs and IRECVs, now we have
    ! to ensure that the data has actually been communicated before moving on to
    ! the next generation

    call MPI_WAITALL(n_request, request, MPI_STATUSES_IGNORE, mpi_err)

    ! Deallocate space for bank_position on the very last generation
    if (current_batch == n_max_batches .and. current_gen == gen_per_batch) &
         deallocate(bank_position)
#else
    source_bank = temp_sites(1:n_particles)
#endif

    call time_bank_sendrecv % stop()

    ! Deallocate space for the temporary source bank on the last generation
    if (current_batch == n_max_batches .and. current_gen == gen_per_batch) &
         deallocate(temp_sites)

  end subroutine synchronize_bank

!===============================================================================
! SHANNON_ENTROPY calculates the Shannon entropy of the fission source
! distribution to assess source convergence
!===============================================================================

  subroutine shannon_entropy()

    integer :: ent_idx        ! entropy index
    integer :: i, j, k        ! index for bank sites
    integer :: n              ! # of boxes in each dimension
    logical :: sites_outside  ! were there sites outside entropy box?
    type(RegularMesh), pointer :: m

    ! Get pointer to entropy mesh
    m => entropy_mesh

    ! On the first pass through this subroutine, we need to determine how big
    ! the entropy mesh should be in each direction and then allocate a
    ! three-dimensional array to store the fraction of source sites in each mesh
    ! box

    if (.not. allocated(entropy_p)) then
      if (.not. allocated(m % dimension)) then
        ! If the user did not specify how many mesh cells are to be used in
        ! each direction, we automatically determine an appropriate number of
        ! cells
        n = ceiling((n_particles/20)**(ONE/THREE))

        ! copy dimensions
        m % n_dimension = 3
        allocate(m % dimension(3))
        m % dimension = n

        ! determine width
        m % width = (m % upper_right - m % lower_left) / m % dimension

      end if

      ! allocate p
      allocate(entropy_p(1, m % dimension(1), m % dimension(2), &
           m % dimension(3)))
    end if

    ! count number of fission sites over mesh
    call count_bank_sites(m, fission_bank, entropy_p, &
         size_bank=n_bank, sites_outside=sites_outside)

    ! display warning message if there were sites outside entropy box
    if (sites_outside) then
      if (master) call warning("Fission source site(s) outside of entropy box.")
    end if

    ! sum values to obtain shannon entropy
    if (master) then
      ! Normalize to total weight of bank sites
      entropy_p = entropy_p / sum(entropy_p)

      ent_idx = current_gen + gen_per_batch*(current_batch - 1)
      entropy(ent_idx) = ZERO
      do i = 1, m % dimension(1)
        do j = 1, m % dimension(2)
          do k = 1, m % dimension(3)
            if (entropy_p(1,i,j,k) > ZERO) then
              entropy(ent_idx) = entropy(ent_idx) - &
                   entropy_p(1,i,j,k) * log(entropy_p(1,i,j,k))/log(TWO)
            end if
          end do
        end do
      end do
    end if

  end subroutine shannon_entropy

!===============================================================================
! CALCULATE_GENERATION_KEFF collects the single-processor tracklength k's onto
! the master processor and normalizes them. This should work whether or not the
! no-reduce method is being used.
!===============================================================================

  subroutine calculate_generation_keff()

    ! Get keff for this generation by subtracting off the starting value
    keff_generation = global_tallies(K_TRACKLENGTH) % value - keff_generation

#ifdef MPI
    ! Combine values across all processors
    call MPI_ALLREDUCE(keff_generation, k_generation(overall_gen), 1, &
         MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpi_err)
#else
    k_generation(overall_gen) = keff_generation
#endif

    ! Normalize single batch estimate of k
    ! TODO: This should be normalized by total_weight, not by n_particles
    k_generation(overall_gen) = k_generation(overall_gen) / n_particles

  end subroutine calculate_generation_keff

!===============================================================================
! CALCULATE_AVERAGE_KEFF calculates the mean and standard deviation of the mean
! of k-effective during active generations and broadcasts the mean to all
! processors
!===============================================================================

  subroutine calculate_average_keff()

    integer :: n        ! number of active generations
    real(8) :: alpha    ! significance level for CI
    real(8) :: t_value  ! t-value for confidence intervals

    ! Determine number of active generations
    n = overall_gen - n_inactive*gen_per_batch

    if (n <= 0) then
      ! For inactive generations, use current generation k as estimate for next
      ! generation
      keff = k_generation(overall_gen)

    else
      ! Sample mean of keff
      k_sum(1) = k_sum(1) + k_generation(overall_gen)
      k_sum(2) = k_sum(2) + k_generation(overall_gen)**2

      ! Determine mean
      keff = k_sum(1) / n

      if (n > 1) then
        if (confidence_intervals) then
          ! Calculate t-value for confidence intervals
          alpha = ONE - CONFIDENCE_LEVEL
          t_value = t_percentile(ONE - alpha/TWO, n - 1)
        else
          t_value = ONE
        end if

        ! Standard deviation of the sample mean of k
        keff_std = t_value * sqrt((k_sum(2)/n - keff**2) / (n - 1))
      end if
    end if

  end subroutine calculate_average_keff

!===============================================================================
! CALCULATE_COMBINED_KEFF calculates a minimum variance estimate of k-effective
! based on a linear combination of the collision, absorption, and tracklength
! estimates. The theory behind this can be found in M. Halperin, "Almost
! linearly-optimum combination of unbiased estimates," J. Am. Stat. Assoc., 56,
! 36-43 (1961), doi:10.1080/01621459.1961.10482088. The implementation here
! follows that described in T. Urbatsch et al., "Estimation and interpretation
! of keff confidence intervals in MCNP," Nucl. Technol., 111, 169-182 (1995).
!===============================================================================

  subroutine calculate_combined_keff()

    integer :: l        ! loop index
    integer :: i, j, k  ! indices referring to collision, absorption, or track
    integer :: n        ! number of realizations
    real(8) :: kv(3)    ! vector of k-effective estimates
    real(8) :: cov(3,3) ! sample covariance matrix
    real(8) :: f        ! weighting factor
    real(8) :: g        ! sum of weighting factors
    real(8) :: S(3)     ! sums used for variance calculation

    ! Make sure we have at least four realizations. Notice that at the end,
    ! there is a N-3 term in a denominator.
    if (n_realizations <= 3) return

    ! Initialize variables
    n = n_realizations
    g = ZERO
    S = ZERO
    k_combined = ZERO

    ! Copy estimates of k-effective and its variance (not variance of the mean)
    kv(1) = global_tallies(K_COLLISION) % sum / n
    kv(2) = global_tallies(K_ABSORPTION) % sum / n
    kv(3) = global_tallies(K_TRACKLENGTH) % sum / n
    cov(1,1) = (global_tallies(K_COLLISION) % sum_sq - &
         n * kv(1) * kv(1)) / (n - 1)
    cov(2,2) = (global_tallies(K_ABSORPTION) % sum_sq - &
         n * kv(2) * kv(2)) / (n - 1)
    cov(3,3) = (global_tallies(K_TRACKLENGTH) % sum_sq - &
         n * kv(3) * kv(3)) / (n - 1)

    ! Calculate covariances based on sums with Bessel's correction
    cov(1,2) = (k_col_abs - n * kv(1) * kv(2))/(n - 1)
    cov(1,3) = (k_col_tra - n * kv(1) * kv(3))/(n - 1)
    cov(2,3) = (k_abs_tra - n * kv(2) * kv(3))/(n - 1)
    cov(2,1) = cov(1,2)
    cov(3,1) = cov(1,3)
    cov(3,2) = cov(2,3)

    do l = 1, 3
      ! Permutations of estimates
      if (l == 1) then
        ! i = collision, j = absorption, k = tracklength
        i = 1
        j = 2
        k = 3
      elseif (l == 2) then
        ! i = absortion, j = tracklength, k = collision
        i = 2
        j = 3
        k = 1
      elseif (l == 3) then
        ! i = tracklength, j = collision, k = absorption
        i = 3
        j = 1
        k = 2
      end if

      ! Calculate weighting
      f = cov(j,j)*(cov(k,k) - cov(i,k)) - cov(k,k)*cov(i,j) + &
           cov(j,k)*(cov(i,j) + cov(i,k) - cov(j,k))

      ! Add to S sums for variance of combined estimate
      S(1) = S(1) + f * cov(1,l)
      S(2) = S(2) + (cov(j,j) + cov(k,k) - TWO*cov(j,k))*kv(l)*kv(l)
      S(3) = S(3) + (cov(k,k) + cov(i,j) - cov(j,k) - cov(i,k))*kv(l)*kv(j)

      ! Add to sum for combined k-effective
      k_combined(1) = k_combined(1) + f * kv(l)
      g = g + f
    end do

    ! Complete calculations of S sums
    S = (n - 1)*S
    S(1) = (n - 1)**2 * S(1)

    ! Calculate combined estimate of k-effective
    k_combined(1) = k_combined(1) / g

    ! Calculate standard deviation of combined estimate
    g = (n - 1)**2 * g
    k_combined(2) = sqrt(S(1)/(g*n*(n-3)) * (ONE + n*((S(2) - TWO*S(3))/g)))

  end subroutine calculate_combined_keff

!===============================================================================
! COUNT_SOURCE_FOR_UFS determines the source fraction in each UFS mesh cell and
! reweights the source bank so that the sum of the weights is equal to
! n_particles. The 'source_frac' variable is used later to bias the production
! of fission sites
!===============================================================================

  subroutine count_source_for_ufs()

    real(8) :: total         ! total weight in source bank
    logical :: sites_outside ! were there sites outside the ufs mesh?
#ifdef MPI
    integer :: n             ! total number of ufs mesh cells
#endif

    if (current_batch == 1 .and. current_gen == 1) then
      ! On the first generation, just assume that the source is already evenly
      ! distributed so that effectively the production of fission sites is not
      ! biased

      source_frac = ufs_mesh % volume_frac

    else
      ! count number of source sites in each ufs mesh cell
      call count_bank_sites(ufs_mesh, source_bank, source_frac, &
           sites_outside=sites_outside, size_bank=work)

      ! Check for sites outside of the mesh
      if (master .and. sites_outside) then
        call fatal_error("Source sites outside of the UFS mesh!")
      end if

#ifdef MPI
      ! Send source fraction to all processors
      n = product(ufs_mesh % dimension)
      call MPI_BCAST(source_frac, n, MPI_REAL8, 0, MPI_COMM_WORLD, mpi_err)
#endif

      ! Normalize to total weight to get fraction of source in each cell
      total = sum(source_frac)
      source_frac = source_frac / total

      ! Since the total starting weight is not equal to n_particles, we need to
      ! renormalize the weight of the source sites

      source_bank % wgt = source_bank % wgt * n_particles / total
    end if

  end subroutine count_source_for_ufs

#ifdef _OPENMP
!===============================================================================
! JOIN_BANK_FROM_THREADS joins threadprivate fission banks into a single fission
! bank that can be sampled. Note that this operation is necessarily sequential
! to preserve the order of the bank when using varying numbers of threads.
!===============================================================================

  subroutine join_bank_from_threads()

    integer(8) :: total ! total number of fission bank sites
    integer    :: i     ! loop index for threads

    ! Initialize the total number of fission bank sites
    total = 0

!$omp parallel

    ! Copy thread fission bank sites to one shared copy
!$omp do ordered schedule(static)
    do i = 1, n_threads
!$omp ordered
      master_fission_bank(total+1:total+n_bank) = fission_bank(1:n_bank)
      total = total + n_bank
!$omp end ordered
    end do
!$omp end do

    ! Make sure all threads have made it to this point
!$omp barrier

    ! Now copy the shared fission bank sites back to the master thread's copy.
    if (thread_id == 0) then
      n_bank = total
      fission_bank(1:n_bank) = master_fission_bank(1:n_bank)
    else
      n_bank = 0
    end if

!$omp end parallel

  end subroutine join_bank_from_threads
#endif

end module eigenvalue
