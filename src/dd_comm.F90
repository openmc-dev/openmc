module dd_comm
  
  use bank_header,      only: Bank
  use constants
  use dd_header,        only: dd_type
  use error,            only: fatal_error
  use global,           only: domain_decomp, MPI_BANK, MPI_PARTICLEBUFFER, n_procs, &
                              rank, mpi_err, work, message, &
                              n_particles, work_index, &
                              source_bank, fission_bank, size_source_bank, &
                              size_fission_bank, &
                              time_dd_sendrecv, time_dd_info, &
                              current_batch, current_stage
  use mesh,             only: get_mesh_bin
  use output,           only: write_message
  use particle_header,  only: buffer_to_particle, particle_to_buffer
  use resize_arr,       only: resize_array, extend_array
  use search,           only: binary_search
  use string,           only: to_str
  
  use mpi
  
  implicit none
  public

contains

!===============================================================================
! DISTRIBUTE_SOURCE loops through all particles in the source bank and sends
! them to the process that is handling the domain they are located on.  This
! routine also receives particles from other processes, and saves them in the
! source bank, overwriting the previous values.
!===============================================================================

  subroutine distribute_source()
  
    integer(8) :: i           ! loop index over initially-sampled source sites
    integer :: n_source_sites ! number of sites that will start in this domain
    integer :: to_domain
    integer :: to_rank        ! global rank of the processor to send to
    integer :: alloc_err      ! allocation error code
    integer :: n_request      ! number of communication requests
    integer :: new_source     ! size of the source_bank after distribution
    integer, allocatable    :: request(:) ! communication requests
    integer, allocatable    :: n_send_domain(:) ! number sent to each domain
    integer, allocatable    :: n_send_rank(:) ! number sent to each process
    type(dd_type), pointer  :: dd => domain_decomp
    type(Bank), allocatable :: buffer(:)

    ! Since remote domains could have more than one process working on them, we
    ! need to keep track of how many particles we've sent to each so that we can
    ! rotate through processes to send to
    allocate(n_send_domain(n_procs))
    n_send_domain = 0

    ! We're going to send sites one at a time, then ALLREDUCE the number of
    ! sent sites so each rank knows how many revcs to post
    allocate(n_send_rank(n_procs))
    n_send_rank = 0

    ! We need a copy of the send data
    allocate(buffer(work), STAT=alloc_err)
    
    ! Check for allocation errors 
    if (alloc_err /= 0) then
      message = "Failed to allocate source bank send buffer during domain &
                &decomposition initialization."
      call fatal_error()
    end if
    
    ! TODO: This may be too conservative an allocation for most use-cases
    allocate(request(work + n_particles), STAT=alloc_err)
    n_request = 0

    ! Check for allocation errors 
    if (alloc_err /= 0) then
      message = "Failed to allocate communication request array during domain &
                &decomposition initialization."
      call fatal_error()
    end if

    ! Copy initially-sampled source sites
    buffer = source_bank
    
    ! We'll save sites that belong in this domain to the source bank as we
    ! loop through and find them, as well as when we receive them from other
    ! processes
    n_source_sites = 0

    ! Loop over all sampled source sites and send/save sites as needed
    do i = 1, work
    
      ! Determine which domain meshbin this should start on
      call get_mesh_bin(dd % mesh, buffer(i) % xyz, to_domain)
    
      if (to_domain == dd % meshbin) then
        ! This site should start on this domain
        
        n_source_sites = n_source_sites + 1
        source_bank(n_source_sites) = buffer(i)
        
      else
        ! This site needs to be sent
        
        to_rank = dd % domain_masters(to_domain) + &
            mod(n_send_domain(to_domain), dd % domain_n_procs(to_domain))
        
!        print *, rank, work_index(rank)+i, dd % meshbin, "sending", to_domain, n_send_domain(to_domain), to_rank, &
!            dd % domain_n_procs(to_domain)
        
        ! Increment send counters
        n_send_rank(to_rank + 1) = n_send_rank(to_rank + 1) + 1
        n_send_domain(to_domain) = n_send_domain(to_domain) + 1

        ! Post the send
        n_request = n_request + 1
        call MPI_ISEND(buffer(i), 1, MPI_BANK, to_rank, rank, MPI_COMM_WORLD, &
            request(n_request), mpi_err)

      end if
    
    end do

    ! Determine how many sites were sent here by reducing n_send_rank
    call MPI_ALLREDUCE(MPI_IN_PLACE, n_send_rank, n_procs, MPI_INTEGER, &
        MPI_SUM, MPI_COMM_WORLD, mpi_err)

    ! Check if we're not going to have enough space in the source_bank (which
    ! we're using as the receive buffer), and resize the array if needed
    new_source = n_send_rank(rank + 1) + n_source_sites
    if (new_source > size_source_bank) &
        call extend_array(source_bank, size_source_bank, int(new_source, 8), &
            alloc_err)
    

    ! Receive sites from other processes
    do i = 1, n_send_rank(rank + 1)
    
      n_source_sites = n_source_sites + 1
      n_request = n_request + 1
      call MPI_IRECV(source_bank(n_source_sites), 1, MPI_BANK, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &
                 request(n_request), mpi_err)
    end do

    ! Wait for all sends/recvs to complete
    call MPI_WAITALL(n_request, request, MPI_STATUSES_IGNORE, mpi_err)

    ! Reset how much work needs to be done on this process
    work = n_source_sites
    
    ! Free memory
    deallocate(n_send_domain)
    deallocate(n_send_rank)
    deallocate(buffer)
  
  end subroutine distribute_source

!===============================================================================
! SYNCHRONIZE_BANK_DD does the same thing as the non-domain-decomposed
! synchronize_bank routine in the eigenvalue module, except it does it for the
! local domain communicator group only. All fission sites are sampled using
! stored prn_seeds, maintaining tally reproducibility. Special handling was
! introduced for cases where no fission particles are banked.
!===============================================================================

  subroutine synchronize_bank_dd()

!    integer    :: i            ! loop indices
!    integer    :: j            ! loop indices
!    integer(8) :: start        ! starting index in global bank
!    integer(8) :: finish       ! ending index in global bank
!    integer(8) :: total        ! total sites in global fission bank
!    integer(8) :: index_temp   ! index in temporary source bank
!    integer(8) :: sites_needed ! # of sites to be sampled
!    real(8)    :: p_sample     ! probability of sampling a site
!    type(Bank), save, allocatable :: &
!         & temp_sites(:)       ! local array of extra sites on each node

!#ifdef MPI
!    integer(8) :: n            ! number of sites to send/recv
!    integer    :: neighbor     ! processor to send/recv data from
!    integer    :: request(20)  ! communication request for send/recving sites
!    integer    :: n_request    ! number of communication requests
!    integer(8) :: index_local  ! index in local source bank
!    integer(8), save, allocatable :: &
!         & bank_position(:)    ! starting positions in global source bank
!#endif

!    ! In order to properly understand the fission bank algorithm, you need to
!    ! think of the fission and source bank as being one global array divided
!    ! over multiple processors. At the start, each processor has a random amount
!    ! of fission bank sites -- each processor needs to know the total number of
!    ! sites in order to figure out the probability for selecting
!    ! sites. Furthermore, each proc also needs to know where in the 'global'
!    ! fission bank its own sites starts in order to ensure reproducibility by
!    ! skipping ahead to the proper seed.

!#ifdef MPI
!    start = 0_8
!    call MPI_EXSCAN(n_bank, start, 1, MPI_INTEGER8, MPI_SUM, & 
!         MPI_COMM_WORLD, mpi_err)

!    ! While we would expect the value of start on rank 0 to be 0, the MPI
!    ! standard says that the receive buffer on rank 0 is undefined and not
!    ! significant
!    if (rank == 0) start = 0_8

!    finish = start + n_bank
!    total = finish
!    call MPI_BCAST(total, 1, MPI_INTEGER8, n_procs - 1, & 
!         MPI_COMM_WORLD, mpi_err)

!#else
!    start  = 0_8
!    finish = n_bank
!    total  = n_bank
!#endif

!    ! If there are not that many particles per generation, it's possible that no
!    ! fission sites were created at all on a single processor. Rather than add
!    ! extra logic to treat this circumstance, we really want to ensure the user
!    ! runs enough particles to avoid this in the first place.

!    if (n_bank == 0) then
!      message = "No fission sites banked on processor " // to_str(rank)
!      call fatal_error()
!    end if

!    ! Make sure all processors start at the same point for random sampling. Then
!    ! skip ahead in the sequence using the starting index in the 'global'
!    ! fission bank for each processor.

!    call set_particle_seed(int((current_batch - 1)*gen_per_batch + &
!         current_gen,8))
!    call prn_skip(start)

!    ! Determine how many fission sites we need to sample from the source bank
!    ! and the probability for selecting a site.

!    if (total < n_particles) then
!      sites_needed = mod(n_particles,total)
!    else
!      sites_needed = n_particles
!    end if
!    p_sample = real(sites_needed,8)/real(total,8)

!    call time_bank_sample % start()

!    ! ==========================================================================
!    ! SAMPLE N_PARTICLES FROM FISSION BANK AND PLACE IN TEMP_SITES

!    ! Allocate temporary source bank
!    index_temp = 0_8
!    if (.not. allocated(temp_sites)) allocate(temp_sites(3*work))

!    do i = 1, int(n_bank,4)

!      ! If there are less than n_particles particles banked, automatically add
!      ! int(n_particles/total) sites to temp_sites. For example, if you need
!      ! 1000 and 300 were banked, this would add 3 source sites per banked site
!      ! and the remaining 100 would be randomly sampled.
!      if (total < n_particles) then
!        do j = 1, int(n_particles/total)
!          index_temp = index_temp + 1
!          temp_sites(index_temp) = fission_bank(i)
!        end do
!      end if

!      ! Randomly sample sites needed
!      if (prn() < p_sample) then
!        index_temp = index_temp + 1
!        temp_sites(index_temp) = fission_bank(i)
!      end if
!    end do

!    ! At this point, the sampling of source sites is done and now we need to
!    ! figure out where to send source sites. Since it is possible that one
!    ! processor's share of the source bank spans more than just the immediate
!    ! neighboring processors, we have to perform an ALLGATHER to determine the
!    ! indices for all processors

!#ifdef MPI
!    ! First do an exclusive scan to get the starting indices for 
!    start = 0_8
!    call MPI_EXSCAN(index_temp, start, 1, MPI_INTEGER8, MPI_SUM, & 
!         MPI_COMM_WORLD, mpi_err)
!    finish = start + index_temp

!    ! Allocate space for bank_position if this hasn't been done yet
!    if (.not. allocated(bank_position)) allocate(bank_position(n_procs))
!    call MPI_ALLGATHER(start, 1, MPI_INTEGER8, bank_position, 1, &
!         MPI_INTEGER8, MPI_COMM_WORLD, mpi_err)
!#else
!    start  = 0_8
!    finish = index_temp
!#endif

!    ! Now that the sampling is complete, we need to ensure that we have exactly
!    ! n_particles source sites. The way this is done in a reproducible manner is
!    ! to adjust only the source sites on the last processor.

!    if (rank == n_procs - 1) then
!      if (finish > n_particles) then
!        ! If we have extra sites sampled, we will simply discard the extra
!        ! ones on the last processor
!        index_temp = n_particles - start

!      elseif (finish < n_particles) then
!        ! If we have too few sites, repeat sites from the very end of the
!        ! fission bank
!        sites_needed = n_particles - finish
!        do i = 1, int(sites_needed,4)
!          index_temp = index_temp + 1
!          temp_sites(index_temp) = fission_bank(n_bank - sites_needed + i)
!        end do
!      end if

!      ! the last processor should not be sending sites to right
!      finish = work_index(rank + 1)
!    end if

!    call time_bank_sample % stop()
!    call time_bank_sendrecv % start()

!#ifdef MPI
!    ! ==========================================================================
!    ! SEND BANK SITES TO NEIGHBORS

!    index_local = 1
!    n_request = 0

!    if (start < n_particles) then
!      ! Determine the index of the processor which has the first part of the
!      ! source_bank for the local processor
!      neighbor = binary_search(work_index, n_procs + 1, start) - 1

!      SEND_SITES: do while (start < finish)
!        ! Determine the number of sites to send
!        n = min(work_index(neighbor + 1), finish) - start

!        ! Initiate an asynchronous send of source sites to the neighboring
!        ! process
!        if (neighbor /= rank) then
!          n_request = n_request + 1
!          call MPI_ISEND(temp_sites(index_local), n, MPI_BANK, neighbor, &
!               rank, MPI_COMM_WORLD, request(n_request), mpi_err)
!        end if

!        ! Increment all indices
!        start       = start       + n
!        index_local = index_local + n
!        neighbor    = neighbor    + 1

!        ! Check for sites out of bounds -- this only happens in the rare
!        ! circumstance that a processor close to the end has so many sites that
!        ! it would exceed the bank on the last processor
!        if (neighbor > n_procs - 1) exit
!      end do SEND_SITES
!    end if

!    ! ==========================================================================
!    ! RECEIVE BANK SITES FROM NEIGHBORS OR TEMPORARY BANK

!    start = work_index(rank)
!    index_local = 1

!    ! Determine what process has the source sites that will need to be stored at
!    ! the beginning of this processor's source bank.

!    if (start >= bank_position(n_procs)) then
!      neighbor = n_procs - 1
!    else
!      neighbor = binary_search(bank_position, n_procs, start) - 1
!    end if

!    RECV_SITES: do while (start < work_index(rank + 1))
!      ! Determine how many sites need to be received
!      if (neighbor == n_procs - 1) then
!        n = work_index(rank + 1) - start
!      else
!        n = min(bank_position(neighbor + 2), work_index(rank + 1)) - start
!      end if

!      if (neighbor /= rank) then
!        ! If the source sites are not on this processor, initiate an
!        ! asynchronous receive for the source sites

!        n_request = n_request + 1
!        call MPI_IRECV(source_bank(index_local), n, MPI_BANK, &
!             neighbor, neighbor, MPI_COMM_WORLD, request(n_request), mpi_err)

!      else
!        ! If the source sites are on this procesor, we can simply copy them
!        ! from the temp_sites bank

!        index_temp = start - bank_position(rank+1) + 1
!        source_bank(index_local:index_local+n-1) = &
!             temp_sites(index_temp:index_temp+n-1)
!      end if

!      ! Increment all indices
!      start       = start       + n
!      index_local = index_local + n
!      neighbor    = neighbor    + 1
!    end do RECV_SITES

!    ! Since we initiated a series of asynchronous ISENDs and IRECVs, now we have
!    ! to ensure that the data has actually been communicated before moving on to
!    ! the next generation

!    call MPI_WAITALL(n_request, request, MPI_STATUSES_IGNORE, mpi_err)

!    ! Deallocate space for bank_position on the very last generation
!    if (current_batch == n_batches .and. current_gen == gen_per_batch) &
!         deallocate(bank_position)
!#else
!    source_bank = temp_sites(1:n_particles)
!#endif

!    call time_bank_sendrecv % stop()

!    ! Deallocate space for the temporary source bank on the last generation
!    if (current_batch == n_batches .and. current_gen == gen_per_batch) &
!         deallocate(temp_sites)

  end subroutine synchronize_bank_dd


!===============================================================================
! SYNCHRONIZE_SCATTER_INFO communicates all required scattering information to
! the local domain neighborhood
!===============================================================================

  subroutine synchronize_transfer_info(dd)
  
    type(dd_type), intent(inout) :: dd
  
  end subroutine synchronize_transfer_info


!===============================================================================
! SYNCHRONIZE_PARTICLES contains the logic to send particles that moved out of
! this domain to the processors that handle the receiving domains, as well as
! recieve particles from other domains.
!===============================================================================

  subroutine synchronize_particles(dd)
  
    type(dd_type), intent(inout) :: dd
  
    integer(8) :: i8
    integer :: i              ! loop index over particle buffer
    integer :: j
    integer :: pr             ! loop index over processors to send to
    integer :: nd             ! loop index over neighbor domains
    integer :: to_bin         ! domain bin to send particles to
    integer :: to_meshbin     ! domain meshbin to send particles to
    integer :: from_bin       ! domain bin to receive particles from
    integer :: from_meshbin   ! domain meshbin to receive particles from
    integer :: neighbor_meshbin
    integer :: to_ddrank      ! rank of the processor to send to inside its domain communicator group
    integer :: to_rank        ! global rank of the processor to send to
    integer :: from_rank      ! global rank of the processor to recv from
    integer :: neighbor_bin   ! domain meshbin of a neighbor to to_bin
    integer :: index_scatt    ! index in the 'global' scatter array
    integer :: index_pbuffer  ! index in the particle buffer
    integer :: domain_start
    integer :: request(10000)   ! communication request for send/recving particles TODO this should be allocated to an appropriate max size depending on the number of procs per domain
    integer :: n_request      ! number of communication requests
    integer :: n_recv
    integer :: n_send
    integer(8) :: start
    integer :: local_start
    integer :: pbuffer_start
    integer :: debug(4)
    integer :: pr_bin
    integer :: size_buff
    integer :: mod_
    integer :: alloc_err
    
    ! This subroutine operates in a manner similar to synchronize_bank for the
    ! fission bank, where the group of processors that operate on a certain
    ! domain should be thought of as sharing arrays ('global' to just this
    ! domain). We want scattered particles to be distributed evenly across all
    ! processors that operate on each domain, so to do that we figure out where
    ! the current particles to transmit are located in this 'global' scattering
    ! array.  The steps followed are:
    !   Determine for this domain how many particles are going to each other one
    !   Synchronize domain-to-domains scatter info in the 2nd order neighborhood
    !   Determine how many particles this process sends to each other process
    !   Synchronize the process-to-process send info
    !   Send/recv particles using the process-to-process send info
    
    ! This barrier isn't necessary since everyone will wait at the first
    ! allreduce, but this allows for a more honest assessment of the
    ! communication costs in the algorithm.  (i.e. if we don't have the barrier,
    ! profilers will show the allreduce as taking the majority of the runtime)
    call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)

    call time_dd_info % start()
    
    !===========================================================================
    ! SYNCHRONIZE SCATTER INFORMATION
    
    ! First find the global number of scatters, using dd_n_scatters_domain as
    ! a temp holder - the order doesn't matter for this global sum
    dd % n_global_scatters = 0
    call MPI_ALLREDUCE(dd % n_scatters_local, dd % n_scatters_domain, &
        N_CARTESIAN_NEIGHBORS, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, mpi_err)
    dd % n_global_scatters = sum(dd % n_scatters_domain)

    ! Now reduce to find how many particles are leaving this domain to each of
    ! the neighbors
    call MPI_ALLREDUCE(dd % n_scatters_local, dd % n_scatters_domain, &
        N_CARTESIAN_NEIGHBORS, MPI_INTEGER, MPI_SUM, dd % comm, mpi_err)

!    ! +++++++++++ DD DEBUG PRINT +++++++++++
!    if (dd % local_master) call write_leakage_info()
!    ! +++++++++++ DD DEBUG PRINT +++++++++++

    ! Communicate domain scatter information to all 2nd neighbors in 2nd order
    ! neighborhood. The first processor in each domain sends this information
    ! to all other processors in the neighborhood.
    
    ! Send scatter info
    n_request = 0
    if (dd % local_master) then
      do to_bin = 1, N_CARTESIAN_NEIGHBORS
        to_meshbin = dd % neighbor_meshbins(to_bin)
        if (to_meshbin == NO_BIN_FOUND) cycle
        
        ! Send to neighbors neighbors
        do nd = 1, N_CARTESIAN_NEIGHBORS
          neighbor_bin = to_bin*N_CARTESIAN_NEIGHBORS + nd
          
          neighbor_meshbin = dd % neighbor_meshbins(neighbor_bin)
          if (neighbor_meshbin == NO_BIN_FOUND) cycle
          if (neighbor_meshbin == dd % meshbin) cycle ! don't send to ourself
          
          ! We need to send this info to each process on the receiving domain
          do pr = 1, dd % domain_n_procs(neighbor_meshbin)
            to_rank = dd % domain_masters(neighbor_meshbin) + pr - 1
          
            n_request = n_request + 1
            call MPI_ISEND(dd % n_scatters_domain(to_bin), 1, MPI_INTEGER, &
                to_rank, to_meshbin, MPI_COMM_WORLD, request(n_request), &
                mpi_err)
          
          end do
            
        end do
      end do
    end if
    
    ! Receive scatter info
    dd % n_scatters_neighborhood = 0
    do to_bin = 1, N_CARTESIAN_NEIGHBORS
      to_meshbin = dd % neighbor_meshbins(to_bin)
      if (to_meshbin == NO_BIN_FOUND) cycle
        
      do nd = 1, N_CARTESIAN_NEIGHBORS
        from_bin = to_bin*N_CARTESIAN_NEIGHBORS + nd
        from_meshbin = dd % neighbor_meshbins(from_bin)
        if (from_meshbin == NO_BIN_FOUND) cycle
        
        if (from_meshbin == dd % meshbin) then
          dd % n_scatters_neighborhood(from_bin, to_bin) = &
              dd % n_scatters_domain(to_bin)
        else
          from_rank = dd % domain_masters(from_meshbin)
          n_request = n_request + 1
          call MPI_IRECV(dd % n_scatters_neighborhood(from_bin, to_bin), 1, &
              MPI_INTEGER, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &
              request(n_request), mpi_err)
         end if
                 
      end do
    end do
    
    ! Wait for the scatter information to be synchronized
    call MPI_WAITALL(n_request, request, MPI_STATUSES_IGNORE, mpi_err)
    n_request = 0

    !===========================================================================
    ! DETERMINE PARTICLE SEND INFORMATION
    
    ! First get processor offsets - we do this here to ensure that every process
    ! calls the exscan in the same order
    dd % scatter_offest = 0
    do to_bin = 1, N_CARTESIAN_NEIGHBORS
      to_meshbin = dd % neighbor_meshbins(to_bin)
      if (to_meshbin == NO_BIN_FOUND) cycle
      call MPI_EXSCAN(dd % n_scatters_local(to_bin), &
          dd % scatter_offest(to_bin), 1, MPI_INTEGER, MPI_SUM, dd % comm, &
          mpi_err)
      if (dd % local_master) dd % scatter_offest(to_bin) = 0
    end do
    
      ! Now figure out where to send the particles
    dd % send_rank_info = 0
    do to_bin = 1, N_CARTESIAN_NEIGHBORS
      ! loop over each of the neighbors we need to send to

      to_meshbin = dd % neighbor_meshbins(to_bin)
      if (to_meshbin == NO_BIN_FOUND) cycle

      if (dd % n_scatters_local(to_bin) == 0) cycle ! no outscat to this domain

      ! To determine which processor (or processors) to send to, and how many to
      ! send to each, we first need to figure out which parts of the 'global'
      ! inscatter array of the target domain are made up of particles on this
      ! processor
      
      domain_start = 0
      do nd = 1, N_CARTESIAN_NEIGHBORS
        neighbor_bin = to_bin*N_CARTESIAN_NEIGHBORS + nd
        neighbor_meshbin = dd % neighbor_meshbins(neighbor_bin)
        if (neighbor_meshbin == NO_BIN_FOUND) cycle
        if (neighbor_meshbin == dd % meshbin) exit
        domain_start = domain_start + &
            dd % n_scatters_neighborhood(neighbor_bin, to_bin)
      end do
      
      ! At this point domain_start is the start of the scatter array for this
      ! domain.  We need it for this processer within the domain, so we do a
      ! scan like what's done in synchronize_bank, but for the domain group
      ! communicator, and add the result to domain_start. (done previously)
      index_scatt = domain_start + dd % scatter_offest(to_bin)

      ! Now we determine which parts of the array are owned by the receiving
      ! processors

      n_recv = sum(dd % n_scatters_neighborhood(:, to_bin))
      mod_ = modulo(n_recv, dd % domain_n_procs(to_meshbin))
      do j = 1, dd % domain_n_procs(to_meshbin)
        dd % proc_finish(j) = j * n_recv/dd % domain_n_procs(to_meshbin)
        if (j-1 < mod_) then
          dd % proc_finish(j) = dd % proc_finish(j) + j
        else
          dd % proc_finish(j) = dd % proc_finish(j) + mod_
        end if
      end do

      ! Determine which receiving processor this starting index corresponds to
      if (index_scatt <= dd % proc_finish(1)) then
        to_ddrank = 0
      else
        to_ddrank = binary_search(dd % proc_finish, &
            dd % domain_n_procs(to_meshbin), index_scatt)
      end if

      ! Now we can loop through particles to send and build up send information
      pr_bin = (to_bin - 1) * dd % max_domain_procs + to_ddrank + 1
      do i = 1, dd % size_particle_buffer
        
        if (dd % particle_buffer(i) % outscatter_destination == to_bin) then

          ! Accumulate send info
          dd % send_rank_info(pr_bin) = dd % send_rank_info(pr_bin) + 1
          index_scatt = index_scatt + 1
          
          ! Increment the processor rank and bin if we're now on the next proc
          if (to_ddrank + 1 < dd % domain_n_procs(to_meshbin) .and. &
              index_scatt + 1 > dd % proc_finish(to_ddrank + 1)) then
            to_ddrank = to_ddrank + 1
            pr_bin = (to_bin - 1) * dd % max_domain_procs + to_ddrank + 1
          end if
            
        end if
        
      end do
      
    end do
    
    !===========================================================================
    ! COMMUNICATE SEND INFORMATION
    
    ! Send particle communication info
    do to_bin = 1, N_CARTESIAN_NEIGHBORS
      ! Loop over each of the neighbors we might need to send to
      to_meshbin = dd % neighbor_meshbins(to_bin)
      if (to_meshbin == NO_BIN_FOUND) cycle

      do pr = 1, dd % domain_n_procs(to_meshbin)
        ! Loop over each processor in the receiving bin
        
        pr_bin = (to_bin - 1) * dd % max_domain_procs + pr
        n_request = n_request + 1
        to_rank = dd % domain_masters(to_meshbin) + pr - 1
        call MPI_ISEND(dd % send_rank_info(pr_bin), 1, MPI_INTEGER, to_rank, &
          rank, MPI_COMM_WORLD, request(n_request), mpi_err)

      end do
    end do
    
    ! Receive particle communication info
    dd % recv_rank_info = 0
    do from_bin = 1, N_CARTESIAN_NEIGHBORS
      ! Loop over each of the neighbors we might need to receive from
      from_meshbin = dd % neighbor_meshbins(from_bin)
      if (from_meshbin == NO_BIN_FOUND) cycle

      do pr = 1, dd % domain_n_procs(from_meshbin)
        ! Loop over each processor in the sending bin
        
        pr_bin = (from_bin - 1) * dd % max_domain_procs + pr
        n_request = n_request + 1
        from_rank = dd % domain_masters(from_meshbin) + pr - 1
        call MPI_IRECV(dd % recv_rank_info(pr_bin), 1, MPI_INTEGER, from_rank, &
          from_rank, MPI_COMM_WORLD, request(n_request), mpi_err)
        
      end do
    end do
    
    ! Wait for the scatter information to be synchronized
    call MPI_WAITALL(n_request, request, MPI_STATUSES_IGNORE, mpi_err)
    n_request = 0


!call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)
!    debug = 0
!    do to_bin = 1, 6
!        to_meshbin = dd % neighbor_meshbins(to_bin)
!        if (to_meshbin == NO_BIN_FOUND) cycle
!        debug(to_meshbin) = dd % n_scatters_local(to_bin)
!    end do
!    print *, dd % meshbin, debug
!call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)

    call time_dd_info % stop()

    !===========================================================================
    ! VERIFY BUFFERS

    ! There are certainly more efficient ways to use memory than to have two
    ! separate buffers like what's done here, but this helps for clarity

    ! First check that we have enough space in the inscatter bank
    n_recv = sum(dd % recv_rank_info)
    if (n_recv > dd % size_recv_buffer) then
    
      ! Resize the receive buffer
      n_recv = ceiling(dble(n_recv) * DD_BUFFER_HEADROOM)
      call resize_array(dd % recv_buffer, dd % size_recv_buffer, n_recv, &
          alloc_err)
      
    end if

    ! Also ensure we have enough space in the send buffer
    n_send = sum(dd % n_scatters_local)
    if (n_send > dd % size_send_buffer) then
    
      ! Resize the send buffer
      n_send = ceiling(dble(n_send) * DD_BUFFER_HEADROOM)
      call resize_array(dd % send_buffer, dd % size_send_buffer, n_send, &
          alloc_err)
      
    end if

    !===========================================================================
    ! RECEIVE PARTICLES
 
    ! We post the receives first so that we can match sends as soon as possible.
    ! Supposedly there is an upper limit on the number of sends that can be
    ! outstanding at once.
   
    call time_dd_sendrecv % start()

    local_start = 1
    do from_bin = 1, N_CARTESIAN_NEIGHBORS
      ! Loop through each neighbor we might need to receive from
      from_meshbin = dd % neighbor_meshbins(from_bin)
      if (from_meshbin == NO_BIN_FOUND) cycle
      
      do pr = 1, dd % domain_n_procs(from_meshbin)
        ! Loop over each processor in the sending bin
        
        pr_bin = (from_bin - 1) * dd % max_domain_procs + pr
        n_recv = dd % recv_rank_info(pr_bin)
        if (n_recv == 0) cycle

        n_request = n_request + 1
        from_rank = dd % domain_masters(from_meshbin) + pr - 1
        call MPI_IRECV(dd % recv_buffer(local_start), n_recv, &
            MPI_PARTICLEBUFFER, from_rank, from_rank, MPI_COMM_WORLD, &
            request(n_request), mpi_err)
        local_start = local_start + n_recv

        print *, current_stage,'recv',from_rank,'-->',rank, n_recv

      end do
    end do

    dd % n_inscatt = local_start - 1

    !===========================================================================
    ! SEND PARTICLES

    index_pbuffer = 1

    do to_bin = 1, N_CARTESIAN_NEIGHBORS
      ! Loop over each of the neighbors we might need to send to
      to_meshbin = dd % neighbor_meshbins(to_bin)
      if (to_meshbin == NO_BIN_FOUND) cycle

      start = 1_8
      do pr = 1, dd % domain_n_procs(to_meshbin)
        ! Loop over each processor in the receiving bin
        
        pr_bin = (to_bin - 1) * dd % max_domain_procs + pr
        if (dd % send_rank_info(pr_bin) == 0) cycle 

        pbuffer_start = index_pbuffer
        n_send = 0
        do i8 = start, work
        
          ! Add outscatter particles to send buffer
          if (dd % particle_buffer(i8) % outscatter_destination == to_bin) then
          
            if (dd % particle_buffer(i8) % outscatter_destination == NO_OUTSCATTER) then
              message = "Trying to send particle that didn't leave domain!"
              call fatal_error()
            end if
            
            call particle_to_buffer(dd % particle_buffer(i8), &
                                    dd % send_buffer(index_pbuffer))
            index_pbuffer = index_pbuffer + 1
            n_send = n_send + 1
          end if
          
          if (n_send == dd % send_rank_info(pr_bin)) then
            start = i8
            exit
          end if
          
        end do

        n_request = n_request + 1
        to_rank = dd % domain_masters(to_meshbin) + pr - 1
        call MPI_ISEND(dd % send_buffer(pbuffer_start), n_send, &
            MPI_PARTICLEBUFFER, to_rank, rank, MPI_COMM_WORLD, &
            request(n_request), mpi_err)

        print *, current_stage,'send',rank,'-->',to_rank, n_send

      end do
      
    end do

    ! Wait for the particles to be sent
    call MPI_WAITALL(n_request, request, MPI_STATUSES_IGNORE, mpi_err)

    call time_dd_sendrecv % stop()

    !===========================================================================
    ! REBUILD PARTICLE BUFFER
    
    ! Grow the particle buffer if we don't have enough space
    if (dd % n_inscatt > dd % size_particle_buffer) then
    
      size_buff = ceiling(dble(dd % n_inscatt) * DD_BUFFER_HEADROOM)
      call resize_array(dd % particle_buffer, dd % size_particle_buffer, &
          size_buff, alloc_err)
      
      ! We also need to resize the source bank and fission bank now that we
      ! have more particles
      call resize_array(source_bank, size_source_bank, int(size_buff, 8), &
          alloc_err)
      call extend_array(fission_bank, size_fission_bank, int(3*size_buff, 8), &
          alloc_err)
      !TODO: fission bank might not be big enough! This only covers the particle
      ! buffer from this stage!
      
    end if

    ! Populate the particle source bank from the receive buffer
    do i = 1, dd % n_inscatt
      call buffer_to_particle(dd % recv_buffer(i), dd % particle_buffer(i))
    end do
    
    work = dd % n_inscatt
    
    print *, rank, current_stage, 'n_scatt', dd % n_inscatt

  end subroutine synchronize_particles


end module dd_comm
