module dd_comm
  
  use bank_header,      only: Bank
  use constants
  use dd_header,        only: dd_type
  use error,            only: fatal_error
  use global,           only: domain_decomp, source_bank, MPI_BANK, n_procs, &
                              rank, mpi_err, work, size_source, message, &
                              n_particles, work_index
  use mesh,             only: get_mesh_bin
  use output,           only: write_message
  use string,           only: to_str
  
  use mpi
  
  implicit none
  private
  
  public :: distribute_source
  public :: synchronize_bank_dd

contains

!===============================================================================
! DISTRIBUTE_SOURCE loops through all particles in the source bank and sends
! them to the process that is handling the domain they are located on.  This
! routine also receives particles from other processes, and saves them in the
! source bank, overwriting the previous values.
!===============================================================================

  subroutine distribute_source()
  
    integer(8) :: i           ! loop index over initially-sampled source sites
    integer :: n_source_sites ! number of sites that will start in this domains
    integer :: to_domain
    integer :: to_rank        ! global rank of the processor to send to
    integer :: alloc_err      ! allocation error code
    integer :: n_request      ! number of communication requests
    integer, allocatable    :: request(:) ! communication requests
    integer, allocatable    :: n_send_domain(:) ! number sent to each domain
    integer, allocatable    :: n_send_rank(:) ! number sent to each process
    type(dd_type), pointer  :: dd => domain_decomp
    type(Bank), allocatable :: buffer(:)
    type(Bank), allocatable :: temp(:)

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
    
    ! The initial size of the source_bank is equal to work.  This may change
    size_source = work

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
    if (n_send_rank(rank + 1) + n_source_sites > size_source) then

      ! Allocate new source_bank array of the required size
      allocate(temp(n_send_rank(rank + 1) + n_source_sites))

      ! Copy sites over to temporary array
      temp(1:n_source_sites) = source_bank(1:n_source_sites)

      ! Move allocation from temporary array
      call move_alloc(FROM=temp, TO=source_bank)
    
    end if

    ! Receive sites from other processes
    do i = 1, n_send_rank(rank + 1)
    
      n_source_sites = n_source_sites + 1
      n_request = n_request + 1
      call MPI_IRECV(source_bank(n_source_sites), 1, MPI_BANK, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &
                 request(n_request), mpi_err)
    end do

    ! Wait for all sends/recvs to complete
    call MPI_WAITALL(n_request, request, MPI_STATUSES_IGNORE, mpi_err)

    ! 
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


end module dd_comm
