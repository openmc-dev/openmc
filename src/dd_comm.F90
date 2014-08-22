module dd_comm
  
  use bank_header,      only: Bank
  use constants
  use dd_header,        only: dd_type
  use error,            only: fatal_error
  use global,           only: domain_decomp, source_bank, MPI_BANK, n_procs, &
                              rank, mpi_err, work, size_source, message, work_index
  use mesh,             only: get_mesh_bin
  use output,           only: write_message
  use string,           only: to_str
  
  use mpi
  
  implicit none
  private
  
  public :: distribute_source

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
    
    ! We'll assume, quite arbitrarily, that each process will never receive more
    ! than twice what was originally allocated from n_particles/n_procs.  In
    ! principle though this could be as large as n_particles if for some reason
    ! one domain should track ALL source sites.  That could only result from a
    ! terrible domain mesh and node distribution.
    allocate(request(3*work), STAT=alloc_err)
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

end module dd_comm
