module intercycle

  use, intrinsic :: ISO_FORTRAN_ENV

  use error,           only: fatal_error, warning
  use global
  use mesh,            only: get_mesh_bin
  use mesh_header,     only: StructuredMesh
  use output,          only: write_message
  use random_lcg,      only: prn, set_particle_seed, prn_skip
  use search,          only: binary_search
  use string,          only: to_str
  use tally,           only: accumulate_cycle_estimate
  use tally_header,    only: TallyObject
  use timing,          only: timer_start, timer_stop

#ifdef MPI
  use mpi
#endif

contains

!===============================================================================
! SYNCHRONIZE_BANK samples source sites from the fission sites that were
! accumulated during the cycle. This routine is what allows this Monte Carlo to
! scale to large numbers of processors where other codes cannot.
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
    integer    :: n            ! number of sites to send/recv
    integer    :: neighbor     ! processor to send/recv data from
    integer    :: request(20)  ! communication request for send/recving sites
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
    finish = start + n_bank
    total = finish
    call MPI_BCAST(total, 1, MPI_INTEGER8, n_procs - 1, & 
         MPI_COMM_WORLD, mpi_err)

#else
    start  = 0_8
    finish = n_bank
    total  = n_bank
#endif

    ! If there are not that many particles per cycle, it's possible that no
    ! fission sites were created at all on a single processor. Rather than add
    ! extra logic to treat this circumstance, we really want to ensure the user
    ! runs enough particles to avoid this in the first place.

    if (n_bank == 0) then
       message = "No fission sites banked on processor " // to_str(rank)
       call fatal_error()
    end if

    ! Make sure all processors start at the same point for random sampling by
    ! using the current_cycle as the starting seed. Then skip ahead in the
    ! sequence using the starting index in the 'global' fission bank for each
    ! processor Skip ahead however many random numbers are needed

    call set_particle_seed(int(current_cycle,8))
    call prn_skip(start)

    ! Determine how many fission sites we need to sample from the source bank
    ! and the probability for selecting a site.

    if (total < n_particles) then
       sites_needed = mod(n_particles,total)
    else
       sites_needed = n_particles
    end if
    p_sample = real(sites_needed,8)/real(total,8)

    call timer_start(time_ic_sample)

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
       finish = bank_last
    end if

    call timer_stop(time_ic_sample)
    call timer_start(time_ic_sendrecv)
    
#ifdef MPI
    ! ==========================================================================
    ! SEND BANK SITES TO NEIGHBORS

    index_local = 1
    n_request = 0

    SEND_SITES: do while (start < finish)
       ! Determine the index of the processor which has the first part of the
       ! source_bank for the local processor
       neighbor = start / maxwork

       ! Determine the number of sites to send
       n = min((neighbor + 1)*maxwork, finish) - start

       ! Initiate an asynchronous send of source sites to the neighboring
       ! process
       if (neighbor /= rank) then
          n_request = n_request + 1
          call MPI_ISEND(temp_sites(index_local), n, MPI_BANK, neighbor, &
               rank, MPI_COMM_WORLD, request(n_request), mpi_err)
       end if

       ! Increment all indices
       start       = start       + n
       index_local = index_local + n
       neighbor    = neighbor    + 1
    end do SEND_SITES

    ! ==========================================================================
    ! RECEIVE BANK SITES FROM NEIGHBORS OR TEMPORARY BANK

    start = bank_first - 1
    index_local = 1

    ! Determine what process has the source sites that will need to be stored at
    ! the beginning of this processor's source bank.

    if (start >= bank_position(n_procs)) then
       neighbor = n_procs - 1
    else
       neighbor = binary_search(bank_position, n_procs, start) - 1
    end if

    RECV_SITES: do while (start < bank_last)
       ! Determine how many sites need to be received
       if (neighbor == n_procs - 1) then
          n = min(n_particles, (rank+1)*maxwork) - start
       else
          n = min(bank_position(neighbor+2), min(n_particles, &
               (rank+1)*maxwork)) - start
       end if

       if (neighbor /= rank) then
          ! If the source sites are not on this processor, initiate an
          ! asynchronous receive for the source sites

          n_request = n_request + 1
          call MPI_IRECV(source_bank(index_local), n, MPI_BANK, &
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

    ! Send we initiated a series of asynchronous ISENDs and IRECVs, now we have
    ! to ensure that the data has actually been communicated before moving on to
    ! the next cycle

    call MPI_WAITALL(n_request, request, MPI_STATUSES_IGNORE, mpi_err)

    ! Deallocate space for bank_position on the last cycle
    if (current_cycle == n_cycles) deallocate(bank_position)
#else
    source_bank = temp_sites(1:n_particles)
#endif

    call timer_stop(time_ic_sendrecv)

    ! Deallocate space for the temporary source bank on the last cycle
    if (current_cycle == n_cycles) deallocate(temp_sites)

  end subroutine synchronize_bank

!===============================================================================
! SHANNON_ENTROPY calculates the Shannon entropy of the fission source
! distribution to assess source convergence
!===============================================================================

  subroutine shannon_entropy()

    integer :: i              ! index for bank sites
    integer :: n              ! # of boxes in each dimension
    integer :: bin            ! index in entropy_p
    integer(8) :: total_bank  ! total # of fission bank sites
    integer, save :: n_box    ! total # of boxes on mesh
    logical :: outside_box    ! were there sites outside entropy box?
    type(StructuredMesh), pointer :: m => null()

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
          n = ceiling((n_particles/20)**(1.0/3.0))

          ! copy dimensions
          m % n_dimension = 3
          allocate(m % dimension(3))
          m % dimension = n
       end if

       ! Determine total number of mesh boxes
       n_box = product(m % dimension)

       ! allocate and determine width
       allocate(m % width(3))
       m % width = (m % upper_right - m % lower_left) / m % dimension

       ! allocate p
       allocate(entropy_p(n_box))
    end if

    ! initialize p
    entropy_p = ZERO
    outside_box = .false.

    ! loop over fission sites and count how many are in each mesh box
    FISSION_SITES: do i = 1, int(n_bank,4)
       ! determine scoring bin for entropy mesh
       call get_mesh_bin(m, fission_bank(i) % xyz, bin)

       ! if outside mesh, skip particle
       if (bin == NO_BIN_FOUND) then
          outside_box = .true.
          cycle
       end if

       ! add to appropriate mesh box
       entropy_p(bin) = entropy_p(bin) + 1
    end do FISSION_SITES

    ! display warning message if there were sites outside entropy box
    if (outside_box) then
       message = "Fission source site(s) outside of entropy box."
       call warning()
    end if

#ifdef MPI
    ! collect values from all processors
    if (master) then
       call MPI_REDUCE(MPI_IN_PLACE, entropy_p, n_box, MPI_REAL8, MPI_SUM, &
            0, MPI_COMM_WORLD, mpi_err)
    else
       call MPI_REDUCE(entropy_p, entropy_p, n_box, MPI_REAL8, MPI_SUM, &
            0, MPI_COMM_WORLD, mpi_err)
    end if

    ! determine total number of bank sites
    call MPI_REDUCE(n_bank, total_bank, 1, MPI_INTEGER8, MPI_SUM, 0, &
         MPI_COMM_WORLD, mpi_err)
#else
    total_bank = n_bank
#endif

    ! sum values to obtain shannon entropy
    if (master) then
       ! Normalize to number of bank sites
       entropy_p = entropy_p / total_bank

       entropy = 0
       do i = 1, n_box
          if (entropy_p(i) > 0) then
             entropy = entropy - entropy_p(i) * log(entropy_p(i))/log(2.0)
          end if
       end do
    end if

  end subroutine shannon_entropy

!===============================================================================
! CALCULATE_KEFF calculates the single cycle estimate of keff as well as the
! mean and standard deviation of the mean for active cycles and displays them
!===============================================================================

  subroutine calculate_keff()

    integer :: n        ! active cycle number
    real(8) :: k_cycle  ! single cycle estimate of keff
    real(8) :: global_temp(N_GLOBAL_TALLIES)

    message = "Calculate cycle keff..."
    call write_message(8)

    ! Since the creation of bank sites was originally weighted by the last
    ! cycle keff, we need to multiply by that keff to get the current cycle's
    ! value
#ifdef MPI
    global_tallies(K_ANALOG) % value = n_bank * keff

    ! Copy global tallies into array to be reduced
    global_temp = global_tallies(:) % value

    if (master) then
       call MPI_REDUCE(MPI_IN_PLACE, global_temp, N_GLOBAL_TALLIES, &
            MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, mpi_err)

       ! Transfer values back to global_tallies on master
       global_tallies(:) % value = global_temp
    else
       call MPI_REDUCE(global_temp, global_temp, N_GLOBAL_TALLIES, &
            MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, mpi_err)
       
       ! Reset value on other processors
       global_tallies(:) % value = ZERO
    end if
#else
    global_tallies(K_ANALOG) % value = n_bank * keff
#endif

    ! Collect statistics and print output
    if (master) then
       k_cycle = global_tallies(K_ANALOG) % value/n_particles

       if (current_cycle > n_inactive) then
          ! Active cycle number
          n = current_cycle - n_inactive

          ! Accumulate single cycle realizations of k
          call accumulate_cycle_estimate(global_tallies)

          ! Determine mean and standard deviation of mean
          keff = global_tallies(K_ANALOG) % sum/n
          keff_std = sqrt((global_tallies(K_ANALOG) % sum_sq/n - keff*keff)/n)

          ! Display output for this cycle
          if (current_cycle > n_inactive + 1) then
             if (entropy_on) then
                write(UNIT=OUTPUT_UNIT, FMT=103) current_cycle, k_cycle, &
                     entropy, keff, keff_std
             else
                write(UNIT=OUTPUT_UNIT, FMT=101) current_cycle, k_cycle, &
                     keff, keff_std
             end if
          else
             if (entropy_on) then
                write(UNIT=OUTPUT_UNIT, FMT=102) current_cycle, k_cycle, entropy
             else
                write(UNIT=OUTPUT_UNIT, FMT=100) current_cycle, k_cycle
             end if
          end if
       else
          ! Display output for inactive cycle
          if (entropy_on) then
             write(UNIT=OUTPUT_UNIT, FMT=102) current_cycle, k_cycle, entropy
          else
             write(UNIT=OUTPUT_UNIT, FMT=100) current_cycle, k_cycle
          end if
          keff = k_cycle
       end if
    end if

#ifdef MPI
    ! Broadcast new keff value to all processors
    call MPI_BCAST(keff, 1, MPI_REAL8, 0, MPI_COMM_WORLD, mpi_err)
#endif

100 format (2X,I5,2X,F8.5)
101 format (2X,I5,2X,F8.5,5X,F8.5," +/-",F8.5)
102 format (2X,I5,2X,F8.5,3X,F8.5)
103 format (2X,I5,2X,F8.5,3X,F8.5,3X,F8.5," +/-",F8.5)

  end subroutine calculate_keff

end module intercycle
