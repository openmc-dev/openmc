module intercycle

  use ISO_FORTRAN_ENV

  use error,           only: fatal_error, warning
  use global
  use output,          only: write_message
  use particle_header, only: Particle, initialize_particle
  use random_lcg,      only: prn, set_particle_seed, prn_skip
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

  subroutine synchronize_bank(i_cycle)

    integer, intent(in) :: i_cycle

    integer    :: i, j, k         ! loop indices
    integer(8) :: start           ! starting index in local fission bank
    integer(8) :: finish          ! ending index in local fission bank
    integer(8) :: total           ! total sites in global fission bank
    integer(8) :: index_local     ! index for source bank
    integer    :: send_to_left    ! # of bank sites to send/recv to or from left
    integer    :: send_to_right   ! # of bank sites to send/recv to or from right
    integer(8) :: sites_needed    ! # of sites to be sampled
    real(8)    :: p_sample        ! probability of sampling a site
    type(Bank), allocatable :: &
         temp_sites(:), & ! local array of extra sites on each node
         left_bank(:),  & ! bank sites to send/recv to or from left node
         right_bank(:)    ! bank sites to send/recv to or fram right node

#ifdef MPI
    integer    :: status(MPI_STATUS_SIZE) ! message status
    integer    :: request         ! communication request for sending sites
    integer    :: request_left    ! communication request for recv sites from left
    integer    :: request_right   ! communication request for recv sites from right
#endif

    message = "Collecting number of fission sites..."
    call write_message(8)

#ifdef MPI
    ! Determine starting index for fission bank and total sites in fission bank
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

    ! Check if there are no fission sites
    if (total == 0) then
       message = "No fission sites banked!"
       call fatal_error()
    end if

    ! Make sure all processors start at the same point for random sampling
    call set_particle_seed(int(i_cycle,8))

    ! Skip ahead however many random numbers are needed
    call prn_skip(start)

    allocate(temp_sites(2*work))
    index_local = 0_8 ! Index for local source_bank

    if (total < n_particles) then
       sites_needed = mod(n_particles,total)
    else
       sites_needed = n_particles
    end if
    p_sample = real(sites_needed,8)/real(total,8)

    message = "Sampling fission sites..."
    call write_message(8)

    call timer_start(time_ic_sample)

    ! ==========================================================================
    ! SAMPLE N_PARTICLES FROM FISSION BANK AND PLACE IN TEMP_SITES
    do i = 1, int(n_bank,4)

       ! If there are less than n_particles particles banked, automatically add
       ! int(n_particles/total) sites to temp_sites. For example, if you need
       ! 1000 and 300 were banked, this would add 3 source sites per banked site
       ! and the remaining 100 would be randomly sampled.
       if (total < n_particles) then
          do j = 1,int(n_particles/total)
             ! If index is within this node's range, add site to source
             index_local = index_local + 1
             temp_sites(index_local) = fission_bank(i)
          end do
       end if

       ! Randomly sample sites needed
       if (prn() < p_sample) then
          index_local = index_local + 1
          temp_sites(index_local) = fission_bank(i)
       end if
    end do

    ! Now that we've sampled sites, check where the boundaries of data are for
    ! the source bank
#ifdef MPI
    start = 0_8
    call MPI_EXSCAN(index_local, start, 1, MPI_INTEGER8, MPI_SUM, & 
         MPI_COMM_WORLD, mpi_err)
    finish = start + index_local
    total = finish
    call MPI_BCAST(total, 1, MPI_INTEGER8, n_procs - 1, & 
         MPI_COMM_WORLD, mpi_err)
#else
    start  = 0_8
    finish = index_local
    total  = index_local
#endif

    ! Determine how many sites to send to adjacent nodes
    send_to_left  = int(bank_first - 1_8 - start, 4)
    send_to_right = int(finish - bank_last, 4)

    ! Check to make sure number of sites is not more than size of bank
    if (abs(send_to_left) > work .or. abs(send_to_right) > work) then
       message = "Tried sending sites to neighboring process greater than " &
            // "the size of the source bank."
       call fatal_error()
    end if

    if (rank == n_procs - 1) then
       if (total > n_particles) then
          ! If we have extra sites sampled, we will simply discard the extra
          ! ones on the last processor
          if (rank == n_procs - 1) then
             index_local = index_local - send_to_right
          end if

       elseif (total < n_particles) then
          ! If we have too few sites, grab sites from the very end of the
          ! fission bank
          sites_needed = n_particles - total
          do i = 1, int(sites_needed,4)
             index_local = index_local + 1
             temp_sites(index_local) = fission_bank(n_bank - sites_needed + i)
          end do
       end if

       ! the last processor should not be sending sites to right
       send_to_right = 0
    end if

    call timer_stop(time_ic_sample)
    call timer_start(time_ic_sendrecv)

#ifdef MPI
    message = "Sending fission sites..."
    call write_message(8)

    ! ==========================================================================
    ! SEND BANK SITES TO NEIGHBORS
    allocate(left_bank(abs(send_to_left)))
    allocate(right_bank(abs(send_to_right)))

    if (send_to_right > 0) then
       i = index_local - send_to_right + 1
       call MPI_ISEND(temp_sites(i), send_to_right, MPI_BANK, rank+1, 0, &
            MPI_COMM_WORLD, request, mpi_err)
    else if (send_to_right < 0) then
       call MPI_IRECV(right_bank, -send_to_right, MPI_BANK, rank+1, 1, &
            MPI_COMM_WORLD, request_right, mpi_err)
    end if

    if (send_to_left < 0) then
       call MPI_IRECV(left_bank, -send_to_left, MPI_BANK, rank-1, 0, &
            MPI_COMM_WORLD, request_left, mpi_err)
    else if (send_to_left > 0) then
       call MPI_ISEND(temp_sites(1), send_to_left, MPI_BANK, rank-1, 1, &
            MPI_COMM_WORLD, request, mpi_err)
    end if
#endif

    call timer_stop(time_ic_sendrecv)
    call timer_start(time_ic_rebuild)

    message = "Constructing source bank..."
    call write_message(8)

    ! ==========================================================================
    ! RECONSTRUCT SOURCE BANK
    if (send_to_left < 0 .and. send_to_right >= 0) then
       i = -send_to_left                      ! size of first block
       j = int(index_local,4) - send_to_right ! size of second block
       call copy_from_bank(temp_sites, i+1, j)
#ifdef MPI
       call MPI_WAIT(request_left, status, mpi_err)
#endif
       call copy_from_bank(left_bank, 1, i)
    else if (send_to_left >= 0 .and. send_to_right < 0) then
       i = int(index_local,4) - send_to_left ! size of first block
       j = -send_to_right                    ! size of second block
       call copy_from_bank(temp_sites(1+send_to_left), 1, i)
#ifdef MPI
       call MPI_WAIT(request_right, status, mpi_err)
#endif
       call copy_from_bank(right_bank, i+1, j)
    else if (send_to_left >= 0 .and. send_to_right >= 0) then
       i = int(index_local,4) - send_to_left - send_to_right
       call copy_from_bank(temp_sites(1+send_to_left), 1, i)
    else if (send_to_left < 0 .and. send_to_right < 0) then
       i = -send_to_left
       j = int(index_local,4)
       k = -send_to_right
       call copy_from_bank(temp_sites, i+1, j)
#ifdef MPI
       call MPI_WAIT(request_left, status, mpi_err)
#endif
       call copy_from_bank(left_bank, 1, i)
#ifdef MPI
       call MPI_WAIT(request_right, status, mpi_err)
#endif
       call copy_from_bank(right_bank, i+j+1, k)
    end if

    ! Reset source index
    source_index = 0_8

    call timer_stop(time_ic_rebuild)
    
#ifdef MPI
    deallocate(left_bank)
    deallocate(right_bank)
#endif
    deallocate(temp_sites)

  end subroutine synchronize_bank

!===============================================================================
! COPY_FROM_BANK
!===============================================================================

  subroutine copy_from_bank(temp_bank, i_start, n_sites)

    integer,    intent(in) :: n_sites  ! # of bank sites to copy
    type(Bank), intent(in) :: temp_bank(n_sites)
    integer,    intent(in) :: i_start  ! starting index in source_bank

    integer :: i        ! index in temp_bank
    integer :: i_source ! index in source_bank
    type(Particle), pointer :: p
    
    do i = 1, n_sites
       i_source = i_start + i - 1
       p => source_bank(i_source)

       ! set defaults
       call initialize_particle(p)

       p % coord % xyz = temp_bank(i) % xyz
       p % coord % uvw = temp_bank(i) % uvw
       p % last_xyz    = temp_bank(i) % xyz
       p % E           = temp_bank(i) % E
       p % last_E      = temp_bank(i) % E

    end do

  end subroutine copy_from_bank

!===============================================================================
! REDUCE_TALLIES collects all the results from tallies onto one processor
!===============================================================================

#ifdef MPI
  subroutine reduce_tallies()

    integer :: i
    integer :: n
    integer :: m
    integer :: n_bins
    real(8), allocatable :: tally_temp(:,:)
    type(TallyObject), pointer :: t

    do i = 1, n_tallies
       t => tallies(i)

       n = t % n_total_bins
       m = t % n_macro_bins
       n_bins = n*m

       allocate(tally_temp(n,m))

       tally_temp = t % scores(:,:) % val_history

       if (master) then
          ! The MPI_IN_PLACE specifier allows the master to copy values into a
          ! receive buffer without having a temporary variable
          call MPI_REDUCE(MPI_IN_PLACE, tally_temp, n_bins, MPI_REAL8, MPI_SUM, &
               0, MPI_COMM_WORLD, mpi_err)

          ! Transfer values to val_history on master
          t % scores(:,:) % val_history = tally_temp
       else
          ! Receive buffer not significant at other processors
          call MPI_REDUCE(tally_temp, tally_temp, n_bins, MPI_REAL8, MPI_SUM, &
               0, MPI_COMM_WORLD, mpi_err)

          ! Reset val_history on other processors
          t % scores(:,:) % val_history = 0
       end if

       deallocate(tally_temp)

    end do

  end subroutine reduce_tallies
#endif

!===============================================================================
! SHANNON_ENTROPY calculates the Shannon entropy of the fission source
! distribution to assess source convergence
!===============================================================================

  subroutine shannon_entropy()

    integer :: i              ! x-index for entropy mesh
    integer :: j              ! y-index for entropy mesh
    integer :: k              ! z-index for entropy mesh
    integer :: m              ! index for bank sites
    integer(8) :: total_bank  ! total # of fission bank sites
    integer, save :: n_box    ! total # of boxes on mesh
    integer, save :: n        ! # of boxes in each dimension
    real(8), save :: width(3) ! width of box in each dimension
    logical :: outside_box    ! were there sites outside entropy box?

    ! On the first pass through this subroutine, we need to determine how big
    ! the entropy mesh should be in each direction and then allocate a
    ! three-dimensional array to store the fraction of source sites in each mesh
    ! box

    if (.not. allocated(entropy_p)) then
       ! determine number of boxes in each direction
       n = ceiling((n_particles/20)**(1.0/3.0))
       n_box = n*n*n

       ! determine width
       width = (entropy_upper_right - entropy_lower_left)/n

       ! allocate p
       allocate(entropy_p(n,n,n))
    end if

    ! initialize p
    entropy_p = ZERO
    outside_box = .false.

    ! loop over fission sites and count how many are in each mesh box
    FISSION_SITES: do m = 1, int(n_bank,4)
       ! determine indices for entropy mesh box
       i = int((fission_bank(m) % xyz(1) - entropy_lower_left(1))/width(1)) + 1
       j = int((fission_bank(m) % xyz(2) - entropy_lower_left(2))/width(2)) + 1 
       k = int((fission_bank(m) % xyz(3) - entropy_lower_left(3))/width(3)) + 1

       ! if outside mesh, skip particle
       if (i < 1 .or. i > n .or. j < 1 .or. &
            j > n .or. k < 1 .or. k > n) then
          outside_box = .true.
          cycle
       end if

       ! add to appropriate mesh box
       entropy_p(i,j,k) = entropy_p(i,j,k) + 1
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
       entropy_p = entropy_p / total_bank
       entropy = -sum(entropy_p * log(entropy_p)/log(2.0), entropy_p > ZERO)
    end if

  end subroutine shannon_entropy

!===============================================================================
! CALCULATE_KEFF calculates the single cycle estimate of keff as well as the
! mean and standard deviation of the mean for active cycles and displays them
!===============================================================================

  subroutine calculate_keff(i_cycle)

    integer, intent(in) :: i_cycle ! index of current cycle

    integer(8)    :: total_bank ! total number of source sites
    integer       :: n          ! active cycle number
    real(8)       :: k_cycle    ! single cycle estimate of keff
    real(8), save :: k_sum      ! accumulated keff
    real(8), save :: k_sum_sq   ! accumulated keff**2

    message = "Calculate cycle keff..."
    call write_message(8)

    ! initialize sum and square of sum at beginning of run
    if (i_cycle == 1) then
       k_sum = ZERO
       k_sum_sq = ZERO
    end if

#ifdef MPI
    ! Collect number bank sites onto master process
    call MPI_REDUCE(n_bank, total_bank, 1, MPI_INTEGER8, MPI_SUM, 0, &
         MPI_COMM_WORLD, mpi_err)
#else
    total_bank = n_bank
#endif

    ! Collect statistics and print output
    if (master) then
       ! Since the creation of bank sites was originally weighted by the last
       ! cycle keff, we need to multiply by that keff to get the current cycle's
       ! value

       k_cycle = real(total_bank)/real(n_particles)*keff

       if (i_cycle > n_inactive) then
          ! Active cycle number
          n = i_cycle - n_inactive

          ! Accumulate cycle estimate of k
          k_sum =    k_sum    + k_cycle
          k_sum_sq = k_sum_sq + k_cycle*k_cycle

          ! Determine mean and standard deviation of mean
          keff = k_sum/n
          keff_std = sqrt((k_sum_sq/n - keff*keff)/n)

          ! Display output for this cycle
          if (i_cycle > n_inactive+1) then
             if (entropy_on) then
                write(UNIT=OUTPUT_UNIT, FMT=103) i_cycle, k_cycle, entropy, &
                     keff, keff_std
             else
                write(UNIT=OUTPUT_UNIT, FMT=101) i_cycle, k_cycle, keff, keff_std
             end if
          else
             if (entropy_on) then
                write(UNIT=OUTPUT_UNIT, FMT=102) i_cycle, k_cycle, entropy
             else
                write(UNIT=OUTPUT_UNIT, FMT=100) i_cycle, k_cycle
             end if
          end if
       else
          ! Display output for inactive cycle
          if (entropy_on) then
             write(UNIT=OUTPUT_UNIT, FMT=102) i_cycle, k_cycle, entropy
          else
             write(UNIT=OUTPUT_UNIT, FMT=100) i_cycle, k_cycle
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
