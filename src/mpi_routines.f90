module mpi_routines

  use constants,       only: MAX_LINE_LEN
  use error,           only: fatal_error
  use global
  use output,          only: write_message
  use particle_header, only: Particle, initialize_particle
  use random_lcg,      only: prn, set_particle_seed, prn_skip
  use tally_header,    only: TallyObject

#ifdef MPI
  use mpi
#endif

  implicit none

  integer    :: MPI_BANK   ! MPI datatype for fission bank
  integer(8) :: bank_index ! Fission bank site unique identifier
  real(8)    :: t_sync(4)  ! synchronization time

contains

!===============================================================================
! SETUP_MPI initilizes the Message Passing Interface (MPI) and determines the
! number of processors the problem is being run with as well as the rank of each
! processor.
!===============================================================================

  subroutine setup_mpi()

#ifdef MPI
    integer        :: i
    integer        :: ierr           ! Error status
    integer        :: bank_blocks(4) ! Count for each datatype
    integer        :: bank_types(4)  ! Datatypes
    integer(MPI_ADDRESS_KIND) :: bank_disp(4)   ! Displacements
    integer(MPI_ADDRESS_KIND) :: base
    type(Bank)     :: b

    mpi_enabled = .true.

    ! Initialize MPI
    call MPI_INIT(ierr)
    if (ierr /= MPI_SUCCESS) then
       message = "Failed to initialize MPI."
       call fatal_error()
    end if

    ! Determine number of processors
    call MPI_COMM_SIZE(MPI_COMM_WORLD, n_procs, ierr)
    if (ierr /= MPI_SUCCESS) then
       message = "Could not determine number of processors."
       call fatal_error()
    end if

    ! Determine rank of each processor
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
    if (ierr /= MPI_SUCCESS) then
       message = "Could not determine MPI rank."
       call fatal_error()
    end if

    ! Determine master
    if (rank == 0) then
       master = .true.
    else
       master = .false.
    end if

    ! Determine displacements for MPI_BANK type
    call MPI_GET_ADDRESS(b % id,  bank_disp(1), ierr)
    call MPI_GET_ADDRESS(b % xyz, bank_disp(2), ierr)
    call MPI_GET_ADDRESS(b % uvw, bank_disp(3), ierr)
    call MPI_GET_ADDRESS(b % E,   bank_disp(4), ierr)

    ! Adjust displacements 
    base = bank_disp(1)
    do i = 1, 4
       bank_disp(i) = bank_disp(i) - base
    end do
    
    ! Define MPI_BANK for fission sites
    bank_blocks = (/ 1, 3, 3, 1 /)
    bank_types = (/ MPI_INTEGER8, MPI_REAL8, MPI_REAL8, MPI_REAL8 /)
    call MPI_TYPE_CREATE_STRUCT(4, bank_blocks, bank_disp, & 
         & bank_types, MPI_BANK, ierr)
    call MPI_TYPE_COMMIT(MPI_BANK, ierr)

    t_sync    = ZERO
#else
    ! if no MPI, set processor to master
    mpi_enabled = .false.
    rank = 0
    n_procs = 1
    master = .true.
#endif

  end subroutine setup_mpi

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
    integer(8) :: count           ! index for source bank
    integer(8) :: index           ! index for id -- accounts for all nodes
    integer    :: send_to_left    ! # of bank sites to send/recv to or from left
    integer    :: send_to_right   ! # of bank sites to send/recv to or from right
    integer(8) :: sites_needed    ! # of sites to be sampled
    real(8)    :: p_sample        ! probability of sampling a site
    type(Bank), allocatable :: &
         & temp_sites(:),      & ! local array of extra sites on each node
         & left_bank(:),       & ! bank sites to send/recv to or from left node
         & right_bank(:)         ! bank sites to send/recv to or fram right node

#ifdef MPI
    integer    :: ierr
    integer    :: status(MPI_STATUS_SIZE) ! message status
    integer    :: request         ! communication request for sending sites
    integer    :: request_left    ! communication request for recv sites from left
    integer    :: request_right   ! communication request for recv sites from right
    real(8)    :: t0, t1, t2, t3, t4
#endif

    message = "Collecting number of fission sites..."
    call write_message(8)

#ifdef MPI
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    t0 = MPI_WTIME()

    ! Determine starting index for fission bank and total sites in fission bank
    start = 0_8
    call MPI_EXSCAN(n_bank, start, 1, MPI_INTEGER8, MPI_SUM, & 
         & MPI_COMM_WORLD, ierr)
    finish = start + n_bank
    total = finish
    call MPI_BCAST(total, 1, MPI_INTEGER8, n_procs - 1, & 
         & MPI_COMM_WORLD, ierr)

    t1 = MPI_WTIME()
    t_sync(1) = t_sync(1) + (t1 - t0)
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
    count = 0_8 ! Index for local source_bank
    index = 0_8 ! Index for global source id -- must account for all nodes

    if (total < n_particles) then
       sites_needed = mod(n_particles,total)
    else
       sites_needed = n_particles
    end if
    p_sample = real(sites_needed,8)/real(total,8)

    message = "Sampling fission sites..."
    call write_message(8)

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
             count = count + 1
             temp_sites(count) = fission_bank(i)
          end do
       end if

       ! Randomly sample sites needed
       if (prn() < p_sample) then
          count = count + 1
          temp_sites(count) = fission_bank(i)
       end if
    end do

    ! Now that we've sampled sites, check where the boundaries of data are for
    ! the source bank
#ifdef MPI
    start = 0_8
    call MPI_EXSCAN(count, start, 1, MPI_INTEGER8, MPI_SUM, & 
         & MPI_COMM_WORLD, ierr)
    finish = start + count
    total = finish
    call MPI_BCAST(total, 1, MPI_INTEGER8, n_procs - 1, & 
         & MPI_COMM_WORLD, ierr)
#else
    start  = 0_8
    finish = count
    total  = count
#endif

    ! Determine how many sites to send to adjacent nodes
    send_to_left  = int(bank_first - 1_8 - start, 4)
    send_to_right = int(finish - bank_last, 4)

    if (rank == n_procs - 1) then
       if (total > n_particles) then
          ! If we have extra sites sampled, we will simply discard the extra
          ! ones on the last processor
          if (rank == n_procs - 1) then
             count = count - send_to_right
          end if

       elseif (total < n_particles) then
          ! If we have too few sites, grab sites from the very end of the
          ! fission bank
          sites_needed = n_particles - total
          do i = 1, int(sites_needed,4)
             count = count + 1
             temp_sites(count) = fission_bank(n_bank - sites_needed + i)
          end do
       end if

       ! the last processor should not be sending sites to right
       send_to_right = 0
    end if

#ifdef MPI
    t2 = MPI_WTIME()
    t_sync(2) = t_sync(2) + (t2 - t1)

    message = "Sending fission sites..."
    call write_message(8)

    ! ==========================================================================
    ! SEND BANK SITES TO NEIGHBORS
    allocate(left_bank(abs(send_to_left)))
    allocate(right_bank(abs(send_to_right)))

    if (send_to_right > 0) then
       i = count - send_to_right + 1
       call MPI_ISEND(temp_sites(i), send_to_right, MPI_BANK, rank+1, 0, &
            & MPI_COMM_WORLD, request, ierr)
    else if (send_to_right < 0) then
       call MPI_IRECV(right_bank, -send_to_right, MPI_BANK, rank+1, 1, &
            & MPI_COMM_WORLD, request_right, ierr)
    end if

    if (send_to_left < 0) then
       call MPI_IRECV(left_bank, -send_to_left, MPI_BANK, rank-1, 0, &
            & MPI_COMM_WORLD, request_left, ierr)
    else if (send_to_left > 0) then
       call MPI_ISEND(temp_sites(1), send_to_left, MPI_BANK, rank-1, 1, &
            & MPI_COMM_WORLD, request, ierr)
    end if

    t3 = MPI_WTIME()
    t_sync(3) = t_sync(3) + (t3 - t2)
#endif

    message = "Constructing source bank..."
    call write_message(8)

    ! ==========================================================================
    ! RECONSTRUCT SOURCE BANK
    if (send_to_left < 0 .and. send_to_right >= 0) then
       i = -send_to_left                ! size of first block
       j = int(count,4) - send_to_right ! size of second block
       call copy_from_bank(temp_sites, i+1, j)
#ifdef MPI
       call MPI_WAIT(request_left, status, ierr)
#endif
       call copy_from_bank(left_bank, 1, i)
    else if (send_to_left >= 0 .and. send_to_right < 0) then
       i = int(count,4) - send_to_left ! size of first block
       j = -send_to_right              ! size of second block
       call copy_from_bank(temp_sites(1+send_to_left), 1, i)
#ifdef MPI
       call MPI_WAIT(request_right, status, ierr)
#endif
       call copy_from_bank(right_bank, i+1, j)
    else if (send_to_left >= 0 .and. send_to_right >= 0) then
       i = int(count,4) - send_to_left - send_to_right
       call copy_from_bank(temp_sites(1+send_to_left), 1, i)
    else if (send_to_left < 0 .and. send_to_right < 0) then
       i = -send_to_left
       j = int(count,4)
       k = -send_to_right
       call copy_from_bank(temp_sites, i+1, j)
#ifdef MPI
       call MPI_WAIT(request_left, status, ierr)
#endif
       call copy_from_bank(left_bank, 1, i)
#ifdef MPI
       call MPI_WAIT(request_right, status, ierr)
#endif
       call copy_from_bank(right_bank, i+j+1, k)
    end if

    ! Reset source index
    source_index = 0_8
    
#ifdef MPI
    t4 = MPI_WTIME()
    t_sync(4) = t_sync(4) + (t4 - t3)

    deallocate(left_bank )
    deallocate(right_bank)
#endif
    deallocate(temp_sites)

  end subroutine synchronize_bank

!===============================================================================
! COPY_FROM_BANK
!===============================================================================

  subroutine copy_from_bank(temp_bank, index, n_sites)

    integer,    intent(in) :: n_sites  ! # of bank sites to copy
    type(Bank), intent(in) :: temp_bank(n_sites)
    integer,    intent(in) :: index    ! starting index in source_bank

    integer :: i            ! index in temp_bank
    integer :: index_source ! index in source_bank
    type(Particle), pointer :: p
    
    do i = 1, n_sites
       index_source = index + i - 1
       p => source_bank(index_source)

       p % xyz       = temp_bank(i) % xyz
       p % xyz_local = temp_bank(i) % xyz
       p % last_xyz  = temp_bank(i) % xyz
       p % uvw       = temp_bank(i) % uvw
       p % E         = temp_bank(i) % E
       p % last_E    = p % E

       ! set defaults
       call initialize_particle(p)

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
    integer :: count
    integer :: ierr
    real(8), allocatable :: tally_temp(:,:)
    type(TallyObject), pointer :: t

    do i = 1, n_tallies
       t => tallies(i)

       n = t % n_total_bins
       m = t % n_macro_bins
       count = n*m

       allocate(tally_temp(n,m))

       tally_temp = t % scores(:,:) % val_history

       if (master) then
          ! Description of MPI_IN_PLANE
          call MPI_REDUCE(MPI_IN_PLACE, tally_temp, count, MPI_REAL8, MPI_SUM, &
               0, MPI_COMM_WORLD, ierr)

          ! Transfer values to val_history on master
          t % scores(:,:) % val_history = tally_temp
       else
          ! Receive buffer not significant at other processors
          call MPI_REDUCE(tally_temp, tally_temp, count, MPI_REAL8, MPI_SUM, &
               0, MPI_COMM_WORLD, ierr)

          ! Reset val_history on other processors
          t % scores(:,:) % val_history = 0
       end if

       deallocate(tally_temp)

    end do

  end subroutine reduce_tallies
#endif

end module mpi_routines
