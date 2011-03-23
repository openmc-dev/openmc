module mpi_routines

  use mpi
  use global
  use output, only: message, error
  use mcnp_random, only: rang
  use source, only: copy_from_bank, source_index

  implicit none

  integer    :: MPI_BANK   ! MPI datatype for fission bank
  integer(8) :: bank_index ! Fission bank site unique identifier

contains

!=====================================================================
! SETUP_MPI initilizes the Message Passing Interface (MPI) and
! determines the number of processors the problem is being run with as
! well as the rank of each processor.
!=====================================================================

  subroutine setup_mpi()

    integer        :: i
    integer        :: ierr           ! Error status
    integer        :: bank_blocks(4) ! Count for each datatype
    integer        :: bank_types(4)  ! Datatypes
    integer(MPI_ADDRESS_KIND) :: bank_disp(4)   ! Displacements
    integer(MPI_ADDRESS_KIND) :: base
    character(250) :: msg            ! Error message
    type(Bank)     :: b

    ! Initialize MPI
    call MPI_INIT(ierr)
    if (ierr /= MPI_SUCCESS) then
       msg = "Failed to initialize MPI."
       call error(msg)
    end if
    mpi_enabled = .true.

    ! Determine number of processors
    call MPI_COMM_SIZE(MPI_COMM_WORLD, n_procs, ierr)
    if (ierr /= MPI_SUCCESS) then
       msg = "Could not determine number of processors."
       call error(msg)
    end if

    ! Determine rank of each processor
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
    if (ierr /= MPI_SUCCESS) then
       msg = "Could not determine MPI rank."
       call error(msg)
    end if

    ! Determine master
    if (rank == 0) then
       master = .true.
    else
       master = .false.
    end if

    ! Determine displacements for MPI_BANK type
    call MPI_GET_ADDRESS(b % uid, bank_disp(1), ierr)
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
    bank_types = (/ MPI_INTEGER, MPI_REAL8, MPI_REAL8, MPI_REAL8 /)
    call MPI_TYPE_CREATE_STRUCT(4, bank_blocks, bank_disp, & 
         & bank_types, MPI_BANK, ierr)
    call MPI_TYPE_COMMIT(MPI_BANK, ierr)

  end subroutine setup_mpi

!=====================================================================
! SYNCHRONIZE_BANK
!=====================================================================

  subroutine synchronize_bank()

    implicit none

    integer :: i, j, k                 ! loop indices
    integer :: status(MPI_STATUS_SIZE) ! message status
    integer :: ierr
    integer :: start             ! starting index in local fission bank
    integer :: finish            ! ending index in local fission bank
    integer :: total             ! total sites in global fission bank
    integer :: count             ! index for source bank
    integer :: index             ! index for uid -- accounts for all nodes
    integer :: send_to_left      ! # of bank sites to send/recv to or from left
    integer :: send_to_right     ! # of bank sites to send/recv to or from right
    integer :: sites_needed      ! # of sites to be sampled
    integer :: sites_remaining   ! # of sites left in fission bank
    type(Bank), allocatable :: &
         & temp_sites(:),      & ! local array of extra sites on each node
         & left_bank(:),       & ! bank sites to send/recv to or from left node
         & right_bank(:)         ! bank sites to send/recv to or fram right node
    character(250) :: msg

    msg = "Synchronizing fission bank..."
    call message(msg, 8)

    ! Determine starting index for fission bank and total sites in
    ! fission bank
    start = 0
    call MPI_EXSCAN(n_bank, start, 1, MPI_INTEGER, MPI_SUM, & 
         & MPI_COMM_WORLD, ierr)
    finish = start + n_bank
    total = finish
    call MPI_BCAST(total, 1, MPI_INTEGER, n_procs - 1, & 
         & MPI_COMM_WORLD, ierr)

    ! Check if there are no fission sites
    if (total == 0) then
       msg = "No fission sites banked!"
       call error(msg)
    end if

    allocate(temp_sites(2*work))
    count = 0 ! Index for source_bank
    index = 0 ! Index for uid -- must account for all nodes

    if (total < n_particles) then
       sites_needed = mod(n_particles,total)
    else
       sites_needed = n_particles
    end if

    ! ================================================================
    ! SAMPLE N_PARTICLES FROM FISSION BANK AND PLACE IN TEMP_SITES
    do i = 1, total

       ! If there are less than n_particles particles banked,
       ! automatically add int(n_particles/total) sites to
       ! temp_sites. For example, if you need 1000 and 300 were
       ! banked, this would add 3 source sites per banked site and the
       ! remaining 100 would be randomly sampled.
       if (total < n_particles) then
          do j = 1,int(n_particles/total)
             index = index + 1
             ! If index is within this node's range, add site to source
             if (i > start .and. i <= finish) then
                count = count + 1
                temp_sites(count) = fission_bank(i-start)
                temp_sites(count) % uid = index
             end if
          end do
       end if

       ! Randomly sample sites needed
       sites_remaining = total - i + 1
       if (sites_needed == sites_remaining .or. &
            rang() < real(sites_needed)/real(sites_remaining)) then
          index = index + 1
          if (i > start .and. i <= finish) then
             count = count + 1
             temp_sites(count) = fission_bank(i-start)
             temp_sites(count) % uid = index
          end if
          sites_needed = sites_needed - 1
       end if
    end do

    ! Determine how many sites to send to adjacent nodes
    send_to_left  = bank_first - temp_sites(1)%uid
    send_to_right = temp_sites(1)%uid + (count-1) - bank_last

    ! ================================================================
    ! SEND BANK SITES TO NEIGHBORS
    allocate(left_bank(abs(send_to_left)))
    allocate(right_bank(abs(send_to_right)))

    if (send_to_right > 0) then
       i = count - send_to_right + 1
       call MPI_SEND(temp_sites(i), send_to_right, MPI_BANK, rank+1, 0, &
            & MPI_COMM_WORLD, ierr)
    else if (send_to_right < 0) then
       call MPI_RECV(right_bank, -send_to_right, MPI_BANK, rank+1, 1, &
            & MPI_COMM_WORLD, status, ierr)
    end if

    if (send_to_left < 0) then
       call MPI_RECV(left_bank, -send_to_left, MPI_BANK, rank-1, 0, &
            & MPI_COMM_WORLD, status, ierr)
    else if (send_to_left > 0) then
       call MPI_SEND(temp_sites(1), send_to_left, MPI_BANK, rank-1, 1, &
            & MPI_COMM_WORLD, ierr)
    end if

    ! ================================================================
    ! RECONSTRUCT SOURCE BANK
    if (send_to_left < 0 .and. send_to_right >= 0) then
       i = -send_to_left         ! size of first block
       j = count - send_to_right ! size of second block
       call copy_from_bank(left_bank, 1, i)
       call copy_from_bank(temp_sites, i+1, j)
    else if (send_to_left >= 0 .and. send_to_right < 0) then
       i = count - send_to_left ! size of first block
       j = -send_to_right       ! size of second block
       call copy_from_bank(temp_sites(1+send_to_left), 1, i)
       call copy_from_bank(right_bank, i+1, j)
    else if (send_to_left >= 0 .and. send_to_right >= 0) then
       i = count - send_to_left - send_to_right
       call copy_from_bank(temp_sites(1+send_to_left), 1, i)
    else if (send_to_left < 0 .and. send_to_right < 0) then
       i = -send_to_left
       j = count
       k = -send_to_right
       call copy_from_bank(left_bank, 1, i)
       call copy_from_bank(temp_sites, i+1, j)
       call copy_from_bank(right_bank, i+j+1, k)
    end if

    ! Reset source index
    source_index = 0
    
    deallocate(left_bank )
    deallocate(right_bank)
    deallocate(temp_sites)

  end subroutine synchronize_bank
  
end module mpi_routines
