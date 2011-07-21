program unittest

  use global
  use types,       only: Bank
  use energy_grid, only: add_grid_points
  use mpi_routines, only: setup_mpi, MPI_BANK
  use mpi

  call test_mpi_bank()

contains

!===============================================================================
! TEST_ENERGY_GRID is a test of the energy grid unionization subroutine on a
! simplified set of data
!===============================================================================

  subroutine test_energy_grid()

    type(ListReal), pointer :: list => null()
    type(ListReal), pointer :: current => null()
    real(8), allocatable :: energy(:)
    integer :: n, i
    integer :: n_list

    n = 10
    allocate(energy(n))
    do i = 1,n
       energy(i) = 1.1248_8*i
    end do

    ! create list
    allocate(list)
    current => list
    current%data = 4.23_8
    allocate(current%next)
    current => current%next
    current%data = 7.23_8

    call add_grid_points(list, energy, n)

    current => list
    do while(associated(current))
       print *,current%data
       current => current%next
    end do

    deallocate(energy)

  end subroutine test_energy_grid

!===============================================================================
! TEST_MPI_BANK is a test of the derived datatype MPI_BANK that is used to send
! fission bank sites during synchronization
!===============================================================================

  subroutine test_mpi_bank()

    integer :: ierr
    type(Bank) :: temp_bank(2)

    call setup_mpi()

    if (rank == 0) then
       temp_bank(1) % uid = 123
       temp_bank(1) % xyz = (/ 13.0_8, 9.0_8, 8.0_8 /)
       temp_bank(1) % uvw = (/ 7.5_8, 6.0_8, 2.0_8 /)
       temp_bank(1) % E   = 4.592834_8

       temp_bank(2) % uid = 321
       temp_bank(2) % xyz = (/ 23.0_8, 29.0_8, 28.0_8 /)
       temp_bank(2) % uvw = (/ 1.0_8, 9.0_8, 11.0_8 /)
       temp_bank(2) % E   = 0.327438_8
       call MPI_SEND(temp_bank, 2, MPI_BANK, 1, 3, MPI_COMM_WORLD, ierr)
    elseif (rank == 1) then
       call MPI_RECV(temp_bank, 2, MPI_BANK, 0, 3, MPI_COMM_WORLD, &
            & MPI_STATUS_IGNORE, ierr)
       print *, temp_bank(1) % uid
       print *, temp_bank(1) % xyz
       print *, temp_bank(1) % uvw
       print *, temp_bank(1) % E
       print *, temp_bank(2) % uid
       print *, temp_bank(2) % xyz
       print *, temp_bank(2) % uvw
       print *, temp_bank(2) % E
    end if

    call MPI_FINALIZE(ierr)
    
  end subroutine test_mpi_bank

end program unittest
