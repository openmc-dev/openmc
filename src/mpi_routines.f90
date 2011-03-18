module mpi_routines

  use mpi
  use global
  use output, only: message, error

  implicit none

  integer    :: MPI_bank   ! MPI datatype for fission bank
  integer(8) :: bank_index ! Fission bank site unique identifier

contains

!=====================================================================
! SETUP_MPI initilizes the Message Passing Interface (MPI) and
! determines the number of processors the problem is being run with as
! well as the rank of each processor.
!=====================================================================

  subroutine setup_mpi()

    integer :: ierr
    character(250) :: msg

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

  end subroutine setup_mpi
  
end module mpi_routines
