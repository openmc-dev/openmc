module message_passing

#ifdef MPI
#ifdef MPIF08
  use mpi_f08
#else
  use mpi
#endif
#endif

  ! The defaults set here for the number of processors, rank, and master and
  ! mpi_enabled flag are for when MPI is not being used at all, i.e. a serial
  ! run. In this case, these variables are still used at times.

  integer :: n_procs     = 1       ! number of processes
  integer :: rank        = 0       ! rank of process
  logical :: master      = .true.  ! master process?
  logical :: mpi_enabled = .false. ! is MPI in use and initialized?
#ifdef MPIF08
  type(MPI_Datatype) :: MPI_BANK   ! MPI datatype for fission bank
  type(MPI_Comm) :: mpi_intracomm  ! MPI intra-communicator
#else
  integer :: MPI_BANK              ! MPI datatype for fission bank
  integer :: mpi_intracomm         ! MPI intra-communicator
#endif

end module message_passing
