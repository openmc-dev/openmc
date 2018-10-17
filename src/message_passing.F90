module message_passing

  use, intrinsic :: ISO_C_BINDING

  ! The defaults set here for the number of processors, rank, and master and are
  ! for when MPI is not being used at all, i.e. a serial run. In this case, these
  ! variables are still used at times.

  integer(C_INT), bind(C, name='openmc_n_procs') :: n_procs = 1  ! number of processes
  integer(C_INT), bind(C, name='openmc_rank')    :: rank    = 0  ! rank of process
  logical(C_BOOL), bind(C, name='openmc_master') :: master  = .true. ! master process?

end module message_passing
