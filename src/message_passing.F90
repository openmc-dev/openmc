module message_passing

#ifdef MPI
#ifdef MPIF08
  use mpi_f08
#else
  use mpi
#endif
#endif

end module message_passing
