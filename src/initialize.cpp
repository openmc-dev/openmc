#include "openmc.h"

#include "message_passing.h"


int
openmc_init(const void* intracomm)
{
#ifdef OPENMC_MPI
  openmc::mpi::intracomm = *static_cast<const MPI_Comm *>(intracomm);

  MPI_Fint fcomm = MPI_Comm_c2f(openmc::mpi::intracomm);
  openmc_init_f(&fcomm);
#else
  openmc_init_f(nullptr);
#endif
  return 0;
}
