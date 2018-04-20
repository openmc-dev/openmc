#include "openmc.h"

#ifdef OPENMC_MPI
MPI_Comm openmc_intracomm;
#endif


int
openmc_init(const void* intracomm)
{
#ifdef OPENMC_MPI
  openmc_intracomm = *static_cast<const MPI_Comm *>(intracomm);

  MPI_Fint fcomm = MPI_Comm_c2f(openmc_intracomm);
  openmc_init_f(&fcomm);
#else
  openmc_init_f(nullptr);
#endif
  return 0;
}
