#include "openmc/finalize.h"

#include "openmc/message_passing.h"

namespace openmc {

void openmc_free_bank()
{
#ifdef OPENMC_MPI
  MPI_Type_free(&mpi::bank);
#endif
}

} // namespace openmc
