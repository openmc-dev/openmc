#include "openmc/message_passing.h"

namespace openmc {
namespace mpi {

int rank {0};
int n_procs {1};

#ifdef OPENMC_MPI
MPI_Comm intracomm;
MPI_Datatype bank;
#endif

} // namespace mpi
} // namespace openmc
