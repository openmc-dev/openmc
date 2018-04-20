#include "message_passing.h"

namespace openmc {
namespace mpi {

int rank;
int n_procs;

#ifdef OPENMC_MPI
MPI_Comm intracomm;
MPI_Datatype bank;
#endif

} // namespace mpi
} // namespace openmc
