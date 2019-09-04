#include "openmc/message_passing.h"

namespace openmc {
namespace mpi {

int rank {0};
int n_procs {1};
bool master {true};

#ifdef OPENMC_MPI
MPI_Comm intracomm {MPI_COMM_NULL};
MPI_Datatype bank {MPI_DATATYPE_NULL};
#endif

extern "C" bool openmc_master() { return mpi::master; }

} // namespace mpi

} // namespace openmc
