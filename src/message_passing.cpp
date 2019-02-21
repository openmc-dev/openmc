#include "openmc/message_passing.h"

namespace openmc {
namespace mpi {

int rank {0};
int n_procs {1};
bool master {true};

#ifdef OPENMC_MPI
MPI_Comm intracomm;
MPI_Datatype bank;
#endif

extern "C" bool openmc_master() { return mpi::master; }

} // namespace mpi

} // namespace openmc
