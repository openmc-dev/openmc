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

} // namespace mpi

//==============================================================================
// Fortran compatibility functions
//==============================================================================

#ifdef OPENMC_MPI
extern "C" void
send_int(void* buffer, int count, int dest, int tag)
{
  MPI_Send(buffer, count, MPI_INTEGER, dest, tag, mpi::intracomm);
}

extern "C" void
recv_int(void* buffer, int count, int source, int tag)
{
  MPI_Recv(buffer, count, MPI_INTEGER, source, tag, mpi::intracomm, MPI_STATUS_IGNORE);
}
#endif

} // namespace openmc
