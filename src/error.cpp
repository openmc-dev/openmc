#include "openmc/error.h"

#include "openmc/message_passing.h"

namespace openmc {

#ifdef OPENMC_MPI
void abort_mpi(int code)
{
    MPI_Abort(mpi::intracomm, code);
}
#endif

} // namespace openmc
