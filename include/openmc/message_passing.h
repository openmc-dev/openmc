#ifndef OPENMC_MESSAGE_PASSING_H
#define OPENMC_MESSAGE_PASSING_H

#include <cstdint>

#ifdef OPENMC_MPI
#include <mpi.h>
#endif

#include "openmc/vector.h"

namespace openmc {
namespace mpi {

extern int rank;
extern int n_procs;
extern bool master;

#ifdef OPENMC_MPI
extern MPI_Datatype source_site;
extern MPI_Comm intracomm;
#endif

// Calculates global indices of the bank particles
// across all ranks using a parallel scan. This is used to write
// the surface source file in parallel runs. It will probably
// be used in the future for other types of bank like particles
// in flight used to kick off transient simulations.
//
// More abstractly, this just takes a number from each MPI rank,
// and returns a vector which is the exclusive parallel scan across
// all of those numbers, having a length of the number of MPI ranks
// plus one.
vector<int64_t> calculate_parallel_index_vector(int64_t size);

} // namespace mpi
} // namespace openmc

#endif // OPENMC_MESSAGE_PASSING_H
