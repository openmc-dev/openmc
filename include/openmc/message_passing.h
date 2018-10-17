#ifndef OPENMC_MESSAGE_PASSING_H
#define OPENMC_MESSAGE_PASSING_H

#ifdef OPENMC_MPI
#include <mpi.h>
#endif

namespace openmc {
namespace mpi {

  extern int rank;
  extern int n_procs;
  extern bool master;

#ifdef OPENMC_MPI
  extern MPI_Datatype bank;
  extern MPI_Comm intracomm;
#endif

} // namespace mpi
} // namespace openmc

#endif // OPENMC_MESSAGE_PASSING_H
