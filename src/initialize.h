#ifndef INITIALIZE_H
#define INITIALIZE_H

#ifdef OPENMC_MPI
#include "mpi.h"
#endif

namespace openmc {

#ifdef OPENMC_MPI
  void initialize_mpi(MPI_Comm intracomm);
#endif

}

#endif // INITIALIZE_H
