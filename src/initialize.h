#ifndef INITIALIZE_H
#define INITIALIZE_H

#include "mpi.h"

namespace openmc {

#ifdef OPENMC_MPI
  void initialize_mpi(MPI_Comm intracomm);
#endif

}

#endif // INITIALIZE_H
