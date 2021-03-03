#ifndef OPENMC_INITIALIZE_H
#define OPENMC_INITIALIZE_H

#ifdef OPENMC_MPI
#include "mpi.h"
#endif

namespace openmc {

int parse_command_line(int argc, char* argv[]);
#ifdef OPENMC_MPI
void initialize_mpi(MPI_Comm intracomm);
#endif
void read_input_xml();

}

#endif // OPENMC_INITIALIZE_H
