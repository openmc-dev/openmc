#ifndef INITIALIZE_H
#define INITIALIZE_H

#ifdef OPENMC_MPI
#include "mpi.h"
#endif

extern "C" void print_usage();
extern "C" void print_version();

namespace openmc {

int parse_command_line(int argc, char* argv[]);
#ifdef OPENMC_MPI
void initialize_mpi(MPI_Comm intracomm);
#endif

}

#endif // INITIALIZE_H
