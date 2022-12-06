#ifndef OPENMC_INITIALIZE_H
#define OPENMC_INITIALIZE_H

#include <string>

#ifdef OPENMC_MPI
#include "mpi.h"
#endif

namespace openmc {

int parse_command_line(int argc, char* argv[]);
#ifdef OPENMC_MPI
void initialize_mpi(MPI_Comm intracomm);
#endif

//! Read material, geometry, settings, and tallies from a single XML file
bool read_model_xml();
//! Read inputs from separate XML files
void read_separate_xml_files();
//! Write some output that occurs right after initialization
void initial_output();

} // namespace openmc

#endif // OPENMC_INITIALIZE_H
