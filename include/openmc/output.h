//! \file output.h
//! Functions for ASCII output.

#ifndef OPENMC_OUTPUT_H
#define OPENMC_OUTPUT_H

#include <string>

#include "openmc/particle.h"

namespace openmc {

//! \brief Display the main title banner as well as information about the
//! program developers, version, and date/time which the problem was run.
void title();

//! Display a header block.
//
//! \param msg The main text of the header
//! \param level The lowest verbosity level at which this header is printed
void header(const char* msg, int level);

//! Retrieve a time stamp.
//
//! \return current time stamp (format: "yyyy-mm-dd hh:mm:ss")
std::string time_stamp();

//! Display the attributes of a particle.
extern "C" void print_particle(Particle* p);

//! Display plot information.
void print_plot();

//! Display information regarding cell overlap checking.
void print_overlap_check();

//! Display information about command line usage of OpenMC
void print_usage();

//! Display current version and copright/license information
void print_version();

//! Display header listing what physical values will displayed
void print_columns();

//! Display information about a generation of neutrons
void print_generation();

//! \brief Display last batch's tallied value of the neutron multiplication
//! factor as well as the average value if we're in active batches
void print_batch_keff();

//! Display time elapsed for various stages of a run
void print_runtime();

//! Display results for global tallies including k-effective estimators
void print_results();

void write_tallies();

} // namespace openmc
#endif // OPENMC_OUTPUT_H
