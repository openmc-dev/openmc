//! \file output.h
//! Functions for ASCII output.

#ifndef OPENMC_OUTPUT_H
#define OPENMC_OUTPUT_H

#include <string>

namespace openmc {

//==============================================================================
//! Display a header block.
//!
//! \param msg The main text of the header
//! \param level The lowest verbosity level at which this header is printed
//==============================================================================

void header(const char* msg, int level);


//==============================================================================
//! Retrieve a time stamp with the format "yyyy-mm-dd hh:mm:ss"
//!
//==============================================================================
 
std::string time_stamp();
 
//==============================================================================
//! Display information regarding cell overlap checking.
//==============================================================================

void print_overlap_check();

extern "C" void title();

} // namespace openmc
#endif // OPENMC_OUTPUT_H
