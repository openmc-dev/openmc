#ifndef OPENMC_MCPL_INTERFACE_H
#define OPENMC_MCPL_INTERFACE_H

#include "openmc/particle_data.h"
#include "openmc/vector.h"

#include <string>

namespace openmc {

//==============================================================================
// Constants
//==============================================================================

extern "C" const bool MCPL_ENABLED;

//==============================================================================
// Functions
//==============================================================================

//! Get a vector of source sites from an MCPL file
//
//! \param[in] path  Path to MCPL file
//! \return  Vector of source sites
vector<SourceSite> mcpl_source_sites(std::string path);

//! Write an MCPL source file
//
//! \param[in] filename  Path to MCPL file
//! \param[in] surf_source_bank  Whether to use the surface source bank
void write_mcpl_source_point(
  const char* filename, bool surf_source_bank = false);

} // namespace openmc

#endif // OPENMC_MCPL_INTERFACE_H
