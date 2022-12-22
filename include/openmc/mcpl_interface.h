#ifndef OPENMC_MCPL_INTERFACE_H
#define OPENMC_MCPL_INTERFACE_H

#include "openmc/particle_data.h"

#include <string>
#include <vector>

#ifdef OPENMC_MCPL
#include <mcpl.h>
#endif

namespace openmc {

extern "C" const bool MCPL_ENABLED;

//! Get a vector of source sites from an MCPL file
//
//! \param[in] path  Path to MCPL file
//! \return  Vector of source sites
vector<SourceSite> mcpl_source_sites(std::string path);

} // namespace openmc

#endif // OPENMC_MCPL_INTERFACE_H
