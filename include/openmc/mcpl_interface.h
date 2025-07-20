#ifndef OPENMC_MCPL_INTERFACE_H
#define OPENMC_MCPL_INTERFACE_H

#include "openmc/particle_data.h"
#include "openmc/span.h"
#include "openmc/vector.h"

#include <string>

namespace openmc {

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
//! \param[in] filename     Path to MCPL file
//! \param[in] source_bank  Vector of SourceSites to write to file for this
//!                         MPI rank.
//! \param[in] collision_track_bank  Vector of CollisionTrackSites to write to file 
//!                         for this MPI rank.
//! \param[in] bank_index   Pointer to vector of site index ranges over all
//!                         MPI ranks.
void write_mcpl_source_point(const char* filename, span<SourceSite> source_bank,
  const vector<int64_t>& bank_index);
} // namespace openmc

#endif // OPENMC_MCPL_INTERFACE_H
