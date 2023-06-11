#ifndef OPENMC_MCPL_INTERFACE_H
#define OPENMC_MCPL_INTERFACE_H

#include "openmc/particle_data.h"
#include "openmc/vector.h"

#include <gsl/gsl-lite.hpp>

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
//! \param[in] filename     Path to MCPL file
//! \param[in] source_bank  Vector of SourceSites to write to file for this
//!                         MPI rank
//! \param[in] bank_indx    Pointer to vector of site index ranges over all
//!                         MPI ranks. This can be computed by calling
//!                         calculate_parallel_index_vector on
//!                         source_bank.size().
void write_mcpl_source_point(const char* filename,
  gsl::span<SourceSite> source_bank, vector<int64_t> const& bank_index);
} // namespace openmc

#endif // OPENMC_MCPL_INTERFACE_H
