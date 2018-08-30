#ifndef OPENMC_EIGENVALUE_H
#define OPENMC_EIGENVALUE_H

#include <cstdint> // for int64_t
#include <vector>

#include "xtensor/xtensor.hpp"

#include "openmc/particle.h"

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

extern std::vector<double> entropy; //!< Shannon entropy at each generation
extern xt::xtensor<double, 1> source_frac; //!< Source fraction for UFS

extern "C" int64_t n_bank;
#pragma omp threadprivate(n_bank)

//==============================================================================
// Non-member functions
//==============================================================================

//! Calculates the Shannon entropy of the fission source distribution to assess
//! source convergence
extern "C" void shannon_entropy();

//! Determines the source fraction in each UFS mesh cell and reweights the
//! source bank so that the sum of the weights is equal to n_particles. The
//! 'source_frac' variable is used later to bias the production of fission sites
extern "C" void ufs_count_sites();

//! Get UFS weight corresponding to particle's location
extern "C" double ufs_get_weight(const Particle* p);

} // namespace openmc

#endif // OPENMC_EIGENVALUE_H
