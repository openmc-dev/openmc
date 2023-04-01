#ifndef OPENMC_CMFD_SOLVER_H
#define OPENMC_CMFD_SOLVER_H

namespace openmc {

//==============================================================================
// Constants
//==============================================================================

// For non-accelerated regions on coarse mesh overlay
constexpr int CMFD_NOACCEL {-1};

//==============================================================================
// Non-member functions
//==============================================================================

void free_memory_cmfd();

} // namespace openmc

#endif // OPENMC_CMFD_SOLVER_H
