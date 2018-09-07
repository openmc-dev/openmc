#ifndef OPENMC_TALLY_FILTER_H
#define OPENMC_TALLY_FILTER_H

#include <cstdint>
#include <vector>


namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

extern "C" int32_t n_filters;

class TallyFilterMatch;
extern std::vector<TallyFilterMatch> filter_matches;
#pragma omp threadprivate(filter_matches)

//==============================================================================
//! Stores bins and weights for filtered tally events
//==============================================================================

class TallyFilterMatch
{
public:
  int i_bin;
  std::vector<int> bins;
  std::vector<double> weights;
  bool bins_present;
};

} // namespace openmc
#endif // OPENMC_TALLY_FILTER
