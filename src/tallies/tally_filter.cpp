#include "openmc/tallies/tally_filter.h"


namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

std::vector<TallyFilterMatch> filter_matches;

//==============================================================================
// Fortran compatibility functions
//==============================================================================

extern "C" {
  TallyFilterMatch* filter_match_pointer(int indx)
  {return &filter_matches[indx];}
}

} // namespace openmc
