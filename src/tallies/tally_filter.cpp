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

  void
  filter_match_bins_push_back(TallyFilterMatch* match, int val)
  {match->bins.push_back(val);}

  void
  filter_match_weights_push_back(TallyFilterMatch* match, double val)
  {match->weights.push_back(val);}

  void
  filter_match_bins_clear(TallyFilterMatch* match)
  {match->bins.clear();}

  void
  filter_match_weights_clear(TallyFilterMatch* match)
  {match->weights.clear();}

  int
  filter_match_bins_size(TallyFilterMatch* match)
  {return match->bins.size();}

  int
  filter_match_bins_data(TallyFilterMatch* match, int indx)
  {return match->bins.at(indx-1);}

  double
  filter_match_weights_data(TallyFilterMatch* match, int indx)
  {return match->weights.at(indx-1);}

  void
  filter_match_bins_set_data(TallyFilterMatch* match, int indx, int val)
  {match->bins.at(indx-1) = val;}
}

} // namespace openmc
