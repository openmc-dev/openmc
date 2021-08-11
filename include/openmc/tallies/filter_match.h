#ifndef OPENMC_TALLIES_FILTERMATCH_H
#define OPENMC_TALLIES_FILTERMATCH_H

namespace openmc {

//==============================================================================
//! Stores bins and weights for filtered tally events.
//==============================================================================

class FilterMatch {
public:
  vector<int> bins_;
  vector<double> weights_;
  int i_bin_;
  bool bins_present_ {false};
};

} // namespace openmc
#endif // OPENMC_TALLIES_FILTERMATCH_H
