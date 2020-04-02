#ifndef OPENMC_TALLIES_FILTERMATCH_H
#define OPENMC_TALLIES_FILTERMATCH_H

#define FILTERMATCH_BINS_WEIGHTS_SIZE 125


namespace openmc {

//==============================================================================
//! Stores bins and weights for filtered tally events.
//==============================================================================

class FilterMatch
{
public:
  //std::vector<int> bins_;
  int bins_[FILTERMATCH_BINS_WEIGHTS_SIZE];
  //std::vector<double> weights_;
  double weights_[FILTERMATCH_BINS_WEIGHTS_SIZE];
  int bins_weights_length_ {0};
  int i_bin_;
  bool bins_present_ {false};
};

class BigFilterMatch
{
public:
  std::vector<int> bins_;
  //int bins_[FILTERMATCH_BINS_WEIGHTS_SIZE];
  std::vector<double> weights_;
  //double weights_[FILTERMATCH_BINS_WEIGHTS_SIZE];
  //int bins_weights_length_ {0};
  int i_bin_;
  bool bins_present_ {false};
};

} // namespace openmc
#endif // OPENMC_TALLIES_FILTERMATCH_H
