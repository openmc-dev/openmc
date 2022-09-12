#ifndef OPENMC_TALLIES_FILTERMATCH_H
#define OPENMC_TALLIES_FILTERMATCH_H

#define FILTERMATCH_BINS_WEIGHTS_SIZE 125

#ifndef DEVICE_PRINTF
#define printf(fmt, ...) (0)
#endif

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

  void push_back(int bin, double weight) {
    if (bins_weights_length_ >= FILTERMATCH_BINS_WEIGHTS_SIZE ) {
      printf("Error: Too many filter matches - tally data will be incorrect. Increase size of FILTERMATCH_BINS_WEIGHTS_SIZE macro.\n");
    } else {
      bins_[bins_weights_length_] = bin;
      weights_[bins_weights_length_] = weight;
      bins_weights_length_++;
    }
  }
  
  void push_back(int bin) {
    if (bins_weights_length_ >= FILTERMATCH_BINS_WEIGHTS_SIZE ) {
      printf("Error: Too many filter matches - tally data will be incorrect. Increase size of FILTERMATCH_BINS_WEIGHTS_SIZE macro.\n");
    } else {
      bins_[bins_weights_length_] = bin;
      bins_weights_length_++;
    }
  }

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
