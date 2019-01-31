#include "openmc/tallies/trigger.h"

#include <cmath>

#include "openmc/capi.h"
#include "openmc/tallies/tally.h"

namespace openmc {

extern "C" void
get_tally_uncertainty(int i_tally, int score_index, int filter_index,
  double* std_dev, double* rel_err)
{
  int n;
  int err = openmc_tally_get_n_realizations(i_tally, &n);

  auto results = tally_results(i_tally);
  //TODO: off-by-one
  auto sum = results(filter_index-1, score_index-1, RESULT_SUM);
  auto sum_sq = results(filter_index-1, score_index-1, RESULT_SUM_SQ);

  auto mean = sum / n;
  *std_dev = std::sqrt((sum_sq/n - mean*mean) / (n-1));

  if (mean > 0) {
    *rel_err = *std_dev / mean;
  } else {
    *rel_err = 0.;
  }
}

} // namespace openmc
