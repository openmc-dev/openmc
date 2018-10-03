#ifndef OPENMC_TALLY_FILTER_CELLFROM_H
#define OPENMC_TALLY_FILTER_CELLFROM_H

#include "openmc/tallies/tally_filter_cell.h"


namespace openmc {

class CellFromFilter : public CellFilter
{
public:
  virtual void
  get_all_bins(Particle* p, int estimator, TallyFilterMatch& match)
  const override
  {
    for (int i = 0; i < p->last_n_coord; i++) {
      auto search = map_.find(p->last_cell[i]);
      if (search != map_.end()) {
        // TODO: off-by-one
        match.bins.push_back(search->second + 1);
        match.weights.push_back(1);
      }
    }
  }
};

} // namespace openmc
#endif // OPENMC_TALLY_FILTER_CELLFROM_H
