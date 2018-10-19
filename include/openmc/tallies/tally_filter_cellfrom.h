#ifndef OPENMC_TALLY_FILTER_CELLFROM_H
#define OPENMC_TALLY_FILTER_CELLFROM_H

#include "openmc/tallies/tally_filter_cell.h"


namespace openmc {

//==============================================================================
//! Specifies which geometric cells particles exit when crossing a surface.
//==============================================================================

class CellFromFilter : public CellFilter
{
public:
  virtual std::string type() const override {return "cellfrom";}

  virtual void
  get_all_bins(Particle* p, int estimator, TallyFilterMatch& match)
  const override
  {
    for (int i = 0; i < p->last_n_coord; i++) {
      auto search = map_.find(p->last_cell[i]);
      if (search != map_.end()) {
        // TODO: off-by-one
        match.bins_.push_back(search->second + 1);
        match.weights_.push_back(1);
      }
    }
  }

  virtual std::string
  text_label(int bin) const override
  {
    return "Cell from " + std::to_string(cells[cells_[bin-1]]->id_);
  }
};

} // namespace openmc
#endif // OPENMC_TALLY_FILTER_CELLFROM_H
