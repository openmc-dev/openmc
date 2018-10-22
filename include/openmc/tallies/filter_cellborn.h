#ifndef OPENMC_TALLIES_FILTER_CELLBORN_H
#define OPENMC_TALLIES_FILTER_CELLBORN_H

#include "openmc/tallies/filter_cell.h"

namespace openmc {

//==============================================================================
//! Specifies which cell the particle was born in.
//==============================================================================

class CellbornFilter : public CellFilter
{
public:
  std::string type() const override {return "cellborn";}

  void
  get_all_bins(Particle* p, int estimator, FilterMatch& match)
  const override
  {
    auto search = map_.find(p->cell_born);
    if (search != map_.end()) {
      // TODO: off-by-one
      match.bins_.push_back(search->second + 1);
      match.weights_.push_back(1);
    }
  }

  std::string
  text_label(int bin) const override
  {
    return "Birth Cell " + std::to_string(cells[cells_[bin-1]]->id_);
  }
};

} // namespace openmc
#endif // OPENMC_TALLIES_FILTER_CELLBORN_H
