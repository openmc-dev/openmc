#include "openmc/tallies/filter.h"

#include "openmc/cell.h"

namespace openmc {

void
Filter::CellFromFilter_get_all_bins(const Particle& p, TallyEstimator estimator,
                             FilterMatch& match) const
{
  for (int i = 0; i < p.n_coord_last_; i++) {
    auto search = map_.find(p.cell_last_[i]);
    if (search != map_.end()) {
    //match.bins_.push_back(search->second);
    //match.weights_.push_back(1.0);
    match.push_back(search->second, 1.0);
    }
  }
}

std::string
Filter::CellFromFilter_text_label(int bin) const
{
  return "Cell from " + std::to_string(model::cells[cells_[bin]].id_);
}

} // namespace openmc
