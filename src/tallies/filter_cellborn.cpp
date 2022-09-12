#include "openmc/tallies/filter.h"

#include "openmc/cell.h"

namespace openmc {

void
Filter::CellbornFilter_get_all_bins(const Particle& p, TallyEstimator estimator,
                             FilterMatch& match) const
{
  auto search = map_.find(p.cell_born_);
  if (search != map_.end()) {
    //match.bins_.push_back(search->second);
    //match.weights_.push_back(1.0);
    match.push_back(search->second, 1.0);
  }
}

std::string
Filter::CellbornFilter_text_label(int bin) const
{
  return "Birth Cell " + std::to_string(model::cells[cells_[bin]].id_);
}

} // namespace openmc
