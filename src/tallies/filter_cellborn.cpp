#include "openmc/tallies/filter_cellborn.h"

#include "openmc/cell.h"

namespace openmc {

void
CellbornFilter::get_all_bins(const Particle* p, int estimator,
                             FilterMatch& match) const
{
  auto search = map_.find(p->cell_born);
  if (search != map_.end()) {
    //TODO: off-by-one
    match.bins_.push_back(search->second + 1);
    match.weights_.push_back(1.0);
  }
}

std::string
CellbornFilter::text_label(int bin) const
{
  return "Birth Cell " + std::to_string(cells[cells_[bin-1]]->id_);
}

} // namespace openmc
