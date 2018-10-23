#include "openmc/tallies/filter_cellborn.h"

#include "openmc/cell.h"

namespace openmc {

void
CellbornFilter::get_all_bins(Particle* p, int estimator, FilterMatch& match)
const
{
  auto search = map_.find(p->cell_born);
  if (search != map_.end()) {
    match.bins_.push_back(search->second + 1);
    match.weights_.push_back(1);
  }
}

std::string
CellbornFilter::text_label(int bin) const
{
  return "Birth Cell " + std::to_string(cells[cells_[bin-1]]->id_);
}

} // namespace openmc
