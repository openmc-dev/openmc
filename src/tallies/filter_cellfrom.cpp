#include "openmc/tallies/filter_cellfrom.h"

#include "openmc/cell.h"

namespace openmc {

void
CellFromFilter::get_all_bins(const Particle* p, int estimator,
                             FilterMatch& match) const
{
  for (int i = 0; i < p->last_n_coord; i++) {
    auto search = map_.find(p->last_cell[i]);
    if (search != map_.end()) {
      //TODO: off-by-one
      match.bins_.push_back(search->second + 1);
      match.weights_.push_back(1.0);
    }
  }
}

std::string
CellFromFilter::text_label(int bin) const
{
  //TODO: off-by-one
  return "Cell from " + std::to_string(cells[cells_[bin-1]]->id_);
}

} // namespace openmc
