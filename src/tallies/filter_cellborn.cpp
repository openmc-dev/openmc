#include "openmc/tallies/filter_cellborn.h"

#include "openmc/cell.h"

namespace openmc {

void CellBornFilter::get_all_bins(
  const Particle& p, TallyEstimator estimator, FilterMatch& match) const
{
  auto search = map_.find(p.cell_born());
  if (search != map_.end()) {
    match.set(search->second);
  }
}

std::string CellBornFilter::text_label(int bin) const
{
  return "Birth Cell " + std::to_string(model::cells[cells_[bin]]->id_);
}

} // namespace openmc
