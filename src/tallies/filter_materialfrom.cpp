#include "openmc/tallies/filter_materialfrom.h"

#include "openmc/cell.h"
#include "openmc/material.h"

namespace openmc {

void MaterialFromFilter::get_all_bins(
  const Particle& p, TallyEstimator estimator, FilterMatch& match) const
{
  auto search = map_.find(p.material_last());
  if (search != map_.end()) {
    match.bins_.push_back(search->second);
    match.weights_.push_back(1.0);
  }
}

std::string MaterialFromFilter::text_label(int bin) const
{
  return "Material from " +
         std::to_string(model::materials[materials_[bin]]->id_);
}

} // namespace openmc
