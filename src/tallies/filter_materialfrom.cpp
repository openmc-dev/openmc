#include "openmc/tallies/filter_materialfrom.h"

#include "openmc/material.h"

namespace openmc {

void MaterialFromFilter::get_all_bins(
  const Particle& p, TallyEstimator estimator, FilterMatch& match) const
{
  for (int i = 0; i < p.n_coord_last(); i++) {
    if (p.cell_last(i).material_.size() == 1) {
      auto search = map_.find(p.cell_last(i).material_[0]);
      if (search != map_.end()) {
        match.bins_.push_back(search->second);
        match.weights_.push_back(1.0);
      }
    }
  }
}

std::string MaterialFromFilter::text_label(int bin) const
{
  return "Material from " +
         std::to_string(model::materials[material_[bin]]->id_);
}

} // namespace openmc
