#include "openmc/tallies/filter_materialfrom.h"

#include "openmc/cell.h"
#include "openmc/constants.h"
#include "openmc/material.h"

namespace openmc {

void MaterialFromFilter::get_all_bins(
  const Particle& p, TallyEstimator estimator, FilterMatch& match) const
{
  int32_t i_cell = p.cell_last(p.n_coord_last() - 1);
  if (i_cell == C_NONE)
    return;

  int32_t i_material = model::cells[i_cell]->material(p.cell_instance_last());

  auto search = map_.find(i_material);
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
