#include "openmc/tallies/filter_materialfrom.h"

#include "openmc/cell.h"
#include "openmc/constants.h"
#include "openmc/material.h"

namespace openmc {

void MaterialFromFilter::get_all_bins(
  const Particle& p, TallyEstimator estimator, FilterMatch& match) const
{
  // Derive the material the particle came from rather than tracking it
  // separately, so that this filter and CellFromFilter can never disagree
  // about where a score originated. The material fill lives on the cell at the
  // lowest coordinate level; the instance is needed because a cell with
  // distributed materials resolves a different material per instance.
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
