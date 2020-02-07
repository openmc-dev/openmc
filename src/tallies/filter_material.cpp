#include "openmc/tallies/filter_material.h"

#include <fmt/core.h>

#include "openmc/capi.h"
#include "openmc/material.h"
#include "openmc/xml_interface.h"

namespace openmc {

void
MaterialFilter::from_xml(pugi::xml_node node)
{
  // Get material IDs and convert to indices in the global materials vector
  auto mats = get_node_array<int32_t>(node, "bins");
  for (auto& m : mats) {
    auto search = model::material_map.find(m);
    if (search == model::material_map.end()) {
      throw std::runtime_error{fmt::format(
        "Could not find material {} specified on tally filter.", m)};
    }
    m = search->second;
  }

  this->set_materials(mats);
}

void
MaterialFilter::set_materials(gsl::span<const int32_t> materials)
{
  // Clear existing materials
  materials_.clear();
  materials_.reserve(materials.size());
  map_.clear();

  // Update materials and mapping
  for (auto& index : materials) {
    Expects(index >= 0);
    Expects(index < model::materials.size());
    materials_.push_back(index);
    map_[index] = materials_.size() - 1;
  }

  n_bins_ = materials_.size();
}

void
MaterialFilter::get_all_bins(const Particle* p, TallyEstimator estimator,
                             FilterMatch& match) const
{
  auto search = map_.find(p->material_);
  if (search != map_.end()) {
    match.bins_.push_back(search->second);
    match.weights_.push_back(1.0);
  }
}

void
MaterialFilter::to_statepoint(hid_t filter_group) const
{
  Filter::to_statepoint(filter_group);
  std::vector<int32_t> material_ids;
  for (auto c : materials_) material_ids.push_back(model::materials[c]->id_);
  write_dataset(filter_group, "bins", material_ids);
}

std::string
MaterialFilter::text_label(int bin) const
{
  return fmt::format("Material {}", model::materials[materials_[bin]]->id_);
}

//==============================================================================
// C-API functions
//==============================================================================

extern "C" int
openmc_material_filter_get_bins(int32_t index, const int32_t** bins, size_t* n)
{
  // Make sure this is a valid index to an allocated filter.
  if (int err = verify_filter(index)) return err;

  // Get a pointer to the filter and downcast.
  const auto& filt_base = model::tally_filters[index].get();
  auto* filt = dynamic_cast<MaterialFilter*>(filt_base);

  // Check the filter type.
  if (!filt) {
    set_errmsg("Tried to get material filter bins on a non-material filter.");
    return OPENMC_E_INVALID_TYPE;
  }

  // Output the bins.
  *bins = filt->materials().data();
  *n = filt->materials().size();
  return 0;
}

extern "C" int
openmc_material_filter_set_bins(int32_t index, size_t n, const int32_t* bins)
{
  // Make sure this is a valid index to an allocated filter.
  if (int err = verify_filter(index)) return err;

  // Get a pointer to the filter and downcast.
  const auto& filt_base = model::tally_filters[index].get();
  auto* filt = dynamic_cast<MaterialFilter*>(filt_base);

  // Check the filter type.
  if (!filt) {
    set_errmsg("Tried to set material filter bins on a non-material filter.");
    return OPENMC_E_INVALID_TYPE;
  }

  // Update the filter.
  filt->set_materials({bins, n});
  return 0;
}

} // namespace openmc
