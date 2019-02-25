#include "openmc/tallies/filter_material.h"

#include <sstream>

#include "openmc/capi.h"
#include "openmc/error.h"
#include "openmc/material.h"
#include "openmc/xml_interface.h"

namespace openmc {

void
MaterialFilter::from_xml(pugi::xml_node node)
{
  materials_ = get_node_array<int32_t>(node, "bins");
  n_bins_ = materials_.size();
}

void
MaterialFilter::initialize()
{
  // Convert material IDs to indices of the global array.
  for (auto& m : materials_) {
    auto search = model::material_map.find(m);
    if (search != model::material_map.end()) {
      m = search->second;
    } else {
      std::stringstream err_msg;
      err_msg << "Could not find material " << m
              << " specified on tally filter.";
      fatal_error(err_msg);
    }
  }

  // Populate the index->bin map.
  for (int i = 0; i < materials_.size(); i++) {
    map_[materials_[i]] = i;
  }
}

void
MaterialFilter::get_all_bins(const Particle* p, int estimator,
                             FilterMatch& match) const
{
  auto search = map_.find(p->material);
  if (search != map_.end()) {
    //TODO: off-by-one
    match.bins_.push_back(search->second + 1);
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
  //TODO: off-by-one
  return "Material " + std::to_string(model::materials[materials_[bin-1]]->id_);
}

//==============================================================================
// C-API functions
//==============================================================================

extern "C" int
openmc_material_filter_get_bins(int32_t index, int32_t** bins, int32_t* n)
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
  *bins = filt->materials_.data();
  *n = filt->materials_.size();
  return 0;
}

extern "C" int
openmc_material_filter_set_bins(int32_t index, int32_t n, const int32_t* bins)
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
  filt->materials_.clear();
  filt->materials_.resize(n);
  for (int i = 0; i < n; i++) filt->materials_[i] = bins[i];
  filt->n_bins_ = filt->materials_.size();
  filt->map_.clear();
  for (int i = 0; i < n; i++) filt->map_[filt->materials_[i]] = i;
  return 0;
}

} // namespace openmc
