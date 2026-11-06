#include "openmc/tallies/filter_importance.h"

#include <fmt/core.h>

#include "openmc/capi.h"
#include "openmc/search.h"
#include "openmc/settings.h"
#include "openmc/mesh.h"
#include "openmc/xml_interface.h"

namespace openmc {

//==============================================================================
// ImportanceFilter implementation
//==============================================================================

void ImportanceFilter::from_xml(pugi::xml_node node)
{
  // Set the importance (adjoint flux)
  auto importance = get_node_array<double>(node, "importance");
  this->set_importance(importance);

  // Set the importance mesh
  auto bins_ = get_node_array<int32_t>(node, "mesh");
  if (bins_.size() != 1) {
    fatal_error("Only one mesh can be specified per " + type_str()
                + " importance filter.");
  }

  auto id = bins_[0];
  auto search = model::mesh_map.find(id);
  if (search != model::mesh_map.end()) {
    set_mesh(search->second);
  } else{
    fatal_error(fmt::format(
      "Could not find mesh {} specified on tally filter.", id));
  }
}

void ImportanceFilter::set_importance(gsl::span<const double> importance)
{
  // Clear existing importance
  importance_.clear();
  importance_.reserve(importance.size());

  // Copy importance
  for (gsl::index i = 0; i < importance.size(); ++i) {
    importance_.push_back(importance[i]);
  }

  n_bins_ = 1; // This is probably redundant, I set this in set_mesh()
}

void ImportanceFilter::get_all_bins(
  const Particle& p, TallyEstimator estimator, FilterMatch& match) const
{
  auto bin = model::meshes[mesh_]->get_bin(p.r());
  if (bin >= 0) {
    match.bins_.push_back(0); // There is only one bin, (is it 0-indexed?)
    match.weights_.push_back(importance_[bin]); // the weight is the importance in the mesh bin
  }
}

void ImportanceFilter::to_statepoint(hid_t filter_group) const
{
  Filter::to_statepoint(filter_group);
  write_dataset(filter_group, "importance", importance_);
  write_dataset(filter_group, "mesh", model::meshes[mesh_]->id_);
}

std::string ImportanceFilter::text_label(int bin) const
{
  // returns the label of the importance mesh
  auto& mesh = *model::meshes[mesh_];
  return mesh.bin_label(bin);
}

void ImportanceFilter::set_mesh(int32_t mesh)
{
  mesh_ = mesh;
  n_bins_ = 1;  // We only have 1 bin in this filter
}

//==============================================================================
// C-API functions
//==============================================================================

extern"C" int
openmc_importance_filter_get_importance(int32_t index, const double** importance, size_t* n)
{
  // Make sure this is a valid index to an allocated filter.
  if (int err = verify_filter(index)) return err;

  // Get a pointer to the filter and downcast.
  const auto& filt_base = model::tally_filters[index].get();
  auto* filt = dynamic_cast<ImportanceFilter*>(filt_base);

  // Check the filter type.
  if (!filt) {
    set_errmsg("Tried to get importance bins on a non-importance filter.");
    return OPENMC_E_INVALID_TYPE;
  }

  // Output the importance.
  *importance = filt->importance().data();
  *n = filt->importance().size();
  return 0;
}

extern "C" int
openmc_importance_filter_set_importance(int32_t index, size_t n, const double* importance)
{
  // Make sure this is a valid index to an allocated filter.
  if (int err = verify_filter(index)) return err;

  // Get a pointer to the filter and downcast.
  const auto& filt_base = model::tally_filters[index].get();
  auto* filt = dynamic_cast<ImportanceFilter*>(filt_base);

  // Check the filter type.
  if (!filt) {
    set_errmsg("Tried to set importance bins on a non-importance filter.");
    return OPENMC_E_INVALID_TYPE;
  }

  // Update the filter.
  filt->set_importance({importance, n});
  return 0;
}

extern "C" int
openmc_importance_filter_get_mesh(int32_t index, int32_t* index_mesh)
{
  // Make sure this is a valid index to an allocated filter.
  if (int err = verify_filter(index)) return err;

  // Get a pointer to the filter and downcast.
  const auto& filt_base = model::tally_filters[index].get();
  auto* filt = dynamic_cast<ImportanceFilter*>(filt_base);

  // Check the filter type.
  if (!filt) {
    set_errmsg("Tried to get importance mesh on a non-importance filter.");
    return OPENMC_E_INVALID_TYPE;
  }

  // Output the mesh.
  *index_mesh = filt->mesh();
  return 0;
}

extern "C" int
openmc_importance_filter_set_mesh(int32_t index, int32_t index_mesh)
{
  // Make sure this is a valid index to an allocated filter.
  if (int err = verify_filter(index)) return err;

  // Get a pointer to the filter and downcast.
  const auto& filt_base = model::tally_filters[index].get();
  auto* filt = dynamic_cast<ImportanceFilter*>(filt_base);

  // Check the filter type.
  if (!filt) {
    set_errmsg("Tried to set importance mesh on a non-importance filter.");
    return OPENMC_E_INVALID_TYPE;
  }

  // Check the mesh index.
  if (index_mesh < 0 || index_mesh >= model::meshes.size()) {
    set_errmsg("Index in 'meshes' array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }

  // Update the filter.
  filt->set_mesh(index_mesh);
  return 0;
}

}// namespace openmc