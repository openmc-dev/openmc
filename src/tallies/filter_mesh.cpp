#include "openmc/tallies/filter_mesh.h"

#include <fmt/core.h>

#include "openmc/capi.h"
#include "openmc/constants.h"
#include "openmc/error.h"
#include "openmc/mesh.h"
#include "openmc/xml_interface.h"

namespace openmc {

void MeshFilter::from_xml(pugi::xml_node node)
{
  auto bins_ = get_node_array<int32_t>(node, "bins");
  if (bins_.size() != 1) {
    fatal_error(
      "Only one mesh can be specified per " + type_str() + " mesh filter.");
  }

  auto id = bins_[0];
  auto search = model::mesh_map.find(id);
  if (search != model::mesh_map.end()) {
    set_mesh(search->second);
  } else {
    fatal_error(
      fmt::format("Could not find mesh {} specified on tally filter.", id));
  }

  if (check_for_node(node, "translation")) {
    set_translation(get_node_array<double>(node, "translation"));
  }
}

void MeshFilter::get_all_bins(
  const Particle& p, TallyEstimator estimator, FilterMatch& match) const
{

  Position last_r = p.r_last();
  Position r = p.r();
  Position u = p.u();

  // apply translation if present
  if (translated_) {
    last_r -= translation();
    r -= translation();
  }

  if (estimator != TallyEstimator::TRACKLENGTH) {
    auto bin = model::meshes[mesh_]->get_bin(r);
    if (bin >= 0) {
      match.bins_.push_back(bin);
      match.weights_.push_back(1.0);
    }
  } else {
    model::meshes[mesh_]->bins_crossed(
      last_r, r, u, match.bins_, match.weights_);
  }
}

void MeshFilter::to_statepoint(hid_t filter_group) const
{
  Filter::to_statepoint(filter_group);
  write_dataset(filter_group, "bins", model::meshes[mesh_]->id_);
  if (translated_) {
    write_dataset(filter_group, "translation", translation_);
  }
}

std::string MeshFilter::text_label(int bin) const
{
  auto& mesh = *model::meshes.at(mesh_);
  std::string label = mesh.bin_label(bin);
  return label;
}

void MeshFilter::set_mesh(int32_t mesh)
{
  // perform any additional preparation for mesh tallies here
  mesh_ = mesh;
  n_bins_ = model::meshes[mesh_]->n_bins();
  model::meshes[mesh_]->prepare_for_point_location();
}

void MeshFilter::set_translation(const Position& translation)
{
  translated_ = true;
  translation_ = translation;
}

void MeshFilter::set_translation(const double translation[3])
{
  this->set_translation({translation[0], translation[1], translation[2]});
}

//==============================================================================
// C-API functions
//==============================================================================

extern "C" int openmc_mesh_filter_get_mesh(int32_t index, int32_t* index_mesh)
{
  if (!index_mesh) {
    set_errmsg("Mesh index argument is a null pointer.");
    return OPENMC_E_INVALID_ARGUMENT;
  }

  // Make sure this is a valid index to an allocated filter.
  if (int err = verify_filter(index))
    return err;

  // Get a pointer to the filter and downcast.
  const auto& filt_base = model::tally_filters[index].get();
  auto* filt = dynamic_cast<MeshFilter*>(filt_base);

  // Check the filter type.
  if (!filt) {
    set_errmsg("Tried to get mesh on a non-mesh filter.");
    return OPENMC_E_INVALID_TYPE;
  }

  // Output the mesh.
  *index_mesh = filt->mesh();
  return 0;
}

extern "C" int openmc_mesh_filter_set_mesh(int32_t index, int32_t index_mesh)
{
  // Make sure this is a valid index to an allocated filter.
  if (int err = verify_filter(index))
    return err;

  // Get a pointer to the filter and downcast.
  const auto& filt_base = model::tally_filters[index].get();
  auto* filt = dynamic_cast<MeshFilter*>(filt_base);

  // Check the filter type.
  if (!filt) {
    set_errmsg("Tried to set mesh on a non-mesh filter.");
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

extern "C" int openmc_mesh_filter_get_translation(
  int32_t index, double translation[3])
{
  // Make sure this is a valid index to an allocated filter
  if (int err = verify_filter(index))
    return err;

  // Check the filter type
  const auto& filter = model::tally_filters[index];
  if (filter->type() != FilterType::MESH &&
      filter->type() != FilterType::MESHBORN &&
      filter->type() != FilterType::MESH_SURFACE) {
    set_errmsg("Tried to get a translation from a non-mesh-based filter.");
    return OPENMC_E_INVALID_TYPE;
  }

  // Get translation from the mesh filter and set value
  auto mesh_filter = dynamic_cast<MeshFilter*>(filter.get());
  const auto& t = mesh_filter->translation();
  for (int i = 0; i < 3; i++) {
    translation[i] = t[i];
  }

  return 0;
}

extern "C" int openmc_mesh_filter_set_translation(
  int32_t index, double translation[3])
{
  // Make sure this is a valid index to an allocated filter
  if (int err = verify_filter(index))
    return err;

  const auto& filter = model::tally_filters[index];
  // Check the filter type
  if (filter->type() != FilterType::MESH &&
      filter->type() != FilterType::MESHBORN &&
      filter->type() != FilterType::MESH_SURFACE) {
    set_errmsg("Tried to set mesh on a non-mesh-based filter.");
    return OPENMC_E_INVALID_TYPE;
  }

  // Get a pointer to the filter and downcast
  auto mesh_filter = dynamic_cast<MeshFilter*>(filter.get());

  // Set the translation
  mesh_filter->set_translation(translation);

  return 0;
}

} // namespace openmc
