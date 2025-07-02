#include "openmc/tallies/filter_meshmaterial.h"

#include <utility> // for move

#include <fmt/core.h>

#include "openmc/capi.h"
#include "openmc/constants.h"
#include "openmc/container_util.h"
#include "openmc/error.h"
#include "openmc/material.h"
#include "openmc/mesh.h"
#include "openmc/xml_interface.h"

namespace openmc {

void MeshMaterialFilter::from_xml(pugi::xml_node node)
{
  // Get mesh ID
  auto mesh = get_node_array<int32_t>(node, "mesh");
  if (mesh.size() != 1) {
    fatal_error(
      "Only one mesh can be specified per " + type_str() + " mesh filter.");
  }

  auto id = mesh[0];
  auto search = model::mesh_map.find(id);
  if (search == model::mesh_map.end()) {
    fatal_error(
      fmt::format("Could not find mesh {} specified on tally filter.", id));
  }
  set_mesh(search->second);

  // Get pairs of (element index, material) and set the bins
  auto bins = get_node_array<int32_t>(node, "bins");
  this->set_bins(bins);

  if (check_for_node(node, "translation")) {
    set_translation(get_node_array<double>(node, "translation"));
  }
}

void MeshMaterialFilter::set_bins(span<int32_t> bins)
{
  if (bins.size() % 2 != 0) {
    fatal_error(
      fmt::format("Size of mesh material bins is not even: {}", bins.size()));
  }

  // Create a vector of ElementMat pairs from the flat vector of bins
  vector<ElementMat> element_mats;
  for (int64_t i = 0; i < bins.size() / 2; ++i) {
    int32_t element = bins[2 * i];
    int32_t mat_id = bins[2 * i + 1];
    auto search = model::material_map.find(mat_id);
    if (search == model::material_map.end()) {
      fatal_error(fmt::format(
        "Could not find material {} specified on tally filter.", mat_id));
    }
    int32_t mat_index = search->second;
    element_mats.push_back({element, mat_index});
  }

  this->set_bins(std::move(element_mats));
}

void MeshMaterialFilter::set_bins(vector<ElementMat>&& bins)
{
  // Swap internal bins_ with the provided vector to avoid copying
  bins_.swap(bins);

  // Clear and update the mapping and vector of materials
  materials_.clear();
  map_.clear();
  for (std::size_t i = 0; i < bins_.size(); ++i) {
    const auto& x = bins_[i];
    assert(x.index_mat >= 0);
    assert(x.index_mat < model::materials.size());
    materials_.insert(x.index_mat);
    map_[x] = i;
  }

  n_bins_ = bins_.size();
}

void MeshMaterialFilter::set_mesh(int32_t mesh)
{
  // perform any additional perparation for mesh tallies here
  mesh_ = mesh;
  model::meshes[mesh_]->prepare_for_point_location();
}

void MeshMaterialFilter::set_translation(const Position& translation)
{
  translated_ = true;
  translation_ = translation;
}

void MeshMaterialFilter::set_translation(const double translation[3])
{
  this->set_translation({translation[0], translation[1], translation[2]});
}

void MeshMaterialFilter::get_all_bins(
  const Particle& p, TallyEstimator estimator, FilterMatch& match) const
{
  // If current material is not in any bins, don't bother checking
  if (!contains(materials_, p.material())) {
    return;
  }

  Position last_r = p.r_last();
  Position r = p.r();
  Position u = p.u();

  // apply translation if present
  if (translated_) {
    last_r -= translation();
    r -= translation();
  }

  if (estimator != TallyEstimator::TRACKLENGTH) {
    int32_t index_element = model::meshes[mesh_]->get_bin(r);
    if (index_element >= 0) {
      auto search = map_.find({index_element, p.material()});
      if (search != map_.end()) {
        match.bins_.push_back(search->second);
        match.weights_.push_back(1.0);
      }
    }
  } else {
    // First determine which elements the particle crosses (may or may not
    // actually match bins so we have to adjust bins_/weight_ after)
    int32_t n_start = match.bins_.size();
    model::meshes[mesh_]->bins_crossed(
      last_r, r, u, match.bins_, match.weights_);
    int32_t n_end = match.bins_.size();

    // Go through bins and weights and check which ones are actually a match
    // based on the (element, material) pair. For matches, overwrite the bin.
    int i = 0;
    for (int j = n_start; j < n_end; ++j) {
      int32_t index_element = match.bins_[j];
      double weight = match.weights_[j];
      auto search = map_.find({index_element, p.material()});
      if (search != map_.end()) {
        match.bins_[n_start + i] = search->second;
        match.weights_[n_start + i] = weight;
        ++i;
      }
    }

    // Resize the vectors to remove the unmatched bins
    match.bins_.resize(n_start + i);
  }
}

void MeshMaterialFilter::to_statepoint(hid_t filter_group) const
{
  Filter::to_statepoint(filter_group);
  write_dataset(filter_group, "mesh", model::meshes[mesh_]->id_);

  size_t n = bins_.size();
  xt::xtensor<size_t, 2> data({n, 2});
  for (int64_t i = 0; i < n; ++i) {
    const auto& x = bins_[i];
    data(i, 0) = x.index_element;
    data(i, 1) = model::materials[x.index_mat]->id_;
  }
  write_dataset(filter_group, "bins", data);

  if (translated_) {
    write_dataset(filter_group, "translation", translation_);
  }
}

std::string MeshMaterialFilter::text_label(int bin) const
{
  auto& x = bins_[bin];
  auto& mesh = *model::meshes.at(mesh_);
  return fmt::format("Mesh {}, {}, Material {}", mesh.id(),
    mesh.bin_label(x.index_element), model::materials[x.index_mat]->id_);
}

} // namespace openmc
