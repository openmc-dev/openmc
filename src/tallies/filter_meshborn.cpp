#include "openmc/tallies/filter_meshborn.h"

#include "openmc/capi.h"
#include "openmc/constants.h"
#include "openmc/error.h"
#include "openmc/mesh.h"

namespace openmc {

void MeshBornFilter::get_all_bins(
  const Particle& p, TallyEstimator estimator, FilterMatch& match) const
{
  Position r_born = p.r_born();

  // apply translation if present
  if (translated_) {
    r_born -= translation();
  }

  auto bin = model::meshes[mesh_]->get_bin(r_born);
  if (bin >= 0) {
    match.bins_.push_back(bin);
    match.weights_.push_back(1.0);
  }
}

std::string MeshBornFilter::text_label(int bin) const
{
  auto& mesh = *model::meshes.at(mesh_);
  return mesh.bin_label(bin) + " (born)";
}

//==============================================================================
// C-API functions
//==============================================================================

extern "C" int openmc_meshborn_filter_get_mesh(
  int32_t index, int32_t* index_mesh)
{
  return openmc_mesh_filter_get_mesh(index, index_mesh);
}

extern "C" int openmc_meshborn_filter_set_mesh(
  int32_t index, int32_t index_mesh)
{
  return openmc_mesh_filter_set_mesh(index, index_mesh);
}

extern "C" int openmc_meshborn_filter_get_translation(
  int32_t index, double translation[3])
{
  return openmc_mesh_filter_get_translation(index, translation);
}

extern "C" int openmc_meshborn_filter_set_translation(
  int32_t index, double translation[3])
{
  return openmc_mesh_filter_set_translation(index, translation);
}

} // namespace openmc
