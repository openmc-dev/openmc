#include "openmc/tallies/filter_meshsurface.h"

#include "openmc/capi.h"
#include "openmc/constants.h"
#include "openmc/error.h"
#include "openmc/mesh.h"

namespace openmc {

void MeshAngularFilter::get_all_bins(
  const Particle& p, TallyEstimator estimator, FilterMatch& match) const
{
  Direction u = p.u();
  if (!rotation_.empty()) {
    u = u.rotate(rotation_);
  }
  auto bin = model::meshes[mesh_]->get_bin(r);
  if (bin >= 0) {
    match.bins_.push_back(bin);
    match.weights_.push_back(1.0);
  }
}

void MeshAngularFilter::to_statepoint(hid_t filter_group) const
{
  Filter::to_statepoint(filter_group);
  write_dataset(filter_group, "bins", model::meshes[mesh_]->id_);
  if (rotated_) {
    write_dataset(filter_group, "rotation", rotation_);
  }
}

//==============================================================================
// C-API functions
//==============================================================================

extern "C" int openmc_meshangular_filter_get_mesh(
  int32_t index, int32_t* index_mesh)
{
  return openmc_mesh_filter_get_mesh(index, index_mesh);
}

extern "C" int openmc_meshangular_filter_set_mesh(
  int32_t index, int32_t index_mesh)
{
  return openmc_mesh_filter_set_mesh(index, index_mesh);
}

} // namespace openmc
