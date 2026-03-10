#include "openmc/tallies/filter_meshsurface.h"

#include "openmc/capi.h"
#include "openmc/constants.h"
#include "openmc/error.h"
#include "openmc/mesh.h"

namespace openmc {

void MeshSurfaceFilter::get_all_bins(
  const Particle& p, TallyEstimator estimator, FilterMatch& match) const
{
  Position r0 = p.r_last_current();
  Position r1 = p.r();
  if (translated_) {
    r0 -= translation();
    r1 -= translation();
  }

  Direction u = p.u();
  model::meshes[mesh_]->surface_bins_crossed(r0, r1, u, match.bins_);
  for (auto b : match.bins_)
    match.weights_.push_back(1.0);
}

std::string MeshSurfaceFilter::text_label(int bin) const
{
  auto mesh = dynamic_cast<StructuredMesh*>(model::meshes[mesh_].get());
  int n_dim = mesh->n_dimension_;

  // Get flattend mesh index and surface index.
  int i_mesh = bin / (4 * n_dim);
  int i_surf = bin % (4 * n_dim);

  std::string out = MeshFilter::text_label(i_mesh);
  out += " ";
  out += mesh->surface_label(i_surf);
  return out;
}

void MeshSurfaceFilter::set_mesh(int32_t mesh)
{
  mesh_ = mesh;
  if (!dynamic_cast<StructuredMesh*>(model::meshes[mesh_].get()))
    fatal_error("Only structured mesh is supported in MeshSurfaceFilter.");
  n_bins_ = model::meshes[mesh_]->n_surface_bins();
}

//==============================================================================
// C-API functions
//==============================================================================

extern "C" int openmc_meshsurface_filter_get_mesh(
  int32_t index, int32_t* index_mesh)
{
  return openmc_mesh_filter_get_mesh(index, index_mesh);
}

extern "C" int openmc_meshsurface_filter_set_mesh(
  int32_t index, int32_t index_mesh)
{
  return openmc_mesh_filter_set_mesh(index, index_mesh);
}

extern "C" int openmc_meshsurface_filter_get_translation(
  int32_t index, double translation[3])
{
  return openmc_mesh_filter_get_translation(index, translation);
}

extern "C" int openmc_meshsurface_filter_set_translation(
  int32_t index, double translation[3])
{
  return openmc_mesh_filter_set_translation(index, translation);
}

} // namespace openmc
