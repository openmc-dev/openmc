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
  auto& mesh = *model::meshes[mesh_];

  int n_dim = mesh.n_dimension_;

  if (mesh.get_mesh_type() == "hexagonal") {
    // Hexagonal mesh with 2 x 4 surfaces
    return hexbin_label(bin);
  } else {
    // The regular behavior for a structured mesh with 2 x ndim surfaces
    return bin_label(n_dim, bin);
  }
}

std::string MeshSurfaceFilter::bin_label(int n_dim, int bin) const
{
  // Get flattend mesh index and surface index.
  int i_mesh = bin / (4 * n_dim);
  MeshDir surf_dir = static_cast<MeshDir>(bin % (4 * n_dim));

  // Get mesh index part of label.
  std::string out = MeshFilter::text_label(i_mesh);

  // Get surface part of label.
  switch (surf_dir) {
  case MeshDir::OUT_LEFT:
    out += " Outgoing, x-min";
    break;
  case MeshDir::IN_LEFT:
    out += " Incoming, x-min";
    break;
  case MeshDir::OUT_RIGHT:
    out += " Outgoing, x-max";
    break;
  case MeshDir::IN_RIGHT:
    out += " Incoming, x-max";
    break;
  case MeshDir::OUT_BACK:
    out += " Outgoing, y-min";
    break;
  case MeshDir::IN_BACK:
    out += " Incoming, y-min";
    break;
  case MeshDir::OUT_FRONT:
    out += " Outgoing, y-max";
    break;
  case MeshDir::IN_FRONT:
    out += " Incoming, y-max";
    break;
  case MeshDir::OUT_BOTTOM:
    out += " Outgoing, z-min";
    break;
  case MeshDir::IN_BOTTOM:
    out += " Incoming, z-min";
    break;
  case MeshDir::OUT_TOP:
    out += " Outgoing, z-max";
    break;
  case MeshDir::IN_TOP:
    out += " Incoming, z-max";
    break;
  }

  return out;
}

std::string MeshSurfaceFilter::hexbin_label(int bin) const
{

  // Get flattend mesh index and surface index.
  int i_mesh = bin / (4 * 4);
  HexMeshDir surf_dir = static_cast<HexMeshDir>(bin % (4 * 4));

  // Get mesh index part of label.
  std::string out = MeshFilter::text_label(i_mesh);

  switch (surf_dir) {
  case HexMeshDir::OUT_SE:
    out += " Outgoing, r_max SE";
    break;
  case HexMeshDir::IN_SE:
    out += " Incoming, r_max SE";
    break;
  case HexMeshDir::OUT_NW:
    out += " OUtgoing, r min NW";
    break;
  case HexMeshDir::IN_NW:
    out += " Incoming, r min NW";
    break;
  case HexMeshDir::OUT_E:
    out += " Outgoing, q max E";
    break;
  case HexMeshDir::IN_E:
    out += " Incoming, q max E";
    break;
  case HexMeshDir::OUT_W:
    out += " Outgoing, q min W";
    break;
  case HexMeshDir::IN_W:
    out += " Incoming, q min W";
    break;
  case HexMeshDir::OUT_NE:
    out += " Outgoing, s max NE";
    break;
  case HexMeshDir::IN_NE:
    out += " Incoming, s max NE";
    break;
  case HexMeshDir::OUT_SW:
    out += " Outgoing, s min SW";
    break;
  case HexMeshDir::IN_SW:
    out += " Incoming, s min SW";
    break;
  case HexMeshDir::OUT_BOTTOM:
    out += " Outgoing, z min";
    break;
  case HexMeshDir::IN_BOTTOM:
    out += " Incoming, z min";
    break;
  case HexMeshDir::OUT_TOP:
    out += " Outgoing, z max";
    break;
  case HexMeshDir::IN_TOP:
    out += " Incoming, z max";
    break;
  }

  return out;
}

void MeshSurfaceFilter::set_mesh(int32_t mesh)
{
  mesh_ = mesh;
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
