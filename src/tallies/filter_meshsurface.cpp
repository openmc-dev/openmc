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

  // Get flattend mesh index and surface index.
  int i_mesh = bin / (4 * n_dim);

  // Get mesh index part of label.
  std::string out = MeshFilter::text_label(i_mesh);
  auto surf_dir = (bin % (4 * n_dim));
  if (mesh->mesh_type != HexagonalMesh::mesh_type) {
    // Get surface part of label.
    switch (surf_dir) {
    case 0:
      out += fmt::format(" Outgoing, {}-min", mesh->labels_[0]);
      break;
    case 1:
      out += fmt::format(" Incoming, {}-min", mesh->labels_[0]);
      break;
    case 2:
      out += fmt::format(" Outgoing, {}-max", mesh->labels_[0]);
      break;
    case 3:
      out += fmt::format(" Incoming, {}-max", mesh->labels_[0]);
      break;
    case 4:
      out += fmt::format(" Outgoing, {}-min", mesh->labels_[1]);
      break;
    case 5:
      out += fmt::format(" Incoming, {}-min", mesh->labels_[1]);
      break;
    case 6:
      out += fmt::format(" Outgoing, {}-max", mesh->labels_[1]);
      break;
    case 7:
      out += fmt::format(" Incoming, {}-max", mesh->labels_[1]);
      break;
    case 8:
      out += fmt::format(" Outgoing, {}-min", mesh->labels_[2]);
      break;
    case 9:
      out += fmt::format(" Incoming, {}-min", mesh->labels_[2]);
      break;
    case 10:
      out += fmt::format(" Outgoing, {}-max", mesh->labels_[2]);
      break;
    case 11:
      out += fmt::format(" Incoming, {}-max", mesh->labels_[2]);
      break;
    }
  } else {
    // Get surface part of label.
    switch (surf_dir) {
    case 0:
      out += fmt::format(
        " Outgoing, {}-min {}-max", mesh->labels_[0], mesh->labels_[1]);
      break;
    case 1:
      out += fmt::format(
        " Incoming, {}-min {}-max", mesh->labels_[0], mesh->labels_[1]);
      break;
    case 2:
      out += fmt::format(
        " Outgoing, {}-max {}-min", mesh->labels_[0], mesh->labels_[1]);
      break;
    case 3:
      out += fmt::format(
        " Incoming, {}-max {}-min", mesh->labels_[0], mesh->labels_[1]);
      break;
    case 4:
      out += fmt::format(
        " Outgoing, {}-min {}-max", mesh->labels_[1], mesh->labels_[2]);
      break;
    case 5:
      out += fmt::format(
        " Incoming, {}-min {}-max", mesh->labels_[1], mesh->labels_[2]);
      break;
    case 6:
      out += fmt::format(
        " Outgoing, {}-max {}-min", mesh->labels_[1], mesh->labels_[2]);
      break;
    case 7:
      out += fmt::format(
        " Incoming, {}-max {}-min", mesh->labels_[1], mesh->labels_[2]);
      break;
    case 8:
      out += fmt::format(
        " Outgoing, {}-min {}-max", mesh->labels_[2], mesh->labels_[0]);
      break;
    case 9:
      out += fmt::format(
        " Incoming, {}-min {}-max", mesh->labels_[2], mesh->labels_[0]);
      break;
    case 10:
      out += fmt::format(
        " Outgoing, {}-max {}-min", mesh->labels_[2], mesh->labels_[0]);
      break;
    case 11:
      out += fmt::format(
        " Incoming, {}-max {}-min", mesh->labels_[2], mesh->labels_[0]);
      break;
    case 12:
      out += fmt::format(" Outgoing, {}-min", mesh->labels_[3]);
      break;
    case 13:
      out += fmt::format(" Incoming, {}-min", mesh->labels_[3]);
      break;
    case 14:
      out += fmt::format(" Outgoing, {}-max", mesh->labels_[3]);
      break;
    case 15:
      out += fmt::format(" Incoming, {}-max", mesh->labels_[3]);
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
