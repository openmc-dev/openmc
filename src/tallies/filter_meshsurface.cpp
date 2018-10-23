#include "openmc/tallies/filter_meshsurface.h"

#include "openmc/capi.h"
#include "openmc/constants.h"
#include "openmc/error.h"
#include "openmc/mesh.h"

namespace openmc {

void
MeshSurfaceFilter::from_xml(pugi::xml_node node)
{
  MeshFilter::from_xml(node);
  n_bins_ = 4 * meshes[mesh_]->n_dimension_;;
  for (auto dim : meshes[mesh_]->shape_) n_bins_ *= dim;
}

void
MeshSurfaceFilter::get_all_bins(Particle* p, int estimator, FilterMatch& match)
const
{
  meshes[mesh_]->surface_bins_crossed(p, match.bins_);
  for (auto b : match.bins_) match.weights_.push_back(1.0);
}

std::string
MeshSurfaceFilter::text_label(int bin) const
{
  auto& mesh = *meshes[mesh_];
  int n_dim = mesh.n_dimension_;

  // Get flattend mesh index and surface index.
  int i_mesh = (bin - 1) / (4 * n_dim) + 1;
  int i_surf = ((bin - 1) % (4 * n_dim)) + 1;

  // Get mesh index part of label.
  std::string out = MeshFilter::text_label(i_mesh);

  // Get surface part of label.
  if (i_surf == OUT_LEFT) {
    out += " Outgoing, x-min";
  } else if (i_surf == IN_LEFT) {
    out += " Incoming, x-min";
  } else if (i_surf == OUT_RIGHT) {
    out += " Outgoing, x-max";
  } else if (i_surf == IN_RIGHT) {
    out += " Incoming, x-max";
  } else if (i_surf == OUT_BACK) {
    out += " Outgoing, y-min";
  } else if (i_surf == IN_BACK) {
    out += " Incoming, y-min";
  } else if (i_surf == OUT_FRONT) {
    out += " Outgoing, y-max";
  } else if (i_surf == IN_FRONT) {
    out += " Incoming, y-max";
  } else if (i_surf == OUT_BOTTOM) {
    out += " Outgoing, z-min";
  } else if (i_surf == IN_BOTTOM) {
    out += " Incoming, z-min";
  } else if (i_surf == OUT_TOP) {
    out += " Outgoing, z-max";
  } else if (i_surf == IN_TOP) {
    out += " Incoming, z-max";
  }

  return out;
}

//==============================================================================
// C-API functions
//==============================================================================

extern"C" int
openmc_meshsurface_filter_get_mesh(int32_t index, int32_t* index_mesh)
{return openmc_mesh_filter_get_mesh(index, index_mesh);}

extern"C" int
openmc_meshsurface_filter_set_mesh(int32_t index, int32_t index_mesh)
{return openmc_mesh_filter_set_mesh(index, index_mesh);}

} // namespace openmc
