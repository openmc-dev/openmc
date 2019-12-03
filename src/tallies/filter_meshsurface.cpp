#include "openmc/tallies/filter_meshsurface.h"

#include "openmc/capi.h"
#include "openmc/constants.h"
#include "openmc/error.h"
#include "openmc/mesh.h"

namespace openmc {

void
MeshSurfaceFilter::get_all_bins(const Particle* p, int estimator,
                                FilterMatch& match) const
{
  model::meshes[mesh_]->surface_bins_crossed(p, match.bins_);
  for (auto b : match.bins_) match.weights_.push_back(1.0);
}

std::string
MeshSurfaceFilter::text_label(int bin) const
{
  auto& mesh = *model::meshes[mesh_];
  int n_dim = mesh.n_dimension_;

  // Get flattend mesh index and surface index.
  int i_mesh = bin / (4 * n_dim);
  int i_surf = (bin % (4 * n_dim)) + 1;

  // Get mesh index part of label.
  std::string out = MeshFilter::text_label(i_mesh);

  // Get surface part of label.
  switch (i_surf) {
    case OUT_LEFT:
      out += " Outgoing, x-min";
      break;
    case IN_LEFT:
      out += " Incoming, x-min";
      break;
    case OUT_RIGHT:
      out += " Outgoing, x-max";
      break;
    case IN_RIGHT:
      out += " Incoming, x-max";
      break;
    case OUT_BACK:
      out += " Outgoing, y-min";
      break;
    case IN_BACK:
      out += " Incoming, y-min";
      break;
    case OUT_FRONT:
      out += " Outgoing, y-max";
      break;
    case IN_FRONT:
      out += " Incoming, y-max";
      break;
    case OUT_BOTTOM:
      out += " Outgoing, z-min";
      break;
    case IN_BOTTOM:
      out += " Incoming, z-min";
      break;
    case OUT_TOP:
      out += " Outgoing, z-max";
      break;
    case IN_TOP:
      out += " Incoming, z-max";
      break;
  }

  return out;
}

void
MeshSurfaceFilter::set_mesh(int32_t mesh)
{
  mesh_ = mesh;
  n_bins_ = model::meshes[mesh_]->n_surface_bins();
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
