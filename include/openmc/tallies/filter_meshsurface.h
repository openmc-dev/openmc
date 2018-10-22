#ifndef OPENMC_TALLIES_FILTER_MESHSURFACE_H
#define OPENMC_TALLIES_FILTER_MESHSURFACE_H

#include <cstdint>
#include <sstream>

#include "openmc/constants.h"
#include "openmc/error.h"
#include "openmc/mesh.h"
#include "openmc/tallies/filter_mesh.h"


namespace openmc {

class MeshSurfaceFilter : public MeshFilter
{
public:
  std::string type() const override {return "meshsurface";}

  void
  from_xml(pugi::xml_node node) override
  {
    MeshFilter::from_xml(node);
    n_bins_ = 4 * meshes[mesh_]->n_dimension_;;
    for (auto dim : meshes[mesh_]->shape_) n_bins_ *= dim;
  }

  void
  get_all_bins(Particle* p, int estimator, FilterMatch& match)
  const override
  {
    meshes[mesh_]->surface_bins_crossed(p, match.bins_);
    for (auto b : match.bins_) match.weights_.push_back(1.0);
  }

  std::string
  text_label(int bin) const override
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
};

} // namespace openmc
#endif // OPENMC_TALLIES_FILTER_MESHSURFACE_H
