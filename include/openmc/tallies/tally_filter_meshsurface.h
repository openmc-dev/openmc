#ifndef OPENMC_TALLY_FILTER_MESHSURFACE_H
#define OPENMC_TALLY_FILTER_MESHSURFACE_H

#include <cstdint>
#include <sstream>

#include "openmc/constants.h"
#include "openmc/error.h"
#include "openmc/mesh.h"
#include "openmc/tallies/tally_filter.h"


namespace openmc {

class MeshSurfaceFilter : public TallyFilter
{
public:
  virtual ~MeshSurfaceFilter() override = default;

  virtual void
  from_xml(pugi::xml_node node) override
  {
    auto bins_ = get_node_array<int32_t>(node, "bins");
    if (bins_.size() != 1)
      fatal_error("Only one mesh can be specified per mesh filter.");

    auto id = bins_[0];
    auto search = mesh_map.find(id);
    if (search != mesh_map.end()) {
      mesh_ = search->second;
    } else{
      std::stringstream err_msg;
      err_msg << "Could not find cell " << id << " specified on tally filter.";
      fatal_error(err_msg);
    }

    n_bins_ = 4 * meshes[mesh_]->n_dimension_;;
    for (auto dim : meshes[mesh_]->shape_) n_bins_ *= dim;
  }

  virtual void
  get_all_bins(Particle* p, int estimator, TallyFilterMatch& match)
  const override
  {
    meshes[mesh_]->surface_bins_crossed(p, match.bins);
    for (auto b : match.bins) match.weights.push_back(1.0);
  }

  virtual void
  to_statepoint(hid_t filter_group) const override
  {
    write_dataset(filter_group, "type", "meshsurface");
    write_dataset(filter_group, "n_bins", n_bins_);
    write_dataset(filter_group, "bins", meshes[mesh_]->id_);
  }

  virtual std::string
  text_label(int bin) const override
  {
    auto& mesh = *meshes[mesh_];
    int n_dim = mesh.n_dimension_;

    // Get flattend mesh index and surface index.
    int i_mesh = (bin - 1) / (4 * n_dim) + 1;
    int i_surf = ((bin - 1) % (4 * n_dim)) + 1;

    // Get mesh index part of label.
    int ijk[n_dim];
    mesh.get_indices_from_bin(i_mesh, ijk);
    std::stringstream out;
    out << "Mesh Index (" << ijk[0];
    if (n_dim > 1) out << ", " << ijk[1];
    if (n_dim > 2) out << ", " << ijk[2];
    out << ")";

    // Get surface part of label.
    if (i_surf == OUT_LEFT) {
      out << " Outgoing, x-min";
    } else if (i_surf == IN_LEFT) {
      out << " Incoming, x-min";
    } else if (i_surf == OUT_RIGHT) {
      out << " Outgoing, x-max";
    } else if (i_surf == IN_RIGHT) {
      out << " Incoming, x-max";
    } else if (i_surf == OUT_BACK) {
      out << " Outgoing, y-min";
    } else if (i_surf == IN_BACK) {
      out << " Incoming, y-min";
    } else if (i_surf == OUT_FRONT) {
      out << " Outgoing, y-max";
    } else if (i_surf == IN_FRONT) {
      out << " Incoming, y-max";
    } else if (i_surf == OUT_BOTTOM) {
      out << " Outgoing, z-min";
    } else if (i_surf == IN_BOTTOM) {
      out << " Incoming, z-min";
    } else if (i_surf == OUT_TOP) {
      out << " Outgoing, z-max";
    } else if (i_surf == IN_TOP) {
      out << " Incoming, z-max";
    }

    return out.str();
  }

  int32_t mesh_;
};

} // namespace openmc
#endif // OPENMC_TALLY_FILTER_MESHSURFACE_H
