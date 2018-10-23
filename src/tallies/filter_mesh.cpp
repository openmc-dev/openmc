#include "openmc/tallies/filter_mesh.h"

#include <sstream>

#include "openmc/capi.h"
#include "openmc/constants.h"
#include "openmc/error.h"
#include "openmc/mesh.h"

namespace openmc {

void
MeshFilter::from_xml(pugi::xml_node node)
{
  auto bins_ = get_node_array<int32_t>(node, "bins");
  if (bins_.size() != 1) {
    fatal_error("Only one mesh can be specified per " + type()
                + " mesh filter.");
  }

  auto id = bins_[0];
  auto search = mesh_map.find(id);
  if (search != mesh_map.end()) {
    mesh_ = search->second;
  } else{
    std::stringstream err_msg;
    err_msg << "Could not find cell " << id << " specified on tally filter.";
    fatal_error(err_msg);
  }

  n_bins_ = 1;
  for (auto dim : meshes[mesh_]->shape_) n_bins_ *= dim;
}

void
MeshFilter::get_all_bins(Particle* p, int estimator, FilterMatch& match)
const
{
  if (estimator != ESTIMATOR_TRACKLENGTH) {
    auto bin = meshes[mesh_]->get_bin(p->coord[0].xyz);
    if (bin >= 0) {
      match.bins_.push_back(bin);
      match.weights_.push_back(1.0);
    }
  } else {
    meshes[mesh_]->bins_crossed(p, match.bins_, match.weights_);
  }
}

void
MeshFilter::to_statepoint(hid_t filter_group) const
{
  Filter::to_statepoint(filter_group);
  write_dataset(filter_group, "bins", meshes[mesh_]->id_);
}

std::string
MeshFilter::text_label(int bin) const
{
  auto& mesh = *meshes[mesh_];
  int n_dim = mesh.n_dimension_;

  int ijk[n_dim];
  mesh.get_indices_from_bin(bin, ijk);

  std::stringstream out;
  out << "Mesh Index (" << ijk[0];
  if (n_dim > 1) out << ", " << ijk[1];
  if (n_dim > 2) out << ", " << ijk[2];
  out << ")";

  return out.str();
}

//==============================================================================
// C-API functions
//==============================================================================

extern "C" int
openmc_mesh_filter_get_mesh(int32_t index, int32_t* index_mesh)
{
  int err = verify_filter(index);
  if (err) return err;

  auto filt = filter_from_f(index);
  if (filt->type() != "mesh" && filt->type() != "meshsurface") {
    set_errmsg("Tried to get mesh on a non-mesh filter.");
    return OPENMC_E_INVALID_TYPE;
  }

  auto mesh_filt = static_cast<MeshFilter*>(filt);
  *index_mesh = mesh_filt->mesh_;
  return 0;
}

extern "C" int
openmc_mesh_filter_set_mesh(int32_t index, int32_t index_mesh)
{
  int err = verify_filter(index);
  if (err) return err;

  auto filt = filter_from_f(index);
  if (filt->type() != "mesh" && filt->type() != "meshsurface") {
    set_errmsg("Tried to set mesh on a non-mesh filter.");
    return OPENMC_E_INVALID_TYPE;
  }

  if (index_mesh < 0 || index_mesh >= meshes.size()) {
    set_errmsg("Index in 'meshes' array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }

  auto mesh_filt = static_cast<MeshFilter*>(filt);
  mesh_filt->mesh_ = index_mesh;
  if (filt->type() == "mesh") {
    mesh_filt->n_bins_ = 1;
  } else {
    filt->n_bins_ = 4 * meshes[index_mesh]->n_dimension_;
  }
  for (auto dim : meshes[index_mesh]->shape_) mesh_filt->n_bins_ *= dim;
  filter_update_n_bins(index);
  return 0;
}

} // namespace openmc
