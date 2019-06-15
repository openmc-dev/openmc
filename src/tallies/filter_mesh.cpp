#include "openmc/tallies/filter_mesh.h"

#include <sstream>

#include "openmc/capi.h"
#include "openmc/constants.h"
#include "openmc/error.h"
#include "openmc/mesh.h"
#include "openmc/xml_interface.h"

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
  auto search = model::mesh_map.find(id);
  if (search != model::mesh_map.end()) {
    set_mesh(search->second);
  } else{
    std::stringstream err_msg;
    err_msg << "Could not find mesh " << id << " specified on tally filter.";
    fatal_error(err_msg);
  }
}

void
MeshFilter::get_all_bins(const Particle* p, int estimator, FilterMatch& match)
const
{
  if (estimator != ESTIMATOR_TRACKLENGTH) {
    auto bin = model::meshes[mesh_]->get_bin(p->r());
    if (bin >= 0) {
      match.bins_.push_back(bin);
      match.weights_.push_back(1.0);
    }
  } else {
    model::meshes[mesh_]->bins_crossed(p, match.bins_, match.weights_);
  }
}

void
MeshFilter::to_statepoint(hid_t filter_group) const
{
  Filter::to_statepoint(filter_group);
  write_dataset(filter_group, "bins", model::meshes[mesh_]->id_);
}

std::string
MeshFilter::text_label(int bin) const
{
  auto& mesh = *model::meshes[mesh_];
  int n_dim = mesh.n_dimension_;

  std::vector<int> ijk(n_dim);
  mesh.get_indices_from_bin(bin, ijk.data());

  std::stringstream out;
  out << "Mesh Index (" << ijk[0];
  if (n_dim > 1) out << ", " << ijk[1];
  if (n_dim > 2) out << ", " << ijk[2];
  out << ")";

  return out.str();
}

void
MeshFilter::set_mesh(int32_t mesh)
{
  mesh_ = mesh;
  n_bins_ = model::meshes[mesh_]->n_bins();
}

//==============================================================================
// C-API functions
//==============================================================================

extern "C" int
openmc_mesh_filter_get_mesh(int32_t index, int32_t* index_mesh)
{
  // Make sure this is a valid index to an allocated filter.
  if (int err = verify_filter(index)) return err;

  // Get a pointer to the filter and downcast.
  const auto& filt_base = model::tally_filters[index].get();
  auto* filt = dynamic_cast<MeshFilter*>(filt_base);

  // Check the filter type.
  if (!filt) {
    set_errmsg("Tried to get mesh on a non-mesh filter.");
    return OPENMC_E_INVALID_TYPE;
  }

  // Output the mesh.
  *index_mesh = filt->mesh();
  return 0;
}

extern "C" int
openmc_mesh_filter_set_mesh(int32_t index, int32_t index_mesh)
{
  // Make sure this is a valid index to an allocated filter.
  if (int err = verify_filter(index)) return err;

  // Get a pointer to the filter and downcast.
  const auto& filt_base = model::tally_filters[index].get();
  auto* filt = dynamic_cast<MeshFilter*>(filt_base);

  // Check the filter type.
  if (!filt) {
    set_errmsg("Tried to set mesh on a non-mesh filter.");
    return OPENMC_E_INVALID_TYPE;
  }

  // Check the mesh index.
  if (index_mesh < 0 || index_mesh >= model::meshes.size()) {
    set_errmsg("Index in 'meshes' array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }

  // Update the filter.
  filt->set_mesh(index_mesh);
  return 0;
}

} // namespace openmc
