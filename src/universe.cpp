#include "openmc/universe.h"

#include <set>

#include "openmc/hdf5_interface.h"
#include "openmc/particle.h"

namespace openmc {

namespace model {

std::unordered_map<int32_t, int32_t> universe_map;
vector<unique_ptr<Universe>> universes;

} // namespace model

//==============================================================================
// Universe implementation
//==============================================================================

void Universe::to_hdf5(hid_t universes_group) const
{
  // Create a group for this universe.
  auto group = create_group(universes_group, fmt::format("universe {}", id_));

  // Write the geometry representation type.
  write_string(group, "geom_type", "csg", false);

  // Write the contained cells.
  if (cells_.size() > 0) {
    vector<int32_t> cell_ids;
    for (auto i_cell : cells_)
      cell_ids.push_back(model::cells[i_cell]->id_);
    write_dataset(group, "cells", cell_ids);
  }

  close_group(group);
}

bool Universe::find_cell(GeometryState& p) const
{
  if (filled_with_triso_base_ != -1) {
    bool found = find_cell_in_virtual_lattice(p);
    if (found) {
      return found;
    }
  }
  const auto& cells {
    !partitioner_ ? cells_ : partitioner_->get_cells(p.r_local(), p.u_local())};

  Position r {p.r_local()};
  Position u {p.u_local()};
  auto surf = p.surface();
  int32_t i_univ = p.lowest_coord().universe();

  for (auto i_cell : cells) {
    if (model::cells[i_cell]->universe_ != i_univ)
      continue;
    // Check if this cell contains the particle
    if (model::cells[i_cell]->contains(r, u, surf)) {
      p.lowest_coord().cell() = i_cell;
      return true;
    }
  }
  return false;
}
bool Universe::find_cell_in_virtual_lattice(GeometryState& p) const
{
  Cell& c {*model::cells[model::cell_map[filled_with_triso_base_]]};
  vector<int> lat_ind(3);
  Position r {p.r_local()};
  lat_ind[0] = std::max<int>(
    std::min<int>(
      floor((r.x - c.vl_lower_left_[0]) / c.vl_pitch_[0]), c.vl_shape_[0] - 1),
    0);
  lat_ind[1] = std::max<int>(
    std::min<int>(
      floor((r.y - c.vl_lower_left_[1]) / c.vl_pitch_[1]), c.vl_shape_[1] - 1),
    0);
  lat_ind[2] = std::max<int>(
    std::min<int>(
      floor((r.z - c.vl_lower_left_[2]) / c.vl_pitch_[2]), c.vl_shape_[2] - 1),
    0);

  int32_t i_univ = p.lowest_coord().universe();
  for (int token :
    c.vl_triso_distribution_[lat_ind[0] + lat_ind[1] * c.vl_shape_[0] +
                             lat_ind[2] * c.vl_shape_[0] * c.vl_shape_[1]]) {
    vector<double> triso_center = model::surfaces[abs(token) - 1]->get_center();
    double triso_radius = model::surfaces[abs(token) - 1]->get_radius();
    if (model::cells
          [model::cell_map[model::surfaces[abs(token) - 1]->triso_base_index_]]
            ->universe_ != i_univ)
      continue;
    if (abs(token) == abs(p.surface())) {
      if (p.surface() < 0) {
        p.lowest_coord().cell() =
          model::cell_map[model::surfaces[abs(token) - 1]
                            ->triso_particle_index_];
        return true;
      } else {
        p.lowest_coord().cell() = model::cell_map[filled_with_triso_base_];
        return true;
      }
    }
    if (pow(r.x - triso_center[0], 2) + pow(r.y - triso_center[1], 2) +
          pow(r.z - triso_center[2], 2) <
        pow(triso_radius, 2)) {
      p.lowest_coord().cell() =
        model::cell_map[model::surfaces[abs(token) - 1]->triso_particle_index_];
      return true;
    }
  }
  if (model::cells[model::cell_map[filled_with_triso_base_]]->universe_ ==
      i_univ) {
    p.lowest_coord().cell() = model::cell_map[filled_with_triso_base_];
    return true;
  }
  return false;
}

BoundingBox Universe::bounding_box() const
{
  BoundingBox bbox = {INFTY, -INFTY, INFTY, -INFTY, INFTY, -INFTY};
  if (cells_.size() == 0) {
    return {};
  } else {
    for (const auto& cell : cells_) {
      auto& c = model::cells[cell];
      bbox |= c->bounding_box();
    }
  }
  return bbox;
}

//==============================================================================
// UniversePartitioner implementation
//==============================================================================

UniversePartitioner::UniversePartitioner(const Universe& univ)
{
  // Define an ordered set of surface indices that point to z-planes.  Use a
  // functor to to order the set by the z0_ values of the corresponding planes.
  struct compare_surfs {
    bool operator()(const int32_t& i_surf, const int32_t& j_surf) const
    {
      const auto* surf = model::surfaces[i_surf].get();
      const auto* zplane = dynamic_cast<const SurfaceZPlane*>(surf);
      double zi = zplane->z0_;
      surf = model::surfaces[j_surf].get();
      zplane = dynamic_cast<const SurfaceZPlane*>(surf);
      double zj = zplane->z0_;
      return zi < zj;
    }
  };
  std::set<int32_t, compare_surfs> surf_set;

  // Find all of the z-planes in this universe.  A set is used here for the
  // O(log(n)) insertions that will ensure entries are not repeated.
  for (auto i_cell : univ.cells_) {
    for (auto token : model::cells[i_cell]->surfaces()) {
      auto i_surf = std::abs(token) - 1;
      const auto* surf = model::surfaces[i_surf].get();
      if (const auto* zplane = dynamic_cast<const SurfaceZPlane*>(surf))
        surf_set.insert(i_surf);
    }
  }

  // Populate the surfs_ vector from the ordered set.
  surfs_.insert(surfs_.begin(), surf_set.begin(), surf_set.end());

  // Populate the partition lists.
  partitions_.resize(surfs_.size() + 1);
  for (auto i_cell : univ.cells_) {
    // It is difficult to determine the bounds of a complex cell, so add complex
    // cells to all partitions.
    if (!model::cells[i_cell]->is_simple()) {
      for (auto& p : partitions_)
        p.push_back(i_cell);
      continue;
    }

    // Find the tokens for bounding z-planes.
    int32_t lower_token = 0, upper_token = 0;
    double min_z, max_z;
    for (auto token : model::cells[i_cell]->surfaces()) {
      const auto* surf = model::surfaces[std::abs(token) - 1].get();
      if (const auto* zplane = dynamic_cast<const SurfaceZPlane*>(surf)) {
        if (lower_token == 0 || zplane->z0_ < min_z) {
          lower_token = token;
          min_z = zplane->z0_;
        }
        if (upper_token == 0 || zplane->z0_ > max_z) {
          upper_token = token;
          max_z = zplane->z0_;
        }
      }
    }

    // If there are no bounding z-planes, add this cell to all partitions.
    if (lower_token == 0) {
      for (auto& p : partitions_)
        p.push_back(i_cell);
      continue;
    }

    // Find the first partition this cell lies in.  If the lower_token indicates
    // a negative halfspace, then the cell is unbounded in the lower direction
    // and it lies in the first partition onward.  Otherwise, it is bounded by
    // the positive halfspace given by the lower_token.
    int first_partition = 0;
    if (lower_token > 0) {
      for (int i = 0; i < surfs_.size(); ++i) {
        if (lower_token == surfs_[i] + 1) {
          first_partition = i + 1;
          break;
        }
      }
    }

    // Find the last partition this cell lies in.  The logic is analogous to the
    // logic for first_partition.
    int last_partition = surfs_.size();
    if (upper_token < 0) {
      for (int i = first_partition; i < surfs_.size(); ++i) {
        if (upper_token == -(surfs_[i] + 1)) {
          last_partition = i;
          break;
        }
      }
    }

    // Add the cell to all relevant partitions.
    for (int i = first_partition; i <= last_partition; ++i) {
      partitions_[i].push_back(i_cell);
    }
  }
}

const vector<int32_t>& UniversePartitioner::get_cells(
  Position r, Direction u) const
{
  // Perform a binary search for the partition containing the given coordinates.
  int left = 0;
  int middle = (surfs_.size() - 1) / 2;
  int right = surfs_.size() - 1;
  while (true) {
    // Check the sense of the coordinates for the current surface.
    const auto& surf = *model::surfaces[surfs_[middle]];
    if (surf.sense(r, u)) {
      // The coordinates lie in the positive halfspace.  Recurse if there are
      // more surfaces to check.  Otherwise, return the cells on the positive
      // side of this surface.
      int right_leaf = right - (right - middle) / 2;
      if (right_leaf != middle) {
        left = middle + 1;
        middle = right_leaf;
      } else {
        return partitions_[middle + 1];
      }

    } else {
      // The coordinates lie in the negative halfspace.  Recurse if there are
      // more surfaces to check.  Otherwise, return the cells on the negative
      // side of this surface.
      int left_leaf = left + (middle - left) / 2;
      if (left_leaf != middle) {
        right = middle - 1;
        middle = left_leaf;
      } else {
        return partitions_[middle];
      }
    }
  }
}

} // namespace openmc
