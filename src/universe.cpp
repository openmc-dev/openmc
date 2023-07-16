#include "openmc/universe.h"
#include "openmc/partitioners.h"
#include "openmc/hdf5_interface.h"
#include "openmc/timer.h"
#include "openmc/error.h"

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

bool Universe::find_cell(Particle& p) const
{

  const auto& cells {
    !partitioner_ ? cells_ : partitioner_->get_cells(p.r_local(), p.u_local())};

  bool result = find_cell_in_list(cells, p);

#ifdef PARTITIONER_FALLBACK_ENABLED
  if(!result && partitioner_) {
    result = find_cell_in_list(
      partitioner_->get_cells_fallback(p.r_local(), p.u_local()), p);
  }
  #endif

  return result;
}

bool Universe::find_cell_in_list(
  const std::vector<int>& cells, Particle& p) const
{
  for (auto it = cells.begin(); it != cells.end(); it++) {
    int32_t i_cell = *it;
    int32_t i_univ = p.coord(p.n_coord() - 1).universe;
    if (model::cells[i_cell]->universe_ != i_univ)
      continue;

    // Check if this cell contains the particle;
    Position r {p.r_local()};
    Direction u {p.u_local()};
    auto surf = p.surface();
    if (model::cells[i_cell]->contains(r, u, surf)) {
      p.coord(p.n_coord() - 1).cell = i_cell;
      return true;
    }
  }

  return false;
}

int Universe::find_cell_for_point(const std::vector<int>& cells_to_search, const Position& r) const {
  for (int i = 0; i < cells_to_search.size(); i++) {
    int32_t i_cell = cells_to_search[i];
    // Check if this cell contains the particle;
    Direction u {0.0, 0.0, 0.0};
    int32_t surf = false;
    if (model::cells[i_cell]->contains(r, u, surf)) {
      return i_cell;
    }
  }

  return -1;
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

double max_ignore_infinities(double orig, double val)
{
  if (val == INFTY || val == -INFTY) {
    return orig;
  } else {
    return std::max(orig, val);
  }
}

double min_ignore_infinities(double orig, double val)
{
  if (val == INFTY || val == -INFTY) {
    return orig;
  } else {
    return std::min(orig, val);
  }
}

AABB Universe::partitioner_bounding_box() const
{
  AABB box;
  if (cells_.size() == 0) {
    box.max_ = {0, 0, 0};
    box.min_ = {0, 0, 0};
  } else {
    for (const auto& cell : cells_) {
      auto bbox = model::cells[cell]->bounding_box();
      // extend on each axis, ignoring any infinities

      box.max_.x = max_ignore_infinities(box.max_.x, bbox.xmax);
      box.max_.y = max_ignore_infinities(box.max_.y, bbox.ymax);
      box.max_.z = max_ignore_infinities(box.max_.z, bbox.zmax);

      box.min_.x = min_ignore_infinities(box.min_.x, bbox.xmin);
      box.min_.y = min_ignore_infinities(box.min_.y, bbox.ymin);
      box.min_.z = min_ignore_infinities(box.min_.z, bbox.zmin);
    }

    write_message(
      "Paritioner bounding box ranges from (" + std::to_string(box.min_.x) +
        ", " + std::to_string(box.min_.y) + ", " + std::to_string(box.min_.z) +
        ") to (" + std::to_string(box.max_.x) + ", " +
        std::to_string(box.max_.y) + ", " + std::to_string(box.max_.z) + ")",
      5);
  }

  return box;
}

//! We do nothing here since it is a dummy destructor
UniversePartitioner::~UniversePartitioner() {}

const std::vector<int32_t>& UniversePartitioner::get_cells_fallback(Position r, Direction u) const {
  fatal_error("Fallback is not enabled for this partitioner. Please recompile with the PARTITIONER_FALLBACK_ENABLED macro disabled in partitioners.h.");
}

void UniversePartitioner::export_to_file(const std::string& path) const {
  warning("Failed to export partitioner to " + path +
          " because current partitioner type currently does not support "
          "exporting as a file.");
}


} // namespace openmc
