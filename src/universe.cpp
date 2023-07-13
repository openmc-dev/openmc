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

  bool result = find_cell_in_list(p, cells);

  #ifdef PARTITIONER_FALLBACK_ENABLED
  if(!result && partitioner_) {
    result = find_cell_in_list(p, partitioner_->get_cells_fallback(p.r_local(), p.u_local()));
  }
  #endif

  return result;
}

bool Universe::find_cell_in_list(Particle& p, const std::vector<int>& cells) const {
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

void Universe::find_cell_in_list(const std::vector<int>& cells_to_search, std::vector<int>& cells_found, std::vector<bool>& skip_cell, Position& r) const
{
  for (int i = 0; i < cells_to_search.size(); i++) {
    if(skip_cell[i]) continue;
    int32_t i_cell = cells_to_search[i];
    // Check if this cell contains the particle;
    Direction u {0.0, 0.0, 0.0};
    int32_t surf = false;
    if (model::cells[i_cell]->contains(r, u, surf)) {
      //std::cout << "Found " << i_cell << "\n";
      cells_found.push_back(i_cell);
      skip_cell[i] = true;
    }
  }
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

int Universe::find_cell_for_point(const Position& r) const {
  return find_cell_for_point(cells_, r);
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

//! We do nothing here since it is a dummy destructor
UniversePartitioner::~UniversePartitioner() {}

const std::vector<int32_t>& UniversePartitioner::get_cells_fallback(Position r, Direction u) const {
  fatal_error("Fallback is not enabled for this partitioner. Please recompile with the PARTITIONER_FALLBACK_ENABLED macro disabled in partitioners.h.");
}

void UniversePartitioner::export_to_file(const std::string& path) const {
  warning("Failed to export partitioner to " + path 
        + " because current partitioner type currently does not support exporting as a file.");
  // no need to abort
}


} // namespace openmc
