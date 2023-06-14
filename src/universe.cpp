#include "openmc/universe.h"

#include <set>

#include "openmc/hdf5_interface.h"
#include "openmc/timer.h"
#include <iostream>
#include <mutex>
#include <time.h>

double finding_time = 0.0;

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

std::mutex finding_time_sum_mutex;

bool Universe::find_cell(Particle& p) const
{
  Timer t;
  t.start();

  const auto& cells {
    !partitioner_ ? cells_ : partitioner_->get_cells(p.r_local(), p.u_local())};

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
      finding_time_sum_mutex.lock();
      finding_time += t.elapsed();
      finding_time_sum_mutex.unlock();
      return true;
    }
  }

  finding_time_sum_mutex.lock();
  finding_time += t.elapsed();
  finding_time_sum_mutex.unlock();
  return false;
}

// not complete
bool Universe::find_cell_in_list(std::vector<int> cells_to_search, std::vector<int> cells_found, Position& r) const
{
  Timer t;
  t.start();

  for (int32_t i_cell : cells_to_search) {
    // Check if this cell contains the particle;
    Position r {p.r_local()};
    Direction u {p.u_local()};
    auto surf = p.surface();
    if (model::cells[i_cell]->contains(r, u, surf)) {
      p.coord(p.n_coord() - 1).cell = i_cell;
      finding_time_sum_mutex.lock();
      finding_time += t.elapsed();
      finding_time_sum_mutex.unlock();
      return true;
    }
  }

  finding_time_sum_mutex.lock();
  finding_time += t.elapsed();
  finding_time_sum_mutex.unlock();
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

//! We do nothing here since it is a dummy destructor
UniversePartitioner::~UniversePartitioner() {}



} // namespace openmc
