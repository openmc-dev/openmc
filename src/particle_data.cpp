#include "openmc/particle_data.h"

#include <sstream>

#include "openmc/cell.h"
#include "openmc/error.h"
#include "openmc/geometry.h"
#include "openmc/material.h"
#include "openmc/nuclide.h"
#include "openmc/photon.h"
#include "openmc/settings.h"
#include "openmc/tallies/derivative.h"
#include "openmc/tallies/filter.h"
#include "openmc/tallies/tally.h"

namespace openmc {

void GeometryState::mark_as_lost(const char* message)
{
  fatal_error(message);
}

void GeometryState::mark_as_lost(const std::string& message)
{
  mark_as_lost(message.c_str());
}

void GeometryState::mark_as_lost(const std::stringstream& message)
{
  mark_as_lost(message.str());
}

void LocalCoord::rotate(const vector<double>& rotation)
{
  r = r.rotate(rotation);
  u = u.rotate(rotation);
  rotated = true;
}

void LocalCoord::reset()
{
  cell = C_NONE;
  universe = C_NONE;
  lattice = C_NONE;
  lattice_i[0] = 0;
  lattice_i[1] = 0;
  lattice_i[2] = 0;
  rotated = false;
}

GeometryState::GeometryState()
{
  // Create and clear coordinate levels
  coord_.resize(model::n_coord_levels);
  cell_last_.resize(model::n_coord_levels);
  clear();
}

void
GeometryState::advance_to_boundary_from_void() {
  auto root_coord = this->coord(0);
  const auto& root_universe = model::universes[model::root_universe];
  boundary().reset();

  for (auto c_i : root_universe->cells_) {
    auto dist =
      model::cells.at(c_i)->distance(root_coord.r, root_coord.u, 0, this);
    if (dist.first < boundary().distance) {
      boundary().distance = dist.first;
      boundary().surface_index = dist.second;
    }
  }

  // if no intersection or near-infinite intersection, reset
  // boundary information
  if (boundary().distance > 1e300) {
    boundary().distance = INFTY;
    boundary().surface_index = -1;
    return;
  }

  // move the particle up to (and just past) the boundary
  move_distance(boundary().distance + TINY_BIT);
}

void GeometryState::move_distance(double length)
{
  for (int j = 0; j < n_coord(); ++j) {
    coord(j).r += length * coord(j).u;
  }
}

ParticleData::ParticleData()
{
  zero_delayed_bank();

  // Every particle starts with no accumulated flux derivative.  Note that in
  // event mode, we construct the particle once up front, so have to run this
  // even if the current batch is inactive.
  if (!model::active_tallies.empty() || settings::event_based) {
    flux_derivs_.resize(model::tally_derivs.size());
    zero_flux_derivs();
  }

  // Allocate space for tally filter matches
  filter_matches_.resize(model::tally_filters.size());

  // Create microscopic cross section caches
  neutron_xs_.resize(data::nuclides.size());
  photon_xs_.resize(data::elements.size());

  // Creates the pulse-height storage for the particle
  if (!model::pulse_height_cells.empty()) {
    pht_storage_.resize(model::pulse_height_cells.size(), 0.0);
  }
}

TrackState ParticleData::get_track_state() const
{
  TrackState state;
  state.r = this->r();
  state.u = this->u();
  state.E = this->E();
  state.time = this->time();
  state.wgt = this->wgt();
  state.cell_id = model::cells[this->lowest_coord().cell]->id_;
  state.cell_instance = this->cell_instance();
  if (this->material() != MATERIAL_VOID) {
    state.material_id = model::materials[material()]->id_;
  }
  return state;
}

} // namespace openmc
