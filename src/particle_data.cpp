#include "openmc/particle_data.h"

#include "openmc/geometry.h"
#include "openmc/nuclide.h"
#include "openmc/photon.h"
#include "openmc/settings.h"
#include "openmc/tallies/derivative.h"
#include "openmc/tallies/filter.h"
#include "openmc/tallies/tally.h"

namespace openmc {

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

ParticleData::ParticleData()
{
  // Create and clear coordinate levels
  coord_.resize(model::n_coord_levels);
  cell_last_.resize(model::n_coord_levels);
  clear();

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
}

} // namespace openmc
