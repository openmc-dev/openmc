#include "openmc/particle_data.h"

#include "openmc/geometry.h"
#include "openmc/nuclide.h"
#include "openmc/photon.h"
#include "openmc/tallies/derivative.h"
#include "openmc/tallies/filter.h"
#include "openmc/tallies/tally.h"

namespace openmc {

ParticleData::ParticleData()
{
  // Create and clear coordinate levels
  coord_.resize(model::n_coord_levels);
  cell_last_.resize(model::n_coord_levels);
  clear();

  zero_delayed_bank();

  // Every particle starts with no accumulated flux derivative.
  if (!model::active_tallies.empty()) {
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
