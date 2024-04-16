#include "openmc/random_ray/random_ray.h"

#include "openmc/geometry.h"
#include "openmc/message_passing.h"
#include "openmc/mgxs_interface.h"
#include "openmc/random_ray/flat_source_domain.h"
#include "openmc/search.h"
#include "openmc/settings.h"
#include "openmc/simulation.h"
#include "openmc/source.h"

namespace openmc {

//==============================================================================
// Non-method functions
//==============================================================================

// returns 1 - exp(-tau)
// Equivalent to -(_expm1f(-tau)), but faster
// Written by Colin Josey.
float cjosey_exponential(float tau)
{
  constexpr float c1n = -1.0000013559236386308f;
  constexpr float c2n = 0.23151368626911062025f;
  constexpr float c3n = -0.061481916409314966140f;
  constexpr float c4n = 0.0098619906458127653020f;
  constexpr float c5n = -0.0012629460503540849940f;
  constexpr float c6n = 0.00010360973791574984608f;
  constexpr float c7n = -0.000013276571933735820960f;

  constexpr float c0d = 1.0f;
  constexpr float c1d = -0.73151337729389001396f;
  constexpr float c2d = 0.26058381273536471371f;
  constexpr float c3d = -0.059892419041316836940f;
  constexpr float c4d = 0.0099070188241094279067f;
  constexpr float c5d = -0.0012623388962473160860f;
  constexpr float c6d = 0.00010361277635498731388f;
  constexpr float c7d = -0.000013276569500666698498f;

  float x = -tau;

  float den = c7d;
  den = den * x + c6d;
  den = den * x + c5d;
  den = den * x + c4d;
  den = den * x + c3d;
  den = den * x + c2d;
  den = den * x + c1d;
  den = den * x + c0d;

  float num = c7n;
  num = num * x + c6n;
  num = num * x + c5n;
  num = num * x + c4n;
  num = num * x + c3n;
  num = num * x + c2n;
  num = num * x + c1n;
  num = num * x;

  return num / den;
}

//==============================================================================
// RandomRay implementation
//==============================================================================

// Static Variable Declarations
double RandomRay::distance_inactive_;
double RandomRay::distance_active_;
unique_ptr<Source> RandomRay::ray_source_;

RandomRay::RandomRay()
  : negroups_(data::mg.num_energy_groups_),
    angular_flux_(data::mg.num_energy_groups_),
    delta_psi_(data::mg.num_energy_groups_)
{}

RandomRay::RandomRay(uint64_t ray_id, FlatSourceDomain* domain) : RandomRay()
{
  initialize_ray(ray_id, domain);
}

// Transports ray until termination criteria are met
uint64_t RandomRay::transport_history_based_single_ray()
{
  using namespace openmc;
  while (alive()) {
    event_advance_ray();
    if (!alive())
      break;
    event_cross_surface();
  }

  return n_event();
}

// Transports ray across a single source region
void RandomRay::event_advance_ray()
{
  // Find the distance to the nearest boundary
  boundary() = distance_to_boundary(*this);
  double distance = boundary().distance;

  if (distance <= 0.0) {
    mark_as_lost("Negative transport distance detected for particle " +
                 std::to_string(id()));
    return;
  }

  if (is_active_) {
    // If the ray is in the active length, need to check if it has
    // reached its maximum termination distance. If so, reduce
    // the ray traced length so that the ray does not overrun the
    // maximum numerical length (so as to avoid numerical bias).
    if (distance_travelled_ + distance >= distance_active_) {
      distance = distance_active_ - distance_travelled_;
      wgt() = 0.0;
    }

    distance_travelled_ += distance;
    attenuate_flux(distance, true);
  } else {
    // If the ray is still in the dead zone, need to check if it
    // has entered the active phase. If so, split into two segments (one
    // representing the final part of the dead zone, the other representing the
    // first part of the active length) and attenuate each. Otherwise, if the
    // full length of the segment is within the dead zone, attenuate as normal.
    if (distance_travelled_ + distance >= distance_inactive_) {
      is_active_ = true;
      double distance_dead = distance_inactive_ - distance_travelled_;
      attenuate_flux(distance_dead, false);

      double distance_alive = distance - distance_dead;

      // Ensure we haven't travelled past the active phase as well
      if (distance_alive > distance_active_) {
        distance_alive = distance_active_;
        wgt() = 0.0;
      }

      attenuate_flux(distance_alive, true);
      distance_travelled_ = distance_alive;
    } else {
      distance_travelled_ += distance;
      attenuate_flux(distance, false);
    }
  }

  // Advance particle
  for (int j = 0; j < n_coord(); ++j) {
    coord(j).r += distance * coord(j).u;
  }
}

// This function forms the inner loop of the random ray transport process.
// It is responsible for several tasks. Based on the incoming angular flux
// of the ray and the source term in the region, the outgoing angular flux
// is computed. The delta psi between the incoming and outgoing fluxes is
// contributed to the estimate of the total scalar flux in the source region.
// Additionally, the contribution of the ray path to the stochastically
// estimated volume is also kept track of. All tasks involving writing
// to the data for the source region are done with a lock over the entire
// source region.  Locks are used instead of atomics as all energy groups
// must be written, such that locking once is typically much more efficient
// than use of many atomic operations corresponding to each energy group
// individually (at least on CPU). Several other bookkeeping tasks are also
// performed when inside the lock.
void RandomRay::attenuate_flux(double distance, bool is_active)
{
  // The number of geometric intersections is counted for reporting purposes
  n_event()++;

  // Determine source region index etc.
  int i_cell = lowest_coord().cell;

  // The source region is the spatial region index
  int64_t source_region =
    domain_->source_region_offsets_[i_cell] + cell_instance();

  // The source element is the energy-specific region index
  int64_t source_element = source_region * negroups_;
  int material = this->material();

  // Temperature and angle indices, if using multiple temperature
  // data sets and/or anisotropic data sets.
  // TODO: Currently assumes we are only using single temp/single
  // angle data.
  const int t = 0;
  const int a = 0;

  // MOC incoming flux attenuation + source contribution/attenuation equation
  for (int g = 0; g < negroups_; g++) {
    float sigma_t = data::mg.macro_xs_[material].get_xs(
      MgxsType::TOTAL, g, NULL, NULL, NULL, t, a);
    float tau = sigma_t * distance;
    float exponential = cjosey_exponential(tau); // exponential = 1 - exp(-tau)
    float new_delta_psi =
      (angular_flux_[g] - domain_->source_[source_element + g]) * exponential;
    delta_psi_[g] = new_delta_psi;
    angular_flux_[g] -= new_delta_psi;
  }

  // If ray is in the active phase (not in dead zone), make contributions to
  // source region bookkeeping
  if (is_active) {

    // Aquire lock for source region
    domain_->lock_[source_region].lock();

    // Accumulate delta psi into new estimate of source region flux for
    // this iteration
    for (int g = 0; g < negroups_; g++) {
      domain_->scalar_flux_new_[source_element + g] += delta_psi_[g];
    }

    // If the source region hasn't been hit yet this iteration,
    // indicate that it now has
    if (domain_->was_hit_[source_region] == 0) {
      domain_->was_hit_[source_region] = 1;
    }

    // Accomulate volume (ray distance) into this iteration's estimate
    // of the source region's volume
    domain_->volume_[source_region] += distance;

    // Tally valid position inside the source region (e.g., midpoint of
    // the ray) if not done already
    if (!domain_->position_recorded_[source_region]) {
      Position midpoint = r() + u() * (distance / 2.0);
      domain_->position_[source_region] = midpoint;
      domain_->position_recorded_[source_region] = 1;
    }

    // Release lock
    domain_->lock_[source_region].unlock();
  }
}

void RandomRay::initialize_ray(uint64_t ray_id, FlatSourceDomain* domain)
{
  domain_ = domain;

  // Reset particle event counter
  n_event() = 0;

  is_active_ = (distance_inactive_ <= 0.0);

  wgt() = 1.0;

  // set identifier for particle
  id() = simulation::work_index[mpi::rank] + ray_id;

  // set random number seed
  int64_t particle_seed =
    (simulation::current_batch - 1) * settings::n_particles + id();
  init_particle_seeds(particle_seed, seeds());
  stream() = STREAM_TRACKING;

  // Sample from ray source distribution
  SourceSite site {ray_source_->sample(current_seed())};
  site.E = lower_bound_index(
    data::mg.rev_energy_bins_.begin(), data::mg.rev_energy_bins_.end(), site.E);
  site.E = negroups_ - site.E - 1.;
  this->from_source(&site);

  // Locate ray
  if (lowest_coord().cell == C_NONE) {
    if (!exhaustive_find_cell(*this)) {
      this->mark_as_lost(
        "Could not find the cell containing particle " + std::to_string(id()));
    }

    // Set birth cell attribute
    if (cell_born() == C_NONE)
      cell_born() = lowest_coord().cell;
  }

  // Initialize ray's starting angular flux to starting location's isotropic
  // source
  int i_cell = lowest_coord().cell;
  int64_t source_region_idx =
    domain_->source_region_offsets_[i_cell] + cell_instance();

  for (int g = 0; g < negroups_; g++) {
    angular_flux_[g] = domain_->source_[source_region_idx * negroups_ + g];
  }
}

} // namespace openmc
