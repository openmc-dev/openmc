#include "openmc/random_ray/ray.h"
#include "openmc/geometry.h"
#include "openmc/settings.h"
#include "openmc/mgxs_interface.h"
#include "openmc/random_ray/source_region.h"
#include "openmc/source.h"

namespace openmc {

//==============================================================================
// Ray implementation
//==============================================================================

Ray::Ray()
{
  angular_flux_.resize(data::mg.num_energy_groups_);
  delta_psi_.resize(data::mg.num_energy_groups_);
}

Ray::Ray(uint64_t index_source, uint64_t nrays, int iter) : Ray::Ray()
{
  initialize_ray(index_source, nrays, iter);
}

uint64_t Ray::transport_history_based_single_ray(double distance_inactive, double distance_active)
{
  using namespace openmc;
  while (true) {
    if (!alive())
      break;
    event_advance_ray(distance_inactive, distance_active);
    if (!alive())
      break;
    event_cross_surface();
  }

  return n_event();
}

void Ray::event_advance_ray(double distance_inactive, double distance_active)
{
  // Find the distance to the nearest boundary
  boundary() = distance_to_boundary(*this);
  double distance = boundary().distance;
  assert(distance > 0.0);

  // Check for final termination
  if (is_active_) {
    if (distance_travelled_ + distance >= distance_active) {
      distance = distance_active - distance_travelled_;
      wgt() = 0.0;
    }
    distance_travelled_ += distance;
    attenuate_flux(distance, true);
  }

  // Check for end of inactive region (dead zone)
  if(!is_active_)
  {
    if(distance_travelled_ + distance >= distance_inactive)
    {
      is_active_ = true;
      double distance_dead = distance_inactive - distance_travelled_;
      attenuate_flux(distance_dead,  false);

      double distance_alive = distance - distance_dead;

      // Ensure we haven't travelled past the active phase as well
      if (distance_alive > distance_active)
      {
        distance_alive = distance_active;
        wgt() = 0.0;
      }
      attenuate_flux(distance_alive, true);

      distance_travelled_ = distance_alive;
    }
    else
    {
      distance_travelled_ += distance;
      attenuate_flux(distance,  false);
    }
  }

  // Advance particle
  for (int j = 0; j < n_coord(); ++j) {
    coord(j).r += distance * coord(j).u;
  }
}


void Ray::attenuate_flux(double distance, bool is_active)
{
  assert(distance > 0.0);

  int negroups = data::mg.num_energy_groups_;

  n_event()++;

  // Determine Cell Index etc.
  int i_cell = lowest_coord().cell;
  int64_t source_region_idx = random_ray::source_region_offsets[i_cell] + cell_instance();
  int64_t source_region_group_idx = source_region_idx * negroups;
  int material = this->material();

  for( int e = 0; e < negroups; e++ )
  {
    float Sigma_t = data::mg.macro_xs_[material].get_xs(MgxsType::TOTAL, e, NULL, NULL, NULL);
    float tau = Sigma_t * distance;
    float exponential = cjosey_exponential(tau);
    float new_delta_psi = (angular_flux_[e] - random_ray::source[source_region_group_idx + e]) * exponential;
    delta_psi_[e] = new_delta_psi;
    angular_flux_[e] -= new_delta_psi;
  }

  if(is_active)
  {
    random_ray::lock[source_region_idx].lock();

    for (int e = 0; e < negroups; e++) {
      random_ray::scalar_flux_new[source_region_group_idx + e] += delta_psi_[e];
    }

    if (random_ray::was_hit[source_region_idx] == 0) {
      random_ray::was_hit[source_region_idx] = 1;
    }

    random_ray::volume[source_region_idx] += distance;

    // Tally position if not done already
    if (!random_ray::position_recorded[source_region_idx]) {
      Position midpoint = r() + u() * (distance/2.0);
      random_ray::position[source_region_idx] = midpoint;
      random_ray::position_recorded[source_region_idx] = 1;
    }

    random_ray::lock[source_region_idx].unlock();
  }
}



void Ray::initialize_ray(uint64_t index_source, uint64_t nrays, int iter)
{
  id() = index_source;

  // Reset particle event counter
  n_event() = 0;

  if( settings::ray_distance_inactive <= 0.0 )
    is_active_ = true;
  else
    is_active_ = false;

  wgt() = 1.0;

  // set random number seed
  int64_t particle_seed = (iter-1) * nrays + id();
  init_particle_seeds(particle_seed, seeds());
  stream() = STREAM_TRACKING;

  // sample from external source distribution (should use box)
  auto site = sample_external_source(current_seed());
  from_source(&site);

  // If the cell hasn't been determined based on the particle's location,
  // initiate a search for the current cell. This generally happens at the
  // beginning of the history and again for any secondary particles
  if (lowest_coord().cell == C_NONE) {
    if (!exhaustive_find_cell(*this)) {
      this->mark_as_lost("Could not find the cell containing particle "
          + std::to_string(id()));
    }

    // Set birth cell attribute
    if (cell_born() == C_NONE) cell_born() = lowest_coord().cell;
  }

  // Initialize ray's starting angular flux to starting location's isotropic source
  int negroups = data::mg.num_energy_groups_;
  int i_cell = lowest_coord().cell;
  int64_t source_region_idx = random_ray::source_region_offsets[i_cell] + cell_instance();

  for (int e = 0; e < negroups; e++) {
    angular_flux_[e] = random_ray::source[source_region_idx * negroups + e];
  }
}

//==============================================================================
// Non-method functions
//==============================================================================

// returns 1 - exp(-tau)
// Equivalent to -(_expm1f(-tau)), but is stable
inline float cjosey_exponential(const float tau)
{
  const float c1n = -1.0000013559236386308f;
  const float c2n = 0.23151368626911062025f;
  const float c3n = -0.061481916409314966140f;
  const float c4n = 0.0098619906458127653020f;
  const float c5n = -0.0012629460503540849940f;
  const float c6n = 0.00010360973791574984608f;
  const float c7n = -0.000013276571933735820960f;

  const float c0d = 1.0f;
  const float c1d = -0.73151337729389001396f;
  const float c2d = 0.26058381273536471371f;
  const float c3d = -0.059892419041316836940f;
  const float c4d = 0.0099070188241094279067f;
  const float c5d = -0.0012623388962473160860f;
  const float c6d = 0.00010361277635498731388f;
  const float c7d = -0.000013276569500666698498f;

  float x = -tau;
  float num, den;

  den = c7d;
  den = den * x + c6d;
  den = den * x + c5d;
  den = den * x + c4d;
  den = den * x + c3d;
  den = den * x + c2d;
  den = den * x + c1d;
  den = den * x + c0d;

  num = c7n;
  num = num * x + c6n;
  num = num * x + c5n;
  num = num * x + c4n;
  num = num * x + c3n;
  num = num * x + c2n;
  num = num * x + c1n;
  num = num * x;

  const float exponential = num / den;
  return exponential;
}

} // namespace openmc
