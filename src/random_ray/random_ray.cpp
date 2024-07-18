#include "openmc/random_ray/random_ray.h"

#include "openmc/constants.h"
#include "openmc/geometry.h"
#include "openmc/message_passing.h"
#include "openmc/mgxs_interface.h"
#include "openmc/random_ray/flat_source_domain.h"
#include "openmc/random_ray/linear_source_domain.h"
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

// The below two functions (exponentialG and exponentialG2) were developed
// by Colin Josey. The implementation of these functions is closely based
// on the OpenMOC versions of these functions. The OpenMOC license is given
// below:

// Copyright (C) 2012-2023 Massachusetts Institute of Technology and OpenMOC
// contributors
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

// Computes y = 1/x-(1-exp(-x))/x**2 using a 5/6th order rational
// approximation. It is accurate to 2e-7 over [0, 1e5]. Developed by Colin
// Josey using Remez's algorithm, with original implementation in OpenMOC at:
// https://github.com/mit-crpg/OpenMOC/blob/develop/src/exponentials.h
float exponentialG(float tau)
{
  // Numerator coefficients in rational approximation for 1/x - (1 - exp(-x)) /
  // x^2
  constexpr float d0n = 0.5f;
  constexpr float d1n = 0.176558112351595f;
  constexpr float d2n = 0.04041584305811143f;
  constexpr float d3n = 0.006178333902037397f;
  constexpr float d4n = 0.0006429894635552992f;
  constexpr float d5n = 0.00006064409107557148f;

  // Denominator coefficients in rational approximation for 1/x - (1 - exp(-x))
  // / x^2
  constexpr float d0d = 1.0f;
  constexpr float d1d = 0.6864462055546078f;
  constexpr float d2d = 0.2263358514260129f;
  constexpr float d3d = 0.04721469893686252f;
  constexpr float d4d = 0.006883236664917246f;
  constexpr float d5d = 0.0007036272419147752f;
  constexpr float d6d = 0.00006064409107557148f;

  float x = tau;

  float num = d5n;
  num = num * x + d4n;
  num = num * x + d3n;
  num = num * x + d2n;
  num = num * x + d1n;
  num = num * x + d0n;

  float den = d6d;
  den = den * x + d5d;
  den = den * x + d4d;
  den = den * x + d3d;
  den = den * x + d2d;
  den = den * x + d1d;
  den = den * x + d0d;

  return num / den;
}

// Computes G2 : y = 2/3 - (1 + 2/x) * (1/x + 0.5 - (1 + 1/x) * (1-exp(-x)) /
// x) using a 5/5th order rational approximation. It is accurate to 1e-6 over
// [0, 1e6]. Developed by Colin Josey using Remez's algorithm, with original
// implementation in OpenMOC at:
// https://github.com/mit-crpg/OpenMOC/blob/develop/src/exponentials.h
float exponentialG2(float tau)
{

  // Coefficients for numerator in rational approximation
  constexpr float g1n = -0.08335775885589858f;
  constexpr float g2n = -0.003603942303847604f;
  constexpr float g3n = 0.0037673183263550827f;
  constexpr float g4n = 0.00001124183494990467f;
  constexpr float g5n = 0.00016837426505799449f;

  // Coefficients for denominator in rational approximation
  constexpr float g1d = 0.7454048371823628f;
  constexpr float g2d = 0.23794300531408347f;
  constexpr float g3d = 0.05367250964303789f;
  constexpr float g4d = 0.006125197988351906f;
  constexpr float g5d = 0.0010102514456857377f;

  float x = tau;

  float num = g5n;
  num = num * x + g4n;
  num = num * x + g3n;
  num = num * x + g2n;
  num = num * x + g1n;
  num = num * x;

  float den = g5d;
  den = den * x + g4d;
  den = den * x + g3d;
  den = den * x + g2d;
  den = den * x + g1d;
  den = den * x + 1.0f;

  return num / den;
}

//==============================================================================
// RandomRay implementation
//==============================================================================

// Static Variable Declarations
double RandomRay::distance_inactive_;
double RandomRay::distance_active_;
unique_ptr<Source> RandomRay::ray_source_;
RandomRaySourceShape RandomRay::source_shape_ {RandomRaySourceShape::FLAT};

RandomRay::RandomRay()
  : angular_flux_(data::mg.num_energy_groups_),
    delta_psi_(data::mg.num_energy_groups_),
    negroups_(data::mg.num_energy_groups_)
{
  if (source_shape_ == RandomRaySourceShape::LINEAR ||
      source_shape_ == RandomRaySourceShape::LINEAR_XY) {
    delta_moments_.resize(negroups_);
  }
}

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

void RandomRay::attenuate_flux(double distance, bool is_active)
{
  switch (source_shape_) {
  case RandomRaySourceShape::FLAT:
    attenuate_flux_flat_source(distance, is_active);
    break;
  case RandomRaySourceShape::LINEAR:
  case RandomRaySourceShape::LINEAR_XY:
    attenuate_flux_linear_source(distance, is_active);
    break;
  default:
    fatal_error("Unknown source shape for random ray transport.");
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
void RandomRay::attenuate_flux_flat_source(double distance, bool is_active)
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

void RandomRay::attenuate_flux_linear_source(double distance, bool is_active)
{
  // Cast domain to LinearSourceDomain
  LinearSourceDomain* domain = dynamic_cast<LinearSourceDomain*>(domain_);
  if (!domain) {
    fatal_error("RandomRay::attenuate_flux_linear_source() called with "
                "non-LinearSourceDomain domain.");
  }

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

  Position& centroid = domain->centroid_[source_region];
  Position midpoint = r() + u() * (distance / 2.0);

  // Determine the local position of the midpoint and the ray origin
  // relative to the source region's centroid
  Position rm_local;
  Position r0_local;

  // In the first few iterations of the simulation, the source region
  // may not yet have had any ray crossings, in which case there will
  // be no estimate of its centroid. We detect this by checking if it has
  // any accumulated volume. If its volume is zero, just use the midpoint
  // of the ray as the region's centroid.
  if (domain->volume_t_[source_region]) {
    rm_local = midpoint - centroid;
    r0_local = r() - centroid;
  } else {
    rm_local = {0.0, 0.0, 0.0};
    r0_local = -u() * 0.5 * distance;
  }
  double distance_2 = distance * distance;

  // Linear Source MOC incoming flux attenuation + source
  // contribution/attenuation equation
  for (int g = 0; g < negroups_; g++) {

    // Compute tau, the optical thickness of the ray segment
    float sigma_t = data::mg.macro_xs_[material].get_xs(
      MgxsType::TOTAL, g, NULL, NULL, NULL, t, a);
    float tau = sigma_t * distance;

    // If tau is very small, set it to zero to avoid numerical issues.
    // The following computations will still work with tau = 0.
    if (tau < 1.0e-8f) {
      tau = 0.0f;
    }

    // Compute linear source terms, spatial and directional (dir),
    // calculated from the source gradients dot product with local centroid
    // and direction, respectively.
    float spatial_source =
      domain_->source_[source_element + g] +
      rm_local.dot(domain->source_gradients_[source_element + g]);
    float dir_source = u().dot(domain->source_gradients_[source_element + g]);

    float gn = exponentialG(tau);
    float f1 = 1.0f - tau * gn;
    float f2 = (2.0f * gn - f1) * distance_2;
    float new_delta_psi = (angular_flux_[g] - spatial_source) * f1 * distance -
                          0.5 * dir_source * f2;

    float h1 = f1 - gn;
    float g1 = 0.5f - h1;
    float g2 = exponentialG2(tau);
    g1 = g1 * spatial_source;
    g2 = g2 * dir_source * distance * 0.5f;
    h1 = h1 * angular_flux_[g];
    h1 = (g1 + g2 + h1) * distance_2;
    spatial_source = spatial_source * distance + new_delta_psi;

    // Store contributions for this group into arrays, so that they can
    // be accumulated into the source region's estimates inside of the locked
    // region.
    delta_psi_[g] = new_delta_psi;
    delta_moments_[g] = r0_local * spatial_source + u() * h1;

    // Update the angular flux for this group
    angular_flux_[g] -= new_delta_psi * sigma_t;

    // If 2D mode is enabled, the z-component of the flux moments is forced
    // to zero
    if (source_shape_ == RandomRaySourceShape::LINEAR_XY) {
      delta_moments_[g].z = 0.0;
    }
  }

  // If ray is in the active phase (not in dead zone), make contributions to
  // source region bookkeeping
  if (is_active) {
    // Compute an estimate of the spatial moments matrix for the source
    // region based on parameters from this ray's crossing
    MomentMatrix moment_matrix_estimate;
    moment_matrix_estimate.compute_spatial_moments_matrix(
      rm_local, u(), distance);

    // Aquire lock for source region
    domain_->lock_[source_region].lock();

    // Accumulate deltas into the new estimate of source region flux for this
    // iteration
    for (int g = 0; g < negroups_; g++) {
      domain_->scalar_flux_new_[source_element + g] += delta_psi_[g];
      domain->flux_moments_new_[source_element + g] += delta_moments_[g];
    }

    // Accumulate the volume (ray segment distance), centroid, and spatial
    // momement estimates into the running totals for the iteration for this
    // source region. The centroid and spatial momements estimates are scaled by
    // the ray segment length as part of length averaging of the estimates.
    domain_->volume_[source_region] += distance;
    domain->centroid_iteration_[source_region] += midpoint * distance;
    moment_matrix_estimate *= distance;
    domain->mom_matrix_[source_region] += moment_matrix_estimate;

    // If the source region hasn't been hit yet this iteration,
    // indicate that it now has
    if (domain_->was_hit_[source_region] == 0) {
      domain_->was_hit_[source_region] = 1;
    }

    // Tally valid position inside the source region (e.g., midpoint of
    // the ray) if not done already
    if (!domain_->position_recorded_[source_region]) {
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
