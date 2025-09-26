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

#include "openmc/distribution_spatial.h"
#include "openmc/random_dist.h"
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

// Implementation of the Fisher-Yates shuffle algorithm.
// Algorithm adapted from:
//    https://en.cppreference.com/w/cpp/algorithm/random_shuffle#Version_3
void fisher_yates_shuffle(vector<int64_t>& arr, uint64_t* seed)
{
  // Loop over the array from the last element down to the second
  for (int i = arr.size() - 1; i > 0; --i) {
    // Generate a random index in the range [0, i]
    int j = uniform_int_distribution(0, i, seed);
    std::swap(arr[i], arr[j]);
  }
}

// Function to generate randomized Halton sequence samples
//
// Algorithm adapted from:
//      A. B. Owen. A randomized halton algorithm in r. Arxiv, 6 2017.
//      URL https://arxiv.org/abs/1706.02808
vector<double> rhalton(int dim, uint64_t* seed, int64_t skip = 0)
{
  if (dim > 10) {
    fatal_error("Halton sampling dimension too large");
  }
  int64_t b, res, dig;
  double b2r, ans;
  const std::array<int64_t, 10> primes = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29};
  vector<double> halton(dim, 0.0);

  vector<int64_t> perm;
  for (int D = 0; D < dim; ++D) {
    b = primes[D];
    perm.resize(b);
    b2r = 1.0 / b;
    res = skip;
    ans = 0.0;

    while ((1.0 - b2r) < 1.0) {
      std::iota(perm.begin(), perm.end(), 0);
      fisher_yates_shuffle(perm, seed);
      dig = res % b;
      ans += perm[dig] * b2r;
      res = (res - dig) / b;
      b2r /= b;
    }

    halton[D] = ans;
  }

  return halton;
}

//==============================================================================
// RandomRay implementation
//==============================================================================

// Static Variable Declarations
double RandomRay::distance_inactive_;
double RandomRay::distance_active_;
unique_ptr<Source> RandomRay::ray_source_;
RandomRaySourceShape RandomRay::source_shape_ {RandomRaySourceShape::FLAT};
RandomRaySampleMethod RandomRay::sample_method_ {RandomRaySampleMethod::PRNG};

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
    // If ray has too many events, display warning and kill it
    if (n_event() >= settings::max_particle_events) {
      warning("Ray " + std::to_string(id()) +
              " underwent maximum number of events, terminating ray.");
      wgt() = 0.0;
    }
  }

  return n_event();
}

// Transports ray across a single source region
void RandomRay::event_advance_ray()
{
  // Find the distance to the nearest boundary
  boundary() = distance_to_boundary(*this);
  double distance = boundary().distance();

  if (distance < 0.0) {
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

      attenuate_flux(distance_alive, true, distance_dead);
      distance_travelled_ = distance_alive;
    } else {
      distance_travelled_ += distance;
      attenuate_flux(distance, false);
    }
  }

  // Advance particle
  for (int j = 0; j < n_coord(); ++j) {
    coord(j).r() += distance * coord(j).u();
  }
}

void RandomRay::attenuate_flux(double distance, bool is_active, double offset)
{
  // Lookup base source region index
  int64_t sr = domain_->lookup_base_source_region_idx(*this);

  // Perform ray tracing across mesh
  // Determine the mesh index for the base source region, if any
  int mesh_idx = domain_->lookup_mesh_idx(sr);

  if (mesh_idx == C_NONE) {
    // If there's no mesh being applied to this cell, then
    // we just attenuate the flux as normal, and set
    // the mesh bin to 0
    attenuate_flux_inner(distance, is_active, sr, 0, r());
  } else {
    // If there is a mesh being applied to this cell, then
    // we loop over all the bin crossings and attenuate
    // separately.
    Mesh* mesh = model::meshes[mesh_idx].get();

    // We adjust the start and end positions of the ray slightly
    // to accomodate for floating point precision issues that tend
    // to occur at mesh boundaries that overlap with geometry lattice
    // boundaries.
    Position start = r() + (offset + TINY_BIT) * u();
    Position end = start + (distance - 2.0 * TINY_BIT) * u();
    double reduced_distance = (end - start).norm();

    // Ray trace through the mesh and record bins and lengths
    mesh_bins_.resize(0);
    mesh_fractional_lengths_.resize(0);
    mesh->bins_crossed(start, end, u(), mesh_bins_, mesh_fractional_lengths_);

    // Loop over all mesh bins and attenuate flux
    for (int b = 0; b < mesh_bins_.size(); b++) {
      double physical_length = reduced_distance * mesh_fractional_lengths_[b];
      attenuate_flux_inner(
        physical_length, is_active, sr, mesh_bins_[b], start);
      start += physical_length * u();
    }
  }
}

void RandomRay::attenuate_flux_inner(
  double distance, bool is_active, int64_t sr, int mesh_bin, Position r)
{
  SourceRegionKey sr_key {sr, mesh_bin};
  SourceRegionHandle srh;
  srh = domain_->get_subdivided_source_region_handle(sr_key, r, u());
  if (srh.is_numerical_fp_artifact_) {
    return;
  }

  switch (source_shape_) {
  case RandomRaySourceShape::FLAT:
    if (srh.material() == MATERIAL_VOID) {
      attenuate_flux_flat_source_void(srh, distance, is_active, r);
    } else {
      attenuate_flux_flat_source(srh, distance, is_active, r);
    }
    break;
  case RandomRaySourceShape::LINEAR:
  case RandomRaySourceShape::LINEAR_XY:
    if (srh.material() == MATERIAL_VOID) {
      attenuate_flux_linear_source_void(srh, distance, is_active, r);
    } else {
      attenuate_flux_linear_source(srh, distance, is_active, r);
    }
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
void RandomRay::attenuate_flux_flat_source(
  SourceRegionHandle& srh, double distance, bool is_active, Position r)
{
  // The number of geometric intersections is counted for reporting purposes
  n_event()++;

  // Get material
  int material = srh.material();

  // MOC incoming flux attenuation + source contribution/attenuation equation
  for (int g = 0; g < negroups_; g++) {
    float sigma_t = domain_->sigma_t_[material * negroups_ + g];
    float tau = sigma_t * distance;
    float exponential = cjosey_exponential(tau); // exponential = 1 - exp(-tau)
    float new_delta_psi = (angular_flux_[g] - srh.source(g)) * exponential;
    delta_psi_[g] = new_delta_psi;
    angular_flux_[g] -= new_delta_psi;
  }

  // If ray is in the active phase (not in dead zone), make contributions to
  // source region bookkeeping

  // Aquire lock for source region
  srh.lock();

  if (is_active) {
    // Accumulate delta psi into new estimate of source region flux for
    // this iteration
    for (int g = 0; g < negroups_; g++) {
      srh.scalar_flux_new(g) += delta_psi_[g];
    }

    // Accomulate volume (ray distance) into this iteration's estimate
    // of the source region's volume
    srh.volume() += distance;

    srh.n_hits() += 1;
  }

  // Tally valid position inside the source region (e.g., midpoint of
  // the ray) if not done already
  if (!srh.position_recorded()) {
    Position midpoint = r + u() * (distance / 2.0);
    srh.position() = midpoint;
    srh.position_recorded() = 1;
  }

  // Release lock
  srh.unlock();
}

// Alternative flux attenuation function for true void regions.
void RandomRay::attenuate_flux_flat_source_void(
  SourceRegionHandle& srh, double distance, bool is_active, Position r)
{
  // The number of geometric intersections is counted for reporting purposes
  n_event()++;

  int material = srh.material();

  // If ray is in the active phase (not in dead zone), make contributions to
  // source region bookkeeping
  if (is_active) {

    // Aquire lock for source region
    srh.lock();

    // Accumulate delta psi into new estimate of source region flux for
    // this iteration
    for (int g = 0; g < negroups_; g++) {
      srh.scalar_flux_new(g) += angular_flux_[g] * distance;
    }

    // Accomulate volume (ray distance) into this iteration's estimate
    // of the source region's volume
    srh.volume() += distance;
    srh.volume_sq() += distance * distance;
    srh.n_hits() += 1;

    // Tally valid position inside the source region (e.g., midpoint of
    // the ray) if not done already
    if (!srh.position_recorded()) {
      Position midpoint = r + u() * (distance / 2.0);
      srh.position() = midpoint;
      srh.position_recorded() = 1;
    }

    // Release lock
    srh.unlock();
  }

  // Add source to incoming angular flux, assuming void region
  if (settings::run_mode == RunMode::FIXED_SOURCE) {
    for (int g = 0; g < negroups_; g++) {
      angular_flux_[g] += srh.external_source(g) * distance;
    }
  }
}

void RandomRay::attenuate_flux_linear_source(
  SourceRegionHandle& srh, double distance, bool is_active, Position r)
{
  // The number of geometric intersections is counted for reporting purposes
  n_event()++;

  int material = srh.material();

  Position& centroid = srh.centroid();
  Position midpoint = r + u() * (distance / 2.0);

  // Determine the local position of the midpoint and the ray origin
  // relative to the source region's centroid
  Position rm_local;
  Position r0_local;

  // In the first few iterations of the simulation, the source region
  // may not yet have had any ray crossings, in which case there will
  // be no estimate of its centroid. We detect this by checking if it has
  // any accumulated volume. If its volume is zero, just use the midpoint
  // of the ray as the region's centroid.
  if (srh.volume_t()) {
    rm_local = midpoint - centroid;
    r0_local = r - centroid;
  } else {
    rm_local = {0.0, 0.0, 0.0};
    r0_local = -u() * 0.5 * distance;
  }
  double distance_2 = distance * distance;

  // Linear Source MOC incoming flux attenuation + source
  // contribution/attenuation equation
  for (int g = 0; g < negroups_; g++) {

    // Compute tau, the optical thickness of the ray segment
    float sigma_t = domain_->sigma_t_[material * negroups_ + g];
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
      srh.source(g) + rm_local.dot(srh.source_gradients(g));
    float dir_source = u().dot(srh.source_gradients(g));

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

  // Compute an estimate of the spatial moments matrix for the source
  // region based on parameters from this ray's crossing
  MomentMatrix moment_matrix_estimate;
  moment_matrix_estimate.compute_spatial_moments_matrix(
    rm_local, u(), distance);

  // Aquire lock for source region
  srh.lock();

  // If ray is in the active phase (not in dead zone), make contributions to
  // source region bookkeeping

  if (is_active) {
    // Accumulate deltas into the new estimate of source region flux for this
    // iteration
    for (int g = 0; g < negroups_; g++) {
      srh.scalar_flux_new(g) += delta_psi_[g];
      srh.flux_moments_new(g) += delta_moments_[g];
    }

    // Accumulate the volume (ray segment distance), centroid, and spatial
    // momement estimates into the running totals for the iteration for this
    // source region. The centroid and spatial momements estimates are scaled
    // by the ray segment length as part of length averaging of the estimates.
    srh.volume() += distance;
    srh.centroid_iteration() += midpoint * distance;
    moment_matrix_estimate *= distance;
    srh.mom_matrix() += moment_matrix_estimate;

    srh.n_hits() += 1;
  }

  // Tally valid position inside the source region (e.g., midpoint of
  // the ray) if not done already
  if (!srh.position_recorded()) {
    srh.position() = midpoint;
    srh.position_recorded() = 1;
  }

  // Release lock
  srh.unlock();
}

// If traveling through a void region, the source term is either zero
// or an external source. As all external sources are currently assumed
// to be flat, we don't really need this function and could instead just call
// the "attenuate_flux_flat_source_void" function and get the same numerical and
// tally results. However, computation of the flux moments in void regions is
// nonetheless useful as this information is still used by the plotter when
// estimating the flux at specific pixel coordinates. Thus, plots will look
// nicer/more accurate if we record flux moments, so this function is useful.
void RandomRay::attenuate_flux_linear_source_void(
  SourceRegionHandle& srh, double distance, bool is_active, Position r)
{
  // The number of geometric intersections is counted for reporting purposes
  n_event()++;

  Position& centroid = srh.centroid();
  Position midpoint = r + u() * (distance / 2.0);

  // Determine the local position of the midpoint and the ray origin
  // relative to the source region's centroid
  Position rm_local;
  Position r0_local;

  // In the first few iterations of the simulation, the source region
  // may not yet have had any ray crossings, in which case there will
  // be no estimate of its centroid. We detect this by checking if it has
  // any accumulated volume. If its volume is zero, just use the midpoint
  // of the ray as the region's centroid.
  if (srh.volume_t()) {
    rm_local = midpoint - centroid;
    r0_local = r - centroid;
  } else {
    rm_local = {0.0, 0.0, 0.0};
    r0_local = -u() * 0.5 * distance;
  }
  double distance_2 = distance * distance;

  // Compared to linear flux attenuation through solid regions,
  // transport through a void region is greatly simplified. Here we
  // compute the updated flux moments.
  for (int g = 0; g < negroups_; g++) {
    float spatial_source = 0.f;
    if (settings::run_mode == RunMode::FIXED_SOURCE) {
      spatial_source = srh.external_source(g);
    }
    float new_delta_psi = (angular_flux_[g] - spatial_source) * distance;
    float h1 = 0.5f;
    h1 = h1 * angular_flux_[g];
    h1 = h1 * distance_2;
    spatial_source = spatial_source * distance + new_delta_psi;

    // Store contributions for this group into arrays, so that they can
    // be accumulated into the source region's estimates inside of the locked
    // region.
    delta_moments_[g] = r0_local * spatial_source + u() * h1;

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
    srh.lock();

    // Accumulate delta psi into new estimate of source region flux for
    // this iteration, and update flux momements
    for (int g = 0; g < negroups_; g++) {
      srh.scalar_flux_new(g) += angular_flux_[g] * distance;
      srh.flux_moments_new(g) += delta_moments_[g];
    }

    // Accumulate the volume (ray segment distance), centroid, and spatial
    // momement estimates into the running totals for the iteration for this
    // source region. The centroid and spatial momements estimates are scaled by
    // the ray segment length as part of length averaging of the estimates.
    srh.volume() += distance;
    srh.volume_sq() += distance_2;
    srh.centroid_iteration() += midpoint * distance;
    moment_matrix_estimate *= distance;
    srh.mom_matrix() += moment_matrix_estimate;

    // Tally valid position inside the source region (e.g., midpoint of
    // the ray) if not done already
    if (!srh.position_recorded()) {
      srh.position() = midpoint;
      srh.position_recorded() = 1;
    }

    srh.n_hits() += 1;

    // Release lock
    srh.unlock();
  }

  // Add source to incoming angular flux, assuming void region
  if (settings::run_mode == RunMode::FIXED_SOURCE) {
    for (int g = 0; g < negroups_; g++) {
      angular_flux_[g] += srh.external_source(g) * distance;
    }
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
  id() = ray_id;

  // generate source site using sample method
  SourceSite site;
  switch (sample_method_) {
  case RandomRaySampleMethod::PRNG:
    site = sample_prng();
    break;
  case RandomRaySampleMethod::HALTON:
    site = sample_halton();
    break;
  default:
    fatal_error("Unknown sample method for random ray transport.");
  }

  site.E = 0.0;
  this->from_source(&site);

  // Locate ray
  if (lowest_coord().cell() == C_NONE) {
    if (!exhaustive_find_cell(*this)) {
      this->mark_as_lost(
        "Could not find the cell containing particle " + std::to_string(id()));
    }

    // Set birth cell attribute
    if (cell_born() == C_NONE)
      cell_born() = lowest_coord().cell();
  }

  SourceRegionKey sr_key = domain_->lookup_source_region_key(*this);
  SourceRegionHandle srh =
    domain_->get_subdivided_source_region_handle(sr_key, r(), u());

  // Initialize ray's starting angular flux to starting location's isotropic
  // source
  if (!srh.is_numerical_fp_artifact_) {
    for (int g = 0; g < negroups_; g++) {
      angular_flux_[g] = srh.source(g);
    }
  }
}

SourceSite RandomRay::sample_prng()
{
  // set random number seed
  int64_t particle_seed =
    (simulation::current_batch - 1) * settings::n_particles + id();
  init_particle_seeds(particle_seed, seeds());
  stream() = STREAM_TRACKING;

  // Sample from ray source distribution
  SourceSite site {ray_source_->sample(current_seed())};

  return site;
}

SourceSite RandomRay::sample_halton()
{
  SourceSite site;

  // Set random number seed
  int64_t batch_seed = (simulation::current_batch - 1) * settings::n_particles;
  int64_t skip = id();
  init_particle_seeds(batch_seed, seeds());
  stream() = STREAM_TRACKING;

  // Calculate next samples in LDS across 5 dimensions
  vector<double> samples = rhalton(5, current_seed(), skip = skip);

  // Get spatial box of ray_source_
  SpatialBox* sb = dynamic_cast<SpatialBox*>(
    dynamic_cast<IndependentSource*>(RandomRay::ray_source_.get())->space());

  // Sample spatial distribution
  Position xi {samples[0], samples[1], samples[2]};
  // make a small shift in position to avoid geometry floating point issues
  Position shift {FP_COINCIDENT, FP_COINCIDENT, FP_COINCIDENT};
  site.r = (sb->lower_left() + shift) +
           xi * ((sb->upper_right() - shift) - (sb->lower_left() + shift));

  // Sample Polar cosine and azimuthal angles
  double mu = 2.0 * samples[3] - 1.0;
  double azi = 2.0 * PI * samples[4];
  // Convert to Cartesian coordinates
  double c = std::sqrt(1.0 - mu * mu);
  site.u.x = mu;
  site.u.y = std::cos(azi) * c;
  site.u.z = std::sin(azi) * c;

  return site;
}

} // namespace openmc
