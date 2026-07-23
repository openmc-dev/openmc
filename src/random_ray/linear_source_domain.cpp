#include "openmc/random_ray/linear_source_domain.h"

#include "openmc/cell.h"
#include "openmc/geometry.h"
#include "openmc/material.h"
#include "openmc/message_passing.h"
#include "openmc/mgxs_interface.h"
#include "openmc/output.h"
#include "openmc/plot.h"
#include "openmc/random_ray/random_ray.h"
#include "openmc/simulation.h"
#include "openmc/tallies/filter.h"
#include "openmc/tallies/tally.h"
#include "openmc/tallies/tally_scoring.h"
#include "openmc/timer.h"

namespace openmc {

//==============================================================================
// LinearSourceDomain implementation
//==============================================================================

void LinearSourceDomain::batch_reset()
{
  FlatSourceDomain::batch_reset();
#pragma omp parallel for
  for (int64_t sr = 0; sr < n_source_regions(); sr++) {
    source_regions_.centroid_iteration(sr) = {0.0, 0.0, 0.0};
    source_regions_.mom_matrix(sr) = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  }
#pragma omp parallel for
  for (int64_t se = 0; se < n_source_elements(); se++) {
    source_regions_.flux_moments_new(se) = {0.0, 0.0, 0.0};
  }
}

void LinearSourceDomain::update_single_neutron_source(SourceRegionHandle& srh)
{
  // Reset all source regions to zero (important for void regions)
  for (int g = 0; g < negroups_; g++) {
    srh.source(g) = 0.0;
  }

  // Add scattering + fission source
  int material = srh.material();
  int temp = srh.temperature_idx();
  double density_mult = srh.density_mult();
  if (material != MATERIAL_VOID) {
    double inverse_k_eff = 1.0 / k_eff_;
    MomentMatrix invM = srh.mom_matrix().inverse();

    for (int g_out = 0; g_out < negroups_; g_out++) {
      double sigma_t =
        sigma_t_[(material * ntemperature_ + temp) * negroups_ + g_out] *
        density_mult;

      double scatter_flat = 0.0f;
      double fission_flat = 0.0f;
      MomentArray scatter_linear = {0.0, 0.0, 0.0};
      MomentArray fission_linear = {0.0, 0.0, 0.0};

      for (int g_in = 0; g_in < negroups_; g_in++) {
        // Handles for the flat and linear components of the flux
        double flux_flat = srh.scalar_flux_old(g_in);
        MomentArray flux_linear = srh.flux_moments_old(g_in);

        // Handles for cross sections
        double sigma_s =
          sigma_s_[((material * ntemperature_ + temp) * negroups_ + g_out) *
                     negroups_ +
                   g_in] *
          density_mult;
        double nu_sigma_f =
          nu_sigma_f_[(material * ntemperature_ + temp) * negroups_ + g_in] *
          density_mult;
        double chi =
          chi_[(material * ntemperature_ + temp) * negroups_ + g_out];

        // Compute source terms for flat and linear components of the flux
        scatter_flat += sigma_s * flux_flat;
        scatter_linear += sigma_s * flux_linear;
        if (settings::create_fission_neutrons) {
          fission_flat += nu_sigma_f * flux_flat * chi;
          fission_linear += nu_sigma_f * flux_linear * chi;
        }
      }

      // Compute the flat source term
      srh.source(g_out) =
        (scatter_flat + fission_flat * inverse_k_eff) / sigma_t;

      // Compute the linear source terms. In the first 10 iterations when the
      // centroids and spatial moments are not well known, we will leave the
      // source gradients as zero so as to avoid causing any numerical
      // instability. If a negative source is encountered, this region must be
      // very small/noisy or have poorly developed spatial moments, so we zero
      // the source gradients (effectively making this a flat source region
      // temporarily), so as to improve stability.
      if (simulation::current_batch > 10 && srh.source(g_out) >= 0.0) {
        srh.source_gradients(g_out) =
          invM * ((scatter_linear + fission_linear * inverse_k_eff) / sigma_t);
      } else {
        srh.source_gradients(g_out) = {0.0, 0.0, 0.0};
      }
    }
  }

  // Add external source if in fixed source mode
  if (settings::run_mode == RunMode::FIXED_SOURCE) {
    for (int g = 0; g < negroups_; g++) {
      srh.source(g) += srh.external_source(g);
    }
  }
}

void LinearSourceDomain::normalize_scalar_flux_and_volumes(
  double total_active_distance_per_iteration)
{
  double normalization_factor = 1.0 / total_active_distance_per_iteration;
  double volume_normalization_factor =
    1.0 / (total_active_distance_per_iteration * simulation::current_batch);

// Normalize flux to total distance travelled by all rays this iteration
#pragma omp parallel for
  for (int64_t se = 0; se < n_source_elements(); se++) {
    source_regions_.scalar_flux_new(se) *= normalization_factor;
    source_regions_.flux_moments_new(se) *= normalization_factor;
  }

// Accumulate cell-wise ray length tallies collected this iteration, then
// update the simulation-averaged cell-wise volume estimates
#pragma omp parallel for
  for (int64_t sr = 0; sr < n_source_regions(); sr++) {
    // Offset of this batch's track centroid from the accumulated centroid
    // the transport sweep evaluated the linear source against (the value
    // held here before this batch is folded in). The flux update uses it to
    // add back the source as actually sampled by this batch's tracks (see
    // set_flux_to_flux_plus_source). Regions with no prior accumulated
    // volume used per-segment midpoints in place of a centroid, for which
    // the batch offset is identically zero.
    if (source_regions_.volume_t(sr) > 0.0 &&
        source_regions_.volume(sr) > 0.0) {
      source_regions_.centroid_offset(sr) =
        source_regions_.centroid_iteration(sr) *
          (1.0 / source_regions_.volume(sr)) -
        source_regions_.centroid(sr);
    } else {
      source_regions_.centroid_offset(sr) = {0.0, 0.0, 0.0};
    }
    source_regions_.centroid_t(sr) += source_regions_.centroid_iteration(sr);
    source_regions_.mom_matrix_t(sr) += source_regions_.mom_matrix(sr);
    source_regions_.volume_t(sr) += source_regions_.volume(sr);
    source_regions_.volume_sq_t(sr) += source_regions_.volume_sq(sr);
    source_regions_.volume_naive(sr) =
      source_regions_.volume(sr) * normalization_factor;
    source_regions_.volume(sr) =
      source_regions_.volume_t(sr) * volume_normalization_factor;
    source_regions_.volume_sq(sr) =
      source_regions_.volume_sq_t(sr) / source_regions_.volume_t(sr);
    if (source_regions_.volume_t(sr) > 0.0) {
      double inv_volume = 1.0 / source_regions_.volume_t(sr);
      source_regions_.centroid(sr) = source_regions_.centroid_t(sr);
      source_regions_.centroid(sr) *= inv_volume;
      source_regions_.mom_matrix(sr) = source_regions_.mom_matrix_t(sr);
      source_regions_.mom_matrix(sr) *= inv_volume;
    }
  }
}

void LinearSourceDomain::set_flux_to_flux_plus_source(
  int64_t sr, double volume, int g)
{
  int material = source_regions_.material(sr);
  if (material == MATERIAL_VOID) {
    FlatSourceDomain::set_flux_to_flux_plus_source(sr, volume, g);
  } else {
    source_regions_.scalar_flux_new(sr, g) /= volume;
    // The source term added back includes the gradient part of the source
    // that this batch's rays actually integrated. The sweep evaluates the
    // linear source against the accumulated centroid, but the batch's tracks
    // average it at their own track-length-weighted centroid: summing the
    // exact chord integrals of a linear source over the batch's tracks gives
    // a mean emission of source + gradient . centroid_offset, entering the
    // transport sum with the batch's tracked volume -- hence the
    // volume_naive / volume weight, which also prevents the correlation
    // between the batch volume and batch centroid (both sampled by the same
    // rays) from biasing the update when the two volumes differ. Adding back
    // the flat source alone instead leaves the difference in the flux
    // estimate as noise proportional to the gradient times the per-batch
    // centroid fluctuation. In near-cancellation regions (scattering ratio
    // near one), where the reduced source approximately equals the flux,
    // that noise can exceed the flux itself at modest per-batch hit counts
    // and ignite self-sustaining negativity through the scattering feedback.
    // When the region is updated with its own batch (naive) volume the
    // weight is one and the update becomes an exact per-batch track
    // identity: the flux estimate is the track-length average of the angular
    // flux, and so inherits its sign -- restoring to linear sources the
    // property that a naive-volume update adds no noise of its own beyond
    // the angular fluxes it averages.
    source_regions_.scalar_flux_new(sr, g) +=
      source_regions_.source(sr, g) +
      (source_regions_.volume_naive(sr) / volume) *
        source_regions_.source_gradients(sr, g)
          .dot(source_regions_.centroid_offset(sr));
  }
  // If a source region is small, then the moments are likely noisy, so we zero
  // them. This is reasonable, given that small regions can get by with a flat
  // source approximation anyhow.
  if (source_regions_.is_small(sr)) {
    source_regions_.flux_moments_new(sr, g) = {0.0, 0.0, 0.0};
  } else {
    source_regions_.flux_moments_new(sr, g) *= (1.0 / volume);
  }
}

void LinearSourceDomain::set_flux_to_old_flux(int64_t sr, int g)
{
  source_regions_.scalar_flux_new(sr, g) =
    source_regions_.scalar_flux_old(sr, g);
  source_regions_.flux_moments_new(sr, g) =
    source_regions_.flux_moments_old(sr, g);
}

void LinearSourceDomain::accumulate_iteration_flux()
{
  // Accumulate scalar flux
  FlatSourceDomain::accumulate_iteration_flux();

  // Accumulate scalar flux moments
#pragma omp parallel for
  for (int64_t se = 0; se < n_source_elements(); se++) {
    source_regions_.flux_moments_t(se) += source_regions_.flux_moments_new(se);
  }
}

double LinearSourceDomain::evaluate_flux_at_point(
  Position r, int64_t sr, int g) const
{
  double phi_flat = FlatSourceDomain::evaluate_flux_at_point(r, sr, g);

  Position local_r = r - source_regions_.centroid(sr);
  MomentArray phi_linear = source_regions_.flux_moments_t(sr, g);
  phi_linear *= 1.0 / (settings::n_batches - settings::n_inactive);

  MomentMatrix invM = source_regions_.mom_matrix(sr).inverse();
  MomentArray phi_solved = invM * phi_linear;

  return phi_flat + phi_solved.dot(local_r);
}

} // namespace openmc
