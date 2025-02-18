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
  for (int64_t sr = 0; sr < n_source_regions_; sr++) {
    source_regions_.centroid_iteration(sr) = {0.0, 0.0, 0.0};
    source_regions_.mom_matrix(sr) = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  }
#pragma omp parallel for
  for (int64_t se = 0; se < n_source_elements_; se++) {
    source_regions_.flux_moments_new(se) = {0.0, 0.0, 0.0};
  }
}

void LinearSourceDomain::update_neutron_source(double k_eff)
{
  simulation::time_update_src.start();

  double inverse_k_eff = 1.0 / k_eff;

#pragma omp parallel for
  for (int64_t sr = 0; sr < n_source_regions_; sr++) {
    int material = source_regions_.material(sr);
    if (material == MATERIAL_VOID) {
      continue;
    }
    MomentMatrix invM = source_regions_.mom_matrix(sr).inverse();

    for (int g_out = 0; g_out < negroups_; g_out++) {
      double sigma_t = sigma_t_[material * negroups_ + g_out];

      double scatter_flat = 0.0f;
      double fission_flat = 0.0f;
      MomentArray scatter_linear = {0.0, 0.0, 0.0};
      MomentArray fission_linear = {0.0, 0.0, 0.0};

      for (int g_in = 0; g_in < negroups_; g_in++) {
        // Handles for the flat and linear components of the flux
        double flux_flat = source_regions_.scalar_flux_old(sr, g_in);
        MomentArray flux_linear = source_regions_.flux_moments_old(sr, g_in);

        // Handles for cross sections
        double sigma_s =
          sigma_s_[material * negroups_ * negroups_ + g_out * negroups_ + g_in];
        double nu_sigma_f = nu_sigma_f_[material * negroups_ + g_in];
        double chi = chi_[material * negroups_ + g_out];

        // Compute source terms for flat and linear components of the flux
        scatter_flat += sigma_s * flux_flat;
        fission_flat += nu_sigma_f * flux_flat * chi;
        scatter_linear += sigma_s * flux_linear;
        fission_linear += nu_sigma_f * flux_linear * chi;
      }

      // Compute the flat source term
      source_regions_.source(sr, g_out) =
        (scatter_flat + fission_flat * inverse_k_eff) / sigma_t;

      // Compute the linear source terms
      // In the first 10 iterations when the centroids and spatial moments
      // are not well known, we will leave the source gradients as zero
      // so as to avoid causing any numerical instability.
      if (simulation::current_batch > 10) {
        source_regions_.source_gradients(sr, g_out) =
          invM * ((scatter_linear + fission_linear * inverse_k_eff) / sigma_t);
      }
    }
  }

  if (settings::run_mode == RunMode::FIXED_SOURCE) {
// Add external source to flat source term if in fixed source mode
#pragma omp parallel for
    for (int64_t se = 0; se < n_source_elements_; se++) {
      source_regions_.source(se) += source_regions_.external_source(se);
    }
  }

  simulation::time_update_src.stop();
}

void LinearSourceDomain::normalize_scalar_flux_and_volumes(
  double total_active_distance_per_iteration)
{
  double normalization_factor = 1.0 / total_active_distance_per_iteration;
  double volume_normalization_factor =
    1.0 / (total_active_distance_per_iteration * simulation::current_batch);

// Normalize flux to total distance travelled by all rays this iteration
#pragma omp parallel for
  for (int64_t se = 0; se < n_source_elements_; se++) {
    source_regions_.scalar_flux_new(se) *= normalization_factor;
    source_regions_.flux_moments_new(se) *= normalization_factor;
  }

// Accumulate cell-wise ray length tallies collected this iteration, then
// update the simulation-averaged cell-wise volume estimates
#pragma omp parallel for
  for (int64_t sr = 0; sr < n_source_regions_; sr++) {
    source_regions_.centroid_t(sr) += source_regions_.centroid_iteration(sr);
    source_regions_.mom_matrix_t(sr) += source_regions_.mom_matrix(sr);
    source_regions_.volume_t(sr) += source_regions_.volume(sr);
    source_regions_.volume_sq_t(sr) += source_regions_.volume_sq(sr);
    source_regions_.volume_naive(sr) =
      source_regions_.volume(sr) * normalization_factor;
    source_regions_.volume(sr) =
      source_regions_.volume_t(sr) * volume_normalization_factor;
    source_regions_.volume_sq(sr) =
      (source_regions_.volume_sq_t(sr) / source_regions_.volume_t(sr)) *
      volume_normalization_factor;
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
    source_regions_.scalar_flux_new(sr, g) += source_regions_.source(sr, g);
  }
  source_regions_.flux_moments_new(sr, g) *= (1.0 / volume);
}

void LinearSourceDomain::set_flux_to_old_flux(int64_t sr, int g)
{
  source_regions_.scalar_flux_new(g) = source_regions_.scalar_flux_old(g);
  source_regions_.flux_moments_new(g) = source_regions_.flux_moments_old(g);
}

void LinearSourceDomain::accumulate_iteration_flux()
{
  // Accumulate scalar flux
  FlatSourceDomain::accumulate_iteration_flux();

  // Accumulate scalar flux moments
#pragma omp parallel for
  for (int64_t se = 0; se < n_source_elements_; se++) {
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
