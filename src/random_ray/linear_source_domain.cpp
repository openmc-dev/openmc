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

#include <cstdio>

namespace openmc {

//==============================================================================
// LinearSourceDomain implementation
//==============================================================================

LinearSourceDomain::LinearSourceDomain() : FlatSourceDomain()
{
  flux_moments_old_.assign(n_source_elements_, {0.0, 0.0, 0.0});
  flux_moments_new_.assign(n_source_elements_, {0.0, 0.0, 0.0});
  flux_moments_t_.assign(n_source_elements_, {0.0, 0.0, 0.0});
  source_moments_.assign(n_source_elements_, {0.0, 0.0, 0.0});

  centroid_.assign(n_source_regions_, {nan(""), nan(""), nan("")});
  centroid_t_.assign(n_source_regions_, {0.0, 0.0, 0.0});
  mom_matrix_.assign(n_source_regions_, {0.0, 0.0, 0.0, 0.0, 0.0, 0.0});
  mom_matrix_t_.assign(n_source_regions_, {0.0, 0.0, 0.0, 0.0, 0.0, 0.0});
}

void LinearSourceDomain::batch_reset()
{
  FlatSourceDomain::batch_reset();
#pragma omp parallel for
  for (auto& m : flux_moments_new_) {
    m = {0.0, 0.0, 0.0};
  }
}

void LinearSourceDomain::update_neutron_source(double k_eff)
{
  simulation::time_update_src.start();

  double inverse_k_eff = 1.0 / k_eff;

  // Temperature and angle indices, if using multiple temperature
  // data sets and/or anisotropic data sets.
  // TODO: Currently assumes we are only using single temp/single
  // angle data.
  const int t = 0;
  const int a = 0;
// double det;
#pragma omp parallel for
  for (int sr = 0; sr < n_source_regions_; sr++) {

    int material = material_[sr];
    SymmetricMatrix invM = mom_matrix_[sr].inverse();

    for (int e_out = 0; e_out < negroups_; e_out++) {
      float sigma_t = data::mg.macro_xs_[material].get_xs(
        MgxsType::TOTAL, e_out, nullptr, nullptr, nullptr, t, a);

      float scatter_flat = 0.0f;
      float fission_flat = 0.0f;
      MomentArray scatter_linear = {0.0, 0.0, 0.0};
      MomentArray fission_linear = {0.0, 0.0, 0.0};

      for (int e_in = 0; e_in < negroups_; e_in++) {
        // Handles for the flat and linear components of the flux
        float flux_flat = scalar_flux_old_[sr * negroups_ + e_in];
        MomentArray flux_linear = flux_moments_old_[sr * negroups_ + e_in];

        // Handles for cross sections
        float sigma_s = data::mg.macro_xs_[material].get_xs(
          MgxsType::NU_SCATTER, e_in, &e_out, nullptr, nullptr, t, a);
        float nu_sigma_f = data::mg.macro_xs_[material].get_xs(
          MgxsType::NU_FISSION, e_in, nullptr, nullptr, nullptr, t, a);
        float chi = data::mg.macro_xs_[material].get_xs(
          MgxsType::CHI_PROMPT, e_in, &e_out, nullptr, nullptr, t, a);

        // Compute source terms for flat and linear components of the flux
        scatter_flat += sigma_s * flux_flat;
        fission_flat += nu_sigma_f * flux_flat * chi;
        scatter_linear += sigma_s * flux_linear;
        fission_linear += nu_sigma_f * flux_linear * chi;
      }

      // Compute the flat source term
      source_[sr * negroups_ + e_out] =
        (scatter_flat + fission_flat * inverse_k_eff) / sigma_t;

      // Compute the linear source terms
      if (simulation::current_batch > 2) {
        source_moments_[sr * negroups_ + e_out] =
          invM * ((scatter_linear + fission_linear * inverse_k_eff) / sigma_t);
      }
    }
  }

  if (settings::run_mode == RunMode::FIXED_SOURCE) {
// Add external source to flat source term if in fixed source mode
#pragma omp parallel for
    for (int se = 0; se < n_source_elements_; se++) {
      source_[se] += external_source_[se];
    }
  }

  simulation::time_update_src.stop();
}

void LinearSourceDomain::normalize_scalar_flux_and_volumes(
  double total_active_distance_per_iteration)
{
  float normalization_factor = 1.0 / total_active_distance_per_iteration;
  double volume_normalization_factor =
    1.0 / (total_active_distance_per_iteration * simulation::current_batch);

// Normalize flux to total distance travelled by all rays this iteration
#pragma omp parallel for
  for (int64_t e = 0; e < scalar_flux_new_.size(); e++) {
    scalar_flux_new_[e] *= normalization_factor;
    flux_moments_new_[e] *= normalization_factor;
  }

// Accumulate cell-wise ray length tallies collected this iteration, then
// update the simulation-averaged cell-wise volume estimates
#pragma omp parallel for
  for (int64_t sr = 0; sr < n_source_regions_; sr++) {
    volume_t_[sr] += volume_[sr];
    volume_[sr] = volume_t_[sr] * volume_normalization_factor;
  }
}

int64_t LinearSourceDomain::add_source_to_scalar_flux()
{
  int64_t n_hits = 0;

  // Temperature and angle indices, if using multiple temperature
  // data sets and/or anisotropic data sets.
  // TODO: Currently assumes we are only using single temp/single
  // angle data.
  const int t = 0;
  const int a = 0;

#pragma omp parallel for reduction(+ : n_hits)
  for (int sr = 0; sr < n_source_regions_; sr++) {

    double volume = volume_[sr];
    double volume_tracks = volume_t_[sr];
    double invvol = 1.0 / volume_tracks;
    int material = material_[sr];

    // Check if this cell was hit this iteration
    int was_cell_hit = was_hit_[sr];
    if (was_cell_hit) {
      n_hits++;
      centroid_[sr] = centroid_t_[sr];
      centroid_[sr] *= invvol;
      mom_matrix_[sr] = mom_matrix_t_[sr];
      mom_matrix_[sr] *= invvol;
    }

    for (int g = 0; g < negroups_; g++) {
      int64_t idx = (sr * negroups_) + g;
      // There are three scenarios we need to consider:
      if (was_cell_hit) {
        // 1. If the FSR was hit this iteration, then the new flux is equal to
        // the flat source from the previous iteration plus the contributions
        // from rays passing through the source region (computed during the
        // transport sweep)
        scalar_flux_new_[idx] /= volume;
        scalar_flux_new_[idx] += source_[idx];
        flux_moments_new_[idx] *= (1.0 / volume);
      } else if (volume > 0.0) {
        // 2. If the FSR was not hit this iteration, but has been hit some
        // previous iteration, then we simply set the new scalar flux to be
        // equal to the contribution from the flat source alone.
        scalar_flux_new_[idx] = source_[idx];
      } else {
        // If the FSR was not hit this iteration, and it has never been hit in
        // any iteration (i.e., volume is zero), then we want to set this to 0
        // to avoid dividing anything by a zero volume.
        scalar_flux_new_[idx] = 0.0f;
        flux_moments_new_[idx] *= 0.0f;
      }
    }
  }

  return n_hits;
}

void LinearSourceDomain::flux_swap()
{
  FlatSourceDomain::flux_swap();
  flux_moments_old_.swap(flux_moments_new_);
}

void LinearSourceDomain::accumulate_iteration_flux()
{
  // Accumulate scalar flux
  FlatSourceDomain::accumulate_iteration_flux();

  // Accumulate scalar flux moments
#pragma omp parallel for
  for (int64_t se = 0; se < n_source_elements_; se++) {
    flux_moments_t_[se] += flux_moments_new_[se];
  }
}

double LinearSourceDomain::evaluate_flux_at_point(
  const Position r, const int64_t sr, const int g) const
{
  float phi_flat = FlatSourceDomain::evaluate_flux_at_point(r, sr, g);

  Position local_r = r - centroid_[sr];
  MomentArray phi_linear = flux_moments_t_[sr * negroups_ + g];
  phi_linear *= 1.0 / (settings::n_batches - settings::n_inactive);
  MomentArray phi_solved = mom_matrix_[sr].solve(phi_linear);
  
  return phi_flat + phi_solved.dot(local_r);
}

} // namespace openmc