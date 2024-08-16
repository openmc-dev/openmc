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

LinearSourceDomain::LinearSourceDomain() : FlatSourceDomain()
{
  // First order spatial moment of the scalar flux
  flux_moments_old_.assign(n_source_elements_, {0.0, 0.0, 0.0});
  flux_moments_new_.assign(n_source_elements_, {0.0, 0.0, 0.0});
  flux_moments_t_.assign(n_source_elements_, {0.0, 0.0, 0.0});
  // Source gradients given by M inverse multiplied by source moments
  source_gradients_.assign(n_source_elements_, {0.0, 0.0, 0.0});

  centroid_.assign(n_source_regions_, {0.0, 0.0, 0.0});
  centroid_iteration_.assign(n_source_regions_, {0.0, 0.0, 0.0});
  centroid_t_.assign(n_source_regions_, {0.0, 0.0, 0.0});
  mom_matrix_.assign(n_source_regions_, {0.0, 0.0, 0.0, 0.0, 0.0, 0.0});
  mom_matrix_t_.assign(n_source_regions_, {0.0, 0.0, 0.0, 0.0, 0.0, 0.0});
}

void LinearSourceDomain::batch_reset()
{
  FlatSourceDomain::batch_reset();
#pragma omp parallel for
  for (int64_t se = 0; se < n_source_elements_; se++) {
    flux_moments_new_[se] = {0.0, 0.0, 0.0};
  }
#pragma omp parallel for
  for (int64_t sr = 0; sr < n_source_regions_; sr++) {
    centroid_iteration_[sr] = {0.0, 0.0, 0.0};
    mom_matrix_[sr] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
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

#pragma omp parallel for
  for (int sr = 0; sr < n_source_regions_; sr++) {

    int material = material_[sr];
    MomentMatrix invM = mom_matrix_[sr].inverse();

    for (int e_out = 0; e_out < negroups_; e_out++) {
      double sigma_t = data::mg.macro_xs_[material].get_xs(
        MgxsType::TOTAL, e_out, nullptr, nullptr, nullptr, t, a);

      double scatter_flat = 0.0f;
      double fission_flat = 0.0f;
      MomentArray scatter_linear = {0.0, 0.0, 0.0};
      MomentArray fission_linear = {0.0, 0.0, 0.0};

      for (int e_in = 0; e_in < negroups_; e_in++) {
        // Handles for the flat and linear components of the flux
        double flux_flat = scalar_flux_old_[sr * negroups_ + e_in];
        MomentArray flux_linear = flux_moments_old_[sr * negroups_ + e_in];

        // Handles for cross sections
        double sigma_s = data::mg.macro_xs_[material].get_xs(
          MgxsType::NU_SCATTER, e_in, &e_out, nullptr, nullptr, t, a);
        double nu_sigma_f = data::mg.macro_xs_[material].get_xs(
          MgxsType::NU_FISSION, e_in, nullptr, nullptr, nullptr, t, a);
        double chi = data::mg.macro_xs_[material].get_xs(
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
      // In the first 10 iterations when the centroids and spatial moments
      // are not well known, we will leave the source gradients as zero
      // so as to avoid causing any numerical instability.
      if (simulation::current_batch > 10) {
        source_gradients_[sr * negroups_ + e_out] =
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
  double normalization_factor = 1.0 / total_active_distance_per_iteration;
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
    centroid_t_[sr] += centroid_iteration_[sr];
    mom_matrix_t_[sr] += mom_matrix_[sr];
    volume_t_[sr] += volume_[sr];
    volume_naive_[sr] = volume_[sr] * normalization_factor;
    volume_[sr] = volume_t_[sr] * volume_normalization_factor;
    if (volume_t_[sr] > 0.0) {
      double inv_volume = 1.0 / volume_t_[sr];
      centroid_[sr] = centroid_t_[sr];
      centroid_[sr] *= inv_volume;
      mom_matrix_[sr] = mom_matrix_t_[sr];
      mom_matrix_[sr] *= inv_volume;
    }
  }
}

void LinearSourceDomain::set_flux_to_flux_plus_source(
  int64_t idx, double volume, int material, int g)
{
  scalar_flux_new_[idx] /= volume;
  scalar_flux_new_[idx] += source_[idx];
  flux_moments_new_[idx] *= (1.0 / volume);
}

void LinearSourceDomain::set_flux_to_old_flux(int64_t idx)
{
  scalar_flux_new_[idx] = scalar_flux_old_[idx];
  flux_moments_new_[idx] = flux_moments_old_[idx];
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

void LinearSourceDomain::all_reduce_replicated_source_regions()
{
#ifdef OPENMC_MPI
  FlatSourceDomain::all_reduce_replicated_source_regions();
  simulation::time_bank_sendrecv.start();

  // We are going to assume we can safely cast Position, MomentArray,
  // and MomentMatrix to contiguous arrays of doubles for the MPI
  // allreduce operation. This is a safe assumption as typically
  // compilers will at most pad to 8 byte boundaries. If a new FP32 MomentArray
  // type is introduced, then there will likely be padding, in which case this
  // function will need to become more complex.
  if (sizeof(MomentArray) != 3 * sizeof(double) ||
      sizeof(MomentMatrix) != 6 * sizeof(double)) {
    fatal_error("Unexpected buffer padding in linear source domain reduction.");
  }

  MPI_Allreduce(MPI_IN_PLACE, static_cast<void*>(flux_moments_new_.data()),
    n_source_elements_ * 3, MPI_DOUBLE, MPI_SUM, mpi::intracomm);
  MPI_Allreduce(MPI_IN_PLACE, static_cast<void*>(mom_matrix_.data()),
    n_source_regions_ * 6, MPI_DOUBLE, MPI_SUM, mpi::intracomm);
  MPI_Allreduce(MPI_IN_PLACE, static_cast<void*>(centroid_iteration_.data()),
    n_source_regions_ * 3, MPI_DOUBLE, MPI_SUM, mpi::intracomm);

  simulation::time_bank_sendrecv.stop();
#endif
}

double LinearSourceDomain::evaluate_flux_at_point(
  Position r, int64_t sr, int g) const
{
  double phi_flat = FlatSourceDomain::evaluate_flux_at_point(r, sr, g);

  Position local_r = r - centroid_[sr];
  MomentArray phi_linear = flux_moments_t_[sr * negroups_ + g];
  phi_linear *= 1.0 / (settings::n_batches - settings::n_inactive);

  MomentMatrix invM = mom_matrix_[sr].inverse();
  MomentArray phi_solved = invM * phi_linear;

  return phi_flat + phi_solved.dot(local_r);
}

} // namespace openmc
