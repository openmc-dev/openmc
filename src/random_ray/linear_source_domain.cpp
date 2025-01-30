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
  for (int sr = 0; sr < n_source_regions_; sr++) {
    SourceRegion& region = source_regions_[sr];
    region.centroid_iteration_ = {0.0, 0.0, 0.0};
    region.mom_matrix_ = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    for (int g = 0; g < negroups_; g++) {
      region.flux_moments_new_[g] = {0.0, 0.0, 0.0};
    }
  }
}

void LinearSourceDomain::update_neutron_source(double k_eff)
{
  simulation::time_update_src.start();

  double inverse_k_eff = 1.0 / k_eff;

#pragma omp parallel for
  for (int sr = 0; sr < n_source_regions_; sr++) {
    SourceRegion& region = source_regions_[sr];

    int material = region.material_;
    MomentMatrix invM = region.mom_matrix_.inverse();

    for (int g_out = 0; g_out < negroups_; g_out++) {
      double sigma_t = sigma_t_[material * negroups_ + g_out];

      double scatter_flat = 0.0f;
      double fission_flat = 0.0f;
      MomentArray scatter_linear = {0.0, 0.0, 0.0};
      MomentArray fission_linear = {0.0, 0.0, 0.0};

      for (int g_in = 0; g_in < negroups_; g_in++) {
        // Handles for the flat and linear components of the flux
        double flux_flat = region.scalar_flux_old_[g_in];
        MomentArray flux_linear = region.flux_moments_old_[g_in];

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
      region.source_[g_out] =
        (scatter_flat + fission_flat * inverse_k_eff) / sigma_t;

      // Compute the linear source terms
      // In the first 10 iterations when the centroids and spatial moments
      // are not well known, we will leave the source gradients as zero
      // so as to avoid causing any numerical instability.
      if (simulation::current_batch > 10) {
        region.source_gradients_[g_out] =
          invM * ((scatter_linear + fission_linear * inverse_k_eff) / sigma_t);
      }
    }
  }

  if (settings::run_mode == RunMode::FIXED_SOURCE) {
// Add external source to flat source term if in fixed source mode
#pragma omp parallel for
    for (int sr = 0; sr < n_source_regions_; sr++) {
      SourceRegion& region = source_regions_[sr];
      for (int g = 0; g < negroups_; g++) {
        region.source_[g] += region.external_source_[g];
      }
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
  for (int sr = 0; sr < n_source_regions_; sr++) {
    SourceRegion& region = source_regions_[sr];
    for (int g = 0; g < negroups_; g++) {
      region.scalar_flux_new_[g] *= normalization_factor;
      region.flux_moments_new_[g] *= normalization_factor;
    }
  }

// Accumulate cell-wise ray length tallies collected this iteration, then
// update the simulation-averaged cell-wise volume estimates
#pragma omp parallel for
  for (int64_t sr = 0; sr < n_source_regions_; sr++) {
    SourceRegion& region = source_regions_[sr];
    region.centroid_t_ += region.centroid_iteration_;
    region.mom_matrix_t_ += region.mom_matrix_;
    region.volume_t_ += region.volume_;
    region.volume_naive_ = region.volume_ * normalization_factor;
    region.volume_ = region.volume_t_ * volume_normalization_factor;
    if (region.volume_t_ > 0.0) {
      double inv_volume = 1.0 / region.volume_t_;
      region.centroid_ = region.centroid_t_;
      region.centroid_ *= inv_volume;
      region.mom_matrix_ = region.mom_matrix_t_;
      region.mom_matrix_ *= inv_volume;
    }
  }
}

void LinearSourceDomain::set_flux_to_flux_plus_source(
  int sr, double volume, int g)
{
  SourceRegion& region = source_regions_[sr];
  region.scalar_flux_new_[g] /= volume;
  region.scalar_flux_new_[g] += region.source_[g];
  region.flux_moments_new_[g] *= (1.0 / volume);
}

void LinearSourceDomain::set_flux_to_old_flux(int sr, int g)
{
  SourceRegion& region = source_regions_[sr];
  region.scalar_flux_new_[g] = region.scalar_flux_old_[g];
  region.flux_moments_new_[g] = region.flux_moments_old_[g];
}

void LinearSourceDomain::flux_swap()
{
  FlatSourceDomain::flux_swap();
#pragma omp parallel for
  for (int sr = 0; sr < n_source_regions_; sr++) {
    SourceRegion& region = source_regions_[sr];
    region.flux_moments_old_.swap(region.flux_moments_new_);
  }
}

void LinearSourceDomain::accumulate_iteration_flux()
{
  // Accumulate scalar flux
  FlatSourceDomain::accumulate_iteration_flux();

  // Accumulate scalar flux moments
#pragma omp parallel for
  for (int sr = 0; sr < n_source_regions_; sr++) {
    SourceRegion& region = source_regions_[sr];
    for (int g = 0; g < negroups_; g++) {
      region.flux_moments_t_[g] += region.flux_moments_new_[g];
    }
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

  // Perform separate all reductions for each source region
  for (int sr = 0; sr < n_source_regions_; sr++) {
    SourceRegion& region = source_regions_[sr];

    MPI_Allreduce(MPI_IN_PLACE,
      static_cast<void*>(region.flux_moments_new_.data()), negroups_ * 3,
      MPI_DOUBLE, MPI_SUM, mpi::intracomm);
    MPI_Allreduce(MPI_IN_PLACE, static_cast<void*>(&region.mom_matrix_), 6,
      MPI_DOUBLE, MPI_SUM, mpi::intracomm);
    MPI_Allreduce(MPI_IN_PLACE, static_cast<void*>(&region.centroid_iteration_),
      3, MPI_DOUBLE, MPI_SUM, mpi::intracomm);
  }

  simulation::time_bank_sendrecv.stop();
#endif
}

double LinearSourceDomain::evaluate_flux_at_point(
  Position r, int64_t sr, int g) const
{
  const SourceRegion& region = source_regions_[sr];
  double phi_flat = FlatSourceDomain::evaluate_flux_at_point(r, sr, g);

  Position local_r = r - region.centroid_;
  MomentArray phi_linear = region.flux_moments_t_[g];
  phi_linear *= 1.0 / (settings::n_batches - settings::n_inactive);

  MomentMatrix invM = region.mom_matrix_.inverse();
  MomentArray phi_solved = invM * phi_linear;

  return phi_flat + phi_solved.dot(local_r);
}

} // namespace openmc
