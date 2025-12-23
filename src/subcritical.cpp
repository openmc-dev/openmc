#include "openmc/subcritical.h"
#include "openmc/eigenvalue.h"
#include "openmc/math_functions.h"
#include "openmc/simulation.h"
#include "openmc/tallies/tally.h"

#include "hdf5.h"
#include <cmath>
#include <utility>

namespace openmc{
    
void convert_to_subcritical_k(double& k, double& k_std) {
    double k_sub = k / (k + 1);
    double k_sub_std = k_sub * sqrt( pow(k_std / k, 2) + pow(k_std / (k + 1), 2) );
    k = k_sub;
    k_std = k_sub_std;
}

double convert_to_subcritical_k(double k) {
    // Used to convert keff estimators in fixed source mode, which are multiplicities, specifically = M-1
    // into the corresponding estimators for subcritical k
    double k_sub = (k/simulation::total_weight) / (k/simulation::total_weight + 1) * simulation::total_weight;
    return k_sub;
}

void calculate_generation_kq()
{
  const auto& gt = simulation::global_tallies_first_gen;

  // Get keff for this generation by subtracting off the starting value
  simulation::kq_generation_val =
    gt(GlobalTally::K_TRACKLENGTH, TallyResult::VALUE) -
    simulation::kq_generation_val;

  double kq_reduced;
#ifdef OPENMC_MPI
  if (settings::solver_type != SolverType::RANDOM_RAY) {
    // Combine values across all processors
    MPI_Allreduce(&simulation::kq_generation_val, &kq_reduced, 1, MPI_DOUBLE,
      MPI_SUM, mpi::intracomm);
  } else {
    // If using random ray, MPI parallelism is provided by domain replication.
    // As such, all fluxes will be reduced at the end of each transport sweep,
    // such that all ranks have identical scalar flux vectors, and will all
    // independently compute the same value of k. Thus, there is no need to
    // perform any additional MPI reduction here.
    kq_reduced = simulation::kq_generation_val;
  }
#else
  kq_reduced = simulation::kq_generation_val;
#endif

  // Normalize single batch estimate of k
  // TODO: This should be normalized by total_weight, not by n_particles
  if (settings::solver_type != SolverType::RANDOM_RAY) {
    kq_reduced /= settings::n_particles;
  }

  simulation::kq_generation.push_back(kq_reduced);
}

void calculate_average_kq()
{
  // Determine overall generation and number of active generations
  int i = overall_generation() - 1;
  int n;
  if (simulation::current_batch > settings::n_inactive) {
    n = settings::gen_per_batch * simulation::n_realizations +
        simulation::current_gen;
  } else {
    n = 0;
  }

  if (n <= 0) {
    // For inactive generations, use current generation k as estimate for next
    // generation
    simulation::kq = simulation::kq_generation[i];
  } else {
    // Sample mean of keff
    simulation::kq_sum[0] += simulation::kq_generation[i];
    simulation::kq_sum[1] += std::pow(simulation::kq_generation[i], 2);
    // Determine mean
    simulation::kq = simulation::kq_sum[0] / n;

    if (n > 1) {
      double t_value;
      if (settings::confidence_intervals) {
        // Calculate t-value for confidence intervals
        double alpha = 1.0 - CONFIDENCE_LEVEL;
        t_value = t_percentile(1.0 - alpha / 2.0, n - 1);
      } else {
        t_value = 1.0;
      }

      // Standard deviation of the sample mean of k
      simulation::keff_std =
        t_value *
        std::sqrt(
          (simulation::kq_sum[1] / n - std::pow(simulation::kq, 2)) / (n - 1));

      // In some cases (such as an infinite medium problem), random ray
      // may estimate k exactly and in an unvarying manner between iterations.
      // In this case, the floating point roundoff between the division and the
      // power operations may cause an extremely small negative value to occur
      // inside the sqrt operation, leading to NaN. If this occurs, we check for
      // it and set the std dev to zero.
      if (!std::isfinite(simulation::kq_std)) {
        simulation::kq_std = 0.0;
      }
    }
  }
}

int openmc_get_subcritical_kq(double* k_combined)
{
  k_combined[0] = 0.0;
  k_combined[1] = 0.0;

  // Special case for n <=3. Notice that at the end,
  // there is a N-3 term in a denominator.
  if (simulation::n_realizations <= 3 ||
      settings::solver_type == SolverType::RANDOM_RAY) {
    k_combined[0] = simulation::kq;
    k_combined[1] = simulation::kq_std;
    if (simulation::n_realizations <= 1) {
      k_combined[1] = std::numeric_limits<double>::infinity();
    }
    return 0;
  }

  // Initialize variables
  int64_t n = simulation::n_realizations;

  // Copy estimates of k-effective and its variance (not variance of the mean)
  const auto& gt = simulation::global_tallies_first_gen;

  array<double, 3> kv {};
  xt::xtensor<double, 2> cov = xt::zeros<double>({3, 3});
  kv[0] = gt(GlobalTally::K_COLLISION, TallyResult::SUM) / n;
  kv[1] = gt(GlobalTally::K_ABSORPTION, TallyResult::SUM) / n;
  kv[2] = gt(GlobalTally::K_TRACKLENGTH, TallyResult::SUM) / n;
  cov(0, 0) =
    (gt(GlobalTally::K_COLLISION, TallyResult::SUM_SQ) - n * kv[0] * kv[0]) /
    (n - 1);
  cov(1, 1) =
    (gt(GlobalTally::K_ABSORPTION, TallyResult::SUM_SQ) - n * kv[1] * kv[1]) /
    (n - 1);
  cov(2, 2) =
    (gt(GlobalTally::K_TRACKLENGTH, TallyResult::SUM_SQ) - n * kv[2] * kv[2]) /
    (n - 1);

  // Calculate covariances based on sums with Bessel's correction
  cov(0, 1) = (simulation::kq_col_abs - n * kv[0] * kv[1]) / (n - 1);
  cov(0, 2) = (simulation::kq_col_tra - n * kv[0] * kv[2]) / (n - 1);
  cov(1, 2) = (simulation::kq_abs_tra - n * kv[1] * kv[2]) / (n - 1);
  cov(1, 0) = cov(0, 1);
  cov(2, 0) = cov(0, 2);
  cov(2, 1) = cov(1, 2);

  // Check to see if two estimators are the same; this is guaranteed to happen
  // in MG-mode with survival biasing when the collision and absorption
  // estimators are the same, but can theoretically happen at anytime.
  // If it does, the standard estimators will produce floating-point
  // exceptions and an expression specifically derived for the combination of
  // two estimators (vice three) should be used instead.

  // First we will identify if there are any matching estimators
  int i, j;
  bool use_three = false;
  if ((std::abs(kv[0] - kv[1]) / kv[0] < FP_REL_PRECISION) &&
      (std::abs(cov(0, 0) - cov(1, 1)) / cov(0, 0) < FP_REL_PRECISION)) {
    // 0 and 1 match, so only use 0 and 2 in our comparisons
    i = 0;
    j = 2;

  } else if ((std::abs(kv[0] - kv[2]) / kv[0] < FP_REL_PRECISION) &&
             (std::abs(cov(0, 0) - cov(2, 2)) / cov(0, 0) < FP_REL_PRECISION)) {
    // 0 and 2 match, so only use 0 and 1 in our comparisons
    i = 0;
    j = 1;

  } else if ((std::abs(kv[1] - kv[2]) / kv[1] < FP_REL_PRECISION) &&
             (std::abs(cov(1, 1) - cov(2, 2)) / cov(1, 1) < FP_REL_PRECISION)) {
    // 1 and 2 match, so only use 0 and 1 in our comparisons
    i = 0;
    j = 1;

  } else {
    // No two estimators match, so set boolean to use all three estimators.
    use_three = true;
  }

  if (use_three) {
    // Use three estimators as derived in the paper by Urbatsch

    // Initialize variables
    double g = 0.0;
    array<double, 3> S {};

    for (int l = 0; l < 3; ++l) {
      // Permutations of estimates
      int k;
      switch (l) {
      case 0:
        // i = collision, j = absorption, k = tracklength
        i = 0;
        j = 1;
        k = 2;
        break;
      case 1:
        // i = absortion, j = tracklength, k = collision
        i = 1;
        j = 2;
        k = 0;
        break;
      case 2:
        // i = tracklength, j = collision, k = absorption
        i = 2;
        j = 0;
        k = 1;
        break;
      }

      // Calculate weighting
      double f = cov(j, j) * (cov(k, k) - cov(i, k)) - cov(k, k) * cov(i, j) +
                 cov(j, k) * (cov(i, j) + cov(i, k) - cov(j, k));

      // Add to S sums for variance of combined estimate
      S[0] += f * cov(0, l);
      S[1] += (cov(j, j) + cov(k, k) - 2.0 * cov(j, k)) * kv[l] * kv[l];
      S[2] += (cov(k, k) + cov(i, j) - cov(j, k) - cov(i, k)) * kv[l] * kv[j];

      // Add to sum for combined k-effective
      k_combined[0] += f * kv[l];
      g += f;
    }

    // Complete calculations of S sums
    for (auto& S_i : S) {
      S_i *= (n - 1);
    }
    S[0] *= (n - 1) * (n - 1);

    // Calculate combined estimate of k-effective
    k_combined[0] /= g;

    // Calculate standard deviation of combined estimate
    g *= (n - 1) * (n - 1);
    k_combined[1] =
      std::sqrt(S[0] / (g * n * (n - 3)) * (1 + n * ((S[1] - 2 * S[2]) / g)));

  } else {
    // Use only two estimators
    // These equations are derived analogously to that done in the paper by
    // Urbatsch, but are simpler than for the three estimators case since the
    // block matrices of the three estimator equations reduces to scalars here

    // Store the commonly used term
    double f = kv[i] - kv[j];
    double g = cov(i, i) + cov(j, j) - 2.0 * cov(i, j);

    // Calculate combined estimate of k-effective
    k_combined[0] = kv[i] - (cov(i, i) - cov(i, j)) / g * f;

    // Calculate standard deviation of combined estimate
    k_combined[1] = (cov(i, i) * cov(j, j) - cov(i, j) * cov(i, j)) *
                    (g + n * f * f) / (n * (n - 2) * g * g);
    k_combined[1] = std::sqrt(k_combined[1]);
  }
  return 0;
}

double calculate_ks(double k, double kq)
{
  // Calculate ks from k = kq/(1 - ks + kq) -> ks = 1 + kq*(k - 1)/k
  double ks = 1 + kq * (k - 1) / k;
  return ks;
}

double calculate_sigma_ks(double k, double k_std, double kq, double kq_std)
{
  double sigma_ks =
    sqrt(pow(kq_std, 2) +
         pow(kq / k, 2) * (pow(kq_std / kq, 2) + pow(k_std / k, 2)));
  return sigma_ks;
}

void write_subcritical_hdf5(hid_t group)
{
  write_dataset(group, "n_inactive", settings::n_inactive);
  write_dataset(group, "generations_per_batch", settings::gen_per_batch);
  write_dataset(group, "k_generation", simulation::k_generation);
  if (settings::entropy_on) {
    write_dataset(group, "entropy", simulation::entropy);
  }
  write_dataset(group, "k_col_abs", simulation::k_col_abs);
  write_dataset(group, "k_col_tra", simulation::k_col_tra);
  write_dataset(group, "k_abs_tra", simulation::k_abs_tra);
  write_dataset(group, "kq_col_abs", simulation::kq_col_abs);
  write_dataset(group, "kq_col_tra", simulation::kq_col_tra);
  write_dataset(group, "kq_abs_tra", simulation::kq_abs_tra);
  array<double, 2> k_combined;
  openmc_get_keff(k_combined.data());
  array<double, 2> kq_combined;
  openmc_get_subcritical_kq(kq_combined.data());
  write_dataset(group, "k_combined", k_combined);
  write_dataset(group, "kq_combined", kq_combined);
  array<double, 2> ks_combined;
  ks_combined[0] = calculate_ks(k_combined[0], kq_combined[0]);
  ks_combined[1] = calculate_sigma_ks(
    k_combined[0], k_combined[1], kq_combined[0], kq_combined[1]);
  write_dataset(group, "ks_combined", ks_combined);
  write_dataset(group, "kq_generation", simulation::kq_generation);
  write_dataset(group, "ks_generation", simulation::ks_generation);
}

} // namespace openmc