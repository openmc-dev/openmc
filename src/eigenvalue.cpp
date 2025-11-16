#include "openmc/eigenvalue.h"

#include "xtensor/xbuilder.hpp"
#include "xtensor/xmath.hpp"
#include "xtensor/xtensor.hpp"
#include "xtensor/xview.hpp"

#include "openmc/array.h"
#include "openmc/bank.h"
#include "openmc/capi.h"
#include "openmc/constants.h"
#include "openmc/error.h"
#include "openmc/hdf5_interface.h"
#include "openmc/ifp.h"
#include "openmc/math_functions.h"
#include "openmc/mesh.h"
#include "openmc/message_passing.h"
#include "openmc/random_lcg.h"
#include "openmc/search.h"
#include "openmc/settings.h"
#include "openmc/simulation.h"
#include "openmc/tallies/tally.h"
#include "openmc/timer.h"

#include <algorithm> // for min
#include <cmath>     // for sqrt, abs, pow
#include <iterator>  // for back_inserter
#include <limits>    //for infinity
#include <string>

#include <fmt/core.h>

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

namespace simulation {

double keff_generation;
array<double, 2> k_sum;
vector<double> entropy;
xt::xtensor<double, 1> source_frac;

// Delayed neutron kinetics parameters
double keff_prompt_generation {0.0};
vector<double> k_prompt;
double keff_prompt {0.0};
double keff_prompt_std {0.0};
double beta_eff {0.0};
double beta_eff_std {0.0};
double alpha_k_based {0.0};
double alpha_k_based_std {0.0};
double alpha_static {0.0};
double alpha_static_std {0.0};
double prompt_gen_time {0.0};
double prompt_gen_time_std {0.0};

// Index of internal kinetics tally (for alpha calculations)
int kinetics_tally_index {-1};

// Alpha eigenvalue calculation (COG static method) - iteration state
double alpha_previous {0.0};          // Previous iteration's alpha value
double pseudo_absorption_sigma {0.0}; // Pseudo-absorption cross section
int alpha_iteration {0};              // Current alpha iteration number
bool alpha_converged {false};         // Alpha convergence flag

} // namespace simulation

//==============================================================================
// Non-member functions
//==============================================================================

void calculate_generation_keff()
{
  const auto& gt = simulation::global_tallies;

  // Get keff for this generation by subtracting off the starting value
  simulation::keff_generation =
    gt(GlobalTally::K_TRACKLENGTH, TallyResult::VALUE) -
    simulation::keff_generation;

  double keff_reduced;
#ifdef OPENMC_MPI
  if (settings::solver_type != SolverType::RANDOM_RAY) {
    // Combine values across all processors
    MPI_Allreduce(&simulation::keff_generation, &keff_reduced, 1, MPI_DOUBLE,
      MPI_SUM, mpi::intracomm);
  } else {
    // If using random ray, MPI parallelism is provided by domain replication.
    // As such, all fluxes will be reduced at the end of each transport sweep,
    // such that all ranks have identical scalar flux vectors, and will all
    // independently compute the same value of k. Thus, there is no need to
    // perform any additional MPI reduction here.
    keff_reduced = simulation::keff_generation;
  }
#else
  keff_reduced = simulation::keff_generation;
#endif

  // Normalize single batch estimate of k
  // TODO: This should be normalized by total_weight, not by n_particles
  if (settings::solver_type != SolverType::RANDOM_RAY) {
    keff_reduced /= settings::n_particles;
  }

  simulation::k_generation.push_back(keff_reduced);
}

void calculate_generation_prompt_keff()
{
  // Only calculate if enabled
  if (!settings::calculate_prompt_k)
    return;

  // Get k_prompt for this generation by subtracting off the starting value
  simulation::keff_prompt_generation =
    global_tally_prompt_tracklength - simulation::keff_prompt_generation;

  double keff_prompt_reduced;
#ifdef OPENMC_MPI
  if (settings::solver_type != SolverType::RANDOM_RAY) {
    // Combine values across all processors
    MPI_Allreduce(&simulation::keff_prompt_generation, &keff_prompt_reduced, 1,
      MPI_DOUBLE, MPI_SUM, mpi::intracomm);
  } else {
    // For random ray, all ranks have identical flux and compute the same k
    keff_prompt_reduced = simulation::keff_prompt_generation;
  }
#else
  keff_prompt_reduced = simulation::keff_prompt_generation;
#endif

  // Normalize single batch estimate of k_prompt
  if (settings::solver_type != SolverType::RANDOM_RAY) {
    keff_prompt_reduced /= settings::n_particles;
  }

  simulation::k_prompt.push_back(keff_prompt_reduced);
}

void synchronize_bank()
{
  simulation::time_bank.start();

  // In order to properly understand the fission bank algorithm, you need to
  // think of the fission and source bank as being one global array divided
  // over multiple processors. At the start, each processor has a random amount
  // of fission bank sites -- each processor needs to know the total number of
  // sites in order to figure out the probability for selecting
  // sites. Furthermore, each proc also needs to know where in the 'global'
  // fission bank its own sites starts in order to ensure reproducibility by
  // skipping ahead to the proper seed.

#ifdef OPENMC_MPI
  int64_t start = 0;
  int64_t n_bank = simulation::fission_bank.size();
  MPI_Exscan(&n_bank, &start, 1, MPI_INT64_T, MPI_SUM, mpi::intracomm);

  // While we would expect the value of start on rank 0 to be 0, the MPI
  // standard says that the receive buffer on rank 0 is undefined and not
  // significant
  if (mpi::rank == 0)
    start = 0;

  int64_t finish = start + simulation::fission_bank.size();
  int64_t total = finish;
  MPI_Bcast(&total, 1, MPI_INT64_T, mpi::n_procs - 1, mpi::intracomm);

#else
  int64_t start = 0;
  int64_t finish = simulation::fission_bank.size();
  int64_t total = finish;
#endif

  // If there are not that many particles per generation, it's possible that no
  // fission sites were created at all on a single processor. Rather than add
  // extra logic to treat this circumstance, we really want to ensure the user
  // runs enough particles to avoid this in the first place.

  if (simulation::fission_bank.size() == 0) {
    fatal_error(
      "No fission sites banked on MPI rank " + std::to_string(mpi::rank));
  }

  simulation::time_bank_sample.start();

  // Allocate temporary source bank -- we don't really know how many fission
  // sites were created, so overallocate by a factor of 3
  int64_t index_temp = 0;

  vector<SourceSite> temp_sites(3 * simulation::work_per_rank);

  // Temporary banks for IFP
  vector<vector<int>> temp_delayed_groups;
  vector<vector<double>> temp_lifetimes;
  if (settings::ifp_on) {
    resize_ifp_data(
      temp_delayed_groups, temp_lifetimes, 3 * simulation::work_per_rank);
  }

  // ==========================================================================
  // SAMPLE N_PARTICLES FROM FISSION BANK AND PLACE IN TEMP_SITES

  // We use Uniform Combing method to exactly get the targeted particle size
  // [https://doi.org/10.1080/00295639.2022.2091906]

  // Make sure all processors use the same random number seed.
  int64_t id = simulation::total_gen + overall_generation();
  uint64_t seed = init_seed(id, STREAM_TRACKING);

  // Comb specification
  double teeth_distance = static_cast<double>(total) / settings::n_particles;
  double teeth_offset = prn(&seed) * teeth_distance;

  // First and last hitting tooth
  int64_t end = start + simulation::fission_bank.size();
  int64_t tooth_start = std::ceil((start - teeth_offset) / teeth_distance);
  int64_t tooth_end = std::floor((end - teeth_offset) / teeth_distance) + 1;

  // Locally comb particles in fission_bank
  double tooth = tooth_start * teeth_distance + teeth_offset;
  for (int64_t i = tooth_start; i < tooth_end; i++) {
    int64_t idx = std::floor(tooth) - start;
    temp_sites[index_temp] = simulation::fission_bank[idx];
    if (settings::ifp_on) {
      copy_ifp_data_from_fission_banks(
        idx, temp_delayed_groups[index_temp], temp_lifetimes[index_temp]);
    }
    ++index_temp;

    // Next tooth
    tooth += teeth_distance;
  }

  // At this point, the sampling of source sites is done and now we need to
  // figure out where to send source sites. Since it is possible that one
  // processor's share of the source bank spans more than just the immediate
  // neighboring processors, we have to perform an ALLGATHER to determine the
  // indices for all processors

#ifdef OPENMC_MPI
  // First do an exclusive scan to get the starting indices for
  start = 0;
  MPI_Exscan(&index_temp, &start, 1, MPI_INT64_T, MPI_SUM, mpi::intracomm);
  finish = start + index_temp;

  // TODO: protect for MPI_Exscan at rank 0

  // Allocate space for bank_position if this hasn't been done yet
  int64_t bank_position[mpi::n_procs];
  MPI_Allgather(
    &start, 1, MPI_INT64_T, bank_position, 1, MPI_INT64_T, mpi::intracomm);
#else
  start = 0;
  finish = index_temp;
#endif

  simulation::time_bank_sample.stop();
  simulation::time_bank_sendrecv.start();

#ifdef OPENMC_MPI
  // ==========================================================================
  // SEND BANK SITES TO NEIGHBORS

  // IFP number of generation
  int ifp_n_generation;
  if (settings::ifp_on) {
    broadcast_ifp_n_generation(
      ifp_n_generation, temp_delayed_groups, temp_lifetimes);
  }

  int64_t index_local = 0;
  vector<MPI_Request> requests;

  // IFP send buffers
  vector<int> send_delayed_groups;
  vector<double> send_lifetimes;

  if (start < settings::n_particles) {
    // Determine the index of the processor which has the first part of the
    // source_bank for the local processor
    int neighbor = upper_bound_index(
      simulation::work_index.begin(), simulation::work_index.end(), start);

    // Resize IFP send buffers
    if (settings::ifp_on && mpi::n_procs > 1) {
      resize_ifp_data(send_delayed_groups, send_lifetimes,
        ifp_n_generation * 3 * simulation::work_per_rank);
    }

    while (start < finish) {
      // Determine the number of sites to send
      int64_t n =
        std::min(simulation::work_index[neighbor + 1], finish) - start;

      // Initiate an asynchronous send of source sites to the neighboring
      // process
      if (neighbor != mpi::rank) {
        requests.emplace_back();
        MPI_Isend(&temp_sites[index_local], static_cast<int>(n),
          mpi::source_site, neighbor, mpi::rank, mpi::intracomm,
          &requests.back());

        if (settings::ifp_on) {
          // Send IFP data
          send_ifp_info(index_local, n, ifp_n_generation, neighbor, requests,
            temp_delayed_groups, send_delayed_groups, temp_lifetimes,
            send_lifetimes);
        }
      }

      // Increment all indices
      start += n;
      index_local += n;
      ++neighbor;

      // Check for sites out of bounds -- this only happens in the rare
      // circumstance that a processor close to the end has so many sites that
      // it would exceed the bank on the last processor
      if (neighbor > mpi::n_procs - 1)
        break;
    }
  }

  // ==========================================================================
  // RECEIVE BANK SITES FROM NEIGHBORS OR TEMPORARY BANK

  start = simulation::work_index[mpi::rank];
  index_local = 0;

  // IFP receive buffers
  vector<int> recv_delayed_groups;
  vector<double> recv_lifetimes;
  vector<DeserializationInfo> deserialization_info;

  // Determine what process has the source sites that will need to be stored at
  // the beginning of this processor's source bank.

  int neighbor;
  if (start >= bank_position[mpi::n_procs - 1]) {
    neighbor = mpi::n_procs - 1;
  } else {
    neighbor =
      upper_bound_index(bank_position, bank_position + mpi::n_procs, start);
  }

  // Resize IFP receive buffers
  if (settings::ifp_on && mpi::n_procs > 1) {
    resize_ifp_data(recv_delayed_groups, recv_lifetimes,
      ifp_n_generation * simulation::work_per_rank);
  }

  while (start < simulation::work_index[mpi::rank + 1]) {
    // Determine how many sites need to be received
    int64_t n;
    if (neighbor == mpi::n_procs - 1) {
      n = simulation::work_index[mpi::rank + 1] - start;
    } else {
      n = std::min(bank_position[neighbor + 1],
            simulation::work_index[mpi::rank + 1]) -
          start;
    }

    if (neighbor != mpi::rank) {
      // If the source sites are not on this processor, initiate an
      // asynchronous receive for the source sites

      requests.emplace_back();
      MPI_Irecv(&simulation::source_bank[index_local], static_cast<int>(n),
        mpi::source_site, neighbor, neighbor, mpi::intracomm, &requests.back());

      if (settings::ifp_on) {
        // Receive IFP data
        receive_ifp_data(index_local, n, ifp_n_generation, neighbor, requests,
          recv_delayed_groups, recv_lifetimes, deserialization_info);
      }

    } else {
      // If the source sites are on this processor, we can simply copy them
      // from the temp_sites bank

      index_temp = start - bank_position[mpi::rank];
      std::copy(&temp_sites[index_temp], &temp_sites[index_temp + n],
        &simulation::source_bank[index_local]);

      if (settings::ifp_on) {
        copy_partial_ifp_data_to_source_banks(
          index_temp, n, index_local, temp_delayed_groups, temp_lifetimes);
      }
    }

    // Increment all indices
    start += n;
    index_local += n;
    ++neighbor;
  }

  // Since we initiated a series of asynchronous ISENDs and IRECVs, now we have
  // to ensure that the data has actually been communicated before moving on to
  // the next generation

  int n_request = requests.size();
  MPI_Waitall(n_request, requests.data(), MPI_STATUSES_IGNORE);

  if (settings::ifp_on) {
    deserialize_ifp_info(ifp_n_generation, deserialization_info,
      recv_delayed_groups, recv_lifetimes);
  }

#else
  std::copy(temp_sites.data(), temp_sites.data() + settings::n_particles,
    simulation::source_bank.begin());
  if (settings::ifp_on) {
    copy_complete_ifp_data_to_source_banks(temp_delayed_groups, temp_lifetimes);
  }
#endif

  simulation::time_bank_sendrecv.stop();
  simulation::time_bank.stop();
}

void calculate_average_keff()
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
    simulation::keff = simulation::k_generation[i];
  } else {
    // Sample mean of keff
    simulation::k_sum[0] += simulation::k_generation[i];
    simulation::k_sum[1] += std::pow(simulation::k_generation[i], 2);

    // Determine mean
    simulation::keff = simulation::k_sum[0] / n;

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
          (simulation::k_sum[1] / n - std::pow(simulation::keff, 2)) / (n - 1));
    }
  }
}

void calculate_kinetics_parameters()
{
  // Skip kinetics calculation during alpha iterations to avoid corrupting
  // the k_prompt and generation time values that were calculated from the
  // main eigenvalue batches
  if (simulation::alpha_iteration > 0)
    return;

  // Only calculate if enabled
  if (!settings::calculate_prompt_k)
    return;

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
    // For inactive generations, use current generation values as estimates
    simulation::keff_prompt = simulation::k_prompt[i];
  } else {
    // Accumulate sums for k_prompt
    static double k_prompt_sum = 0.0;
    static double k_prompt_sum_sq = 0.0;
    k_prompt_sum += simulation::k_prompt[i];
    k_prompt_sum_sq += std::pow(simulation::k_prompt[i], 2);

    // Calculate mean k_prompt
    simulation::keff_prompt = k_prompt_sum / n;

    // Calculate standard deviation if we have enough samples
    if (n > 1) {
      double t_value;
      if (settings::confidence_intervals) {
        double alpha = 1.0 - CONFIDENCE_LEVEL;
        t_value = t_percentile(1.0 - alpha / 2.0, n - 1);
      } else {
        t_value = 1.0;
      }
      simulation::keff_prompt_std =
        t_value *
        std::sqrt((k_prompt_sum_sq / n - std::pow(simulation::keff_prompt, 2)) /
                  (n - 1));
    }

    // Calculate beta_eff = (k_eff - k_prompt) / k_eff
    if (simulation::keff > 0.0) {
      simulation::beta_eff =
        (simulation::keff - simulation::keff_prompt) / simulation::keff;

      // Estimate standard deviation of beta_eff using error propagation
      // For β = (k - k_p) / k:
      // σ_β² ≈ (1/k)² σ_kp² + (k_p/k²)² σ_k²
      if (n > 1) {
        double term1 = std::pow(1.0 / simulation::keff, 2) *
                       std::pow(simulation::keff_prompt_std, 2);
        double term2 =
          std::pow(simulation::keff_prompt / std::pow(simulation::keff, 2), 2) *
          std::pow(simulation::keff_std, 2);
        simulation::beta_eff_std = std::sqrt(term1 + term2);
      }
    }

    // Calculate alpha eigenvalues if enabled and tally exists
    if (settings::calculate_alpha && simulation::kinetics_tally_index >= 0) {
      auto& tally = *model::tallies[simulation::kinetics_tally_index];
      const auto& results = tally.results();

      // Extract tally results (shape: [filter_bins, scores*nuclides,
      // result_types]) No filters (n_filter_bins=1), no nuclides, so just
      // access scores Score indices: 0=gen_time_num, 1=gen_time_denom,
      // 2=nu_fission_rate,
      //                3=absorption_rate, 4=leakage_rate, 5=population
      // Result type indices: 0=VALUE, 1=SUM, 2=SUM_SQ
      //
      // NOTE: Must use TallyResult::SUM (index 1), not VALUE (index 0)!
      // The SUM is already normalized per particle and accumulated over
      // batches. Divide by n_realizations to get the average per batch.

      int sum_idx = static_cast<int>(TallyResult::SUM);
      int sum_sq_idx = static_cast<int>(TallyResult::SUM_SQ);

      double gen_time_num = results(0, 0, sum_idx) / n;
      double gen_time_denom = results(0, 1, sum_idx) / n;
      double nu_fission_rate = results(0, 2, sum_idx) / n;
      double absorption_rate = results(0, 3, sum_idx) / n;
      double leakage_rate = results(0, 4, sum_idx) / n;
      double population = results(0, 5, sum_idx) / n;

      // Calculate standard deviations for tally scores (if n > 1)
      double gen_time_num_std = 0.0;
      double gen_time_denom_std = 0.0;
      double nu_fission_rate_std = 0.0;
      double absorption_rate_std = 0.0;
      double leakage_rate_std = 0.0;
      double population_std = 0.0;

      if (n > 1) {
        // Standard deviation of mean: σ = sqrt((SUM_SQ/n - mean²)/(n-1))
        auto calc_std = [&](int score_idx) {
          double mean = results(0, score_idx, sum_idx) / n;
          double sum_sq = results(0, score_idx, sum_sq_idx) / n;
          double variance = (sum_sq - mean * mean) / (n - 1);
          return (variance > 0.0) ? std::sqrt(variance) : 0.0;
        };

        gen_time_num_std = calc_std(0);
        gen_time_denom_std = calc_std(1);
        nu_fission_rate_std = calc_std(2);
        absorption_rate_std = calc_std(3);
        leakage_rate_std = calc_std(4);
        population_std = calc_std(5);
      }

      // Calculate prompt generation time: Λ_prompt = num / (k_prompt × denom)
      if (gen_time_denom > 0.0 && simulation::keff_prompt > 0.0) {
        simulation::prompt_gen_time =
          gen_time_num / (simulation::keff_prompt * gen_time_denom);

        // Error propagation for prompt generation time
        // For Λ = num / (k_p × denom):
        // σ_Λ² ≈ (∂Λ/∂num)² σ_num² + (∂Λ/∂denom)² σ_denom² + (∂Λ/∂k_p)² σ_kp²
        if (n > 1 && simulation::prompt_gen_time > 0.0) {
          double dLambda_dnum =
            1.0 / (simulation::keff_prompt * gen_time_denom);
          double dLambda_ddenom =
            -gen_time_num /
            (simulation::keff_prompt * gen_time_denom * gen_time_denom);
          double dLambda_dkp =
            -gen_time_num / (simulation::keff_prompt * simulation::keff_prompt *
                              gen_time_denom);

          double var_Lambda =
            dLambda_dnum * dLambda_dnum * gen_time_num_std * gen_time_num_std +
            dLambda_ddenom * dLambda_ddenom * gen_time_denom_std *
              gen_time_denom_std +
            dLambda_dkp * dLambda_dkp * simulation::keff_prompt_std *
              simulation::keff_prompt_std;

          simulation::prompt_gen_time_std = std::sqrt(var_Lambda);
        }

        // Calculate alpha (k-based): α = (k_prompt - 1) / Λ_prompt
        if (simulation::prompt_gen_time > 0.0) {
          simulation::alpha_k_based =
            (simulation::keff_prompt - 1.0) / simulation::prompt_gen_time;

          // Error propagation for alpha_k_based
          // For α = (k_p - 1) / Λ: σ_α² ≈ (1/Λ)² σ_kp² + ((k_p-1)/Λ²)² σ_Λ²
          if (n > 1 && simulation::prompt_gen_time > 0.0) {
            double dAlpha_dkp = 1.0 / simulation::prompt_gen_time;
            double dAlpha_dLambda =
              -(simulation::keff_prompt - 1.0) /
              (simulation::prompt_gen_time * simulation::prompt_gen_time);

            double var_alpha = dAlpha_dkp * dAlpha_dkp *
                                 simulation::keff_prompt_std *
                                 simulation::keff_prompt_std +
                               dAlpha_dLambda * dAlpha_dLambda *
                                 simulation::prompt_gen_time_std *
                                 simulation::prompt_gen_time_std;

            simulation::alpha_k_based_std = std::sqrt(var_alpha);
          }
        }
      }

      // ========================================================================
      // ALPHA EIGENVALUE CALCULATION (COG STATIC METHOD)
      // ========================================================================
      //
      // This method is implemented in run_alpha_iterations() which is called
      // after normal eigenvalue batches complete, based on the COG static
      // method which uses iterative refinement with pseudo-absorption.
      //
      // The method:
      //   1. Initializes α₀ = (k_prompt - 1) / Λ_prompt
      //   2. Adds pseudo-absorption σ_α = α / v to cross sections
      //   3. Runs batches and computes K'
      //   4. Updates α and iterates until K' → 1.0
      //
      // Initialize to NaN; will be set by run_alpha_iterations() if enabled
      simulation::alpha_static = std::numeric_limits<double>::quiet_NaN();
      simulation::alpha_static_std = std::numeric_limits<double>::quiet_NaN();
    }
  }
}

int openmc_get_keff(double* k_combined)
{
  k_combined[0] = 0.0;
  k_combined[1] = 0.0;

  // Special case for n <=3. Notice that at the end,
  // there is a N-3 term in a denominator.
  if (simulation::n_realizations <= 3 ||
      settings::solver_type == SolverType::RANDOM_RAY) {
    k_combined[0] = simulation::keff;
    k_combined[1] = simulation::keff_std;
    if (simulation::n_realizations <= 1) {
      k_combined[1] = std::numeric_limits<double>::infinity();
    }
    return 0;
  }

  // Initialize variables
  int64_t n = simulation::n_realizations;

  // Copy estimates of k-effective and its variance (not variance of the mean)
  const auto& gt = simulation::global_tallies;

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
  cov(0, 1) = (simulation::k_col_abs - n * kv[0] * kv[1]) / (n - 1);
  cov(0, 2) = (simulation::k_col_tra - n * kv[0] * kv[2]) / (n - 1);
  cov(1, 2) = (simulation::k_abs_tra - n * kv[1] * kv[2]) / (n - 1);
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

void shannon_entropy()
{
  // Get source weight in each mesh bin
  bool sites_outside;
  xt::xtensor<double, 1> p =
    simulation::entropy_mesh->count_sites(simulation::fission_bank.data(),
      simulation::fission_bank.size(), &sites_outside);

  // display warning message if there were sites outside entropy box
  if (sites_outside) {
    if (mpi::master)
      warning("Fission source site(s) outside of entropy box.");
  }

  if (mpi::master) {
    // Normalize to total weight of bank sites
    p /= xt::sum(p);

    // Sum values to obtain Shannon entropy
    double H = 0.0;
    for (auto p_i : p) {
      if (p_i > 0.0) {
        H -= p_i * std::log2(p_i);
      }
    }

    // Add value to vector
    simulation::entropy.push_back(H);
  }
}

void ufs_count_sites()
{
  if (simulation::current_batch == 1 && simulation::current_gen == 1) {
    // On the first generation, just assume that the source is already evenly
    // distributed so that effectively the production of fission sites is not
    // biased

    std::size_t n = simulation::ufs_mesh->n_bins();
    double vol_frac = simulation::ufs_mesh->volume_frac_;
    simulation::source_frac = xt::xtensor<double, 1>({n}, vol_frac);

  } else {
    // count number of source sites in each ufs mesh cell
    bool sites_outside;
    simulation::source_frac =
      simulation::ufs_mesh->count_sites(simulation::source_bank.data(),
        simulation::source_bank.size(), &sites_outside);

    // Check for sites outside of the mesh
    if (mpi::master && sites_outside) {
      fatal_error("Source sites outside of the UFS mesh!");
    }

#ifdef OPENMC_MPI
    // Send source fraction to all processors
    int n_bins = simulation::ufs_mesh->n_bins();
    MPI_Bcast(
      simulation::source_frac.data(), n_bins, MPI_DOUBLE, 0, mpi::intracomm);
#endif

    // Normalize to total weight to get fraction of source in each cell
    double total = xt::sum(simulation::source_frac)();
    simulation::source_frac /= total;

    // Since the total starting weight is not equal to n_particles, we need to
    // renormalize the weight of the source sites
    for (int i = 0; i < simulation::work_per_rank; ++i) {
      simulation::source_bank[i].wgt *= settings::n_particles / total;
    }
  }
}

double ufs_get_weight(const Particle& p)
{
  // Determine indices on ufs mesh for current location
  int mesh_bin = simulation::ufs_mesh->get_bin(p.r());
  if (mesh_bin < 0) {
    p.write_restart();
    fatal_error("Source site outside UFS mesh!");
  }

  if (simulation::source_frac(mesh_bin) != 0.0) {
    return simulation::ufs_mesh->volume_frac_ /
           simulation::source_frac(mesh_bin);
  } else {
    return 1.0;
  }
}

void write_eigenvalue_hdf5(hid_t group)
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
  array<double, 2> k_combined;
  openmc_get_keff(k_combined.data());
  write_dataset(group, "k_combined", k_combined);

  // Write delayed neutron kinetics parameters if calculated
  if (settings::calculate_prompt_k) {
    write_dataset(group, "k_prompt_generation", simulation::k_prompt);
    array<double, 2> k_prompt_vals {
      simulation::keff_prompt, simulation::keff_prompt_std};
    write_dataset(group, "k_prompt", k_prompt_vals);
    array<double, 2> beta_eff_vals {
      simulation::beta_eff, simulation::beta_eff_std};
    write_dataset(group, "beta_eff", beta_eff_vals);

    // Write alpha eigenvalues if calculated
    if (settings::calculate_alpha) {
      array<double, 2> prompt_gen_time_vals {
        simulation::prompt_gen_time, simulation::prompt_gen_time_std};
      write_dataset(group, "prompt_gen_time", prompt_gen_time_vals);
      array<double, 2> alpha_k_vals {
        simulation::alpha_k_based, simulation::alpha_k_based_std};
      write_dataset(group, "alpha_k_based", alpha_k_vals);
    }
  }
}

void read_eigenvalue_hdf5(hid_t group)
{
  read_dataset(group, "generations_per_batch", settings::gen_per_batch);
  int n = simulation::restart_batch * settings::gen_per_batch;
  simulation::k_generation.resize(n);
  read_dataset(group, "k_generation", simulation::k_generation);
  if (settings::entropy_on) {
    read_dataset(group, "entropy", simulation::entropy);
  }
  read_dataset(group, "k_col_abs", simulation::k_col_abs);
  read_dataset(group, "k_col_tra", simulation::k_col_tra);
  read_dataset(group, "k_abs_tra", simulation::k_abs_tra);

  // Read delayed neutron kinetics parameters if they exist
  if (settings::calculate_prompt_k && object_exists(group, "k_prompt")) {
    simulation::k_prompt.resize(n);
    read_dataset(group, "k_prompt_generation", simulation::k_prompt);
    array<double, 2> k_prompt_vals;
    read_dataset(group, "k_prompt", k_prompt_vals);
    simulation::keff_prompt = k_prompt_vals[0];
    simulation::keff_prompt_std = k_prompt_vals[1];
    array<double, 2> beta_eff_vals;
    read_dataset(group, "beta_eff", beta_eff_vals);
    simulation::beta_eff = beta_eff_vals[0];
    simulation::beta_eff_std = beta_eff_vals[1];

    // Read alpha eigenvalues if they exist
    if (settings::calculate_alpha && object_exists(group, "alpha_k_based")) {
      array<double, 2> prompt_gen_time_vals;
      read_dataset(group, "prompt_gen_time", prompt_gen_time_vals);
      simulation::prompt_gen_time = prompt_gen_time_vals[0];
      simulation::prompt_gen_time_std = prompt_gen_time_vals[1];
      array<double, 2> alpha_k_vals;
      read_dataset(group, "alpha_k_based", alpha_k_vals);
      simulation::alpha_k_based = alpha_k_vals[0];
      simulation::alpha_k_based_std = alpha_k_vals[1];
    }
  }
}

void setup_kinetics_tallies()
{
  // Only create tallies if alpha calculations are enabled
  if (!settings::calculate_alpha)
    return;

  // Create internal tally for kinetics parameters
  auto* tally = Tally::create();
  simulation::kinetics_tally_index = tally->index();
  tally->set_writable(false); // Don't write to tallies.out

  // Set scores for alpha eigenvalue calculations
  vector<std::string> scores;

  // Scores for k-based alpha: α = (k_prompt - 1) / Λ_prompt
  scores.push_back("prompt-chain-gen-time-num");
  scores.push_back("prompt-chain-gen-time-denom");

  // Scores for rate-based alpha: α = (R_prod - R_removal) / N_prompt
  scores.push_back("prompt-chain-nu-fission-rate");
  scores.push_back("prompt-chain-absorption-rate");
  scores.push_back("prompt-chain-leakage-rate");
  scores.push_back("prompt-chain-population");

  tally->set_scores(scores);

  // No filters - tally over entire geometry
  tally->set_filters({});
}

void run_alpha_iterations()
{
  using namespace openmc;

  // ============================================================================
  // ALPHA EIGENVALUE CALCULATION (COG STATIC METHOD)
  // ============================================================================
  //
  // This implements the alpha eigenvalue calculation using the COG static
  // method, which uses iterative refinement with pseudo-absorption to find
  // the alpha eigenvalue.
  //
  // Method:
  //   1. Initialize: α₀ = (k_prompt - 1) / Λ_prompt (k-based estimate)
  //   2. Add pseudo-absorption: σ_α(E) = α / v(E) to material cross sections
  //   3. Run eigenvalue batch with modified cross sections → get K'
  //   4. Update: α_new = α_old + (K' - 1) / Λ
  //   5. Iterate until convergence: |K' - 1.0| < tolerance
  //
  // The converged α satisfies K'(α) = 1, where K' is the eigenvalue of the
  // transport equation with pseudo-absorption included.
  //
  // This method runs AFTER normal eigenvalue batches complete, using the
  // converged source distribution. The alpha_iteration counter (> 0) signals
  // that pseudo-absorption should be added during transport and that kinetics
  // parameter calculation should be skipped to avoid corrupting k_prompt/Λ.
  //
  // Reference: COG User Manual, "Alpha Eigenvalue" section
  // ============================================================================

  // Only run if calculate_alpha is enabled and we're in eigenvalue mode
  if (!settings::calculate_alpha || settings::run_mode != RunMode::EIGENVALUE)
    return;

  // Verify that normal eigenvalue calculation has completed and we have
  // valid kinetics parameters from the main simulation
  if (simulation::keff_prompt <= 0.0 || simulation::prompt_gen_time <= 0.0) {
    if (mpi::master) {
      warning(
        "Cannot run alpha iterations: invalid k_prompt or generation time");
    }
    return;
  }

  // Additional safety check: ensure we have completed some active batches
  if (simulation::current_batch < settings::n_inactive + 1) {
    if (mpi::master) {
      warning("Cannot run alpha iterations: no active batches completed");
    }
    return;
  }

  if (mpi::master) {
    header("ALPHA EIGENVALUE CALCULATION (COG STATIC METHOD)", 3);
    fmt::print("\n");
    fmt::print(" Initial k_prompt        = {:.6f}\n", simulation::keff_prompt);
    fmt::print(
      " Initial gen time        = {:.6e} s\n", simulation::prompt_gen_time);
  }

  // Initialize alpha using k-based estimate: α = (k_prompt - 1) / Λ
  simulation::alpha_previous =
    (simulation::keff_prompt - 1.0) / simulation::prompt_gen_time;

  if (mpi::master) {
    fmt::print(
      " Initial alpha estimate  = {:.6e} 1/s\n\n", simulation::alpha_previous);
    fmt::print(" Iteration    Alpha (1/s)         K'          |K'-1|\n");
    fmt::print(" ---------    ------------      --------      --------\n");
  }

  // Alpha iteration loop
  simulation::alpha_converged = false;
  double k_prime = 0.0;

  for (simulation::alpha_iteration = 1;
       simulation::alpha_iteration <= settings::max_alpha_iterations;
       ++simulation::alpha_iteration) {

    // Run a single batch with pseudo-absorption enabled
    // The pseudo-absorption σ_α = α/v is added in particle.cpp
    // Note: This runs additional batches beyond n_max_batches
    int status = 0;
    openmc_next_batch(&status);

    // Get K' from this batch (track-length estimator)
    k_prime = simulation::keff_generation;

    // Check convergence: |K' - 1.0| < tolerance
    double k_error = std::abs(k_prime - 1.0);

    if (mpi::master) {
      fmt::print(" {:4d}       {:.6e}    {:.6f}    {:.6e}\n",
        simulation::alpha_iteration, simulation::alpha_previous, k_prime,
        k_error);
    }

    // Check for convergence
    if (k_error < settings::alpha_tolerance) {
      simulation::alpha_converged = true;
      if (mpi::master) {
        fmt::print(
          "\n *** CONVERGED: K' = {:.6f} (target: 1.0) ***\n", k_prime);
      }
      break;
    }

    // Update alpha for next iteration
    // Using COG's Order1 method: α_new = α_old + (K' - 1) / Λ
    double delta_alpha = (k_prime - 1.0) / simulation::prompt_gen_time;
    simulation::alpha_previous += delta_alpha;

    // Check for divergence
    if (std::abs(simulation::alpha_previous) > 1.0e10) {
      if (mpi::master) {
        warning("Alpha iteration diverging, stopping iterations");
      }
      break;
    }
  }

  // Report final result
  if (mpi::master) {
    if (!simulation::alpha_converged &&
        simulation::alpha_iteration > settings::max_alpha_iterations) {
      fmt::print(
        "\n *** WARNING: Maximum iterations reached without convergence ***\n");
    }
    fmt::print("\n Final alpha (COG static) = {:.6e} +/- N/A 1/s\n",
      simulation::alpha_previous);
    fmt::print(
      " Converged in {} iterations\n\n", simulation::alpha_iteration - 1);
  }

  // Store the converged alpha value from the COG static method
  simulation::alpha_static = simulation::alpha_previous;
  simulation::alpha_static_std = std::numeric_limits<double>::quiet_NaN();

  // Reset iteration counter to indicate we're done with alpha iterations
  simulation::alpha_iteration = 0;
}

} // namespace openmc
