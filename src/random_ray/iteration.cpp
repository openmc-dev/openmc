#include "openmc/random_ray/iteration.h"
#include "openmc/output.h"
#include "openmc/geometry.h"
#include "openmc/tallies/tally.h"
#include "openmc/tallies/filter.h"
#include "openmc/tallies/tally_scoring.h"
#include "openmc/simulation.h"
#include "openmc/eigenvalue.h"
#include "openmc/timer.h"
#include "openmc/mgxs_interface.h"
#include "openmc/message_passing.h"

namespace openmc {

void random_ray_tally_energy_integrated()
{
  openmc::simulation::time_tallies.start();

  int negroups = data::mg.num_energy_groups_;

  #pragma omp parallel for
  for (int sr = 0; sr < random_ray::n_source_regions; sr++) {

    double volume = random_ray::volume[sr];
    int was_cell_hit = random_ray::was_hit[sr];
    int material = random_ray::material[sr]; 

    Particle p;
    p.r() = random_ray::position[sr];  
    bool found = exhaustive_find_cell(p);
    assert(found);
    double score_flux = 0.0;
    double score_fission = 0.0;
    for (int e = 0; e < negroups; e++) {
      int64_t idx = (sr * negroups) + e;
      float flux = random_ray::scalar_flux_new[idx];
      double Sigma_f = data::mg.macro_xs_[material].get_xs(MgxsType::FISSION, e, NULL, NULL, NULL);
      score_flux += flux * volume;
      score_fission += flux * volume * Sigma_f;
    }

    for (auto i_tally : model::active_tallies) {
      Tally& tally {*model::tallies[i_tally]};

      // Initialize an iterator over valid filter bin combinations.  If there are
      // no valid combinations, use a continue statement to ensure we skip the
      // assume_separate break below.
      auto filter_iter = FilterBinIter(tally, p);
      auto end = FilterBinIter(tally, true, &p.filter_matches());
      if (filter_iter == end)
        continue;

      // Loop over filter bins.
      for (; filter_iter != end; ++filter_iter) {
        auto filter_index = filter_iter.index_;
        auto filter_weight = filter_iter.weight_;

        // Loop over scores
        for (auto score_index = 0; score_index < tally.scores_.size(); score_index++) {
          auto score_bin = tally.scores_[score_index];
          double score;

          if (score_bin == SCORE_FLUX) {
            score = score_flux;
          } else {
            score = score_fission;
          }

          #pragma omp atomic
          tally.results_(filter_index, score_index, TallyResult::VALUE) += score * filter_weight;
        }
      }
    }
    // Reset all the filter matches for the next tally event.
    for (auto& match : p.filter_matches())
      match.bins_present_ = false;
  }
  openmc::simulation::time_tallies.stop();
}

void random_ray_tally()
{
  openmc::simulation::time_tallies.start();

  int negroups = data::mg.num_energy_groups_;

  #pragma omp parallel for
  for (int sr = 0; sr < random_ray::n_source_regions; sr++) {

    double volume = random_ray::volume[sr];
    int was_cell_hit = random_ray::was_hit[sr];
    int material = random_ray::material[sr]; 

    Particle p;
    p.r() = random_ray::position[sr];  
    bool found = exhaustive_find_cell(p);
    assert(found);
    for (int e = 0; e < negroups; e++) {
      p.g() = e;
      p.g_last() = e;

      int64_t idx = (sr * negroups) + e;
      float flux = random_ray::scalar_flux_new[idx];
      double Sigma_f = data::mg.macro_xs_[material].get_xs(MgxsType::FISSION, e, NULL, NULL, NULL);
      double score_flux = flux * volume;
      double score_fission = flux * volume * Sigma_f;

      for (auto i_tally : model::active_tallies) {
        Tally& tally {*model::tallies[i_tally]};

        // Initialize an iterator over valid filter bin combinations.  If there are
        // no valid combinations, use a continue statement to ensure we skip the
        // assume_separate break below.
        auto filter_iter = FilterBinIter(tally, p);
        auto end = FilterBinIter(tally, true, &p.filter_matches());
        if (filter_iter == end)
          continue;

        // Loop over filter bins.
        for (; filter_iter != end; ++filter_iter) {
          auto filter_index = filter_iter.index_;
          auto filter_weight = filter_iter.weight_;

          // Loop over scores
          for (auto score_index = 0; score_index < tally.scores_.size(); score_index++) {
            auto score_bin = tally.scores_[score_index];
            double score;

            if (score_bin == SCORE_FLUX) {
              score = score_flux;
            } else {
              score = score_fission;
            }

            #pragma omp atomic
            tally.results_(filter_index, score_index, TallyResult::VALUE) += score * filter_weight;
          }
        }
      }
      // Reset all the filter matches for the next tally event.
      for (auto& match : p.filter_matches())
        match.bins_present_ = false;
    }
  }
  openmc::simulation::time_tallies.stop();
}

int openmc_run_random_ray(void)
{
  // Display header
  header("RANDOM RAY K EIGENVALUE SIMULATION", 3);
  print_columns();

  // Allocate tally results arrays if they're not allocated yet
  for (auto& t : model::tallies) {
    t->set_strides();
    t->init_results();
  }

  // Reset global variables -- this is done before loading state point (as that
  // will potentially populate k_generation and entropy)
  simulation::current_batch = 0;
  simulation::current_gen = 1;
  simulation::n_realizations = 0;
  simulation::k_generation.clear();
  simulation::entropy.clear();
  openmc_reset();

  //openmc_simulation_init();

  // Enable all tallies, and enforce
  for (auto& tally : model::tallies) {
    for (auto score_bin : tally->scores_) {
      if (score_bin != SCORE_FLUX && score_bin != SCORE_FISSION) {
        fatal_error("Only flux and fission scores are supported in random ray mode");
      }
    }
    tally->active_ = true;
  }

  setup_active_tallies();

  double k_eff = 1.0;

  // Intialize Cell (FSR) data
  initialize_source_regions();

  uint64_t total_geometric_intersections = 0;

  int n_iters_total = settings::n_batches;
  int n_iters_inactive = settings::n_inactive;
  int n_iters_active = n_iters_total - n_iters_inactive;

  int nrays = settings::n_particles;
  double distance_active = settings::ray_distance_active;
  double distance_inactive = settings::ray_distance_inactive;
  double total_active_distance_per_iteration = distance_active * nrays;

  openmc::simulation::time_total.start();

  // Power Iteration Loop
  for (int iter = 1; iter <= n_iters_total; iter++) {
    // Increment current batch
    simulation::current_batch++;
    //initialize_batch();

    // Update neutron source
    update_neutron_source(k_eff);

    // Reset scalar fluxes and iteration volume tallies to zero
    std::fill(random_ray::scalar_flux_new.begin(), random_ray::scalar_flux_new.end(), 0.0);
    std::fill(random_ray::volume.begin(), random_ray::volume.end(), 0.0);

    // Start timer for transport
    simulation::time_transport.start();

    // Transport Sweep
    #pragma omp parallel for schedule(runtime) reduction(+:total_geometric_intersections)
    for (int i = 0; i < nrays; i++)
    {
      Ray r;
      r.initialize_ray(i, nrays, iter);
      total_geometric_intersections += r.transport_history_based_single_ray(distance_inactive, distance_active);
    }

    // Stop timer for transport
    simulation::time_transport.stop();

    // Normalize scalar flux and update volumes
    normalize_scalar_flux_and_volumes(total_active_distance_per_iteration, iter);

    // Add source to scalar flux
    int64_t n_hits = add_source_to_scalar_flux();

    // Compute k-eff
    k_eff = compute_k_eff(k_eff);
    simulation::k_generation.push_back(k_eff);
    calculate_average_keff();

    // Output status data
    if (settings::verbosity >= 7) {
      print_generation();
    }

    // Tally fission rates
    if (iter > settings::n_inactive) {
      //random_ray_tally();
      random_ray_tally_energy_integrated();
      accumulate_tallies();
    }

    // Set phi_old = phi_new
    random_ray::scalar_flux_old.swap(random_ray::scalar_flux_new);

    // Report if source region miss rate is higher than 0.01%
    double percent_missed = (1.0 - (static_cast<double>(n_hits) / static_cast<double>(random_ray::n_source_regions))) * 100.0;
    if( percent_missed > 0.01 )
      printf(" High FSR miss rate detected (%.4lf%%)! Consider increasing ray density by adding more particles and/or active distance.\n", percent_missed);

    if (k_eff > 2.0 || k_eff < 0.25 || !(std::isfinite(k_eff))) {
      fatal_error("Instability detected");
    }

    //finalize_batch();
  }
  openmc::simulation::time_total.stop();

  // Write tally results to tallies.out
  if (settings::output_tallies) write_tallies();

  print_results_random_ray(total_geometric_intersections);

  //openmc_simulation_finalize();

  return 0;
}

void update_neutron_source(double k_eff)
{
  simulation::time_update_src.start();

  double inverse_k_eff = 1.0 / k_eff;
  int negroups = data::mg.num_energy_groups_;

  #pragma omp parallel for
  for (int sr = 0; sr < random_ray::n_source_regions; sr++) {
    int material = random_ray::material[sr]; 

    for (int energy_group_out = 0; energy_group_out < negroups; energy_group_out++) {
      float Sigma_t = data::mg.macro_xs_[material].get_xs(MgxsType::TOTAL, energy_group_out, NULL, NULL, NULL);
      float scatter_source = 0.0f;
      float fission_source = 0.0f;

      for (int energy_group_in = 0; energy_group_in < negroups; energy_group_in++) {
        float scalar_flux = random_ray::scalar_flux_old[sr * negroups + energy_group_in];
        float Sigma_s = data::mg.macro_xs_[material].get_xs(MgxsType::NU_SCATTER, energy_group_in, &energy_group_out, NULL, NULL);
        float nu_Sigma_f = data::mg.macro_xs_[material].get_xs(MgxsType::NU_FISSION, energy_group_in, NULL, NULL, NULL);
        float Chi = data::mg.macro_xs_[material].get_xs(MgxsType::CHI_PROMPT, energy_group_in, &energy_group_out, NULL, NULL);
        scatter_source += Sigma_s    * scalar_flux;
        fission_source += nu_Sigma_f * scalar_flux * Chi;
      }

      fission_source *= inverse_k_eff;
      float new_isotropic_source = (scatter_source + fission_source)  / Sigma_t;
      random_ray::source[sr * negroups + energy_group_out] = new_isotropic_source;
    }

  }

  simulation::time_update_src.stop();
}

void normalize_scalar_flux_and_volumes(double total_active_distance_per_iteration, int iter)
{
  int negroups = data::mg.num_energy_groups_;

  float  normalization_factor =        1.0 /  total_active_distance_per_iteration;
  double volume_normalization_factor = 1.0 / (total_active_distance_per_iteration * iter);

  // Normalize Scalar flux to total distance travelled by all rays this iteration
  #pragma omp parallel for
  for (auto& element : random_ray::scalar_flux_new) {
    element *= normalization_factor;
  }

  // Accumulate cell-wise ray length tallies collected this iteration, then
  // update the simulation-averaged cell-wise volume estimates
  #pragma omp parallel for
  for (int sr = 0; sr < random_ray::n_source_regions; sr++) {
    random_ray::volume_t[sr] += random_ray::volume[sr];
    random_ray::volume[sr] = random_ray::volume_t[sr] * volume_normalization_factor;
  }
}

int64_t add_source_to_scalar_flux(void)
{
  int negroups = data::mg.num_energy_groups_;

  int64_t n_hits = 0;

  #pragma omp parallel for reduction(+:n_hits)
  for (int sr = 0; sr < random_ray::n_source_regions; sr++) {
    double volume = random_ray::volume[sr];
    int was_cell_hit = random_ray::was_hit[sr];
    int material = random_ray::material[sr]; 
    for (int e = 0; e < negroups; e++) {
      int64_t idx = (sr * negroups) + e;

      // There are three scenarios we need to consider:
      if (was_cell_hit) {
        // If it was hit, then finish computing the scalar flux estimate for the iteration
        float Sigma_t = data::mg.macro_xs_[material].get_xs(MgxsType::TOTAL, e, NULL, NULL, NULL);
        random_ray::scalar_flux_new[idx] /= (Sigma_t * volume);
        random_ray::scalar_flux_new[idx] += random_ray::source[idx];
        n_hits++;
      } else if (volume == 0.0) {
        // If it was not hit this iteration, and the simulation averaged volume is still 0,
        // then it has never been hit. In this case, just set it to 0
        random_ray::scalar_flux_new[idx] = 0.f;
      } else {
        // If it was not hit this iteration, but the simulation averaged volume is nonzero,
        // then we can reuse the previous iteration's estimate of the scalar flux. This
        // induces a small amount of error, but may be beneficial for maintaining stability
        // when the ray density is extremely low.
        random_ray::scalar_flux_new[idx] = random_ray::scalar_flux_old[idx];
      }
    }
  }

  // Return the number of source regions that were hit this iteration
  return n_hits;
}

double compute_k_eff(double k_eff_old)
{
  int negroups = data::mg.num_energy_groups_;
  double fission_rate_old = 0;
  double fission_rate_new = 0;

  #pragma omp parallel for reduction(+:fission_rate_old, fission_rate_new)
  for (int sr = 0; sr < random_ray::n_source_regions; sr++) {

    // If simulation averaged volume is zero, don't include this cell
    double volume = random_ray::volume[sr];
    if (volume == 0.0) {
      continue;
    }

    int material = random_ray::material[sr]; 

    double sr_fission_source_old = 0;
    double sr_fission_source_new = 0;

    for (int e = 0; e < negroups; e++) {
      int64_t idx = (sr * negroups) + e;
      double nu_Sigma_f = data::mg.macro_xs_[material].get_xs(MgxsType::NU_FISSION, e, NULL, NULL, NULL);
      sr_fission_source_old += nu_Sigma_f * random_ray::scalar_flux_old[idx];
      sr_fission_source_new += nu_Sigma_f * random_ray::scalar_flux_new[idx];
    }

    fission_rate_old += sr_fission_source_old * volume;
    fission_rate_new += sr_fission_source_new * volume;
  }


  double k_eff_new = k_eff_old * (fission_rate_new / fission_rate_old);

  return k_eff_new;
}

} // namespace openmc
