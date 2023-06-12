#include "openmc/random_ray/iteration.h"

namespace openmc {

int openmc_run_random_ray(void)
{
  print_inputs();

  // Display header
  header("RANDOM RAY K EIGENVALUE SIMULATION", 3);
  print_columns();

  // Allocate tally results arrays if they're not allocated yet
  for (auto& t : model::tallies) {
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

  // Enable all tallies, and enforce
  // Note: Currently, only tallies of mesh type that score fission are allowed
  for( int i = 0; i < model::tallies.size(); i++ )
  {
    auto& tally {*model::tallies[i]};
    assert(tally.scores_.size() == 1 && "Only a single fission score per mesh tally is supported in random ray mode.");
    assert(tally.scores_[0] == SCORE_FISSION && "Only fission scores are supported in random ray mode.");
    assert(tally.filters().size() == 1 && "Only a single mesh filter per tally is supported in random ray mode.");
    uint64_t tally_filter_idx = tally.filters(0);
    const auto& filt_base = model::tally_filters[tally_filter_idx].get();
    auto* filt = dynamic_cast<MeshFilter*>(filt_base);
    assert(filt != NULL && "Only mesh filter types are supported in random ray mode.");

    tally.active_ = true;
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

  // Serialize Material XS data
  //prep_xs();

  openmc::simulation::time_total.start();

  // Power Iteration Loop
  for (int iter = 1; iter <= n_iters_total; iter++) {
    // Increment current batch
    simulation::current_batch++;

    // Update neutron source
    update_neutron_source(k_eff);

    // Reset scalar and volumes flux to zero
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
    if (mpi::master && settings::verbosity >= 7) {
      print_generation();
    }

    // Tally fission rates
    simulation::time_tally_fission_rates.start();
    if( iter > settings::n_inactive)
      tally_fission_rates();
    simulation::time_tally_fission_rates.stop();

    // Set phi_old = phi_new
    random_ray::scalar_flux_old.swap(random_ray::scalar_flux_new);

    // Report if source region miss rate is higher than 0.01%
    double percent_missed = (1.0 - (static_cast<double>(n_hits) / static_cast<double>(random_ray::n_source_regions))) * 100.0;
    if( percent_missed > 0.01 )
      printf(" High FSR miss rate detected (%.4lf%%)! Consider increasing ray density by adding more particles and/or active distance.\n", percent_missed);

    if (k_eff > 2.0 || k_eff < 0.25 || !(std::isfinite(k_eff))) {
      fatal_error("Instability detected");
    }
  }
}

void update_neutron_source(double k_eff)
{
  simulation::time_update_src.start();

  double inverse_k_eff = 1.0 / k_eff;
  int negroups = data::mg.num_energy_groups_;

  #pragma omp parallel for
  for (int sr = 0; sr < random_ray::n_source_regions; sr++) {
    int material = random_ray::materials[sr]; 

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
    random_ray::volume[c] = random_ray::volume_t[sr] * volume_normalization_factor;
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
    int material = random_ray::materials[sr]; 
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

    int material = random_ray::materials[sr]; 

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

void tally_fission_rates(void)
{
  return;
}

} // namespace openmc
