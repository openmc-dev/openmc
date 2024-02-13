#ifndef OPENMC_RANDOM_RAY_ITERATION_H
#define OPENMC_RANDOM_RAY_ITERATION_H

#include "openmc/random_ray/random_ray.h"
#include "openmc/random_ray/source_region.h"

namespace openmc {

void update_neutron_source(double k_eff);
double compute_k_eff(double k_eff_old);
void normalize_scalar_flux_and_volumes();
int64_t add_source_to_scalar_flux();
double calculate_miss_rate();
void openmc_run_random_ray();
void instability_check(int64_t n_hits, double k_eff, double& avg_miss_rate);
void validate_random_ray_inputs();
void all_reduce_random_ray_batch_results(bool mapped_all_tallies);
int get_random_ray_sampling_source_index();

template<typename T>
void parallel_fill(std::vector<T>& arr, T value)
{
#pragma omp parallel for schedule(static)
  for (int i = 0; i < arr.size(); i++) {
    arr[i] = value;
  }
}

} // namespace openmc

#endif // OPENMC_RANDOM_RAY_ITERATION_H
