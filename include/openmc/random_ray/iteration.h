#ifndef OPENMC_RANDOM_RAY_ITERATION_H
#define OPENMC_RANDOM_RAY_ITERATION_H

#include "openmc/random_ray/ray.h"
#include "openmc/random_ray/source_region.h"

namespace openmc {

void update_neutron_source(double k_eff);
double compute_k_eff(double k_eff_old);
void normalize_scalar_flux_and_volumes(double total_active_distance_per_iteration, int iter);
int64_t add_source_to_scalar_flux(void);
void tally_fission_rates(void);
double calculate_miss_rate(void);
int openmc_run_random_ray(void);

} // namespace openmc

#endif // OPENMC_RANDOM_RAY_ITERATION_H
