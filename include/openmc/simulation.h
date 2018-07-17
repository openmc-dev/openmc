#ifndef OPENMC_SIMULATION_H
#define OPENMC_SIMULATION_H

#include <cstdint>
#include <vector>

namespace openmc {

extern "C" int openmc_current_batch;
extern "C" int openmc_current_gen;
extern "C" int64_t openmc_current_work;
extern "C" int openmc_n_lost_particles;
extern "C" int openmc_total_gen;

extern std::vector<int64_t> work_index;

#pragma omp threadprivate(openmc_current_work)

extern "C" void openmc_simulation_init_c();
void calculate_work();

} // namespace openmc

#endif // OPENMC_SIMULATION_H
