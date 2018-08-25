#ifndef OPENMC_SIMULATION_H
#define OPENMC_SIMULATION_H

#include <cstdint>

extern "C" int openmc_current_batch;
extern "C" int openmc_current_gen;
extern "C" int64_t openmc_current_work;
extern "C" int openmc_n_lost_particles;

#pragma omp threadprivate(openmc_current_work)

#endif // OPENMC_SIMULATION_H
