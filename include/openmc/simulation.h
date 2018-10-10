//! \file simulation.h
//! \brief Variables/functions related to a running simulation

#ifndef OPENMC_SIMULATION_H
#define OPENMC_SIMULATION_H

#include <cstdint>
#include <vector>

namespace openmc {

//==============================================================================
// Global variable declarations
//==============================================================================

namespace simulation {

extern "C" int current_batch;
extern "C" int current_gen;
extern "C" int64_t current_work;
extern "C" int n_lost_particles;
extern "C" int total_gen;
extern "C" bool trace;
extern "C" int64_t work;

#ifdef _OPENMP
extern "C" int n_threads;
#endif

extern std::vector<int64_t> work_index;

#pragma omp threadprivate(current_work, trace)

} // namespace simulation

//==============================================================================
// Functions
//==============================================================================

//! Initialize simulation
extern "C" void openmc_simulation_init_c();

//! Determine number of particles to transport per process
void calculate_work();

} // namespace openmc

#endif // OPENMC_SIMULATION_H
