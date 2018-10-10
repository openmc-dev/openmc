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

extern "C" int current_batch;    //!< current batch
extern "C" int current_gen;      //!< current fission generation
extern "C" int64_t current_work; //!< index in source back of current particle
extern "C" double keff;          //!< average k over batches
extern "C" double keff_std;      //!< standard deviation of average k
extern "C" double k_col_abs;     //!< sum over batches of k_collision * k_absorption
extern "C" double k_col_tra;     //!< sum over batches of k_collision * k_tracklength
extern "C" double k_abs_tra;     //!< sum over batches of k_absorption * k_tracklength
extern "C" double log_spacing;   //!< lethargy spacing for energy grid searches
extern "C" int n_lost_particles; //!< cumulative number of lost particles
extern "C" bool need_depletion_rx; //!< need to calculate depletion rx?
extern "C" bool restart_batch;   //!< batch at which a restart job resumed
extern "C" bool satisfy_triggers; //!< have tally triggers been satisfied?
extern "C" bool simulation_initialized; //!< has simulation been initialized?
extern "C" int total_gen;        //!< total number of generations simulated
extern "C" int64_t work;         //!< number of particles per process

extern std::vector<int64_t> work_index;

// Threadprivate variables
extern "C" bool trace;     //!< flag to show debug information
#ifdef _OPENMP
extern "C" int n_threads;  //!< number of OpenMP threads
extern "C" int thread_id;  //!< ID of a given thread
#endif

#pragma omp threadprivate(current_work, thread_id, trace)

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
