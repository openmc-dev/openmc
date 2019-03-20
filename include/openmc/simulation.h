//! \file simulation.h
//! \brief Variables/functions related to a running simulation

#ifndef OPENMC_SIMULATION_H
#define OPENMC_SIMULATION_H

#include "openmc/particle.h"

#include <cstdint>
#include <vector>

namespace openmc {

constexpr int STATUS_EXIT_NORMAL {0};
constexpr int STATUS_EXIT_MAX_BATCH {1};
constexpr int STATUS_EXIT_ON_TRIGGER {2};

//==============================================================================
// Global variable declarations
//==============================================================================

namespace simulation {

extern "C" int current_batch;    //!< current batch
extern "C" int current_gen;      //!< current fission generation
extern "C" int64_t current_work; //!< index in source back of current particle
extern "C" bool initialized;     //!< has simulation been initialized?
extern "C" double keff;          //!< average k over batches
extern "C" double keff_std;      //!< standard deviation of average k
extern "C" double k_col_abs;     //!< sum over batches of k_collision * k_absorption
extern "C" double k_col_tra;     //!< sum over batches of k_collision * k_tracklength
extern "C" double k_abs_tra;     //!< sum over batches of k_absorption * k_tracklength
extern double log_spacing;       //!< lethargy spacing for energy grid searches
extern "C" int n_lost_particles; //!< cumulative number of lost particles
extern "C" bool need_depletion_rx; //!< need to calculate depletion rx?
extern "C" int restart_batch;   //!< batch at which a restart job resumed
extern "C" bool satisfy_triggers; //!< have tally triggers been satisfied?
extern "C" int total_gen;        //!< total number of generations simulated
extern double total_weight;  //!< Total source weight in a batch
extern int64_t work_per_rank;         //!< number of particles per MPI rank

extern std::vector<double> k_generation;
extern std::vector<int64_t> work_index;

// Threadprivate variables
extern "C" bool trace;     //!< flag to show debug information

#pragma omp threadprivate(current_work, trace)

} // namespace simulation

//==============================================================================
// Functions
//==============================================================================

//! Allocate space for source and fission banks
void allocate_banks();

//! Determine number of particles to transport per process
void calculate_work();

//! Initialize a batch
void initialize_batch();

//! Initialize a fission generation
void initialize_generation();

void initialize_history(Particle* p, int64_t index_source);

//! Finalize a batch
//!
//! Handles synchronization and accumulation of tallies, calculation of Shannon
//! entropy, getting single-batch estimate of keff, and turning on tallies when
//! appropriate
void finalize_batch();

//! Finalize a fission generation
void finalize_generation();

//! Determine overall generation number
extern "C" int overall_generation();

#ifdef OPENMC_MPI
void broadcast_results();
#endif

void free_memory_simulation();

} // namespace openmc

#endif // OPENMC_SIMULATION_H
