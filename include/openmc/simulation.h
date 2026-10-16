//! \file simulation.h
//! \brief Variables/functions related to a running simulation

#ifndef OPENMC_SIMULATION_H
#define OPENMC_SIMULATION_H

#include "openmc/mesh.h"
#include "openmc/particle.h"
#include "openmc/shared_array.h"
#include "openmc/vector.h"

#include <cstdint>

namespace openmc {

constexpr int STATUS_EXIT_NORMAL {0};
constexpr int STATUS_EXIT_MAX_BATCH {1};
constexpr int STATUS_EXIT_ON_TRIGGER {2};

//==============================================================================
// Global variable declarations
//==============================================================================

namespace simulation {

extern int ct_current_file; //!< current collision track file index
extern int current_batch;   //!< current batch
extern int current_gen;     //!< current fission generation
extern bool initialized;    //!< has simulation been initialized?
extern double keff;         //!< average k over batches
extern double keff_std;     //!< standard deviation of average k
extern double k_col_abs;    //!< sum over batches of k_collision * k_absorption
extern double k_col_tra;    //!< sum over batches of k_collision * k_tracklength
extern double k_abs_tra;   //!< sum over batches of k_absorption * k_tracklength
extern double log_spacing; //!< lethargy spacing for energy grid searches
extern int n_lost_particles;   //!< cumulative number of lost particles
extern bool need_depletion_rx; //!< need to calculate depletion rx?
extern int restart_batch;      //!< batch at which a restart job resumed
extern bool satisfy_triggers;  //!< have tally triggers been satisfied?
extern int ssw_current_file;   //!< current surface source file
extern int total_gen;          //!< total number of generations simulated
extern double total_weight;    //!< Total source weight in a batch
extern int64_t work_per_rank;  //!< number of particles per MPI rank

extern const RegularMesh* entropy_mesh;
extern const RegularMesh* ufs_mesh;

extern vector<double> k_generation;
extern vector<int64_t> work_index;

extern int64_t
  simulation_tracks_completed; //!< Number of tracks completed on this rank

} // namespace simulation

//==============================================================================
// Functions
//==============================================================================

//! Allocate space for source and fission banks
void allocate_banks();

//! Determine number of particles to transport per process
void calculate_work(int64_t n_particles);

//! First primary index owned by a rank under the phase-1 partition
//!
//! Recomputes what calculate_work(settings::n_particles) would produce, so
//! that the primary partition is available after work_index has been
//! overwritten for a secondary generation. Valid for rank in [0, n_procs],
//! where n_procs returns the total primary count.
//!
//! \param rank MPI rank
//! \return Index of that rank's first primary
int64_t phase1_first_root(int rank);

//! Rank owning a given root index under the phase-1 partition
//!
//! \param root Root index in [0, n_particles)
//! \return Rank that transported that primary
int phase1_owner_of_root(int64_t root);

//! Replace the placement key of every site in a collected bank with the root
//! of its history
//!
//! Must be called once per generation, after the sites have been placed and
//! before any MPI migration, since the placement key indexes a bank that is
//! local to the rank that produced the sites.
//!
//! \param sites Bank whose sites were just collected
//! \param parents Bank the parents were transported from, or nullptr when the
//!   parents were the primaries
void resolve_root_indices(
  SharedArray<SourceSite>& sites, const SharedArray<SourceSite>* parents);

//! Initialize nuclear data before a simulation
void initialize_data();

//! Initialize a batch
void initialize_batch();

//! Initialize a fission generation
void initialize_generation();

//! Full initialization of a particle track
void initialize_particle_track(
  Particle& p, int64_t index_source, bool is_secondary);

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

//! Compute unique particle ID from a 1-based source index
//! \param index_source 1-based source index within this rank's work
//! \return globally unique particle ID
int64_t compute_particle_id(int64_t index_source);

//! Compute the transport RNG seed from a particle ID
//! \param particle_id the particle's globally unique ID
//! \return seed value passed to init_particle_seeds()
int64_t compute_transport_seed(int64_t particle_id);

//! Simulate a single particle history from birth to death, inclusive of any
//! secondary particles. In shared secondary mode, only a single track is
//! transported and secondaries are deposited into a shared bank instead.
void transport_history_based_single_particle(Particle& p);

//! Simulate all particle histories using history-based parallelism
void transport_history_based();

//! Simulate all particles using history-based parallelism, with a shared
//! secondary bank
void transport_history_based_shared_secondary();

//! Simulate all particle histories using event-based parallelism
void transport_event_based();

//! Simulate all particles using event-based parallelism, with a shared
//! secondary bank
void transport_event_based_shared_secondary();

} // namespace openmc

#endif // OPENMC_SIMULATION_H
