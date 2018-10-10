#include "openmc/simulation.h"

#include "openmc/capi.h"
#include "openmc/message_passing.h"
#include "openmc/settings.h"

// OPENMC_RUN encompasses all the main logic where iterations are performed
// over the batches, generations, and histories in a fixed source or k-eigenvalue
// calculation.

int openmc_run() {
  openmc_simulation_init();

  int err = 0;
  int status = 0;
  while (status == 0 && err == 0) {
    err = openmc_next_batch(&status);
  }

  openmc_simulation_finalize();
  return err;
}

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

namespace simulation {

int current_batch;
int current_gen;
int64_t current_work;
double keff {1.0};
double keff_std;
double k_col_abs {0.0};
double k_col_tra {0.0};
double k_abs_tra {0.0};
double log_spacing;
int n_lost_particles {0};
bool need_depletion_rx {false};
int restart_batch;
bool satisfy_triggers {false};
bool simulation_initialized {false};
int total_gen {0};
int64_t work;

std::vector<double> k_generation;
std::vector<int64_t> work_index;

// Threadprivate variables
bool trace;     //!< flag to show debug information
#ifdef _OPENMP
int n_threads {-1};  //!< number of OpenMP threads
int thread_id;  //!< ID of a given thread
#endif

} // namespace simulation

//==============================================================================
// Functions
//==============================================================================

void openmc_simulation_init_c()
{
  // Determine how much work each process should do
  calculate_work();
}

int overall_generation()
{
  using namespace simulation;
  return settings::gen_per_batch*(current_batch - 1) + current_gen;
}

void calculate_work()
{
  // Determine minimum amount of particles to simulate on each processor
  int64_t min_work = settings::n_particles / mpi::n_procs;

  // Determine number of processors that have one extra particle
  int64_t remainder = settings::n_particles % mpi::n_procs;

  int64_t i_bank = 0;
  simulation::work_index.reserve(mpi::n_procs);
  simulation::work_index.push_back(0);
  for (int i = 0; i < mpi::n_procs; ++i) {
    // Number of particles for rank i
    int64_t work_i = i < remainder ? min_work + 1 : min_work;

    // Set number of particles
    if (mpi::rank == i) simulation::work = work_i;

    // Set index into source bank for rank i
    i_bank += work_i;
    simulation::work_index.push_back(i_bank);
  }
}

} // namespace openmc
