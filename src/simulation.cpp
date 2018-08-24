#include "openmc/simulation.h"

#include "openmc/capi.h"
#include "openmc/message_passing.h"

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

std::vector<int64_t> work_index;

//==============================================================================
// Functions
//==============================================================================

void openmc_simulation_init_c()
{
  // Determine how much work each process should do
  calculate_work();
}

void calculate_work()
{
  // Determine minimum amount of particles to simulate on each processor
  int64_t min_work = n_particles/mpi::n_procs;

  // Determine number of processors that have one extra particle
  int64_t remainder = n_particles % mpi::n_procs;

  int64_t i_bank = 0;
  work_index.reserve(mpi::n_procs);
  work_index.push_back(0);
  for (int i = 0; i < mpi::n_procs; ++i) {
    // Number of particles for rank i
    int64_t work_i = i < remainder ? min_work + 1 : min_work;

    // Set number of particles
    if (mpi::rank == i) openmc_work = work_i;

    // Set index into source bank for rank i
    i_bank += work_i;
    work_index.push_back(i_bank);
  }
}

} // namespace openmc
