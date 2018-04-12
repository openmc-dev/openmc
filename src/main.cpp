#include "openmc.h"

#ifdef MPI
#include <mpi.h>
#else
#define MPI_COMM_WORLD nullptr
#endif


int main(int argc, char** argv) {
  int err;

  // Initialize run -- when run with MPI, pass communicator
  openmc_init(MPI_COMM_WORLD);

  // start problem based on mode
  switch (openmc_run_mode) {
  case RUN_MODE_FIXEDSOURCE:
  case RUN_MODE_EIGENVALUE:
    err = openmc_run();
    break;
  case RUN_MODE_PLOTTING:
    openmc_plot_geometry();
    break;
  case RUN_MODE_PARTICLE:
    if (openmc_master) err = openmc_particle_restart();
    break;
  case RUN_MODE_VOLUME:
    openmc_calculate_volumes();
    break;
  }

  // Finalize and free up memory
  openmc_finalize();

  // If MPI is in use and enabled, terminate it
#ifdef MPI
  err = MPI_Finalize();
#endif
}
