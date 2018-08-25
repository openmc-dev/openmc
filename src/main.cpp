#ifdef OPENMC_MPI
#include "mpi.h"
#endif
#include "openmc/capi.h"
#include "openmc/error.h"


int main(int argc, char* argv[]) {
  int err;

  // Initialize run -- when run with MPI, pass communicator
#ifdef OPENMC_MPI
  MPI_Comm world {MPI_COMM_WORLD};
  err = openmc_init(argc, argv, &world);
#else
  err = openmc_init(argc, argv, nullptr);
#endif
  if (err == -1) {
    // This happens for the -h and -v flags
    return 0;
  } else if (err) {
    openmc::fatal_error(openmc_err_msg);
  }

  // start problem based on mode
  switch (openmc_run_mode) {
    case RUN_MODE_FIXEDSOURCE:
    case RUN_MODE_EIGENVALUE:
      err = openmc_run();
      break;
    case RUN_MODE_PLOTTING:
      err = openmc_plot_geometry();
      break;
    case RUN_MODE_PARTICLE:
      if (openmc_master) err = openmc_particle_restart();
      break;
    case RUN_MODE_VOLUME:
      err = openmc_calculate_volumes();
      break;
  }
  if (err) openmc::fatal_error(openmc_err_msg);

  // Finalize and free up memory
  err = openmc_finalize();
  if (err) openmc::fatal_error(openmc_err_msg);

  // If MPI is in use and enabled, terminate it
#ifdef OPENMC_MPI
  MPI_Finalize();
#endif
}
