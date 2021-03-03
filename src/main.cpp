#ifdef OPENMC_MPI
#include <mpi.h>
#endif
#include "openmc/capi.h"
#include "openmc/constants.h"
#include "openmc/error.h"
#include "openmc/message_passing.h"
#include "openmc/particle_restart.h"
#include "openmc/settings.h"


int main(int argc, char* argv[]) {
  using namespace openmc;
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
    fatal_error(openmc_err_msg);
  }

  // start problem based on mode
  switch (settings::run_mode) {
    case RunMode::FIXED_SOURCE:
    case RunMode::EIGENVALUE:
      err = openmc_run();
      break;
    case RunMode::PLOTTING:
      err = openmc_plot_geometry();
      break;
    case RunMode::PARTICLE:
      if (mpi::master) run_particle_restart();
      err = 0;
      break;
    case RunMode::VOLUME:
      err = openmc_calculate_volumes();
      break;
    default:
      break;
  }
  if (err) fatal_error(openmc_err_msg);

  // Finalize and free up memory
  err = openmc_finalize();
  if (err) fatal_error(openmc_err_msg);

  // If MPI is in use and enabled, terminate it
#ifdef OPENMC_MPI
  MPI_Finalize();
#endif
}
