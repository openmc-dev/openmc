#include "openmc/capi.h"

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
