#include <algorithm> // for copy
#include <cstdint>
#include <iostream>

#include "xtensor/xarray.hpp"
#include "xtensor/xio.hpp"

#include "openmc/bank.h"
#include "openmc/capi.h"
#include "openmc/mesh.h"
#include "openmc/message_passing.h"
#include "openmc/simulation.h"
#include "openmc/settings.h"

namespace openmc {


extern "C" void
cmfd_populate_sourcecounts(int n_energy, const double* energies,
  double* source_counts, bool* outside)
{
  // Get source counts in each mesh bin / energy bin
  auto& m = model::meshes.at(settings::index_cmfd_mesh);
  xt::xarray<double> counts = m->count_sites(simulation::work,
    simulation::source_bank.data(), n_energy, energies, outside);

  // Copy data from the xarray into the source counts array
  std::copy(counts.begin(), counts.end(), source_counts);
}

#ifdef OPENMC_MPI
extern "C" void
cmfd_broadcast(int n, double* buffer)
{
  MPI_Bcast(buffer, n, MPI_DOUBLE, 0, mpi::intracomm);
}
#endif


} // namespace openmc
