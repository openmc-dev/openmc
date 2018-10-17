#include <algorithm> // for copy
#include <cstdint>
#include <iostream>

#include "xtensor/xarray.hpp"
#include "xtensor/xio.hpp"

#include "openmc/capi.h"
#include "openmc/mesh.h"
#include "openmc/message_passing.h"
#include "openmc/simulation.h"

namespace openmc {

extern "C" int index_cmfd_mesh;

extern "C" void
cmfd_populate_sourcecounts(int n_energy, const double* energies,
  double* source_counts, bool* outside)
{
  // Get pointer to source bank
  Bank* source_bank;
  int64_t n;
  openmc_source_bank(&source_bank, &n);

  // Get source counts in each mesh bin / energy bin
  auto& m = meshes.at(index_cmfd_mesh);
  xt::xarray<double> counts = m->count_sites(simulation::work, source_bank, n_energy, energies, outside);

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
