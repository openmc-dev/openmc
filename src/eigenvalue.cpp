#include "openmc/eigenvalue.h"

#include "xtensor/xmath.hpp"
#include "xtensor/xtensor.hpp"
#include "xtensor/xview.hpp"

#include "openmc/capi.h"
#include "openmc/error.h"
#include "openmc/hdf5_interface.h"
#include "openmc/mesh.h"
#include "openmc/message_passing.h"
#include "openmc/settings.h"
#include "openmc/simulation.h"

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

std::vector<double> entropy;
xt::xtensor<double, 1> source_frac;

//==============================================================================
// Non-member functions
//==============================================================================

void shannon_entropy()
{
  // Get pointer to entropy mesh
  auto& m = meshes[settings::index_entropy_mesh];

  // Get pointer to fission bank
  Bank* fission_bank;
  int64_t n;
  openmc_fission_bank(&fission_bank, &n);

  // Get source weight in each mesh bin
  bool sites_outside;
  xt::xtensor<double, 1> p = m->count_sites(
    n_bank, fission_bank, 0, nullptr, &sites_outside);

  // display warning message if there were sites outside entropy box
  if (sites_outside) {
    if (mpi::master) warning("Fission source site(s) outside of entropy box.");
  }

  // sum values to obtain shannon entropy
  if (mpi::master) {
    // Normalize to total weight of bank sites
    p /= xt::sum(p);

    double H = 0.0;
    for (auto p_i : p) {
      if (p_i > 0.0) {
        H -= p_i * std::log(p_i)/std::log(2.0);
      }
    }

    // Add value to vector
    entropy.push_back(H);
  }
}

void ufs_count_sites()
{
  auto &m = meshes[settings::index_ufs_mesh];

  if (openmc_current_batch == 1 && openmc_current_gen == 1) {
    // On the first generation, just assume that the source is already evenly
    // distributed so that effectively the production of fission sites is not
    // biased

    auto s = xt::view(source_frac, xt::all());
    s = m->volume_frac_;

  } else {
    // Get pointer to source bank
    Bank* source_bank;
    int64_t n;
    openmc_source_bank(&source_bank, &n);


    // count number of source sites in each ufs mesh cell
    bool sites_outside;
    source_frac = m->count_sites(openmc_work, source_bank, 0, nullptr,
      &sites_outside);

    // Check for sites outside of the mesh
    if (mpi::master && sites_outside) {
      fatal_error("Source sites outside of the UFS mesh!");
    }

#ifdef OPENMC_MPI
    // Send source fraction to all processors
    int n = xt::prod(m->shape_)();
    MPI_Bcast(source_frac.data(), n, MPI_DOUBLE, 0, mpi::intracomm)
#endif

    // Normalize to total weight to get fraction of source in each cell
    double total = xt::sum(source_frac)();
    source_frac /= total;

    // Since the total starting weight is not equal to n_particles, we need to
    // renormalize the weight of the source sites
    for (int i = 0; i < openmc_work; ++i) {
      source_bank[i].wgt *= settings::n_particles / total;
    }
  }
}

double ufs_get_weight(const Particle* p)
{
  auto& m = meshes[settings::index_entropy_mesh];

  // Determine indices on ufs mesh for current location
  // TODO: off by one
  int mesh_bin = m->get_bin({p->coord[0].xyz}) - 1;
  if (mesh_bin < 0) {
    p->write_restart();
    fatal_error("Source site outside UFS mesh!");
  }

  if (source_frac(mesh_bin) != 0.0) {
    return m->volume_frac_ / source_frac(mesh_bin);
  } else {
    return 1.0;
  }
}

extern "C" void entropy_to_hdf5(hid_t group)
{
  if (settings::entropy_on) {
    write_dataset(group, "entropy", entropy);
  }
}

extern "C" void entropy_from_hdf5(hid_t group)
{
  if (settings::entropy_on) {
    read_dataset(group, "entropy", entropy);
  }
}

extern "C" double entropy_c(int i)
{
  return entropy.at(i - 1);
}

extern "C" double entropy_clear()
{
  entropy.clear();
}


} // namespace openmc
