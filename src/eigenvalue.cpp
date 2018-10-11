#include "openmc/eigenvalue.h"

#include "xtensor/xmath.hpp"
#include "xtensor/xtensor.hpp"
#include "xtensor/xview.hpp"

#include "openmc/capi.h"
#include "openmc/error.h"
#include "openmc/hdf5_interface.h"
#include "openmc/mesh.h"
#include "openmc/message_passing.h"
#include "openmc/random_lcg.h"
#include "openmc/search.h"
#include "openmc/settings.h"
#include "openmc/simulation.h"

#include <algorithm> // for min
#include <string>

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

std::vector<double> entropy;
xt::xtensor<double, 1> source_frac;

//==============================================================================
// Non-member functions
//==============================================================================

void synchronize_bank()
{
  // Get pointers to source/fission bank
  Bank* source_bank;
  Bank* fission_bank;
  int64_t n;
  openmc_source_bank(&source_bank, &n);
  openmc_fission_bank(&fission_bank, &n);

  // In order to properly understand the fission bank algorithm, you need to
  // think of the fission and source bank as being one global array divided
  // over multiple processors. At the start, each processor has a random amount
  // of fission bank sites -- each processor needs to know the total number of
  // sites in order to figure out the probability for selecting
  // sites. Furthermore, each proc also needs to know where in the 'global'
  // fission bank its own sites starts in order to ensure reproducibility by
  // skipping ahead to the proper seed.

#ifdef OPENMC_MPI
  int64_t start = 0;
  MPI_Exscan(&n_bank, &start, 1, MPI_INT64_T, MPI_SUM, mpi::intracomm);

  // While we would expect the value of start on rank 0 to be 0, the MPI
  // standard says that the receive buffer on rank 0 is undefined and not
  // significant
  if (mpi::rank == 0) start = 0;

  int64_t finish = start + n_bank;
  int64_t total = finish;
  MPI_Bcast(&total, 1, MPI_INT64_T, mpi::n_procs - 1, mpi::intracomm);

#else
  int64_t start  = 0;
  int64_t finish = n_bank;
  int64_t total  = n_bank;
#endif

  // If there are not that many particles per generation, it's possible that no
  // fission sites were created at all on a single processor. Rather than add
  // extra logic to treat this circumstance, we really want to ensure the user
  // runs enough particles to avoid this in the first place.

  if (n_bank == 0) {
    fatal_error("No fission sites banked on MPI rank " + std::to_string(mpi::rank));
  }

  // Make sure all processors start at the same point for random sampling. Then
  // skip ahead in the sequence using the starting index in the 'global'
  // fission bank for each processor.

  set_particle_seed(simulation::total_gen + overall_generation());
  advance_prn_seed(start);

  // Determine how many fission sites we need to sample from the source bank
  // and the probability for selecting a site.

  int64_t sites_needed;
  if (total < settings::n_particles) {
    sites_needed = settings::n_particles % total;
  } else {
    sites_needed = settings::n_particles;
  }
  double p_sample = static_cast<double>(sites_needed) / total;

  //time_bank_sample % start()

  // ==========================================================================
  // SAMPLE N_PARTICLES FROM FISSION BANK AND PLACE IN TEMP_SITES

  // Allocate temporary source bank
  int64_t index_temp = 0;
  Bank temp_sites[3*simulation::work];

  for (int64_t i = 0; i < n_bank; ++i) {
    // If there are less than n_particles particles banked, automatically add
    // int(n_particles/total) sites to temp_sites. For example, if you need
    // 1000 and 300 were banked, this would add 3 source sites per banked site
    // and the remaining 100 would be randomly sampled.
    if (total < settings::n_particles) {
      for (int64_t j = 1; j <= settings::n_particles / total; ++j) {
        temp_sites[index_temp] = fission_bank[i];
        ++index_temp;
      }
    }

    // Randomly sample sites needed
    if (prn() < p_sample) {
      temp_sites[index_temp] = fission_bank[i];
      ++index_temp;
    }
  }

  // At this point, the sampling of source sites is done and now we need to
  // figure out where to send source sites. Since it is possible that one
  // processor's share of the source bank spans more than just the immediate
  // neighboring processors, we have to perform an ALLGATHER to determine the
  // indices for all processors

#ifdef OPENMC_MPI
  // First do an exclusive scan to get the starting indices for
  start = 0;
  MPI_Exscan(&index_temp, &start, 1, MPI_INT64_T, MPI_SUM, mpi::intracomm);
  finish = start + index_temp;

  // Allocate space for bank_position if this hasn't been done yet
  int64_t bank_position[mpi::n_procs];
  MPI_Allgather(&start, 1, MPI_INT64_T, bank_position, 1,
    MPI_INT64_T, mpi::intracomm);
#else
  start = 0;
  finish = index_temp;
#endif

  // Now that the sampling is complete, we need to ensure that we have exactly
  // n_particles source sites. The way this is done in a reproducible manner is
  // to adjust only the source sites on the last processor.

  if (mpi::rank == mpi::n_procs - 1) {
    if (finish > settings::n_particles) {
      // If we have extra sites sampled, we will simply discard the extra
      // ones on the last processor
      index_temp = settings::n_particles - start;

    } else if (finish < settings::n_particles) {
      // If we have too few sites, repeat sites from the very end of the
      // fission bank
      sites_needed = settings::n_particles - finish;
      for (int i = 0; i < sites_needed; ++i) {
        temp_sites[index_temp] = fission_bank[n_bank - sites_needed + i];
        ++index_temp;
      }
    }

    // the last processor should not be sending sites to right
    finish = simulation::work_index[mpi::rank + 1];
  }

  //time_bank_sample % stop()
  //time_bank_sendrecv % start()

#ifdef OPENMC_MPI
  // ==========================================================================
  // SEND BANK SITES TO NEIGHBORS

  int64_t index_local = 0;
  std::vector<MPI_Request> requests;

  if (start < settings::n_particles) {
    // Determine the index of the processor which has the first part of the
    // source_bank for the local processor
    int neighbor = upper_bound_index(simulation::work_index.begin(),
      simulation::work_index.end(), start);

    while (start < finish) {
      // Determine the number of sites to send
      int64_t n = std::min(simulation::work_index[neighbor + 1], finish) - start;

      // Initiate an asynchronous send of source sites to the neighboring
      // process
      if (neighbor != mpi::rank) {
        requests.emplace_back();
        MPI_Isend(&temp_sites[index_local], static_cast<int>(n), mpi::bank,
          neighbor, mpi::rank, mpi::intracomm, &requests.back());
      }

      // Increment all indices
      start += n;
      index_local += n;
      ++neighbor;

      // Check for sites out of bounds -- this only happens in the rare
      // circumstance that a processor close to the end has so many sites that
      // it would exceed the bank on the last processor
      if (neighbor > mpi::n_procs - 1) break;
    }
  }

  // ==========================================================================
  // RECEIVE BANK SITES FROM NEIGHBORS OR TEMPORARY BANK

  start = simulation::work_index[mpi::rank];
  index_local = 0;

  // Determine what process has the source sites that will need to be stored at
  // the beginning of this processor's source bank.

  int neighbor;
  if (start >= bank_position[mpi::n_procs - 1]) {
    neighbor = mpi::n_procs - 1;
  } else {
    neighbor = upper_bound_index(bank_position, bank_position + mpi::n_procs, start);
  }

  while (start < simulation::work_index[mpi::rank + 1]) {
    // Determine how many sites need to be received
    int64_t n;
    if (neighbor == mpi::n_procs - 1) {
      n = simulation::work_index[mpi::rank + 1] - start;
    } else {
      n = std::min(bank_position[neighbor + 1], simulation::work_index[mpi::rank + 1]) - start;
    }

    if (neighbor != mpi::rank) {
      // If the source sites are not on this processor, initiate an
      // asynchronous receive for the source sites

      requests.emplace_back();
      MPI_Irecv(&source_bank[index_local], static_cast<int>(n), mpi::bank,
            neighbor, neighbor, mpi::intracomm, &requests.back());

    } else {
      // If the source sites are on this procesor, we can simply copy them
      // from the temp_sites bank

      index_temp = start - bank_position[mpi::rank];
      std::copy(&temp_sites[index_temp], &temp_sites[index_temp + n],
        &source_bank[index_local]);
    }

    // Increment all indices
    start += n;
    index_local += n;
    ++neighbor;
  }

  // Since we initiated a series of asynchronous ISENDs and IRECVs, now we have
  // to ensure that the data has actually been communicated before moving on to
  // the next generation

  int n_request = requests.size();
  MPI_Waitall(n_request, requests.data(), MPI_STATUSES_IGNORE);

#else
  std::copy(temp_sites, temp_sites + settings::n_particles, source_bank);
#endif

  //time_bank_sendrecv % stop()
}

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

  if (simulation::current_batch == 1 && simulation::current_gen == 1) {
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
    source_frac = m->count_sites(simulation::work, source_bank, 0, nullptr,
      &sites_outside);

    // Check for sites outside of the mesh
    if (mpi::master && sites_outside) {
      fatal_error("Source sites outside of the UFS mesh!");
    }

#ifdef OPENMC_MPI
    // Send source fraction to all processors
    int n_bins = xt::prod(m->shape_)();
    MPI_Bcast(source_frac.data(), n_bins, MPI_DOUBLE, 0, mpi::intracomm);
#endif

    // Normalize to total weight to get fraction of source in each cell
    double total = xt::sum(source_frac)();
    source_frac /= total;

    // Since the total starting weight is not equal to n_particles, we need to
    // renormalize the weight of the source sites
    for (int i = 0; i < simulation::work; ++i) {
      source_bank[i].wgt *= settings::n_particles / total;
    }
  }
}

double ufs_get_weight(const Particle* p)
{
  auto& m = meshes[settings::index_ufs_mesh];

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

extern "C" void entropy_clear()
{
  entropy.clear();
}


} // namespace openmc
