#include "openmc/tallies/pulse_height.h"

#include <algorithm> // sort
#include <cstddef>
#include <numeric> // iota
#include <utility> // make_pair

#include "openmc/message_passing.h"
#include "openmc/openmp_interface.h"
#include "openmc/particle.h"
#include "openmc/settings.h"
#include "openmc/simulation.h"
#include "openmc/tallies/tally.h"
#include "openmc/tallies/tally_scoring.h"

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

namespace simulation {

vector<vector<PulseHeightContribution>> pht_thread_buffers;

} // namespace simulation

//==============================================================================
// Non-member functions
//==============================================================================

void init_pulse_height_buffers()
{
  simulation::pht_thread_buffers.resize(num_threads());
  for (auto& buffer : simulation::pht_thread_buffers) {
    buffer.clear();
  }
}

void free_memory_pulse_height()
{
  simulation::pht_thread_buffers.clear();
  simulation::pht_thread_buffers.shrink_to_fit();
}

void stage_pulse_height(
  int64_t root_index, int64_t track_id, const vector<double>& pht)
{
  // A particle whose root was never assigned cannot be attributed to a
  // history. This should not happen, but dropping the fragment is safer than
  // adding it to an arbitrary history.
  if (root_index < 0)
    return;

  // Defensive: staging is only reachable from the shared-secondary drivers,
  // which call init_pulse_height_buffers() before transporting anything.
  if (simulation::pht_thread_buffers.empty())
    return;

  // Histories that deposit nothing still have to be scored, but they are
  // recovered from the full root range in finalize_pulse_height_tallies()
  // rather than from staged entries, so an all-zero fragment carries no
  // information and is not worth moving between ranks.
  bool nonzero = false;
  for (double e : pht) {
    if (e != 0.0) {
      nonzero = true;
      break;
    }
  }
  if (!nonzero)
    return;

  PulseHeightContribution contribution;
  contribution.root_index = root_index;
  contribution.track_id = track_id;
  contribution.energy = pht;
  simulation::pht_thread_buffers[thread_num()].push_back(
    std::move(contribution));
}

void finalize_pulse_height_tallies()
{
  int n_cells = model::pulse_height_cells.size();
  if (n_cells == 0)
    return;

  // Range of root indices owned by this rank
  int64_t first_root = phase1_first_root(mpi::rank);
  int64_t last_root = phase1_first_root(mpi::rank + 1);
  int64_t n_owned = last_root - first_root;

  // Per-history, per-cell deposited energy for the histories owned here.
  // Entries left at zero correspond to histories whose tree deposited nothing
  // in any pulse-height cell; those are still scored below, matching the
  // behaviour of the non-shared path where every primary is scored at death.
  vector<double> totals(n_owned * n_cells, 0.0);

  // Fragments are collected rather than summed as they are found, so that the
  // summation order below can be fixed by track id instead of left to arrival
  // order, which depends on thread scheduling and on where a descendant
  // landed. A track id is the particle's global slot within its generation plus
  // the tracks completed in earlier generations, both global quantities, so the
  // resulting order is the same for any thread or rank count.
  // Root and track id travel together as a pair, so one exchange carries both
  vector<int64_t> own_keys;
  vector<double> own_energy;

  // Flatten the per-thread staging buffers, keeping everything already destined
  // for this rank and packing the rest by destination.
#ifdef OPENMC_MPI
  vector<int> send_counts(mpi::n_procs, 0);
  vector<int64_t> send_keys;
  vector<double> send_energy;
  vector<vector<int64_t>> keys_by_rank(mpi::n_procs);
  vector<vector<double>> energy_by_rank(mpi::n_procs);
#endif

  for (auto& buffer : simulation::pht_thread_buffers) {
    for (auto& contribution : buffer) {
      int64_t root = contribution.root_index;
      int owner = phase1_owner_of_root(root);
      if (owner == mpi::rank) {
        own_keys.push_back(root);
        own_keys.push_back(contribution.track_id);
        own_energy.insert(own_energy.end(), contribution.energy.begin(),
          contribution.energy.end());
      } else {
#ifdef OPENMC_MPI
        keys_by_rank[owner].push_back(root);
        keys_by_rank[owner].push_back(contribution.track_id);
        energy_by_rank[owner].insert(energy_by_rank[owner].end(),
          contribution.energy.begin(), contribution.energy.end());
        send_counts[owner]++;
#endif
      }
    }
    buffer.clear();
  }

#ifdef OPENMC_MPI
  if (mpi::n_procs > 1) {
    // Concatenate the per-destination buffers into contiguous send buffers
    vector<int> send_displs(mpi::n_procs, 0);
    int total_send = 0;
    for (int r = 0; r < mpi::n_procs; ++r) {
      send_displs[r] = total_send;
      total_send += send_counts[r];
    }
    send_keys.reserve(static_cast<size_t>(total_send) * 2);
    send_energy.reserve(static_cast<size_t>(total_send) * n_cells);
    for (int r = 0; r < mpi::n_procs; ++r) {
      send_keys.insert(
        send_keys.end(), keys_by_rank[r].begin(), keys_by_rank[r].end());
      send_energy.insert(
        send_energy.end(), energy_by_rank[r].begin(), energy_by_rank[r].end());
      keys_by_rank[r].clear();
      keys_by_rank[r].shrink_to_fit();
      energy_by_rank[r].clear();
      energy_by_rank[r].shrink_to_fit();
    }

    // Exchange how many contributions each rank is sending to each other rank
    vector<int> recv_counts(mpi::n_procs, 0);
    MPI_Alltoall(send_counts.data(), 1, MPI_INT, recv_counts.data(), 1, MPI_INT,
      mpi::intracomm);

    vector<int> recv_displs(mpi::n_procs, 0);
    int total_recv = 0;
    for (int r = 0; r < mpi::n_procs; ++r) {
      recv_displs[r] = total_recv;
      total_recv += recv_counts[r];
    }

    // Both payloads are a fixed number of items per contribution, so their
    // counts and displacements are the contribution ones scaled
    auto scaled = [&](const vector<int>& v, int factor) {
      vector<int> out(v.size());
      for (size_t i = 0; i < v.size(); ++i)
        out[i] = v[i] * factor;
      return out;
    };

    vector<int64_t> recv_keys(static_cast<size_t>(total_recv) * 2);
    MPI_Alltoallv(send_keys.data(), scaled(send_counts, 2).data(),
      scaled(send_displs, 2).data(), MPI_INT64_T, recv_keys.data(),
      scaled(recv_counts, 2).data(), scaled(recv_displs, 2).data(), MPI_INT64_T,
      mpi::intracomm);

    vector<double> recv_energy(static_cast<size_t>(total_recv) * n_cells);
    MPI_Alltoallv(send_energy.data(), scaled(send_counts, n_cells).data(),
      scaled(send_displs, n_cells).data(), MPI_DOUBLE, recv_energy.data(),
      scaled(recv_counts, n_cells).data(), scaled(recv_displs, n_cells).data(),
      MPI_DOUBLE, mpi::intracomm);

    // Received fragments join the local ones, to be ordered together below
    own_keys.insert(own_keys.end(), recv_keys.begin(), recv_keys.end());
    own_energy.insert(own_energy.end(), recv_energy.begin(), recv_energy.end());
  }
#endif

  // Sum each history's fragments in track id order. Sorting by root first is
  // not needed for the result, since fragments of different histories land in
  // disjoint accumulators, but it keeps the accumulation local in memory.
  vector<int64_t> order(own_keys.size() / 2);
  std::iota(order.begin(), order.end(), int64_t {0});
  std::sort(order.begin(), order.end(), [&](int64_t a, int64_t b) {
    return std::make_pair(own_keys[2 * a], own_keys[2 * a + 1]) <
           std::make_pair(own_keys[2 * b], own_keys[2 * b + 1]);
  });

  for (int64_t i : order) {
    int64_t offset = (own_keys[2 * i] - first_root) * n_cells;
    for (int c = 0; c < n_cells; ++c) {
      totals[offset + c] += own_energy[static_cast<size_t>(i) * n_cells + c];
    }
  }

  // Score one pulse per owned history. score_pulse_height_tally() drives filter
  // matching off a Particle, so give each thread a default-constructed one; its
  // cell and E_last are overwritten and restored inside the call.
#pragma omp parallel
  {
    Particle p;
    vector<double> pht(n_cells);

#pragma omp for schedule(static)
    for (int64_t i = 0; i < n_owned; ++i) {
      for (int c = 0; c < n_cells; ++c) {
        pht[c] = totals[i * n_cells + c];
      }
      score_pulse_height_tally(p, pht, model::active_pulse_height_tallies);
    }
  }
}

} // namespace openmc
