#include "openmc/tallies/pulse_height.h"

#include <algorithm> // upper_bound
#include <cstddef>

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
  simulation::phase1_work_index.clear();
  simulation::phase1_work_index.shrink_to_fit();
}

void stage_pulse_height(int64_t root_index, const vector<double>& pht)
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
  contribution.energy = pht;
  simulation::pht_thread_buffers[thread_num()].push_back(
    std::move(contribution));
}

namespace {

//! Rank that owns a given root index, from the phase-1 primary partition.
int owner_of_root(int64_t root_index)
{
  const auto& index = simulation::phase1_work_index;
  auto it = std::upper_bound(index.begin(), index.end(), root_index);
  return static_cast<int>(std::distance(index.begin(), it)) - 1;
}

} // namespace

void finalize_pulse_height_tallies()
{
  int n_cells = model::pulse_height_cells.size();
  if (n_cells == 0)
    return;

  // Range of root indices owned by this rank
  int64_t first_root = simulation::phase1_work_index[mpi::rank];
  int64_t last_root = simulation::phase1_work_index[mpi::rank + 1];
  int64_t n_owned = last_root - first_root;

  // Per-history, per-cell deposited energy for the histories owned here.
  // Entries left at zero correspond to histories whose tree deposited nothing
  // in any pulse-height cell; those are still scored below, matching the
  // behaviour of the non-shared path where every primary is scored at death.
  vector<double> totals(n_owned * n_cells, 0.0);

  // Flatten the per-thread staging buffers, folding in everything already
  // destined for this rank and packing the rest by destination.
#ifdef OPENMC_MPI
  vector<int> send_counts(mpi::n_procs, 0);
  vector<int64_t> send_roots;
  vector<double> send_energy;
  vector<vector<int64_t>> roots_by_rank(mpi::n_procs);
  vector<vector<double>> energy_by_rank(mpi::n_procs);
#endif

  for (auto& buffer : simulation::pht_thread_buffers) {
    for (auto& contribution : buffer) {
      int64_t root = contribution.root_index;
      int owner = owner_of_root(root);
      if (owner == mpi::rank) {
        int64_t offset = (root - first_root) * n_cells;
        for (int c = 0; c < n_cells; ++c) {
          totals[offset + c] += contribution.energy[c];
        }
      } else {
#ifdef OPENMC_MPI
        roots_by_rank[owner].push_back(root);
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
    send_roots.reserve(total_send);
    send_energy.reserve(static_cast<size_t>(total_send) * n_cells);
    for (int r = 0; r < mpi::n_procs; ++r) {
      send_roots.insert(
        send_roots.end(), roots_by_rank[r].begin(), roots_by_rank[r].end());
      send_energy.insert(
        send_energy.end(), energy_by_rank[r].begin(), energy_by_rank[r].end());
      roots_by_rank[r].clear();
      roots_by_rank[r].shrink_to_fit();
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

    // Root indices, one per contribution
    vector<int64_t> recv_roots(total_recv);
    MPI_Alltoallv(send_roots.data(), send_counts.data(), send_displs.data(),
      MPI_INT64_T, recv_roots.data(), recv_counts.data(), recv_displs.data(),
      MPI_INT64_T, mpi::intracomm);

    // Energies, n_cells per contribution
    vector<int> send_counts_e(mpi::n_procs);
    vector<int> send_displs_e(mpi::n_procs);
    vector<int> recv_counts_e(mpi::n_procs);
    vector<int> recv_displs_e(mpi::n_procs);
    for (int r = 0; r < mpi::n_procs; ++r) {
      send_counts_e[r] = send_counts[r] * n_cells;
      send_displs_e[r] = send_displs[r] * n_cells;
      recv_counts_e[r] = recv_counts[r] * n_cells;
      recv_displs_e[r] = recv_displs[r] * n_cells;
    }
    vector<double> recv_energy(static_cast<size_t>(total_recv) * n_cells);
    MPI_Alltoallv(send_energy.data(), send_counts_e.data(),
      send_displs_e.data(), MPI_DOUBLE, recv_energy.data(),
      recv_counts_e.data(), recv_displs_e.data(), MPI_DOUBLE, mpi::intracomm);

    for (int i = 0; i < total_recv; ++i) {
      int64_t offset = (recv_roots[i] - first_root) * n_cells;
      for (int c = 0; c < n_cells; ++c) {
        totals[offset + c] += recv_energy[static_cast<size_t>(i) * n_cells + c];
      }
    }
  }
#endif

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
