#include "openmc/bank.h"
#include "openmc/capi.h"
#include "openmc/error.h"
#include "openmc/ifp.h"
#include "openmc/message_passing.h"
#include "openmc/simulation.h"
#include "openmc/vector.h"

#include <cstdint>

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

namespace simulation {

vector<SourceSite> source_bank;

SharedArray<SourceSite> surf_source_bank;

SharedArray<CollisionTrackSite> collision_track_bank;

// The fission bank is allocated as a SharedArray, rather than a vector, as it
// will be shared by all threads in the simulation. It will be allocated to a
// fixed maximum capacity in the init_fission_bank() function. Then, Elements
// will be added to it by using SharedArray's special thread_safe_append()
// function.
SharedArray<SourceSite> fission_bank;

vector<vector<int>> ifp_source_delayed_group_bank;

vector<vector<double>> ifp_source_lifetime_bank;

vector<vector<int>> ifp_fission_delayed_group_bank;

vector<vector<double>> ifp_fission_lifetime_bank;

// Each entry in this vector corresponds to the number of progeny produced
// this generation for the particle located at that index. This vector is
// used to efficiently sort the fission bank after each iteration.
vector<int64_t> progeny_per_particle;

// A shared bank for secondary particles generated during transport
vector<SourceSite> global_secondary_bank;

} // namespace simulation

//==============================================================================
// Non-member functions
//==============================================================================

void free_memory_bank()
{
  simulation::source_bank.clear();
  simulation::surf_source_bank.clear();
  simulation::collision_track_bank.clear();
  simulation::fission_bank.clear();
  simulation::progeny_per_particle.clear();
  simulation::ifp_source_delayed_group_bank.clear();
  simulation::ifp_source_lifetime_bank.clear();
  simulation::ifp_fission_delayed_group_bank.clear();
  simulation::ifp_fission_lifetime_bank.clear();
}

void init_fission_bank(int64_t max)
{
  simulation::fission_bank.reserve(max);
  simulation::progeny_per_particle.resize(simulation::work_per_rank);
}

// Performs an O(n) sort on a fission or secondary bank, by leveraging
// the parent_id and progeny_id fields of banked particles. See the following
// paper for more details:
// "Reproducibility and Monte Carlo Eigenvalue Calculations," F.B. Brown and
// T.M. Sutton, 1992 ANS Annual Meeting, Transactions of the American Nuclear
// Society, Volume 65, Page 235.
void sort_bank(SharedArray<SourceSite>& bank, bool is_fission_bank)
{
  // Debugging Sanity Check
  // TODO: delete this after debugging
  int64_t n_progeny = 0;
  for (int i = 0; i < simulation::progeny_per_particle.size(); i++) {
    n_progeny += simulation::progeny_per_particle[i];
  }
  int64_t bank_size = bank.size();
  if (n_progeny != bank_size) {
    fatal_error(fmt::format(
      "Discrepancy detected in sort_bank: total progeny count ({}) does not "
      "match bank size ({}). This may indicate a bug.",
      n_progeny, bank_size));
  }

  // Ensure we don't read off the end of the array if we ran with 0 particles
  if (simulation::progeny_per_particle.size() == 0) {
    return;
  }

  // Perform exclusive scan summation to determine starting indices in fission
  // bank for each parent particle id
  std::exclusive_scan(simulation::progeny_per_particle.begin(),
    simulation::progeny_per_particle.end(),
    simulation::progeny_per_particle.begin(), 0);

  // We need a scratch vector to make permutation of the fission bank into
  // sorted order easy. Under normal usage conditions, the fission bank is
  // over provisioned, so we can use that as scratch space.
  SourceSite* sorted_bank;
  vector<SourceSite> sorted_bank_holder;
  vector<vector<int>> sorted_ifp_delayed_group_bank;
  vector<vector<double>> sorted_ifp_lifetime_bank;

  // If there is not enough space, allocate a temporary vector and point to it
  if (bank.size() > bank.capacity() / 2) {
    sorted_bank_holder.resize(bank.size());
    sorted_bank = sorted_bank_holder.data();
  } else { // otherwise, point sorted_bank to unused portion of the fission bank
    sorted_bank = &bank[bank.size()];
  }

  if (settings::ifp_on && is_fission_bank) {
    allocate_temporary_vector_ifp(
      sorted_ifp_delayed_group_bank, sorted_ifp_lifetime_bank);
  }

  // Use parent and progeny indices to sort fission bank
  for (int64_t i = 0; i < bank.size(); i++) {
    const auto& site = bank[i];
    int64_t idx =
      simulation::progeny_per_particle[site.parent_id] + site.progeny_id;
    sorted_bank[idx] = site;
    if (settings::ifp_on && is_fission_bank) {
      copy_ifp_data_from_fission_banks(
        i, sorted_ifp_delayed_group_bank[idx], sorted_ifp_lifetime_bank[idx]);
    }
  }

  // Copy sorted bank into the fission bank
  std::copy(sorted_bank, sorted_bank + bank.size(), bank.data());
  if (settings::ifp_on && is_fission_bank) {
    copy_ifp_data_to_fission_banks(
      sorted_ifp_delayed_group_bank.data(), sorted_ifp_lifetime_bank.data());
  }

  // Set the parent_id fields to something larger per rank
  // TODO: DELETE THIS AFTER DEBUGGING
  /*
  for (int64_t i = 0; i < bank.size(); i++) {
    bank[i].parent_id += mpi::rank * 10000000;
  }
    */
}

// This function redistributes SourceSite particles across MPI ranks to
// achieve load balancing while preserving the global ordering of particles.
//
// GUARANTEES:
// -----------
// 1. Global Order Preservation: After redistribution, each rank holds a
//    contiguous slice of the original global ordering. For example, if the
//    input across 3 ranks was:
//      - Rank 0: IDs 0-4
//      - Rank 1: IDs 5-6
//      - Rank 2: IDs 7-200
//    Then after redistribution (assuming ~67 particles per rank):
//      - Rank 0: IDs 0-66   (contiguous)
//      - Rank 1: IDs 67-133 (contiguous)
//      - Rank 2: IDs 134-200 (contiguous)
//    The global ordering is always preserved - no rank will ever hold
//    non-contiguous ID ranges like "0-4 and 100-200".
//
// 2. Even Load Balancing: Particles are distributed as evenly as possible.
//    If total % n_procs != 0, the first 'remainder' ranks each get one extra
//    particle (i.e., floor division with remainder distributed to lower
//    ranks). This follows the same logic as calculate_work().
//
// HOW IT WORKS:
// -------------
// The algorithm uses overlap-based redistribution:
// 1. Each rank's current data occupies a range [cumulative_before[rank],
//    cumulative_before[rank+1]) in the global index space.
// 2. Each rank's target data should occupy [cumulative_target[rank],
//    cumulative_target[rank+1]) in the same global index space.
// 3. For each pair of (source_rank, dest_rank), we calculate the overlap
//    between what source_rank currently has and what dest_rank needs.
// 4. MPI_Alltoallv transfers exactly these overlapping regions, with
//    displacements ensuring data lands at the correct position in the
//    receiving buffer.
//
// EDGE CASES HANDLED:
// -------------------
// - Single rank (n_procs == 1): Returns immediately with local size, no MPI.
// - Empty total (all ranks have 0 particles): Returns 0 immediately.
// - Imbalanced input (e.g., one rank has all particles): Works correctly;
//   that rank will send portions to all other ranks based on target ranges.
// - Non-divisible totals: First 'remainder' ranks get one extra particle.
int64_t synchronize_global_secondary_bank(
  SharedArray<SourceSite>& shared_secondary_bank)
{
  // Get current size of local bank
  int64_t local_size = shared_secondary_bank.size();

  if (mpi::n_procs == 1) {
    return local_size;
  }

#ifdef OPENMC_MPI
  // Gather all sizes to all ranks
  vector<int64_t> all_sizes(mpi::n_procs);
  MPI_Allgather(&local_size, 1, MPI_INT64_T, all_sizes.data(), 1, MPI_INT64_T,
    mpi::intracomm);

  // Calculate total and check for empty case
  int64_t total = 0;
  for (int64_t size : all_sizes) {
    total += size;
  }

  // If we don't have any items to distribute, return
  if (total == 0) {
    return total;
  }

  int64_t base_count = total / mpi::n_procs;
  int64_t remainder = total % mpi::n_procs;

  // Calculate target size for each rank
  // First 'remainder' ranks get base_count + 1, rest get base_count
  SharedArray<int64_t> target_sizes(mpi::n_procs);
  for (int i = 0; i < mpi::n_procs; ++i) {
    target_sizes[i] = base_count + (i < remainder ? 1 : 0);
  }

  // Calculate send and receive counts in terms of SourceSite objects
  // (not bytes)
  vector<int> send_counts(mpi::n_procs, 0);
  vector<int> recv_counts(mpi::n_procs, 0);
  vector<int> send_displs(mpi::n_procs, 0);
  vector<int> recv_displs(mpi::n_procs, 0);

  // Calculate cumulative positions (starting index for each rank in the
  // global array)
  vector<int64_t> cumulative_before(mpi::n_procs + 1, 0);
  vector<int64_t> cumulative_target(mpi::n_procs + 1, 0);
  for (int i = 0; i < mpi::n_procs; ++i) {
    cumulative_before[i + 1] = cumulative_before[i] + all_sizes[i];
    cumulative_target[i + 1] = cumulative_target[i] + target_sizes[i];
  }

  // Determine send amounts from this rank to others
  int64_t my_start = cumulative_before[mpi::rank];
  int64_t my_end = cumulative_before[mpi::rank + 1];

  for (int dest = 0; dest < mpi::n_procs; ++dest) {
    int64_t dest_start = cumulative_target[dest];
    int64_t dest_end = cumulative_target[dest + 1];

    // Calculate overlap between my current range and destination's
    // target range
    int64_t overlap_start = std::max(my_start, dest_start);
    int64_t overlap_end = std::min(my_end, dest_end);

    if (overlap_start < overlap_end) {
      int64_t count = overlap_end - overlap_start;
      send_counts[dest] =
        static_cast<int>(count); // Count of SourceSite objects
      send_displs[dest] = static_cast<int>(
        overlap_start - my_start); // Displacement in SourceSite objects
    }
  }
  // Determine receive amounts from other ranks
  int64_t my_target_start = cumulative_target[mpi::rank];
  int64_t my_target_end = cumulative_target[mpi::rank + 1];

  for (int src = 0; src < mpi::n_procs; ++src) {
    int64_t src_start = cumulative_before[src];
    int64_t src_end = cumulative_before[src + 1];

    // Calculate overlap between source's current range and my target
    // range
    int64_t overlap_start = std::max(src_start, my_target_start);
    int64_t overlap_end = std::min(src_end, my_target_end);

    if (overlap_start < overlap_end) {
      int64_t count = overlap_end - overlap_start;
      recv_counts[src] = static_cast<int>(count); // Count of SourceSite objects
      recv_displs[src] = static_cast<int>(
        overlap_start - my_target_start); // Displacement in SourceSite objects
    }
  }

  // Prepare receive buffer with target size
  SharedArray<SourceSite> new_bank(target_sizes[mpi::rank]);

  // Perform all-to-all redistribution using the custom MPI type
  MPI_Alltoallv(shared_secondary_bank.data(), send_counts.data(),
    send_displs.data(), mpi::source_site, new_bank.data(), recv_counts.data(),
    recv_displs.data(), mpi::source_site, mpi::intracomm);

  // Replace old bank with redistributed data
  shared_secondary_bank = std::move(new_bank);

  return total;
#endif
}

//==============================================================================
// C API
//==============================================================================

extern "C" int openmc_source_bank(void** ptr, int64_t* n)
{
  if (!ptr || !n) {
    set_errmsg("Received null pointer.");
    return OPENMC_E_INVALID_ARGUMENT;
  }

  if (simulation::source_bank.size() == 0) {
    set_errmsg("Source bank has not been allocated.");
    return OPENMC_E_ALLOCATE;
  } else {
    *ptr = simulation::source_bank.data();
    *n = simulation::source_bank.size();
    return 0;
  }
}

extern "C" int openmc_fission_bank(void** ptr, int64_t* n)
{
  if (!ptr || !n) {
    set_errmsg("Received null pointer.");
    return OPENMC_E_INVALID_ARGUMENT;
  }

  if (simulation::fission_bank.size() == 0) {
    set_errmsg("Fission bank has not been allocated.");
    return OPENMC_E_ALLOCATE;
  } else {
    *ptr = simulation::fission_bank.data();
    *n = simulation::fission_bank.size();
    return 0;
  }
}

} // namespace openmc
