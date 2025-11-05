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

} // namespace simulation

//==============================================================================
// Non-member functions
//==============================================================================

void free_memory_bank()
{
  simulation::source_bank.clear();
  simulation::surf_source_bank.clear();
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

// Performs an O(n) sort on the fission bank, by leveraging
// the parent_id and progeny_id fields of banked particles. See the following
// paper for more details:
// "Reproducibility and Monte Carlo Eigenvalue Calculations," F.B. Brown and
// T.M. Sutton, 1992 ANS Annual Meeting, Transactions of the American Nuclear
// Society, Volume 65, Page 235.
void sort_fission_bank()
{
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
  if (simulation::fission_bank.size() >
      simulation::fission_bank.capacity() / 2) {
    sorted_bank_holder.resize(simulation::fission_bank.size());
    sorted_bank = sorted_bank_holder.data();
  } else { // otherwise, point sorted_bank to unused portion of the fission bank
    sorted_bank = &simulation::fission_bank[simulation::fission_bank.size()];
  }

  if (settings::ifp_on) {
    allocate_temporary_vector_ifp(
      sorted_ifp_delayed_group_bank, sorted_ifp_lifetime_bank);
  }

  // Use parent and progeny indices to sort fission bank
  for (int64_t i = 0; i < simulation::fission_bank.size(); i++) {
    const auto& site = simulation::fission_bank[i];
    int64_t offset = site.parent_id - 1 - simulation::work_index[mpi::rank];
    int64_t idx = simulation::progeny_per_particle[offset] + site.progeny_id;
    if (idx >= simulation::fission_bank.size()) {
      fatal_error("Mismatch detected between sum of all particle progeny and "
                  "shared fission bank size.");
    }
    sorted_bank[idx] = site;
    if (settings::ifp_on) {
      copy_ifp_data_from_fission_banks(
        i, sorted_ifp_delayed_group_bank[idx], sorted_ifp_lifetime_bank[idx]);
    }
  }

  // Copy sorted bank into the fission bank
  std::copy(sorted_bank, sorted_bank + simulation::fission_bank.size(),
    simulation::fission_bank.data());
  if (settings::ifp_on) {
    copy_ifp_data_to_fission_banks(
      sorted_ifp_delayed_group_bank.data(), sorted_ifp_lifetime_bank.data());
  }
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
