#ifndef OPENMC_IFP_H
#define OPENMC_IFP_H

#include "openmc/message_passing.h"
#include "openmc/particle.h"
#include "openmc/particle_data.h"
#include "openmc/settings.h"

namespace openmc {

// ---------------------------------------------------------
// Common helpers
// ---------------------------------------------------------

bool is_beta_effective_or_both();

bool is_generation_time_or_both();

// ---------------------------------------------------------
// physics.cpp functions
// ---------------------------------------------------------

template<typename T>
vector<T> _ifp(const T& value, const vector<T>& data)
{
  vector<T> updated_data;
  size_t ifp_idx = data.size();
  if (ifp_idx < settings::ifp_n_generation) {
    updated_data.resize(ifp_idx + 1);
    for (size_t i = 0; i < ifp_idx; i++) {
      updated_data[i] = data[i];
    }
    updated_data[ifp_idx] = value;
  } else if (ifp_idx == settings::ifp_n_generation) {
    updated_data.resize(ifp_idx);
    for (size_t i = 0; i < ifp_idx - 1; i++) {
      updated_data[i] = data[i + 1];
    }
    updated_data[ifp_idx - 1] = value;
  }
  return updated_data;
}

//! Iterated Fission Probability (IFP) method
//!
//! Needs to be done after the delayed group is found.
//!
//! Add the IFP information in the IFP banks using the same index
//! as the one used to append the fission site to the fission bank.
//! Multithreading protection is guaranteed by the index returned by the
//! thread_safe_append call in physics.cpp.
void ifp(const Particle& p, const SourceSite& site, int64_t idx);

// ---------------------------------------------------------
// simulation.cpp functions
// ---------------------------------------------------------

//! Resize the IFP banks used in the simulation
void resize_simulation_ifp_banks();

// ---------------------------------------------------------
// eigenvalue.cpp functions
// ---------------------------------------------------------

//! Resize IFP vectors
template<typename T, typename U>
void resize_ifp_data(
  vector<T>& delayed_groups, vector<U>& lifetimes, int64_t n_dg, int64_t n_l)
{
  if (is_beta_effective_or_both()) {
    delayed_groups.resize(n_dg);
  }
  if (is_generation_time_or_both()) {
    lifetimes.resize(n_l);
  }
}

void initialize_ifp_pointers(int64_t i_bank,
  const vector<int>*& delayed_groups_ptr, const vector<double>*& lifetimes_ptr);

void add_ifp_data(int64_t idx, vector<vector<int>>& delayed_groups,
  const vector<int>* delayed_groups_ptr, vector<vector<double>>& lifetimes,
  const vector<double>* lifetimes_ptr);

void retrieve_ifp_data_from_fission_banks(int64_t idx, int i_bank,
  vector<vector<int>>& delayed_groups, vector<vector<double>>& lifetimes);

#ifdef OPENMC_MPI

//! Deserialization info for IFP transfer via MPI
struct DeserializationInfo {
  int64_t index_local;
  int64_t n;
};

//! Broadcast the number of generation from the size of the first element
void broadcast_ifp_n_generation(int& n_generation,
  const vector<vector<int>>& delayed_groups,
  const vector<vector<double>>& lifetimes);

//! Send IFP info via MPI
void send_ifp_info(int64_t idx, int64_t n, int n_generation, int neighbor,
  vector<MPI_Request>& requests, const vector<vector<int>> delayed_groups,
  vector<int> send_delayed_groups, const vector<vector<double>> lifetimes,
  vector<double> send_lifetimes);

//! Receive IFP info through MPI
void receive_ifp_data(int64_t idx, int64_t n, int n_generation, int neighbor,
  vector<MPI_Request>& requests, vector<int>& delayed_groups,
  vector<double>& lifetimes, vector<DeserializationInfo>& deserialization);

void copy_partial_ifp_data_to_source_banks(int64_t idx, int n, int64_t i_bank,
  const vector<vector<int>>& delayed_groups,
  const vector<vector<double>>& lifetimes);

//! Deserialize IFP information
void deserialize_ifp_info(int n_generation,
  const vector<DeserializationInfo>& deserialization,
  const vector<int>& delayed_groups, const vector<double>& lifetimes);

#endif

//! Copy IFP temporary vector to source banks
void copy_complete_ifp_data_to_source_banks(
  const vector<vector<int>>& delayed_groups,
  const vector<vector<double>>& lifetimes);

// ---------------------------------------------------------
// bank.cpp functions
// ---------------------------------------------------------

void allocate_temporary_vector_ifp(
  vector<vector<int>>& delayed_group_bank_holder,
  vector<int>*& delayed_group_bank,
  vector<vector<double>>& lifetime_bank_holder, vector<double>*& lifetime_bank);

void sort_ifp_data_from_fission_banks(int64_t i_bank, int64_t idx,
  vector<int>*& delayed_group_bank, vector<double>*& lifetime_bank);

void copy_ifp_data_to_fission_banks(
  vector<int>* const& delayed_group_bank, vector<double>* const& lifetime_bank);

} // namespace openmc

#endif // OPENMC_IFP_H
