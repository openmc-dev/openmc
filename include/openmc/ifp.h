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
vector<T> _ifp(const T& site_value, const vector<T>& data) {
  vector<T> updated_data;
  size_t ifp_idx = data.size();
  if (ifp_idx < settings::ifp_n_generation) {
    updated_data.resize(ifp_idx + 1);
    for (size_t i = 0; i < ifp_idx; i++) {
      updated_data[i] = data[i];
    }
    updated_data[ifp_idx] = site_value;
  } else if (ifp_idx == settings::ifp_n_generation) {
    updated_data.resize(ifp_idx);
    for (size_t i = 0; i < ifp_idx - 1; i++) {
      updated_data[i] = data[i + 1];
    }
    updated_data[ifp_idx - 1] = site_value;
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
void ifp(Particle& p, SourceSite& site, int64_t idx);

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

void initialize_ifp_pointers(int64_t i, const vector<int>*& delayed_groups_ptr,
  const vector<double>*& lifetimes_ptr);

void add_ifp_data_to_temp(int64_t index_temp,
  vector<vector<int>>& temp_delayed_groups,
  const vector<int>* delayed_groups_ptr, vector<vector<double>>& temp_lifetimes,
  const vector<double>* lifetimes_ptr);

void copy_ifp_data_to_temp_from_fission_banks(int64_t index_temp, int i_bank,
  vector<vector<int>>& temp_delayed_groups,
  vector<vector<double>>& temp_lifetimes);

#ifdef OPENMC_MPI

//! Deserialization info for IFP transfer via MPI
struct DeserializationInfo {
  int64_t index_local;
  int64_t n;
};

//! Broadcast the number of generation from the size of the first element
void broadcast_ifp_n_generation(int& ifp_n_generation,
  vector<vector<int>>& temp_delayed_groups,
  vector<vector<double>>& temp_lifetimes);

//! Send IFP info via MPI
void send_ifp_info(int64_t index_local, int64_t n, int ifp_n_generation,
  int neighbor, vector<MPI_Request>& requests,
  vector<vector<int>> temp_delayed_groups, vector<int> send_delayed_groups,
  vector<vector<double>> temp_lifetimes, vector<double> send_lifetimes);

//! Receive IFP info through MPI
void receive_ifp_data(int64_t index_local, int64_t n, int ifp_n_generation,
  int neighbor, vector<MPI_Request>& requests, vector<int>& recv_delayed_groups,
  vector<double>& recv_lifetimes,
  vector<DeserializationInfo>& deserialization_info);

void copy_ifp_temp_to_source_banks_partial(int64_t index_temp, int n,
  int64_t index_local, const vector<vector<int>>& temp_delayed_groups,
  const vector<vector<double>>& temp_lifetimes);

//! Deserialize IFP information
void deserialize_ifp_info(int ifp_n_generation,
  vector<DeserializationInfo>& deserialization_info,
  vector<int>& recv_delayed_groups, vector<double>& recv_lifetimes);

#endif

//! Copy IFP temporary vector to source banks
void copy_ifp_temp_to_source_banks(vector<vector<int>>& temp_delayed_groups,
  vector<vector<double>>& temp_lifetimes);

// ---------------------------------------------------------
// bank.cpp functions
// ---------------------------------------------------------

void allocate_temporary_vector_ifp(
  vector<vector<int>>& delayed_group_bank_holder,
  vector<int>*& delayed_group_bank,
  vector<vector<double>>& lifetime_bank_holder, vector<double>*& lifetime_bank);

void sort_ifp_banks(int64_t i, int64_t idx, vector<int>*& delayed_group_bank,
  vector<double>*& lifetime_bank);

void copy_ifp_banks_to_fission_banks(
  vector<int>*& delayed_group_bank, vector<double>*& lifetime_bank);

} // namespace openmc

#endif // OPENMC_IFP_H
