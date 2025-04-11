#ifndef OPENMC_IFP_H
#define OPENMC_IFP_H

#include "openmc/message_passing.h"
#include "openmc/particle.h"
#include "openmc/particle_data.h"
#include "openmc/settings.h"

namespace openmc {

//! Check the value of the IFP parameter for beta effective or both.
//!
//! \return true if "BetaEffective" or "Both", false otherwise.
bool is_beta_effective_or_both();

//! Check the value of the IFP parameter for generation time or both.
//!
//! \return true if "GenerationTime" or "Both", false otherwise.
bool is_generation_time_or_both();

//! Resize IFP vectors
//!
//! \param[in,out] delayed_groups List of delayed group numbers
//! \param[in,out] lifetimes List of lifetimes
//! \param[in] n  Dimension to resize vectors
template<typename T, typename U>
void resize_ifp_data(vector<T>& delayed_groups, vector<U>& lifetimes, int64_t n)
{
  if (is_beta_effective_or_both()) {
    delayed_groups.resize(n);
  }
  if (is_generation_time_or_both()) {
    lifetimes.resize(n);
  }
}

//! Update a list of values by adding a new value if the size
//! of the list can accomodate the new value or by shifting all
//! values to the left (removing the first value of the list
//! and adding the new value at the end of the list).
//!
//! \param[in] value Value to add to the list
//! \param[in] data Initial version of the list
//! \return Updated list
template<typename T>
vector<T> _ifp(const T& value, const vector<T>& data)
{
  vector<T> updated;
  size_t source_idx = data.size();

  if (source_idx < settings::ifp_n_generation) {
    updated.resize(source_idx + 1);
    for (size_t i = 0; i < source_idx; i++) {
      updated[i] = data[i];
    }
    updated[source_idx] = value;
  } else if (source_idx == settings::ifp_n_generation) {
    updated.resize(source_idx);
    for (size_t i = 0; i < source_idx - 1; i++) {
      updated[i] = data[i + 1];
    }
    updated[source_idx - 1] = value;
  }
  return updated;
}

//! \brief Iterated Fission Probability (IFP) method.
//!
//! Add the IFP information in the IFP banks using the same index
//! as the one used to append the fission site to the fission bank.
//! Multithreading protection is guaranteed by the index returned by the
//! thread_safe_append call in physics.cpp.
//!
//! Needs to be done after the delayed group is found.
//!
//! \param[in] p Particle
//! \param[in] site Fission site
//! \param[in] idx Bank index from the thread_safe_append call in physics.cpp
void ifp(const Particle& p, const SourceSite& site, int64_t idx);

//! Resize the IFP banks used in the simulation
void resize_simulation_ifp_banks();

//! Retrieve IFP data from the IFP fission banks.
//!
//! \param[in] i_bank Index in the fission banks
//! \param[in,out] delayed_groups Delayed group numbers
//! \param[in,out] lifetimes Lifetimes lists
void copy_ifp_data_from_fission_banks(
  int i_bank, vector<int>& delayed_groups, vector<double>& lifetimes);

#ifdef OPENMC_MPI

//! Deserialization information for transfer of IFP data using MPI
struct DeserializationInfo {
  int64_t index_local; //!< local index
  int64_t n;           //!< number of sites sent
};

//! Broadcast the number of generation determined by the size of the first
//! element on the first processor.
//!
//! \param[in] n_generation Number of generations
//! \param[in] delayed_groups List of delayed group numbers lists
//! \param[in] lifetimes List of lifetimes lists
void broadcast_ifp_n_generation(int& n_generation,
  const vector<vector<int>>& delayed_groups,
  const vector<vector<double>>& lifetimes);

//! Send IFP data using MPI.
//!
//! \param[in] idx Index of the first site
//! \param[in] n Number of sites to send
//! \param[in] n_generation Number of generations
//! \param[in] neighbor Index of the neighboring processor
//! \param[in] requests MPI requests
//! \param[in] delayed_groups List of delayed group numbers lists
//! \param[out] send_delayed_groups Delayed group numbers buffer
//! \param[in] lifetimes List of lifetimes lists
//! \param[out] send_lifetimes Lifetimes buffer
void send_ifp_info(int64_t idx, int64_t n, int n_generation, int neighbor,
  vector<MPI_Request>& requests, const vector<vector<int>>& delayed_groups,
  vector<int>& send_delayed_groups, const vector<vector<double>>& lifetimes,
  vector<double>& send_lifetimes);

//! Receive IFP data using MPI.
//!
//! \param[in] idx Index of the first site
//! \param[in] n Number of sites to receive
//! \param[in] n_generation Number of generations
//! \param[in] neighbor Index of the neighboring processor
//! \param[in] requests MPI requests
//! \param[in] delayed_groups List of delayed group numbers
//! \param[in] lifetimes List of lifetimes
//! \param[out] deserialization Information to deserialize the received data
void receive_ifp_data(int64_t idx, int64_t n, int n_generation, int neighbor,
  vector<MPI_Request>& requests, vector<int>& delayed_groups,
  vector<double>& lifetimes, vector<DeserializationInfo>& deserialization);

//! Copy partial IFP data from local lists to source banks.
//!
//! \param[in] idx Index of the first site
//! \param[in] n Number of sites to copy
//! \param[in] i_bank Index in the IFP source banks
//! \param[in] delayed_groups List of delayed group numbers lists
//! \param[in] lifetimes List of lifetimes lists
void copy_partial_ifp_data_to_source_banks(int64_t idx, int n, int64_t i_bank,
  const vector<vector<int>>& delayed_groups,
  const vector<vector<double>>& lifetimes);

//! Deserialize IFP information received using MPI and store it in
//! the IFP source banks.
//!
//! \param[in] n_generation Number of generations
//! \param[out] deserialization Information to deserialize the received data
//! \param[in] delayed_groups List of delayed group numbers
//! \param[in] lifetimes List of lifetimes
void deserialize_ifp_info(int n_generation,
  const vector<DeserializationInfo>& deserialization,
  const vector<int>& delayed_groups, const vector<double>& lifetimes);

#endif

//! Copy IFP temporary vectors to source banks.
//!
//! \param[in] delayed_groups List of delayed group numbers lists
//! \param[in] lifetimes List of lifetimes lists
void copy_complete_ifp_data_to_source_banks(
  const vector<vector<int>>& delayed_groups,
  const vector<vector<double>>& lifetimes);

//! Allocate temporary vectors for IFP data.
//!
//! \param[in,out] delayed_groups List of delayed group numbers lists
//! \param[in,out] lifetimes List of delayed group numbers lists
void allocate_temporary_vector_ifp(
  vector<vector<int>>& delayed_groups, vector<vector<double>>& lifetimes);

//! Copy local IFP data to IFP fission banks.
//!
//! \param[in] delayed_groups_ptr Pointer to delayed group numbers
//! \param[in] lifetimes_ptr Pointer to lifetimes
void copy_ifp_data_to_fission_banks(
  const vector<int>* delayed_groups_ptr, const vector<double>* lifetimes_ptr);

} // namespace openmc

#endif // OPENMC_IFP_H
