#ifndef OPENMC_IFP_H
#define OPENMC_IFP_H

#include "openmc/message_passing.h"
#include "openmc/particle.h"
#include "openmc/particle_data.h"
#include "openmc/settings.h"

#include <algorithm> // for copy

namespace openmc {

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

//==============================================================================
//! One stream of Iterated Fission Probability history data.
//
//! IFP tracks two quantities with identical bookkeeping: the delayed group
//! number of the ancestor neutron (an int) and its lifetime (a double). Each is
//! a per-particle list of the last ifp_n_generation values, held in a source
//! bank and a fission bank, and each is maintained only if a tally asked for
//! it.
//!
//! Bundling the flag with the two banks it governs means every operation is a
//! no-op when the stream is off, so callers do not branch, and the banks cannot
//! be resized independently of the flag.
//==============================================================================

template<typename T>
class IFPStream {
public:
  //! Whether this stream is being maintained
  bool enabled() const { return enabled_; }

  //! Begin maintaining this stream, when a tally requests a dependent score
  void enable() { enabled_ = true; }

  //! Clear the stream, including the flag. Derived from the tallies present in
  //! the model, so it must not persist into the next model in this process.
  void reset()
  {
    enabled_ = false;
    source_bank_.clear();
    fission_bank_.clear();
  }

  vector<vector<T>>& source_bank() { return source_bank_; }
  const vector<vector<T>>& source_bank() const { return source_bank_; }
  vector<vector<T>>& fission_bank() { return fission_bank_; }

  //! Resize both banks
  void resize_banks(int64_t n_source, int64_t n_fission)
  {
    if (!enabled_)
      return;
    source_bank_.resize(n_source);
    fission_bank_.resize(n_fission);
  }

  //! Append a value to the ancestor history of a new fission site
  void store(int64_t i_source, int64_t i_fission, const T& value)
  {
    if (!enabled_)
      return;
    fission_bank_[i_fission] = _ifp(value, source_bank_[i_source]);
  }

  //! Resize a caller-owned buffer for this stream, if it is enabled.
  //
  //! Templated on the container because callers pass both the per-particle
  //! banks (`vector<vector<T>>`) and the flat MPI serialization buffers
  //! (`vector<T>`).
  template<typename V>
  void resize_temp(V& temp, int64_t n) const
  {
    if (enabled_)
      temp.resize(n);
  }

  //! Copy one entry of the fission bank into a caller-owned buffer
  //
  //! \param[in] i_fission Index in this stream's fission bank
  //! \param[in] i_temp Index in the destination buffer
  //! \param[out] out Destination buffer
  void copy_from_fission(
    int64_t i_fission, int64_t i_temp, vector<vector<T>>& out) const
  {
    if (enabled_)
      out[i_temp] = fission_bank_[i_fission];
  }

  //! Copy a temporary buffer into the fission bank
  void copy_temp_to_fission(const vector<T>* temp, size_t n)
  {
    if (enabled_)
      std::copy(temp, temp + n, fission_bank_.data());
  }

  //! Copy a run of a temporary buffer into the source bank
  void copy_temp_to_source(
    int64_t i_temp, int64_t n, int64_t i_source, const vector<vector<T>>& temp)
  {
    if (enabled_)
      std::copy(&temp[i_temp], &temp[i_temp + n], &source_bank_[i_source]);
  }

  //! Copy a whole temporary buffer into the source bank
  void copy_all_temp_to_source(const vector<vector<T>>& temp, int64_t n)
  {
    if (enabled_)
      std::copy(temp.data(), temp.data() + n, source_bank_.begin());
  }

private:
  bool enabled_ {false};
  vector<vector<T>> source_bank_;
  vector<vector<T>> fission_bank_;
};

namespace simulation {

extern IFPStream<int> ifp_delayed_group; //!< Delayed group numbers
extern IFPStream<double> ifp_lifetime;   //!< Neutron lifetimes

} // namespace simulation

//! Whether Iterated Fission Probability is in use at all
inline bool ifp_on()
{
  return simulation::ifp_delayed_group.enabled() ||
         simulation::ifp_lifetime.enabled();
}

//! \brief Iterated Fission Probability (IFP) method.
//!
//! Add the IFP information in the IFP banks using the same index
//! as the one used to append the fission site to the fission bank.
//! Multithreading protection is guaranteed by the index returned by the
//! thread_safe_append call in physics.cpp.
//!
//! \param[in] p Particle
//! \param[in] idx Bank index from the thread_safe_append call in physics.cpp
void ifp(const Particle& p, int64_t idx);

//! Resize the IFP banks used in the simulation
void resize_simulation_ifp_banks();

//! Clear both streams, including their flags
void reset_ifp_streams();

#ifdef OPENMC_MPI

//! Deserialization information for transfer of IFP data using MPI
struct DeserializationInfo {
  int64_t index_local; //!< local index
  int64_t n;           //!< number of sites sent
};

//! Broadcast the number of generations, determined from the first element on
//! the first processor of whichever stream is enabled.
//!
//! \param[in,out] n_generation Number of generations
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
//! \param[in] data List of data lists
//! \param[out] send_data data buffer
template<typename T>
void send_ifp_info(int64_t idx, int64_t n, int n_generation, int neighbor,
  vector<MPI_Request>& requests, const vector<vector<T>>& data,
  vector<T>& send_data)
{
  // Copy data in buffer
  for (int i = idx; i < idx + n; i++) {
    std::copy(
      data[i].begin(), data[i].end(), send_data.begin() + i * n_generation);
  }

  // Send data
  requests.emplace_back();
  MPI_Datatype datatype = mpi::MPITypeMap<T>::mpi_type;
  MPI_Isend(&send_data[n_generation * idx], n_generation * static_cast<int>(n),
    datatype, neighbor, mpi::rank, mpi::intracomm, &requests.back());
}

//! Receive IFP data using MPI.
//!
//! \param[in] idx Index of the first site
//! \param[in] n Number of sites to receive
//! \param[in] n_generation Number of generations
//! \param[in] neighbor Index of the neighboring processor
//! \param[in] requests MPI requests
//! \param[in] data data buffer
//! \param[out] deserialization Information to deserialize the received data
template<typename T>
void receive_ifp_data(int64_t idx, int64_t n, int n_generation, int neighbor,
  vector<MPI_Request>& requests, vector<T>& data,
  vector<DeserializationInfo>& deserialization)
{
  requests.emplace_back();
  MPI_Datatype datatype = mpi::MPITypeMap<T>::mpi_type;
  MPI_Irecv(&data[n_generation * idx], n_generation * static_cast<int>(n),
    datatype, neighbor, neighbor, mpi::intracomm, &requests.back());

  // Deserialization info to reconstruct data later
  DeserializationInfo info = {idx, n};
  deserialization.push_back(info);
}

//! Deserialize IFP information received using MPI into a stream's source bank.
//!
//! \param[in] n_generation Number of generations
//! \param[in] data data to deserialize
//! \param[in,out] stream Stream whose source bank receives the data
//! \param[in] deserialization Information to deserialize the received data
template<typename T>
void deserialize_ifp_info(int n_generation, const vector<T>& data,
  IFPStream<T>& stream, const vector<DeserializationInfo>& deserialization)
{
  if (!stream.enabled())
    return;

  auto& bank = stream.source_bank();
  for (auto info : deserialization) {
    int64_t index_local = info.index_local;
    int64_t n = info.n;

    for (int i = index_local; i < index_local + n; i++) {
      vector<T> data_received(
        data.begin() + n_generation * i, data.begin() + n_generation * (i + 1));
      bank[i] = data_received;
    }
  }
}

#endif

} // namespace openmc

#endif // OPENMC_IFP_H
