#include "openmc/ifp.h"

#include "openmc/bank.h"
#include "openmc/message_passing.h"
#include "openmc/particle.h"
#include "openmc/particle_data.h"
#include "openmc/settings.h"
#include "openmc/simulation.h"
#include "openmc/vector.h"

namespace openmc {

bool is_beta_effective_or_both()
{
  if (settings::ifp_parameter == IFPParameter::BetaEffective ||
      settings::ifp_parameter == IFPParameter::Both) {
    return true;
  }
  return false;
}

bool is_generation_time_or_both()
{
  if (settings::ifp_parameter == IFPParameter::GenerationTime ||
      settings::ifp_parameter == IFPParameter::Both) {
    return true;
  }
  return false;
}

void ifp(const Particle& p, const SourceSite& site, int64_t idx)
{
  if (is_beta_effective_or_both()) {
    const auto& delayed_groups =
      simulation::ifp_source_delayed_group_bank[p.current_work() - 1];
    simulation::ifp_fission_delayed_group_bank[idx] =
      _ifp(site.delayed_group, delayed_groups);
  }
  if (is_generation_time_or_both()) {
    const auto& lifetimes =
      simulation::ifp_source_lifetime_bank[p.current_work() - 1];
    simulation::ifp_fission_lifetime_bank[idx] = _ifp(p.lifetime(), lifetimes);
  }
}

void resize_simulation_ifp_banks()
{
  resize_ifp_data(simulation::ifp_source_delayed_group_bank,
    simulation::ifp_source_lifetime_bank, simulation::work_per_rank);
  resize_ifp_data(simulation::ifp_fission_delayed_group_bank,
    simulation::ifp_fission_lifetime_bank, 3 * simulation::work_per_rank);
}

void copy_ifp_data_from_fission_banks(
  int i_bank, vector<int>& delayed_groups, vector<double>& lifetimes)
{
  if (is_beta_effective_or_both()) {
    delayed_groups = simulation::ifp_fission_delayed_group_bank[i_bank];
  }
  if (is_generation_time_or_both()) {
    lifetimes = simulation::ifp_fission_lifetime_bank[i_bank];
  }
}

#ifdef OPENMC_MPI

void broadcast_ifp_n_generation(int& n_generation,
  const vector<vector<int>>& delayed_groups,
  const vector<vector<double>>& lifetimes)
{
  if (mpi::rank == 0) {
    if (is_beta_effective_or_both()) {
      n_generation = static_cast<int>(delayed_groups[0].size());
    } else {
      n_generation = static_cast<int>(lifetimes[0].size());
    }
  }
  MPI_Bcast(&n_generation, 1, MPI_INT, 0, mpi::intracomm);
}

void send_ifp_info(int64_t idx, int64_t n, int n_generation, int neighbor,
  vector<MPI_Request>& requests, const vector<vector<int>>& delayed_groups,
  vector<int>& send_delayed_groups, const vector<vector<double>>& lifetimes,
  vector<double>& send_lifetimes)
{
  // Copy data in send buffers
  for (int i = idx; i < idx + n; i++) {
    if (is_beta_effective_or_both()) {
      std::copy(delayed_groups[i].begin(), delayed_groups[i].end(),
        send_delayed_groups.begin() + i * n_generation);
    }
    if (is_generation_time_or_both()) {
      std::copy(lifetimes[i].begin(), lifetimes[i].end(),
        send_lifetimes.begin() + i * n_generation);
    }
  }
  // Send delayed groups
  if (is_beta_effective_or_both()) {
    requests.emplace_back();
    MPI_Isend(&send_delayed_groups[n_generation * idx],
      n_generation * static_cast<int>(n), MPI_INT, neighbor, mpi::rank,
      mpi::intracomm, &requests.back());
  }
  // Send lifetimes
  if (is_generation_time_or_both()) {
    requests.emplace_back();
    MPI_Isend(&send_lifetimes[n_generation * idx],
      n_generation * static_cast<int>(n), MPI_DOUBLE, neighbor, mpi::rank,
      mpi::intracomm, &requests.back());
  }
}

void receive_ifp_data(int64_t idx, int64_t n, int n_generation, int neighbor,
  vector<MPI_Request>& requests, vector<int>& delayed_groups,
  vector<double>& lifetimes, vector<DeserializationInfo>& deserialization)
{
  // Receive delayed groups
  if (is_beta_effective_or_both()) {
    requests.emplace_back();
    MPI_Irecv(&delayed_groups[n_generation * idx],
      n_generation * static_cast<int>(n), MPI_INT, neighbor, neighbor,
      mpi::intracomm, &requests.back());
  }
  // Receive lifetimes
  if (is_generation_time_or_both()) {
    requests.emplace_back();
    MPI_Irecv(&lifetimes[n_generation * idx],
      n_generation * static_cast<int>(n), MPI_DOUBLE, neighbor, neighbor,
      mpi::intracomm, &requests.back());
  }
  // Deserialization info to reconstruct data later
  DeserializationInfo info = {idx, n};
  deserialization.push_back(info);
}

void copy_partial_ifp_data_to_source_banks(int64_t idx, int n, int64_t i_bank,
  const vector<vector<int>>& delayed_groups,
  const vector<vector<double>>& lifetimes)
{
  if (is_beta_effective_or_both()) {
    std::copy(&delayed_groups[idx], &delayed_groups[idx + n],
      &simulation::ifp_source_delayed_group_bank[i_bank]);
  }
  if (is_generation_time_or_both()) {
    std::copy(&lifetimes[idx], &lifetimes[idx + n],
      &simulation::ifp_source_lifetime_bank[i_bank]);
  }
}

void deserialize_ifp_info(int n_generation,
  const vector<DeserializationInfo>& deserialization,
  const vector<int>& delayed_groups, const vector<double>& lifetimes)
{
  for (auto info : deserialization) {
    int64_t index_local = info.index_local;
    int64_t n = info.n;

    for (int i = index_local; i < index_local + n; i++) {
      if (is_beta_effective_or_both()) {
        vector<int> delayed_groups_received(
          delayed_groups.begin() + n_generation * i,
          delayed_groups.begin() + n_generation * (i + 1));
        simulation::ifp_source_delayed_group_bank[i] = delayed_groups_received;
      }
      if (is_generation_time_or_both()) {
        vector<double> lifetimes_received(lifetimes.begin() + n_generation * i,
          lifetimes.begin() + n_generation * (i + 1));
        simulation::ifp_source_lifetime_bank[i] = lifetimes_received;
      }
    }
  }
}

#endif

void copy_complete_ifp_data_to_source_banks(
  const vector<vector<int>>& delayed_groups,
  const vector<vector<double>>& lifetimes)
{
  if (is_beta_effective_or_both()) {
    std::copy(delayed_groups.data(),
      delayed_groups.data() + settings::n_particles,
      simulation::ifp_source_delayed_group_bank.begin());
  }
  if (is_generation_time_or_both()) {
    std::copy(lifetimes.data(), lifetimes.data() + settings::n_particles,
      simulation::ifp_source_lifetime_bank.begin());
  }
}

void allocate_temporary_vector_ifp(
  vector<vector<int>>& delayed_groups, vector<vector<double>>& lifetimes)
{
  if (is_beta_effective_or_both()) {
    delayed_groups.resize(simulation::fission_bank.size());
  }
  if (is_generation_time_or_both()) {
    lifetimes.resize(simulation::fission_bank.size());
  }
}

void copy_ifp_data_to_fission_banks(const vector<int>* const delayed_groups_ptr,
  const vector<double>* lifetimes_ptr)
{
  if (is_beta_effective_or_both()) {
    std::copy(delayed_groups_ptr,
      delayed_groups_ptr + simulation::fission_bank.size(),
      simulation::ifp_fission_delayed_group_bank.data());
  }
  if (is_generation_time_or_both()) {
    std::copy(lifetimes_ptr, lifetimes_ptr + simulation::fission_bank.size(),
      simulation::ifp_fission_lifetime_bank.data());
  }
}

} // namespace openmc
