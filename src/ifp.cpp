#include "openmc/ifp.h"

#include "openmc/bank.h"
#include "openmc/message_passing.h"
#include "openmc/particle.h"
#include "openmc/particle_data.h"
#include "openmc/settings.h"
#include "openmc/simulation.h"
#include "openmc/vector.h"

namespace openmc {

// ---------------------------------------------------------
// Common helpers
// ---------------------------------------------------------

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

// ---------------------------------------------------------
// physics.cpp functions
// ---------------------------------------------------------

void ifp(Particle& p, SourceSite& site, int64_t idx)
{
  // Beta effective
  if (is_beta_effective_or_both()) {
    const auto& ifp_delayed_groups =
      simulation::ifp_source_delayed_group_bank[p.current_work() - 1];
    vector<int> updated_ifp_delayed_groups =
      _ifp(site.delayed_group, ifp_delayed_groups);
    simulation::ifp_fission_delayed_group_bank[idx] =
      updated_ifp_delayed_groups;
  }
  // Generation time
  if (is_generation_time_or_both()) {
    const auto& ifp_lifetimes =
      simulation::ifp_source_lifetime_bank[p.current_work() - 1];
    vector<double> updated_ifp_lifetimes = _ifp(p.lifetime(), ifp_lifetimes);
    simulation::ifp_fission_lifetime_bank[idx] = updated_ifp_lifetimes;
  }
}

// ---------------------------------------------------------
// simulation.cpp functions
// ---------------------------------------------------------

void resize_simulation_ifp_banks() {
  if (is_beta_effective_or_both()) {
    simulation::ifp_source_delayed_group_bank.resize(
      simulation::work_per_rank);
    simulation::ifp_fission_delayed_group_bank.resize(
      3 * simulation::work_per_rank);
  }
  if (is_generation_time_or_both()) {
    simulation::ifp_source_lifetime_bank.resize(simulation::work_per_rank);
    simulation::ifp_fission_lifetime_bank.resize(
      3 * simulation::work_per_rank);
  }
}

// ---------------------------------------------------------
// eigenvalue.cpp functions
// ---------------------------------------------------------

void initialize_ifp_pointers(int64_t i, const vector<int>*& delayed_groups_ptr,
  const vector<double>*& lifetimes_ptr)
{
  if (is_beta_effective_or_both()) {
    delayed_groups_ptr = &simulation::ifp_fission_delayed_group_bank[i];
  }
  if (is_generation_time_or_both()) {
    lifetimes_ptr = &simulation::ifp_fission_lifetime_bank[i];
  }
}

void add_ifp_data_to_temp(int64_t index_temp,
  vector<vector<int>>& temp_delayed_groups,
  const vector<int>* delayed_groups_ptr, vector<vector<double>>& temp_lifetimes,
  const vector<double>* lifetimes_ptr)
{
  if (is_beta_effective_or_both()) {
    temp_delayed_groups[index_temp] = *delayed_groups_ptr;
  }
  if (is_generation_time_or_both()) {
    temp_lifetimes[index_temp] = *lifetimes_ptr;
  }
}

void copy_ifp_data_to_temp_from_fission_banks(int64_t index_temp, int i_bank,
  vector<vector<int>>& temp_delayed_groups,
  vector<vector<double>>& temp_lifetimes)
{
  if (is_beta_effective_or_both()) {
    temp_delayed_groups[index_temp] =
      simulation::ifp_fission_delayed_group_bank[i_bank];
  }
  if (is_generation_time_or_both()) {
    temp_lifetimes[index_temp] = simulation::ifp_fission_lifetime_bank[i_bank];
  }
}

#ifdef OPENMC_MPI

void broadcast_ifp_n_generation(int& ifp_n_generation,
  vector<vector<int>>& temp_delayed_groups,
  vector<vector<double>>& temp_lifetimes)
{
  if (mpi::rank == 0) {
    if (is_beta_effective_or_both()) {
      ifp_n_generation = static_cast<int>(temp_delayed_groups[0].size());
    } else {
      ifp_n_generation = static_cast<int>(temp_lifetimes[0].size());
    }
  }
  MPI_Bcast(&ifp_n_generation, 1, MPI_INT, 0, mpi::intracomm);
}

void send_ifp_info(int64_t index_local, int64_t n, int ifp_n_generation,
  int neighbor, vector<MPI_Request>& requests,
  vector<vector<int>> temp_delayed_groups, vector<int> send_delayed_groups,
  vector<vector<double>> temp_lifetimes, vector<double> send_lifetimes)
{
  for (int i = index_local; i < index_local + n; i++) {
    if (is_beta_effective_or_both()) {
      std::copy(temp_delayed_groups[i].begin(), temp_delayed_groups[i].end(),
        send_delayed_groups.begin() + i * ifp_n_generation);
    }
    if (is_generation_time_or_both()) {
      std::copy(temp_lifetimes[i].begin(), temp_lifetimes[i].end(),
        send_lifetimes.begin() + i * ifp_n_generation);
    }
  }

  // Send delayed groups
  if (is_beta_effective_or_both()) {
    requests.emplace_back();
    MPI_Isend(&send_delayed_groups[ifp_n_generation * index_local],
      ifp_n_generation * static_cast<int>(n), MPI_INT, neighbor, mpi::rank,
      mpi::intracomm, &requests.back());
  }

  // Send lifetimes
  if (is_generation_time_or_both()) {
    requests.emplace_back();
    MPI_Isend(&send_lifetimes[ifp_n_generation * index_local],
      ifp_n_generation * static_cast<int>(n), MPI_DOUBLE, neighbor, mpi::rank,
      mpi::intracomm, &requests.back());
  }
}

void receive_ifp_data(int64_t index_local, int64_t n, int ifp_n_generation,
  int neighbor, vector<MPI_Request>& requests, vector<int>& recv_delayed_groups,
  vector<double>& recv_lifetimes,
  vector<DeserializationInfo>& deserialization_info)
{
  // Receive delayed groups
  if (is_beta_effective_or_both()) {
    requests.emplace_back();
    MPI_Irecv(&recv_delayed_groups[ifp_n_generation * index_local],
      ifp_n_generation * static_cast<int>(n), MPI_INT, neighbor, neighbor,
      mpi::intracomm, &requests.back());
  }

  // Receive lifetimes
  if (is_generation_time_or_both()) {
    requests.emplace_back();
    MPI_Irecv(&recv_lifetimes[ifp_n_generation * index_local],
      ifp_n_generation * static_cast<int>(n), MPI_DOUBLE, neighbor, neighbor,
      mpi::intracomm, &requests.back());
  }

  // Deserialization info to reconstruct data later
  DeserializationInfo info = {index_local, n};
  deserialization_info.push_back(info);
}

void copy_ifp_temp_to_source_banks_partial(int64_t index_temp, int n,
  int64_t index_local, const vector<vector<int>>& temp_delayed_groups,
  const vector<vector<double>>& temp_lifetimes)
{
  if (is_beta_effective_or_both()) {
    std::copy(&temp_delayed_groups[index_temp],
      &temp_delayed_groups[index_temp + n],
      &simulation::ifp_source_delayed_group_bank[index_local]);
  }
  if (is_generation_time_or_both()) {
    std::copy(&temp_lifetimes[index_temp], &temp_lifetimes[index_temp + n],
      &simulation::ifp_source_lifetime_bank[index_local]);
  }
}

void deserialize_ifp_info(int ifp_n_generation,
  vector<DeserializationInfo>& deserialization_info,
  vector<int>& recv_delayed_groups, vector<double>& recv_lifetimes)
{
  int64_t n;
  for (auto info : deserialization_info) {
    int64_t index_local = info.index_local;
    n = info.n;

    for (int i = index_local; i < index_local + n; i++) {
      if (is_beta_effective_or_both()) {
        vector<int> delayed_groups_received(
          recv_delayed_groups.begin() + ifp_n_generation * i,
          recv_delayed_groups.begin() + ifp_n_generation * (i + 1));
        simulation::ifp_source_delayed_group_bank[i] = delayed_groups_received;
      }
      if (is_generation_time_or_both()) {
        vector<double> lifetimes_received(
          recv_lifetimes.begin() + ifp_n_generation * i,
          recv_lifetimes.begin() + ifp_n_generation * (i + 1));
        simulation::ifp_source_lifetime_bank[i] = lifetimes_received;
      }
    }
  }
}

#endif

void copy_ifp_temp_to_source_banks(vector<vector<int>>& temp_delayed_groups,
  vector<vector<double>>& temp_lifetimes)
{
  if (is_beta_effective_or_both()) {
    std::copy(temp_delayed_groups.data(),
      temp_delayed_groups.data() + settings::n_particles,
      simulation::ifp_source_delayed_group_bank.begin());
  }
  if (is_generation_time_or_both()) {
    std::copy(temp_lifetimes.data(),
      temp_lifetimes.data() + settings::n_particles,
      simulation::ifp_source_lifetime_bank.begin());
  }
}

// ---------------------------------------------------------
// bank.cpp functions
// ---------------------------------------------------------

void allocate_temporary_vector_ifp(
  vector<vector<int>>& delayed_group_bank_holder,
  vector<int>*& delayed_group_bank,
  vector<vector<double>>& lifetime_bank_holder, vector<double>*& lifetime_bank)
{
  if (is_beta_effective_or_both()) {
    delayed_group_bank_holder.resize(simulation::fission_bank.size());
    delayed_group_bank = delayed_group_bank_holder.data();
  }
  if (is_generation_time_or_both()) {
    lifetime_bank_holder.resize(simulation::fission_bank.size());
    lifetime_bank = lifetime_bank_holder.data();
  }
}

void sort_ifp_banks(int64_t i, int64_t idx, vector<int>*& delayed_group_bank,
  vector<double>*& lifetime_bank)
{
  if (is_beta_effective_or_both()) {
    const auto& delayed_groups = simulation::ifp_fission_delayed_group_bank[i];
    delayed_group_bank[idx] = delayed_groups;
  }
  if (is_generation_time_or_both()) {
    const auto& lifetimes = simulation::ifp_fission_lifetime_bank[i];
    lifetime_bank[idx] = lifetimes;
  }
}

void copy_ifp_banks_to_fission_banks(
  vector<int>*& delayed_group_bank, vector<double>*& lifetime_bank)
{
  if (is_beta_effective_or_both()) {
    std::copy(delayed_group_bank,
      delayed_group_bank + simulation::fission_bank.size(),
      simulation::ifp_fission_delayed_group_bank.data());
  }
  if (is_generation_time_or_both()) {
    std::copy(lifetime_bank, lifetime_bank + simulation::fission_bank.size(),
      simulation::ifp_fission_lifetime_bank.data());
  }
}

} // namespace openmc
