#include "openmc/ifp.h"

#include "openmc/bank.h"
#include "openmc/particle.h"
#include "openmc/particle_data.h"
#include "openmc/settings.h"
#include "openmc/simulation.h"
#include "openmc/vector.h"

namespace openmc {

void ifp(Particle& p, SourceSite& site, int64_t idx)
{
  // Iterated Fission Probability (IFP) method

  // Needs to be done after the delayed group is found.

  // Add the IFP information in the IFP banks using the same index
  // as the one used to append the fission site to the fission bank.
  // Multithreading protection is guaranteed by the index returned by the
  // thread_safe_append call in physics.cpp.

  // Beta effective
  if (settings::ifp_parameter == IFPParameter::BetaEffective ||
      settings::ifp_parameter == IFPParameter::Both) {

    vector<int> updated_ifp_delayed_groups;

    const auto& ifp_delayed_groups =
      simulation::ifp_source_delayed_group_bank[p.current_work() - 1];
    size_t ifp_idx = ifp_delayed_groups.size();

    if (ifp_idx < settings::ifp_n_generation) {
      updated_ifp_delayed_groups.resize(ifp_idx + 1);
      for (size_t i = 0; i < ifp_idx; i++) {
        updated_ifp_delayed_groups[i] = ifp_delayed_groups[i];
      }
      updated_ifp_delayed_groups[ifp_idx] = site.delayed_group;
    } else if (ifp_idx == settings::ifp_n_generation) {
      updated_ifp_delayed_groups.resize(ifp_idx);
      for (size_t i = 0; i < ifp_idx - 1; i++) {
        updated_ifp_delayed_groups[i] = ifp_delayed_groups[i + 1];
      }
      updated_ifp_delayed_groups[ifp_idx - 1] = site.delayed_group;
    }

    simulation::ifp_fission_delayed_group_bank[idx] =
      updated_ifp_delayed_groups;
  }

  // Generation time
  if (settings::ifp_parameter == IFPParameter::GenerationTime ||
      settings::ifp_parameter == IFPParameter::Both) {

    vector<double> updated_ifp_lifetimes;

    const auto& ifp_lifetimes =
      simulation::ifp_source_lifetime_bank[p.current_work() - 1];
    size_t ifp_idx = ifp_lifetimes.size();

    if (ifp_idx < settings::ifp_n_generation) {
      updated_ifp_lifetimes.resize(ifp_idx + 1);
      for (size_t i = 0; i < ifp_idx; i++) {
        updated_ifp_lifetimes[i] = ifp_lifetimes[i];
      }
      updated_ifp_lifetimes[ifp_idx] = p.lifetime();
    } else if (ifp_idx == settings::ifp_n_generation) {
      updated_ifp_lifetimes.resize(ifp_idx);
      for (size_t i = 0; i < ifp_idx - 1; i++) {
        updated_ifp_lifetimes[i] = ifp_lifetimes[i + 1];
      }
      updated_ifp_lifetimes[ifp_idx - 1] = p.lifetime();
    }

    simulation::ifp_fission_lifetime_bank[idx] = updated_ifp_lifetimes;
  }
}

} // namespace openmc
