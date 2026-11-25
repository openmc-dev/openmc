#include "openmc/ifp.h"

#include "openmc/bank.h"
#include "openmc/message_passing.h"
#include "openmc/particle.h"
#include "openmc/particle_data.h"
#include "openmc/settings.h"
#include "openmc/simulation.h"
#include "openmc/vector.h"

namespace openmc {

namespace simulation {

IFPStream<int> ifp_delayed_group;
IFPStream<double> ifp_lifetime;

} // namespace simulation

void ifp(const Particle& p, int64_t idx)
{
  simulation::ifp_delayed_group.store(p.current_work(), idx, p.delayed_group());
  simulation::ifp_lifetime.store(p.current_work(), idx, p.lifetime());
}

void resize_simulation_ifp_banks()
{
  int64_t n = simulation::work_per_rank;
  simulation::ifp_delayed_group.resize_banks(n, 3 * n);
  simulation::ifp_lifetime.resize_banks(n, 3 * n);
}

void reset_ifp_streams()
{
  simulation::ifp_delayed_group.reset();
  simulation::ifp_lifetime.reset();
}

#ifdef OPENMC_MPI
void broadcast_ifp_n_generation(int& n_generation,
  const vector<vector<int>>& delayed_groups,
  const vector<vector<double>>& lifetimes)
{
  if (mpi::rank == 0) {
    n_generation = simulation::ifp_delayed_group.enabled()
                     ? static_cast<int>(delayed_groups[0].size())
                     : static_cast<int>(lifetimes[0].size());
  }
  MPI_Bcast(&n_generation, 1, MPI_INT, 0, mpi::intracomm);
}
#endif

} // namespace openmc
