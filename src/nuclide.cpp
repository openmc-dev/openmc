#include "openmc/nuclide.h"

namespace openmc {

// Order corresponds to ParticleType enum
std::array<double, 2> energy_min {0.0, 0.0};
std::array<double, 2> energy_max {INFTY, INFTY};

extern "C" void
set_particle_energy_bounds(int particle, double E_min, double E_max)
{
  energy_min[particle - 1] = E_min;
  energy_max[particle - 1] = E_max;
}

} // namespace openmc
