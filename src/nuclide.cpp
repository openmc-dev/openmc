#include "openmc/nuclide.h"

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

namespace data {

std::array<double, 2> energy_min {0.0, 0.0};
std::array<double, 2> energy_max {INFTY, INFTY};

} // namespace data

//==============================================================================
// Fortran compatibility functions
//==============================================================================

extern "C" void
set_particle_energy_bounds(int particle, double E_min, double E_max)
{
  data::energy_min[particle - 1] = E_min;
  data::energy_max[particle - 1] = E_max;
}

} // namespace openmc
