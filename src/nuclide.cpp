#include "openmc/nuclide.h"

#include "openmc/hdf5_interface.h"

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

namespace data {

std::array<double, 2> energy_min {0.0, 0.0};
std::array<double, 2> energy_max {INFTY, INFTY};

std::vector<Nuclide> nuclides;

} // namespace data

//==============================================================================
// Nuclide implementation
//==============================================================================

Nuclide::Nuclide(hid_t group)
{
  // Get name of nuclide from group, removing leading '/'
  name_ = object_name(group).substr(1);

  read_attribute(group, "Z", Z_);
  read_attribute(group, "A", A_);
  read_attribute(group, "metastable", metastable_);
  read_attribute(group, "atomic_weight_ratio", awr_);
}

//==============================================================================
// Fortran compatibility functions
//==============================================================================

extern "C" void
set_particle_energy_bounds(int particle, double E_min, double E_max)
{
  data::energy_min[particle - 1] = E_min;
  data::energy_max[particle - 1] = E_max;
}

extern "C" void nuclide_from_hdf5_c(hid_t group)
{
  data::nuclides.emplace_back(group);
}

} // namespace openmc
