#include "openmc/secondary_uncorrelated.h"

#include <string> // for string

#include <fmt/core.h>

#include "openmc/error.h"
#include "openmc/hdf5_interface.h"
#include "openmc/nuclide.h"
#include "openmc/particle.h"
#include "openmc/random_dist.h"

namespace openmc {

//==============================================================================
// UncorrelatedAngleEnergy implementation
//==============================================================================

UncorrelatedAngleEnergy::UncorrelatedAngleEnergy(hid_t group)
{
  // Check if angle group is present & read
  if (object_exists(group, "angle")) {
    hid_t angle_group = open_group(group, "angle");
    angle_ = AngleDistribution {angle_group};
    close_group(angle_group);
  }

  // Check if energy group is present & read
  if (object_exists(group, "energy")) {
    hid_t energy_group = open_group(group, "energy");

    std::string type;
    read_attribute(energy_group, "type", type);
    using UPtrEDist = unique_ptr<EnergyDistribution>;
    if (type == "discrete_photon") {
      energy_ = UPtrEDist {new DiscretePhoton {energy_group}};
    } else if (type == "level") {
      energy_ = UPtrEDist {new LevelInelastic {energy_group}};
    } else if (type == "continuous") {
      energy_ = UPtrEDist {new ContinuousTabular {energy_group}};
    } else if (type == "maxwell") {
      energy_ = UPtrEDist {new MaxwellEnergy {energy_group}};
    } else if (type == "evaporation") {
      energy_ = UPtrEDist {new Evaporation {energy_group}};
    } else if (type == "watt") {
      energy_ = UPtrEDist {new WattEnergy {energy_group}};
    } else {
      warning(
        fmt::format("Energy distribution type '{}' not implemented.", type));
    }
    close_group(energy_group);
  }
}

void UncorrelatedAngleEnergy::sample(
  double E_in, double& E_out, double& mu, uint64_t* seed) const
{
  // Sample cosine of scattering angle
  if (!angle_.empty()) {
    mu = angle_.sample(E_in, seed);
  } else {
    // no angle distribution given => assume isotropic for all energies
    mu = uniform_distribution(-1., 1., seed);
  }

  // Sample outgoing energy
  E_out = energy_->sample(E_in, seed);
}

double UncorrelatedAngleEnergy::get_pdf(
  double E_in, double mu, double& E_out, uint64_t* seed) const
{
  // Sample outgoing energy
  E_out = energy_->sample(E_in, seed);
}

} // namespace openmc
