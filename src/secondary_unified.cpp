#include "openmc/secondary_unified.h"

namespace openmc {

UnifiedAngleEnergy::UnifiedAngleEnergy(AngleEnergyType type, std::unique_ptr<uint8_t[]> data)
  : type_(type), data_(std::move(data)) { }

void
UnifiedAngleEnergy::sample(double E_in, double& E_out, double& mu, uint64_t* seed) const
{
  switch (type_) {
  case AngleEnergyType::UNCORRELATED:
    break;
  case AngleEnergyType::KALBACH_MANN:
    break;
  case AngleEnergyType::CORRELATED:
    break;
  case AngleEnergyType::NBODY:
    break;
  }
}

} // namespace openmc
