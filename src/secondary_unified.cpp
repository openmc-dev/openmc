#include "openmc/secondary_unified.h"

#include "openmc/secondary_kalbach.h"
#include "openmc/secondary_nbody.h"
#include "openmc/secondary_uncorrelated.h"

namespace openmc {

UnifiedAngleEnergy::UnifiedAngleEnergy(AngleEnergyType type, DataBuffer buffer)
  : type_(type), buffer_(std::move(buffer)) { }

void
UnifiedAngleEnergy::sample(double E_in, double& E_out, double& mu, uint64_t* seed) const
{
  switch (type_) {
  case AngleEnergyType::UNCORRELATED:
    {
      UncorrelatedAngleEnergyFlat dist(this->data());
      dist.sample(E_in, E_out, mu, seed);
    }
    break;
  case AngleEnergyType::KALBACH_MANN:
    {
      KalbachMannFlat dist(this->data());
      dist.sample(E_in, E_out, mu, seed);
    }
    break;
  case AngleEnergyType::CORRELATED:
    break;
  case AngleEnergyType::NBODY:
    {
      NBodyPhaseSpaceFlat dist(this->data());
      dist.sample(E_in, E_out, mu, seed);
    }
    break;
  }
}

} // namespace openmc
