#include "openmc/secondary_flat.h"

#include "openmc/secondary_correlated.h"
#include "openmc/secondary_kalbach.h"
#include "openmc/secondary_nbody.h"
#include "openmc/secondary_thermal.h"
#include "openmc/secondary_uncorrelated.h"

namespace openmc {

void
AngleEnergyFlat::sample(double E_in, double& E_out, double& mu, uint64_t* seed) const
{
  switch (this->type()) {
  case AngleEnergyType::UNCORRELATED:
    {
      UncorrelatedAngleEnergyFlat dist(data_);
      dist.sample(E_in, E_out, mu, seed);
    }
    break;
  case AngleEnergyType::KALBACH_MANN:
    {
      KalbachMannFlat dist(data_);
      dist.sample(E_in, E_out, mu, seed);
    }
    break;
  case AngleEnergyType::CORRELATED:
    {
      CorrelatedAngleEnergyFlat dist(data_);
      dist.sample(E_in, E_out, mu, seed);
    }
    break;
  case AngleEnergyType::NBODY:
    {
      NBodyPhaseSpaceFlat dist(data_);
      dist.sample(E_in, E_out, mu, seed);
    }
    break;
  case AngleEnergyType::COHERENT_ELASTIC:
    {
      CoherentElasticAEFlat dist(data_);
      dist.sample(E_in, E_out, mu, seed);
    }
    break;
  case AngleEnergyType::INCOHERENT_ELASTIC:
    {
      IncoherentElasticAEFlat dist(data_);
      dist.sample(E_in, E_out, mu, seed);
    }
    break;
  case AngleEnergyType::INCOHERENT_ELASTIC_DISCRETE:
    {
      IncoherentElasticAEDiscreteFlat dist(data_);
      dist.sample(E_in, E_out, mu, seed);
    }
    break;
  case AngleEnergyType::INCOHERENT_INELASTIC:
    {
      IncoherentInelasticAEFlat dist(data_);
      dist.sample(E_in, E_out, mu, seed);
    }
    break;
  case AngleEnergyType::INCOHERENT_INELASTIC_DISCRETE:
    {
      IncoherentInelasticAEDiscreteFlat dist(data_);
      dist.sample(E_in, E_out, mu, seed);
    }
    break;
  }
}

AngleEnergyType AngleEnergyFlat::type() const
{
  int value = *reinterpret_cast<const int*>(data_);
  return static_cast<AngleEnergyType>(value);
}

AngleEnergyFlatContainer::AngleEnergyFlatContainer(const AngleEnergy& dist)
{
  // Determine number of bytes needed and create allocation
  size_t n = buffer_nbytes(dist);

  // Write into buffer
  buffer_.reserve(n);
  dist.serialize(buffer_);
  Ensures(n == buffer_.size());
}

void AngleEnergyFlatContainer::sample(double E_in, double& E_out, double& mu, uint64_t* seed) const
{
  this->dist().sample(E_in, E_out, mu, seed);
}

void AngleEnergyFlatContainer::copy_to_device() const
{
  buffer_.copy_to_device();
}

void AngleEnergyFlatContainer::release_device() const
{
  buffer_.release_device();
}

} // namespace openmc
