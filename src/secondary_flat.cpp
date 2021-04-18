#include "openmc/secondary_flat.h"

#include "openmc/secondary_correlated.h"
#include "openmc/secondary_kalbach.h"
#include "openmc/secondary_nbody.h"
#include "openmc/secondary_thermal.h"
#include "openmc/secondary_uncorrelated.h"

namespace openmc {

AngleEnergyFlat::AngleEnergyFlat(const AngleEnergy& dist)
{
  // Determine number of bytes needed and create allocation
  size_t n = buffer_nbytes(dist);

  // Write into buffer
  buffer_.reserve(n);
  dist.serialize(buffer_);
  Ensures(n == buffer_.size());

  int value = *reinterpret_cast<const int*>(this->data());
  type_ = static_cast<AngleEnergyType>(value);
}

AngleEnergyFlat::AngleEnergyFlat(DataBuffer buffer)
  : buffer_(std::move(buffer))
{
  int value = *reinterpret_cast<const int*>(this->data());
  type_ = static_cast<AngleEnergyType>(value);
}

void AngleEnergyFlat::copy_to_device() const
{
  buffer_.copy_to_device();
}

void AngleEnergyFlat::release_device() const
{
  buffer_.release_device();
}

void
AngleEnergyFlat::sample(double E_in, double& E_out, double& mu, uint64_t* seed) const
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
    {
      CorrelatedAngleEnergyFlat dist(this->data());
      dist.sample(E_in, E_out, mu, seed);
    }
    break;
  case AngleEnergyType::NBODY:
    {
      NBodyPhaseSpaceFlat dist(this->data());
      dist.sample(E_in, E_out, mu, seed);
    }
    break;
  case AngleEnergyType::COHERENT_ELASTIC:
    {
      CoherentElasticAEFlat dist(this->data());
      dist.sample(E_in, E_out, mu, seed);
    }
    break;
  case AngleEnergyType::INCOHERENT_ELASTIC:
    {
      IncoherentElasticAEFlat dist(this->data());
      dist.sample(E_in, E_out, mu, seed);
    }
    break;
  case AngleEnergyType::INCOHERENT_ELASTIC_DISCRETE:
    {
      IncoherentElasticAEDiscreteFlat dist(this->data());
      dist.sample(E_in, E_out, mu, seed);
    }
    break;
  case AngleEnergyType::INCOHERENT_INELASTIC:
    {
      IncoherentInelasticAEFlat dist(this->data());
      dist.sample(E_in, E_out, mu, seed);
    }
    break;
  case AngleEnergyType::INCOHERENT_INELASTIC_DISCRETE:
    {
      IncoherentInelasticAEDiscreteFlat dist(this->data());
      dist.sample(E_in, E_out, mu, seed);
    }
    break;
  }
}

} // namespace openmc
