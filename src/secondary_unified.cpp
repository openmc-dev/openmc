#include "openmc/secondary_unified.h"

namespace openmc {

DataBuffer::DataBuffer(size_t n)
{
  data_ = std::make_unique<uint8_t[]>(n);
}

void DataBuffer::add(int value)
{
  auto data_int = reinterpret_cast<int*>(data_.get() + offset_);
  *data_int = value;
  offset_ += sizeof(int);
}

void DataBuffer::add(double value)
{
  auto data_int = reinterpret_cast<double*>(data_.get() + offset_);
  *data_int = value;
  offset_ += sizeof(double);
}

UnifiedAngleEnergy::UnifiedAngleEnergy(AngleEnergyType type, DataBuffer buffer)
  : type_(type), buffer_(std::move(buffer)) { }

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
