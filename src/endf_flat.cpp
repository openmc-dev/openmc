#include "openmc/endf_flat.h"

#include "openmc/endf.h"
#include "openmc/error.h"

namespace openmc {

double Function1DFlat::operator()(double x) const
{
  switch (this->type()) {
  case FunctionType::TABULATED:
    {
      Tabulated1DFlat dist(data_);
      return dist(x);
    }
  case FunctionType::POLYNOMIAL:
    {
      PolynomialFlat dist(data_);
      return dist(x);
    }
  case FunctionType::COHERENT_ELASTIC:
    {
      CoherentElasticXSFlat dist(data_);
      return dist(x);
    }
  case FunctionType::INCOHERENT_ELASTIC:
    {
      IncoherentElasticXSFlat dist(data_);
      return dist(x);
    }
  default:
    UNREACHABLE();
  }
}

FunctionType Function1DFlat::type() const
{
  int value = *reinterpret_cast<const int*>(data_);
  return static_cast<FunctionType>(value);
}


Function1DFlatContainer::Function1DFlatContainer(const Function1D& func)
{
  // Determine number of bytes needed and create allocation
  size_t n = buffer_nbytes(func);

  // Write into buffer
  buffer_.reserve(n);
  func.serialize(buffer_);
  Ensures(n == buffer_.size());
}

double Function1DFlatContainer::operator()(double x) const
{
  return this->func()(x);
}

void Function1DFlatContainer::copy_to_device()
{
  buffer_.copy_to_device();
}

void Function1DFlatContainer::release_from_device()
{
  buffer_.release_device();
}

} // namespace openmc
