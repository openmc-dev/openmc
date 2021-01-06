#include "openmc/endf_flat.h"

#include "openmc/endf.h"
#include "openmc/error.h"

namespace openmc {

Function1DFlat::Function1DFlat(DataBuffer buffer)
  : buffer_(std::move(buffer))
{
  int value = *reinterpret_cast<const int*>(this->data());
  type_ = static_cast<FunctionType>(value);
}

double Function1DFlat::operator()(double x) const
{
  switch (type_) {
  case FunctionType::TABULATED:
    {
      Tabulated1DFlat dist(this->data());
      return dist(x);
    }
  case FunctionType::POLYNOMIAL:
    {
      PolynomialFlat dist(this->data());
      return dist(x);
    }
  case FunctionType::COHERENT_ELASTIC:
    {
      CoherentElasticXSFlat dist(this->data());
      return dist(x);
    }
  case FunctionType::INCOHERENT_ELASTIC:
    {
      IncoherentElasticXSFlat dist(this->data());
      return dist(x);
    }
  default:
    UNREACHABLE();
  }
}

} // namespace openmc
