//! \file endf_flat.h
//! Flattened classes and functions related to the ENDF-6 format

#ifndef OPENMC_ENDF_FLAT_H
#define OPENMC_ENDF_FLAT_H

#include "openmc/endf.h"
#include "openmc/serialize.h"

namespace openmc {

enum FunctionType {
  POLYNOMIAL,
  TABULATED,
  COHERENT_ELASTIC,
  INCOHERENT_ELASTIC
};

class Function1DFlat : public Function1D {
public:
  Function1DFlat(DataBuffer buffer);

  double operator()(double x) const override;

  const uint8_t* data() const { return buffer_.data_.get(); }
  FunctionType type() const { return type_; }
private:
  // Data members
  FunctionType type_;
  DataBuffer buffer_;
};

} // namespace openmc

#endif // OPENMC_ENDF_FLAT_H
