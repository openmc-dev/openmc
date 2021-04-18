//! \file endf_flat.h
//! Flattened classes and functions related to the ENDF-6 format

#ifndef OPENMC_ENDF_FLAT_H
#define OPENMC_ENDF_FLAT_H

#include "openmc/endf.h"
#include "openmc/serialize.h"

namespace openmc {

enum class FunctionType {
  POLYNOMIAL,
  TABULATED,
  COHERENT_ELASTIC,
  INCOHERENT_ELASTIC
};

class Function1DFlat {
public:
  explicit Function1DFlat(const Function1D& func);
  explicit Function1DFlat(DataBuffer buffer);

  double operator()(double x) const;

  const uint8_t* data() const { return buffer_.data_; }
  FunctionType type() const { return type_; }
private:
  // Data members
  FunctionType type_;
  DataBuffer buffer_;
};

} // namespace openmc

#endif // OPENMC_ENDF_FLAT_H
