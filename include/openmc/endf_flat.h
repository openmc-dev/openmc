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
  // Constructors
  explicit Function1DFlat(const uint8_t* data) : data_(data) { }

  #pragma omp declare target
  double operator()(double x) const;
  #pragma omp end declare target

  FunctionType type() const;
private:
  // Data members
  const uint8_t* data_;
};

class Function1DFlatContainer {
public:
  explicit Function1DFlatContainer(const Function1D& func);

  double operator()(double x) const;

  const uint8_t* data() const { return buffer_.data_; }
  FunctionType type() const { return this->func().type(); }
  Function1DFlat func() const { return Function1DFlat(buffer_.data_); }

private:
  DataBuffer buffer_;
};

} // namespace openmc

#endif // OPENMC_ENDF_FLAT_H
