//! \file serialize.h
//! Unified angle-energy distribution

#ifndef OPENMC_SERIALIZE_H
#define OPENMC_SERIALIZE_H

#include <cstdint> // for uint8_t
#include <memory> // for unique_ptr
#include <vector>

#include <xtensor/xtensor.hpp>

namespace openmc {

class DataBuffer {
public:
  explicit DataBuffer(size_t n);

  void add(int value);
  void add(double value);
  void add(size_t value);
  void add(const std::vector<double>& value);
  void add(const std::vector<int>& value);
  void add(const xt::xtensor<double, 1>& value);

  std::unique_ptr<uint8_t[]> data_;
  size_t offset_{0};
};

} // namespace openmc

#endif // OPENMC_SERIALIZE_H
