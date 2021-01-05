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
  enum Mode {
    count,
    write
  };

  DataBuffer() = default;
  explicit DataBuffer(size_t n);

  void reserve(size_t n);
  size_t size() const { return offset_; }

  template<typename T> std::enable_if_t<std::is_scalar<std::decay_t<T>>::value>
  add(T value);

  template<typename T> void add(const std::vector<T>& value);
  template<typename T> void add(const xt::xtensor<T, 1>& value);

  std::unique_ptr<uint8_t[]> data_;
  size_t offset_{0};
  Mode mode_{Mode::write};
};

template<typename T> inline
std::enable_if_t<std::is_scalar<std::decay_t<T>>::value>
DataBuffer::add(T value)
{
  if (mode_ == Mode::write) {
    auto ptr = reinterpret_cast<T*>(data_.get() + offset_);
    *ptr = value;
  }
  offset_ += sizeof(T);
}

template<typename T> inline
void DataBuffer::add(const std::vector<T>& value)
{
  if (mode_ == Mode::write) {
    std::memcpy(data_.get() + offset_, value.data(), sizeof(T)*value.size());
  }
  offset_ += sizeof(T)*value.size();
}

template<typename T> inline
void DataBuffer::add(const xt::xtensor<T, 1>& value)
{
  if (mode_ == Mode::write) {
    std::memcpy(data_.get() + offset_, value.data(), sizeof(T)*value.size());
  }
  offset_ += sizeof(T)*value.size();
}

template<typename T>
size_t buffer_nbytes(const T& obj)
{
  DataBuffer buffer;
  buffer.mode_ = DataBuffer::Mode::count;
  obj.serialize(buffer);
  return buffer.size();
}

} // namespace openmc

#endif // OPENMC_SERIALIZE_H
