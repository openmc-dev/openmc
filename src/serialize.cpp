#include "openmc/serialize.h"

#include <cstring> // for memcpy

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
  auto data_double = reinterpret_cast<double*>(data_.get() + offset_);
  *data_double = value;
  offset_ += sizeof(double);
}

void DataBuffer::add(size_t value)
{
  auto data_cast = reinterpret_cast<size_t*>(data_.get() + offset_);
  *data_cast = value;
  offset_ += sizeof(size_t);

}

void DataBuffer::add(const std::vector<int>& value)
{
  auto data_int = reinterpret_cast<int*>(data_.get() + offset_);
  std::memcpy(data_int, value.data(), value.size());
}

void DataBuffer::add(const std::vector<double>& value)
{
  auto data_double = reinterpret_cast<int*>(data_.get() + offset_);
  std::memcpy(data_double, value.data(), value.size());
}

void DataBuffer::add(const xt::xtensor<double, 1>& value)
{
  auto data_double = reinterpret_cast<int*>(data_.get() + offset_);
  std::memcpy(data_double, value.data(), value.size());
}


} // end namespace openmc
