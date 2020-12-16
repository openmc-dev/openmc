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
  std::memcpy(data_.get() + offset_, value.data(), 4*value.size());
  offset_ += 4*value.size();
}

void DataBuffer::add(const std::vector<double>& value)
{
  std::memcpy(data_.get() + offset_, value.data(), 8*value.size());
  offset_ += 8*value.size();
}

void DataBuffer::add(const xt::xtensor<double, 1>& value)
{
  std::memcpy(data_.get() + offset_, value.data(), 8*value.size());
  offset_ += 8*value.size();
}


} // end namespace openmc
