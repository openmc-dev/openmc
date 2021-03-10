#include "openmc/serialize.h"

#include <cstring> // for memcpy

namespace openmc {

DataBuffer::DataBuffer(size_t n)
{
  this->reserve(n);
  mode_ = Mode::write;
}

DataBuffer::~DataBuffer()
{
  if (data_) delete[] data_;
}

void DataBuffer::reserve(size_t n)
{
  if (data_) delete[] data_;
  data_ = new uint8_t[n];
  offset_ = 0;
}

} // end namespace openmc
