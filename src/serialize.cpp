#include "openmc/serialize.h"

#include <cstring> // for memcpy

namespace openmc {

DataBuffer::DataBuffer(size_t n)
{
  this->reserve(n);
  mode_ = Mode::write;
}

void DataBuffer::reserve(size_t n)
{
  data_ = std::make_unique<uint8_t[]>(n);
  offset_ = 0;
}

} // end namespace openmc
