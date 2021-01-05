#include "openmc/serialize.h"

#include <cstring> // for memcpy

namespace openmc {

DataBuffer::DataBuffer(size_t n)
{
  data_ = std::make_unique<uint8_t[]>(n);
}

} // end namespace openmc
