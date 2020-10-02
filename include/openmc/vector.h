#ifndef OPENMC_VECTOR_H
#define OPENMC_VECTOR_H

/*
 * In an implementation of OpenMC that offloads computations to an accelerator,
 * we may need to provide replacements for standard library containers and
 * algorithms that have no native implementations on the device of interest.
 * Because some developers are currently in the process of creating such code,
 * introducing the below typedef lessens the amount of rebase conflicts that
 * happen as they rebase their code on OpenMC's develop branch.
 *
 */

#include "openmc/memory.h"
#include <algorithm>
#include <cstddef>
#include <type_traits>
#include <utility>
#include <vector>

namespace openmc {

/**
 * This is an "implementation" of a C++ standard vector,
 * except for how it's actually not.
 *
 * Importantly, this uses CUDA unified memory. It's set up so
 * that all memory control should be done from the host, and
 * provides methods for data access that work on the device.
 * I'm not aware of any CUDA libraries which provide something
 * like this. It makes OpenMC GPU compatibility fit smoothly.
 */
#ifdef __CUDACC__

template<class T, class Alloc = UnifiedAllocator<T>>
class vector {
private:
  T* begin;
  std::size_t length;
  std::size_t capacity;
  Alloc alloc;

  __host__ void grow()
  {
    T* old_begin = begin;
    auto old_capacity = capacity;
    capacity *= 2;
    begin = alloc.allocate(length);
    memcpy(begin, old_begin, length * sizeof(T));
    alloc.deallocate(old_begin, old_capacity);
  }

public:
  __host__ vector() : begin(nullptr), length(0), capacity(0) {}

  // Construct, default-initializing elements in a vector of
  // length new_length.
  __host__ vector(std::size_t new_length, const T& value = T())
  {
    begin = alloc.allocate(new_length);
    capacity = new_length;
    for (std::size_t i = 0; i < new_length; ++i)
      new (begin + i) T(value);
    length = new_length;
  }

  __host__ void reserve(std::size_t new_length)
  {
    T* old_begin = begin;
    auto old_capacity = capacity;
    capacity = new_length;
    begin = alloc.allocate(capacity);
    memcpy(begin, old_begin, length * sizeof(T));
    alloc.deallocate(old_begin, old_capacity);
  }

  __host__ void push_back(const T& value)
  {
    if (length == capacity)
      grow();
    begin[length++] = value;
  }

  __host__ __device__ T& operator[](std::size_t n) { return *(begin + n); }
  __host__ __device__ T const& operator[](std::size_t n) const
  {
    return *(begin + n);
  }

  __host__ __device__ bool empty() { return length > 0; }
  __host__ __device__ T* data() { return begin; }
  __host__ __device__ std::size_t size() { return length; }
};

#else

template<class T, class Alloc>
using vector = std::vector<T, Alloc>;

#endif

} // namespace openmc
#endif // OPENMC_VECTOR_H
