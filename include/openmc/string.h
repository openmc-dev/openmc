#pragma once

#include "openmc/vector.h"
#include <string>

namespace openmc {

// Use a custom-made standard library string with GPU-compatible methods if in
// CUDA
#ifdef __CUDACC__

// By far, this is not a very efficient or perhaps even nonconforming
// implementation of a standard library string. However, string-handling in
// OpenMC is probably 0.0000001% of the work done, and we really just need
// something that compiles and works in this case.
class string {
private:
  vector<char> data_;

public:
  // Let default constructor work on host or device
  __host__ __device__ string() {}

  // Plain old rvalue constructor
  __host__ __device__ string& operator=(string&& copy_from)
  {
    data_ = std::move(copy_from.data_);
    return *this;
  }
  __host__ __device__ string(string&& copy_from)
    : data_(std::move(copy_from.data_))
  {}

  // Convert from std::string to CUDA-compat string
  __host__ string& operator=(std::string const& copy_from)
  {
    data_.resize(copy_from.size());
    for (unsigned i = 0; i < copy_from.size(); i++)
      data_[i] = copy_from[i];
    return *this;
  }

  // Convert from CUDA-compat string to std::string
  __host__ operator std::string() const
  {
    std::string tmp(data_.begin(), data_.end());
    return tmp;
  }

  __host__ __device__ char* begin() const { return data_.begin(); }
  __host__ __device__ char* end() const { return data_.end(); }
  __host__ __device__ bool empty() const { return data_.empty(); }
  __host__ __device__ std::size_t length() const { return data_.size(); }
  __host__ __device__ char* c_str() { return data_.begin(); }
};

#else

using std::string;

#endif
} // namespace openmc
