#pragma once

#include "openmc/vector.h"
#include <algorithm> // lexicographical_compare
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

  // Construction from std::string
  __host__ string(std::string const& copy_from)
    : data_(copy_from.c_str(), copy_from.c_str() + copy_from.size())
  {}
  __host__ string(string const& copy_from) : data_(copy_from.data_) {}
  __host__ string& operator=(string const& copy_from)
  {
    data_ = copy_from.data_;
    return *this;
  }

  // Convert from std::string to CUDA-compat string
  __host__ string& operator=(std::string const& copy_from)
  {
    data_.resize(copy_from.size());
    for (unsigned i = 0; i < copy_from.size(); i++)
      data_[i] = copy_from[i];
    return *this;
  }

  // Assignment from C string
  __host__ string& operator=(const char* data)
  {
    std::size_t j = 0;
    while (data[j] != '\0')
      j++;
    data_.resize(j + 1); // include null char
    for (std::size_t i = 0; i <= j; ++i)
      data_[i] = data[i];
    return *this;
  }

  // Substring constructor
  __host__ string(const char* data, std::size_t length)
    : data_(data, data + length)
  {}

  // Null-terminated constructor
  __host__ string(const char* data)
  {
    // Loop through to resize memory
    std::size_t j = 0;
    while (data[j] != '\0')
      j++;
    data_.resize(j + 1); // include null char
    for (std::size_t i = 0; i <= j; ++i)
      data_[i] = data[i];
  }

  // Convert from CUDA-compat string to std::string. this allows seamless
  // integration with other parts of the codebase written to use std::string.
  __host__ operator std::string() const
  {
    std::string tmp(data_.begin(), data_.end());
    return tmp;
  }

  __host__ __device__ char* begin() const { return data_.begin(); }
  __host__ __device__ char* end() const { return data_.end(); }
  __host__ __device__ bool empty() const { return data_.empty(); }
  __host__ __device__ std::size_t size() const { return data_.size(); }
  __host__ __device__ const char* c_str() const { return data_.begin(); }
  __host__ __device__ char* data() const { return data_.begin(); }
};

inline __host__ bool operator<(string const& left, string const& right)
{
  return std::lexicographical_compare(
    left.begin(), left.end(), right.begin(), right.end());
}
inline __host__ bool operator==(string const& left, string const& right)
{
  if (left.size() != right.size())
    return false;
  // sure hope these are null terminated!!
  return std::strcmp(left.begin(), right.begin()) == 0;
}
inline __host__ bool operator==(string const& left, const char* right)
{
  return std::strcmp(left.begin(), right) == 0;
}
inline __host__ bool operator!=(string const& left, const char* right)
{
  return std::strcmp(left.begin(), right) != 0;
}

template<class CharT, class Traits>
inline std::basic_istream<CharT, Traits>& operator>>(
  std::basic_istream<CharT, Traits>& is, string& str)
{
  // it aint pretty, but it works
  std::string tmp;
  is >> tmp;
  str = tmp;
  return is;
}

template<class CharT, class Traits>
inline std::basic_ostream<CharT, Traits>& operator<<(
  std::basic_ostream<CharT, Traits>& os, string& str)
{
  // it aint pretty, but it works
  std::string tmp(str);
  os << tmp;
  return os;
}

#else

using std::string;

#endif
} // namespace openmc
