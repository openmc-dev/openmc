#pragma once
#include <array>

#ifdef __CUDACC__
template<typename T, unsigned long Size>
class array {
public:
  array() = default;
  __host__ __device__ T& operator[](unsigned indx) { return data[indx]; }
  __host__ __device__ T const& operator[](unsigned indx) const
  {
    return data[indx];
  }

  __host__ __device__ array<T, Size>& operator=(array<T, Size> const& other)
  {
    for (unsigned i = 0; i < Size; ++i) {
      data[i] = other.data[i];
    }
    return *this;
  }
  __host__ __device__ array(array<T, Size> const& other)
  {
    for (unsigned i = 0; i < Size; ++i) {
      data[i] = other.data[i];
    }
  }

  __host__ __device__ array(std::initializer_list<T> list)
  {
    unsigned i = 0;
    for (auto const& val : list)
      data[i++] = val;
  }

  __host__ operator std::array<T, Size>() const
  {
    std::array<T, Size> result;
    for (unsigned i = 0; i < Size; ++i) {
      result[i] = data[i];
    }
    return result;
  }

private:
  T data[Size];
};
#else
using std::array;
#endif
