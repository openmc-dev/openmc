#pragma once
#include <array>

#ifdef __CUDACC__
template<typename T, unsigned long Size>
class array {
public:
  array() = default;
  using value_type = T;
  __host__ __device__ T& operator[](unsigned indx) { return data_[indx]; }
  __host__ __device__ T const& operator[](unsigned indx) const
  {
    return data_[indx];
  }

  __host__ __device__ array<T, Size>& operator=(array<T, Size> const& other)
  {
    for (unsigned i = 0; i < Size; ++i) {
      data_[i] = other.data_[i];
    }
    return *this;
  }
  __host__ __device__ array(array<T, Size> const& other)
  {
    for (unsigned i = 0; i < Size; ++i) {
      data_[i] = other.data_[i];
    }
  }

  __host__ __device__ void fill(T const& val)
  {
    for (unsigned i = 0; i < Size; ++i) {
      data_[i] = val;
    }
  }

  __host__ __device__ array(std::initializer_list<T> list)
  {
    unsigned i = 0;
    for (auto const& val : list)
      data_[i++] = val;
  }

  __host__ operator std::array<T, Size>() const
  {
    std::array<T, Size> result;
    for (unsigned i = 0; i < Size; ++i) {
      result[i] = data_[i];
    }
    return result;
  }

  __host__ __device__ T* data()
  {
    T* start = &data_[0];
    return start;
  }

private:
  T data_[Size];
};
#else
using std::array;
#endif
